"""CLI entrance point."""

import argparse
from typing import Sequence

import numpy as np

from ._config_module import Config
from ._plotting_objects import Molecule
from .examples.get_example_files import all_examples
from .plotter import Plotter
from .tabulator import GridType, Tabulator


def _parse_orbital_expression(expr: str | None, total_mos: int, occupations: Sequence[int]) -> list[int]:
    """Parse a comma-separated orbital expression supporting ranges."""

    if not expr:
        if total_mos == 0:
            return []

        occupied = [i for i, occ in enumerate(occupations) if occ > 0]
        if not occupied:
            return [0]

        homo = occupied[-1]
        if homo + 1 < total_mos:
            return [homo, homo + 1]
        return [homo]

    selections: set[int] = set()
    cleaned = expr.replace(' ', '')
    for chunk in cleaned.split(','):
        if not chunk:
            continue
        if '-' in chunk or ':' in chunk:
            separator = '-' if '-' in chunk else ':'
            start_str, end_str = chunk.split(separator, 1)
            start = int(start_str)
            end = int(end_str)
            if start > end:
                start, end = end, start
            selections.update(range(start, end + 1))
        else:
            selections.add(int(chunk))

    indices = sorted(selections)
    if not indices:
        raise ValueError('No valid orbital indices were parsed from the expression.')
    if indices[0] < 0 or indices[-1] >= total_mos:
        raise ValueError(f'Orbital selection out of bounds. Valid range is 0-{total_mos - 1}.')
    return indices


def _parse_resolution(arg: str | None, defaults: Sequence[int]) -> tuple[int, int, int]:
    """Parse a resolution string of the form ``Nx,Ny,Nz``."""

    if not arg:
        return tuple(int(value) for value in defaults)

    parts = [part.strip() for part in arg.split(',') if part.strip()]
    if len(parts) != 3:
        raise ValueError('Resolution must provide exactly three comma-separated integers.')

    values = tuple(int(part) for part in parts)
    if any(value <= 0 for value in values):
        raise ValueError('Resolution values must be positive integers.')
    return values


def _build_grid(
    tabulator: Tabulator,
    grid_type: GridType,
    resolution: tuple[int, int, int],
    extent: float,
) -> None:
    """Configure the tabulator grid according to CLI options."""

    if grid_type is GridType.CARTESIAN:
        nx, ny, nz = resolution
        axis = np.linspace(-extent, extent, nx)
        tabulator.cartesian_grid(axis, np.linspace(-extent, extent, ny), np.linspace(-extent, extent, nz))
    elif grid_type is GridType.SPHERICAL:
        nr, ntheta, nphi = resolution
        tabulator.spherical_grid(
            np.linspace(0, extent, nr),
            np.linspace(0, np.pi, ntheta),
            np.linspace(0, 2 * np.pi, nphi),
        )
    else:  # pragma: no cover - defensive branch
        raise ValueError(f'Unsupported grid type: {grid_type}')


def main() -> None:
    """Entry point for the moldenViz command-line interface.

    Parses command line arguments and launches the plotter with the specified
    molden file or example molecule. Supports options to plot only the molecule
    structure without molecular orbitals.
    """
    parser = argparse.ArgumentParser(prog='moldenViz')
    source = parser.add_mutually_exclusive_group(required=True)

    source.add_argument('file', nargs='?', default=None, help='Optional molden file path', type=str)
    parser.add_argument('-m', '--only_molecule', action='store_true', help='Only plots the molecule')
    source.add_argument(
        '-e',
        '--example',
        type=str,
        metavar='molecule',
        choices=all_examples.keys(),
        help='Load example %(metavar)s. Options are: %(choices)s',
    )

    parser.add_argument('--export-vtk', metavar='PATH', help='Export selected orbitals to VTK structured grids.')
    parser.add_argument('--export-cube', metavar='PATH', help='Export selected orbitals to Gaussian cube files.')
    parser.add_argument(
        '--orbitals',
        metavar='INDICES',
        help='Comma-separated list or range (e.g. 15,16 or 12-18) specifying orbitals to export.',
    )
    parser.add_argument(
        '--grid',
        choices=[grid_type.value for grid_type in (GridType.SPHERICAL, GridType.CARTESIAN)],
        help='Grid type to use when exporting volumetric data.',
    )
    parser.add_argument(
        '--resolution',
        metavar='Nx,Ny,Nz',
        help='Override the grid resolution used during export.',
    )

    args = parser.parse_args()

    source_payload = args.file or all_examples[args.example]

    exports_requested = bool(args.export_vtk or args.export_cube)

    tabulator = None
    if exports_requested:
        if args.only_molecule:
            raise SystemExit('Orbital export is unavailable when --only-molecule is set.')

        tabulator = Tabulator(source_payload)
        config = Config()
        molecule = Molecule(tabulator._parser.atoms)  # noqa: SLF001

        extent = max(config.grid.max_radius_multiplier * molecule.max_radius, config.grid.min_radius)

        grid_choice = args.grid or ('cartesian' if args.export_cube else GridType.SPHERICAL.value)
        try:
            grid_value = GridType(grid_choice)
        except ValueError as exc:  # pragma: no cover - parser enforces choices but guard for safety
            raise ValueError(f'Unsupported grid choice: {grid_choice}') from exc

        if args.export_cube and grid_value is not GridType.CARTESIAN:
            grid_value = GridType.CARTESIAN

        if grid_value is GridType.CARTESIAN:
            defaults = (
                config.grid.cartesian.num_x_points,
                config.grid.cartesian.num_y_points,
                config.grid.cartesian.num_z_points,
            )
        else:
            defaults = (
                config.grid.spherical.num_r_points,
                config.grid.spherical.num_theta_points,
                config.grid.spherical.num_phi_points,
            )

        try:
            resolution = _parse_resolution(args.resolution, defaults)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc
        _build_grid(tabulator, grid_value, resolution, extent)

        total_mos = len(tabulator._parser.mos)  # noqa: SLF001
        occupations = [mo.occ for mo in tabulator._parser.mos]  # noqa: SLF001
        try:
            indices = _parse_orbital_expression(args.orbitals, total_mos, occupations)
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc

        if args.export_vtk:
            tabulator.export_vtk(args.export_vtk, indices)
        if args.export_cube:
            tabulator.export_cube(args.export_cube, indices)

    Plotter(
        source_payload,
        only_molecule=args.only_molecule,
        tabulator=tabulator,
    )


if __name__ == '__main__':
    main()
