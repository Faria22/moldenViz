"""Generate PySCF Molden fixtures from the bundled molecule geometries.

Run from the repository root with::

    uv run --with pyscf==2.14.0 python scripts/generate_pyscf_molden.py

The generated files use cc-pVQZ because every bundled molecule then contains
shells through ``l=4``. The orbitals diagonalize PySCF's one-electron core
Hamiltonian; an SCF calculation is intentionally unnecessary for parser and
tabulator fixtures. Only a window around the occupied/virtual boundary is
written to keep the files reasonably small; every retained orbital still has a
coefficient for every atomic orbital.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from pyscf import gto, lib, scf
from pyscf.tools import molden

from moldenViz import Parser

if TYPE_CHECKING:
    from collections.abc import Sequence

BASIS = 'cc-pvqz'
DEFAULT_OUTPUT_DIR = Path('tests/fixtures/pyscf')
EXAMPLE_DIR = Path('src/moldenViz/examples/molden_files')
LMAX = 4
O2_TRIPLET_SPIN = 2
REPRESENTATIONS = ('spherical', 'cartesian')


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'molecules',
        nargs='*',
        help='Example names to generate (default: every bundled example).',
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f'Output directory (default: {DEFAULT_OUTPUT_DIR}).',
    )
    parser.add_argument(
        '--representation',
        choices=(*REPRESENTATIONS, 'both'),
        default='both',
        help='Basis representation to write (default: both).',
    )
    parser.add_argument(
        '--max-orbitals',
        type=int,
        default=12,
        help='Maximum frontier orbitals written per file (default: 12).',
    )
    parser.add_argument(
        '--overwrite',
        action='store_true',
        help='Replace existing generated files.',
    )
    return parser.parse_args()


def _example_paths(names: Sequence[str]) -> list[Path]:
    available = {path.stem: path for path in sorted(EXAMPLE_DIR.glob('*.inp'))}
    if not names:
        return list(available.values())

    unknown = sorted(set(names) - available.keys())
    if unknown:
        choices = ', '.join(available)
        raise ValueError(f'Unknown examples: {", ".join(unknown)}. Choices: {choices}')
    return [available[name] for name in names]


def _build_molecule(example_path: Path, *, cartesian: bool) -> gto.Mole:
    parsed = Parser(str(example_path), only_molecule=True)
    molecule = gto.M(
        atom=[(atom.label, atom.position) for atom in parsed.atoms],
        basis=BASIS,
        cart=cartesian,
        spin=O2_TRIPLET_SPIN if example_path.stem == 'o2' else 0,
        unit='Bohr',
        verbose=lib.logger.NOTE,
    )

    angular_momenta = {molecule.bas_angular(index) for index in range(molecule.nbas)}
    if angular_momenta != set(range(LMAX + 1)):
        raise RuntimeError(
            f'{example_path.stem} with {BASIS} has angular momenta '
            f'{sorted(angular_momenta)}, expected 0 through {LMAX}.',
        )
    return molecule


def _frontier_indices(occupations: np.ndarray, max_orbitals: int) -> np.ndarray:
    if max_orbitals < 1:
        raise ValueError('--max-orbitals must be positive.')
    if len(occupations) <= max_orbitals:
        return np.arange(len(occupations))

    occupied = np.flatnonzero(occupations > 0)
    center = occupied[-1] + 1 if len(occupied) else 0
    start = max(0, center - max_orbitals // 2)
    start = min(start, len(occupations) - max_orbitals)
    return np.arange(start, start + max_orbitals)


def _write_molden(
    molecule: gto.Mole,
    coefficients: np.ndarray,
    energies: np.ndarray,
    occupations: np.ndarray,
    output_path: Path,
    *,
    max_orbitals: int,
) -> None:
    indices = _frontier_indices(occupations, max_orbitals)

    molden.from_mo(
        molecule,
        str(output_path),
        coefficients[:, indices],
        ene=energies[indices],
        occ=occupations[indices],
        ignore_h=False,
    )


def generate(
    example_path: Path,
    output_dir: Path,
    representation: str,
    *,
    max_orbitals: int,
    overwrite: bool,
) -> Path:
    """Generate one spherical- or Cartesian-basis Molden fixture.

    Returns
    -------
    Path
        Path to the generated or existing fixture.
    """
    cartesian = representation == 'cartesian'
    output_path = output_dir / f'{example_path.stem}-{BASIS}-{representation}.molden'
    if output_path.exists() and not overwrite:
        print(f'skip {output_path} (already exists)')  # ruff: ignore[print]
        return output_path

    molecule = _build_molecule(example_path, cartesian=cartesian)
    mean_field = scf.ROHF(molecule) if molecule.spin else scf.RHF(molecule)
    energies, coefficients = mean_field.eig(mean_field.get_hcore(), mean_field.get_ovlp())
    occupations = mean_field.get_occ(energies, coefficients)

    output_dir.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(f'{output_path.suffix}.tmp')
    _write_molden(
        molecule,
        np.asarray(coefficients),
        np.asarray(energies),
        np.asarray(occupations),
        temporary_path,
        max_orbitals=max_orbitals,
    )
    temporary_path.replace(output_path)
    print(  # ruff: ignore[print]
        f'wrote {output_path} (atoms={molecule.natm}, aos={molecule.nao_nr()}, lmax={LMAX})',
    )
    return output_path


def main() -> None:
    """Generate the requested PySCF Molden fixtures."""
    args = _parse_args()
    representations = REPRESENTATIONS if args.representation == 'both' else (args.representation,)
    for example_path in _example_paths(args.molecules):
        for representation in representations:
            generate(
                example_path,
                args.output_dir,
                representation,
                max_orbitals=args.max_orbitals,
                overwrite=args.overwrite,
            )


if __name__ == '__main__':
    main()
