"""Tabulator module for creating grids and tabulating Gaussian-type orbitals (GTOs) from Molden files."""

import logging
import os
from concurrent.futures import ThreadPoolExecutor
from enum import Enum
from math import factorial
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.special import assoc_legendre_p_all as s_plm

from .models import Atom, MolecularOrbital
from .parser import Parser

__all__ = ['GridType', 'Tabulator']

logger = logging.getLogger(__name__)

_MOIndices = NDArray[np.integer] | list[int] | tuple[int, ...] | range


def _grid_creation_with_only_molecule_error() -> RuntimeError:
    """Return a consistent error for grid creation when only the molecule is parsed.

    Returns
    -------
    RuntimeError
        The error indicating that grid creation is not allowed when `only_molecule` is set to `True`.

    """
    return RuntimeError('Grid creation is not allowed when `only_molecule` is set to `True`.')


class GridType(Enum):
    """Grid types allowed."""

    SPHERICAL = 'spherical'
    CARTESIAN = 'cartesian'
    UNKNOWN = 'unknown'


class Tabulator:
    """Extends Parser, create grids and tabulates Gaussian-type orbitals (GTOs) from Molden files.

    Parameters
    ----------
    source : str | list[str]
        The path to the molden file, or the lines from the file.
    only_molecule : bool, optional
        Only parse the atoms and skip molecular orbitals.
        Default is ``False``.

    Attributes
    ----------
    grid : NDArray[np.floating]
        The grid points where GTOs and MOs are tabulated.
    gtos : NDArray[np.floating]
        The cached Gaussian-type orbitals (GTOs) on the grid. Access raises
        ``RuntimeError`` when no cache is available.
    has_gtos : bool
        Whether GTO values are currently cached.
    """

    def __init__(
        self,
        source: str | list[str],
        only_molecule: bool = False,
    ) -> None:
        """Initialize the Tabulator with a Molden file or its content."""
        self._parser = Parser(source, only_molecule)

        self._only_molecule = only_molecule

        self._grid: NDArray[np.floating]
        self._grid_type = GridType.UNKNOWN
        self._grid_dimensions: tuple[int, int, int] = (0, 0, 0)

        # Used for when exporting to cube format
        self._grid_axes: tuple[NDArray[np.floating], ...] | None = None

    @property
    def grid(self) -> NDArray[np.floating]:
        """The 2D array of Cartesian grid points used for tabulation."""
        return self._grid

    def set_grid(self, new_grid: Any) -> None:
        """Set an arbitrary Cartesian point grid after validating its structure.

        Parameters
        ----------
        new_grid : Any
            Candidate 2D array whose rows are XYZ coordinates.

        Raises
        ------
        TypeError
            If ``new_grid`` is not a NumPy array.
        ValueError
            If the array does not have three columns or contains no rows.
        """
        min_num_rows = 1
        num_cols = 3
        num_dim = 2

        if not isinstance(new_grid, np.ndarray):
            raise TypeError(f"Expected a NumPy array for 'grid', but got {type(new_grid).__name__}.")

        if new_grid.ndim != num_dim:
            raise ValueError(f"'grid' must be a 2D array, but got shape {new_grid.shape}.")

        if new_grid.shape[0] < min_num_rows:
            raise ValueError("'grid' must have at least one row (one point in space).")

        if new_grid.shape[1] != num_cols:
            raise ValueError(f"'grid' must have exactly 3 columns, but got {new_grid.shape[1]} columns.")

        if self._only_molecule:
            raise _grid_creation_with_only_molecule_error()

        self._grid = new_grid
        self._grid_type = GridType.UNKNOWN
        self._grid_dimensions = (0, 0, 0)
        self._grid_axes = None
        self.clear_gtos()

    @property
    def has_gtos(self) -> bool:
        """Whether Gaussian-type orbitals are currently cached."""
        return hasattr(self, '_gtos')

    @property
    def gtos(self) -> NDArray[np.floating]:
        """The tabulated Gaussian-type orbitals (GTOs) on the grid.

        Raises
        ------
        RuntimeError
            If GTOs have not been tabulated or the cache was cleared.
        """
        if not self.has_gtos:
            raise RuntimeError('GTOs are not available. Call tabulate_gtos() first.')
        return self._gtos

    def clear_gtos(self) -> None:
        """Release any cached Gaussian-type orbital values."""
        if self.has_gtos:
            del self._gtos

    def set_gtos(self, gtos: NDArray[np.floating]) -> None:
        """Install GTO values produced for the current grid.

        Parameters
        ----------
        gtos : NDArray[np.floating]
            Two-dimensional GTO values with one row per grid point.

        Raises
        ------
        ValueError
            If the values are not two-dimensional or do not match the grid.
        """
        expected_dimensions = 2
        if gtos.ndim != expected_dimensions:
            raise ValueError('GTO data must be a two-dimensional array.')
        if gtos.shape[0] != self.grid.shape[0]:
            raise ValueError('GTO data must contain one row per grid point.')
        self._gtos = gtos

    @property
    def atoms(self) -> list[Atom]:
        """Atoms parsed from the Molden source."""
        return self._parser.atoms

    @property
    def molecular_orbitals(self) -> list[MolecularOrbital]:
        """Molecular-orbital metadata parsed from the Molden source."""
        return self._parser.mos

    @property
    def grid_type(self) -> GridType:
        """The coordinate grid type."""
        return self._grid_type

    @property
    def grid_dimensions(self) -> tuple[int, int, int]:
        """The three-dimensional grid shape used by structured exporters."""
        return self._grid_dimensions

    @property
    def grid_axes(self) -> tuple[NDArray[np.floating], ...] | None:
        """The original coordinate axes, if the grid was built from axes."""
        return self._grid_axes

    @staticmethod
    def _axis_spacing(axis: NDArray[np.floating], name: str) -> float:
        if axis.size <= 1:
            return 0.0

        diffs = np.diff(axis)
        if np.any(diffs <= 0):
            raise ValueError(f'{name}-axis values must be strictly increasing.')
        if not np.allclose(diffs, diffs[0]):
            raise ValueError(f'{name}-axis must be evenly spaced.')
        return float(diffs[0])

    @staticmethod
    def spherical_to_cartesian(
        r: NDArray[np.floating],
        theta: NDArray[np.floating],
        phi: NDArray[np.floating],
    ) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
        """Convert spherical coordinates to Cartesian coordinates.

        Returns
        -------
        tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]
            Tuple containing the Cartesian coordinates ``(x, y, z)``.
        """
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)
        return x, y, z

    @staticmethod
    def cartesian_to_spherical(
        x: NDArray[np.floating],
        y: NDArray[np.floating],
        z: NDArray[np.floating],
    ) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
        """Convert Cartesian coordinates to spherical coordinates.

        Returns
        -------
        tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]
            Tuple containing spherical coordinates ``(r, theta, phi)``.
        """
        r = np.sqrt(x * x + y * y + z * z)

        eps = np.finfo(float).eps
        safe_r = np.clip(r, eps, None)
        if np.isclose(safe_r, 0.0).all():
            return r, np.zeros_like(r), np.zeros_like(r)

        safe_ratio = Tabulator._check_bounds(z / safe_r)
        theta = np.arccos(safe_ratio)
        phi = np.arctan2(y, x)
        return r, theta, phi

    _spherical_to_cartesian = spherical_to_cartesian
    _cartesian_to_spherical = cartesian_to_spherical

    def _set_grid(
        self,
        x: NDArray[np.floating],
        y: NDArray[np.floating],
        z: NDArray[np.floating],
        grid_type: GridType,
        tabulate_gtos: bool = True,
    ) -> None:
        r"""Create grid from x, y, z (or r, theta, phi) arrays and tabulate GTOs.

        Parameters
        ----------
        x : NDArray[np.floating]
            1D array of x (or r) coordinates.
        y : NDArray[np.floating]
            1D array of y (or theta) coordinates.
        z : NDArray[np.floating]
            1D array of z (or phi) coordinates.
        grid_type : GridType
            What type of grid. Determines if x, y, z are actual
            x, y, z or r, theta, phi.
        tabulate_gtos : bool, optional
            Whether to tabulate Gaussian-type orbitals (GTOs) after creating the grid.
            Defaults to True.
        """
        if self._only_molecule:
            raise _grid_creation_with_only_molecule_error()

        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        if grid_type == GridType.SPHERICAL:
            xx, yy, zz = self.spherical_to_cartesian(xx, yy, zz)

        self.set_grid(np.column_stack((xx.ravel(), yy.ravel(), zz.ravel())))
        self._grid_type = grid_type
        self._grid_dimensions = (len(x), len(y), len(z))
        self._grid_axes = (x, y, z)
        logger.info(
            'Created %s grid with %d points (%dx%dx%d).',
            grid_type.value,
            self._grid.shape[0],
            len(x),
            len(y),
            len(z),
        )

        if tabulate_gtos:
            self._gtos = self.tabulate_gtos()

    def cartesian_grid(
        self,
        x: NDArray[np.floating],
        y: NDArray[np.floating],
        z: NDArray[np.floating],
        tabulate_gtos: bool = True,
    ) -> None:
        r"""Create cartesian grid from x, y, z arrays and tabulate GTOs.

        Parameters
        ----------
        x : NDArray[np.floating]
            1D array of x coordinates.
        y : NDArray[np.floating]
            1D array of y coordinates.
        z : NDArray[np.floating]
            1D array of z coordinates.
        tabulate_gtos : bool, optional
            Whether to tabulate Gaussian-type orbitals (GTOs) after creating the grid.
            Defaults to True.
        """
        logger.debug('Setting cartesian grid axes with lengths x=%d, y=%d, z=%d.', len(x), len(y), len(z))
        self._set_grid(x, y, z, GridType.CARTESIAN, tabulate_gtos)

    def spherical_grid(
        self,
        r: NDArray[np.floating],
        theta: NDArray[np.floating],
        phi: NDArray[np.floating],
        tabulate_gtos: bool = True,
    ) -> None:
        r"""Create spherical grid from r, theta, phi arrays and tabulate GTOs.

        Parameters
        ----------
        r : NDArray[np.floating]
            1D array of radial coordinates.
        theta : NDArray[np.floating]
            1D array of polar angles (radians).
        phi : NDArray[np.floating]
            1D array of azimuthal angles (radians).
        tabulate_gtos : bool, optional
            Whether to tabulate Gaussian-type orbitals (GTOs) after creating the grid.
            Defaults to True.

        Notes
        -----
        Grid points are converted to Cartesian coordinates.

        """
        logger.debug('Setting spherical grid axes with lengths r=%d, theta=%d, phi=%d.', len(r), len(theta), len(phi))
        self._set_grid(r, theta, phi, GridType.SPHERICAL, tabulate_gtos)

    def tabulate_gtos(self) -> NDArray[np.floating]:
        """Tabulate Gaussian-type orbitals (GTOs) on the current grid.

        Returns
        -------
        NDArray[np.floating]
            Array containing the tabulated GTOs data.

        Raises
        ------
        RuntimeError
            If the grid is not defined before tabulating GTOs,
            or if the `only_molecule` flag is set to `True`.
        """
        if self._only_molecule:
            raise RuntimeError('Grid creation is not allowed when `only_molecule` is set to `True`.')

        if not hasattr(self, 'grid'):
            raise RuntimeError('Grid is not defined. Please create a grid before tabulating GTOs.')

        total_points = self._grid.shape[0]
        total_coeffs = self._parser.mo_coeffs.shape[1]
        logger.info('Tabulating GTOs on %d grid points.', total_points)

        # Having a predefined array makes it faster to fill the data
        gto_data = np.empty((total_points, total_coeffs))
        atom_tasks: list[tuple[Any, slice]] = []
        idx_shell_start = 0

        # Calculate the slices for each atom's shells
        for atom in self._parser.atoms:
            num_gtos_in_shell = sum(2 * shell.l + 1 for shell in atom.shells)
            atom_slice = slice(idx_shell_start, idx_shell_start + num_gtos_in_shell)
            atom_tasks.append((atom, atom_slice))
            idx_shell_start += num_gtos_in_shell

        max_workers = min(len(atom_tasks), os.cpu_count() or 1)
        if max_workers <= 1:
            for atom, atom_slice in atom_tasks:
                self._tabulate_atom(atom, atom_slice, gto_data)
        else:
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [
                    executor.submit(self._tabulate_atom, atom, atom_slice, gto_data) for atom, atom_slice in atom_tasks
                ]
                for future in futures:
                    future.result()

        logger.debug('GTO data shape: %s', gto_data.shape)

        self.set_gtos(gto_data)
        return gto_data

    def _tabulate_atom(self, atom: Any, atom_slice: slice, gto_data: NDArray[np.floating]) -> None:
        """Tabulate all shells for a single atom into the shared GTO array."""
        centered_grid = self._grid - atom.position
        max_l = atom.shells[-1].l
        r_sq = np.einsum('ij,ij->i', centered_grid, centered_grid)
        solid_harmonics = self._tabulate_real_solid_harmonics(centered_grid, max_l)
        atom_block = gto_data[:, atom_slice]
        block_cursor = 0

        for shell in atom.shells:
            l = shell.l
            num_m = 2 * l + 1
            m_inds = np.arange(-l, l + 1)
            inner_slice = slice(block_cursor, block_cursor + num_m)

            contraction = shell._prefactor @ np.exp(  # ruff:ignore[private-member-access]
                -shell._gto_exps[:, None] * r_sq[None, :],  # ruff:ignore[private-member-access]
            )

            atom_block[:, inner_slice] = contraction[:, None] * solid_harmonics[l, m_inds, ...].T
            block_cursor += num_m

    def tabulate_mos(self, mo_inds: int | _MOIndices | None = None) -> NDArray[np.floating]:
        """Tabulate molecular orbitals (MOs) on the current grid.

        Parameters
        ----------
        mo_inds : int, array-like, or None, optional
            Indices of the MOs to tabulate. If None, all MOs are tabulated.

        Returns
        -------
        NDArray[np.floating]
            Array containing the tabulated MOs data.

            If an integer is provided, it will tabulate only that MO.
            If an array-like is provided, it will tabulate the MOs at those indices.

        Raises
        ------
        RuntimeError
            If the grid is not defined before tabulating MOs.
        RuntimeError
            If GTOs are not tabulated before tabulating MOs.
        ValueError
            If provided mo_inds is invalid.
        """
        if not hasattr(self, '_grid'):
            raise RuntimeError('Grid is not defined. Please create a grid before tabulating MOs.')
        if not self.has_gtos:
            raise RuntimeError('GTOs are not tabulated. Please tabulate GTOs before tabulating MOs.')

        if mo_inds is None:
            mo_inds = list(range(len(self._parser.mos)))

        num_requested = 1 if isinstance(mo_inds, int) else len(mo_inds)
        logger.info('Tabulating %d molecular orbital(s).', num_requested)

        if isinstance(mo_inds, range):
            mo_inds = list(mo_inds)

        if not isinstance(mo_inds, int) and not mo_inds:
            raise ValueError('Provided mo_inds is empty. Please provide valid indices.')

        if isinstance(mo_inds, int):
            if mo_inds < 0 or mo_inds >= len(self._parser.mos):
                raise ValueError('Provided mo_index is invalid. Please provide valid index.')
        elif any(mo_ind < 0 or mo_ind >= len(self._parser.mos) for mo_ind in mo_inds):
            raise ValueError('Provided mo_inds contains invalid indices. Please provide valid indices.')

        if isinstance(mo_inds, int):
            mo_data = self.gtos @ self._parser.mo_coeffs[mo_inds]
        else:
            mo_coeffs = self._parser.mo_coeffs[mo_inds]
            mo_data = self.gtos @ mo_coeffs.T
            logger.debug('MO data shape: %s', mo_data.shape)

        return mo_data

    def export(self, path: str | Path, *, mo_index: int | None = None) -> None:
        """Export the current grid data to a VTK-based or cube file.

        Parameters
        ----------
        path : str | Path
            Target path for the exported data. The file extension should
            match the desired exporter (``.vtk`` for VTK,
            ``.cube`` for cube files).
        mo_index : int | None, optional
            Molecular orbital index to export.
            The parameter is optional for VTK exports. If none is given
            then all the molecular orbitals will be exported.
            Required for cube exports.

        Raises
        ------
        RuntimeError
            If a grid has not been generated or only the molecular geometry
            was parsed.
        ValueError
            If an unsupported ``filetype`` is provided, or if ``mo_index`` is
            missing when exporting cube files.
        """
        if not hasattr(self, '_grid'):
            raise RuntimeError('Grid is not defined. Please create a grid before exporting.')

        if any(dim <= 0 for dim in self._grid_dimensions):
            raise RuntimeError('Grid dimensions are not defined. Create a grid before exporting.')

        if mo_index is not None and (mo_index < 0 or mo_index >= len(self._parser.mos)):
            raise ValueError('Provided molecular orbital index is out of range.')

        if self._only_molecule:
            raise RuntimeError('Orbital exports are unavailable when only the molecule was parsed.')

        destination = Path(path)
        filetype = destination.suffix
        logger.info('Exporting data to %s (format: %s).', destination, filetype or 'unknown')

        if filetype == '.vtk':
            self.export_vtk(destination, mo_index)
        elif filetype == '.cube':
            if mo_index is None:
                raise ValueError('Cube exports require a molecular orbital index.')
            self.export_cube(destination, mo_index)
        else:
            raise ValueError("Unsupported export format. Use '.vtk' or '.cube'.")

    def export_vtk(self, destination: Path, mo_index: int | None = None) -> None:
        """Write orbital data to a VTK structured-grid dataset."""
        # Import lazily so tabulator-only workflows do not require PyVista/VTK at import time.
        import pyvista as pv  # ruff:ignore[import-outside-top-level]

        if not self.has_gtos:
            self.tabulate_gtos()

        mo_data = self.tabulate_mos(mo_index)
        dims = self._grid_dimensions[::-1]
        logger.debug('Writing VTK file %s with grid dimensions %s.', destination, dims)

        struct_grid = pv.StructuredGrid()
        struct_grid.points = self._grid.copy()
        struct_grid.dimensions = dims

        if mo_index is None:
            for mo_ind in range(mo_data.shape[1]):
                struct_grid.point_data[f'mo_{mo_ind}'] = mo_data[:, mo_ind]
        else:
            struct_grid.point_data[f'mo_{mo_index}'] = mo_data

        struct_grid.save(str(destination))

    def export_cube(self, destination: Path, mo_index: int) -> None:
        """Write a single molecular orbital to a Gaussian cube file."""
        cube_values_per_line = 6
        if self._grid_type != GridType.CARTESIAN or self._grid_axes is None:
            raise RuntimeError('Cube exports are only supported for Cartesian grids.')

        mo_values = self.tabulate_mos(mo_index)
        logger.debug('Writing cube file %s for MO %d.', destination, mo_index)

        x, y, z = self._grid_axes
        dx = self._axis_spacing(x, 'x')
        dy = self._axis_spacing(y, 'y')
        dz = self._axis_spacing(z, 'z')

        nx = len(x)
        ny = len(y)
        nz = len(z)

        with destination.open('w', encoding='ascii') as cube_file:
            cube_file.write('Generated by moldenViz Tabulator\n')
            cube_file.write(f'Molecular orbital {mo_index}\n')
            cube_file.write(f'{len(self._parser.atoms):5d} {x[0]:13.6f} {y[0]:13.6f} {z[0]:13.6f}\n')
            cube_file.write(f'{nx:5d} {dx:13.6f} {0.0:13.6f} {0.0:13.6f}\n')
            cube_file.write(f'{ny:5d} {0.0:13.6f} {dy:13.6f} {0.0:13.6f}\n')
            cube_file.write(f'{nz:5d} {0.0:13.6f} {0.0:13.6f} {dz:13.6f}\n')

            for atom in self._parser.atoms:
                cube_file.write(
                    f'{atom.atomic_number:5d} {0.0:13.6f} '
                    f'{atom.position[0]:13.6f} {atom.position[1]:13.6f} {atom.position[2]:13.6f}\n',
                )

            data_3d = mo_values.reshape(self._grid_dimensions, order='C')
            for ix in range(nx):
                for iy in range(ny):
                    for iz in range(nz):
                        cube_file.write(f'{data_3d[ix, iy, iz]:13.5e} ')
                        if iz % cube_values_per_line == (cube_values_per_line - 1):
                            cube_file.write('\n')
                    cube_file.write('\n')

    @staticmethod
    def _tabulate_xlms(theta: NDArray[np.floating], phi: NDArray[np.floating], lmax: int) -> NDArray[np.floating]:
        r"""Tabulate the real spherical harmonics for given theta and phi values.

        We define the real spherical harmonics, Xlms
        (see eq.6, M.A. Blanco et al./Journal of Molecular Structure (Theochem) 419 (1997) 19-27), as:

        Xlms = sqrt(2)*Pl|m|s(\theta)*sin(|m|\phi), m<0
        Xlms = sqrt(2)*Plms(\theta)*cos(m\phi), m>0
        Xlms =         Plms             , m=0

        Note: Above, the Plms are normalized, i.e, \Theta_{lm}(\theta) in eq 1 of the paper.

        Parameters
        ----------
        theta : NDArray[np.floating]
            Array of theta values.
        phi : NDArray[np.floating]
            Array of phi values.
        lmax : int
            Maximum angular momentum quantum number.

        Returns
        -------
        NDArray[np.floating]
            Tabulated real spherical harmonics.

        Raises
        ------
        ValueError
            If input arrays are not 1D or of the same size, or if lmax is negative.
        """
        if theta.ndim != 1 or phi.ndim != 1 or theta.size != phi.size or lmax < 0:
            raise ValueError('Invalid input: theta and phi must be 1D arrays of the same size.')
        if theta.size == 0 or phi.size == 0:
            raise ValueError('Input arrays theta and phi must not be empty.')
        if lmax < 0:
            raise ValueError('lmax must be a non-negative integer.')

        plms = s_plm(lmax, lmax, Tabulator._check_bounds(np.cos(theta)), norm=True)[0] / np.sqrt(2 * np.pi)

        xlms = np.empty_like(plms, dtype=float)
        xlms[:, 0, :] = plms[:, 0, :]

        for m in range(1, lmax + 1):
            factor = -1 if (m % 2) else 1  # Condon-Shortley phase
            xlms[:, -m, :] = factor * np.sqrt(2) * plms[:, m, :] * np.sin(m * phi)
            xlms[:, m, :] = factor * np.sqrt(2) * plms[:, m, :] * np.cos(m * phi)

        return xlms

    @staticmethod
    def _tabulate_real_solid_harmonics(
        centered_grid: NDArray[np.floating],
        lmax: int,
    ) -> NDArray[np.floating]:
        r"""Tabulate normalized real solid harmonics directly in Cartesian coordinates.

        The polynomial is obtained by differentiating the finite power series
        for the Legendre polynomial, following the Rodrigues formulas in
        NIST DLMF sections 14.30 and 18.5:

        .. math::

            r^l X_{lm} =
            N_{lm} Q_m(x, y)
            \sum_{k=0}^{\lfloor(l-m)/2\rfloor}
            \frac{(-1)^k(2l-2k)!}
                 {2^l k!(l-k)!(l-2k-m)!}
            z^{l-m-2k}(x^2+y^2+z^2)^k,

        where :math:`Q_m` is the real or imaginary part of
        :math:`(x+iy)^m`. The normalization ``N_lm`` preserves the spherical
        kernel's normalization and Condon--Shortley phase convention.

        Parameters
        ----------
        centered_grid : NDArray[np.floating]
            Cartesian coordinates relative to an atom, shaped ``(n, 3)``.
        lmax : int
            Maximum angular momentum quantum number.

        Returns
        -------
        NDArray[np.floating]
            Solid harmonics indexed by ``[l, m, point]``. Negative ``m``
            values use NumPy's negative indexing convention.

        Raises
        ------
        ValueError
            If the grid is not a non-empty ``(n, 3)`` array or ``lmax`` is
            negative.
        """
        expected_dimensions = 2
        cartesian_dimensions = 3
        if (
            centered_grid.ndim != expected_dimensions
            or centered_grid.shape[1] != cartesian_dimensions
            or centered_grid.shape[0] == 0
        ):
            raise ValueError('centered_grid must be a non-empty array shaped (n, 3).')
        if lmax < 0:
            raise ValueError('lmax must be a non-negative integer.')

        x, y, z = np.asarray(centered_grid, dtype=float).T
        r_sq = x * x + y * y + z * z
        solid_harmonics = np.zeros((lmax + 1, 2 * lmax + 1, x.size), dtype=float)

        xy_real = np.ones_like(x)
        xy_imag = np.zeros_like(x)
        for m in range(lmax + 1):
            if m:
                xy_real, xy_imag = xy_real * x - xy_imag * y, xy_real * y + xy_imag * x

            for l in range(m, lmax + 1):
                polynomial = np.zeros_like(x)
                for k in range((l - m) // 2 + 1):
                    coefficient = (
                        (-1) ** k
                        * factorial(2 * l - 2 * k)
                        / (2**l * factorial(k) * factorial(l - k) * factorial(l - 2 * k - m))
                    )
                    polynomial += coefficient * z ** (l - m - 2 * k) * r_sq**k

                normalization = np.sqrt(
                    (2 * l + 1) * factorial(l - m) / (4 * np.pi * factorial(l + m)),
                )
                if m == 0:
                    solid_harmonics[l, 0, :] = normalization * polynomial
                else:
                    scale = np.sqrt(2) * normalization
                    solid_harmonics[l, m, :] = scale * xy_real * polynomial
                    solid_harmonics[l, -m, :] = scale * xy_imag * polynomial

        return solid_harmonics

    @staticmethod
    def _check_bounds(x: np.ndarray) -> np.ndarray:
        """Ensure that x is within the open interval (-1, 1) for scipy's associated Legendre polynomial.

        Parameters
        ----------
        x : np.ndarray
            Input array of x values.

        Returns
        -------
        np.ndarray
            Clipped array with values in the open interval (-1, 1).
        """
        min_x = np.nextafter(-1.0, 0.0)
        max_x = np.nextafter(1.0, 0.0)
        return np.clip(x, min_x, max_x)
