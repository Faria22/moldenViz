"""Tabulator module for creating grids and tabulating Gaussian-type orbitals (GTOs) from Molden files."""

import logging
from enum import Enum
from functools import lru_cache
from math import factorial
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pyvista as pv
from numpy.typing import NDArray

from .parser import Parser

logger = logging.getLogger(__name__)

array_like_type = NDArray[np.integer] | list[int] | tuple[int, ...] | range


def _grid_creation_with_only_molecule_error() -> RuntimeError:
    """Return a consistent error for grid creation when only the molecule is parsed.

    Returns
    -------
    RuntimeError
        The error indicating that grid creation is not allowed when `only_molecule` is set to `True`.

    """
    return RuntimeError('Grid creation is not allowed when `only_molecule` is set to `True`.')


def _spherical_to_cartesian(
    r: NDArray[np.floating],
    theta: NDArray[np.floating],
    phi: NDArray[np.floating],
) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """Convert spherical coordinates to Cartesian coordinates.

    Parameters
    ----------
    r : NDArray[np.floating]
        1D array of radial coordinates.
    theta : NDArray[np.floating]
        1D array of polar angles (radians).
    phi : NDArray[np.floating]
        1D array of azimuthal angles (radians).

    Returns
    -------
    tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]
        Tuple containing the Cartesian coordinates ``(x, y, z)``.
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    return x, y, z


def _cartesian_to_spherical(
    x: NDArray[np.floating],
    y: NDArray[np.floating],
    z: NDArray[np.floating],
) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """Convert Cartesian coordinates to spherical coordinates.

    Parameters
    ----------
    x : NDArray[np.floating]
        1D array of x coordinates.
    y : NDArray[np.floating]
        1D array of y coordinates.
    z : NDArray[np.floating]
        1D array of z coordinates.

    Returns
    -------
    tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]
        Tuple containing spherical coordinates ``(r, theta, phi)``.
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)

    return r, theta, phi


@lru_cache(maxsize=None)
def _cached_factorial(n: int) -> int:
    """Compute factorial with caching for non-negative integers.

    Returns
    -------
    int
        The factorial of n.
    """
    return factorial(n)


@lru_cache(maxsize=None)
def _binomial(r: float, k: int) -> float:
    """Calculate the generalized binomial coefficient (r over k).

    This function supports real or complex 'r' and non-negative integer 'k'.
    The formula is: C(r, k) = r * (r-1) * ... * (r-k+1) / k!

    Parameters
    ----------
    r : float
        A real or complex number.
    k : int
        A non-negative integer.

    Returns
    -------
    float
        The value of the generalized binomial coefficient as a float or complex number.

    Raises
    ------
    ValueError
        If k is a negative integer.
    """
    if not isinstance(k, int) or k < 0:
        raise ValueError('k must be a non-negative integer.')

    if k == 0:
        return 1
    if k == 1:
        return r

    numerator = 1
    for i in range(k):
        numerator *= r - i

    return numerator / _cached_factorial(k)


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
        The tabulated Gaussian-type orbitals (GTOs) on the grid.
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
        self._gtos: NDArray[np.floating]

        # Used for when exporting to cube format
        self.original_axes: tuple[NDArray[np.floating], ...] | None = None

    @property
    def grid(self) -> NDArray[np.floating]:
        """Return the 2D array of Cartesian grid points used for tabulation."""
        return self._grid

    @grid.setter
    def grid(self, new_grid: Any) -> None:
        """Set the tabulation grid after validating its structure.

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

        del self.grid
        self._grid = new_grid

    @grid.deleter
    def grid(self) -> None:
        """Delete the cached grid and mark its type as unknown."""
        del self._grid
        self._grid_type = GridType.UNKNOWN
        self._grid_dimensions = (0, 0, 0)
        self.original_axes = None

    @property
    def gtos(self) -> NDArray[np.floating]:
        """Get the tabulated Gaussian-type orbitals (GTOs) on the grid."""
        return self._gtos

    @gtos.deleter
    def gtos(self) -> None:
        del self._gtos

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

        self.original_axes = (x, y, z)

        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        if grid_type == GridType.SPHERICAL:
            xx, yy, zz = _spherical_to_cartesian(xx, yy, zz)

        self._grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        self._grid_type = grid_type
        self._grid_dimensions = (len(x), len(y), len(z))

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

        # Having a predefined array makes it faster to fill the data
        gto_data = np.empty((self._grid.shape[0], self._parser.mo_coeffs.shape[1]))
        self._atom_gto_slices = []
        ind = 0
        for atom in self._parser.atoms:
            atom_start = ind
            centered_grid = self._grid - atom.position
            max_l = atom.shells[-1].l

            r, theta, phi = _cartesian_to_spherical(*centered_grid.T)  # pyright: ignore[reportArgumentType]
            xlms = self._tabulate_xlms(theta, phi, max_l)

            for shell in atom.shells:
                l = shell.l
                m_inds = np.arange(-l, l + 1)
                gto_inds = ind + l + m_inds

                radial = shell.norm * r**l * sum(gto.norm * gto.coeff * np.exp(-gto.exp * r**2) for gto in shell.gtos)

                gto_data[:, gto_inds] = radial[:, None] * xlms[l, m_inds, ...].T

                ind += 2 * l + 1

            self._atom_gto_slices.append(slice(atom_start, ind))

        logger.debug('GTO data shape: %s', gto_data.shape)

        self._gtos = gto_data
        return gto_data

    def tabulate_mos(self, mo_inds: Optional[int | array_like_type] = None) -> NDArray[np.floating]:
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
        if not hasattr(self, 'gtos'):
            raise RuntimeError('GTOs are not tabulated. Please tabulate GTOs before tabulating MOs.')

        if mo_inds is None:
            mo_inds = list(range(len(self._parser.mos)))

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
            mo_data = np.sum(self.gtos * self._parser.mo_coeffs[mo_inds][None, :], axis=1)
        else:
            # Use direct slicing of mo_coeffs array
            mo_coeffs = self._parser.mo_coeffs[mo_inds]

            mo_data = np.sum(self.gtos[:, None, :] * mo_coeffs[None, ...], axis=2)
            logger.debug('MO data shape: %s', mo_data.shape)

        return mo_data

    def export(self, path: str | Path, *, mo_index: Optional[int] = None) -> None:
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

        if filetype == '.vtk':
            self._export_vtk(destination, mo_index)
        elif filetype == '.cube':
            if mo_index is None:
                raise ValueError('Cube exports require a molecular orbital index.')
            self._export_cube(destination, mo_index)
        else:
            raise ValueError("Unsupported export format. Use '.vtk' or '.cube'.")

    def _export_vtk(self, destination: Path, mo_index: Optional[int] = None) -> None:
        """Write orbital data to a VTK multiblock dataset."""
        if not hasattr(self, 'gtos'):
            self.tabulate_gtos()

        mo_data = self.tabulate_mos(mo_index)
        dims = self._grid_dimensions[::-1]

        struct_grid = pv.StructuredGrid()
        struct_grid.points = self._grid.copy()
        struct_grid.dimensions = dims

        if mo_index is None:
            for mo_ind in range(mo_data.shape[1]):
                struct_grid.point_data[f'mo_{mo_ind}'] = mo_data[:, mo_ind]
        else:
            struct_grid.point_data[f'mo_{mo_index}'] = mo_data

        struct_grid.save(destination)

    def _export_cube(self, destination: Path, mo_index: int) -> None:
        """Write a single molecular orbital to a Gaussian cube file."""
        if self._grid_type != GridType.CARTESIAN or self.original_axes is None:
            raise RuntimeError('Cube exports are only supported for Cartesian grids.')

        mo_values = self.tabulate_mos(mo_index)

        x, y, z = self.original_axes
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
                        if iz % 6 == 5:  # noqa: PLR2004
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

        plms = Tabulator._tabulate_plms(np.cos(theta), lmax)

        xlms = np.empty((lmax + 1, 2 * lmax + 1, theta.size), dtype=float)
        xlms[:, 0, :] = plms[:, 0, :]

        for m in range(1, lmax + 1):
            xlms[:, -m, :] = np.sqrt(2) * plms[:, m, :] * np.sin(m * phi)
            xlms[:, m, :] = np.sqrt(2) * plms[:, m, :] * np.cos(m * phi)

        return xlms

    @staticmethod
    def _tabulate_plms(x: NDArray[np.floating], lmax: int) -> NDArray[np.floating]:
        """Tabulate normalized associated Legendre polynomials (without Condon-Shortley phase).

        Returns a 3D array with axes `(l, m, index)` where:
            The first axis enumerates ``l`` in ``[0, lmax]``.
            The second axis enumerates ``m`` in ``[0, lmax]``.
            The third axis spans the samples of ``x``.

        Using closed form outlined here:
        https://en.m.wikipedia.org/wiki/Associated_Legendre_polynomials#Closed_Form

        Parameters
        ----------
        x : NDArray[np.floating]
            Array of x values.
        lmax : int
            Maximum angular momentum quantum number.

        Returns
        -------
        NDArray[np.floating]
            Tabulated associated Legendre polynomials.
        """
        plms = np.empty((lmax + 1, lmax + 1, x.size), dtype=float)

        for l in range(lmax + 1):
            for m in range(l + 1):
                plms[l, m, :] = (
                    np.sqrt((2 * l + 1) * _cached_factorial(l - m) / _cached_factorial(l + m) / 4 / np.pi)
                    * 2**l
                    * (1 - x**2) ** (m / 2)
                    * sum(
                        _cached_factorial(k)
                        * x ** (k - m)
                        * _binomial(l, k)
                        * _binomial((l + k - 1) / 2, l)
                        / _cached_factorial(k - m)
                        for k in range(m, l + 1)
                    )
                )

        return plms
