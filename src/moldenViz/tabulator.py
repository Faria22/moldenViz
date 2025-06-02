"""Tabulator module for creating grids and tabulating Gaussian-type orbitals (GTOs) from Molden files."""

import logging
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from scipy.special import assoc_legendre_p_all as leg_p_all

from .parser import Parser

logger = logging.getLogger(__name__)

array_like_type = NDArray[np.integer] | list[int] | tuple[int, ...] | range


def _spherical_to_cartersian(
    r: NDArray[np.floating],
    theta: NDArray[np.floating],
    phi: NDArray[np.floating],
) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """Convert spherical coordinates to Cartesian coordinates.

    Args:
    ----
        r (NDArray[np.floating]): Radial distances.
        theta (NDArray[np.floating]): Polar angles (in radians).
        phi (NDArray[np.floating]): Azimuthal angles (in radians).

    Returns:
    -------
        tuple: Arrays of x, y, z Cartesian coordinates.

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

    Args:
    ----
        x (NDArray[np.floating]): X coordinates.
        y (NDArray[np.floating]): Y coordinates.
        z (NDArray[np.floating]): Z coordinates.

    Returns:
    -------
        tuple: Arrays of r (radius), theta (polar angle), phi (azimuthal angle).

    """
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)

    return r, theta, phi


class Tabulator(Parser):
    """Tabulator class for creating grids and tabulating Gaussian-type orbitals (GTOs) from Molden files.

    Inherits from:
        Parser: Parses Molden files to extract molecular and orbital information.
    """

    def __init__(
        self,
        filename: Optional[str] = None,
        molden_lines: Optional[list[str]] = None,
    ) -> None:
        """Initialize the Tabulator with a Molden file or its content.

        Args:
            filename (Optional[str]): Path to the Molden file.
            molden_lines (Optional[list[str]]): Lines of a Molden file.


        Note: If both `filename` and `molden_lines` are provided, or if neither is provided,
              the class will raise a ValueError. Only one of them should be provided. See the
              `Parser` class for more details on how it handles these parameters.

        """
        super().__init__(filename, molden_lines)

        self.grid: NDArray[np.floating] | None = None
        self.gtos_data: NDArray[np.floating] | None = None
        self.mos_data: NDArray[np.floating] | None = None

    def cartesian_grid(
        self,
        x: NDArray[np.floating],
        y: NDArray[np.floating],
        z: NDArray[np.floating],
        tabulate_gtos: bool = True,
    ) -> None:
        r"""Create cartesian grid from x, y, z arrays and tabulate GTOs.

        Args:
        ----
            x (NDArray[np.floating]): Array of x coordinates.
            y (NDArray[np.floating]): Array of y coordinates.
            z (NDArray[np.floating]): Array of z coordinates.
            tabulate_gtos (bool, optional): Whether to tabulate Gaussian-type orbitals (GTOs) after creating the grid.
                Defaults to True.

        """
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        self.grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))

        if tabulate_gtos:
            self.gtos_data = self.tabulate_gtos()

    def spherical_grid(
        self,
        r: NDArray[np.floating],
        theta: NDArray[np.floating],
        phi: NDArray[np.floating],
        tabulate_gtos: bool = True,
    ) -> None:
        r"""Create spherical grid from r, theta, phi arrays and tabulate GTOs.

        Note: Grid points are converted to Cartesian coordinates.

        Args:
        ----
            r (NDArray[np.floating]): Array of radial coordinates.
            theta (NDArray[np.floating]): Array of polar angles (in radians).
            phi (NDArray[np.floating]): Array of azimuthal angles (in radians).
            tabulate_gtos (bool, optional): Whether to tabulate Gaussian-type orbitals (GTOs) after creating the grid.
                Defaults to True.

        """
        rr, tt, pp = np.meshgrid(r, theta, phi, indexing='ij')
        xx, yy, zz = _spherical_to_cartersian(rr, tt, pp)
        self.grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))

        if tabulate_gtos:
            self.gtos_data = self.tabulate_gtos()

    def tabulate_gtos(self) -> NDArray[np.floating]:
        """Tabulate Gaussian-type orbitals (GTOs) on the current grid.

        Returns
        -------
            NDArray[np.floating]: Array containing the tabulated GTOs data.

        Raises
        ------
            ValueError: If the grid is not defined before tabulating GTOs.

        """
        if self.grid is None:
            raise ValueError('Grid is not defined. Please create a grid before tabulating GTOs.')

        # Having a predefined array makes it faster to fill the data
        gto_data = np.empty((self.grid.shape[0], len(self.mos[0].coeffs)))
        ind = 0
        for atom in self.atoms:
            centered_grid = self.grid - atom.position
            max_l = atom.gtos[-1].l

            r, theta, phi = _cartesian_to_spherical(*centered_grid.T)
            xlms = self._tabulate_xlms(theta, phi, max_l)

            for gto in atom.gtos:
                l = gto.l
                m_inds = np.arange(-l, l + 1)
                gto_inds = ind + l + m_inds

                radial = gto.norm * r**l * sum(prim.norm * prim.coeff * np.exp(-prim.exp * r**2) for prim in gto.prims)

                gto_data[:, gto_inds] = radial[:, None] * xlms[l, m_inds, ...].T

                ind += 2 * l + 1

        logger.debug('GTO data shape: %s', gto_data.shape)
        return gto_data

    def tabulate_mos(self, mo_inds: Optional[int | array_like_type] = None) -> NDArray[np.floating]:
        """Tabulate molecular orbitals (MOs) on the current grid.

        Args:
        ----
            mo_inds (Optional[int | ArrayLike]): Indices of the MOs to tabulate. If None, all MOs are tabulated.

        Returns:
        -------
            NDArray[np.floating]: Array containing the tabulated MOs data.

            If an integer is provided, it will tabulate only that MO.
            If an array-like is provided, it will tabulate the MOs at those indices.

        Raises:
        ------
            ValueError: If the grid is not defined before tabulating MOs.
            ValueError: If GTOs are not tabulated before tabulating MOs.
            ValueError: If provided mo_inds invalid.

        """
        if self.grid is None:
            raise ValueError('Grid is not defined. Please create a grid before tabulating MOs.')

        if self.gtos_data is None:
            raise ValueError('GTOs are not tabulated. Please tabulate GTOs before tabulating MOs.')

        if mo_inds is None:
            mo_inds = list(range(len(self.mos)))

        if isinstance(mo_inds, range):
            mo_inds = list(mo_inds)

        if not isinstance(mo_inds, int) and not mo_inds:
            raise ValueError('Provided mo_inds is empty. Please provide valid indices.')

        if isinstance(mo_inds, int):
            if mo_inds < 0 or mo_inds >= len(self.mos):
                raise ValueError('Provided mo_index is invalid. Please provide valid index.')
        elif any(mo_ind < 0 or mo_ind >= len(self.mos) for mo_ind in mo_inds):
            raise ValueError('Provided mo_inds contains invalid indices. Please provide valid indices.')

        if isinstance(mo_inds, int):
            mo_data = np.sum(self.gtos_data * self.mos[mo_inds].coeffs[None, :], axis=1)
        else:
            mo_coeffs = np.stack([self.mos[i].coeffs for i in mo_inds])

            mo_data = np.sum(self.gtos_data[:, None, :] * mo_coeffs[None, ...], axis=2)
            logger.debug('MO data shape: %s', mo_data.shape)

        return mo_data

    @staticmethod
    def _tabulate_xlms(theta: NDArray[np.floating], phi: NDArray[np.floating], lmax: int) -> NDArray[np.floating]:
        r"""Tabulate the real spherical harmonics for given theta and phi values.

        We define the real spherical harmonics, Xlms
        (see eq.6, M.A. Blanco et al./Journal of Molecular Structure (Theochem) 419 (1997) 19-27), as:

        Xlms = sqrt(2)*Pl|m|(\theta)*sin(|m|\phi), m<0
        Xlms = sqrt(2)*Plms(\theta)*cos(|m|\phi), m>0
        Xlms =         Plms             , m=0

        Note: Above, the Plms are normalized, i.e, \Theta_{lm}(\theta) in eq 1 of the paper.

        Args:
        ----
            theta (NDArray[np.floating]): Array of theta values.
            phi (NDArray[np.floating]): Array of phi values.
            lmax (int): Maximum angular momentum quantum number.

        Returns:
        -------
            NDArray[np.floating]: Tabulated real spherical harmonics.

        """
        # leg_p_all first dimension has always size = 1
        xlms = leg_p_all(lmax, lmax, np.cos(theta), norm=True)[0, ...] / np.sqrt(np.pi)
        xlms[:, 0, :] /= np.sqrt(2)

        for m in range(1, lmax + 1):
            xlms[:, -m, :] = (-1) ** m * xlms[:, m, :] * np.sin(m * phi)
            xlms[:, m, :] *= (-1) ** m * np.cos(m * phi)

        return xlms
