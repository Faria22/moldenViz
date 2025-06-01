"""Tabulator module for creating grids and tabulating Gaussian-type orbitals (GTOs) from Molden files."""

import logging
from parser import Parser
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from scipy.special import assoc_legendre_p_all as leg_p_all

logger = logging.getLogger(__name__)


def _spherical_to_cartersian(
    r: NDArray[np.floating],
    theta: NDArray[np.floating],
    phi: NDArray[np.floating],
) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """Convert spherical coordinates to Cartesian coordinates.

    Args:
        r (NDArray[np.floating]): Radial distances.
        theta (NDArray[np.floating]): Polar angles (in radians).
        phi (NDArray[np.floating]): Azimuthal angles (in radians).

    Returns:
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
        x (NDArray[np.floating]): X coordinates.
        y (NDArray[np.floating]): Y coordinates.
        z (NDArray[np.floating]): Z coordinates.

    Returns:
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

        self.grid: Optional[NDArray[np.floating]] = None

    def carterian_grid(
        self,
        x: NDArray[np.floating],
        y: NDArray[np.floating],
        z: NDArray[np.floating],
        tabulate_gtos: bool = True,
    ) -> None:
        r"""Create cartersian grid from x, y, z arrays.

        Args:
            x (NDArray[np.floating]): Array of x coordinates.
            y (NDArray[np.floating]): Array of y coordinates.
            z (NDArray[np.floating]): Array of z coordinates.
            tabulate_gtos (bool, optional): Whether to tabulate Gaussian-type orbitals (GTOs) after creating the grid.
                Defaults to True.

        """
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        self.grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))

        if tabulate_gtos:
            self.tabulate_gtos()

    def spherical_grid(
        self,
        r: NDArray[np.floating],
        theta: NDArray[np.floating],
        phi: NDArray[np.floating],
        tabulate_gtos: bool = True,
    ) -> None:
        r"""Create spherical grid from r, theta, phi arrays. Grid points are converted to Cartesian coordinates.

        Args:
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
            self.tabulate_gtos()

    def tabulate_gtos(self) -> None:
        """Tabulate Gaussian-type orbitals (GTOs) on the current grid.

        Raises:
            ValueError: If the grid is not defined before tabulating GTOs.

        """
        if self.grid is None:
            raise ValueError('Grid is not defined. Please create a grid before tabulating GTOs.')

        # Having a predefined array makes it faster to fill the data
        gto_data = np.empty((self.grid.shape[0], len(self.mo_coeffs[0].coeffs)))
        ind = 0
        for atom in self.atoms:
            centered_grid = self.grid - atom.position
            max_l = max(gto.l for gto in atom.gtos)

            r, theta, phi = _cartesian_to_spherical(*centered_grid.T)
            xlms = self._tabulate_xlms(theta, phi, max_l)
            for gto in atom.gtos:
                l = gto.l
                inds = ind + np.arange(2 * l + 1)

                data = gto.norm * r**l * sum(prim.coeff * np.exp(-prim.exp * r**2) for prim in gto.prims)

                gto_data[:, inds] = data[:, None] * xlms[l, range(-l, l + 1), ...].T

                ind += 2 * l + 1

        logger.debug('GTO data shape: %s', gto_data.shape)

    @staticmethod
    def _tabulate_xlms(theta: NDArray[np.floating], phi: NDArray[np.floating], lmax: int) -> NDArray[np.floating]:
        r"""Tabulate the real spherical harmonics for given theta and phi values.

        We define the real spherical harmonics, Xlms
        (see eq.6, M.A. Blanco et al./Journal of Molecular Structure (Theochem) 419 (1997) 19-27), as:

        Xlms = sqrt(2)*Plms*sin(|m|\phi), m<0
        Xlms = sqrt(2)*Plms*cos(|m|\phi), m>0
        Xlms =         Plms             , m=0

        Note: Above, the Plms are normalized

        Args:
            theta (NDArray[np.floating]): Array of theta values.
            phi (NDArray[np.floating]): Array of phi values.
            lmax (int): Maximum angular momentum quantum number.

        Returns:
            NDArray[np.floating]: Tabulated real spherical harmonics.

        """
        # leg_p_all first dimension has always size = 1
        xlms = leg_p_all(lmax, lmax, np.cos(theta), norm=True)[0, ...] / np.sqrt(2 * np.pi)

        for m in range(-lmax, lmax + 1):
            if m < 0:
                xlms[:, m, :] *= np.sqrt(2) * np.sin(-m * phi)
            elif m > 0:
                xlms[:, m, :] *= np.sqrt(2) * np.cos(m * phi)

        return xlms
