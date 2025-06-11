"""Tabulator module for creating grids and tabulating Gaussian-type orbitals (GTOs) from Molden files."""

import logging
from enum import Enum
from math import factorial
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from .parser import Parser

logger = logging.getLogger(__name__)

array_like_type = NDArray[np.integer] | list[int] | tuple[int, ...] | range


def _spherical_to_cartesian(
    r: NDArray[np.floating],
    theta: NDArray[np.floating],
    phi: NDArray[np.floating],
) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """Convert spherical coordinates to Cartesian coordinates.

    Args
    ----
        r (NDArray[np.floating]): Radial distances.
        theta (NDArray[np.floating]): Polar angles (in radians).
        phi (NDArray[np.floating]): Azimuthal angles (in radians).

    Returns
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

    Args
    ----
        x (NDArray[np.floating]): X coordinates.
        y (NDArray[np.floating]): Y coordinates.
        z (NDArray[np.floating]): Z coordinates.

    Returns
    -------
        tuple: Arrays of r (radius), theta (polar angle), phi (azimuthal angle).

    """
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)

    return r, theta, phi


def _binomial(r: float, k: int) -> float:
    """
    Calculate the generalized binomial coefficient (r over k).

    This function supports real or complex 'r' and non-negative integer 'k'.
    The formula is: C(r, k) = r * (r-1) * ... * (r-k+1) / k!

    Args:
        r: A real or complex number.
        k: A non-negative integer.

    Returns
    -------
        The value of the generalized binomial coefficient as a float or complex number.

    Raises
    ------
        ValueError: If k is a negative integer.
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

    denominator = factorial(k)

    return numerator / denominator


class GridType(Enum):
    """Grid types allowed."""

    SPHERICAL = 'spherical'
    CARTESIAN = 'cartesian'
    UNKNOWN = 'unknown'


class Tabulator(Parser):
    """Extends Parses, create grids and tabulates Gaussian-type orbitals (GTOs) from Molden files.

    Args
    ----
        source: str | list[str]
            The path to the molden file, or the lines from the file.

        only_molecule: bool, optional
            Only parse the atoms and skip molecular orbitals.
            Default is `False`.
    """

    def __init__(
        self,
        source: str | list[str],
        only_molecule: bool = False,
    ) -> None:
        """Initialize the Tabulator with a Molden file or its content."""
        super().__init__(source, only_molecule)

        self.only_molecule = only_molecule
        self.grid_type = GridType.UNKNOWN

        self.grid: NDArray[np.floating]
        self.grid_dimensions: tuple[int, int, int]
        self.gtos_data: NDArray[np.floating]
        self.mos_data: NDArray[np.floating]

    def cartesian_grid(
        self,
        x: NDArray[np.floating],
        y: NDArray[np.floating],
        z: NDArray[np.floating],
        tabulate_gtos: bool = True,
    ) -> None:
        r"""Create cartesian grid from x, y, z arrays and tabulate GTOs.

        Args
        ----
            x: NDArray[np.floating]
                Array of x coordinates.
            y: NDArray[np.floating]
                Array of y coordinates.
            z: NDArray[np.floating]
                Array of z coordinates.
            tabulate_gtos: bool, optional
                Whether to tabulate Gaussian-type orbitals (GTOs) after creating the grid.
                Defaults to True.

        """
        if self.only_molecule:
            raise RuntimeError('Grid creation is not allowed when `only_molecule` is set to `True`.')

        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        self.grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        self.grid_type = GridType.CARTESIAN
        self.grid_dimensions = (len(x), len(y), len(z))

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

        Args
        ----
            r: NDArray[np.floating]
                Array of radial coordinates.
            theta: NDArray[np.floating]
                Array of polar angles (in radians).
            phi: NDArray[np.floating]
                Array of azimuthal angles (in radians).
            tabulate_gtos: bool, optional
                Whether to tabulate Gaussian-type orbitals (GTOs) after creating the grid.
                Defaults to True.

        Note
        ----
            Grid points are converted to Cartesian coordinates.

        """
        if self.only_molecule:
            raise RuntimeError('Grid creation is not allowed when `only_molecule` is set to `True`.')

        rr, tt, pp = np.meshgrid(r, theta, phi, indexing='ij')
        xx, yy, zz = _spherical_to_cartesian(rr, tt, pp)
        self.grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        self.grid_type = GridType.SPHERICAL
        self.grid_dimensions = (len(r), len(theta), len(phi))

        if tabulate_gtos:
            self.gtos_data = self.tabulate_gtos()

    def tabulate_gtos(self) -> NDArray[np.floating]:
        """Tabulate Gaussian-type orbitals (GTOs) on the current grid.

        Returns
        -------
            NDArray[np.floating]:
                Array containing the tabulated GTOs data.

        Raises
        ------
            ValueError:
                If the grid is not defined before tabulating GTOs.

        """
        if self.only_molecule:
            raise RuntimeError('Grid creation is not allowed when `only_molecule` is set to `True`.')

        if not hasattr(self, 'grid'):
            raise ValueError('Grid is not defined. Please create a grid before tabulating GTOs.')

        # Having a predefined array makes it faster to fill the data
        gto_data = np.empty((self.grid.shape[0], len(self.mos[0].coeffs)))
        ind = 0
        for atom in self.atoms:
            centered_grid = self.grid - atom.position
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

        logger.debug('GTO data shape: %s', gto_data.shape)
        return gto_data

    def tabulate_mos(self, mo_inds: Optional[int | array_like_type] = None) -> NDArray[np.floating]:
        """Tabulate molecular orbitals (MOs) on the current grid.

        Args
        ----
            mo_inds : int, ArrayLike, None
                Indices of the MOs to tabulate. If None, all MOs are tabulated.

        Returns
        -------
            NDArray[np.floating]:
                Array containing the tabulated MOs data.

                If an integer is provided, it will tabulate only that MO.
                If an array-like is provided, it will tabulate the MOs at those indices.

        Raises
        ------
            ValueError:
                If the grid is not defined before tabulating MOs.
            ValueError:
                If GTOs are not tabulated before tabulating MOs.
            ValueError:
                If provided mo_inds invalid.

        """
        if not hasattr(self, 'grid'):
            raise ValueError('Grid is not defined. Please create a grid before tabulating MOs.')

        if not hasattr(self, 'gtos_data'):
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

        Xlms = sqrt(2)*Pl|m|s(\theta)*sin(|m|\phi), m<0
        Xlms = sqrt(2)*Plms(\theta)*cos(m\phi), m>0
        Xlms =         Plms             , m=0

        Note: Above, the Plms are normalized, i.e, \Theta_{lm}(\theta) in eq 1 of the paper.

        Args
        ----
            theta: NDArray[np.floating]
                Array of theta values.
            phi: NDArray[np.floating]
                Array of phi values.
            lmax: int
                Maximum angular momentum quantum number.

        Returns
        -------
            NDArray[np.floating]:
                Tabulated real spherical harmonics.

        """
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

        Returns an array of shape (`lmax+1`, `lmax+1`, `x.size`), where:
            The first index is 0 <= l <= lmax
            The second index is 0 <= m <= lmax
            The third index goes over the x points

        Using closed form outlined here:
        https://en.m.wikipedia.org/wiki/Associated_Legendre_polynomials#Closed_Form

        Args
        ----
            x: NDArray[np.floating]
                Array of x values.
            lmax: int
                Maximum angular momentum quantum number.

        Returns
        -------
            NDArray[np.floating]:
                Tabulated associated Legendre polynomials.


        """
        plms = np.zeros((lmax + 1, lmax + 1, x.size), dtype=float)

        for l in range(lmax + 1):
            for m in range(l + 1):
                plms[l, m, :] = (
                    np.sqrt((2 * l + 1) * factorial(l - m) / factorial(l + m) / 4 / np.pi)
                    * 2**l
                    * (1 - x**2) ** (m / 2)
                    * sum(
                        factorial(k) * x ** (k - m) * _binomial(l, k) * _binomial((l + k - 1) / 2, l) / factorial(k - m)
                        for k in range(m, l + 1)
                    )
                )

        return plms
