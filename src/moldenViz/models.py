"""Public data models used by moldenViz."""

from __future__ import annotations

from dataclasses import dataclass
from math import gamma
from typing import TYPE_CHECKING

import numpy as np
from pydantic import BaseModel, Field

if TYPE_CHECKING:
    from numpy.typing import NDArray

__all__ = ['Atom', 'AtomType', 'GaussianPrimitive', 'MolecularOrbital', 'Shell']


@dataclass
class Atom:
    """A parsed atom and its basis-function shells.

    Parameters
    ----------
    label : str
        Atomic label from the Molden file.
    atomic_number : int
        Atomic number of the element.
    position : NDArray[np.floating]
        Three-dimensional position in Bohr.
    shells : list[Shell]
        Electron shells associated with the atom.
    """

    label: str
    atomic_number: int
    position: NDArray[np.floating]
    shells: list[Shell]


@dataclass
class MolecularOrbital:
    """Metadata for a parsed molecular orbital."""

    sym: str
    energy: float
    spin: str
    occ: int


class GaussianPrimitive:
    """A Gaussian primitive with an exponent and contraction coefficient."""

    def __init__(self, exp: float, coeff: float) -> None:
        self.exp = exp
        self.coeff = coeff
        self._norm = 0.0

    def _normalize(self, l: int) -> None:
        """Calculate and cache the primitive normalization for angular momentum ``l``."""
        self._norm = np.sqrt(2 * (2 * self.exp) ** (l + 1.5) / gamma(l + 1.5))


class Shell:
    """An electron shell containing Gaussian primitives."""

    def __init__(self, l: int, gtos: list[GaussianPrimitive]) -> None:
        self.l = l
        self.gtos = gtos

        self._norm = 0.0
        self._gto_norms = np.empty(len(gtos), dtype=float)
        self._gto_exps = np.array([gto.exp for gto in gtos], dtype=float)
        self._gto_coeffs = np.array([gto.coeff for gto in gtos], dtype=float)
        self._prefactor = np.empty(len(gtos), dtype=float)

    def _normalize(self) -> None:
        """Calculate and cache shell and primitive normalization values."""
        for idx, gto in enumerate(self.gtos):
            gto._normalize(self.l)  # ruff:ignore[private-member-access]
            self._gto_norms[idx] = gto._norm  # ruff:ignore[private-member-access]

        overlap = 0.0
        for i_gto in self.gtos:
            for j_gto in self.gtos:
                overlap += (
                    i_gto.coeff
                    * j_gto.coeff
                    * (2 * np.sqrt(i_gto.exp * j_gto.exp) / (i_gto.exp + j_gto.exp)) ** (self.l + 1.5)
                )

        self._norm = 1 / np.sqrt(overlap)
        self._prefactor = self._norm * self._gto_norms * self._gto_coeffs


class AtomType(BaseModel):
    """Validated visualization properties for an atomic element."""

    name: str = Field(..., min_length=1, max_length=3, description='Atom symbol')
    color: str = Field(..., pattern=r'^[0-9A-Fa-f]{6}$', description='Hex color code without #')
    radius: float = Field(..., gt=0, description='Atom radius (must be positive)')
    max_num_bonds: int = Field(..., ge=0, description='Maximum number of bonds')
