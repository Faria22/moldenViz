"""Shared inputs and helpers for the benchmark suite."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np

from moldenViz import examples
from moldenViz.tabulator import Tabulator

if TYPE_CHECKING:
    from numpy.typing import NDArray

EXAMPLE_NAMES = (
    'acrolein',
    'benzene',
    'co',
    'co2',
    'furan',
    'h2o',
    'o2',
    'prismane',
    'pyridine',
)
GRID_EDGES = (10, 25, 50, 100)
POINT_CHUNK_SIZES = (8_192, 32_768, 65_536, None)
MO_SELECTIONS = ('single', 'several', 'all')
REPRESENTATIVE_EXAMPLES = ('h2o', 'furan', 'benzene')
WORKER_COUNTS = (1, 4)
PYSCF_CO_SPHERICAL = Path(__file__).parents[1] / 'tests/fixtures/pyscf/co-cc-pvqz-spherical.molden'

MOSelection = Literal['single', 'several', 'all']


class GenericSolidHarmonicTabulator(Tabulator):
    """Tabulator using the finite-polynomial solid-harmonic oracle."""

    @staticmethod
    def _tabulate_real_solid_harmonics(
        centered_grid: NDArray[np.floating],
        lmax: int,
    ) -> NDArray[np.floating]:
        return Tabulator._tabulate_real_solid_harmonics_generic(  # ruff:ignore[private-member-access]
            centered_grid,
            lmax,
        )


def example_source(name: str) -> list[str]:
    """Return the bundled Molden source named by ``name``.

    Returns
    -------
    list[str]
        Molden input lines for the requested example.
    """
    source = getattr(examples, name)
    if not isinstance(source, list):
        raise TypeError(f'Example {name!r} did not provide Molden input lines.')
    return source


def grid_axis(edge_size: int) -> np.ndarray:
    """Return one axis for an ``edge_size ** 3`` Cartesian grid.

    Returns
    -------
    np.ndarray
        Evenly spaced grid coordinates.
    """
    return np.linspace(-3.0, 3.0, edge_size)


def mo_indices(tabulator: Tabulator, selection: MOSelection) -> int | list[int] | None:
    """Map a benchmark selection label to public ``tabulate_mos`` input.

    Returns
    -------
    int | list[int] | None
        Indices representing one, several, or all molecular orbitals.
    """
    if selection == 'single':
        return 0
    if selection == 'several':
        return list(range(min(5, len(tabulator.molecular_orbitals))))
    return None
