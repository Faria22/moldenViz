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
PYSCF_SPHERICAL_EXAMPLES = ('co', 'h2o', 'benzene')
PYSCF_FIXTURE_DIR = Path(__file__).parents[1] / 'tests/fixtures/pyscf'

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


def example_content(name: str) -> str:
    """Return the bundled Molden content named by ``name``.

    Returns
    -------
    str
        Complete Molden input for the requested example.
    """
    content = getattr(examples, name)
    if not isinstance(content, str):
        raise TypeError(f'Example {name!r} did not provide Molden content.')
    return content


def grid_axis(edge_size: int) -> np.ndarray:
    """Return one axis for an ``edge_size ** 3`` Cartesian grid.

    Returns
    -------
    np.ndarray
        Evenly spaced grid coordinates.
    """
    return np.linspace(-3.0, 3.0, edge_size)


def pyscf_spherical_path(name: str) -> Path:
    """Return the committed spherical cc-pVQZ fixture for ``name``.

    Returns
    -------
    Path
        Path to the requested PySCF fixture.
    """
    return PYSCF_FIXTURE_DIR / f'{name}-cc-pvqz-spherical.molden'


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
