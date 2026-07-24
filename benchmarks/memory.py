"""Peak-memory benchmarks for representative tabulation workloads."""

import numpy as np

from moldenViz.tabulator import Tabulator

from ._shared import (
    POINT_CHUNK_SIZES,
    PYSCF_CO_SPHERICAL,
    REPRESENTATIVE_EXAMPLES,
    WORKER_COUNTS,
    GenericSolidHarmonicTabulator,
    example_source,
    grid_axis,
)


class PeakMemorySolidHarmonics:
    """Compare generated and generic solid-harmonic temporary memory."""

    params = ((125_000, 1_000_000), (2, 4), ('generated', 'generic'))
    param_names = ['num_points', 'lmax', 'implementation']
    timeout = 180

    def setup(self, num_points: int, lmax: int, implementation: str) -> None:
        """Create deterministic Cartesian points."""
        rng = np.random.default_rng(seed=116)
        self.points = rng.uniform(-6.0, 6.0, size=(num_points, 3))

    def peakmem_solid_harmonics(self, num_points: int, lmax: int, implementation: str) -> None:
        """Measure peak RSS for the selected harmonic implementation."""
        if implementation == 'generated':
            Tabulator._tabulate_real_solid_harmonics(  # ruff:ignore[private-member-access]
                self.points,
                lmax,
            )
        else:
            Tabulator._tabulate_real_solid_harmonics_generic(  # ruff:ignore[private-member-access]
                self.points,
                lmax,
            )


class PeakMemoryHighAngularMomentumGTOTabulation:
    """Compare full cc-pVQZ GTO-tabulation peak memory."""

    params = ((25, 50), ('generated', 'generic'))
    param_names = ['edge_size', 'implementation']
    timeout = 180

    def setup(self, edge_size: int, implementation: str) -> None:
        """Create a sequential cc-pVQZ tabulator with an uncomputed grid."""
        tabulator_class = GenericSolidHarmonicTabulator if implementation == 'generic' else Tabulator
        self.tabulator = tabulator_class(str(PYSCF_CO_SPHERICAL), max_workers=1)
        axis = grid_axis(edge_size)
        self.tabulator.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    def peakmem_tabulate_gtos(self, edge_size: int, implementation: str) -> None:
        """Measure peak RSS while tabulating every cc-pVQZ basis function."""
        self.tabulator.tabulate_gtos()


class PeakMemoryGTOTabulation:
    """Measure process peak RSS while producing large GTO arrays."""

    params = (REPRESENTATIVE_EXAMPLES, (50, 75))
    param_names = ['molecule', 'edge_size']
    timeout = 180

    def setup(self, molecule: str, edge_size: int) -> None:
        """Create a tabulator with an uncomputed Cartesian grid."""
        self.tabulator = Tabulator(example_source(molecule))
        axis = grid_axis(edge_size)
        self.tabulator.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    def peakmem_tabulate_gtos(self, molecule: str, edge_size: int) -> None:
        """Tabulate all GTOs while ASV samples peak resident memory."""
        self.tabulator.tabulate_gtos()


class PeakMemoryGTOChunkSizes:
    """Measure peak RSS across bounded point chunk sizes."""

    params = (REPRESENTATIVE_EXAMPLES, (75,), POINT_CHUNK_SIZES)
    param_names = ['molecule', 'edge_size', 'point_chunk_size']
    timeout = 180

    def setup(self, molecule: str, edge_size: int, point_chunk_size: int | None) -> None:
        """Create a tabulator with an uncomputed Cartesian grid."""
        self.tabulator = Tabulator(example_source(molecule))
        axis = grid_axis(edge_size)
        self.tabulator.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    def peakmem_tabulate_gtos(self, molecule: str, edge_size: int, point_chunk_size: int | None) -> None:
        """Tabulate GTOs with the selected point-chunk bound."""
        self.tabulator.tabulate_gtos(point_chunk_size=point_chunk_size)


class PeakMemoryGTOWorkerScaling:
    """Compare peak RSS for sequential and bounded-parallel GTO work."""

    params = (('benzene',), (75,), WORKER_COUNTS)
    param_names = ['molecule', 'edge_size', 'max_workers']
    timeout = 180

    def setup(self, molecule: str, edge_size: int, max_workers: int) -> None:
        """Create a representative large grid with an explicit worker limit."""
        self.tabulator = Tabulator(example_source(molecule), max_workers=max_workers)
        axis = grid_axis(edge_size)
        self.tabulator.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    def peakmem_tabulate_gtos(self, molecule: str, edge_size: int, max_workers: int) -> None:
        """Tabulate GTOs while ASV samples concurrency-dependent peak RSS."""
        self.tabulator.tabulate_gtos()


class PeakMemoryAllMOContraction:
    """Measure peak RSS while contracting every MO on a large grid."""

    params = (REPRESENTATIVE_EXAMPLES, (50, 75))
    param_names = ['molecule', 'edge_size']
    timeout = 180

    def setup(self, molecule: str, edge_size: int) -> None:
        """Precompute the GTO input for the all-MO contraction."""
        self.tabulator = Tabulator(example_source(molecule))
        axis = grid_axis(edge_size)
        self.tabulator.cartesian_grid(axis, axis, axis)

    def peakmem_tabulate_all_mos(self, molecule: str, edge_size: int) -> None:
        """Contract all MOs while ASV samples peak resident memory."""
        self.tabulator.tabulate_mos(None)
