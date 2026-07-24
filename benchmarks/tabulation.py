"""Runtime benchmarks for GTO and MO tabulation."""

from moldenViz.tabulator import Tabulator

from ._shared import (
    EXAMPLE_NAMES,
    GRID_EDGES,
    MO_SELECTIONS,
    POINT_CHUNK_SIZES,
    PYSCF_CO_SPHERICAL,
    REPRESENTATIVE_EXAMPLES,
    WORKER_COUNTS,
    GenericSolidHarmonicTabulator,
    MOSelection,
    example_source,
    grid_axis,
    mo_indices,
)


class TimeHighAngularMomentumGTOTabulation:
    """Compare generated and generic kernels in a cc-pVQZ tabulation."""

    params = ((25, 50), ('generated', 'generic'))
    param_names = ['edge_size', 'implementation']
    number = 1
    repeat = (3, 5, 1.0)
    timeout = 180

    def setup(self, edge_size: int, implementation: str) -> None:
        """Create a sequential cc-pVQZ tabulator with an uncomputed grid."""
        tabulator_class = GenericSolidHarmonicTabulator if implementation == 'generic' else Tabulator
        self.tabulator = tabulator_class(str(PYSCF_CO_SPHERICAL), max_workers=1)
        axis = grid_axis(edge_size)
        self.tabulator.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    def time_tabulate_gtos(self, edge_size: int, implementation: str) -> None:
        """Tabulate every cc-pVQZ basis function."""
        self.tabulator.tabulate_gtos()


class TimeGTOTabulation:
    """Measure GTO scaling across every example molecule and grid size."""

    params = (EXAMPLE_NAMES, GRID_EDGES)
    param_names = ['molecule', 'edge_size']
    number = 1
    repeat = (3, 5, 1.0)
    timeout = 180

    def setup(self, molecule: str, edge_size: int) -> None:
        """Create a tabulator with an uncomputed Cartesian grid."""
        self.tabulator = Tabulator(example_source(molecule))
        axis = grid_axis(edge_size)
        self.tabulator.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    def time_tabulate_gtos(self, molecule: str, edge_size: int) -> None:
        """Tabulate every basis function on the selected grid."""
        self.tabulator.tabulate_gtos()


class TimeGTOChunkSizes:
    """Measure the runtime tradeoff across bounded point chunk sizes."""

    params = (REPRESENTATIVE_EXAMPLES, (50, 100), POINT_CHUNK_SIZES)
    param_names = ['molecule', 'edge_size', 'point_chunk_size']
    number = 1
    repeat = (3, 5, 1.0)
    timeout = 180

    def setup(self, molecule: str, edge_size: int, point_chunk_size: int | None) -> None:
        """Create a tabulator with an uncomputed Cartesian grid."""
        self.tabulator = Tabulator(example_source(molecule))
        axis = grid_axis(edge_size)
        self.tabulator.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    def time_tabulate_gtos(self, molecule: str, edge_size: int, point_chunk_size: int | None) -> None:
        """Tabulate GTOs with the selected point-chunk bound."""
        self.tabulator.tabulate_gtos(point_chunk_size=point_chunk_size)


class TimeGTOWorkerScaling:
    """Compare sequential and bounded-parallel GTO tabulation."""

    params = (REPRESENTATIVE_EXAMPLES, (50, 100), WORKER_COUNTS)
    param_names = ['molecule', 'edge_size', 'max_workers']
    number = 1
    repeat = (3, 5, 1.0)
    timeout = 180

    def setup(self, molecule: str, edge_size: int, max_workers: int) -> None:
        """Create an uncomputed grid with an explicit worker limit."""
        self.tabulator = Tabulator(example_source(molecule), max_workers=max_workers)
        axis = grid_axis(edge_size)
        self.tabulator.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    def time_tabulate_gtos(self, molecule: str, edge_size: int, max_workers: int) -> None:
        """Tabulate GTOs with the selected concurrency."""
        self.tabulator.tabulate_gtos()


class TimeMOContraction:
    """Measure single, several, and all-MO contractions."""

    params = (EXAMPLE_NAMES, GRID_EDGES, MO_SELECTIONS)
    param_names = ['molecule', 'edge_size', 'mo_selection']
    number = 1
    repeat = (3, 5, 1.0)
    timeout = 180

    def setup(self, molecule: str, edge_size: int, mo_selection: MOSelection) -> None:
        """Precompute GTOs so the timed region contains only MO contraction."""
        self.tabulator = Tabulator(example_source(molecule))
        axis = grid_axis(edge_size)
        self.tabulator.cartesian_grid(axis, axis, axis)
        self.indices = mo_indices(self.tabulator, mo_selection)

    def time_tabulate_mos(self, molecule: str, edge_size: int, mo_selection: MOSelection) -> None:
        """Contract the requested MO group from cached GTO values."""
        self.tabulator.tabulate_mos(self.indices)
