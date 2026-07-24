"""Peak-memory benchmarks for representative tabulation workloads."""

from moldenViz.tabulator import Tabulator

from ._shared import POINT_CHUNK_SIZES, REPRESENTATIVE_EXAMPLES, example_source, grid_axis


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
