"""Benchmarks for structured Cartesian grid creation."""
# ruff:file-ignore[import-private-name]

import numpy as np
import pyvista as pv

from moldenViz import examples
from moldenViz._adaptive_grid import crossed_cell_ids, refined_grid
from moldenViz.tabulator import Tabulator


class TimeGridCreation:
    """Measure coordinate-grid construction independently of GTO work."""

    params = (10, 25, 50, 100)
    param_names = ['edge_size']
    number = 1
    repeat = (3, 10, 1.0)

    def setup(self, edge_size: int) -> None:
        """Create reusable parser state and grid axes."""
        self.tabulator = Tabulator(examples.co)
        self.axis = np.linspace(-3.0, 3.0, edge_size)

    def time_create_cartesian_grid(self, edge_size: int) -> None:
        """Create an ``edge_size ** 3`` Cartesian grid."""
        self.tabulator.cartesian_grid(
            self.axis,
            self.axis,
            self.axis,
            tabulate_gtos=False,
        )


class TimeInitialGridTabulation:
    """Compare initial uniform and adaptive-coarse GTO tabulation."""

    params = ('uniform-100', 'adaptive-coarse-21')
    param_names = ['mode']
    number = 1
    repeat = (3, 5, 1.0)

    def setup(self, mode: str) -> None:
        """Create the requested initial Cartesian grid."""
        edge_size = 100 if mode == 'uniform-100' else 21
        self.tabulator = Tabulator(examples.co)
        axis = np.linspace(-5.0, 5.0, edge_size)
        self.tabulator.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    def time_initial_gto_tabulation(self, mode: str) -> None:
        """Tabulate GTOs for the initial plotter grid."""
        self.tabulator.compute_gtos(self.tabulator.grid)


class TimeAdaptiveGridPreparation:
    """Measure the full all-MO adaptive-cache preparation stage."""

    number = 1
    repeat = (3, 5, 1.0)

    def setup(self) -> None:
        """Prepare a default coarse CO grid and its orbital values."""
        self.tabulator = Tabulator(examples.co)
        axis = np.linspace(-5.0, 5.0, 21)
        self.tabulator.cartesian_grid(axis, axis, axis)
        self.coarse_grid = pv.StructuredGrid()
        self.coarse_grid.points = self.tabulator.grid
        self.coarse_grid.dimensions = self.tabulator.grid_dimensions[::-1]
        self.coarse_mos = self.tabulator.tabulate_mos()

    def time_prepare_adaptive_cache(self) -> None:
        """Find all crossed cells, refine them, and tabulate fine GTOs."""
        cell_ids = crossed_cell_ids(self.coarse_grid, self.coarse_mos, 0.1)
        fine_grid = refined_grid(self.coarse_grid, cell_ids, 5.0)
        self.tabulator.compute_gtos(np.asarray(fine_grid.points))
