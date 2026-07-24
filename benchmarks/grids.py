"""Benchmarks for structured Cartesian grid creation."""

import numpy as np

from moldenViz import examples
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
