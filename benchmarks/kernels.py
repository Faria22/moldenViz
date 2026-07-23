"""Benchmarks for computational kernels used during tabulation."""

import numpy as np

from moldenViz.tabulator import Tabulator


class TimeSolidHarmonics:
    """Measure Cartesian real solid-harmonic scaling."""

    params = ((10_000, 125_000, 1_000_000), (0, 2, 4))
    param_names = ['num_points', 'lmax']
    number = 1
    repeat = (3, 5, 1.0)
    timeout = 180

    def setup(self, num_points: int, lmax: int) -> None:
        """Create deterministic Cartesian points."""
        rng = np.random.default_rng(seed=8300)
        self.points = rng.uniform(-6.0, 6.0, size=(num_points, 3))

    def time_real_solid_harmonics(self, num_points: int, lmax: int) -> None:
        """Evaluate every real solid harmonic through ``lmax``."""
        Tabulator._tabulate_real_solid_harmonics(self.points, lmax)  # ruff:ignore[private-member-access]
