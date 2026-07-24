"""Benchmarks for parsing bundled Molden inputs."""

from moldenViz.tabulator import Tabulator

from ._shared import EXAMPLE_NAMES, example_source


class TimeParsing:
    """Measure parser and model construction for every bundled molecule."""

    params = EXAMPLE_NAMES
    param_names = ['molecule']
    number = 1
    repeat = (3, 10, 1.0)

    def time_parse_example(self, molecule: str) -> None:
        """Parse a bundled Molden input."""
        Tabulator(example_source(molecule))
