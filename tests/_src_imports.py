"""Utilities to import project modules directly from ``src`` during testing."""

from __future__ import annotations

import sys
from importlib import import_module
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

parser_module = import_module("moldenViz.parser")
plotter_module = import_module("moldenViz.plotter")
tabulator_module = import_module("moldenViz.tabulator")

# parser exports
Parser = parser_module.Parser
_GTO = parser_module._GTO
_Shell = parser_module._Shell

# plotter exports
Plotter = plotter_module.Plotter

# tabulator exports
GridType = tabulator_module.GridType
Tabulator = tabulator_module.Tabulator
_cartesian_to_spherical = tabulator_module._cartesian_to_spherical
_spherical_to_cartesian = tabulator_module._spherical_to_cartesian
array_like_type = tabulator_module.array_like_type

__all__ = [
    "Parser",
    "Plotter",
    "Tabulator",
    "GridType",
    "_GTO",
    "_Shell",
    "_cartesian_to_spherical",
    "_spherical_to_cartesian",
    "parser_module",
    "plotter_module",
    "tabulator_module",
    "array_like_type",
]
