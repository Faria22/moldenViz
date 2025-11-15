"""molden_viz - A package for visualizing and analysing Molden files."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any

from .parser import Parser
from .tabulator import Tabulator

__all__ = ['Parser', 'Plotter', 'Tabulator']

if TYPE_CHECKING:  # pragma: no cover - type checking helper
    from .plotter import Plotter as Plotter
else:
    Plotter: Any | None = None


def __getattr__(name: str) -> Any:
    """Lazily import heavy modules such as Plotter.

    Parameters
    ----------
    name : str
        Attribute requested from the package namespace.

    Returns
    -------
    Any
        The requested attribute from the package.

    Raises
    ------
    AttributeError
        If the attribute is not defined.
    """
    if name == 'Plotter':
        module = import_module('moldenViz.plotter')
        plotter_cls = module.Plotter
        globals()['Plotter'] = plotter_cls
        return plotter_cls
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
