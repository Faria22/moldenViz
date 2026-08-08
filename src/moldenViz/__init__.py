"""molden_viz - A package for visualizing and analysing Molden files."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any

from .__about__ import __version__

__all__ = [
    'Atom',
    'AtomType',
    'GaussianPrimitive',
    'GridType',
    'MolecularOrbital',
    'Parser',
    'Plotter',
    'Shell',
    'Tabulator',
    '__version__',
]

if TYPE_CHECKING:  # pragma: no cover - type checking helper
    from ._config_module import AtomType as AtomType
    from .models import Atom as Atom
    from .models import GaussianPrimitive as GaussianPrimitive
    from .models import MolecularOrbital as MolecularOrbital
    from .models import Shell as Shell
    from .parser import Parser as Parser
    from .plotter import Plotter as Plotter
    from .tabulator import GridType as GridType
    from .tabulator import Tabulator as Tabulator

_LAZY_IMPORTS = {  # ruff:ignore[non-empty-init-module]
    'Atom': ('moldenViz.models', 'Atom'),
    'AtomType': ('moldenViz._config_module', 'AtomType'),
    'GaussianPrimitive': ('moldenViz.models', 'GaussianPrimitive'),
    'GridType': ('moldenViz.tabulator', 'GridType'),
    'MolecularOrbital': ('moldenViz.models', 'MolecularOrbital'),
    'Parser': ('moldenViz.parser', 'Parser'),
    'Plotter': ('moldenViz.plotter', 'Plotter'),
    'Shell': ('moldenViz.models', 'Shell'),
    'Tabulator': ('moldenViz.tabulator', 'Tabulator'),
}


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
    target = _LAZY_IMPORTS.get(name)
    if target is not None:
        module_name, attribute_name = target
        value = getattr(import_module(module_name), attribute_name)
        globals()[name] = value
        return value
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')


def __dir__() -> list[str]:
    """Return module attributes, including exports that have not loaded yet.

    Returns
    -------
    list[str]
        Available package attributes and public lazy exports.
    """
    return sorted({*globals(), *__all__})
