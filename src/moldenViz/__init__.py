"""molden_viz - A package for visualizing and analysing Molden files."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING, Any

from .__about__ import __version__
from .models import Atom, GaussianPrimitive, MolecularOrbital, Shell
from .parser import Parser
from .tabulator import GridType, Tabulator

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
    from .plotter import Plotter as Plotter


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
    if name == 'AtomType':
        module = import_module('moldenViz._config_module')
        atom_type_cls = module.AtomType
        globals()['AtomType'] = atom_type_cls
        return atom_type_cls
    if name == 'Plotter':
        module = import_module('moldenViz.plotter')
        plotter_cls = module.Plotter
        globals()['Plotter'] = plotter_cls
        return plotter_cls
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
