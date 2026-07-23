"""Tests for the supported moldenViz API boundary."""

from __future__ import annotations

import os
import subprocess
import sys
from typing import TYPE_CHECKING

import moldenViz
import moldenViz.models as models_module
import moldenViz.parser as parser_module
import moldenViz.tabulator as tabulator_module
from moldenViz import Atom, AtomType, GaussianPrimitive, MolecularOrbital, Shell, examples

if TYPE_CHECKING:
    from pathlib import Path

EXPECTED_ROOT_API = [
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


def test_root_all_defines_supported_api() -> None:
    """The root export list should be explicit and stable."""
    assert moldenViz.__all__ == EXPECTED_ROOT_API


def test_public_models_are_parser_result_types() -> None:
    """Parser-facing models should be available from supported import paths."""
    assert parser_module.Atom is Atom
    assert parser_module.GaussianPrimitive is GaussianPrimitive
    assert parser_module.MolecularOrbital is MolecularOrbital
    assert parser_module.Shell is Shell


def test_removed_v1_names_are_not_public() -> None:
    """The v2 API should not retain compatibility aliases for removed names."""
    for name in ('_Atom', '_GTO', '_MolecularOrbital', '_Shell'):
        assert not hasattr(parser_module, name)
    assert not hasattr(tabulator_module, 'array_like_type')
    assert not hasattr(examples, 'all_examples')


def test_public_module_all_values_are_explicit() -> None:
    """Public modules should export only supported project-owned names."""
    assert models_module.__all__ == ['Atom', 'GaussianPrimitive', 'MolecularOrbital', 'Shell']
    assert parser_module.__all__ == [
        'BOHR_PER_ANGSTROM',
        'Atom',
        'GaussianPrimitive',
        'MolecularOrbital',
        'Parser',
        'Shell',
    ]
    assert tabulator_module.__all__ == ['GridType', 'Tabulator']
    assert examples.__all__ == ['acrolein', 'benzene', 'co', 'co2', 'furan', 'h2o', 'o2', 'prismane', 'pyridine']


def test_root_import_is_lightweight_and_read_only(tmp_path: Path) -> None:
    """Importing the root models must not load GUI modules or create user config."""
    script = """
import importlib.abc
import pathlib
import sys

gui_modules = frozenset({
    'matplotlib',
    'pyvista',
    'pyvistaqt',
    'PySide6',
    'pydantic',
    'qtpy',
    'toml',
})

class CoreOnlyImportBlocker(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.partition('.')[0] in gui_modules:
            raise ModuleNotFoundError(f'blocked GUI dependency: {fullname}')
        return None

sys.meta_path.insert(0, CoreOnlyImportBlocker())

from moldenViz import Parser, Tabulator
from moldenViz.models import Atom

assert Parser.__module__ == 'moldenViz.parser'
assert Tabulator.__module__ == 'moldenViz.tabulator'
assert Atom.__module__ == 'moldenViz.models'
assert 'moldenViz.plotter' not in sys.modules
assert 'moldenViz._plotter_ui' not in sys.modules
assert 'pyvista' not in sys.modules
assert 'pydantic' not in sys.modules
assert not (pathlib.Path.home() / '.config' / 'moldenViz').exists()
"""
    env = os.environ.copy()
    env['HOME'] = str(tmp_path)
    subprocess.run([sys.executable, '-c', script], check=True, env=env)


def test_atom_type_is_public_only_from_package_root() -> None:
    """The GUI model should not be exposed alongside parser result models."""
    assert AtomType.__module__ == 'moldenViz._config_module'
    assert not hasattr(models_module, 'AtomType')
