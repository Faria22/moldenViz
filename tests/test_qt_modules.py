"""Compatibility checks for the focused Qt implementation modules."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from unittest.mock import Mock

import pytest

from moldenViz.qt.export import QFileDialog, export_data_dialog, export_image_dialog


def test_qt_public_imports_preserve_lazy_config_and_backend_loading() -> None:
    """The package facade should retain exports without initializing rendering."""
    script = """
import sys
import moldenViz.qt as qt
from moldenViz.qt import OrbitalViewer

assert qt.__all__ == ['OrbitalViewer', 'ViewerConfig']
assert qt.OrbitalViewer is OrbitalViewer
assert 'moldenViz._config_module' not in sys.modules
for module in ('pyvista', 'pyvistaqt', 'pydantic', 'matplotlib'):
    assert module not in sys.modules, module

from moldenViz.qt import ViewerConfig
from moldenViz._config_module import MainConfig
assert ViewerConfig is MainConfig
assert qt.ViewerConfig is ViewerConfig

from moldenViz.qt import *
assert OrbitalViewer is qt.OrbitalViewer
assert ViewerConfig is qt.ViewerConfig

try:
    qt.missing_attribute
except AttributeError as exc:
    assert "module 'moldenViz.qt' has no attribute 'missing_attribute'" == str(exc)
else:
    raise AssertionError('Unknown facade attributes must raise AttributeError')
"""
    subprocess.run([sys.executable, '-c', script], check=True)


@pytest.mark.parametrize('kind', ['data', 'image'])
def test_cancelled_export_dialog_does_not_export(kind: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """Cancelling either native dialog should leave viewer exports untouched."""
    parent = Mock()
    viewer = Mock(current_orbital_index=0)
    on_error = Mock()
    choose_path = Mock(return_value=('', ''))
    monkeypatch.setattr(QFileDialog, 'getSaveFileName', choose_path)

    dialog = export_data_dialog if kind == 'data' else export_image_dialog
    dialog(parent, viewer, {}, on_error)

    assert choose_path.call_args.args[0] is parent
    viewer.export_data.assert_not_called()
    viewer.export_image.assert_not_called()
    on_error.assert_not_called()


def test_data_dialog_preserves_scope_and_destination(monkeypatch: pytest.MonkeyPatch) -> None:
    """All-orbital exports should remain available without a selected orbital."""
    parent = Mock()
    viewer = Mock(current_orbital_index=-1)
    on_error = Mock()
    choose_path = Mock(return_value=('/tmp/orbitals.vtk', ''))
    monkeypatch.setattr(QFileDialog, 'getSaveFileName', choose_path)

    export_data_dialog(parent, viewer, {'format': 'vtk', 'scope': 'all'}, on_error)

    choose_path.assert_called_once_with(parent, 'Export orbital data', 'orbitals_all.vtk', 'VTK files (*.vtk)')
    viewer.export_data.assert_called_once_with('/tmp/orbitals.vtk', file_format='vtk', scope='all')
    on_error.assert_not_called()


def test_data_dialog_requires_selection_before_asking_for_destination(monkeypatch: pytest.MonkeyPatch) -> None:
    """Current-orbital exports should report missing selection through the host."""
    viewer = Mock(current_orbital_index=-1)
    on_error = Mock()
    choose_path = Mock()
    monkeypatch.setattr(QFileDialog, 'getSaveFileName', choose_path)

    export_data_dialog(Mock(), viewer, {}, on_error)

    choose_path.assert_not_called()
    viewer.export_data.assert_not_called()
    assert on_error.call_args.args[0] == 'Export failed'
    assert str(on_error.call_args.args[1]) == 'No orbital is currently selected.'


def test_image_dialog_preserves_transparency(monkeypatch: pytest.MonkeyPatch) -> None:
    """Image options should reach the viewer's explicit-path export API."""
    viewer = Mock()
    on_error = Mock()
    monkeypatch.setattr(QFileDialog, 'getSaveFileName', Mock(return_value=('/tmp/scene.png', '')))

    export_image_dialog(Mock(), viewer, {'format': 'png', 'transparent': True}, on_error)

    viewer.export_image.assert_called_once_with(Path('/tmp/scene.png'), file_format='png', transparent=True)
    on_error.assert_not_called()


@pytest.mark.parametrize('kind', ['data', 'image'])
def test_export_failure_is_reported_to_host(kind: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """Native dialogs should deliver writer failures to the owning host."""
    error = OSError('Destination is not writable')
    viewer = Mock(current_orbital_index=0)
    viewer.export_data.side_effect = error
    viewer.export_image.side_effect = error
    on_error = Mock()
    monkeypatch.setattr(QFileDialog, 'getSaveFileName', Mock(return_value=('/tmp/export', '')))

    dialog = export_data_dialog if kind == 'data' else export_image_dialog
    dialog(Mock(), viewer, {}, on_error)

    on_error.assert_called_once_with('Export failed', error)
