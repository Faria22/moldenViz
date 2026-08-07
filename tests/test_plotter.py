"""Tests for the Qt-native viewer and standalone facade."""
# ruff:file-ignore[no-self-use, undocumented-public-function, undocumented-public-method]

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from typing import cast
from unittest.mock import Mock

import pytest
from PySide6.QtWidgets import QApplication, QWidget

import moldenViz.qt as qt_module
from moldenViz.plotter import Plotter
from moldenViz.qt import OrbitalViewer

MOLDEN_PATH = Path(__file__).parent / 'sample_molden.inp'


class FakeActor:
    """Small VTK actor stand-in."""

    def __init__(self) -> None:
        self.visible = True
        self.opacity = 1.0

    def GetVisibility(self) -> bool:  # ruff: ignore[invalid-function-name]
        return self.visible

    def SetVisibility(self, visible: bool) -> None:  # ruff: ignore[invalid-function-name]
        self.visible = visible

    def GetProperty(self) -> FakeActor:  # ruff: ignore[invalid-function-name]
        return self

    def SetOpacity(self, opacity: float) -> None:  # ruff: ignore[invalid-function-name]
        self.opacity = opacity


class FakeInteractor(QWidget):
    """Headless ``QtInteractor`` replacement for widget unit tests."""

    instances: list[FakeInteractor] = []

    def __init__(self, parent: QWidget | None = None, **_kwargs: object) -> None:
        super().__init__(parent)
        self.background = ''
        self.closed_count = 0
        self.actors: list[FakeActor] = []
        self.saved_graphic: Path | None = None
        self.saved_screenshot: tuple[Path, bool] | None = None
        self.reset_count = 0
        self.__class__.instances.append(self)

    def set_background(self, color: str) -> None:
        self.background = color

    def show_axes(self) -> None:
        return

    def add_mesh(self, _mesh: object, **_kwargs: object) -> FakeActor:
        actor = FakeActor()
        self.actors.append(actor)
        return actor

    def remove_actor(self, actor: FakeActor) -> None:
        if actor in self.actors:
            self.actors.remove(actor)

    def update(self) -> None:
        return

    def close(self) -> None:  # type: ignore[override]
        self.closed_count += 1

    def save_graphic(self, path: Path) -> None:
        self.saved_graphic = path

    def screenshot(self, path: Path, *, transparent_background: bool) -> None:
        self.saved_screenshot = path, transparent_background

    def reset_camera(self) -> None:
        self.reset_count += 1


@pytest.fixture(scope='session')
def qapplication() -> QApplication:
    """Provide the host-owned Qt application required by widgets.

    Returns
    -------
    QApplication
        Process-wide application used by viewer tests.
    """
    application = QApplication.instance()
    return QApplication([]) if application is None else cast(QApplication, application)


@pytest.fixture(autouse=True)
def fake_interactor(monkeypatch: pytest.MonkeyPatch, qapplication: QApplication) -> None:
    """Avoid creating native VTK windows in unit tests."""
    del qapplication
    FakeInteractor.instances.clear()
    monkeypatch.setattr(qt_module, 'QtInteractor', FakeInteractor)


def test_viewer_requires_existing_qapplication(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(qt_module.QApplication, 'instance', lambda: None)
    with pytest.raises(RuntimeError, match='existing QApplication'):
        OrbitalViewer()


def test_viewer_is_parentable_and_does_not_show_itself(qapplication: QApplication) -> None:
    parent = QWidget()
    viewer = OrbitalViewer(parent=parent)

    assert viewer.parent() is parent
    assert viewer.isAncestorOf(viewer.interactor)
    assert not viewer.isVisible()
    assert QApplication.instance() is qapplication

    viewer.close()


def test_viewers_have_isolated_configuration() -> None:
    first = OrbitalViewer(config={'background_color': '#112233'})
    second = OrbitalViewer(config={'background_color': '#abcdef'})

    first.set_background_color('red')

    assert first.config.background_color == 'red'
    assert second.config.background_color == '#abcdef'
    assert first.interactor.background == 'red'
    first.close()
    second.close()


def test_deferred_molecule_source_loading() -> None:
    viewer = OrbitalViewer(config={'molecule': {'bond': {'show': False}}})
    ready_calls: list[bool] = []
    viewer.source_ready.connect(lambda: ready_calls.append(True))

    viewer.set_source(str(MOLDEN_PATH), only_molecule=True)

    assert viewer.tabulator.atoms
    assert viewer.gtos_ready
    assert ready_calls == [True]
    assert not viewer.controls.tabs.isTabVisible(0)
    viewer.close()


def test_orbital_source_tabulation_and_selection() -> None:
    viewer = OrbitalViewer(
        str(MOLDEN_PATH),
        config={
            'grid': {
                'spherical': {'num_r_points': 3, 'num_theta_points': 3, 'num_phi_points': 3},
            },
        },
    )
    viewer.wait_for_gtos(timeout=5)
    changes: list[int] = []
    viewer.orbital_changed.connect(changes.append)

    viewer.show_orbital(0)
    viewer.show_orbital(-1)

    assert viewer.gtos_ready
    assert changes == [0, -1]
    assert viewer.controls.current_mo_ind == -1
    viewer.close()


def test_set_source_rejects_closed_viewer() -> None:
    viewer = OrbitalViewer()
    viewer.close()

    with pytest.raises(RuntimeError, match='closed'):
        viewer.set_source(str(MOLDEN_PATH), only_molecule=True)


def test_close_releases_interactor_once() -> None:
    viewer = OrbitalViewer()
    interactor = viewer.interactor

    viewer.close()
    viewer.close()

    assert interactor.closed_count == 1


def test_errors_are_emitted_without_dialogs() -> None:
    viewer = OrbitalViewer()
    errors: list[tuple[str, object]] = []
    viewer.error_occurred.connect(lambda title, exc: errors.append((title, exc)))

    viewer.report_error('Example failure', ValueError('bad input'))

    assert errors
    assert errors[0][0] == 'Example failure'
    assert isinstance(errors[0][1], ValueError)
    viewer.close()


def test_exports_use_explicit_paths(tmp_path: Path) -> None:
    viewer = OrbitalViewer()
    viewer.tabulator = Mock()
    viewer.controls.current_mo_ind = 2

    viewer.export_data(tmp_path / 'orbital', file_format='vtk', scope='current')
    viewer.export_image(tmp_path / 'scene', file_format='png', transparent=True)

    viewer.tabulator.export.assert_called_once_with(tmp_path / 'orbital.vtk', mo_index=2)
    assert viewer.interactor.saved_screenshot == (tmp_path / 'scene.png', True)
    viewer.close()


def test_cube_export_rejects_all_orbitals(tmp_path: Path) -> None:
    viewer = OrbitalViewer()
    viewer.tabulator = Mock()

    with pytest.raises(ValueError, match='only supports one'):
        viewer.export_data(tmp_path / 'all.cube', file_format='cube', scope='all')

    viewer.close()


def test_plotter_returns_inside_existing_application(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(qt_module, 'QtInteractor', FakeInteractor)
    window = Plotter(str(MOLDEN_PATH), only_molecule=True)

    assert window.viewer.isVisible()
    assert not window._owns_application  # ruff: ignore[private-member-access]
    window.close()


def test_plotter_owns_event_loop_when_it_creates_application() -> None:
    script = f"""
from PySide6.QtCore import QTimer
import moldenViz.qt as qt_module
from moldenViz.plotter import Plotter
from tests.test_plotter import FakeInteractor

qt_module.QtInteractor = FakeInteractor
original_show = Plotter.show
def show_then_close(window):
    original_show(window)
    QTimer.singleShot(0, window.close)
Plotter.show = show_then_close
window = Plotter({str(MOLDEN_PATH)!r}, only_molecule=True)
assert window._owns_application
"""
    environment = os.environ.copy()
    environment['QT_QPA_PLATFORM'] = 'offscreen'
    subprocess.run([sys.executable, '-c', script], check=True, env=environment)


def test_qt_gui_import_path_does_not_load_tkinter() -> None:
    assert 'tkinter' not in sys.modules
