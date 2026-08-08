"""Tests for the Qt-native viewer and standalone facade."""
# ruff:file-ignore[undocumented-public-function]

from __future__ import annotations

import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING, cast
from unittest.mock import Mock

import numpy as np
import pytest
from PySide6.QtGui import QAction
from PySide6.QtWidgets import QAbstractSpinBox, QApplication, QWidget

import moldenViz.qt as qt_module
from moldenViz.plotter import Plotter
from moldenViz.qt import OrbitalViewer
from moldenViz.testing import NullInteractor, without_rendering

if TYPE_CHECKING:
    from collections.abc import Iterator

MOLDEN_PATH = Path(__file__).parent / 'sample_molden.inp'


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
def null_interactor(qapplication: QApplication) -> Iterator[None]:
    """Use the supported non-rendering context for widget tests."""
    del qapplication
    with without_rendering():
        yield


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
    assert isinstance(viewer.interactor, NullInteractor)

    viewer.close()


def test_viewer_can_hide_and_restore_builtin_controls() -> None:
    viewer = OrbitalViewer(show_controls=False)

    assert not viewer.controls_visible
    assert viewer.controls.isHidden()

    viewer.set_controls_visible(True)

    assert viewer.controls_visible
    assert not viewer.controls.isHidden()
    viewer.close()


def test_orbital_columns_fit_their_contents_with_padding() -> None:
    viewer = OrbitalViewer(str(MOLDEN_PATH), only_molecule=False)
    table = viewer.controls.orbital_table

    assert all(table.columnWidth(column) > table.sizeHintForColumn(column) for column in range(table.columnCount()))
    viewer.close()


def test_orbital_navigation_is_above_table() -> None:
    viewer = OrbitalViewer()
    layout = viewer.controls.orbitals_tab.layout()

    assert layout is not None
    navigation_item = layout.itemAt(1)
    table_item = layout.itemAt(2)
    assert navigation_item is not None
    assert table_item is not None
    navigation = navigation_item.layout()
    assert navigation is not None
    previous_item = navigation.itemAt(0)
    assert previous_item is not None
    assert previous_item.widget() is viewer.controls.previous_button
    assert table_item.widget() is viewer.controls.orbital_table
    viewer.close()


def test_control_tabs_follow_workflow_order() -> None:
    viewer = OrbitalViewer()

    assert [viewer.controls.tabs.tabText(index) for index in range(viewer.controls.tabs.count())] == [
        'Orbitals',
        'Appearance',
        'Grid',
        'Export',
    ]
    viewer.close()


def test_cartesian_grid_uses_column_headers_and_axis_row_labels() -> None:
    viewer = OrbitalViewer()
    controls = viewer.controls

    assert [label.text() for label in controls.cartesian_column_labels] == ['Min', 'Max', 'Num points']
    assert [controls.cartesian_axis_labels[axis].text() for axis in 'xyz'] == ['X', 'Y', 'Z']

    controls.grid_type.setCurrentText('cartesian')
    assert not controls.cartesian_grid.isHidden()
    assert all(not controls.cartesian_axis_labels[axis].isHidden() for axis in 'xyz')
    viewer.close()


def test_appearance_control_steps_and_contour_buttons() -> None:
    viewer = OrbitalViewer()
    controls = viewer.controls

    assert controls.mo_opacity.singleStep() == pytest.approx(0.1)
    assert controls.molecule_opacity.singleStep() == pytest.approx(0.1)
    assert controls.bond_radius.singleStep() == pytest.approx(0.05)
    assert controls.contour.buttonSymbols() is QAbstractSpinBox.ButtonSymbols.NoButtons
    viewer.close()


def test_color_fields_follow_appearance_selections() -> None:
    viewer = OrbitalViewer()
    controls = viewer.controls

    controls.color_scheme.setCurrentText('coolwarm')
    assert controls.negative_color.text() == '#3b4cc0'
    assert controls.positive_color.text() == '#b40426'
    assert not controls.negative_color.isEnabled()
    assert not controls.positive_color.isEnabled()

    controls.color_scheme.setCurrentText('custom')
    assert controls.negative_color.isEnabled()
    assert controls.positive_color.isEnabled()

    controls.bond_color_type.setCurrentText('split')
    assert controls.bond_color_label.isHidden()
    assert controls.bond_color.isHidden()
    controls.bond_color_type.setCurrentText('uniform')
    assert not controls.bond_color_label.isHidden()
    assert not controls.bond_color.isHidden()
    viewer.close()


def test_background_color_presets_and_custom_field() -> None:
    viewer = OrbitalViewer()
    controls = viewer.controls

    assert controls.background_color_choice.currentData() == 'white'
    assert controls.background_color.isHidden()

    controls.background_color_choice.setCurrentText('Light gray')
    assert controls.selected_background_color == '#A0A0A0'

    controls.background_color_choice.setCurrentText('Dark gray')
    assert controls.selected_background_color == '#202124'

    controls.background_color_choice.setCurrentText('Custom')
    controls.background_color.setText('#123456')
    assert not controls.background_color.isHidden()
    assert controls.selected_background_color == '#123456'

    viewer.set_background_color('navy')
    assert controls.background_color_choice.currentText() == 'Custom'
    assert controls.background_color.text() == 'navy'
    viewer.close()


def test_export_scope_uses_descriptive_labels_and_stable_values() -> None:
    viewer = OrbitalViewer()
    scope = viewer.controls.data_scope

    assert [scope.itemText(index) for index in range(scope.count())] == [
        'current orbital',
        'all orbitals',
    ]
    assert [scope.itemData(index) for index in range(scope.count())] == ['current', 'all']

    scope.setCurrentIndex(scope.findData('all'))
    viewer.controls.data_format.setCurrentText('cube')
    assert scope.currentData() == 'current'
    viewer.close()


def test_export_buttons_warn_without_host_handler(caplog: pytest.LogCaptureFixture) -> None:
    viewer = OrbitalViewer()

    assert not viewer.has_export_handler
    with caplog.at_level(logging.WARNING, logger='moldenViz.qt'):
        viewer.controls.image_export_button.click()
    assert 'without an export_requested receiver' in caplog.text

    requests: list[tuple[str, object]] = []
    viewer.export_requested.connect(lambda kind, options: requests.append((kind, options)))
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger='moldenViz.qt'):
        viewer.controls.image_export_button.click()

    assert viewer.has_export_handler
    assert requests
    assert requests[-1][0] == 'image'
    assert 'without an export_requested receiver' not in caplog.text
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
    assert viewer.current_orbital_index == -1
    assert viewer.molecular_orbitals == tuple(viewer.tabulator.molecular_orbitals)
    assert viewer.controls.current_mo_ind == -1
    viewer.close()


def test_partial_appearance_updates_do_not_require_controls() -> None:
    viewer = OrbitalViewer(show_controls=False)
    contour = 0.25
    opacity = 0.7

    viewer.update_appearance(contour=contour, mo_opacity=opacity, background_color='#123456')

    assert viewer.config.mo.contour == pytest.approx(contour)
    assert viewer.config.mo.opacity == pytest.approx(opacity)
    assert viewer.config.molecule.opacity == pytest.approx(1.0)
    assert viewer.interactor.background == '#123456'
    assert viewer.controls.contour.value() == pytest.approx(contour)
    with pytest.raises(ValueError, match='between 0 and 1'):
        viewer.update_appearance(mo_opacity=2.0)
    viewer.close()


def test_dashboard_grid_convenience_methods(monkeypatch: pytest.MonkeyPatch) -> None:
    viewer = OrbitalViewer(show_controls=False)
    viewer.tabulator = Mock()
    update_grid = Mock()
    x_points = 3
    monkeypatch.setattr(viewer, 'update_grid', update_grid)

    viewer.set_spherical_grid(radius=4.0, radial_points=x_points, theta_points=4, phi_points=5)

    spherical_axes, spherical_type = update_grid.call_args.args
    assert spherical_type == qt_module.GridType.SPHERICAL
    np.testing.assert_allclose(spherical_axes[0], [0.0, 2.0, 4.0])
    assert tuple(map(len, spherical_axes)) == (3, 4, 5)
    assert viewer.config.grid.default_type == 'spherical'
    assert viewer.config.grid.spherical.num_r_points == x_points

    viewer.set_cartesian_grid(
        x_range=(-1.0, 1.0),
        y_range=(-2.0, 2.0),
        z_range=(-3.0, 3.0),
        x_points=x_points,
        y_points=4,
        z_points=5,
    )

    cartesian_axes, cartesian_type = update_grid.call_args.args
    assert cartesian_type == qt_module.GridType.CARTESIAN
    np.testing.assert_allclose(cartesian_axes[0], [-1.0, 0.0, 1.0])
    assert tuple(map(len, cartesian_axes)) == (3, 4, 5)
    assert viewer.config.grid.default_type == 'cartesian'
    assert viewer.config.grid.cartesian.num_x_points == x_points
    with pytest.raises(ValueError, match='minimum'):
        viewer.set_cartesian_grid(
            x_range=(1.0, -1.0),
            y_range=(-1.0, 1.0),
            z_range=(-1.0, 1.0),
            x_points=3,
            y_points=3,
            z_points=3,
        )
    viewer.close()


def test_grid_update_requires_a_source() -> None:
    viewer = OrbitalViewer(show_controls=False)

    with pytest.raises(RuntimeError, match='Load a source'):
        viewer.set_spherical_grid(radius=4.0, radial_points=3, theta_points=4, phi_points=5)

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
    viewer = OrbitalViewer(
        str(MOLDEN_PATH),
        config={
            'grid': {
                'spherical': {'num_r_points': 3, 'num_theta_points': 3, 'num_phi_points': 3},
            },
        },
        show_controls=False,
    )
    viewer.wait_for_gtos(timeout=5)
    viewer.show_orbital(0)
    viewer.tabulator = Mock()

    viewer.export_data(tmp_path / 'orbital', file_format='vtk', scope='current')
    viewer.export_image(tmp_path / 'scene', file_format='png', transparent=True)

    viewer.tabulator.export.assert_called_once_with(tmp_path / 'orbital.vtk', mo_index=0)
    assert viewer.interactor.saved_screenshot == (tmp_path / 'scene.png', True)
    viewer.close()


def test_cube_export_rejects_all_orbitals(tmp_path: Path) -> None:
    viewer = OrbitalViewer()
    viewer.tabulator = Mock()

    with pytest.raises(ValueError, match='only supports one'):
        viewer.export_data(tmp_path / 'all.cube', file_format='cube', scope='all')

    viewer.close()


def test_plotter_returns_inside_existing_application() -> None:
    window = Plotter(str(MOLDEN_PATH), only_molecule=True)

    assert window.viewer.isVisible()
    assert not window._owns_application  # ruff: ignore[private-member-access]
    window.close()


def test_plotter_menus_follow_reordered_tabs_and_export_values(monkeypatch: pytest.MonkeyPatch) -> None:
    window = Plotter(str(MOLDEN_PATH), only_molecule=True)
    actions = {action.text(): action for action in window.findChildren(QAction)}

    actions['Grid'].setEnabled(True)
    actions['Grid'].trigger()
    assert window.viewer.controls.tabs.currentWidget() is window.viewer.controls.grid_tab
    actions['Appearance'].trigger()
    assert window.viewer.controls.tabs.currentWidget() is window.viewer.controls.appearance_tab

    window.viewer.controls.data_scope.setCurrentIndex(window.viewer.controls.data_scope.findData('all'))
    handle_export_request = Mock()
    monkeypatch.setattr(window, '_handle_export_request', handle_export_request)
    actions['Orbital data…'].setEnabled(True)
    actions['Orbital data…'].trigger()

    assert handle_export_request.call_args.args[1]['scope'] == 'all'
    window.close()


def test_plotter_owns_event_loop_when_it_creates_application() -> None:
    script = f"""
from PySide6.QtCore import QTimer
from moldenViz.plotter import Plotter
from moldenViz.testing import without_rendering

original_initialize = Plotter._initialize_viewer
def initialize_then_close(window, *args, **kwargs):
    assert window.isVisible()
    assert window.centralWidget().objectName() == 'moldenVizLoadingPlaceholder'
    original_initialize(window, *args, **kwargs)
    QTimer.singleShot(0, window.close)
Plotter._initialize_viewer = initialize_then_close
with without_rendering():
    window = Plotter({str(MOLDEN_PATH)!r}, only_molecule=True)
assert window._owns_application
assert window.viewer is not None
"""
    environment = os.environ.copy()
    environment['QT_QPA_PLATFORM'] = 'offscreen'
    subprocess.run([sys.executable, '-c', script], check=True, env=environment)


def test_offscreen_viewer_fails_cleanly_without_testing_context() -> None:
    script = """
from PySide6.QtWidgets import QApplication
from moldenViz.qt import OrbitalViewer

application = QApplication([])
try:
    OrbitalViewer()
except RuntimeError as exc:
    assert 'without_rendering' in str(exc)
else:
    raise AssertionError('Expected the offscreen safety guard')
"""
    environment = os.environ.copy()
    environment['QT_QPA_PLATFORM'] = 'offscreen'

    result = subprocess.run([sys.executable, '-c', script], env=environment, check=False)

    assert result.returncode == 0


def test_qt_gui_import_path_does_not_load_tkinter() -> None:
    assert 'tkinter' not in sys.modules
