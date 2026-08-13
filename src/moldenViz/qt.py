"""Qt-native embeddable molecular-orbital viewer."""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np
from PySide6.QtCore import QEvent, QObject, Qt, QTimer, Signal, Slot
from PySide6.QtWidgets import (
    QAbstractItemView,
    QAbstractSpinBox,
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from ._adaptive_grid import (
    AdaptiveScale,
    crossed_cell_ids,
    format_scale,
    normalize_scale,
    parse_scale,
    refined_grid,
)
from ._plotter_jobs import BackgroundJob
from ._plotter_rendering import _PlotterRendering
from .parser import _validate_molden_input
from .tabulator import GridType, Tabulator

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping
    from concurrent.futures import Future
    from os import PathLike

    import pyvista as pv
    from numpy.typing import NDArray
    from PySide6.QtGui import QCloseEvent, QShowEvent
    from pyvistaqt import QtInteractor

    from ._config_module import Config, MainConfig
    from .models import MolecularOrbital

    ViewerConfig = MainConfig

logger = logging.getLogger(__name__)

__all__ = ['OrbitalViewer', 'ViewerConfig']

_GTO_EXECUTOR = ThreadPoolExecutor(max_workers=1)
_MO_COLOR_SCHEMES = ['bwr', 'RdBu', 'seismic', 'coolwarm', 'PiYG']
_BACKGROUND_COLORS = (
    ('White', 'white'),
    ('Black', 'black'),
    ('Light gray', '#A0A0A0'),
    ('Dark gray', '#202124'),
)
_ORBITAL_COLUMN_PADDING = 8
_CUSTOM_COLOR_COUNT = 2
_MIN_ADAPTIVE_POINTS = 2


def __getattr__(name: str) -> Any:
    """Load the public configuration model only when it is requested.

    Returns
    -------
    Any
        Requested module attribute.
    """
    if name == 'ViewerConfig':
        viewer_config = import_module('moldenViz._config_module').MainConfig
        globals()[name] = viewer_config
        return viewer_config
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')


def _load_qt_interactor() -> type[QtInteractor]:
    """Load the VTK-backed Qt widget when a viewer is constructed.

    Returns
    -------
    type[QtInteractor]
        PyVista's Qt interactor class.
    """
    testing = import_module('moldenViz.testing')
    override = testing._get_interactor_override()  # ruff: ignore[private-member-access]
    if override is not None:
        return cast('type[QtInteractor]', override)
    application = QApplication.instance()
    if isinstance(application, QApplication) and application.platformName() == 'offscreen':
        raise RuntimeError(
            'OrbitalViewer cannot create a VTK QtInteractor with QT_QPA_PLATFORM=offscreen. '
            'Use moldenViz.testing.without_rendering() for UI smoke tests.',
        )
    return import_module('pyvistaqt').QtInteractor


def _is_color_like(color: object) -> bool:
    """Return whether Matplotlib accepts a color specification.

    Returns
    -------
    bool
        Whether ``color`` is valid.
    """
    return bool(import_module('matplotlib.colors').is_color_like(color))


def _has_colormap(name: str) -> bool:
    """Return whether Matplotlib provides a named colormap.

    Returns
    -------
    bool
        Whether ``name`` identifies an available colormap.
    """
    return name in import_module('matplotlib').colormaps


@dataclass(frozen=True)
class _GTOResult:
    """GTO values and the structured grid snapshot they describe."""

    grid: NDArray[np.floating]
    axes: tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]
    grid_type: GridType
    gtos: NDArray[np.floating]


@dataclass(frozen=True)
class _AdaptiveResult:
    """Refined contour grid and its reusable GTO cache."""

    mesh: pv.UnstructuredGrid
    gtos: NDArray[np.floating]
    crossed_cells: int


class _CompletionDispatcher(QObject):
    """Transfer callbacks from worker threads to the Qt object thread."""

    callback_ready = Signal(object)

    def __init__(self, parent: QObject) -> None:
        super().__init__(parent)
        self.callback_ready.connect(self._run_callback, Qt.ConnectionType.QueuedConnection)

    def dispatch(self, callback: Callable[[], None]) -> None:
        """Queue a callback on the owning Qt thread."""
        try:
            self.callback_ready.emit(callback)
        except RuntimeError:
            # The parent widget may have been deleted while non-cancellable
            # numerical work was still finishing.
            return

    @Slot(object)
    # Qt invokes this through a bound queued connection.
    # ruff: ignore[no-self-use]
    def _run_callback(self, callback: Callable[[], None]) -> None:
        callback()


class OrbitalControlPanel(QWidget):
    """Qt controls associated with an :class:`OrbitalViewer`."""

    def __init__(self, viewer: OrbitalViewer) -> None:
        super().__init__(viewer)
        self._viewer = viewer
        self.current_mo_ind = -1
        self.setMinimumWidth(310)
        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        self.tabs = QTabWidget(self)
        layout.addWidget(self.tabs)
        self._build_orbitals_tab()
        self._build_appearance_tab()
        self._build_grid_tab()
        self._build_export_tab()

    def _build_orbitals_tab(self) -> None:
        tab = QWidget(self)
        self.orbitals_tab = tab
        layout = QVBoxLayout(tab)
        self.loading_label = QLabel('', tab)
        self.loading_label.hide()
        layout.addWidget(self.loading_label)

        nav = QHBoxLayout()
        self.previous_button = QPushButton('Previous', tab)
        self.clear_button = QPushButton('Clear orbital', tab)
        self.next_button = QPushButton('Next', tab)
        self.previous_button.clicked.connect(self._previous_orbital)
        self.clear_button.clicked.connect(lambda: self._viewer.show_orbital(-1))
        self.next_button.clicked.connect(self._next_orbital)
        nav.addWidget(self.previous_button)
        nav.addWidget(self.clear_button)
        nav.addWidget(self.next_button)
        layout.addLayout(nav)

        self.orbital_table = QTableWidget(0, 5, tab)
        self.orbital_table.setHorizontalHeaderLabels(['#', 'Sym', 'Spin', 'Occ', 'Energy'])
        self.orbital_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.orbital_table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.orbital_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.orbital_table.verticalHeader().hide()
        self.orbital_table.itemSelectionChanged.connect(self._on_orbital_selected)
        layout.addWidget(self.orbital_table)
        self._fit_orbital_columns()

        self.tabs.addTab(tab, 'Orbitals')
        self.update_nav_button_states()

    def _build_grid_tab(self) -> None:
        tab = QWidget(self)
        self.grid_tab = tab
        layout = QVBoxLayout(tab)
        form = QFormLayout()
        self.grid_type = QComboBox(tab)
        self.grid_type.addItems(['spherical', 'cartesian', 'adaptive'])
        self.grid_type.currentTextChanged.connect(self._update_grid_field_visibility)
        form.addRow('Grid type', self.grid_type)

        self.radius = self._double_spin(0.001, 1_000.0, 5.0)
        self.radius_points = self._int_spin(1, 10_000, 100)
        self.theta_points = self._int_spin(1, 10_000, 60)
        self.phi_points = self._int_spin(1, 10_000, 120)
        spherical_fields = [
            ('Radius', self.radius),
            ('Radius points', self.radius_points),
            ('Theta points', self.theta_points),
            ('Phi points', self.phi_points),
        ]
        self._spherical_rows: list[tuple[QLabel, QWidget]] = []
        for label_text, widget in spherical_fields:
            label = QLabel(label_text, tab)
            form.addRow(label, widget)
            self._spherical_rows.append((label, widget))

        self.cartesian_grid = QWidget(tab)
        cartesian_layout = QGridLayout(self.cartesian_grid)
        cartesian_layout.setContentsMargins(0, 0, 0, 0)
        cartesian_layout.setColumnStretch(0, 0)
        for column in range(1, 4):
            cartesian_layout.setColumnStretch(column, 1)

        self.cartesian_column_labels = [QLabel(text, self.cartesian_grid) for text in ('Min', 'Max', 'Num points')]
        for column, label in enumerate(self.cartesian_column_labels, start=1):
            label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            cartesian_layout.addWidget(label, 0, column)

        self.cartesian_fields: dict[str, tuple[QDoubleSpinBox, QDoubleSpinBox, QSpinBox]] = {}
        self.cartesian_axis_labels: dict[str, QLabel] = {}
        for row, axis in enumerate('xyz', start=1):
            minimum = self._double_spin(-1_000.0, 1_000.0, -5.0)
            maximum = self._double_spin(-1_000.0, 1_000.0, 5.0)
            points = self._int_spin(1, 10_000, 100)
            label = QLabel(axis.upper(), self.cartesian_grid)
            label.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            cartesian_layout.addWidget(label, row, 0)
            cartesian_layout.addWidget(minimum, row, 1)
            cartesian_layout.addWidget(maximum, row, 2)
            cartesian_layout.addWidget(points, row, 3)
            self.cartesian_fields[axis] = (minimum, maximum, points)
            self.cartesian_axis_labels[axis] = label
        form.addRow(self.cartesian_grid)

        self.adaptive_scale_label = QLabel('Fine scale', tab)
        self.adaptive_scale = QLineEdit('5.0', tab)
        form.addRow(self.adaptive_scale_label, self.adaptive_scale)

        layout.addLayout(form)
        apply_button = QPushButton('Apply grid', tab)
        apply_button.clicked.connect(self._apply_grid)
        layout.addWidget(apply_button)
        layout.addStretch()
        self.tabs.addTab(tab, 'Grid')
        self._update_grid_field_visibility(self.grid_type.currentText())

    def _build_appearance_tab(self) -> None:
        tab = QWidget(self)
        self.appearance_tab = tab
        layout = QVBoxLayout(tab)

        self.orbital_group = QGroupBox('Molecular orbital', tab)
        orbital_form = QFormLayout(self.orbital_group)
        self.contour = self._double_spin(1e-6, 1e6, 0.1, decimals=6, show_buttons=False)
        self.mo_opacity = self._double_spin(0.0, 1.0, 1.0, step=0.1)
        self.color_scheme = QComboBox(tab)
        self.color_scheme.addItems([*_MO_COLOR_SCHEMES, 'custom'])
        self.negative_color = QLineEdit('blue', tab)
        self.positive_color = QLineEdit('red', tab)
        orbital_form.addRow('Contour', self.contour)
        orbital_form.addRow('Opacity', self.mo_opacity)
        orbital_form.addRow('Color scheme', self.color_scheme)
        orbital_form.addRow('Negative color', self.negative_color)
        orbital_form.addRow('Positive color', self.positive_color)
        self.color_scheme.currentTextChanged.connect(self._update_custom_color_visibility)
        layout.addWidget(self.orbital_group)

        molecule_group = QGroupBox('Molecule', tab)
        molecule_form = QFormLayout(molecule_group)
        self.molecule_opacity = self._double_spin(0.0, 1.0, 1.0, step=0.1)
        self.show_atoms = QCheckBox('Show atoms', tab)
        self.show_bonds = QCheckBox('Show bonds', tab)
        self.show_axes = QCheckBox('Show axes', tab)
        self.show_axes.toggled.connect(self._viewer.set_axes_visible)
        self.bond_max_length = self._double_spin(0.001, 1_000.0, 4.0)
        self.bond_radius = self._double_spin(0.001, 100.0, 0.15, step=0.05)
        self.bond_color_type = QComboBox(tab)
        self.bond_color_type.addItems(['uniform', 'split'])
        self.bond_color = QLineEdit('grey', tab)
        self.bond_color_type.currentTextChanged.connect(self._update_bond_color_visibility)
        molecule_form.addRow('Opacity', self.molecule_opacity)
        molecule_form.addRow(self.show_atoms)
        molecule_form.addRow(self.show_bonds)
        molecule_form.addRow(self.show_axes)
        molecule_form.addRow('Bond max length', self.bond_max_length)
        molecule_form.addRow('Bond radius', self.bond_radius)
        molecule_form.addRow('Bond colors', self.bond_color_type)
        self.bond_color_label = QLabel('Uniform color', tab)
        molecule_form.addRow(self.bond_color_label, self.bond_color)
        layout.addWidget(molecule_group)

        background_form = QFormLayout()
        self.background_color_choice = QComboBox(tab)
        for label, color in _BACKGROUND_COLORS:
            self.background_color_choice.addItem(label, color)
        self.background_color_choice.addItem('Custom', None)
        self.background_color = QLineEdit('white', tab)
        self.background_color_label = QLabel('Custom background', tab)
        background_form.addRow('Background color', self.background_color_choice)
        background_form.addRow(self.background_color_label, self.background_color)
        self.background_color_choice.currentIndexChanged.connect(self._update_background_color_visibility)
        layout.addLayout(background_form)

        buttons = QHBoxLayout()
        apply_button = QPushButton('Apply', tab)
        reset_button = QPushButton('Reset', tab)
        save_button = QPushButton('Save as defaults', tab)
        apply_button.clicked.connect(self._apply_appearance)
        reset_button.clicked.connect(self.sync_from_viewer)
        save_button.clicked.connect(self._save_settings)
        buttons.addWidget(apply_button)
        buttons.addWidget(reset_button)
        buttons.addWidget(save_button)
        layout.addLayout(buttons)
        layout.addStretch()
        self.tabs.addTab(tab, 'Appearance')
        self._update_bond_color_visibility(self.bond_color_type.currentText())
        self._update_background_color_visibility()

    def _build_export_tab(self) -> None:
        tab = QWidget(self)
        layout = QFormLayout(tab)
        self.data_format = QComboBox(tab)
        self.data_format.addItems(['vtk', 'cube'])
        self.data_scope = QComboBox(tab)
        self.data_scope.addItem('current orbital', 'current')
        self.data_scope.addItem('all orbitals', 'all')
        self.data_format.currentTextChanged.connect(self._update_export_scope)
        self.data_export_button = QPushButton('Export orbital data…', tab)
        self.data_export_button.clicked.connect(self._request_data_export)
        self.image_format = QComboBox(tab)
        self.image_format.addItems(['png', 'jpeg', 'svg', 'pdf'])
        self.transparent_background = QCheckBox('Transparent PNG background', tab)
        self.image_format.currentTextChanged.connect(self._update_transparency_option)
        self.image_export_button = QPushButton('Export image…', tab)
        self.image_export_button.clicked.connect(self._request_image_export)
        layout.addRow('Data format', self.data_format)
        layout.addRow('Scope', self.data_scope)
        layout.addRow(self.data_export_button)
        layout.addRow('Image format', self.image_format)
        layout.addRow(self.transparent_background)
        layout.addRow(self.image_export_button)
        self.tabs.addTab(tab, 'Export')

    @staticmethod
    def _double_spin(
        minimum: float,
        maximum: float,
        value: float,
        *,
        decimals: int = 3,
        step: float | None = None,
        show_buttons: bool = True,
    ) -> QDoubleSpinBox:
        widget = QDoubleSpinBox()
        widget.setRange(minimum, maximum)
        widget.setDecimals(decimals)
        if step is not None:
            widget.setSingleStep(step)
        if not show_buttons:
            widget.setButtonSymbols(QAbstractSpinBox.ButtonSymbols.NoButtons)
        widget.setValue(value)
        return widget

    @staticmethod
    def _int_spin(minimum: int, maximum: int, value: int) -> QSpinBox:
        widget = QSpinBox()
        widget.setRange(minimum, maximum)
        widget.setValue(value)
        return widget

    def set_molecular_orbitals(self, orbitals: list[MolecularOrbital]) -> None:
        """Populate the orbital table from parsed orbitals."""
        self.orbital_table.blockSignals(True)
        self.orbital_table.setRowCount(len(orbitals))
        for row, orbital in enumerate(orbitals):
            values = (str(row + 1), orbital.sym, orbital.spin, f'{orbital.occ:g}', f'{orbital.energy:.6f}')
            for column, value in enumerate(values):
                self.orbital_table.setItem(row, column, QTableWidgetItem(value))
        self.orbital_table.blockSignals(False)
        self._fit_orbital_columns()
        self.current_mo_ind = -1
        self.update_nav_button_states()

    def _fit_orbital_columns(self) -> None:
        self.orbital_table.resizeColumnsToContents()
        for column in range(self.orbital_table.columnCount()):
            width = self.orbital_table.columnWidth(column)
            self.orbital_table.setColumnWidth(column, width + _ORBITAL_COLUMN_PADDING)

    def sync_from_viewer(self) -> None:
        """Refresh fields from the viewer's instance configuration."""
        config = self._viewer.config
        self.grid_type.setCurrentText(config.grid.default_type)
        self.radius_points.setValue(config.grid.spherical.num_r_points)
        self.theta_points.setValue(config.grid.spherical.num_theta_points)
        self.phi_points.setValue(config.grid.spherical.num_phi_points)
        cartesian_config = config.grid.adaptive if config.grid.default_type == 'adaptive' else config.grid.cartesian
        for axis, field_name in zip('xyz', ('num_x_points', 'num_y_points', 'num_z_points'), strict=True):
            self.cartesian_fields[axis][2].setValue(getattr(cartesian_config, field_name))
        self.adaptive_scale.setText(format_scale(config.grid.adaptive.scale))
        self.contour.setValue(config.mo.contour)
        self.mo_opacity.setValue(config.mo.opacity)
        if config.mo.custom_colors:
            self.color_scheme.setCurrentText('custom')
            self.negative_color.setText(config.mo.custom_colors[0])
            self.positive_color.setText(config.mo.custom_colors[1])
        else:
            if self.color_scheme.findText(config.mo.color_scheme) < 0:
                self.color_scheme.insertItem(0, config.mo.color_scheme)
            self.color_scheme.setCurrentText(config.mo.color_scheme)
        self.molecule_opacity.setValue(config.molecule.opacity)
        self.show_atoms.setChecked(config.molecule.atom.show)
        self.show_bonds.setChecked(config.molecule.bond.show)
        self.show_axes.setChecked(config.show_axes)
        self.bond_max_length.setValue(config.molecule.bond.max_length)
        self.bond_radius.setValue(config.molecule.bond.radius)
        self.bond_color_type.setCurrentText(config.molecule.bond.color_type)
        self.bond_color.setText(config.molecule.bond.color)
        self.set_background_color(config.background_color)
        self._update_custom_color_visibility(self.color_scheme.currentText())
        self._update_bond_color_visibility(self.bond_color_type.currentText())

    def set_molecule_only(self, only_molecule: bool) -> None:
        """Hide controls that require molecular orbitals."""
        self.tabs.setTabVisible(self.tabs.indexOf(self.orbitals_tab), not only_molecule)
        self.tabs.setTabVisible(self.tabs.indexOf(self.grid_tab), not only_molecule)
        self.orbital_group.setVisible(not only_molecule)
        self.data_format.setEnabled(not only_molecule)
        self.data_scope.setEnabled(not only_molecule)
        self.data_export_button.setEnabled(not only_molecule)
        if not only_molecule:
            self._update_export_scope(self.data_format.currentText())

    def _on_orbital_selected(self) -> None:
        rows = self.orbital_table.selectionModel().selectedRows()
        if rows and self._viewer.gtos_ready:
            self._viewer.show_orbital(rows[0].row())

    def _next_orbital(self) -> None:
        if self.current_mo_ind + 1 < self.orbital_table.rowCount():
            self._viewer.show_orbital(self.current_mo_ind + 1)

    def _previous_orbital(self) -> None:
        if self.current_mo_ind > 0:
            self._viewer.show_orbital(self.current_mo_ind - 1)

    def highlight_orbital(self, index: int) -> None:
        self.orbital_table.blockSignals(True)
        if index < 0:
            self.orbital_table.clearSelection()
        else:
            self.orbital_table.selectRow(index)
        self.orbital_table.blockSignals(False)

    def set_loading_state(self, loading: bool, message: str = 'Tabulating orbitals…') -> None:
        self.loading_label.setText(message)
        self.loading_label.setVisible(loading)
        self.orbital_table.setEnabled(not loading)
        self.update_nav_button_states()

    def on_gtos_ready(self) -> None:
        self.set_loading_state(False)

    def update_nav_button_states(self) -> None:
        ready = self._viewer.gtos_ready if hasattr(self, '_viewer') else False
        self.previous_button.setEnabled(ready and self.current_mo_ind > 0)
        self.next_button.setEnabled(ready and self.current_mo_ind + 1 < self.orbital_table.rowCount())
        self.clear_button.setEnabled(ready and self.current_mo_ind >= 0)

    def _update_grid_field_visibility(self, grid_type: str) -> None:
        spherical = grid_type == 'spherical'
        for label, widget in self._spherical_rows:
            label.setVisible(spherical)
            widget.setVisible(spherical)
        self.cartesian_grid.setVisible(not spherical)
        self.adaptive_scale_label.setVisible(grid_type == 'adaptive')
        self.adaptive_scale.setVisible(grid_type == 'adaptive')
        if grid_type in {'cartesian', 'adaptive'}:
            grid_config = getattr(self._viewer.config.grid, grid_type)
            for axis, field_name in zip('xyz', ('num_x_points', 'num_y_points', 'num_z_points'), strict=True):
                self.cartesian_fields[axis][2].setValue(getattr(grid_config, field_name))

    def _apply_grid(self) -> None:
        try:
            axes, grid_type = self._grid_values()
            mode = self.grid_type.currentText()
            scale = parse_scale(self.adaptive_scale.text()) if mode == 'adaptive' else None
            self._viewer.update_grid(axes, grid_type, mode=mode, adaptive_scale=scale)
        except (RuntimeError, ValueError) as exc:
            self._viewer.report_error('Grid update failed', exc)

    def _grid_values(
        self,
    ) -> tuple[tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]], GridType]:
        if self.grid_type.currentText() == 'spherical':
            return (
                (
                    np.linspace(0, self.radius.value(), self.radius_points.value()),
                    np.linspace(0, np.pi, self.theta_points.value()),
                    np.linspace(0, 2 * np.pi, self.phi_points.value()),
                ),
                GridType.SPHERICAL,
            )
        values = [self.cartesian_fields[axis] for axis in 'xyz']
        if any(minimum.value() >= maximum.value() for minimum, maximum, _points in values):
            raise ValueError('Each cartesian minimum must be smaller than its maximum.')
        axes = tuple(
            np.linspace(minimum.value(), maximum.value(), points.value()) for minimum, maximum, points in values
        )
        if self.grid_type.currentText() == 'adaptive' and any(
            points.value() < _MIN_ADAPTIVE_POINTS for _min, _max, points in values
        ):
            raise ValueError('Adaptive grid point counts must be at least 2.')
        return (axes[0], axes[1], axes[2]), GridType.CARTESIAN

    def _apply_appearance(self) -> None:
        try:
            self._viewer.apply_appearance(**self._appearance_values())
        except (RuntimeError, ValueError) as exc:
            self._viewer.report_error('Appearance update failed', exc)

    def _appearance_values(self) -> dict[str, Any]:
        custom_colors = [self.negative_color.text(), self.positive_color.text()]
        if self.color_scheme.currentText() == 'custom' and not all(map(_is_color_like, custom_colors)):
            raise ValueError('Both custom molecular-orbital colors must be valid colors.')
        background_color = self.selected_background_color
        if not _is_color_like(background_color):
            raise ValueError('Background color must be a valid color.')
        if self.bond_color_type.currentText() == 'uniform' and not _is_color_like(self.bond_color.text()):
            raise ValueError('Uniform bond color must be a valid color.')
        return {
            'contour': self.contour.value(),
            'mo_opacity': self.mo_opacity.value(),
            'color_scheme': self.color_scheme.currentText(),
            'custom_colors': custom_colors,
            'molecule_opacity': self.molecule_opacity.value(),
            'show_atoms': self.show_atoms.isChecked(),
            'show_bonds': self.show_bonds.isChecked(),
            'bond_max_length': self.bond_max_length.value(),
            'bond_radius': self.bond_radius.value(),
            'bond_color_type': self.bond_color_type.currentText(),
            'bond_color': self.bond_color.text(),
            'background_color': background_color,
        }

    def _save_settings(self) -> None:
        try:
            self._viewer.apply_appearance(**self._appearance_values())
            self._viewer.save_settings()
        except (OSError, RuntimeError, ValueError) as exc:
            self._viewer.report_error('Saving settings failed', exc)

    def _update_custom_color_visibility(self, scheme: str) -> None:
        custom = scheme == 'custom'
        if not custom:
            colors = import_module('matplotlib.colors')
            colormap = import_module('matplotlib').colormaps[scheme]
            self.negative_color.setText(colors.to_hex(colormap(0.0)))
            self.positive_color.setText(colors.to_hex(colormap(1.0)))
        self.negative_color.setEnabled(custom)
        self.positive_color.setEnabled(custom)

    def _update_bond_color_visibility(self, color_type: str) -> None:
        uniform = color_type == 'uniform'
        self.bond_color_label.setVisible(uniform)
        self.bond_color.setVisible(uniform)

    def _update_background_color_visibility(self) -> None:
        custom = self.background_color_choice.currentData() is None
        self.background_color_label.setVisible(custom)
        self.background_color.setVisible(custom)

    def set_background_color(self, color: str) -> None:
        """Synchronize the background preset and custom color controls."""
        index = self.background_color_choice.findData(color)
        if index < 0:
            index = self.background_color_choice.count() - 1
            self.background_color.setText(color)
        self.background_color_choice.setCurrentIndex(index)
        self._update_background_color_visibility()

    @property
    def selected_background_color(self) -> str:
        """Background color selected by the controls."""
        selected = self.background_color_choice.currentData()
        return self.background_color.text() if selected is None else str(selected)

    def _update_export_scope(self, file_format: str) -> None:
        enabled = file_format != 'cube'
        if not enabled and self.data_scope.currentData() == 'all':
            self.data_scope.setCurrentIndex(self.data_scope.findData('current'))
        self.data_scope.setEnabled(enabled)

    def _update_transparency_option(self, file_format: str) -> None:
        enabled = file_format == 'png'
        if not enabled:
            self.transparent_background.setChecked(False)
        self.transparent_background.setEnabled(enabled)

    def _request_data_export(self) -> None:
        scope = str(self.data_scope.currentData())
        if scope == 'current' and self.current_mo_ind < 0:
            self._viewer.report_error('Export failed', ValueError('No orbital is currently selected.'))
            return
        self._warn_if_export_unhandled()
        self._viewer.export_requested.emit(
            'data',
            {'format': self.data_format.currentText(), 'scope': scope},
        )

    def _request_image_export(self) -> None:
        self._warn_if_export_unhandled()
        self._viewer.export_requested.emit(
            'image',
            {
                'format': self.image_format.currentText(),
                'transparent': self.transparent_background.isChecked(),
            },
        )

    def _warn_if_export_unhandled(self) -> None:
        if not self._viewer.has_export_handler:
            logger.warning(
                'Export button clicked without an export_requested receiver; '
                'connect the signal or call export_data/export_image directly.',
            )


class OrbitalViewer(QWidget, _PlotterRendering):
    """Embeddable Qt widget for molecules and molecular orbitals.

    The host must create and run :class:`QApplication`. Constructing this
    widget never starts an event loop or shows a top-level window.

    Parameters
    ----------
    filename : str | os.PathLike[str] | None, optional
        Path to the Molden file.
    content : str | None, optional
        Complete contents of a Molden file.
    only_molecule : bool, optional
        Skip molecular-orbital parsing and controls.
    tabulator : Tabulator, optional
        Existing structured-grid tabulator with cached GTO values.
    config : Config | MainConfig | Mapping, optional
        Per-viewer configuration overrides.
    show_controls : bool, optional
        Show the built-in control panel. Disable it when the host supplies
        its own dashboard.
    parent : QWidget, optional
        Parent widget supplied by the host application.
    """

    loading_changed = Signal(bool)
    input_ready = Signal()
    orbital_changed = Signal(int)
    error_occurred = Signal(str, object)
    export_requested = Signal(str, object)

    def __init__(
        self,
        *,
        filename: str | PathLike[str] | None = None,
        content: str | None = None,
        only_molecule: bool = False,
        tabulator: Tabulator | None = None,
        config: Config | MainConfig | Mapping[str, Any] | None = None,
        show_controls: bool = True,
        parent: QWidget | None = None,
    ) -> None:
        application = QApplication.instance()
        if application is None or not application.inherits('QApplication'):
            raise RuntimeError('OrbitalViewer requires an existing QApplication.')
        super().__init__(parent)
        config_class = import_module('moldenViz._config_module').Config
        if isinstance(config, config_class):
            current_config = cast(Any, config)
            self._config = config_class(
                current_config._pydantic_config.model_dump(by_alias=True),  # ruff:ignore[private-member-access]
            )
        else:
            self._config = config_class(config)

        self._on_screen = True
        self._closed = False
        self._screen_watch = False
        self._only_molecule = only_molecule
        self._gtos_ready = only_molecule
        self._grid_mode = self._config.grid.default_type
        self._adaptive_scale: AdaptiveScale = self._config.grid.adaptive.scale
        self._adaptive_ready = False
        self._adaptive_mesh = import_module('pyvista').UnstructuredGrid()
        self._adaptive_gtos: NDArray[np.floating] | None = None
        self._current_orbital_index = -1
        self._controls_visible = show_controls
        self._selection_screen: OrbitalControlPanel | None = None
        self._orb_actor: pv.Actor | None = None
        self._molecule_actors: list[Any] = []
        self._atom_actors: list[Any] = []
        self._bond_actors: list[Any] = []
        self._dispatcher = _CompletionDispatcher(self)
        self._gto_job: BackgroundJob[_GTOResult] = BackgroundJob(_GTO_EXECUTOR, self._dispatcher.dispatch)
        self._adaptive_job: BackgroundJob[_AdaptiveResult] = BackgroundJob(_GTO_EXECUTOR, self._dispatcher.dispatch)

        interactor_class = _load_qt_interactor()
        self.interactor = interactor_class(self, auto_update=5.0)
        self._pv_plotter = self.interactor
        self.interactor.set_background(self._config.background_color)
        self.set_axes_visible(self._config.show_axes)
        self.controls = OrbitalControlPanel(self)
        self._selection_screen = self.controls
        splitter = QSplitter(Qt.Orientation.Horizontal, self)
        splitter.addWidget(self.controls)
        splitter.addWidget(self.interactor)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        self._splitter = splitter
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(splitter)
        self.controls.sync_from_viewer()
        self.set_controls_visible(show_controls)

        if filename is not None or content is not None or tabulator is not None:
            self.set_input(
                filename=filename,
                content=content,
                tabulator=tabulator,
                only_molecule=only_molecule,
            )

    @property
    def config(self) -> Config:
        """Per-instance visualization configuration."""
        return self._config

    @property
    def gtos_ready(self) -> bool:
        """Whether orbital data can be rendered."""
        return self._gtos_ready and (self._grid_mode != 'adaptive' or self._adaptive_ready)

    @property
    def axes_visible(self) -> bool:
        """Whether the orientation axes are visible."""
        return bool(self._config.show_axes)

    def set_axes_visible(self, visible: bool) -> None:
        """Show or hide the orientation axes.

        Parameters
        ----------
        visible : bool
            Whether the orientation axes should be visible.
        """
        self._config.config.show_axes = visible
        if visible:
            self.interactor.show_axes()
        else:
            self.interactor.hide_axes()
        if hasattr(self, 'controls'):
            self.controls.show_axes.setChecked(visible)

    def showEvent(self, event: QShowEvent) -> None:  # ruff: ignore[invalid-function-name]
        """Watch for screen changes once the native window exists."""
        super().showEvent(event)
        if self._screen_watch:
            return
        handle = self.window().windowHandle()
        if handle is not None:
            handle.installEventFilter(self)
            self._screen_watch = True

    def eventFilter(self, watched: QObject, event: QEvent) -> bool:  # ruff: ignore[invalid-function-name]
        """Repaint native widgets after a device-pixel-ratio change.

        Returns
        -------
        bool
            Whether the event should stop propagating.
        """
        if event.type() == QEvent.Type.DevicePixelRatioChange:
            QTimer.singleShot(0, self._repaint_screen_widgets)
        return super().eventFilter(watched, event)

    def _repaint_screen_widgets(self) -> None:
        """Repaint widgets promoted to native windows."""
        self.repaint()
        self.controls.repaint()
        for index in range(1, self._splitter.count()):
            self._splitter.handle(index).repaint()

    @property
    def controls_visible(self) -> bool:
        """Whether the built-in control panel is enabled for display."""
        return self._controls_visible

    @property
    def current_orbital_index(self) -> int:
        """Currently displayed orbital index, or ``-1`` when cleared."""
        return self._current_orbital_index

    @property
    def has_export_handler(self) -> bool:
        """Whether a host has connected a receiver to ``export_requested``."""
        signal_index = self.metaObject().indexOfSignal('export_requested(QString,PyObject)')
        return signal_index >= 0 and self.isSignalConnected(self.metaObject().method(signal_index))

    @property
    def molecular_orbitals(self) -> tuple[MolecularOrbital, ...]:
        """Molecular orbitals available to host-provided controllers."""
        if not hasattr(self, 'tabulator'):
            return ()
        return tuple(self.tabulator.molecular_orbitals)

    @property
    def _gto_future(self) -> Future[_GTOResult] | None:
        """Compatibility view of the pending background future."""
        return self._gto_job.future

    def set_input(
        self,
        *,
        filename: str | PathLike[str] | None = None,
        content: str | None = None,
        tabulator: Tabulator | None = None,
        only_molecule: bool | None = None,
    ) -> None:
        """Load or replace the Molden input displayed by this widget."""
        if self._closed:
            raise RuntimeError('Cannot load input into a closed OrbitalViewer.')
        if tabulator is not None:
            if filename is not None or content is not None:
                raise ValueError('filename and content must not be provided with tabulator.')
        else:
            _validate_molden_input(filename, content)
        self._cancel_gto_future()
        self._clear_scene()
        self._invalidate_adaptive_grid()
        self._current_orbital_index = -1
        if only_molecule is not None:
            self._only_molecule = only_molecule

        if tabulator is not None:
            if not hasattr(tabulator, 'grid'):
                raise ValueError('Tabulator does not have grid attribute.')
            if tabulator.grid_type == GridType.UNKNOWN:
                raise ValueError('The viewer only supports spherical and cartesian grids.')
            if not tabulator.has_gtos and not self._only_molecule:
                raise ValueError('Tabulator does not have tabulated GTOs.')
            self.tabulator = tabulator
            self._grid_mode = tabulator.grid_type.value
        else:
            self.tabulator = Tabulator(
                filename=filename,
                content=content,
                only_molecule=self._only_molecule,
            )
            self._grid_mode = self._config.grid.default_type
        self._adaptive_scale = self._config.grid.adaptive.scale

        self._gtos_ready = self._only_molecule or self.tabulator.has_gtos
        self._molecule_opacity = self._config.molecule.opacity
        self._load_molecule(self._config)

        if not self._only_molecule and tabulator is None:
            self._create_default_grid()
            self._gtos_ready = False

        if not self._only_molecule:
            self._orb_mesh = self._create_mo_mesh()
            self._contour = self._config.mo.contour
            self._opacity = self._config.mo.opacity
            self._cmap = (
                self._custom_cmap_from_colors(self._config.mo.custom_colors)
                if self._config.mo.custom_colors
                else self._config.mo.color_scheme
            )
            self.controls.set_molecular_orbitals(self.tabulator.molecular_orbitals)
        else:
            self.controls.set_molecular_orbitals([])

        self.controls.set_molecule_only(self._only_molecule)
        self.controls.set_loading_state(not self._gtos_ready)
        if not self._gtos_ready:
            self._schedule_gto_tabulation()
        self.input_ready.emit()

    def _create_default_grid(self) -> None:
        radius = max(
            self._config.grid.max_radius_multiplier * self._molecule.max_radius,
            self._config.grid.min_radius,
        )
        self.controls.radius.setValue(float(radius))
        if self._config.grid.default_type == 'spherical':
            spherical = self._config.grid.spherical
            self.tabulator.spherical_grid(
                np.linspace(0, radius, spherical.num_r_points),
                np.linspace(0, np.pi, spherical.num_theta_points),
                np.linspace(0, 2 * np.pi, spherical.num_phi_points),
                tabulate_gtos=False,
            )
        else:
            cartesian = (
                self._config.grid.adaptive
                if self._config.grid.default_type == 'adaptive'
                else self._config.grid.cartesian
            )
            axes = []
            for axis, points in zip(
                'xyz',
                (cartesian.num_x_points, cartesian.num_y_points, cartesian.num_z_points),
                strict=True,
            ):
                minimum, maximum, point_widget = self.controls.cartesian_fields[axis]
                minimum.setValue(-radius)
                maximum.setValue(radius)
                point_widget.setValue(points)
                axes.append(np.linspace(-radius, radius, points))
            self.tabulator.cartesian_grid(*axes, tabulate_gtos=False)

    def show_orbital(self, index: int) -> None:
        """Show one orbital, or clear it when ``index`` is ``-1``."""
        if index < -1 or index >= len(self.molecular_orbitals):
            raise IndexError(f'Orbital index out of range: {index}')
        self.plot_orbital(index)
        self._current_orbital_index = index
        self.controls.current_mo_ind = index
        self.controls.highlight_orbital(index)
        self.controls.update_nav_button_states()
        self.orbital_changed.emit(index)

    def set_controls_visible(self, visible: bool) -> None:
        """Show or hide the built-in control panel."""
        self._controls_visible = visible
        self.controls.setVisible(visible)

    def update_grid(
        self,
        axes: tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]],
        grid_type: GridType,
        *,
        mode: str | None = None,
        adaptive_scale: AdaptiveScale | None = None,
    ) -> None:
        """Replace the structured grid and schedule fresh GTO tabulation."""
        self._ensure_grid_update_supported()
        resolved_mode = grid_type.value if mode is None else mode
        if resolved_mode not in {'spherical', 'cartesian', 'adaptive'}:
            raise ValueError(f'Unknown grid mode: {resolved_mode}')
        if resolved_mode == 'adaptive' and grid_type != GridType.CARTESIAN:
            raise ValueError('Adaptive mode requires Cartesian grid axes.')
        self._grid_mode = resolved_mode
        if adaptive_scale is not None:
            normalize_scale(adaptive_scale)
            self._adaptive_scale = adaptive_scale
        self._config.config.grid.default_type = resolved_mode
        if resolved_mode == 'adaptive':
            adaptive = self._config.config.grid.adaptive
            adaptive.num_x_points, adaptive.num_y_points, adaptive.num_z_points = map(len, axes)
            adaptive.scale = self._adaptive_scale
        self._update_mesh(*axes, grid_type)

    def _ensure_grid_update_supported(self) -> None:
        if self._only_molecule:
            raise RuntimeError('Molecule-only viewers do not have an orbital grid.')
        if not hasattr(self, 'tabulator'):
            raise RuntimeError('Load input before updating the orbital grid.')

    def set_spherical_grid(
        self,
        *,
        radius: float,
        radial_points: int,
        theta_points: int,
        phi_points: int,
    ) -> None:
        """Build and apply a spherical grid from dashboard-friendly values."""
        if radius <= 0:
            raise ValueError('radius must be greater than zero.')
        self._validate_grid_points(
            radial_points=radial_points,
            theta_points=theta_points,
            phi_points=phi_points,
        )
        self._ensure_grid_update_supported()
        self.controls.grid_type.setCurrentText('spherical')
        self.controls.radius.setValue(radius)
        self.controls.radius_points.setValue(radial_points)
        self.controls.theta_points.setValue(theta_points)
        self.controls.phi_points.setValue(phi_points)
        self._config.config.grid.default_type = 'spherical'
        self._config.config.grid.spherical.num_r_points = radial_points
        self._config.config.grid.spherical.num_theta_points = theta_points
        self._config.config.grid.spherical.num_phi_points = phi_points
        self.update_grid(
            (
                np.linspace(0, radius, radial_points),
                np.linspace(0, np.pi, theta_points),
                np.linspace(0, 2 * np.pi, phi_points),
            ),
            GridType.SPHERICAL,
        )

    def set_cartesian_grid(
        self,
        *,
        x_range: tuple[float, float],
        y_range: tuple[float, float],
        z_range: tuple[float, float],
        x_points: int,
        y_points: int,
        z_points: int,
    ) -> None:
        """Build and apply a Cartesian grid from bounds and point counts."""
        self._validate_grid_points(x_points=x_points, y_points=y_points, z_points=z_points)
        ranges = {'x': x_range, 'y': y_range, 'z': z_range}
        for axis, bounds in ranges.items():
            if bounds[0] >= bounds[1]:
                raise ValueError(f'{axis}_range minimum must be smaller than its maximum.')
        self._ensure_grid_update_supported()

        self.controls.grid_type.setCurrentText('cartesian')
        for axis, bounds, points in zip(
            'xyz',
            (x_range, y_range, z_range),
            (x_points, y_points, z_points),
            strict=True,
        ):
            minimum, maximum, point_widget = self.controls.cartesian_fields[axis]
            minimum.setValue(bounds[0])
            maximum.setValue(bounds[1])
            point_widget.setValue(points)
        self._config.config.grid.default_type = 'cartesian'
        self._config.config.grid.cartesian.num_x_points = x_points
        self._config.config.grid.cartesian.num_y_points = y_points
        self._config.config.grid.cartesian.num_z_points = z_points
        self.update_grid(
            (
                np.linspace(*x_range, x_points),
                np.linspace(*y_range, y_points),
                np.linspace(*z_range, z_points),
            ),
            GridType.CARTESIAN,
        )

    def set_adaptive_grid(
        self,
        *,
        x_range: tuple[float, float],
        y_range: tuple[float, float],
        z_range: tuple[float, float],
        x_points: int,
        y_points: int,
        z_points: int,
        scale: AdaptiveScale,
    ) -> None:
        """Build and apply an adaptive Cartesian grid from explicit values."""
        self._validate_grid_points(x_points=x_points, y_points=y_points, z_points=z_points)
        if min(x_points, y_points, z_points) < _MIN_ADAPTIVE_POINTS:
            raise ValueError('Adaptive grid point counts must be at least 2.')
        ranges = {'x': x_range, 'y': y_range, 'z': z_range}
        for axis, bounds in ranges.items():
            if bounds[0] >= bounds[1]:
                raise ValueError(f'{axis}_range minimum must be smaller than its maximum.')
        normalize_scale(scale)
        self._ensure_grid_update_supported()

        self.controls.grid_type.setCurrentText('adaptive')
        self.controls.adaptive_scale.setText(format_scale(scale))
        for axis, bounds, points in zip(
            'xyz',
            (x_range, y_range, z_range),
            (x_points, y_points, z_points),
            strict=True,
        ):
            minimum, maximum, point_widget = self.controls.cartesian_fields[axis]
            minimum.setValue(bounds[0])
            maximum.setValue(bounds[1])
            point_widget.setValue(points)
        self.update_grid(
            (
                np.linspace(*x_range, x_points),
                np.linspace(*y_range, y_points),
                np.linspace(*z_range, z_points),
            ),
            GridType.CARTESIAN,
            mode='adaptive',
            adaptive_scale=scale,
        )

    @staticmethod
    def _validate_grid_points(**point_counts: int) -> None:
        for name, count in point_counts.items():
            if count < 1:
                raise ValueError(f'{name} must be at least 1.')

    def apply_appearance(
        self,
        *,
        contour: float,
        mo_opacity: float,
        color_scheme: str,
        custom_colors: list[str] | None,
        molecule_opacity: float,
        show_atoms: bool,
        show_bonds: bool,
        bond_max_length: float,
        bond_radius: float,
        bond_color_type: str,
        bond_color: str,
        background_color: str,
    ) -> None:
        """Apply validated appearance settings to this viewer instance."""
        self._validate_appearance(
            contour=contour,
            mo_opacity=mo_opacity,
            color_scheme=color_scheme,
            custom_colors=custom_colors,
            molecule_opacity=molecule_opacity,
            bond_max_length=bond_max_length,
            bond_radius=bond_radius,
            bond_color_type=bond_color_type,
            bond_color=bond_color,
            background_color=background_color,
        )
        contour_changed = hasattr(self, '_contour') and contour != self._contour
        self._config.config.mo.contour = contour
        self._config.config.mo.opacity = mo_opacity
        self._config.config.mo.color_scheme = color_scheme if color_scheme != 'custom' else 'bwr'
        self._config.config.mo.custom_colors = custom_colors if color_scheme == 'custom' else None
        self._config.config.molecule.opacity = molecule_opacity
        self._config.config.molecule.atom.show = show_atoms
        self._config.config.molecule.bond.show = show_bonds
        self._config.config.molecule.bond.max_length = bond_max_length
        self._config.config.molecule.bond.radius = bond_radius
        self._config.config.molecule.bond.color_type = bond_color_type
        self._config.config.molecule.bond.color = bond_color
        self._config.config.background_color = background_color
        self._contour = contour
        self._opacity = mo_opacity
        self._molecule_opacity = molecule_opacity
        self._cmap = (
            self._custom_cmap_from_colors(custom_colors) if custom_colors and color_scheme == 'custom' else color_scheme
        )
        self.interactor.set_background(background_color)
        if hasattr(self, 'tabulator'):
            self._load_molecule(self._config)
            if contour_changed and self._grid_mode == 'adaptive':
                self._invalidate_adaptive_grid(rebuild=True)
                self.interactor.update()
                return
            if not self._only_molecule and self._current_orbital_index >= 0 and self._gtos_ready:
                self.plot_orbital(self._current_orbital_index)
        self.interactor.update()

    def update_appearance(
        self,
        *,
        contour: float | None = None,
        mo_opacity: float | None = None,
        color_scheme: str | None = None,
        custom_colors: list[str] | None = None,
        molecule_opacity: float | None = None,
        show_atoms: bool | None = None,
        show_bonds: bool | None = None,
        bond_max_length: float | None = None,
        bond_radius: float | None = None,
        bond_color_type: str | None = None,
        bond_color: str | None = None,
        background_color: str | None = None,
    ) -> None:
        """Apply a partial appearance update for host-provided controllers."""
        config = self._config
        resolved_custom_colors = custom_colors if custom_colors is not None else config.mo.custom_colors
        resolved_color_scheme = color_scheme
        if resolved_color_scheme is None:
            resolved_color_scheme = (
                'custom' if custom_colors is not None or config.mo.custom_colors else config.mo.color_scheme
            )
        self.apply_appearance(
            contour=config.mo.contour if contour is None else contour,
            mo_opacity=config.mo.opacity if mo_opacity is None else mo_opacity,
            color_scheme=resolved_color_scheme,
            custom_colors=resolved_custom_colors,
            molecule_opacity=config.molecule.opacity if molecule_opacity is None else molecule_opacity,
            show_atoms=config.molecule.atom.show if show_atoms is None else show_atoms,
            show_bonds=config.molecule.bond.show if show_bonds is None else show_bonds,
            bond_max_length=config.molecule.bond.max_length if bond_max_length is None else bond_max_length,
            bond_radius=config.molecule.bond.radius if bond_radius is None else bond_radius,
            bond_color_type=config.molecule.bond.color_type if bond_color_type is None else bond_color_type,
            bond_color=config.molecule.bond.color if bond_color is None else bond_color,
            background_color=config.background_color if background_color is None else background_color,
        )
        self.controls.sync_from_viewer()

    @staticmethod
    def _validate_appearance(
        *,
        contour: float,
        mo_opacity: float,
        color_scheme: str,
        custom_colors: list[str] | None,
        molecule_opacity: float,
        bond_max_length: float,
        bond_radius: float,
        bond_color_type: str,
        bond_color: str,
        background_color: str,
    ) -> None:
        if contour <= 0:
            raise ValueError('contour must be greater than zero.')
        if not 0 <= mo_opacity <= 1 or not 0 <= molecule_opacity <= 1:
            raise ValueError('opacity values must be between 0 and 1.')
        if color_scheme == 'custom':
            valid_custom_colors = (
                custom_colors is not None
                and len(custom_colors) == _CUSTOM_COLOR_COUNT
                and all(map(_is_color_like, custom_colors))
            )
            if not valid_custom_colors:
                raise ValueError('custom_colors must contain exactly two valid colors.')
        elif not _has_colormap(color_scheme):
            raise ValueError(f'Unknown color scheme: {color_scheme}')
        if bond_max_length <= 0 or bond_radius <= 0:
            raise ValueError('bond dimensions must be greater than zero.')
        if bond_color_type not in {'uniform', 'split'}:
            raise ValueError(f'Unknown bond color type: {bond_color_type}')
        if bond_color_type == 'uniform' and not _is_color_like(bond_color):
            raise ValueError(f'Invalid bond color: {bond_color}')
        if not _is_color_like(background_color):
            raise ValueError(f'Invalid background color: {background_color}')

    def set_background_color(self, color: str) -> None:
        """Set this viewer's render background color."""
        if not _is_color_like(color):
            raise ValueError(f'Invalid background color: {color}')
        self._config.config.background_color = color
        self.controls.set_background_color(color)
        self.interactor.set_background(color)

    def export_data(self, path: str | Path, *, file_format: str, scope: str = 'current') -> None:
        """Export orbital data to an explicit destination without dialogs."""
        if file_format not in {'vtk', 'cube'}:
            raise ValueError(f'Unsupported orbital export format: {file_format}')
        if scope not in {'current', 'all'}:
            raise ValueError(f'Unsupported orbital export scope: {scope}')
        if file_format == 'cube' and scope == 'all':
            raise ValueError('Cube format only supports one orbital.')
        index = self._current_orbital_index
        if scope == 'current' and index < 0:
            raise ValueError('No orbital is currently selected.')
        destination = Path(path)
        expected_suffix = f'.{file_format}'
        if destination.suffix.lower() != expected_suffix:
            destination = destination.with_suffix(expected_suffix)
        self.tabulator.export(destination, mo_index=index if scope == 'current' else None)

    def export_image(
        self,
        path: str | Path,
        *,
        file_format: str,
        transparent: bool = False,
    ) -> None:
        """Export the current render to an explicit destination without dialogs."""
        if file_format not in {'png', 'jpeg', 'svg', 'pdf'}:
            raise ValueError(f'Unsupported image format: {file_format}')
        suffix = '.jpg' if file_format == 'jpeg' else f'.{file_format}'
        destination = Path(path)
        if destination.suffix.lower() not in ({'.jpg', '.jpeg'} if file_format == 'jpeg' else {suffix}):
            destination = destination.with_suffix(suffix)
        if file_format in {'svg', 'pdf'}:
            self.interactor.save_graphic(destination)
        else:
            self.interactor.screenshot(
                destination,
                transparent_background=transparent if file_format == 'png' else False,
            )

    def save_settings(self) -> None:
        """Persist this viewer's current settings to the user TOML file."""
        self._config._save_current_config()  # ruff:ignore[private-member-access]

    def wait_for_gtos(self, timeout: float | None = None) -> None:
        """Block until background GTO tabulation finishes."""
        if self._gtos_ready:
            return
        if self._gto_job.future is None:
            raise RuntimeError('GTO tabulation has not been scheduled.')
        result = self._gto_job.wait(timeout=timeout)
        if not self._gtos_ready:
            self._apply_gtos_ready(result, 0.0)

    def wait_for_adaptive_grid(self, timeout: float | None = None) -> None:
        """Block until adaptive-grid tabulation finishes."""
        if self._grid_mode != 'adaptive' or self._adaptive_ready:
            return
        if self._adaptive_job.future is None:
            raise RuntimeError('Adaptive-grid tabulation has not been scheduled.')
        result = self._adaptive_job.wait(timeout=timeout)
        if not self._adaptive_ready:
            self._apply_adaptive_ready(result, 0.0)

    def _schedule_gto_tabulation(
        self,
        axes: tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]] | None = None,
        grid_type: GridType | None = None,
    ) -> None:
        if self._only_molecule or self._gtos_ready or self._gto_job.pending:
            return
        if axes is None:
            current_axes = self.tabulator.grid_axes
            if current_axes is None:
                raise RuntimeError('Structured grid axes are not available.')
            resolved_axes = current_axes
            resolved_grid_type = self.tabulator.grid_type
            current_grid = self.tabulator.grid.copy()
        elif grid_type is None:
            raise ValueError('Grid type is required when scheduling new axes.')
        else:
            resolved_axes = axes
            resolved_grid_type = grid_type
            current_grid = None

        frozen_axes = (resolved_axes[0].copy(), resolved_axes[1].copy(), resolved_axes[2].copy())
        for axis in frozen_axes:
            axis.setflags(write=False)
        tabulator = self.tabulator

        def build_and_tabulate() -> _GTOResult:
            grid = current_grid
            if grid is None:
                grid = Tabulator._build_grid(  # ruff:ignore[private-member-access]
                    frozen_axes[0],
                    frozen_axes[1],
                    frozen_axes[2],
                    resolved_grid_type,
                )
            grid.setflags(write=False)
            return _GTOResult(
                grid=grid,
                axes=frozen_axes,
                grid_type=resolved_grid_type,
                gtos=tabulator.compute_gtos(grid),
            )

        self.loading_changed.emit(True)
        self._gto_job.start(build_and_tabulate, on_success=self._apply_gtos_ready, on_error=self._handle_gto_error)

    def _apply_gtos_ready(self, result: _GTOResult, elapsed: float) -> None:
        if not self._on_screen:
            return
        self.tabulator._set_structured_grid(result.grid, result.axes, result.grid_type)  # ruff:ignore[private-member-access]
        self.tabulator.set_gtos(result.gtos)
        self._gtos_ready = True
        self._orb_mesh = self._create_mo_mesh()
        logger.info('GTO tabulation completed in %.2fs.', elapsed)
        if self._grid_mode == 'adaptive':
            self._schedule_adaptive_tabulation()
            return
        self.controls.on_gtos_ready()
        if self._current_orbital_index >= 0:
            self.plot_orbital(self._current_orbital_index)
        self.loading_changed.emit(False)

    def _schedule_adaptive_tabulation(self) -> None:
        """Build the union of all coarse MO contours and tabulate it once."""
        if self._only_molecule or self._grid_mode != 'adaptive' or not self._gtos_ready:
            return
        if self._adaptive_job.pending:
            return
        if self.tabulator.grid_type != GridType.CARTESIAN:
            raise RuntimeError('Adaptive grids require a Cartesian coarse grid.')

        coarse_grid = self._create_mo_mesh()
        coarse_mos = self.tabulator.tabulate_mos()
        contour = self._contour
        scale = normalize_scale(self._adaptive_scale)
        tabulator = self.tabulator

        def build_and_tabulate() -> _AdaptiveResult:
            cell_ids = crossed_cell_ids(coarse_grid, coarse_mos, contour)
            mesh = refined_grid(coarse_grid, cell_ids, scale)
            gtos = tabulator.compute_gtos(np.asarray(mesh.points))
            return _AdaptiveResult(mesh=mesh, gtos=gtos, crossed_cells=len(cell_ids))

        self._adaptive_ready = False
        self.controls.set_loading_state(True, 'Refining contour grid…')
        self.loading_changed.emit(True)
        logger.info('Starting adaptive-grid tabulation for all molecular orbitals...')
        self._adaptive_job.start(
            build_and_tabulate,
            on_success=self._apply_adaptive_ready,
            on_error=self._handle_adaptive_error,
        )

    def _apply_adaptive_ready(self, result: _AdaptiveResult, elapsed: float) -> None:
        """Install the reusable adaptive grid and GTO cache."""
        if not self._on_screen or self._grid_mode != 'adaptive':
            return
        self._adaptive_mesh = result.mesh
        self._adaptive_gtos = result.gtos
        self._adaptive_ready = True
        logger.info(
            'Adaptive-grid tabulation completed in %.2fs from %d crossed coarse cells (%d fine points).',
            elapsed,
            result.crossed_cells,
            result.mesh.n_points,
        )
        self.controls.on_gtos_ready()
        if self._current_orbital_index >= 0:
            self.plot_orbital(self._current_orbital_index)
        self.loading_changed.emit(False)

    def _handle_adaptive_error(self, exc: Exception) -> None:
        """Keep the previous actor and report failed adaptive tabulation."""
        self._adaptive_ready = False
        self.controls.on_gtos_ready()
        self.loading_changed.emit(False)
        self.report_error('Adaptive grid failed', exc)

    def _invalidate_adaptive_grid(self, *, rebuild: bool = False) -> None:
        """Discard adaptive cache and optionally rebuild it from the coarse grid."""
        self._adaptive_job.cancel()
        self._adaptive_ready = False
        self._adaptive_mesh = import_module('pyvista').UnstructuredGrid()
        self._adaptive_gtos = None
        if rebuild:
            self._schedule_adaptive_tabulation()

    def _handle_gto_error(self, exc: Exception) -> None:
        self._gtos_ready = self.tabulator.has_gtos
        self.controls.on_gtos_ready()
        self.loading_changed.emit(False)
        self.report_error('Orbital tabulation failed', exc)

    def _dispatch_gto_completion(self, callback: Callable[[], None]) -> None:
        if self._on_screen:
            self._dispatcher.dispatch(callback)

    def _ensure_gtos_ready(self) -> bool:
        return self._gtos_ready and (self._grid_mode != 'adaptive' or self._adaptive_ready)

    def _update_settings_button_states(self) -> None:
        """Mirror actor visibility in the public control widgets."""
        self.controls.show_atoms.setChecked(self.are_atoms_visible())
        self.controls.show_bonds.setChecked(self.are_bonds_visible())

    def _cancel_gto_future(self) -> None:
        self._gto_job.cancel()
        self._adaptive_job.cancel()

    def _clear_scene(self) -> None:
        for actor in [*self._molecule_actors, self._orb_actor]:
            if actor is not None:
                self.interactor.remove_actor(actor)
        self._molecule_actors = []
        self._atom_actors = []
        self._bond_actors = []
        self._orb_actor = None

    def report_error(self, title: str, exc: Exception) -> None:
        """Log an exception and notify the host without opening a dialog."""
        logger.error(title, exc_info=(type(exc), exc, exc.__traceback__))
        self.error_occurred.emit(title, exc)

    def _shutdown(self) -> None:
        if self._closed:
            return
        self._closed = True
        self._on_screen = False
        self._cancel_gto_future()
        self._dispatcher.callback_ready.disconnect()
        self.interactor.close()

    def close(self) -> bool:
        """Close the viewer and explicitly release its VTK render window.

        Returns
        -------
        bool
            Whether Qt accepted the close request.
        """
        self._shutdown()
        return super().close()

    def closeEvent(self, event: QCloseEvent) -> None:  # ruff: ignore[invalid-function-name]
        """Release rendering resources during parent-driven teardown."""
        self._shutdown()
        event.accept()
