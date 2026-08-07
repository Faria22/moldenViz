"""Qt-native embeddable molecular-orbital viewer."""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import matplotlib.colors as mcolors
import numpy as np
from PySide6.QtCore import QObject, Qt, Signal, Slot
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
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
from pyvistaqt import QtInteractor

from ._config_module import Config, MainConfig
from ._plotter_jobs import BackgroundJob
from ._plotter_rendering import _PlotterRendering
from .tabulator import GridType, Tabulator

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping
    from concurrent.futures import Future

    import pyvista as pv
    from numpy.typing import NDArray
    from PySide6.QtGui import QCloseEvent

    from .models import MolecularOrbital

logger = logging.getLogger(__name__)

__all__ = ['OrbitalViewer', 'ViewerConfig']

ViewerConfig = MainConfig
_GTO_EXECUTOR = ThreadPoolExecutor(max_workers=1)
_MO_COLOR_SCHEMES = ['bwr', 'RdBu', 'seismic', 'coolwarm', 'PiYG']
_ORBITAL_COLUMN_PADDING = 8


@dataclass(frozen=True)
class _GTOResult:
    """GTO values and the structured grid snapshot they describe."""

    grid: NDArray[np.floating]
    axes: tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]
    grid_type: GridType
    gtos: NDArray[np.floating]


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
        self._build_grid_tab()
        self._build_appearance_tab()
        self._build_export_tab()

    def _build_orbitals_tab(self) -> None:
        tab = QWidget(self)
        self.orbitals_tab = tab
        layout = QVBoxLayout(tab)
        self.loading_label = QLabel('', tab)
        self.loading_label.hide()
        layout.addWidget(self.loading_label)

        self.orbital_table = QTableWidget(0, 5, tab)
        self.orbital_table.setHorizontalHeaderLabels(['#', 'Sym', 'Spin', 'Occ', 'Energy'])
        self.orbital_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.orbital_table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.orbital_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.orbital_table.verticalHeader().hide()
        self.orbital_table.itemSelectionChanged.connect(self._on_orbital_selected)
        layout.addWidget(self.orbital_table)
        self._fit_orbital_columns()

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
        self.tabs.addTab(tab, 'Orbitals')
        self.update_nav_button_states()

    def _build_grid_tab(self) -> None:
        tab = QWidget(self)
        layout = QVBoxLayout(tab)
        form = QFormLayout()
        self.grid_type = QComboBox(tab)
        self.grid_type.addItems(['spherical', 'cartesian'])
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

        self.cartesian_fields: dict[str, tuple[QDoubleSpinBox, QDoubleSpinBox, QSpinBox]] = {}
        self._cartesian_rows: list[tuple[QLabel, QWidget]] = []
        for axis in 'xyz':
            row = QWidget(tab)
            row_layout = QHBoxLayout(row)
            row_layout.setContentsMargins(0, 0, 0, 0)
            minimum = self._double_spin(-1_000.0, 1_000.0, -5.0)
            maximum = self._double_spin(-1_000.0, 1_000.0, 5.0)
            points = self._int_spin(1, 10_000, 100)
            row_layout.addWidget(minimum)
            row_layout.addWidget(maximum)
            row_layout.addWidget(points)
            label = QLabel(f'{axis.upper()} min / max / points', tab)
            form.addRow(label, row)
            self.cartesian_fields[axis] = (minimum, maximum, points)
            self._cartesian_rows.append((label, row))

        layout.addLayout(form)
        apply_button = QPushButton('Apply grid', tab)
        apply_button.clicked.connect(self._apply_grid)
        layout.addWidget(apply_button)
        layout.addStretch()
        self.tabs.addTab(tab, 'Grid')
        self._update_grid_field_visibility(self.grid_type.currentText())

    def _build_appearance_tab(self) -> None:
        tab = QWidget(self)
        layout = QVBoxLayout(tab)

        self.orbital_group = QGroupBox('Molecular orbital', tab)
        orbital_form = QFormLayout(self.orbital_group)
        self.contour = self._double_spin(1e-6, 1e6, 0.1, decimals=6)
        self.mo_opacity = self._double_spin(0.0, 1.0, 1.0)
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
        self.molecule_opacity = self._double_spin(0.0, 1.0, 1.0)
        self.show_atoms = QCheckBox('Show atoms', tab)
        self.show_bonds = QCheckBox('Show bonds', tab)
        self.bond_max_length = self._double_spin(0.001, 1_000.0, 4.0)
        self.bond_radius = self._double_spin(0.001, 100.0, 0.15)
        self.bond_color_type = QComboBox(tab)
        self.bond_color_type.addItems(['uniform', 'split'])
        self.bond_color = QLineEdit('grey', tab)
        self.bond_color_type.currentTextChanged.connect(
            lambda color_type: self.bond_color.setEnabled(color_type == 'uniform'),
        )
        molecule_form.addRow('Opacity', self.molecule_opacity)
        molecule_form.addRow(self.show_atoms)
        molecule_form.addRow(self.show_bonds)
        molecule_form.addRow('Bond max length', self.bond_max_length)
        molecule_form.addRow('Bond radius', self.bond_radius)
        molecule_form.addRow('Bond colors', self.bond_color_type)
        molecule_form.addRow('Uniform color', self.bond_color)
        layout.addWidget(molecule_group)

        background_row = QHBoxLayout()
        self.background_color = QLineEdit('white', tab)
        background_row.addWidget(self.background_color)
        layout.addWidget(QLabel('Background color', tab))
        layout.addLayout(background_row)

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

    def _build_export_tab(self) -> None:
        tab = QWidget(self)
        layout = QFormLayout(tab)
        self.data_format = QComboBox(tab)
        self.data_format.addItems(['vtk', 'cube'])
        self.data_scope = QComboBox(tab)
        self.data_scope.addItems(['current', 'all'])
        self.data_format.currentTextChanged.connect(self._update_export_scope)
        self.data_export_button = QPushButton('Export orbital data…', tab)
        self.data_export_button.clicked.connect(self._request_data_export)
        self.image_format = QComboBox(tab)
        self.image_format.addItems(['png', 'jpeg', 'svg', 'pdf'])
        self.transparent_background = QCheckBox('Transparent PNG background', tab)
        self.image_format.currentTextChanged.connect(self._update_transparency_option)
        image_button = QPushButton('Export image…', tab)
        image_button.clicked.connect(self._request_image_export)
        layout.addRow('Data format', self.data_format)
        layout.addRow('Scope', self.data_scope)
        layout.addRow(self.data_export_button)
        layout.addRow('Image format', self.image_format)
        layout.addRow(self.transparent_background)
        layout.addRow(image_button)
        self.tabs.addTab(tab, 'Export')

    @staticmethod
    def _double_spin(minimum: float, maximum: float, value: float, *, decimals: int = 3) -> QDoubleSpinBox:
        widget = QDoubleSpinBox()
        widget.setRange(minimum, maximum)
        widget.setDecimals(decimals)
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
        for axis, field_name in zip('xyz', ('num_x_points', 'num_y_points', 'num_z_points'), strict=True):
            self.cartesian_fields[axis][2].setValue(getattr(config.grid.cartesian, field_name))
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
        self.bond_max_length.setValue(config.molecule.bond.max_length)
        self.bond_radius.setValue(config.molecule.bond.radius)
        self.bond_color_type.setCurrentText(config.molecule.bond.color_type)
        self.bond_color.setText(config.molecule.bond.color)
        self.background_color.setText(config.background_color)
        self._update_custom_color_visibility(self.color_scheme.currentText())

    def set_molecule_only(self, only_molecule: bool) -> None:
        """Hide controls that require molecular orbitals."""
        self.tabs.setTabVisible(self.tabs.indexOf(self.orbitals_tab), not only_molecule)
        self.tabs.setTabVisible(1, not only_molecule)
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
        for label, row in self._cartesian_rows:
            label.setVisible(not spherical)
            row.setVisible(not spherical)

    def _apply_grid(self) -> None:
        try:
            axes, grid_type = self._grid_values()
            self._viewer.update_grid(axes, grid_type)
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
        return (axes[0], axes[1], axes[2]), GridType.CARTESIAN

    def _apply_appearance(self) -> None:
        try:
            self._viewer.apply_appearance(**self._appearance_values())
        except (RuntimeError, ValueError) as exc:
            self._viewer.report_error('Appearance update failed', exc)

    def _appearance_values(self) -> dict[str, Any]:
        custom_colors = [self.negative_color.text(), self.positive_color.text()]
        if self.color_scheme.currentText() == 'custom' and not all(map(mcolors.is_color_like, custom_colors)):
            raise ValueError('Both custom molecular-orbital colors must be valid colors.')
        if not mcolors.is_color_like(self.background_color.text()):
            raise ValueError('Background color must be a valid color.')
        if self.bond_color_type.currentText() == 'uniform' and not mcolors.is_color_like(self.bond_color.text()):
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
            'background_color': self.background_color.text(),
        }

    def _save_settings(self) -> None:
        try:
            self._viewer.apply_appearance(**self._appearance_values())
            self._viewer.save_settings()
        except (OSError, RuntimeError, ValueError) as exc:
            self._viewer.report_error('Saving settings failed', exc)

    def _update_custom_color_visibility(self, scheme: str) -> None:
        custom = scheme == 'custom'
        self.negative_color.setEnabled(custom)
        self.positive_color.setEnabled(custom)

    def _update_export_scope(self, file_format: str) -> None:
        enabled = file_format != 'cube'
        if not enabled and self.data_scope.currentText() == 'all':
            self.data_scope.setCurrentText('current')
        self.data_scope.setEnabled(enabled)

    def _update_transparency_option(self, file_format: str) -> None:
        enabled = file_format == 'png'
        if not enabled:
            self.transparent_background.setChecked(False)
        self.transparent_background.setEnabled(enabled)

    def _request_data_export(self) -> None:
        if self.data_scope.currentText() == 'current' and self.current_mo_ind < 0:
            self._viewer.report_error('Export failed', ValueError('No orbital is currently selected.'))
            return
        self._viewer.export_requested.emit(
            'data',
            {'format': self.data_format.currentText(), 'scope': self.data_scope.currentText()},
        )

    def _request_image_export(self) -> None:
        self._viewer.export_requested.emit(
            'image',
            {
                'format': self.image_format.currentText(),
                'transparent': self.transparent_background.isChecked(),
            },
        )


class OrbitalViewer(QWidget, _PlotterRendering):
    """Embeddable Qt widget for molecules and molecular orbitals.

    The host must create and run :class:`QApplication`. Constructing this
    widget never starts an event loop or shows a top-level window.
    """

    loading_changed = Signal(bool)
    source_ready = Signal()
    orbital_changed = Signal(int)
    error_occurred = Signal(str, object)
    export_requested = Signal(str, object)

    def __init__(
        self,
        source: str | list[str] | None = None,
        *,
        only_molecule: bool = False,
        tabulator: Tabulator | None = None,
        config: Config | MainConfig | Mapping[str, Any] | None = None,
        parent: QWidget | None = None,
    ) -> None:
        application = QApplication.instance()
        if application is None or not application.inherits('QApplication'):
            raise RuntimeError('OrbitalViewer requires an existing QApplication.')
        super().__init__(parent)
        if isinstance(config, Config):
            self._config = Config(config._pydantic_config.model_dump(by_alias=True))  # ruff:ignore[private-member-access]
        else:
            self._config = Config(config)

        self._on_screen = True
        self._closed = False
        self._only_molecule = only_molecule
        self._gtos_ready = only_molecule
        self._selection_screen: OrbitalControlPanel | None = None
        self._orb_actor: pv.Actor | None = None
        self._molecule_actors: list[Any] = []
        self._atom_actors: list[Any] = []
        self._bond_actors: list[Any] = []
        self._dispatcher = _CompletionDispatcher(self)
        self._gto_job: BackgroundJob[_GTOResult] = BackgroundJob(_GTO_EXECUTOR, self._dispatcher.dispatch)

        self.interactor = QtInteractor(self, auto_update=5.0)
        self._pv_plotter = self.interactor
        self.interactor.set_background(self._config.background_color)
        self.interactor.show_axes()
        self.controls = OrbitalControlPanel(self)
        self._selection_screen = self.controls
        splitter = QSplitter(Qt.Orientation.Horizontal, self)
        splitter.addWidget(self.controls)
        splitter.addWidget(self.interactor)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(splitter)
        self.controls.sync_from_viewer()

        if source is not None or tabulator is not None:
            self.set_source(source, tabulator=tabulator, only_molecule=only_molecule)

    @property
    def config(self) -> Config:
        """Per-instance visualization configuration."""
        return self._config

    @property
    def gtos_ready(self) -> bool:
        """Whether orbital data can be rendered."""
        return self._gtos_ready

    @property
    def _gto_future(self) -> Future[_GTOResult] | None:
        """Compatibility view of the pending background future."""
        return self._gto_job.future

    def set_source(
        self,
        source: str | list[str] | None,
        *,
        tabulator: Tabulator | None = None,
        only_molecule: bool | None = None,
    ) -> None:
        """Load or replace the Molden source displayed by this widget."""
        if self._closed:
            raise RuntimeError('Cannot load a source into a closed OrbitalViewer.')
        if source is None and tabulator is None:
            raise ValueError('source is required when tabulator is not provided.')
        self._cancel_gto_future()
        self._clear_scene()
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
        else:
            assert source is not None
            self.tabulator = Tabulator(source, only_molecule=self._only_molecule)

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
        self.source_ready.emit()

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
            cartesian = self._config.grid.cartesian
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
        if index < -1 or index >= len(self.tabulator.molecular_orbitals):
            raise IndexError(f'Orbital index out of range: {index}')
        self.plot_orbital(index)
        self.controls.current_mo_ind = index
        self.controls.highlight_orbital(index)
        self.controls.update_nav_button_states()
        self.orbital_changed.emit(index)

    def update_grid(
        self,
        axes: tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]],
        grid_type: GridType,
    ) -> None:
        """Replace the structured grid and schedule fresh GTO tabulation."""
        if self._only_molecule:
            raise RuntimeError('Molecule-only viewers do not have an orbital grid.')
        self._update_mesh(*axes, grid_type)

    def apply_appearance(
        self,
        *,
        contour: float,
        mo_opacity: float,
        color_scheme: str,
        custom_colors: list[str],
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
        self._cmap = self._custom_cmap_from_colors(custom_colors) if color_scheme == 'custom' else color_scheme
        self.interactor.set_background(background_color)
        if hasattr(self, 'tabulator'):
            self._load_molecule(self._config)
            if not self._only_molecule and self.controls.current_mo_ind >= 0 and self._gtos_ready:
                self.plot_orbital(self.controls.current_mo_ind)
        self.interactor.update()

    def set_background_color(self, color: str) -> None:
        """Set this viewer's render background color."""
        if not mcolors.is_color_like(color):
            raise ValueError(f'Invalid background color: {color}')
        self._config.config.background_color = color
        self.controls.background_color.setText(color)
        self.interactor.set_background(color)

    def export_data(self, path: str | Path, *, file_format: str, scope: str = 'current') -> None:
        """Export orbital data to an explicit destination without dialogs."""
        if file_format not in {'vtk', 'cube'}:
            raise ValueError(f'Unsupported orbital export format: {file_format}')
        if scope not in {'current', 'all'}:
            raise ValueError(f'Unsupported orbital export scope: {scope}')
        if file_format == 'cube' and scope == 'all':
            raise ValueError('Cube format only supports one orbital.')
        index = self.controls.current_mo_ind
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
        self.controls.on_gtos_ready()
        if self.controls.current_mo_ind >= 0:
            self.plot_orbital(self.controls.current_mo_ind)
        self.loading_changed.emit(False)
        logger.info('GTO tabulation completed in %.2fs.', elapsed)

    def _handle_gto_error(self, exc: Exception) -> None:
        self._gtos_ready = self.tabulator.has_gtos
        self.controls.on_gtos_ready()
        self.loading_changed.emit(False)
        self.report_error('Orbital tabulation failed', exc)

    def _dispatch_gto_completion(self, callback: Callable[[], None]) -> None:
        if self._on_screen:
            self._dispatcher.dispatch(callback)

    def _ensure_gtos_ready(self) -> bool:
        return self._gtos_ready

    def _update_settings_button_states(self) -> None:
        """Mirror actor visibility in the public control widgets."""
        self.controls.show_atoms.setChecked(self.are_atoms_visible())
        self.controls.show_bonds.setChecked(self.are_bonds_visible())

    def _cancel_gto_future(self) -> None:
        self._gto_job.cancel()

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
