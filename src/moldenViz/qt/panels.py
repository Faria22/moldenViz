"""Orbital selection and viewer control panels."""

from __future__ import annotations

import logging
from importlib import import_module
from typing import TYPE_CHECKING, Any

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QAbstractSpinBox,
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
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from moldenViz._adaptive_grid import (
    format_scale,
    parse_scale,
)
from moldenViz.tabulator import GridType

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from moldenViz.models import MolecularOrbital

    from .viewer import OrbitalViewer

logger = logging.getLogger('moldenViz.qt')

_MO_COLOR_SCHEMES = ['bwr', 'RdBu', 'seismic', 'coolwarm', 'PiYG']
_BACKGROUND_COLORS = (
    ('White', 'white'),
    ('Black', 'black'),
    ('Light gray', '#A0A0A0'),
    ('Dark gray', '#202124'),
)
_ORBITAL_COLUMN_PADDING = 8
_MIN_ADAPTIVE_POINTS = 2


def _is_color_like(color: object) -> bool:
    """Return whether Matplotlib accepts a color specification.

    Returns
    -------
    bool
        Whether ``color`` is valid.
    """
    return bool(import_module('matplotlib.colors').is_color_like(color))


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
        self.reset_camera_button = QPushButton('Reset camera', self)
        self.reset_camera_button.clicked.connect(self._reset_camera)
        layout.addWidget(self.reset_camera_button)

    def _reset_camera(self) -> None:
        """Restore the viewer's default camera orientation and framing."""
        self._viewer.reset_camera()

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
        self.show_axes.toggled.connect(lambda visible: self._viewer.set_axes_visible(visible=visible))
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
        self.orbital_table.blockSignals(True)  # ruff: ignore[boolean-positional-value-in-call]
        self.orbital_table.setRowCount(len(orbitals))
        for row, orbital in enumerate(orbitals):
            values = (str(row + 1), orbital.sym, orbital.spin, f'{orbital.occ:g}', f'{orbital.energy:.6f}')
            for column, value in enumerate(values):
                self.orbital_table.setItem(row, column, QTableWidgetItem(value))
        self.orbital_table.blockSignals(False)  # ruff: ignore[boolean-positional-value-in-call]
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

    def set_molecule_only(self, *, only_molecule: bool) -> None:
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
        """Select the table row matching the viewer orbital."""
        self.orbital_table.blockSignals(True)  # ruff: ignore[boolean-positional-value-in-call]
        if index < 0:
            self.orbital_table.clearSelection()
        else:
            self.orbital_table.selectRow(index)
        self.orbital_table.blockSignals(False)  # ruff: ignore[boolean-positional-value-in-call]

    def set_loading_state(self, *, loading: bool, message: str = 'Tabulating orbitals…') -> None:
        """Show computation progress and refresh navigation availability."""
        self.loading_label.setText(message)
        self.loading_label.setVisible(loading)
        self.orbital_table.setEnabled(not loading)
        self.update_nav_button_states()

    def on_gtos_ready(self) -> None:
        """Restore controls after background computation finishes."""
        self.set_loading_state(loading=False)

    def update_nav_button_states(self) -> None:
        """Enable orbital navigation when the selected data is ready."""
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
