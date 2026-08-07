"""Standalone Qt window for the embeddable orbital viewer."""

from __future__ import annotations

import logging
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

from PySide6.QtGui import QAction, QCloseEvent
from PySide6.QtWidgets import QApplication, QFileDialog, QMainWindow, QMessageBox, QWidget

from .qt import OrbitalViewer

if TYPE_CHECKING:
    from concurrent.futures import Future

    from ._config_module import Config, MainConfig
    from .qt import _GTOResult
    from .tabulator import Tabulator

logger = logging.getLogger(__name__)

__all__ = ['Plotter']

# Qt does not retain an unparented top-level Python wrapper. Keep windows made
# inside an existing application alive until Qt emits ``destroyed``.
_OPEN_WINDOWS: set[Plotter] = set()


class Plotter(QMainWindow):
    """Free-floating Qt window containing an :class:`OrbitalViewer`.

    If no :class:`QApplication` exists, ``Plotter`` creates one and blocks in
    its event loop until the window closes. Inside an existing Qt application,
    construction returns immediately and the host keeps ownership of the loop.

    Parameters
    ----------
    source : str | list[str]
        Molden path or raw Molden lines.
    only_molecule : bool, optional
        Skip molecular-orbital parsing and controls.
    tabulator : Tabulator, optional
        Existing structured-grid tabulator with cached GTO values.
    config : Config | MainConfig | Mapping, optional
        Per-window configuration overrides.
    parent : QWidget, optional
        Optional Qt parent for the free-floating window.
    """

    def __init__(
        self,
        source: str | list[str],
        only_molecule: bool = False,
        tabulator: Tabulator | None = None,
        *,
        config: Config | MainConfig | Mapping[str, Any] | None = None,
        parent: QWidget | None = None,
    ) -> None:
        application = QApplication.instance()
        self._owns_application = application is None
        if application is None:
            application = QApplication(sys.argv[:1])
        elif not application.inherits('QApplication'):
            raise RuntimeError('Plotter requires QApplication, but a non-GUI QCoreApplication already exists.')
        self._application = cast(QApplication, application)

        super().__init__(parent)
        self.setWindowTitle('moldenViz')
        self.resize(1200, 760)
        self.viewer = OrbitalViewer(
            source,
            only_molecule=only_molecule,
            tabulator=tabulator,
            config=config,
            parent=self,
        )
        self.setCentralWidget(self.viewer)
        self._pv_plotter = self.viewer.interactor
        self.viewer.error_occurred.connect(self._show_error)
        self.viewer.export_requested.connect(self._handle_export_request)
        self._build_menus(only_molecule)
        _OPEN_WINDOWS.add(self)
        self.destroyed.connect(lambda: _OPEN_WINDOWS.discard(self))
        self.show()

        if self._owns_application:
            self._application.exec()

    @property
    def tabulator(self) -> Tabulator:
        """Tabulator used by the contained viewer."""
        return self.viewer.tabulator

    @property
    def _gto_future(self) -> Future[_GTOResult] | None:
        """Compatibility view of pending GTO computation."""
        return self.viewer._gto_future  # ruff:ignore[private-member-access]

    def _build_menus(self, only_molecule: bool) -> None:
        view_menu = self.menuBar().addMenu('View')
        clear_action = QAction('Clear orbital', self)
        clear_action.setEnabled(not only_molecule)
        clear_action.triggered.connect(lambda: self.viewer.show_orbital(-1))
        view_menu.addAction(clear_action)
        reset_camera = QAction('Reset camera', self)
        reset_camera.triggered.connect(self.viewer.interactor.reset_camera)
        view_menu.addAction(reset_camera)

        settings_menu = self.menuBar().addMenu('Settings')
        appearance_tab = 2
        for label, tab_index in (('Grid', 1), ('Appearance', appearance_tab)):
            action = QAction(label, self)
            action.setEnabled(not only_molecule or tab_index == appearance_tab)
            action.triggered.connect(
                lambda _checked=False, index=tab_index: self.viewer.controls.tabs.setCurrentIndex(index),
            )
            settings_menu.addAction(action)
        save_action = QAction('Save settings', self)
        save_action.triggered.connect(self._save_settings)
        settings_menu.addAction(save_action)

        export_menu = self.menuBar().addMenu('Export')
        data_action = QAction('Orbital data…', self)
        data_action.setEnabled(not only_molecule)
        data_action.triggered.connect(
            lambda: self._handle_export_request(
                'data',
                {
                    'format': self.viewer.controls.data_format.currentText(),
                    'scope': self.viewer.controls.data_scope.currentText(),
                },
            ),
        )
        export_menu.addAction(data_action)
        image_action = QAction('Image…', self)
        image_action.triggered.connect(
            lambda: self._handle_export_request(
                'image',
                {
                    'format': self.viewer.controls.image_format.currentText(),
                    'transparent': self.viewer.controls.transparent_background.isChecked(),
                },
            ),
        )
        export_menu.addAction(image_action)

    def _save_settings(self) -> None:
        try:
            self.viewer.save_settings()
        except OSError as exc:
            self._show_error('Saving settings failed', exc)
        else:
            QMessageBox.information(self, 'Settings saved', 'Configuration saved successfully.')

    def _handle_export_request(self, kind: str, options: object) -> None:
        values = dict(options) if isinstance(options, Mapping) else {}
        if kind == 'data':
            self._export_data(values)
        elif kind == 'image':
            self._export_image(values)
        else:
            self._show_error('Export failed', ValueError(f'Unknown export kind: {kind}'))

    def _export_data(self, options: Mapping[str, Any]) -> None:
        file_format = str(options.get('format', 'vtk'))
        scope = str(options.get('scope', 'current'))
        if scope == 'current' and self.viewer.controls.current_mo_ind < 0:
            self._show_error('Export failed', ValueError('No orbital is currently selected.'))
            return
        suffix = f'.{file_format}'
        default_name = f'orbital{suffix}' if scope == 'current' else f'orbitals_all{suffix}'
        filters = 'VTK files (*.vtk)' if file_format == 'vtk' else 'Gaussian Cube files (*.cube)'
        destination, _selected = QFileDialog.getSaveFileName(self, 'Export orbital data', default_name, filters)
        if not destination:
            return
        try:
            self.viewer.export_data(destination, file_format=file_format, scope=scope)
        except (OSError, RuntimeError, ValueError) as exc:
            self._show_error('Export failed', exc)

    def _export_image(self, options: Mapping[str, Any]) -> None:
        file_format = str(options.get('format', 'png'))
        transparent = bool(options.get('transparent', False))
        suffix = '.jpg' if file_format == 'jpeg' else f'.{file_format}'
        filters = {
            'png': 'PNG files (*.png)',
            'jpeg': 'JPEG files (*.jpg *.jpeg)',
            'svg': 'SVG files (*.svg)',
            'pdf': 'PDF files (*.pdf)',
        }
        destination, _selected = QFileDialog.getSaveFileName(
            self,
            'Export image',
            f'moldenviz_export{suffix}',
            filters[file_format],
        )
        if not destination:
            return
        try:
            self.viewer.export_image(
                Path(destination),
                file_format=file_format,
                transparent=transparent,
            )
        except (OSError, RuntimeError, ValueError) as exc:
            self._show_error('Export failed', exc)

    def _show_error(self, title: str, exc: object) -> None:
        QMessageBox.critical(self, title, str(exc))

    def show_orbital(self, index: int) -> None:
        """Delegate orbital selection to the contained viewer."""
        self.viewer.show_orbital(index)

    def plot_orbital(self, index: int) -> None:
        """Compatibility alias for :meth:`show_orbital`."""
        self.viewer.show_orbital(index)

    def wait_for_gtos(self, timeout: float | None = None) -> None:
        """Wait for background GTO tabulation."""
        self.viewer.wait_for_gtos(timeout)

    def toggle_molecule(self) -> None:
        """Toggle all molecule actors."""
        self.viewer.toggle_molecule()

    def toggle_atoms(self) -> None:
        """Toggle atom actors."""
        self.viewer.toggle_atoms()

    def toggle_bonds(self) -> None:
        """Toggle bond actors."""
        self.viewer.toggle_bonds()

    def is_molecule_visible(self) -> bool:
        """Return whether molecule actors are visible.

        Returns
        -------
        bool
            Whether at least one molecule actor is visible.
        """
        return self.viewer.is_molecule_visible()

    def are_atoms_visible(self) -> bool:
        """Return whether atom actors are visible.

        Returns
        -------
        bool
            Whether at least one atom actor is visible.
        """
        return self.viewer.are_atoms_visible()

    def are_bonds_visible(self) -> bool:
        """Return whether bond actors are visible.

        Returns
        -------
        bool
            Whether at least one bond actor is visible.
        """
        return self.viewer.are_bonds_visible()

    def closeEvent(self, event: QCloseEvent) -> None:  # ruff: ignore[invalid-function-name]
        """Close VTK resources before the top-level Qt window is destroyed."""
        self.viewer.close()
        _OPEN_WINDOWS.discard(self)
        event.accept()
