"""Native export dialogs invoked by the standalone window."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from PySide6.QtWidgets import QFileDialog

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping

    from PySide6.QtWidgets import QWidget

    from .viewer import OrbitalViewer


def export_data_dialog(
    parent: QWidget,
    viewer: OrbitalViewer,
    options: Mapping[str, Any],
    on_error: Callable[[str, object], None],
) -> None:
    """Choose a volume destination and export through the host viewer."""
    file_format = str(options.get('format', 'vtk'))
    scope = str(options.get('scope', 'current'))
    if scope == 'current' and viewer.current_orbital_index < 0:
        on_error('Export failed', ValueError('No orbital is currently selected.'))
        return
    suffix = f'.{file_format}'
    default_name = f'orbital{suffix}' if scope == 'current' else f'orbitals_all{suffix}'
    filters = 'VTK files (*.vtk)' if file_format == 'vtk' else 'Gaussian Cube files (*.cube)'
    destination, _selected = QFileDialog.getSaveFileName(parent, 'Export orbital data', default_name, filters)
    if not destination:
        return
    try:
        viewer.export_data(destination, file_format=file_format, scope=scope)
    except (OSError, RuntimeError, ValueError) as exc:
        on_error('Export failed', exc)


def export_image_dialog(
    parent: QWidget,
    viewer: OrbitalViewer,
    options: Mapping[str, Any],
    on_error: Callable[[str, object], None],
) -> None:
    """Choose an image destination and export through the host viewer."""
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
        parent,
        'Export image',
        f'moldenviz_export{suffix}',
        filters[file_format],
    )
    if not destination:
        return
    try:
        viewer.export_image(
            Path(destination),
            file_format=file_format,
            transparent=transparent,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        on_error('Export failed', exc)
