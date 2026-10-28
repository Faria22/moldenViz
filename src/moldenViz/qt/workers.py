"""Background computation and Qt completion delivery."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from PySide6.QtCore import QObject, Qt, Signal, Slot

from moldenViz._adaptive_grid import (
    crossed_cell_ids,
    refined_grid,
)
from moldenViz.tabulator import GridType, Tabulator

if TYPE_CHECKING:
    from collections.abc import Callable

    import pyvista as pv
    from numpy.typing import NDArray


_GTO_EXECUTOR = ThreadPoolExecutor(max_workers=1)


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


def tabulate_grid(
    tabulator: Tabulator,
    frozen_axes: tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]],
    grid_type: GridType,
    grid: NDArray[np.floating] | None,
) -> _GTOResult:
    """Compute GTO values using a structured grid snapshot.

    Returns
    -------
    _GTOResult
        Grid, axes, and GTO values for installation on the GUI thread.
    """
    if grid is None:
        grid = Tabulator._build_grid(  # ruff:ignore[private-member-access]
            frozen_axes[0],
            frozen_axes[1],
            frozen_axes[2],
            grid_type,
        )
    grid.setflags(write=False)
    return _GTOResult(
        grid=grid,
        axes=frozen_axes,
        grid_type=grid_type,
        gtos=tabulator.compute_gtos(grid),
    )


def tabulate_adaptive_grid(
    coarse_grid: pv.StructuredGrid,
    coarse_mos: NDArray[np.floating],
    contour: float,
    scale: tuple[float, float, float],
    tabulator: Tabulator,
) -> _AdaptiveResult:
    """Refine coarse contours and compute their reusable GTO values.

    Returns
    -------
    _AdaptiveResult
        Refined contour mesh and GTO values for GUI-thread installation.
    """
    cell_ids = crossed_cell_ids(coarse_grid, coarse_mos, contour)
    mesh = refined_grid(coarse_grid, cell_ids, scale)
    gtos = tabulator.compute_gtos(np.asarray(mesh.points))
    return _AdaptiveResult(mesh=mesh, gtos=gtos, crossed_cells=len(cell_ids))
