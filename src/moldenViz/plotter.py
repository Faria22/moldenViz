"""Plotter module for creating plots of the molecule and it's orbitals."""

from __future__ import annotations

import logging
import sys
import tkinter as tk
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from queue import SimpleQueue
from tkinter import messagebox
from typing import TYPE_CHECKING

import numpy as np
import pyvista as pv
from pyvistaqt import BackgroundPlotter

from ._adaptive_grid import AdaptiveScale, crossed_cell_ids, normalize_scale, refined_grid
from ._config_module import Config
from ._plotter_jobs import BackgroundJob
from ._plotter_rendering import _PlotterRendering
from ._plotter_ui import _OrbitalSelectionScreen, _PlotterUI
from .tabulator import GridType, Tabulator

if TYPE_CHECKING:
    from collections.abc import Callable
    from concurrent.futures import Future

    from numpy.typing import NDArray


def _describe_source(source: str | list[str]) -> str:
    """Return a human readable description of the data source.

    Parameters
    ----------
    source : str | list[str]
        Path to a Molden file or the raw lines read from one.

    Returns
    -------
    str
        Description suitable for logging output.
    """
    if isinstance(source, str):
        return source
    return f'{len(source)} molden lines'


logger = logging.getLogger(__name__)

__all__ = ['Plotter']

config = Config()
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


class Plotter(_PlotterUI, _PlotterRendering):
    """
    Handles the 3D visualization of molecules and molecular orbitals.

    This class uses PyVista for 3D rendering and Tkinter for the user interface
    to control plotting parameters and select orbitals.

    Parameters
    ----------
    source : str | list[str]
        The path to the molden file, or the lines from the file.
    only_molecule : bool, optional
        Only parse the atoms and skip molecular orbitals.
        Default is `False`.
    tabulator : Tabulator, optional
        If `None`, `Plotter` creates a `Tabulator` and tabulates the GTOs and MOs with a default grid.
        A `Tabulator` can be passed to reuse a predetermined grid. When `only_molecule` is `False`,
        the supplied `Tabulator` must already have tabulated GTOs available through `tabulator.gtos`.

        Note: `Tabulator` grid must be spherical or cartesian. Custom grids are not allowed.
    tk_root : tk.Tk, optional
        If user is using the plotter inside a Tk app, `tk_root` can be passed
        to avoid creating a new Tk instance. The caller retains ownership of a
        supplied root and must keep its event loop running while background GTO
        tabulation is active. Plotter does not quit or destroy a supplied root.

    Attributes
    ----------
    tabulator : Tabulator
        The Tabulator object used for tabulating GTOs and MOs.

    Raises
    ------
    ValueError
        If the provided tabulator is invalid
        (e.g., missing grid or GTO data when `only_molecule` is `False`, or has an UNKNOWN grid type).
    """

    _SPHERICAL_GRID_SETTINGS_WINDOW_SIZE = '400x350'
    _CARTESIAN_GRID_SETTINGS_WINDOW_SIZE = '650x400'
    _GTO_COMPLETION_POLL_MS = 10
    _TK_UPDATE_MS = 10

    def __init__(
        self,
        source: str | list[str],
        only_molecule: bool = False,
        tabulator: Tabulator | None = None,
        tk_root: tk.Tk | None = None,
    ) -> None:
        logger.info('Initialising Plotter (only_molecule=%s)', only_molecule)

        self._on_screen = True
        self._only_molecule = only_molecule
        self._selection_screen: _OrbitalSelectionScreen | None = None
        self._gtos_ready = only_molecule
        self._grid_mode = tabulator.grid_type.value if tabulator is not None else config.grid.default_type
        self._adaptive_scale: AdaptiveScale = config.grid.adaptive.scale
        self._adaptive_ready = False
        self._adaptive_mesh = pv.UnstructuredGrid()
        self._adaptive_gtos: NDArray[np.floating] | None = None
        self._contour = config.mo.contour

        self._tk_root = tk_root
        self._no_prev_tk_root = self._tk_root is None
        if self._tk_root is None:
            self._tk_root = tk.Tk()
            self._tk_root.withdraw()  # Hides window
            logger.debug('Created internal Tk root window for Plotter UI.')

        self._gto_completions: SimpleQueue[Callable[[], None]] = SimpleQueue()
        self._gto_completion_poll_id: str | None = None
        self._gto_job: BackgroundJob[_GTOResult] = BackgroundJob(
            _GTO_EXECUTOR,
            self._dispatch_gto_completion,
        )
        self._adaptive_job: BackgroundJob[_AdaptiveResult] = BackgroundJob(
            _GTO_EXECUTOR,
            self._dispatch_gto_completion,
        )
        self._schedule_gto_completion_poll()

        self._pv_plotter = BackgroundPlotter(editor=False)
        self._pv_plotter.set_background(config.background_color)
        self._pv_plotter.show_axes()
        logger.debug('Configured PyVista plotter background colour to %s', config.background_color)

        self._add_orbital_menus_to_pv_plotter()
        self._connect_pv_plotter_close_signal()
        self._override_clear_all_button()

        if tabulator:
            logger.info('Using provided Tabulator instance with grid type %s', tabulator.grid_type.value)
            if not hasattr(tabulator, 'grid'):
                raise ValueError('Tabulator does not have grid attribute.')

            if not tabulator.has_gtos and not only_molecule:
                raise ValueError('Tabulator does not have tabulated GTOs.')

            if tabulator.grid_type == GridType.UNKNOWN:
                raise ValueError('The plotter only supports spherical and cartesian grids.')

            self.tabulator = tabulator
        else:
            logger.info('Creating Tabulator for source %s', _describe_source(source))
            self.tabulator = Tabulator(source, only_molecule=only_molecule)
        self._gtos_ready = self._only_molecule or self.tabulator.has_gtos

        self._molecule_opacity = config.molecule.opacity
        self._load_molecule(config)

        # If no tabulator was passed, create default grid
        if not only_molecule and not tabulator:
            if self._grid_mode == 'spherical':
                logger.info(
                    'Generating default spherical grid with %dx%dx%d samples.',
                    config.grid.spherical.num_r_points,
                    config.grid.spherical.num_theta_points,
                    config.grid.spherical.num_phi_points,
                )
                self.tabulator.spherical_grid(
                    np.linspace(
                        0,
                        max(config.grid.max_radius_multiplier * self._molecule.max_radius, config.grid.min_radius),
                        config.grid.spherical.num_r_points,
                    ),
                    np.linspace(0, np.pi, config.grid.spherical.num_theta_points),
                    np.linspace(0, 2 * np.pi, config.grid.spherical.num_phi_points),
                    tabulate_gtos=False,
                )
            else:  # cartesian or adaptive
                r = max(config.grid.max_radius_multiplier * self._molecule.max_radius, config.grid.min_radius)
                grid_config = config.grid.adaptive if self._grid_mode == 'adaptive' else config.grid.cartesian
                logger.info(
                    'Generating default %s Cartesian grid spanning ±%.2f with %dx%dx%d samples.',
                    self._grid_mode,
                    r,
                    grid_config.num_x_points,
                    grid_config.num_y_points,
                    grid_config.num_z_points,
                )
                self.tabulator.cartesian_grid(
                    np.linspace(-r, r, grid_config.num_x_points),
                    np.linspace(-r, r, grid_config.num_y_points),
                    np.linspace(-r, r, grid_config.num_z_points),
                    tabulate_gtos=False,
                )
            self._gtos_ready = False
            self._schedule_gto_tabulation()

        # If we want to have the molecular orbitals, we need to initiate Tk before Qt
        # That is why we have this weird if statement separated this way
        if only_molecule:
            logger.info('Running in molecule-only mode; skipping orbital mesh creation.')
            if self._no_prev_tk_root:
                self._run_internal_event_loop()
            return

        self._orb_mesh = self._create_mo_mesh()
        self._orb_actor: pv.Actor | None = None

        # Values for MO, not the molecule
        self._opacity = config.mo.opacity

        # Set colormap based on configuration
        if config.mo.custom_colors:
            # Create custom colormap from two colors
            self._cmap = self._custom_cmap_from_colors(config.mo.custom_colors)
        else:
            self._cmap = config.mo.color_scheme

        if not self._only_molecule:
            self._selection_screen = _OrbitalSelectionScreen(self)
            logger.debug('Orbital selection screen initialised.')
            if not self._gtos_ready:
                self._selection_screen._set_loading_state(True)  # ruff:ignore[private-member-access]

        if self._no_prev_tk_root:
            self._run_internal_event_loop()

    def _run_internal_event_loop(self) -> None:
        """Run the native GUI loop for a Plotter-owned Tk root."""
        if sys.platform != 'darwin':
            logger.debug('Entering Tk main loop.')
            self._tk_root.mainloop()
            return

        logger.debug('Entering Qt main loop and polling Tk events on macOS.')
        self._pv_plotter.add_callback(
            self._tk_root.update,
            interval=self._TK_UPDATE_MS,
        )
        self._pv_plotter.app.exec()

    def wait_for_gtos(self, timeout: float | None = None) -> None:
        """Block until the background GTO tabulation finishes."""
        if self._gtos_ready:
            return
        if self._gto_job.future is None:
            raise RuntimeError('GTO tabulation has not been scheduled.')
        try:
            result = self._gto_job.wait(timeout=timeout)
        except RuntimeError:
            if self._gtos_ready:
                return
            raise
        if not self._gtos_ready:
            self._apply_gtos_ready(result, 0.0)

    def wait_for_adaptive_grid(self, timeout: float | None = None) -> None:
        """Block until adaptive-grid tabulation finishes."""
        if self._grid_mode != 'adaptive' or self._adaptive_ready:
            return
        future = self._adaptive_job.future
        if future is None:
            raise RuntimeError('Adaptive-grid tabulation has not been scheduled.')
        result = self._adaptive_job.wait(timeout=timeout)
        if not self._adaptive_ready:
            self._apply_adaptive_ready(result, 0.0)

    @property
    def _gto_future(self) -> Future[_GTOResult] | None:
        """Compatibility view of the pending background future."""
        return self._gto_job.future

    def _dispatch_gto_completion(self, callback: Callable[[], None]) -> None:
        """Publish a completion callback without interacting with Tk."""
        if self._on_screen:
            self._gto_completions.put(callback)

    def _schedule_gto_completion_poll(self) -> None:
        """Schedule completion polling from the Tk-owning thread."""
        if self._tk_root is None or not self._on_screen:
            return
        self._gto_completion_poll_id = self._tk_root.after(
            self._GTO_COMPLETION_POLL_MS,
            self._poll_gto_completions,
        )

    def _poll_gto_completions(self) -> None:
        """Deliver queued background completions on the Tk-owning thread."""
        self._gto_completion_poll_id = None
        if not self._on_screen:
            return
        while not self._gto_completions.empty():
            callback = self._gto_completions.get_nowait()
            callback()
            if not self._on_screen:
                return
        self._schedule_gto_completion_poll()

    def _stop_gto_completion_poll(self) -> None:
        """Stop polling and discard callbacks after the UI closes."""
        if self._tk_root is not None and self._gto_completion_poll_id is not None:
            self._tk_root.after_cancel(self._gto_completion_poll_id)
            self._gto_completion_poll_id = None
        while not self._gto_completions.empty():
            self._gto_completions.get_nowait()

    def _schedule_gto_tabulation(
        self,
        axes: tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]] | None = None,
        grid_type: GridType | None = None,
    ) -> None:
        """Build a grid and tabulate its GTOs in the background."""
        if self._only_molecule or self._gtos_ready or self._gto_job.pending:
            return

        if axes is None:
            current_axes = self.tabulator.grid_axes
            if current_axes is None:
                raise RuntimeError('Structured grid axes are not available.')
            resolved_axes = (current_axes[0], current_axes[1], current_axes[2])
            resolved_grid_type = self.tabulator.grid_type
            current_grid = self.tabulator.grid.copy()
        elif grid_type is None:
            raise ValueError('Grid type is required when scheduling new axes.')
        else:
            resolved_axes = axes
            resolved_grid_type = grid_type
            current_grid = None

        frozen_axes = (
            resolved_axes[0].copy(),
            resolved_axes[1].copy(),
            resolved_axes[2].copy(),
        )
        for axis in frozen_axes:
            axis.setflags(write=False)

        def build_and_tabulate() -> _GTOResult:
            grid = current_grid
            if grid is None:
                grid = Tabulator._build_grid(  # ruff:ignore[private-member-access]
                    *frozen_axes,
                    resolved_grid_type,
                )
            grid.setflags(write=False)
            return _GTOResult(
                grid=grid,
                axes=frozen_axes,
                grid_type=resolved_grid_type,
                gtos=self.tabulator.compute_gtos(grid),
            )

        logger.info('Starting background GTO tabulation...')
        self._gto_job.start(
            build_and_tabulate,
            on_success=self._apply_gtos_ready,
            on_error=self._handle_gto_error,
        )

    def _handle_gto_error(self, exc: Exception) -> None:
        """Restore usable UI state and report a failed GTO job."""
        self._gtos_ready = self.tabulator.has_gtos
        if self._selection_screen:
            self._selection_screen._on_gtos_ready()  # ruff:ignore[private-member-access]
        logger.error(
            'Background GTO tabulation failed.',
            exc_info=(type(exc), exc, exc.__traceback__),
        )
        messagebox.showerror('Orbital Tabulation Failed', f'Failed to tabulate orbitals:\n\n{exc!s}')

    def _apply_gtos_ready(self, result: _GTOResult, elapsed: float) -> None:
        """Store computed GTOs and update UI state."""
        if not self._on_screen:
            return
        self.tabulator._set_structured_grid(  # ruff:ignore[private-member-access]
            result.grid,
            result.axes,
            result.grid_type,
        )
        self.tabulator.set_gtos(result.gtos)
        self._gtos_ready = True
        logger.info('GTO tabulation completed in %.2fs.', elapsed)
        self._orb_mesh = self._create_mo_mesh()
        if self._grid_mode == 'adaptive':
            self._schedule_adaptive_tabulation()
            return
        if self._selection_screen:
            self._selection_screen._on_gtos_ready()  # ruff:ignore[private-member-access]
            if self._selection_screen.current_mo_ind >= 0:
                self.plot_orbital(self._selection_screen.current_mo_ind)

    def _schedule_adaptive_tabulation(self) -> None:
        """Build the union of all coarse MO contours and tabulate it once."""
        if self._only_molecule or self._grid_mode != 'adaptive' or not self._gtos_ready or self._adaptive_job.pending:
            return
        if self.tabulator.grid_type != GridType.CARTESIAN:
            raise RuntimeError('Adaptive grids require a Cartesian coarse grid.')

        coarse_grid = self._create_mo_mesh()
        coarse_mos = self.tabulator.tabulate_mos()
        contour = self._contour
        scale = normalize_scale(self._adaptive_scale)

        def build_and_tabulate() -> _AdaptiveResult:
            cell_ids = crossed_cell_ids(coarse_grid, coarse_mos, contour)
            mesh = refined_grid(coarse_grid, cell_ids, scale)
            gtos = self.tabulator.compute_gtos(np.asarray(mesh.points))
            return _AdaptiveResult(mesh=mesh, gtos=gtos, crossed_cells=len(cell_ids))

        self._adaptive_ready = False
        if self._selection_screen:
            self._selection_screen._set_loading_state(  # ruff:ignore[private-member-access]
                True,
                'Refining contour grid...',
            )
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
        if self._selection_screen:
            self._selection_screen._on_gtos_ready()  # ruff:ignore[private-member-access]
            if self._selection_screen.current_mo_ind >= 0:
                self.plot_orbital(self._selection_screen.current_mo_ind)

    def _handle_adaptive_error(self, exc: Exception) -> None:
        """Keep the previous actor and report failed adaptive tabulation."""
        self._adaptive_ready = False
        if self._selection_screen:
            self._selection_screen._on_gtos_ready()  # ruff:ignore[private-member-access]
        logger.error('Adaptive-grid tabulation failed.', exc_info=(type(exc), exc, exc.__traceback__))
        messagebox.showerror('Adaptive Grid Failed', f'Failed to refine orbital grid:\n\n{exc!s}')

    def _invalidate_adaptive_grid(self, *, rebuild: bool = False) -> None:
        """Discard adaptive cache and optionally rebuild it from the coarse grid."""
        self._adaptive_job.cancel()
        self._adaptive_ready = False
        self._adaptive_mesh = pv.UnstructuredGrid()
        self._adaptive_gtos = None
        if rebuild:
            self._schedule_adaptive_tabulation()

    def _ensure_gtos_ready(self) -> bool:
        """Return True if GTO data are ready for orbital operations.

        Returns
        -------
        bool
            True when orbital plots can be rendered immediately.
        """
        if self._gtos_ready and (self._grid_mode != 'adaptive' or self._adaptive_ready):
            return True
        logger.debug('Ignoring orbital request while GTOs are loading.')
        return False

    def _cancel_gto_future(self) -> None:
        """Cancel any pending GTO computation."""
        if not self._on_screen:
            self._stop_gto_completion_poll()
        for label, job in (('GTO', self._gto_job), ('adaptive-grid', self._adaptive_job)):
            future = job.future
            if future is not None and not future.done():
                logger.info('Cancelling pending %s tabulation job.', label)
            job.cancel()
