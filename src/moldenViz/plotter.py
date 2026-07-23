"""Plotter module for creating plots of the molecule and it's orbitals."""

from __future__ import annotations

import logging
import tkinter as tk
from concurrent.futures import ThreadPoolExecutor
from tkinter import messagebox
from typing import TYPE_CHECKING

import numpy as np
from pyvistaqt import BackgroundPlotter

from ._config_module import Config
from ._plotter_jobs import BackgroundJob
from ._plotter_rendering import _PlotterRendering
from ._plotter_ui import _OrbitalSelectionScreen, _PlotterUI
from .tabulator import GridType, Tabulator

if TYPE_CHECKING:
    from collections.abc import Callable
    from concurrent.futures import Future

    import pyvista as pv
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
        If user is using the plotter inside a tk app, `tk_root` can be passed to not create a new tk instance.

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

        self._tk_root = tk_root
        self._no_prev_tk_root = self._tk_root is None
        if self._tk_root is None:
            self._tk_root = tk.Tk()
            self._tk_root.withdraw()  # Hides window
            logger.debug('Created internal Tk root window for Plotter UI.')

        self._gto_job: BackgroundJob[NDArray[np.floating]] = BackgroundJob(
            _GTO_EXECUTOR,
            self._dispatch_gto_completion,
        )

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
            if config.grid.default_type == 'spherical':
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
            else:  # cartesian
                r = max(config.grid.max_radius_multiplier * self._molecule.max_radius, config.grid.min_radius)
                logger.info(
                    'Generating default cartesian grid spanning ±%.2f with %dx%dx%d samples.',
                    r,
                    config.grid.cartesian.num_x_points,
                    config.grid.cartesian.num_y_points,
                    config.grid.cartesian.num_z_points,
                )
                self.tabulator.cartesian_grid(
                    np.linspace(-r, r, config.grid.cartesian.num_x_points),
                    np.linspace(-r, r, config.grid.cartesian.num_y_points),
                    np.linspace(-r, r, config.grid.cartesian.num_z_points),
                    tabulate_gtos=False,
                )
            self._gtos_ready = False
            self._schedule_gto_tabulation()

        # If we want to have the molecular orbitals, we need to initiate Tk before Qt
        # That is why we have this weird if statement separated this way
        if only_molecule:
            logger.info('Running in molecule-only mode; skipping orbital mesh creation.')
            if self._no_prev_tk_root:
                logger.debug('Entering Tk main loop for molecule-only display.')
                self._tk_root.mainloop()
            return

        self._orb_mesh = self._create_mo_mesh()
        self._orb_actor: pv.Actor | None = None

        # Values for MO, not the molecule
        self._contour = config.mo.contour
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
            logger.debug('Entering Tk main loop for full Plotter UI.')
            self._tk_root.mainloop()

    def wait_for_gtos(self, timeout: float | None = None) -> None:
        """Block until the background GTO tabulation finishes."""
        if self._gtos_ready:
            return
        if self._gto_job.future is None:
            raise RuntimeError('GTO tabulation has not been scheduled.')
        try:
            gtos = self._gto_job.wait(timeout=timeout)
        except RuntimeError:
            if self._gtos_ready:
                return
            raise
        if not self._gtos_ready:
            self._apply_gtos_ready(gtos, 0.0)

    @property
    def _gto_future(self) -> Future[NDArray[np.floating]] | None:
        """Compatibility view of the pending background future."""
        return self._gto_job.future

    def _dispatch_gto_completion(self, callback: Callable[[], None]) -> None:
        """Schedule a completion callback on the owning Tk event loop."""
        if self._tk_root is not None:
            self._tk_root.after_idle(callback)
        else:
            callback()

    def _schedule_gto_tabulation(self) -> None:
        """Submit background GTO tabulation work."""
        if self._only_molecule or self._gtos_ready or self._gto_job.pending:
            return
        grid = self.tabulator.grid.copy()
        grid.setflags(write=False)
        logger.info('Starting background GTO tabulation...')
        self._gto_job.start(
            lambda: self.tabulator.compute_gtos(grid),
            on_success=self._apply_gtos_ready,
            on_error=self._handle_gto_error,
        )

    @staticmethod
    def _handle_gto_error(exc: Exception) -> None:
        """Report a failed GTO job after it returns to the UI thread."""
        logger.error(
            'Background GTO tabulation failed.',
            exc_info=(type(exc), exc, exc.__traceback__),
        )
        messagebox.showerror('Orbital Tabulation Failed', f'Failed to tabulate orbitals:\n\n{exc!s}')

    def _apply_gtos_ready(self, gtos: NDArray[np.floating], elapsed: float) -> None:
        """Store computed GTOs and update UI state."""
        if not self._on_screen:
            return
        self.tabulator.set_gtos(gtos)
        self._gtos_ready = True
        logger.info('GTO tabulation completed in %.2fs.', elapsed)
        self._orb_mesh = self._create_mo_mesh()
        if self._selection_screen:
            self._selection_screen._on_gtos_ready()  # ruff:ignore[private-member-access]
            if self._selection_screen.current_mo_ind >= 0:
                self.plot_orbital(self._selection_screen.current_mo_ind)

    def _ensure_gtos_ready(self) -> bool:
        """Return True if GTO data are ready for orbital operations.

        Returns
        -------
        bool
            True when orbital plots can be rendered immediately.
        """
        if self._gtos_ready:
            return True
        logger.debug('Ignoring orbital request while GTOs are loading.')
        return False

    def _cancel_gto_future(self) -> None:
        """Cancel any pending GTO computation."""
        if not self._gto_job.pending:
            return
        future = self._gto_job.future
        if future is not None and not future.done():
            logger.info('Cancelling pending GTO tabulation job.')
        self._gto_job.cancel()
