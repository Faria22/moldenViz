"""Plotter module for creating plots of the molecule and it's orbitals."""

import logging
import time
import tkinter as tk
from concurrent.futures import Future, ThreadPoolExecutor
from tkinter import messagebox

import numpy as np
import pyvista as pv
from matplotlib.colors import LinearSegmentedColormap
from numpy.typing import NDArray
from pyvistaqt import BackgroundPlotter

from ._config_module import Config
from ._plotter_ui import _OrbitalSelectionScreen, _PlotterUI
from ._plotting_objects import Molecule
from .tabulator import GridType, Tabulator


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


class Plotter(_PlotterUI):
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
        self._gto_future: Future[NDArray[np.floating]] | None = None
        self._gtos_ready = only_molecule
        self._gto_start_time: float | None = None
        self._active_gto_job_id: int | None = None
        self._gto_job_counter = 0

        self._tk_root = tk_root
        self._no_prev_tk_root = self._tk_root is None
        if self._tk_root is None:
            self._tk_root = tk.Tk()
            self._tk_root.withdraw()  # Hides window
            logger.debug('Created internal Tk root window for Plotter UI.')

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

        self._molecule: Molecule
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
        if config.mo.custom_colors is not None:
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
        if self._gto_future is None:
            raise RuntimeError('GTO tabulation has not been scheduled.')
        gtos = self._gto_future.result(timeout=timeout)
        if not self._gtos_ready:
            self._apply_gtos_ready(gtos)

    @staticmethod
    def _custom_cmap_from_colors(colors: list[str]) -> LinearSegmentedColormap:
        """Create a custom colormap from a list of colors.

        Parameters
        ----------
        colors : list[str]
            List of color names or hex codes. Must contain exactly two colors.

        Returns
        -------
        LinearSegmentedColormap
            The resulting custom colormap.
        """
        return LinearSegmentedColormap.from_list('custom_mo', colors)

    def _schedule_gto_tabulation(self) -> None:
        """Submit background GTO tabulation work."""
        if self._only_molecule or self._gtos_ready or self._gto_future is not None:
            return
        logger.info('Starting background GTO tabulation...')
        self._gto_start_time = time.perf_counter()
        self._gto_job_counter += 1
        job_id = self._gto_job_counter
        self._active_gto_job_id = job_id
        self._gto_future = _GTO_EXECUTOR.submit(self.tabulator.tabulate_gtos)
        self._gto_future.add_done_callback(lambda fut, job_id=job_id: self._on_gtos_ready(fut, job_id))

    def _on_gtos_ready(self, future: Future[NDArray[np.floating]], job_id: int) -> None:
        """Handle completion of a background tabulation future."""

        def _finish_on_main_thread() -> None:
            if self._active_gto_job_id != job_id:
                return
            self._active_gto_job_id = None
            if self._gtos_ready:
                self._gto_future = None
                return
            if future.cancelled():
                logger.info('Background GTO tabulation cancelled.')
                self._gto_future = None
                self._gto_start_time = None
                return
            try:
                gtos = future.result()
            except Exception as exc:  # pragma: no cover - defensive logging
                self._gto_future = None
                self._gto_start_time = None
                logger.exception('Background GTO tabulation failed.')
                messagebox.showerror('Orbital Tabulation Failed', f'Failed to tabulate orbitals:\n\n{exc!s}')
                return
            self._apply_gtos_ready(gtos)

        if self._tk_root is not None:
            self._tk_root.after_idle(_finish_on_main_thread)
        else:
            _finish_on_main_thread()

    def _apply_gtos_ready(self, gtos: NDArray[np.floating]) -> None:
        """Store computed GTOs and update UI state."""
        self.tabulator._gtos = gtos  # ruff:ignore[private-member-access]
        self._gtos_ready = True
        self._gto_future = None
        self._active_gto_job_id = None
        if self._gto_start_time is not None:
            elapsed = time.perf_counter() - self._gto_start_time
            logger.info('GTO tabulation completed in %.2fs.', elapsed)
            self._gto_start_time = None
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
        if self._gto_future is None:
            return
        if not self._gto_future.done():
            logger.info('Cancelling pending GTO tabulation job.')
            self._gto_future.cancel()
        self._gto_future = None
        self._gto_start_time = None
        self._active_gto_job_id = None

    def _load_molecule(self, config: Config) -> None:
        """Reload the molecule from the parser data."""
        self._molecule = Molecule(self.tabulator._parser.atoms, config)  # ruff:ignore[private-member-access]
        logger.info('Loaded molecule with %d atoms.', len(self._molecule.atoms))

        for actor in self._molecule_actors if hasattr(self, '_molecule_actors') else []:
            self._pv_plotter.remove_actor(actor)

        self._molecule_actors, self._atom_actors, self._bond_actors = self._molecule._add_meshes(  # ruff:ignore[private-member-access]
            self._pv_plotter,
            self._molecule_opacity,
        )
        logger.debug('Added %d molecule actors to the scene.', len(self._molecule_actors))

    def plot_orbital(self, orb_ind: int) -> None:
        """Render the selected orbital isosurface in the PyVista plotter."""
        if not self._ensure_gtos_ready():
            return
        if self._orb_actor:
            self._pv_plotter.remove_actor(self._orb_actor)
            self._orb_actor = None
        if self._selection_screen:
            self._selection_screen.current_mo_ind = orb_ind

        if orb_ind == -1:
            if self._selection_screen:
                self._selection_screen._update_nav_button_states()  # ruff:ignore[private-member-access]
            logger.info('Clearing molecular orbital from scene.')
            return

        mo = self.tabulator._parser.mos[orb_ind]  # ruff:ignore[private-member-access]
        logger.info(
            'Displaying molecular orbital #%d (%s, spin=%s, occ=%s, energy=%.6f au).',
            orb_ind + 1,
            mo.sym,
            mo.spin,
            mo.occ,
            mo.energy,
        )

        self._orb_mesh['orbital'] = self.tabulator.tabulate_mos(orb_ind)

        contour_mesh = self._orb_mesh.contour([-self._contour, self._contour])

        self._orb_actor = self._pv_plotter.add_mesh(
            contour_mesh,
            clim=[-self._contour, self._contour],
            opacity=self._opacity,
            show_scalar_bar=False,
            cmap=self._cmap,
            smooth_shading=True,
        )
        if self._selection_screen:
            self._selection_screen._update_nav_button_states()  # ruff:ignore[private-member-access]

    def _connect_pv_plotter_close_signal(self) -> None:
        """Connect the PyVista plotter close signal to handle closing both windows."""

        def on_pv_plotter_close() -> None:
            """Handle PyVista plotter close event by closing the selection screen and quitting."""
            if self._on_screen:
                self._on_screen = False
                self._cancel_gto_future()
                if self._selection_screen and self._selection_screen.winfo_exists():
                    self._selection_screen.destroy()
                if self._tk_root and self._no_prev_tk_root:
                    self._tk_root.quit()

        self._pv_plotter.app_window.signal_close.connect(on_pv_plotter_close)

    def _clear_all(self) -> None:
        """Clear all actors from the plotter, including molecule and orbitals."""
        if self._molecule_actors:
            for actor in self._molecule_actors:
                actor.SetVisibility(False)

        if self._orb_actor:
            self._pv_plotter.remove_actor(self._orb_actor)
            self._orb_actor = None
            if self._selection_screen:
                self._selection_screen.current_mo_ind = -1
                self._selection_screen._update_nav_button_states()  # ruff:ignore[private-member-access]

    def toggle_molecule(self) -> None:
        """Toggle the visibility of the molecule."""
        if not self._molecule_actors:
            return

        if self.are_bonds_visible() != self.are_atoms_visible():
            if self.are_bonds_visible():
                self.toggle_atoms()
            else:
                self.toggle_bonds()
        else:
            for actor in self._molecule_actors:
                actor.SetVisibility(not actor.GetVisibility())
            self._pv_plotter.update()

        self._update_settings_button_states()

    def toggle_atoms(self) -> None:
        """Toggle the visibility of the molecule."""
        if self._atom_actors:
            for actor in self._atom_actors:
                actor.SetVisibility(not actor.GetVisibility())
            self._pv_plotter.update()

    def toggle_bonds(self) -> None:
        """Toggle the visibility of the molecule."""
        if self._bond_actors:
            for actor in self._bond_actors:
                actor.SetVisibility(not actor.GetVisibility())
            self._pv_plotter.update()

    def is_molecule_visible(self) -> bool:
        """Check if the molecule is currently visible in the plotter.

        Returns
        -------
        bool
            `True` if the molecule is visible, `False` otherwise.
        """
        if self._molecule_actors:
            return bool(self._molecule_actors[0].GetVisibility())  # Check visibility of the first actor
        return False

    def are_atoms_visible(self) -> bool:
        """Check if the atoms are currently visible in the plotter.

        Returns
        -------
        bool
            `True` if the atoms are visible, `False` otherwise.
        """
        if self._atom_actors:
            return bool(self._atom_actors[0].GetVisibility())  # Check visibility of the first actor
        return False

    def are_bonds_visible(self) -> bool:
        """Check if the bonds are currently visible in the plotter.

        Returns
        -------
        bool
            `True` if the bonds are visible, `False` otherwise.
        """
        if self._bond_actors:
            return bool(self._bond_actors[0].GetVisibility())  # Check visibility of the first actor
        return False

    def _create_mo_mesh(self) -> pv.StructuredGrid:
        """Create a mesh for the orbitals.

        Returns
        -------
            pv.StructuredGrid:
                The mesh object for MO visualization.

        """
        mesh = pv.StructuredGrid()
        mesh.points = pv.pyvista_ndarray(self.tabulator.grid)  # pyright: ignore[reportCallIssue]

        # Pyvista needs the dimensions backwards
        # in other words, (phi, theta, r) or (z, y, x)
        mesh.dimensions = self.tabulator.grid_dimensions[::-1]

        return mesh

    def _update_mesh(
        self,
        i_points: NDArray[np.floating],
        j_points: NDArray[np.floating],
        k_points: NDArray[np.floating],
        grid_type: GridType,
    ) -> None:
        """Update the tabulator grid and rebuild the orbital mesh.

        Parameters
        ----------
        i_points : NDArray[np.floating]
            1D array defining the first dimension (radius or x).
        j_points : NDArray[np.floating]
            1D array defining the second dimension (theta or y).
        k_points : NDArray[np.floating]
            1D array defining the third dimension (phi or z).
        grid_type : GridType
            Target grid type to regenerate (`GridType.SPHERICAL` or
            `GridType.CARTESIAN`).

        Raises
        ------
        ValueError
            If ``grid_type`` is not supported.
        """
        self._cancel_gto_future()
        self._gtos_ready = False
        if self._selection_screen:
            self._selection_screen._set_loading_state(  # ruff:ignore[private-member-access]
                True,
                'Updating grid...',
            )
        if grid_type == GridType.CARTESIAN:
            self.tabulator.cartesian_grid(i_points, j_points, k_points)
        elif grid_type == GridType.SPHERICAL:
            self.tabulator.spherical_grid(i_points, j_points, k_points)
        else:
            raise ValueError('The plotter only supports spherical and cartesian grids.')

        dimensions = 'x'.join(str(val) for val in self.tabulator.grid_dimensions)
        logger.info(
            'Rebuilt %s grid with %d points (dimensions %s).',
            grid_type.value,
            self.tabulator.grid.shape[0],
            dimensions,
        )

        self._orb_mesh = self._create_mo_mesh()
        self._gtos_ready = True
        if self._selection_screen:
            self._selection_screen._on_gtos_ready()  # ruff:ignore[private-member-access]
            if self._selection_screen.current_mo_ind >= 0:
                self.plot_orbital(self._selection_screen.current_mo_ind)
