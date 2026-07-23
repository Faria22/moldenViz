"""PyVista scene and molecular-orbital rendering for :class:`Plotter`."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import pyvista as pv
from matplotlib.colors import LinearSegmentedColormap

from ._plotting_objects import Molecule
from .tabulator import GridType

if TYPE_CHECKING:
    import tkinter as tk

    import numpy as np
    from numpy.typing import NDArray

    from ._config_module import Config
    from .tabulator import Tabulator

logger = logging.getLogger(__name__)


class _PlotterRendering:
    """Mixin responsible for PyVista scene and orbital rendering."""

    if TYPE_CHECKING:
        _atom_actors: list[Any]
        _bond_actors: list[Any]
        _cmap: Any
        _contour: float
        _gtos_ready: bool
        _molecule: Molecule
        _molecule_actors: list[Any]
        _molecule_opacity: float
        _no_prev_tk_root: bool
        _on_screen: bool
        _opacity: float
        _orb_actor: Any | None
        _orb_mesh: pv.StructuredGrid
        _only_molecule: bool
        _pv_plotter: Any
        _selection_screen: Any | None
        _tk_root: tk.Misc | None
        tabulator: Tabulator

        def _cancel_gto_future(self) -> None: ...
        def _ensure_gtos_ready(self) -> bool: ...
        def _update_settings_button_states(self) -> None: ...

    @staticmethod
    def _custom_cmap_from_colors(colors: list[str]) -> LinearSegmentedColormap:
        """Create a custom colormap from two endpoint colors.

        Returns
        -------
        LinearSegmentedColormap
            Colormap interpolating between the supplied colors.
        """
        return LinearSegmentedColormap.from_list('custom_mo', colors)

    def _load_molecule(self, current_config: Config) -> None:
        """Reload the molecule from parsed atom data."""
        self._molecule = Molecule(self.tabulator.atoms, current_config)
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

        mo = self.tabulator.molecular_orbitals[orb_ind]
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
        """Connect the PyVista close signal to the Plotter lifecycle."""

        def on_pv_plotter_close() -> None:
            if self._on_screen:
                self._on_screen = False
                self._cancel_gto_future()
                if self._selection_screen and self._selection_screen.winfo_exists():
                    self._selection_screen.destroy()
                if self._tk_root and self._no_prev_tk_root:
                    self._tk_root.quit()

        self._pv_plotter.app_window.signal_close.connect(on_pv_plotter_close)

    def _clear_all(self) -> None:
        """Clear all molecule and orbital actors."""
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
        """Toggle visibility of all molecule actors."""
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
        """Toggle atom visibility."""
        if self._atom_actors:
            for actor in self._atom_actors:
                actor.SetVisibility(not actor.GetVisibility())
            self._pv_plotter.update()

    def toggle_bonds(self) -> None:
        """Toggle bond visibility."""
        if self._bond_actors:
            for actor in self._bond_actors:
                actor.SetVisibility(not actor.GetVisibility())
            self._pv_plotter.update()

    def is_molecule_visible(self) -> bool:
        """Return whether the molecule is visible.

        Returns
        -------
        bool
            Whether the first molecule actor is visible.
        """
        if self._molecule_actors:
            return bool(self._molecule_actors[0].GetVisibility())
        return False

    def are_atoms_visible(self) -> bool:
        """Return whether atoms are visible.

        Returns
        -------
        bool
            Whether the first atom actor is visible.
        """
        if self._atom_actors:
            return bool(self._atom_actors[0].GetVisibility())
        return False

    def are_bonds_visible(self) -> bool:
        """Return whether bonds are visible.

        Returns
        -------
        bool
            Whether the first bond actor is visible.
        """
        if self._bond_actors:
            return bool(self._bond_actors[0].GetVisibility())
        return False

    def _create_mo_mesh(self) -> pv.StructuredGrid:
        """Create a structured mesh for orbital scalar data.

        Returns
        -------
        pv.StructuredGrid
            Mesh configured from the current Tabulator grid.
        """
        mesh = pv.StructuredGrid()
        mesh.points = pv.pyvista_ndarray(self.tabulator.grid)  # pyright: ignore[reportCallIssue]
        mesh.dimensions = self.tabulator.grid_dimensions[::-1]
        return mesh

    def _update_mesh(
        self,
        i_points: NDArray[np.floating],
        j_points: NDArray[np.floating],
        k_points: NDArray[np.floating],
        grid_type: GridType,
    ) -> None:
        """Update the Tabulator grid and rebuild the orbital mesh."""
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
