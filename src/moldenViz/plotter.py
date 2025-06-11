"""Plotter module for creating plots of the molecule and it's orbitals."""

import logging
import tkinter as tk
from tkinter import messagebox, ttk
from typing import Optional

import numpy as np
import pyvista as pv
from numpy.typing import NDArray
from pyvistaqt import BackgroundPlotter

from ._plotting_objects import Molecule
from .parser import _MolecularOrbital
from .tabulator import GridType, Tabulator, _cartesian_to_spherical, _spherical_to_cartesian

logger = logging.getLogger(__name__)


class Plotter:
    """
    Handles the 3D visualization of molecules and molecular orbitals.

    This class uses PyVista for 3D rendering and Tkinter for the user interface
    to control plotting parameters and select orbitals.

    Args
    ----
        source: str | list[str]
            The path to the molden file, or the lines from the file.

        only_molecule: bool, optional
            Only parse the atoms and skip molecular orbitals.
            Default is `False`.

        tabulator: Tabulator, optional
            If `None`, `Plotter` creates a `Tabulator` and tabulates the GTOs and MOs with a default grid.
            A `Tabulator` can be passed as to tabulate the GTOs in a predetermined grid.

            Note: `Tabulator` grid must be spherical or cartesian. Custom grids are not allowed.

        tk_rook: tk.Tk, optional
            If user is using the plotter inside a tk app, `tk_root` can be passed as to not create a new tk instance.

    Raises
    ------
        ValueError:
            If the provided tabulator is invalid
            (e.g., missing grid or GTO data when `only_molecule` is `False`, or has an UNKNOWN grid type).
    """

    # Default values
    CONTOUR = 0.1
    OPACITY = 1.0
    MOLECULE_OPACITY = 1.0

    MIN_RADIUS = 5.0

    NUM_RADIUS_POINTS = 100
    NUM_THETA_POINTS = 60
    NUM_PHI_POINTS = 120

    NUM_X_POINTS = 100
    NUM_Y_POINTS = 100
    NUM_Z_POINTS = 100

    def __init__(
        self,
        source: str | list[str],
        only_molecule: bool = False,
        tabulator: Optional[Tabulator] = None,
        tk_root: Optional[tk.Tk] = None,
    ) -> None:
        self.on_screen = True

        if tabulator:
            if not hasattr(tabulator, 'grid'):
                raise ValueError('Tabulator does not have grid attribute.')

            if not hasattr(tabulator, 'gto_data') and not only_molecule:
                raise ValueError('Tabulator does not have tabulated GTOs.')

            if tabulator.grid_type == GridType.UNKNOWN:
                raise ValueError('The plotter only supports spherical and cartesian grids.')

            self.tab = tabulator
        else:
            self.tab = Tabulator(source, only_molecule=only_molecule)

        self.molecule = Molecule(self.tab.atoms)
        self.molecule_opacity = self.MOLECULE_OPACITY

        if not only_molecule:
            self.tk_root = tk_root
            self.no_prev_root = self.tk_root is None
            if self.no_prev_root:
                self.tk_root = tk.Tk()
                self.tk_root.withdraw()  # Hides window

        self.pv_plotter = BackgroundPlotter(editor=False)
        self.pv_plotter.show_axes()
        self.molecule_actors = self.molecule.add_meshes(self.pv_plotter, self.molecule_opacity)

        # If we want to have the molecular orbitals, we need to initiate Tk before Qt
        # That is why we have this weird if statement separated this way
        if only_molecule:
            self.pv_plotter.app.exec_()  # pyright: ignore[reportAttributeAccessIssue]
            return

        assert self.tk_root is not None  # To help type hinters

        if not tabulator:
            # Default is a spherical grid
            self.tab.spherical_grid(
                np.linspace(0, min(2 * self.molecule.max_radius, self.MIN_RADIUS), self.NUM_RADIUS_POINTS),
                np.linspace(0, np.pi, self.NUM_THETA_POINTS),
                np.linspace(0, 2 * np.pi, self.NUM_PHI_POINTS),
            )

        self.orb_mesh = self._create_mo_mesh()
        self.orb_actor: pv.Actor | None = None

        # Values for MO, not the molecule
        self.contour = self.CONTOUR
        self.opacity = self.OPACITY

        _OrbitalSelectionScreen(self, self.tk_root)

        self.tk_root.mainloop()

    def toggle_molecule(self) -> None:
        """Toggle the visibility of the molecule."""
        if self.molecule_actors:
            for actor in self.molecule_actors:
                actor.SetVisibility(not actor.GetVisibility())
            self.pv_plotter.update()

    def _create_mo_mesh(self) -> pv.StructuredGrid:
        """Create a mesh for the orbitals.

        Returns
        -------
            pv.StructuredGrid:
                The mesh object for MO visualization.

        """
        mesh = pv.StructuredGrid()
        mesh.points = pv.pyvista_ndarray(self.tab.grid)

        # Pyvista needs the dimensions backwards
        # in other words, (phi, theta, r) or (z, y, x)
        mesh.dimensions = self.tab.grid_dimensions[::-1]

        return mesh

    def update_mesh(
        self,
        i_points: NDArray[np.floating],
        j_points: NDArray[np.floating],
        k_points: NDArray[np.floating],
        grid_type: GridType,
    ) -> None:
        """Update the grid in the Tabulator and recreates the MO mesh.

        Args:
            i_points: NDArray[np.floating]:
                Points for the first dimension (radius or x).
            j_points: NDArray[np.floating]:
                Points for the second dimension (theta or y).
            k_points: NDArray[np.floating]:
                Points for the third dimension (phi or z).
            grid_type: GridType:
                The type of grid to create (SPHERICAL or CARTESIAN).

        Raises
        ------
            ValueError: If the grid_type is not SPHERICAL or CARTESIAN.
        """
        if grid_type == GridType.CARTESIAN:
            self.tab.cartesian_grid(i_points, j_points, k_points)
        elif grid_type == GridType.SPHERICAL:
            self.tab.spherical_grid(i_points, j_points, k_points)
        else:
            raise ValueError('The plotter only supports spherical and cartesian grids.')

        self.orb_mesh = self._create_mo_mesh()


class _OrbitalSelectionScreen(tk.Toplevel):
    SPHERICAL_GRID_SETTINGS_WINDOW_SIZE = '600x500'
    CARTESIAN_GRID_SETTINGS_WINDOW_SIZE = '800x500'

    def __init__(self, plotter: Plotter, tk_master: tk.Tk) -> None:
        super().__init__(tk_master)
        self.title('Orbitals')
        self.geometry('350x500')

        self.protocols()

        self.plotter = plotter
        self.current_orb_ind = -1  # Start with no orbital shown

        nav_frame = ttk.Frame(self)
        nav_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.prev_button = ttk.Button(nav_frame, text='<< Previous', command=self.prev_plot)
        self.prev_button.pack(side=tk.LEFT, padx=5, pady=10)

        self.next_button = ttk.Button(nav_frame, text='Next >>', command=self.next_plot)
        self.next_button.pack(side=tk.RIGHT, padx=5, pady=10)

        self.update_button_states()  # Update buttons for initial state

        self.settings_button = ttk.Button(self, text='Settings', command=self.settings_screen)
        self.settings_button.pack(expand=True, padx=5, pady=10)

        self.orb_tv = _OrbitalsTreeview(self)
        self.orb_tv.populate_tree(self.plotter.tab.mos)
        self.orb_tv.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    def on_close(self) -> None:
        self.plotter.on_screen = False
        self.plotter.pv_plotter.close()
        self.destroy()
        if self.plotter.tk_root and self.plotter.no_prev_root:
            self.plotter.tk_root.destroy()

    def protocols(self) -> None:
        self.protocol('WM_DELETE_WINDOW', self.on_close)
        self.bind('<Command-q>', lambda _event: self.on_close())
        self.bind('<Command-w>', lambda _event: self.on_close())
        self.bind('<Control-q>', lambda _event: self.on_close())
        self.bind('<Control-w>', lambda _event: self.on_close())

    def settings_screen(self) -> None:
        self.settings_window = tk.Toplevel(self)
        self.settings_window.title('Settings')

        settings_frame = ttk.Frame(self.settings_window)
        settings_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Contour level
        contour_label = ttk.Label(settings_frame, text='Molecular Orbital Contour:')
        contour_label.grid(row=0, column=0, padx=5, pady=5)
        self.contour_entry = ttk.Entry(settings_frame)
        self.contour_entry.insert(0, str(self.plotter.contour))
        self.contour_entry.grid(row=1, column=0, padx=5, pady=5)

        # Opacity
        opacity_label = ttk.Label(settings_frame)
        opacity_label.grid(row=2, column=0, padx=5, pady=5)
        self.opacity_scale = ttk.Scale(
            settings_frame,
            length=150,
            command=lambda val: opacity_label.config(text=f'Molecular Orbital Opacity: {float(val):.2f}'),
        )
        self.opacity_scale.set(self.plotter.opacity)
        self.opacity_scale.grid(row=3, column=0, padx=5, pady=5)

        # Molecule Opacity
        molecule_opacity_label = ttk.Label(settings_frame)
        molecule_opacity_label.grid(row=4, column=0, padx=5, pady=5)
        self.molecule_opacity_scale = ttk.Scale(
            settings_frame,
            length=150,
            command=lambda val: molecule_opacity_label.config(text=f'Molecule Opacity: {float(val):.2f}'),
        )
        self.molecule_opacity_scale.set(self.plotter.molecule_opacity)
        self.molecule_opacity_scale.grid(row=5, column=0, padx=5, pady=5)

        # Toggle molecule visibility
        toggle_mol_button = ttk.Button(settings_frame, text='Toggle Molecule', command=self.plotter.toggle_molecule)
        toggle_mol_button.grid(row=6, column=0, padx=5, pady=5)

        # Grid parameters
        ttk.Label(settings_frame, text='MO Grid parameters').grid(row=0, column=1, padx=5, pady=5, columnspan=4)

        self.grid_type_radio_var = tk.StringVar()
        self.grid_type_radio_var.set(self.plotter.tab.grid_type.value)

        ttk.Label(settings_frame, text='Spherical grid:').grid(row=1, column=1, padx=5, pady=5)
        sph_grid_type_button = ttk.Radiobutton(
            settings_frame,
            variable=self.grid_type_radio_var,
            value=GridType.SPHERICAL.value,
            command=self.place_grid_params_frame,
        )

        ttk.Label(settings_frame, text='Cartesian grid:').grid(row=1, column=3, padx=5, pady=5)
        cart_grid_type_button = ttk.Radiobutton(
            settings_frame,
            variable=self.grid_type_radio_var,
            value=GridType.CARTESIAN.value,
            command=self.place_grid_params_frame,
        )

        sph_grid_type_button.grid(row=1, column=2, padx=5, pady=5)
        cart_grid_type_button.grid(row=1, column=4, padx=5, pady=5)

        self.sph_grid_params_frame = self.sph_grid_params_frame_widgets(settings_frame)
        self.cart_grid_params_frame = self.cart_grid_params_frame_widgets(settings_frame)

        self.place_grid_params_frame()

        # Reset button
        reset_button = ttk.Button(settings_frame, text='Reset', command=self.reset_settings)
        reset_button.grid(row=8, column=0, padx=5, pady=5, columnspan=5)

        # Save settings button
        save_button = ttk.Button(settings_frame, text='Apply', command=self.apply_settings)
        save_button.grid(row=9, column=0, padx=5, pady=5, columnspan=5)

    def place_grid_params_frame(self) -> None:
        if self.grid_type_radio_var.get() == GridType.SPHERICAL.value:
            self.settings_window.geometry(self.SPHERICAL_GRID_SETTINGS_WINDOW_SIZE)
            self.cart_grid_params_frame.grid_forget()
            self.settings_window.geometry()
            self.sph_grid_params_frame.grid(row=2, column=1, padx=5, pady=5, rowspan=6, columnspan=4)
            self.sph_grid_params_frame_setup()
        else:
            self.settings_window.geometry(self.CARTESIAN_GRID_SETTINGS_WINDOW_SIZE)
            self.sph_grid_params_frame.grid_forget()
            self.cart_grid_params_frame.grid(row=2, column=1, padx=5, pady=5, rowspan=6, columnspan=4)
            self.cart_grid_params_frame_setup()

    def sph_grid_params_frame_widgets(self, master: ttk.Frame) -> ttk.Frame:
        grid_params_frame = ttk.Frame(master)

        # Radius
        ttk.Label(grid_params_frame, text='Radius:').grid(row=0, column=0, padx=5, pady=5)
        self.radius_entry = ttk.Entry(grid_params_frame)
        self.radius_entry.grid(row=0, column=1, padx=5, pady=5)

        # Radius points
        radius_points_label = ttk.Label(grid_params_frame, text='Number of Radius Points:')
        radius_points_label.grid(row=1, column=0, padx=5, pady=5)
        self.radius_points_entry = ttk.Entry(grid_params_frame)
        self.radius_points_entry.grid(row=1, column=1, padx=5, pady=5)

        # Theta points
        theta_points_label = ttk.Label(grid_params_frame, text='Number of Theta Points:')
        theta_points_label.grid(row=2, column=0, padx=5, pady=5)
        self.theta_points_entry = ttk.Entry(grid_params_frame)
        self.theta_points_entry.grid(row=2, column=1, padx=5, pady=5)

        # Phi points
        phi_points_label = ttk.Label(grid_params_frame, text='Number of Phi Points:')
        phi_points_label.grid(row=3, column=0, padx=5, pady=5)
        self.phi_points_entry = ttk.Entry(grid_params_frame)
        self.phi_points_entry.grid(row=3, column=1, padx=5, pady=5)

        return grid_params_frame

    def cart_grid_params_frame_widgets(self, master: ttk.Frame) -> ttk.Frame:
        grid_params_frame = ttk.Frame(master)

        # X
        ttk.Label(grid_params_frame, text='Min x:').grid(row=0, column=0, padx=5, pady=5)
        ttk.Label(grid_params_frame, text='Max x:').grid(row=0, column=1, padx=5, pady=5)
        ttk.Label(grid_params_frame, text='Num x points:').grid(row=0, column=2, padx=5, pady=5)

        self.x_min_entry = ttk.Entry(grid_params_frame)
        self.x_max_entry = ttk.Entry(grid_params_frame)
        self.x_num_points_entry = ttk.Entry(grid_params_frame)

        self.x_min_entry.grid(row=1, column=0, padx=5, pady=5)
        self.x_max_entry.grid(row=1, column=1, padx=5, pady=5)
        self.x_num_points_entry.grid(row=1, column=2, padx=5, pady=5)

        # Y
        ttk.Label(grid_params_frame, text='Min y:').grid(row=2, column=0, padx=5, pady=5)
        ttk.Label(grid_params_frame, text='Max y:').grid(row=2, column=1, padx=5, pady=5)
        ttk.Label(grid_params_frame, text='Num y points:').grid(row=2, column=2, padx=5, pady=5)

        self.y_min_entry = ttk.Entry(grid_params_frame)
        self.y_max_entry = ttk.Entry(grid_params_frame)
        self.y_num_points_entry = ttk.Entry(grid_params_frame)

        self.y_min_entry.grid(row=3, column=0, padx=5, pady=5)
        self.y_max_entry.grid(row=3, column=1, padx=5, pady=5)
        self.y_num_points_entry.grid(row=3, column=2, padx=5, pady=5)

        # Z
        ttk.Label(grid_params_frame, text='Min z:').grid(row=4, column=0, padx=5, pady=5)
        ttk.Label(grid_params_frame, text='Max z:').grid(row=4, column=1, padx=5, pady=5)
        ttk.Label(grid_params_frame, text='Num z points:').grid(row=4, column=2, padx=5, pady=5)

        self.z_min_entry = ttk.Entry(grid_params_frame)
        self.z_max_entry = ttk.Entry(grid_params_frame)
        self.z_num_points_entry = ttk.Entry(grid_params_frame)

        self.z_min_entry.grid(row=5, column=0, padx=5, pady=5)
        self.z_max_entry.grid(row=5, column=1, padx=5, pady=5)
        self.z_num_points_entry.grid(row=5, column=2, padx=5, pady=5)

        return grid_params_frame

    def sph_grid_params_frame_setup(self) -> None:
        self.radius_entry.delete(0, tk.END)
        self.radius_points_entry.delete(0, tk.END)
        self.theta_points_entry.delete(0, tk.END)
        self.phi_points_entry.delete(0, tk.END)

        # Previous grid was cartesian, so use default values
        if self.plotter.tab.grid_type == GridType.CARTESIAN:
            self.radius_entry.insert(0, str(min(self.plotter.molecule.max_radius * 2, self.plotter.MIN_RADIUS)))
            self.radius_points_entry.insert(0, str(self.plotter.NUM_RADIUS_POINTS))
            self.theta_points_entry.insert(0, str(self.plotter.NUM_THETA_POINTS))
            self.phi_points_entry.insert(0, str(self.plotter.NUM_PHI_POINTS))
            return

        num_r, num_theta, num_phi = self.plotter.tab.grid_dimensions

        # The last point of the grid for sure has the largest r
        r, _, _ = _cartesian_to_spherical(*self.plotter.tab.grid[-1, :])  # pyright: ignore[reportArgumentType]

        self.radius_entry.insert(0, str(r))
        self.radius_points_entry.insert(0, str(num_r))
        self.theta_points_entry.insert(0, str(num_theta))
        self.phi_points_entry.insert(0, str(num_phi))

    def cart_grid_params_frame_setup(self) -> None:
        self.x_min_entry.delete(0, tk.END)
        self.x_max_entry.delete(0, tk.END)
        self.x_num_points_entry.delete(0, tk.END)

        self.y_min_entry.delete(0, tk.END)
        self.y_max_entry.delete(0, tk.END)
        self.y_num_points_entry.delete(0, tk.END)

        self.z_min_entry.delete(0, tk.END)
        self.z_max_entry.delete(0, tk.END)
        self.z_num_points_entry.delete(0, tk.END)

        # Previous grid was sphesical, so use adapted default values
        if self.plotter.tab.grid_type == GridType.SPHERICAL:
            r = min(2 * self.plotter.molecule.max_radius, self.plotter.MIN_RADIUS)

            self.x_min_entry.insert(0, str(-r))
            self.y_min_entry.insert(0, str(-r))
            self.z_min_entry.insert(0, str(-r))

            self.x_max_entry.insert(0, str(r))
            self.y_max_entry.insert(0, str(r))
            self.z_max_entry.insert(0, str(r))

            self.x_num_points_entry.insert(0, str(self.plotter.NUM_X_POINTS))
            self.y_num_points_entry.insert(0, str(self.plotter.NUM_Y_POINTS))
            self.z_num_points_entry.insert(0, str(self.plotter.NUM_Z_POINTS))
            return

        x_num, y_num, z_num = self.plotter.tab.grid_dimensions
        x_min, y_min, z_min = self.plotter.tab.grid[0, :]
        x_max, y_max, z_max = self.plotter.tab.grid[-1, :]

        self.x_min_entry.insert(0, str(x_min))
        self.x_max_entry.insert(0, str(x_max))
        self.x_num_points_entry.insert(0, str(x_num))

        self.y_min_entry.insert(0, str(y_min))
        self.y_max_entry.insert(0, str(y_max))
        self.y_num_points_entry.insert(0, str(y_num))

        self.z_min_entry.insert(0, str(z_min))
        self.z_max_entry.insert(0, str(z_max))
        self.z_num_points_entry.insert(0, str(z_num))

    def reset_settings(self) -> None:
        """Reset settings to default values."""
        self.contour_entry.delete(0, tk.END)
        self.contour_entry.insert(0, str(self.plotter.CONTOUR))

        self.opacity_scale.set(Plotter.OPACITY)

        self.molecule_opacity_scale.set(self.plotter.MOLECULE_OPACITY)

        self.grid_type_radio_var.set(GridType.SPHERICAL.value)

        self.radius_entry.delete(0, tk.END)
        self.radius_entry.insert(0, str(min(self.plotter.molecule.max_radius * 2, self.plotter.MIN_RADIUS)))

        self.radius_points_entry.delete(0, tk.END)
        self.radius_points_entry.insert(0, str(self.plotter.NUM_RADIUS_POINTS))

        self.theta_points_entry.delete(0, tk.END)
        self.theta_points_entry.insert(0, str(self.plotter.NUM_THETA_POINTS))

        self.phi_points_entry.delete(0, tk.END)
        self.phi_points_entry.insert(0, str(self.plotter.NUM_PHI_POINTS))

    def apply_settings(self) -> None:
        self.plotter.molecule_opacity = round(self.molecule_opacity_scale.get(), 2)
        for actor in self.plotter.molecule_actors:
            actor.GetProperty().SetOpacity(self.plotter.molecule_opacity)

        if self.grid_type_radio_var.get() == GridType.SPHERICAL.value:
            radius = float(self.radius_entry.get())
            if radius <= 0:
                messagebox.showerror('Invalid input', 'Radius must be greater than zero.')
                return

            num_r_points = int(self.radius_points_entry.get())
            num_theta_points = int(self.theta_points_entry.get())
            num_phi_points = int(self.phi_points_entry.get())

            if num_r_points <= 0 or num_theta_points <= 0 or num_phi_points <= 0:
                messagebox.showerror('Invalid input', 'Number of points must be greater than zero.')
                return

            r = np.linspace(0, radius, num_r_points)
            theta = np.linspace(0, np.pi, num_theta_points)
            phi = np.linspace(0, 2 * np.pi, num_phi_points)

            rr, tt, pp = np.meshgrid(r, theta, phi, indexing='ij')
            xx, yy, zz = _spherical_to_cartesian(rr, tt, pp)

            # Update the mesh with new points if needed
            new_grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
            if not np.array_equal(new_grid, self.plotter.tab.grid):
                self.plotter.update_mesh(r, theta, phi, GridType.SPHERICAL)

        else:
            x_min = float(self.x_min_entry.get())
            x_max = float(self.x_max_entry.get())
            x_num = int(self.x_num_points_entry.get())

            y_min = float(self.y_min_entry.get())
            y_max = float(self.y_max_entry.get())
            y_num = int(self.y_num_points_entry.get())

            z_min = float(self.z_min_entry.get())
            z_max = float(self.z_max_entry.get())
            z_num = int(self.z_num_points_entry.get())

            if x_num <= 0 or y_num <= 0 or z_num <= 0:
                messagebox.showerror('Invalid input', 'Number of points must be greater than zero.')
                return

            x = np.linspace(x_min, x_max, x_num)
            y = np.linspace(y_min, y_max, y_num)
            z = np.linspace(z_min, z_max, z_num)

            xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')

            # Update the mesh with new points if needed
            new_grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
            if not np.array_equal(new_grid, self.plotter.tab.grid):
                self.plotter.update_mesh(x, y, z, GridType.CARTESIAN)

        self.plotter.contour = float(self.contour_entry.get().strip())
        self.plot_orbital(self.current_orb_ind)

        self.plotter.opacity = round(self.opacity_scale.get(), 2)
        if self.plotter.orb_actor:
            self.plotter.orb_actor.GetProperty().SetOpacity(self.plotter.opacity)

    def next_plot(self) -> None:
        """Go to the next plot."""
        self.current_orb_ind += 1
        self.update_button_states()
        self.orb_tv.highlight_orbital(self.current_orb_ind)
        self.plot_orbital(self.current_orb_ind)

    def prev_plot(self) -> None:
        """Go to the previous plot."""
        self.current_orb_ind -= 1
        self.orb_tv.highlight_orbital(self.current_orb_ind)
        self.update_button_states()
        self.plot_orbital(self.current_orb_ind)

    def update_button_states(self) -> None:
        """Update the enabled/disabled state of nav buttons."""
        can_go_prev = self.current_orb_ind > 0
        can_go_next = self.current_orb_ind < len(self.plotter.tab.mos) - 1
        self.prev_button.config(state=tk.NORMAL if can_go_prev else tk.DISABLED)
        self.next_button.config(state=tk.NORMAL if can_go_next else tk.DISABLED)

    def plot_orbital(self, orb_ind: int) -> None:
        if self.plotter.orb_actor:
            self.plotter.pv_plotter.remove_actor(self.plotter.orb_actor)

        if orb_ind != -1:
            self.plotter.orb_mesh['orbital'] = self.plotter.tab.tabulate_mos(orb_ind)

            contour_mesh = self.plotter.orb_mesh.contour([-self.plotter.contour, self.plotter.contour])

            self.plotter.orb_actor = self.plotter.pv_plotter.add_mesh(
                contour_mesh,
                clim=[-self.plotter.contour, self.plotter.contour],
                opacity=self.plotter.opacity,
                show_scalar_bar=False,
                cmap='bwr',
            )


class _OrbitalsTreeview(ttk.Treeview):
    def __init__(self, selection_screen: _OrbitalSelectionScreen) -> None:
        columns = ['Index', 'Symmetry', 'Occupation', 'Energy [au]']
        widths = [20, 50, 50, 120]

        super().__init__(selection_screen, columns=columns, show='headings', height=20)

        for col, w in zip(columns, widths):
            self.heading(col, text=col)
            self.column(col, width=w)

        self.selection_screen = selection_screen

        self.current_orb_ind = -1  # Start with no orbital shown

        # Configure tag
        self.tag_configure('highlight', background='lightblue')

        self.bind('<<TreeviewSelect>>', self.on_select)

    def highlight_orbital(self, orb_ind: int) -> None:
        """Highlight the selected orbital."""
        if self.current_orb_ind != -1:
            self.item(self.current_orb_ind, tags=('!hightlight',))

        self.current_orb_ind = orb_ind
        self.item(orb_ind, tags=('highlight',))
        self.see(orb_ind)  # Scroll to the selected item

    def erase(self) -> None:
        for item in self.get_children():
            self.delete(item)

    def populate_tree(self, mos: list[_MolecularOrbital]) -> None:
        self.erase()

        # Counts the number of MOs with a given symmetry
        mo_syms = list({mo.sym for mo in mos})
        mo_sym_count: dict[str, int] = dict.fromkeys(mo_syms, 0)
        for ind, mo in enumerate(mos):
            mo_sym_count[mo.sym] += 1
            self.insert('', 'end', iid=ind, values=(ind + 1, f'{mo.sym}.{mo_sym_count[mo.sym]}', mo.occ, mo.energy))

    def on_select(self, _event: tk.Event) -> None:
        """Handle selection of an orbital."""
        selected_item = self.selection()
        self.selection_remove(selected_item)
        if selected_item:
            orb_ind = int(selected_item[0])
            self.highlight_orbital(orb_ind)
            self.selection_screen.current_orb_ind = orb_ind
            self.selection_screen.plot_orbital(orb_ind)
            self.selection_screen.update_button_states()
