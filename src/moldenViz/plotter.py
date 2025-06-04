import logging
import tkinter as tk
from tkinter import ttk
from typing import Optional

import numpy as np
import pyvista as pv
from pyvistaqt import BackgroundPlotter

# from ._plotting_objects import Molecule
from .tabulator import Tabulator

logger = logging.getLogger(__name__)


class Plotter(tk.Toplevel):
    CONTOUR = 0.1
    OPACITY = 1
    MOLECULE_OPACITY = 1
    NUM_RADIUS_POINTS = 100
    NUM_THETA_POINTS = 60
    NUM_PHI_POINTS = 120

    def __init__(
        self,
        molden_file: Optional[str],
        molden_lines: Optional[list[str]],
        tabulator: Optional[Tabulator],
        tk_master: Optional[tk.Tk],
        only_molecule: bool = False,
    ) -> None:
        super().__init__(tk_master)
        self.title('Orbitals')
        self.geometry('350x500')

        self.showing = True

        if tabulator:
            self.tab = tabulator
        else:
            self.tab = Tabulator(filename=molden_file, molden_lines=molden_lines)

        self.molecule = molecule
        self.molecule_opacity = self.MOLECULE_OPACITY

        # Default values
        self.r_points = np.linspace(
            0,
            self.molecule.max_radius * 2,
            self.NUM_RADIUS_POINTS,
        )
        self.theta_points = np.linspace(0, np.pi, self.NUM_THETA_POINTS)
        self.phi_points = np.linspace(0, 2 * np.pi, self.NUM_PHI_POINTS)

        self.tab_obj = TargettClass(
            molden_lines,
            self.r_points,
            self.theta_points,
            self.phi_points,
        )
        self.orb_mesh = self.create_mesh(
            self.r_points,
            self.theta_points,
            self.phi_points,
        )

        self.orb_data = orb_data

        self.plotter = BackgroundPlotter(editor=False)
        self.plotter.show_axes()
        self.molecule_actors = self.molecule.add_meshes(
            self.plotter,
            self.molecule_opacity,
        )

        self.protocol('WM_DELETE_WINDOW', self.on_close)
        self.bind('<Control-q>', lambda _event: self.on_close())
        self.bind('<Command-w>', lambda _event: self.on_close())

        self.contour = self.CONTOUR
        self.opacity = self.OPACITY

        self.orbital_selection_screen()
        self.orb_actor = None
        self.plot_orbital(-1)

    def toggle_molecule(self) -> None:
        """Toggle the visibility of the molecule."""
        if self.molecule_actors:
            for actor in self.molecule_actors:
                actor.SetVisibility(not actor.GetVisibility())
            self.plotter.update()

    def on_close(self) -> None:
        self.showing = False
        self.plotter.close()
        self.destroy()

    def orbital_selection_screen(self) -> None:
        self.current_orb_ind = -1  # Start with no orbital shown

        nav_frame = ttk.Frame(self)
        nav_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.prev_button = ttk.Button(
            nav_frame,
            text='<< Previous',
            command=self.prev_plot,
        )
        self.prev_button.pack(side=tk.LEFT, padx=5, pady=10)

        self.next_button = ttk.Button(nav_frame, text='Next >>', command=self.next_plot)
        self.next_button.pack(side=tk.RIGHT, padx=5, pady=10)

        self.update_button_states()  # Update buttons for initial state

        self.settings_button = ttk.Button(self, text='Settings', command=self.settings_screen)
        self.settings_button.pack(expand=True, padx=5, pady=10)

        self.orb_tv = _OrbitsTreeview(self)
        self.orb_tv.populate_tree(self.orb_data)
        self.orb_tv.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    def settings_screen(self) -> None:
        settings_window = tk.Toplevel(self)
        settings_window.title('Settings')
        settings_window.geometry('500x500')

        settings_frame = ttk.Frame(settings_window)
        settings_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Contour level
        contour_label = ttk.Label(settings_frame, text='Molecular Orbital Contour:')
        contour_label.grid(row=0, column=0, padx=5, pady=5)
        self.contour_entry = ttk.Entry(settings_frame)
        self.contour_entry.insert(0, str(self.contour))
        self.contour_entry.grid(row=1, column=0, padx=5, pady=5)

        # Opacity
        opacity_label = ttk.Label(settings_frame)
        opacity_label.grid(row=2, column=0, padx=5, pady=5)
        self.opacity_scale = ttk.Scale(
            settings_frame,
            length=150,
            command=lambda val: opacity_label.config(text=f'Molecular Orbital Opacity: {float(val):.2f}'),
        )
        self.opacity_scale.set(self.opacity)
        self.opacity_scale.grid(row=3, column=0, padx=5, pady=5)

        # Molecule Opacity
        molecule_opacity_label = ttk.Label(settings_frame)
        molecule_opacity_label.grid(row=4, column=0, padx=5, pady=5)
        self.molecule_opacity_scale = ttk.Scale(
            settings_frame,
            length=150,
            command=lambda val: molecule_opacity_label.config(text=f'Molecule Opacity: {float(val):.2f}'),
        )
        self.molecule_opacity_scale.set(self.molecule_opacity)
        self.molecule_opacity_scale.grid(row=5, column=0, padx=5, pady=5)

        # Toggle molecule visibility
        toggle_molecule_button = ttk.Button(
            settings_frame,
            text='Toggle Molecule',
            command=self.toggle_molecule,
        )
        toggle_molecule_button.grid(row=6, column=0, padx=5, pady=5)

        # Radius
        # FIX: move hover_text
        radius_label = self.hover_widget(
            ttk.Label,
            settings_frame,
            text='Maximum Radius:',
            hover_text='Maximum radius for molecular orbital. Default is 2 * max_radius of the molecule.',
        )
        radius_label.grid(row=0, column=1, padx=5, pady=5)
        self.radius_entry = ttk.Entry(settings_frame)
        self.radius_entry.insert(0, str(self.r_points[-1]))
        self.radius_entry.grid(row=1, column=1, padx=5, pady=5)

        # Radius points
        radius_points_label = ttk.Label(settings_frame, text='Number of Radius Points:')
        radius_points_label.grid(row=2, column=1, padx=5, pady=5)
        self.radius_points_entry = ttk.Entry(settings_frame)
        self.radius_points_entry.insert(0, str(len(self.r_points)))
        self.radius_points_entry.grid(row=3, column=1, padx=5, pady=5)

        # Theta points
        theta_points_label = ttk.Label(settings_frame, text='Number of Theta Points:')
        theta_points_label.grid(row=4, column=1, padx=5, pady=5)
        self.theta_points_entry = ttk.Entry(settings_frame)
        self.theta_points_entry.insert(0, str(len(self.theta_points)))
        self.theta_points_entry.grid(row=5, column=1, padx=5, pady=5)

        # Phi points
        phi_points_label = ttk.Label(settings_frame, text='Number of Phi Points:')
        phi_points_label.grid(row=6, column=1, padx=5, pady=5)
        self.phi_points_entry = ttk.Entry(settings_frame)
        self.phi_points_entry.insert(0, str(len(self.phi_points)))
        self.phi_points_entry.grid(row=7, column=1, padx=5, pady=5)

        # Reset button
        reset_button = ttk.Button(
            settings_frame,
            text='Reset',
            command=self.reset_settings,
        )
        reset_button.grid(row=8, column=0, padx=5, pady=5, columnspan=2)

        # Save settings button
        save_button = ttk.Button(
            settings_frame,
            text='Apply',
            command=self.apply_settings,
        )
        save_button.grid(row=9, column=0, padx=5, pady=5, columnspan=2)

    def reset_settings(self) -> None:
        """Reset settings to default values."""
        self.contour_entry.delete(0, tk.END)
        self.contour_entry.insert(0, str(self.CONTOUR))

        self.opacity_scale.set(self.OPACITY)

        self.molecule_opacity_scale.set(self.MOLECULE_OPACITY)

        self.radius_entry.delete(0, tk.END)
        self.radius_entry.insert(0, str(self.molecule.max_radius * 2))

        self.radius_points_entry.delete(0, tk.END)
        self.radius_points_entry.insert(0, str(self.NUM_RADIUS_POINTS))

        self.theta_points_entry.delete(0, tk.END)
        self.theta_points_entry.insert(0, str(self.NUM_THETA_POINTS))

        self.phi_points_entry.delete(0, tk.END)
        self.phi_points_entry.insert(0, str(self.NUM_PHI_POINTS))

    def apply_settings(self) -> None:
        self.molecule_opacity = round(self.molecule_opacity_scale.get(), 2)
        for actor in self.molecule_actors:
            actor.GetProperty().SetOpacity(self.molecule_opacity)

        radius = float(self.radius_entry.get())
        if radius <= 0:
            InvalidInputPopup('Invalid radius value: Must be greater than zero.')
            return

        num_r_points = int(self.radius_points_entry.get())
        num_theta_points = int(self.theta_points_entry.get())
        num_phi_points = int(self.phi_points_entry.get())

        if num_r_points <= 0 or num_theta_points <= 0 or num_phi_points <= 0:
            InvalidInputPopup('Invalid number of points: Must be greater than zero.')
            return

        if (
            radius != self.r_points[-1]
            or num_r_points != len(self.r_points)
            or num_theta_points != len(self.theta_points)
            or num_phi_points != len(self.phi_points)
        ):
            self.r_points = np.linspace(0, radius, num_r_points)
            self.theta_points = np.linspace(0, np.pi, num_theta_points)
            self.phi_points = np.linspace(0, 2 * np.pi, num_phi_points)

            # Update the mesh with new points
            self.update_points()

        self.contour = float(self.contour_entry.get().strip())
        self.plot_orbital(self.current_orb_ind)

        self.opacity = round(self.opacity_scale.get(), 2)
        if self.orb_actor:
            self.orb_actor.GetProperty().SetOpacity(self.opacity)

    def update_points(self) -> None:
        self.tab_obj.update_points(self.r_points, self.theta_points, self.phi_points)
        self.orb_mesh = self.create_mesh(
            self.r_points,
            self.theta_points,
            self.phi_points,
        )

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
        can_go_next = self.current_orb_ind < len(self.orb_data) - 1
        self.prev_button.config(state=tk.NORMAL if can_go_prev else tk.DISABLED)
        self.next_button.config(state=tk.NORMAL if can_go_next else tk.DISABLED)

    def create_mesh(
        self,
        r_points: np.ndarray,
        theta_points: np.ndarray,
        phi_points: np.ndarray,
    ) -> pv.StructuredGrid:
        """Create a mesh for the orbitals."""
        xyz_grid = self.tab_obj.spherical_to_cartesian(
            r_points,
            theta_points,
            phi_points,
        )
        dimensions = (len(phi_points), len(theta_points), len(r_points))

        mesh = pv.StructuredGrid()
        mesh.points = xyz_grid
        mesh.dimensions = dimensions

        return mesh

    def plot_orbital(self, orb_ind: int) -> None:
        if self.orb_actor:
            self.plotter.remove_actor(self.orb_actor)

        if orb_ind != -1:
            self.orb_mesh['orbital'] = self.tab_obj.tab_mo(orb_ind)

            contour_mesh = self.orb_mesh.contour(
                [-self.contour, self.contour],  # pyright: ignore[reportArgumentType]
            )

            self.orb_actor = self.plotter.add_mesh(
                contour_mesh,
                clim=[-self.contour, self.contour],
                opacity=self.opacity,
                show_scalar_bar=False,
                cmap='bwr',
            )


class _OrbitsTreeview(ttk.Treeview):
    def __init__(self, plot_window: OrbitalPlotWindow) -> None:
        columns = ['Symmetry', 'Energy [au]']
        widths = [100, 120]

        super().__init__(plot_window, columns=columns, show='headings', height=20)

        for col, w in zip(columns, widths):
            self.heading(col, text=col)
            self.column(col, width=w)

        self.plot_window = plot_window

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

    def populate_tree(self, orb_data: np.ndarray) -> None:
        self.erase()
        for ind, orbtital in enumerate(orb_data):
            self.insert('', 'end', iid=ind, values=tuple(orbtital))

    def on_select(self, _event: tk.Event) -> None:
        """Handle selection of an orbital."""
        selected_item = self.selection()
        self.selection_remove(selected_item)
        if selected_item:
            orb_ind = int(selected_item[0])
            self.highlight_orbital(orb_ind)
            self.plot_window.current_orb_ind = orb_ind
            self.plot_window.plot_orbital(orb_ind)
            self.plot_window.update_button_states()
