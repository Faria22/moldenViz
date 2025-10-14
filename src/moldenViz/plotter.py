"""Plotter module for creating plots of the molecule and it's orbitals."""

import logging
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib.colors as mcolors
import numpy as np
import pyvista as pv
from matplotlib.colors import LinearSegmentedColormap
from numpy.typing import NDArray
from pyvistaqt import BackgroundPlotter
from qtpy.QtWidgets import QAction, QMenu  # pyright: ignore[reportPrivateImportUsage]
from shiboken6 import isValid

from ._config_module import Config
from ._plotting_objects import Molecule
from .parser import _MolecularOrbital
from .tabulator import GridType, Tabulator, _cartesian_to_spherical, _spherical_to_cartesian

logger = logging.getLogger(__name__)

config = Config()


class Plotter:
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
        A `Tabulator` can be passed to tabulate the GTOs in a predetermined grid.

        Note: `Tabulator` grid must be spherical or cartesian. Custom grids are not allowed.
    tk_root : tk.Tk, optional
        If user is using the plotter inside a tk app, `tk_root` can be passed to not create a new tk instance.

    Attributes
    ----------
    on_screen : bool
        Indicates if the plotter window is currently open.
    tabulator : Tabulator
        The Tabulator object used for tabulating GTOs and MOs.
    molecule : Molecule
        The Molecule object representing the molecular structure.
    molecule_opacity : float
        The opacity of the molecule in the visualization.
    molecule_actors : list[pv.Actor]
        List of PyVista actors representing the molecule.
    atom_actors : list[pv.Actor]
        List of PyVista actors representing the atoms.
    bond_actors : list[pv.Actor]
        List of PyVista actors representing the bonds.
    tk_root : tk.Tk | None
        The Tkinter root window.
    pv_plotter : BackgroundPlotter
        The PyVista BackgroundPlotter for 3D rendering.
    molecule_actors : list[pv.Actor]
        List of PyVista actors representing the molecule.
    orb_mesh : pv.StructuredGrid
        The mesh used for visualizing molecular orbitals.
    orb_actor : pv.Actor | None
        The PyVista actor for the currently displayed molecular orbital.
    contour : float
        The contour level for molecular orbital visualization.
    opacity : float
        The opacity of the molecular orbital in the visualization.
    cmap : str | LinearSegmentedColormap
        The colormap used for molecular orbital visualization.

    Raises
    ------
    ValueError
        If the provided tabulator is invalid
        (e.g., missing grid or GTO data when `only_molecule` is `False`, or has an UNKNOWN grid type).
    """

    SPHERICAL_GRID_SETTINGS_WINDOW_SIZE = '400x350'
    CARTESIAN_GRID_SETTINGS_WINDOW_SIZE = '650x400'

    def __init__(
        self,
        source: str | list[str],
        only_molecule: bool = False,
        tabulator: Tabulator | None = None,
        tk_root: tk.Tk | None = None,
    ) -> None:
        self.on_screen = True
        self.only_molecule = only_molecule
        self.selection_screen: _OrbitalSelectionScreen | None = None

        self.tk_root = tk_root
        self._no_prev_tk_root = self.tk_root is None
        if self.tk_root is None:
            self.tk_root = tk.Tk()
            self.tk_root.withdraw()  # Hides window

        self.pv_plotter = BackgroundPlotter(editor=False)
        self.pv_plotter.set_background(config.background_color)
        self.pv_plotter.show_axes()

        self._add_orbital_menus_to_pv_plotter()
        self._connect_pv_plotter_close_signal()
        self._override_clear_all_button()

        if tabulator:
            if not hasattr(tabulator, 'grid'):
                raise ValueError('Tabulator does not have grid attribute.')

            if not hasattr(tabulator, 'gto_data') and not only_molecule:
                raise ValueError('Tabulator does not have tabulated GTOs.')

            if tabulator._grid_type == GridType.UNKNOWN:  # noqa: SLF001
                raise ValueError('The plotter only supports spherical and cartesian grids.')

            # Check if grid is uniform (PyVista requires uniform grids)
            if tabulator.original_axes is not None:
                Tabulator._axis_spacing(tabulator.original_axes[0], 'x')  # noqa: SLF001
                Tabulator._axis_spacing(tabulator.original_axes[1], 'y')  # noqa: SLF001
                Tabulator._axis_spacing(tabulator.original_axes[2], 'z')  # noqa: SLF001

            self.tabulator = tabulator
        else:
            self.tabulator = Tabulator(source, only_molecule=only_molecule)

        self.molecule: Molecule
        self.molecule_opacity = config.molecule.opacity
        self.load_molecule(config)

        # It no tabulator was passed, create default grid
        if not only_molecule and not tabulator:
            if config.grid.default_type == 'spherical':
                self.tabulator.spherical_grid(
                    np.linspace(
                        0,
                        max(config.grid.max_radius_multiplier * self.molecule.max_radius, config.grid.min_radius),
                        config.grid.spherical.num_r_points,
                    ),
                    np.linspace(0, np.pi, config.grid.spherical.num_theta_points),
                    np.linspace(0, 2 * np.pi, config.grid.spherical.num_phi_points),
                )
            else:  # cartesian
                r = max(config.grid.max_radius_multiplier * self.molecule.max_radius, config.grid.min_radius)
                self.tabulator.cartesian_grid(
                    np.linspace(-r, r, config.grid.cartesian.num_x_points),
                    np.linspace(-r, r, config.grid.cartesian.num_y_points),
                    np.linspace(-r, r, config.grid.cartesian.num_z_points),
                )

        # If we want to have the molecular orbitals, we need to initiate Tk before Qt
        # That is why we have this weird if statement separated this way
        if only_molecule:
            if self._no_prev_tk_root:
                self.tk_root.mainloop()
            return

        self.orb_mesh = self._create_mo_mesh()
        self.orb_actor: pv.Actor | None = None

        # Values for MO, not the molecule
        self.contour = config.mo.contour
        self.opacity = config.mo.opacity

        # Set colormap based on configuration
        if config.mo.custom_colors is not None:
            # Create custom colormap from two colors
            self.cmap = self.custom_cmap_from_colors(config.mo.custom_colors)
        else:
            self.cmap = config.mo.color_scheme

        if not self.only_molecule:
            self.selection_screen = _OrbitalSelectionScreen(self)

        if self._no_prev_tk_root:
            self.tk_root.mainloop()

    @staticmethod
    def custom_cmap_from_colors(colors: list[str]) -> LinearSegmentedColormap:
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

    def load_molecule(self, config: Config) -> None:
        """Reload the molecule from the parser data."""
        self.molecule = Molecule(self.tabulator._parser.atoms, config)  # noqa: SLF001

        for actor in self.molecule_actors if hasattr(self, 'molecule_actors') else []:
            self.pv_plotter.remove_actor(actor)

        self.molecule_actors, self.atom_actors, self.bond_actors = self.molecule.add_meshes(
            self.pv_plotter,
            self.molecule_opacity,
        )

    def _override_clear_all_button(self) -> None:
        """Override the default "Clear All" action in the PyVista plotter's View menu."""
        view_menu = None
        for action in self.pv_plotter.main_menu.actions():
            if action.text() == 'View':
                view_menu = action.menu()
                break

        if view_menu is None:
            raise RuntimeError('Could not find View menu in PyVista plotter.')

        for action in view_menu.actions():  # pyright: ignore[reportAttributeAccessIssue]
            if action is not None and isValid(action) and action.text().lower() == 'clear all':
                while action.triggered.disconnect():
                    pass
                action.triggered.connect(self._clear_all)
                break

    def _add_orbital_menus_to_pv_plotter(self) -> None:
        """Add Settings and Export menus to the PyVista plotter's main menu."""
        # Create Settings menu with dropdown
        settings_menu = QMenu('Settings', self.pv_plotter.app_window)

        # Add Settings submenu items
        if not self.only_molecule:
            grid_settings_action = QAction('Grid Settings', self.pv_plotter.app_window)
            grid_settings_action.triggered.connect(self.grid_settings_screen)
            settings_menu.addAction(grid_settings_action)

            mo_settings_action = QAction('MO Settings', self.pv_plotter.app_window)
            mo_settings_action.triggered.connect(self.mo_settings_screen)
            settings_menu.addAction(mo_settings_action)

        molecule_settings_action = QAction('Molecule Settings', self.pv_plotter.app_window)
        molecule_settings_action.triggered.connect(self.molecule_settings_screen)
        settings_menu.addAction(molecule_settings_action)

        color_settings_action = QAction('Color Settings', self.pv_plotter.app_window)
        color_settings_action.triggered.connect(self.color_settings_screen)
        settings_menu.addAction(color_settings_action)

        settings_menu.addSeparator()

        save_settings_action = QAction('Save Settings', self.pv_plotter.app_window)
        save_settings_action.triggered.connect(self.save_settings)
        settings_menu.addAction(save_settings_action)

        # Create Export menu with dropdown
        export_menu = QMenu('Export', self.pv_plotter.app_window)

        # Add Export submenu items
        if not self.only_molecule:
            export_data_action = QAction('Data', self.pv_plotter.app_window)
            export_data_action.triggered.connect(self.export_orbitals_dialog)
            export_menu.addAction(export_data_action)

        export_image_action = QAction('Image', self.pv_plotter.app_window)
        export_image_action.triggered.connect(self.export_image_dialog)
        export_menu.addAction(export_image_action)

        # Add menus to main menu bar
        self.pv_plotter.main_menu.addMenu(settings_menu)
        self.pv_plotter.main_menu.addMenu(export_menu)

    def _do_export(self, export_window: tk.Toplevel, format_var: tk.StringVar, scope_var: tk.StringVar) -> None:
        """Execute the export operation.

        Parameters
        ----------
        export_window : tk.Toplevel
            The export dialog window to close on success.
        format_var : tk.StringVar
            Variable holding the selected export format ('vtk' or 'cube').
        scope_var : tk.StringVar
            Variable holding the selected scope ('current' or 'all').
        """
        assert self.selection_screen is not None

        file_format = format_var.get()
        scope = scope_var.get()

        # Validate selection
        if scope == 'current' and self.selection_screen.current_mo_ind < 0:
            messagebox.showerror('Export Error', 'No orbital is currently selected.')
            return

        if file_format == 'cube' and scope == 'all':
            messagebox.showerror(
                'Export Error',
                'Cube format only supports exporting a single orbital.\n\n'
                'Please select "Current orbital" or choose VTK format.',
            )
            return

        # Determine file extension and default name
        ext = '.vtk' if file_format == 'vtk' else '.cube'
        default_name = (
            f'orbitals_all{ext}' if scope == 'all' else f'orbital_{self.selection_screen.current_mo_ind}{ext}'
        )

        # Show file save dialog
        file_path = filedialog.asksaveasfilename(
            parent=export_window,
            title='Save Orbital Export',
            defaultextension=ext,
            initialfile=default_name,
            filetypes=[('VTK Files', '*.vtk'), ('Gaussian Cube Files', '*.cube'), ('All Files', '*.*')],
        )

        if not file_path:
            return  # User cancelled

        # Perform the export
        try:
            mo_index = self.selection_screen.current_mo_ind if scope == 'current' else None
            self.tabulator.export(file_path, mo_index=mo_index)
            messagebox.showinfo('Export Successful', f'Orbital(s) exported successfully to:\n{file_path}')
            export_window.destroy()
        except (RuntimeError, ValueError) as e:
            messagebox.showerror('Export Failed', f'Failed to export orbital(s):\n\n{e!s}')

    def export_orbitals_dialog(self) -> None:
        """Open a dialog to configure and export molecular orbitals."""
        assert self.selection_screen is not None

        export_window = tk.Toplevel(self.tk_root)
        export_window.title('Export Orbitals')
        export_window.geometry('400x300')

        main_frame = ttk.Frame(export_window, padding=20)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # File format selection
        ttk.Label(main_frame, text='Export Format:', font=('TkDefaultFont', 10, 'bold')).grid(
            row=0,
            column=0,
            columnspan=2,
            sticky=tk.W,
            pady=(0, 10),
        )

        format_var = tk.StringVar(value='vtk')
        ttk.Radiobutton(main_frame, text='VTK (.vtk) - All orbitals or single', variable=format_var, value='vtk').grid(
            row=1,
            column=0,
            columnspan=2,
            sticky=tk.W,
            padx=20,
        )
        ttk.Radiobutton(
            main_frame,
            text='Gaussian Cube (.cube) - Single orbital only',
            variable=format_var,
            value='cube',
        ).grid(row=2, column=0, columnspan=2, sticky=tk.W, padx=20, pady=(5, 15))

        # Orbital selection
        ttk.Label(main_frame, text='Orbital Selection:', font=('TkDefaultFont', 10, 'bold')).grid(
            row=3,
            column=0,
            columnspan=2,
            sticky=tk.W,
            pady=(0, 10),
        )

        scope_var = tk.StringVar(value='current')
        # Use 1-based indexing for display (add 1 to current_mo_ind)
        orbital_display = (
            self.selection_screen.current_mo_ind + 1 if self.selection_screen.current_mo_ind >= 0 else 'None'
        )
        current_orb_radio = ttk.Radiobutton(
            main_frame,
            text=f'Current orbital (#{orbital_display})',
            variable=scope_var,
            value='current',
        )
        current_orb_radio.grid(row=4, column=0, columnspan=2, sticky=tk.W, padx=20)
        if self.selection_screen.current_mo_ind < 0:
            current_orb_radio.config(state=tk.DISABLED)

        all_orb_radio = ttk.Radiobutton(main_frame, text='All orbitals', variable=scope_var, value='all')
        all_orb_radio.grid(row=5, column=0, columnspan=2, sticky=tk.W, padx=20, pady=(5, 0))

        # Store references for updating the label dynamically
        self._export_window = export_window
        self._export_current_orb_radio = current_orb_radio
        self._export_all_orb_radio = all_orb_radio

        def update_scope_options(*_args: object) -> None:
            """Adjust which export scopes are available based on the format."""
            if self._export_all_orb_radio is None:
                return

            if format_var.get() == 'cube':
                self._export_all_orb_radio.config(state=tk.DISABLED)
                if scope_var.get() == 'all':
                    scope_var.set('current')
            else:
                self._export_all_orb_radio.config(state=tk.NORMAL)

        format_var.trace_add('write', update_scope_options)
        update_scope_options()

        # Clean up references when window is closed
        def on_close() -> None:
            self._export_window = None
            self._export_current_orb_radio = None
            self._export_all_orb_radio = None
            export_window.destroy()

        # Buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=6, column=0, columnspan=2, pady=(20, 0))

        ttk.Button(
            button_frame,
            text='Export',
            command=lambda: self._do_export(export_window, format_var, scope_var),
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text='Cancel', command=on_close).pack(side=tk.LEFT, padx=5)
        export_window.protocol('WM_DELETE_WINDOW', on_close)

    def export_image_dialog(self) -> None:
        """Open a dialog to export the current visualization as an image."""
        export_window = tk.Toplevel(self.tk_root)
        export_window.title('Export Image')
        export_window.geometry('400x250')

        main_frame = ttk.Frame(export_window, padding=20)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # File format selection
        ttk.Label(main_frame, text='Image Format:', font=('TkDefaultFont', 10, 'bold')).grid(
            row=0,
            column=0,
            columnspan=2,
            sticky=tk.W,
            pady=(0, 10),
        )

        format_var = tk.StringVar(value='png')
        ttk.Radiobutton(main_frame, text='PNG (.png) - Raster format', variable=format_var, value='png').grid(
            row=1,
            column=0,
            columnspan=2,
            sticky=tk.W,
            padx=20,
        )
        ttk.Radiobutton(main_frame, text='JPEG (.jpg) - Raster format', variable=format_var, value='jpeg').grid(
            row=2,
            column=0,
            columnspan=2,
            sticky=tk.W,
            padx=20,
            pady=(5, 0),
        )
        ttk.Radiobutton(main_frame, text='SVG (.svg) - Vector format', variable=format_var, value='svg').grid(
            row=3,
            column=0,
            columnspan=2,
            sticky=tk.W,
            padx=20,
            pady=(5, 0),
        )
        ttk.Radiobutton(main_frame, text='PDF (.pdf) - Vector format', variable=format_var, value='pdf').grid(
            row=4,
            column=0,
            columnspan=2,
            sticky=tk.W,
            padx=20,
            pady=(5, 15),
        )

        # Transparent background option (only for PNG)
        transparent_var = tk.BooleanVar(value=False)
        transparent_check = ttk.Checkbutton(
            main_frame,
            text='Transparent background (PNG only)',
            variable=transparent_var,
        )
        transparent_check.grid(row=5, column=0, columnspan=2, sticky=tk.W, padx=20, pady=(0, 15))

        def update_transparent_option(*_args: object) -> None:
            """Enable/disable transparent option based on format."""
            if format_var.get() == 'png':
                transparent_check.config(state=tk.NORMAL)
            else:
                transparent_check.config(state=tk.DISABLED)

        format_var.trace_add('write', update_transparent_option)
        update_transparent_option()

        # Buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.grid(row=6, column=0, columnspan=2, pady=(20, 0))

        ttk.Button(
            button_frame,
            text='Export',
            command=lambda: self._do_image_export(export_window, format_var, transparent_var),
        ).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text='Cancel', command=export_window.destroy).pack(side=tk.LEFT, padx=5)
        export_window.protocol('WM_DELETE_WINDOW', export_window.destroy)

    def _do_image_export(
        self,
        export_window: tk.Toplevel,
        format_var: tk.StringVar,
        transparent_var: tk.BooleanVar,
    ) -> None:
        """Execute the image export operation.

        Parameters
        ----------
        export_window : tk.Toplevel
            The export dialog window to close on success.
        format_var : tk.StringVar
            Variable holding the selected export format ('png', 'jpeg', 'svg', or 'pdf').
            Note: JPEG files are saved with .jpg extension (the standard).
        transparent_var : tk.BooleanVar
            Variable indicating whether to use a transparent background (PNG only).
        """
        file_format = format_var.get()
        transparent = transparent_var.get()

        # Determine file extension and default name
        # Note: JPEG format uses .jpg as the standard extension
        ext_map = {'png': '.png', 'jpeg': '.jpg', 'svg': '.svg', 'pdf': '.pdf'}
        ext = ext_map[file_format]
        default_name = f'moldenviz_export{ext}'

        # Define file types for dialog
        file_types = {
            'png': ('PNG Files', '*.png'),
            'jpeg': ('JPEG Files', '*.jpg *.jpeg'),
            'svg': ('SVG Files', '*.svg'),
            'pdf': ('PDF Files', '*.pdf'),
        }

        # Show file save dialog
        file_path = filedialog.asksaveasfilename(
            parent=export_window,
            title='Save Image Export',
            defaultextension=ext,
            initialfile=default_name,
            filetypes=[file_types[file_format], ('All Files', '*.*')],
        )

        if not file_path:
            return  # User cancelled

        # Perform the export
        try:
            if file_format in {'svg', 'pdf'}:
                # Use save_graphic for vector formats
                self.pv_plotter.save_graphic(file_path)
            else:
                # Use screenshot for raster formats
                self.pv_plotter.screenshot(
                    file_path,
                    transparent_background=transparent if file_format == 'png' else False,
                )

            messagebox.showinfo('Export Successful', f'Image exported successfully to:\n{file_path}')
            export_window.destroy()
        except (RuntimeError, OSError, ValueError) as e:
            messagebox.showerror('Export Failed', f'Failed to export image:\n\n{e!s}')

    def _settings_parent(self) -> tk.Misc:
        """Return the appropriate parent widget for settings dialogs.

        Returns
        -------
        tk.Misc
            The parent widget for settings dialogs.
        """
        parent = self.selection_screen if self.selection_screen is not None else self.tk_root
        if parent is None:
            raise RuntimeError('No Tk root available to host settings dialogs.')
        return parent

    def _get_current_mo_index(self) -> int:
        """Return the currently selected molecular orbital index.

        Returns
        -------
        int
            The index of the currently selected molecular orbital, or -1 if none is selected.

        """
        if self.selection_screen:
            return self.selection_screen.current_mo_ind
        return -1

    def grid_settings_screen(self) -> None:
        """Open the grid settings window."""
        parent = self._settings_parent()
        self.grid_settings_window = tk.Toplevel(parent)
        self.grid_settings_window.title('Grid Settings')

        settings_frame = ttk.Frame(self.grid_settings_window)
        settings_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Grid parameters
        ttk.Label(settings_frame, text='MO Grid parameters').grid(row=0, column=0, padx=5, pady=5, columnspan=5)

        self.grid_type_radio_var = tk.StringVar()
        self.grid_type_radio_var.set(self.tabulator._grid_type.value)  # noqa: SLF001

        ttk.Label(settings_frame, text='Spherical grid:').grid(row=1, column=0, padx=5, pady=5)
        sph_grid_type_button = ttk.Radiobutton(
            settings_frame,
            variable=self.grid_type_radio_var,
            value=GridType.SPHERICAL.value,
            command=self.place_grid_params_frame,
        )

        ttk.Label(settings_frame, text='Cartesian grid:').grid(row=1, column=2, padx=5, pady=5)
        cart_grid_type_button = ttk.Radiobutton(
            settings_frame,
            variable=self.grid_type_radio_var,
            value=GridType.CARTESIAN.value,
            command=self.place_grid_params_frame,
        )

        sph_grid_type_button.grid(row=1, column=1, padx=5, pady=5)
        cart_grid_type_button.grid(row=1, column=3, padx=5, pady=5)

        self.sph_grid_params_frame = self.sph_grid_params_frame_widgets(settings_frame)
        self.cart_grid_params_frame = self.cart_grid_params_frame_widgets(settings_frame)

        self.place_grid_params_frame()

        # Reset button
        reset_button = ttk.Button(settings_frame, text='Reset', command=self.reset_grid_settings)
        reset_button.grid(row=8, column=0, padx=5, pady=5, columnspan=5)

        # Apply settings button
        apply_button = ttk.Button(settings_frame, text='Apply', command=self.apply_grid_settings)
        apply_button.grid(row=9, column=0, padx=5, pady=5, columnspan=5)

    def mo_settings_screen(self) -> None:
        """Open the molecular orbital settings window."""
        parent = self._settings_parent()
        self.mo_settings_window = tk.Toplevel(parent)
        self.mo_settings_window.title('MO Settings')

        settings_frame = ttk.Frame(self.mo_settings_window)
        settings_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Contour level
        contour_label = ttk.Label(settings_frame, text='Molecular Orbital Contour:')
        contour_label.grid(row=0, column=0, padx=5, pady=5, sticky='w')
        self.contour_entry = ttk.Entry(settings_frame)
        self.contour_entry.insert(0, str(self.contour))
        self.contour_entry.grid(row=1, column=0, padx=5, pady=5, sticky='ew')

        # Bind to apply changes on Enter key or focus out
        self.contour_entry.bind('<Return>', lambda _e: self.apply_mo_contour())
        self.contour_entry.bind('<FocusOut>', lambda _e: self.apply_mo_contour())

        # Opacity
        opacity_label = ttk.Label(settings_frame)
        opacity_label.grid(row=2, column=0, padx=5, pady=5, sticky='w')
        self.opacity_scale = ttk.Scale(
            settings_frame,
            length=200,
            command=self.on_opacity_change,
        )
        self.opacity_scale.set(self.opacity)
        self.opacity_scale.grid(row=3, column=0, padx=5, pady=5, sticky='ew')
        # Initialize label
        opacity_label.config(text=f'Molecular Orbital Opacity: {self.opacity:.2f}')

        # Configure grid column weight for proper resizing
        settings_frame.columnconfigure(0, weight=1)

        # Reset button
        reset_button = ttk.Button(settings_frame, text='Reset', command=self.reset_mo_settings)
        reset_button.grid(row=4, column=0, padx=5, pady=5, sticky='ew')

    def on_opacity_change(self, val: str) -> None:
        """Handle opacity slider changes and apply immediately."""
        opacity = round(float(val), 2)
        self.opacity = opacity
        if self.orb_actor:
            self.orb_actor.GetProperty().SetOpacity(opacity)

        # Update label
        for widget in self.mo_settings_window.winfo_children():
            for child in widget.winfo_children():
                if isinstance(child, ttk.Label) and 'Molecular Orbital Opacity:' in child.cget('text'):
                    child.config(text=f'Molecular Orbital Opacity: {opacity:.2f}')

    def apply_mo_contour(self) -> None:
        """Apply contour changes immediately."""
        try:
            new_contour = float(self.contour_entry.get().strip())
            self.contour = new_contour
            # Replot the current orbital with the new contour
            idx = self._get_current_mo_index()
            if idx >= 0:
                self.plot_orbital(idx)
        except ValueError:
            pass  # Ignore invalid input

    def update_settings_button_states(self) -> None:
        """Update the state of the settings buttons based on current plotter state."""
        if hasattr(self, 'show_atoms_var'):
            self.show_atoms_var.set(self.are_atoms_visible())
        if hasattr(self, 'show_bonds_var'):
            self.show_bonds_var.set(self.are_bonds_visible())

        config.molecule.atom.show = self.are_atoms_visible()
        config.molecule.bond.show = self.are_bonds_visible()

    def molecule_settings_screen(self) -> None:
        """Open the molecule settings window."""
        parent = self._settings_parent()
        self.molecule_settings_window = tk.Toplevel(parent)
        self.molecule_settings_window.title('Molecule Settings')

        settings_frame = ttk.Frame(self.molecule_settings_window)
        settings_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Molecule Opacity
        molecule_opacity_label = ttk.Label(settings_frame)
        molecule_opacity_label.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky='w')
        self.molecule_opacity_scale = ttk.Scale(
            settings_frame,
            length=100,
            command=self.on_molecule_opacity_change,
        )
        self.molecule_opacity_scale.set(self.molecule_opacity)
        self.molecule_opacity_scale.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky='ew')
        # Initialize label
        molecule_opacity_label.config(text=f'Molecule Opacity: {self.molecule_opacity:.2f}')

        # Toggle molecule visibility
        toggle_mol_button = ttk.Button(
            settings_frame,
            text='Toggle Molecule',
            command=self.toggle_molecule,
            width=20,
        )
        toggle_mol_button.grid(row=2, column=0, columnspan=2, padx=5, pady=5, sticky='ew')

        # Separator
        ttk.Separator(settings_frame, orient='horizontal').grid(row=3, column=0, columnspan=2, sticky='ew', pady=10)

        # Show atoms checkbox
        self.show_atoms_var = tk.BooleanVar(value=config.molecule.atom.show)
        show_atoms_check = ttk.Checkbutton(
            settings_frame,
            text='Show Atoms',
            variable=self.show_atoms_var,
            command=self.toggle_atoms,
        )
        show_atoms_check.grid(row=5, column=0, columnspan=2, padx=5, pady=5, sticky='w')

        # Show bonds checkbox
        self.show_bonds_var = tk.BooleanVar(value=config.molecule.bond.show)
        show_bonds_check = ttk.Checkbutton(
            settings_frame,
            text='Show Bonds',
            variable=self.show_bonds_var,
            command=self.toggle_bonds,
        )
        show_bonds_check.grid(row=8, column=0, columnspan=2, padx=5, pady=5, sticky='w')

        # Separator
        ttk.Separator(settings_frame, orient='horizontal').grid(row=9, column=0, columnspan=2, sticky='ew', pady=10)

        # Bond max length
        ttk.Label(settings_frame, text='Max Bond Length:').grid(row=10, column=0, padx=5, pady=5, sticky='w')
        self.bond_max_length_entry = ttk.Entry(settings_frame, width=15)
        self.bond_max_length_entry.insert(0, str(config.molecule.bond.max_length))
        self.bond_max_length_entry.grid(row=10, column=1, padx=5, pady=5, sticky='w')

        # Bond radius
        ttk.Label(settings_frame, text='Bond Radius:').grid(row=12, column=0, padx=5, pady=5, sticky='w')
        self.bond_radius_entry = ttk.Entry(settings_frame, width=15)
        self.bond_radius_entry.insert(0, str(config.molecule.bond.radius))
        self.bond_radius_entry.grid(row=12, column=1, padx=5, pady=5, sticky='w')

        # Configure grid column weights for proper resizing
        settings_frame.columnconfigure(0, weight=1)
        settings_frame.columnconfigure(1, weight=1)

        # Reset button
        reset_button = ttk.Button(settings_frame, text='Reset', command=self.reset_molecule_settings)
        reset_button.grid(row=13, column=0, columnspan=2, padx=5, pady=5, sticky='ew')

        # Apply settings button
        apply_button = ttk.Button(settings_frame, text='Apply', command=self.apply_molecule_settings)
        apply_button.grid(row=14, column=0, columnspan=2, padx=5, pady=5, sticky='ew')

    def on_molecule_opacity_change(self, val: str) -> None:
        """Handle molecule opacity slider changes and apply immediately."""
        opacity = round(float(val), 2)
        self.molecule_opacity = opacity
        for actor in self.molecule_actors:
            actor.GetProperty().SetOpacity(opacity)
        # Update label
        for widget in self.molecule_settings_window.winfo_children():
            for child in widget.winfo_children():
                if isinstance(child, ttk.Label) and 'Molecule Opacity:' in child.cget('text'):
                    child.config(text=f'Molecule Opacity: {opacity:.2f}')

    def apply_background_color(self) -> None:
        """Apply background color changes immediately."""
        try:
            color = self.background_color_entry.get().strip()
            if mcolors.is_color_like(color):
                self.pv_plotter.set_background(color)
            else:
                messagebox.showerror('Invalid Input', f'"{color}" is not a valid color.')
        except (ValueError, RuntimeError) as e:
            messagebox.showerror('Error', f'Failed to set background color: {e!s}')

    def on_mo_color_scheme_change(self, _event: tk.Event) -> None:
        """Handle MO color scheme dropdown change to show/hide custom color entries."""
        if self.mo_color_scheme_var.get() == 'custom':
            for widget in self.mo_custom_color_widgets:
                widget.grid()
        else:
            for widget in self.mo_custom_color_widgets:
                widget.grid_remove()

    def on_bond_color_type_change(self) -> None:
        """Handle bond color type change to show/hide bond color entry."""
        if self.bond_color_type_var.get() == 'uniform':
            self.bond_color_label.grid()
            self.bond_color_entry.grid()
        else:
            self.bond_color_label.grid_remove()
            self.bond_color_entry.grid_remove()

        self.apply_bond_color_settings()

    def color_settings_screen(self) -> None:
        """Open the color settings window."""
        parent = self._settings_parent()
        self.color_settings_window = tk.Toplevel(parent)
        self.color_settings_window.title('Color Settings')

        settings_frame = ttk.Frame(self.color_settings_window)
        settings_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Background Color section
        ttk.Label(settings_frame, text='Background Color', font=('TkDefaultFont', 10, 'bold')).grid(
            row=0,
            column=0,
            columnspan=2,
            padx=5,
            pady=5,
            sticky='w',
        )

        ttk.Label(settings_frame, text='Background Color:').grid(row=1, column=0, padx=5, pady=5, sticky='w')
        self.background_color_entry = ttk.Entry(settings_frame, width=15)
        self.background_color_entry.insert(0, str(config.background_color))
        self.background_color_entry.grid(row=1, column=1, padx=5, pady=5, sticky='w')
        # Bind to apply changes on Enter key or focus out
        self.background_color_entry.bind('<Return>', lambda _e: self.apply_background_color())
        self.background_color_entry.bind('<FocusOut>', lambda _e: self.apply_background_color())

        # Separator
        ttk.Separator(settings_frame, orient='horizontal').grid(row=2, column=0, columnspan=2, sticky='ew', pady=10)

        if not self.only_molecule:
            # MO Color section
            ttk.Label(settings_frame, text='Molecular Orbital Colors', font=('TkDefaultFont', 10, 'bold')).grid(
                row=3,
                column=0,
                columnspan=2,
                padx=5,
                pady=5,
                sticky='w',
            )

            ttk.Label(settings_frame, text='Color Scheme:').grid(row=4, column=0, padx=5, pady=5, sticky='w')
            # Create dropdown with predefined color schemes
            predefined_schemes = ['bwr', 'RdBu', 'seismic', 'coolwarm', 'PiYG']

            # Check if user has a custom scheme in config that's not in predefined list
            if config.mo.color_scheme not in predefined_schemes and config.mo.color_scheme != 'custom':
                # Add the user's scheme as the first item
                color_schemes = [config.mo.color_scheme, *predefined_schemes, 'custom']
                default_scheme = config.mo.color_scheme
            else:
                color_schemes = [*predefined_schemes, 'custom']
                default_scheme = config.mo.color_scheme if config.mo.color_scheme in predefined_schemes else 'custom'

            self.mo_color_scheme_var = tk.StringVar(value=default_scheme)
            self.mo_color_scheme_dropdown = ttk.Combobox(
                settings_frame,
                state='readonly',
                textvariable=self.mo_color_scheme_var,
                values=color_schemes,
            )
            self.mo_color_scheme_dropdown.grid(row=4, column=1, padx=5, pady=5, sticky='w')
            self.mo_color_scheme_dropdown.bind('<<ComboboxSelected>>', self.on_mo_color_scheme_change)
            self.mo_color_scheme_dropdown.bind('<<ComboboxSelected>>', lambda _e: self.apply_mo_color_settings())

            # Custom color entries (initially hidden unless 'custom' is selected)
            negative_color_label = ttk.Label(settings_frame, text='Negative Color:')
            negative_color_label.grid(row=5, column=0, padx=5, pady=5, sticky='w')
            self.mo_negative_color_entry = ttk.Entry(settings_frame, width=15)
            if config.mo.custom_colors and len(config.mo.custom_colors) > 0:
                self.mo_negative_color_entry.insert(0, config.mo.custom_colors[0])
            self.mo_negative_color_entry.grid(row=5, column=1, padx=5, pady=5, sticky='w')

            positive_color_label = ttk.Label(settings_frame, text='Positive Color:')
            positive_color_label.grid(row=6, column=0, padx=5, pady=5, sticky='w')
            self.mo_positive_color_entry = ttk.Entry(settings_frame, width=15)
            if config.mo.custom_colors and len(config.mo.custom_colors) > 1:
                self.mo_positive_color_entry.insert(0, config.mo.custom_colors[1])
            self.mo_positive_color_entry.grid(row=6, column=1, padx=5, pady=5, sticky='w')

            # Store references to custom color widgets for show/hide
            self.mo_custom_color_widgets = [
                self.mo_negative_color_entry,
                self.mo_positive_color_entry,
                negative_color_label,
                positive_color_label,
            ]

            # Hide custom color entries if predefined scheme is selected
            if self.mo_color_scheme_var.get() != 'custom':
                for widget in self.mo_custom_color_widgets:
                    widget.grid_remove()

            # Separator
            ttk.Separator(settings_frame, orient='horizontal').grid(row=7, column=0, columnspan=2, sticky='ew', pady=10)

        # Bond Color section
        ttk.Label(settings_frame, text='Bond Colors', font=('TkDefaultFont', 10, 'bold')).grid(
            row=8,
            column=0,
            columnspan=2,
            padx=5,
            pady=5,
            sticky='w',
        )

        # Bond color type
        ttk.Label(settings_frame, text='Bond Color Type:').grid(row=9, column=0, padx=5, pady=5, sticky='w')
        self.bond_color_type_var = tk.StringVar(value=config.molecule.bond.color_type)
        bond_color_frame = ttk.Frame(settings_frame)
        bond_color_frame.grid(row=9, column=1, padx=5, pady=5, sticky='w')
        ttk.Radiobutton(
            bond_color_frame,
            text='Uniform',
            variable=self.bond_color_type_var,
            value='uniform',
            command=self.on_bond_color_type_change,
        ).pack(side=tk.LEFT)
        ttk.Radiobutton(
            bond_color_frame,
            text='Split',
            variable=self.bond_color_type_var,
            value='split',
            command=self.on_bond_color_type_change,
        ).pack(side=tk.LEFT)

        # Bond color (for uniform type only)
        self.bond_color_label = ttk.Label(settings_frame, text='Bond Color:')
        self.bond_color_label.grid(row=10, column=0, padx=5, pady=5, sticky='w')
        self.bond_color_entry = ttk.Entry(settings_frame, width=15)
        self.bond_color_entry.insert(0, str(config.molecule.bond.color))
        self.bond_color_entry.grid(row=10, column=1, padx=5, pady=5, sticky='w')
        self.bond_color_entry.bind('<Return>', lambda _e: self.apply_bond_color_settings())

        # Hide bond color entry if split is selected
        if self.bond_color_type_var.get() == 'split':
            self.bond_color_label.grid_remove()
            self.bond_color_entry.grid_remove()

        # Configure grid column weights for proper resizing
        settings_frame.columnconfigure(0, weight=1)
        settings_frame.columnconfigure(1, weight=1)

        # Reset button
        reset_button = ttk.Button(settings_frame, text='Reset', command=self.reset_color_settings)
        reset_button.grid(row=11, column=0, columnspan=2, padx=5, pady=5, sticky='ew')

        # Apply settings button
        apply_button = ttk.Button(settings_frame, text='Apply', command=self.apply_color_settings)
        apply_button.grid(row=13, column=0, columnspan=2, padx=5, pady=5, sticky='ew')

    def place_grid_params_frame(self) -> None:
        """Render the parameter frame that matches the selected grid type."""
        if self.grid_type_radio_var.get() == GridType.SPHERICAL.value:
            self.grid_settings_window.geometry(self.SPHERICAL_GRID_SETTINGS_WINDOW_SIZE)
            self.cart_grid_params_frame.grid_forget()
            self.sph_grid_params_frame.grid(row=2, column=0, padx=5, pady=5, rowspan=6, columnspan=4)
            self.sph_grid_params_frame_setup()
        else:
            self.grid_settings_window.geometry(self.CARTESIAN_GRID_SETTINGS_WINDOW_SIZE)
            self.sph_grid_params_frame.grid_forget()
            self.cart_grid_params_frame.grid(row=2, column=0, padx=5, pady=5, rowspan=6, columnspan=4)
            self.cart_grid_params_frame_setup()

    def sph_grid_params_frame_widgets(self, master: ttk.Frame) -> ttk.Frame:
        """Build widgets that capture spherical grid parameters.

        Parameters
        ----------
        master : ttk.Frame
            The parent frame to contain the spherical grid parameter widgets.

        Returns
        -------
        ttk.Frame
            The frame containing the spherical grid parameter widgets.

        """
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
        ttk.Label(grid_params_frame, text='Number of Theta Points:').grid(row=2, column=0, padx=5, pady=5)
        self.theta_points_entry = ttk.Entry(grid_params_frame)
        self.theta_points_entry.grid(row=2, column=1, padx=5, pady=5)

        # Phi points
        ttk.Label(grid_params_frame, text='Number of Phi Points:').grid(row=3, column=0, padx=5, pady=5)
        self.phi_points_entry = ttk.Entry(grid_params_frame)
        self.phi_points_entry.grid(row=3, column=1, padx=5, pady=5)

        return grid_params_frame

    def cart_grid_params_frame_widgets(self, master: ttk.Frame) -> ttk.Frame:
        """Build widgets that capture cartesian grid parameters.

        Parameters
        ----------
        master : ttk.Frame
            The parent frame to contain the cartesian grid parameter widgets.

        Returns
        -------
        ttk.Frame
            The frame containing the cartesian grid parameter widgets.
        """
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
        """Populate the spherical grid widgets with defaults or existing values."""
        self.radius_entry.delete(0, tk.END)
        self.radius_points_entry.delete(0, tk.END)
        self.theta_points_entry.delete(0, tk.END)
        self.phi_points_entry.delete(0, tk.END)

        # Previous grid was cartesian, so use default values
        if self.tabulator._grid_type == GridType.CARTESIAN:  # noqa: SLF001
            self.radius_entry.insert(
                0,
                str(max(config.grid.max_radius_multiplier * self.molecule.max_radius, config.grid.min_radius)),
            )
            self.radius_points_entry.insert(0, str(config.grid.spherical.num_r_points))
            self.theta_points_entry.insert(0, str(config.grid.spherical.num_theta_points))
            self.phi_points_entry.insert(0, str(config.grid.spherical.num_phi_points))
            return

        num_r, num_theta, num_phi = self.tabulator._grid_dimensions  # noqa: SLF001

        # The last point of the grid for sure has the largest r
        r, _, _ = _cartesian_to_spherical(*self.tabulator.grid[-1, :])  # pyright: ignore[reportArgumentType]

        self.radius_entry.insert(0, str(r))
        self.radius_points_entry.insert(0, str(num_r))
        self.theta_points_entry.insert(0, str(num_theta))
        self.phi_points_entry.insert(0, str(num_phi))

    def cart_grid_params_frame_setup(self) -> None:
        """Populate the Cartesian grid widgets with defaults or existing values."""
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
        if self.tabulator._grid_type == GridType.SPHERICAL:  # noqa: SLF001
            r = max(config.grid.max_radius_multiplier * self.molecule.max_radius, config.grid.min_radius)

            self.x_min_entry.insert(0, str(-r))
            self.y_min_entry.insert(0, str(-r))
            self.z_min_entry.insert(0, str(-r))

            self.x_max_entry.insert(0, str(r))
            self.y_max_entry.insert(0, str(r))
            self.z_max_entry.insert(0, str(r))

            self.x_num_points_entry.insert(0, str(config.grid.cartesian.num_x_points))
            self.y_num_points_entry.insert(0, str(config.grid.cartesian.num_y_points))
            self.z_num_points_entry.insert(0, str(config.grid.cartesian.num_z_points))
            return

        x_num, y_num, z_num = self.tabulator._grid_dimensions  # noqa: SLF001
        x_min, y_min, z_min = self.tabulator.grid[0, :]
        x_max, y_max, z_max = self.tabulator.grid[-1, :]

        self.x_min_entry.insert(0, str(x_min))
        self.x_max_entry.insert(0, str(x_max))
        self.x_num_points_entry.insert(0, str(x_num))

        self.y_min_entry.insert(0, str(y_min))
        self.y_max_entry.insert(0, str(y_max))
        self.y_num_points_entry.insert(0, str(y_num))

        self.z_min_entry.insert(0, str(z_min))
        self.z_max_entry.insert(0, str(z_max))
        self.z_num_points_entry.insert(0, str(z_num))

    def reset_grid_settings(self) -> None:
        """Restore grid settings widgets back to configuration defaults."""
        self.grid_type_radio_var.set(config.grid.default_type)

        self.radius_entry.delete(0, tk.END)
        self.radius_entry.insert(
            0,
            str(max(config.grid.max_radius_multiplier * self.molecule.max_radius, config.grid.min_radius)),
        )

        self.radius_points_entry.delete(0, tk.END)
        self.radius_points_entry.insert(0, str(config.grid.spherical.num_r_points))

        self.theta_points_entry.delete(0, tk.END)
        self.theta_points_entry.insert(0, str(config.grid.spherical.num_theta_points))

        self.phi_points_entry.delete(0, tk.END)
        self.phi_points_entry.insert(0, str(config.grid.spherical.num_phi_points))

        self.place_grid_params_frame()

    def reset_mo_settings(self) -> None:
        """Restore MO settings widgets back to configuration defaults."""
        self.contour_entry.delete(0, tk.END)
        self.contour_entry.insert(0, str(config.mo.contour))

        self.opacity_scale.set(config.mo.opacity)

        self.apply_mo_contour()  # Reapply contour with new value

    def reset_molecule_settings(self) -> None:
        """Restore molecule settings widgets back to configuration defaults."""
        config = Config()  # Reload config to discard unsaved changes

        self.molecule_opacity_scale.set(config.molecule.opacity)

        self.show_atoms_var.set(config.molecule.atom.show)
        self.show_bonds_var.set(config.molecule.bond.show)

        if not self.are_atoms_visible():
            self.toggle_atoms()
        if not self.are_bonds_visible():
            self.toggle_bonds()

        self.bond_max_length_entry.delete(0, tk.END)
        self.bond_max_length_entry.insert(0, str(config.molecule.bond.max_length))

        self.bond_radius_entry.delete(0, tk.END)
        self.bond_radius_entry.insert(0, str(config.molecule.bond.radius))

        self.apply_molecule_settings()  # Reapply molecule settings with new values

    def reset_color_settings(self) -> None:
        """Restore color settings widgets back to configuration defaults."""
        config = Config()  # Reload config to discard unsaved changes

        self.background_color_entry.delete(0, tk.END)
        self.background_color_entry.insert(0, str(config.background_color))

        # Reset MO color scheme dropdown
        predefined_schemes = ['bwr', 'RdBu', 'seismic', 'coolwarm', 'PiYG']

        if config.mo.color_scheme not in predefined_schemes and config.mo.color_scheme != 'custom':
            color_schemes = [config.mo.color_scheme, *predefined_schemes, 'custom']
            self.mo_color_scheme_dropdown['values'] = color_schemes
            self.mo_color_scheme_var.set(config.mo.color_scheme)
        else:
            color_schemes = [*predefined_schemes, 'custom']
            self.mo_color_scheme_dropdown['values'] = color_schemes
            if config.mo.color_scheme in predefined_schemes:
                self.mo_color_scheme_var.set(config.mo.color_scheme)
            else:
                self.mo_color_scheme_var.set('custom')

        # Reset custom color entries
        self.mo_negative_color_entry.delete(0, tk.END)
        self.mo_positive_color_entry.delete(0, tk.END)
        if config.mo.custom_colors:
            if len(config.mo.custom_colors) > 0:
                self.mo_negative_color_entry.insert(0, config.mo.custom_colors[0])
            if len(config.mo.custom_colors) > 1:
                self.mo_positive_color_entry.insert(0, config.mo.custom_colors[1])

        if self.mo_color_scheme_var.get() == 'custom':
            for widget in self.mo_custom_color_widgets:
                widget.grid()
        else:
            for widget in self.mo_custom_color_widgets:
                widget.grid_remove()

        # Reset bond color type
        self.bond_color_type_var.set(config.molecule.bond.color_type)

        # Reset bond color entry
        self.bond_color_entry.delete(0, tk.END)
        self.bond_color_entry.insert(0, str(config.molecule.bond.color))

        # Show/hide bond color entry based on type
        if self.bond_color_type_var.get() == 'uniform':
            self.bond_color_label.grid()
            self.bond_color_entry.grid()
        else:
            self.bond_color_label.grid_remove()
            self.bond_color_entry.grid_remove()

        self.apply_background_color()  # Reapply background color with new value
        self.apply_color_settings()  # Reapply MO and bond color settings with new values

    def apply_grid_settings(self) -> None:
        """Validate UI inputs and apply the chosen grid parameters."""
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

            new_grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
            if not np.array_equal(new_grid, self.tabulator.grid):
                self.update_mesh(r, theta, phi, GridType.SPHERICAL)

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

            new_grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
            if not np.array_equal(new_grid, self.tabulator.grid):
                self.update_mesh(x, y, z, GridType.CARTESIAN)

        # Replot the current orbital with the new grid
        idx = self._get_current_mo_index()
        if idx >= 0:
            self.plot_orbital(idx)

    def apply_molecule_settings(self) -> None:
        """Validate UI inputs and apply the chosen molecule rendering parameters."""
        redraw_molecule = False

        try:
            bond_max_length = float(self.bond_max_length_entry.get())
            if bond_max_length != config.molecule.bond.max_length:
                config.molecule.bond.max_length = bond_max_length
                redraw_molecule = True
        except ValueError:
            messagebox.showerror('Invalid Input', 'Bond Max Length must be a valid number.')
            return

        try:
            bond_radius = float(self.bond_radius_entry.get())
            if bond_radius != config.molecule.bond.radius:
                config.molecule.bond.radius = bond_radius
                redraw_molecule = True
        except ValueError:
            messagebox.showerror('Invalid Input', 'Bond Radius must be a valid number.')
            return

        if redraw_molecule:
            self.load_molecule(config)

    def apply_color_settings(self) -> None:
        """Apply both MO and bond color settings."""
        self.apply_mo_color_settings()
        self.apply_custom_mo_color_settings()
        self.apply_bond_color_settings()

    def apply_mo_color_settings(self) -> None:
        """Validate UI inputs and apply the chosen MO color settings."""
        self.on_mo_color_scheme_change(tk.Event())  # Update visibility of custom color entries

        mo_color_scheme = self.mo_color_scheme_var.get().strip()
        if mo_color_scheme == 'custom':
            return

        if mo_color_scheme != config.mo.color_scheme:
            self.cmap = mo_color_scheme
            config.mo.color_scheme = mo_color_scheme
            config.mo.custom_colors = []

            idx = self._get_current_mo_index()
            if idx >= 0:
                self.plot_orbital(idx)

    def apply_custom_mo_color_settings(self) -> None:
        """Validate UI inputs and apply the chosen MO color settings."""
        if self.mo_color_scheme_var.get().strip() != 'custom':
            return

        custom_colors = [self.mo_negative_color_entry.get().strip(), self.mo_positive_color_entry.get().strip()]

        if any(not mcolors.is_color_like(c) for c in custom_colors):
            messagebox.showerror('Invalid Input', 'One or more custom colors are not valid.')
            return

        if custom_colors != config.mo.custom_colors:
            config.mo.custom_colors = custom_colors
            config.mo.color_scheme = 'custom'
            self.cmap = self.custom_cmap_from_colors(custom_colors)

            idx = self._get_current_mo_index()
            if idx >= 0:
                self.plot_orbital(idx)

    def apply_bond_color_settings(self) -> None:
        """Validate UI inputs and apply the chosen color settings."""
        redraw_molecule = False

        bond_color_type = self.bond_color_type_var.get().strip()
        if bond_color_type != config.molecule.bond.color_type:
            config.molecule.bond.color_type = bond_color_type
            redraw_molecule = True

        bond_color = self.bond_color_entry.get().strip()
        if self.bond_color_type_var.get() == 'uniform' and bond_color != config.molecule.bond.color:
            config.molecule.bond.color = bond_color
            redraw_molecule = True

        if redraw_molecule:
            self.load_molecule(config)

    @staticmethod
    def save_settings() -> None:
        """Save current configuration to the user's custom config file."""
        try:
            config.save_current_config()
            messagebox.showinfo('Settings Saved', 'Configuration saved successfully to ~/.config/moldenViz/config.toml')
        except (OSError, ValueError) as e:
            messagebox.showerror('Save Error', f'Failed to save configuration: {e!s}')

    def plot_orbital(self, orb_ind: int) -> None:
        """Render the selected orbital isosurface in the PyVista plotter."""
        if self.orb_actor:
            self.pv_plotter.remove_actor(self.orb_actor)
            self.orb_actor = None
        if self.selection_screen:
            self.selection_screen.current_mo_ind = orb_ind

        if orb_ind == -1:
            if self.selection_screen:
                self.selection_screen.update_nav_button_states()
            return

        self.orb_mesh['orbital'] = self.tabulator.tabulate_mos(orb_ind)

        contour_mesh = self.orb_mesh.contour([-self.contour, self.contour])

        self.orb_actor = self.pv_plotter.add_mesh(
            contour_mesh,
            clim=[-self.contour, self.contour],
            opacity=self.opacity,
            show_scalar_bar=False,
            cmap=self.cmap,
            smooth_shading=True,
        )
        if self.selection_screen:
            self.selection_screen.update_nav_button_states()

    def _connect_pv_plotter_close_signal(self) -> None:
        """Connect the PyVista plotter close signal to handle closing both windows."""

        def on_pv_plotter_close() -> None:
            """Handle PyVista plotter close event by closing the selection screen and quitting."""
            if self.on_screen:
                self.on_screen = False
                if self.selection_screen and self.selection_screen.winfo_exists():
                    self.selection_screen.destroy()
                if self.tk_root and self._no_prev_tk_root:
                    self.tk_root.quit()

        self.pv_plotter.app_window.signal_close.connect(on_pv_plotter_close)

    def _clear_all(self) -> None:
        """Clear all actors from the plotter, including molecule and orbitals."""
        if self.molecule_actors:
            for actor in self.molecule_actors:
                actor.SetVisibility(False)

        if self.orb_actor:
            self.pv_plotter.remove_actor(self.orb_actor)
            self.orb_actor = None
            if self.selection_screen:
                self.selection_screen.current_mo_ind = -1
                self.selection_screen.update_nav_button_states()

    def toggle_molecule(self) -> None:
        """Toggle the visibility of the molecule."""
        if not self.molecule_actors:
            return

        if self.are_bonds_visible() != self.are_atoms_visible():
            if self.are_bonds_visible():
                self.toggle_atoms()
            else:
                self.toggle_bonds()
        else:
            for actor in self.molecule_actors:
                actor.SetVisibility(not actor.GetVisibility())
            self.pv_plotter.update()

        self.update_settings_button_states()

    def toggle_atoms(self) -> None:
        """Toggle the visibility of the molecule."""
        if self.atom_actors:
            for actor in self.atom_actors:
                actor.SetVisibility(not actor.GetVisibility())
            self.pv_plotter.update()

    def toggle_bonds(self) -> None:
        """Toggle the visibility of the molecule."""
        if self.bond_actors:
            for actor in self.bond_actors:
                actor.SetVisibility(not actor.GetVisibility())
            self.pv_plotter.update()

    def is_molecule_visible(self) -> bool:
        """Check if the molecule is currently visible in the plotter.

        Returns
        -------
        bool
            `True` if the molecule is visible, `False` otherwise.
        """
        if self.molecule_actors:
            return bool(self.molecule_actors[0].GetVisibility())  # Check visibility of the first actor
        return False

    def are_atoms_visible(self) -> bool:
        """Check if the atoms are currently visible in the plotter.

        Returns
        -------
        bool
            `True` if the atoms are visible, `False` otherwise.
        """
        if self.atom_actors:
            return bool(self.atom_actors[0].GetVisibility())  # Check visibility of the first actor
        return False

    def are_bonds_visible(self) -> bool:
        """Check if the bonds are currently visible in the plotter.

        Returns
        -------
        bool
            `True` if the bonds are visible, `False` otherwise.
        """
        if self.bond_actors:
            return bool(self.bond_actors[0].GetVisibility())  # Check visibility of the first actor
        return False

    def _create_mo_mesh(self) -> pv.StructuredGrid:
        """Create a mesh for the orbitals.

        Returns
        -------
            pv.StructuredGrid:
                The mesh object for MO visualization.

        """
        mesh = pv.StructuredGrid()
        mesh.points = pv.pyvista_ndarray(self.tabulator.grid)

        # Pyvista needs the dimensions backwards
        # in other words, (phi, theta, r) or (z, y, x)
        mesh.dimensions = self.tabulator._grid_dimensions[::-1]  # noqa: SLF001

        return mesh

    def update_mesh(
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
        if grid_type == GridType.CARTESIAN:
            self.tabulator.cartesian_grid(i_points, j_points, k_points)
        elif grid_type == GridType.SPHERICAL:
            self.tabulator.spherical_grid(i_points, j_points, k_points)
        else:
            raise ValueError('The plotter only supports spherical and cartesian grids.')

        self.orb_mesh = self._create_mo_mesh()


class _OrbitalSelectionScreen(tk.Toplevel):
    """Modal dialog that lets users browse and configure molecular orbitals."""

    SPHERICAL_GRID_SETTINGS_WINDOW_SIZE = '400x350'
    CARTESIAN_GRID_SETTINGS_WINDOW_SIZE = '650x400'

    def __init__(self, plotter: Plotter) -> None:
        """Create the orbital selection dialog for a plotter instance.

        Parameters
        ----------
        plotter : Plotter
            Active plotter that supplies molecular orbital data.
        tk_master : tk.Tk
            Tk root or parent window that owns this dialog.
        """
        super().__init__(plotter.tk_root)
        self.title('Orbitals')
        self.geometry('350x500')

        self.protocols()

        self.plotter = plotter
        self.current_mo_ind = -1  # Start with no orbital shown

        # Initialize export window attributes
        self._export_window = None
        self._export_current_orb_radio = None
        self._export_all_orb_radio = None

        nav_frame = ttk.Frame(self)
        nav_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.prev_button = ttk.Button(nav_frame, text='<< Previous', command=self.prev_plot)
        self.prev_button.pack(side=tk.LEFT, padx=5, pady=10)

        self.next_button = ttk.Button(nav_frame, text='Next >>', command=self.next_plot)
        self.next_button.pack(side=tk.RIGHT, padx=5, pady=10)

        self.update_nav_button_states()  # Update buttons for initial state

        self.orb_tv = _OrbitalsTreeview(self)
        self.orb_tv.populate_tree(self.plotter.tabulator._parser.mos)  # noqa: SLF001
        self.orb_tv.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    def on_close(self) -> None:
        """Close the selection dialog and release GUI resources."""
        self.plotter.on_screen = False
        self.plotter.pv_plotter.close()
        self.destroy()
        if self.plotter.tk_root and self.plotter._no_prev_tk_root:  # noqa: SLF001
            self.plotter.tk_root.quit()
            self.plotter.tk_root.destroy()

    def protocols(self) -> None:
        """Attach standard close shortcuts to the dialog window."""
        self.protocol('WM_DELETE_WINDOW', self.on_close)
        self.bind('<Command-q>', lambda _event: self.on_close())
        self.bind('<Command-w>', lambda _event: self.on_close())
        self.bind('<Control-q>', lambda _event: self.on_close())
        self.bind('<Control-w>', lambda _event: self.on_close())

    def next_plot(self) -> None:
        """Advance to the next molecular orbital."""
        max_index = len(self.plotter.tabulator._parser.mos) - 1  # noqa: SLF001
        if max_index < 0:
            return
        current = self.current_mo_ind
        new_index = 0 if current < 0 else min(current + 1, max_index)
        self.plotter.plot_orbital(new_index)
        if self.current_mo_ind >= 0:
            self.orb_tv.highlight_orbital(self.current_mo_ind)

    def prev_plot(self) -> None:
        """Return to the previous molecular orbital."""
        if self.current_mo_ind <= 0:
            return
        new_index = self.current_mo_ind - 1
        self.plotter.plot_orbital(new_index)
        self.orb_tv.highlight_orbital(self.current_mo_ind)

    def update_nav_button_states(self) -> None:
        """Synchronize navigation button state with the current orbital index."""
        total = len(self.plotter.tabulator._parser.mos)  # noqa: SLF001
        can_go_prev = self.current_mo_ind > 0
        can_go_next = total > 0 and self.current_mo_ind < total - 1
        self.prev_button.config(state=tk.NORMAL if can_go_prev else tk.DISABLED)
        self.next_button.config(state=tk.NORMAL if can_go_next else tk.DISABLED)
        self._update_export_dialog_label()

    def _update_export_dialog_label(self) -> None:
        """Update the export dialog label to reflect the current orbital index."""
        if self._export_current_orb_radio is not None:
            # Use 1-based indexing for display (add 1 to current_mo_ind)
            orbital_display = self.current_mo_ind + 1 if self.current_mo_ind >= 0 else 'None'
            self._export_current_orb_radio.config(text=f'Current orbital (#{orbital_display})')
            # Update the state based on whether an orbital is selected
            if self.current_mo_ind < 0:
                self._export_current_orb_radio.config(state=tk.DISABLED)
            else:
                self._export_current_orb_radio.config(state=tk.NORMAL)

    def plot_orbital(self, orb_ind: int) -> None:
        """Render the selected orbital isosurface in the PyVista plotter.

        Parameters
        ----------
        orb_ind : int
            Index of the orbital to display; ``-1`` clears the current mesh.
        """
        self.plotter.plot_orbital(orb_ind)
        self.current_mo_ind = orb_ind


class _OrbitalsTreeview(ttk.Treeview):
    def __init__(self, selection_screen: _OrbitalSelectionScreen) -> None:
        """Initialise the tree view that lists available molecular orbitals.

        Parameters
        ----------
        selection_screen : _OrbitalSelectionScreen
            Parent dialog that handles selection changes.
        """
        columns = ['Index', 'Symmetry', 'Occupation', 'Energy [au]']
        widths = [20, 50, 50, 120]

        super().__init__(selection_screen, columns=columns, show='headings', height=20)

        for col, w in zip(columns, widths, strict=False):
            self.heading(col, text=col)
            self.column(col, width=w)

        self.selection_screen = selection_screen

        self.current_mo_ind = -1  # Start with no orbital shown

        # Configure tag
        self.tag_configure('highlight', background='lightblue')

        self.bind('<<TreeviewSelect>>', self.on_select)

    def highlight_orbital(self, orb_ind: int) -> None:
        """Highlight the given orbital within the tree view.

        Parameters
        ----------
        orb_ind : int
            Index to highlight.
        """
        if self.current_mo_ind != -1:
            self.item(self.current_mo_ind, tags=('!hightlight',))

        self.current_mo_ind = orb_ind
        self.item(orb_ind, tags=('highlight',))
        self.see(orb_ind)  # Scroll to the selected item

    def erase(self) -> None:
        """Remove all orbital entries from the tree view."""
        for item in self.get_children():
            self.delete(item)

    def populate_tree(self, mos: list[_MolecularOrbital]) -> None:
        """Populate the tree view with molecular orbital metadata.

        Parameters
        ----------
        mos : list[_MolecularOrbital]
            Orbitals sourced from the parser.
        """
        self.erase()

        # Counts the number of MOs with a given symmetry
        mo_syms = list({mo.sym for mo in mos})
        mo_sym_count: dict[str, int] = dict.fromkeys(mo_syms, 0)
        for ind, mo in enumerate(mos):
            mo_sym_count[mo.sym] += 1
            self.insert('', 'end', iid=ind, values=(ind + 1, f'{mo.sym}.{mo_sym_count[mo.sym]}', mo.occ, mo.energy))

    def on_select(self, _event: tk.Event) -> None:
        """Handle user selection events raised by the tree view.

        Parameters
        ----------
        _event : tk.Event
            Tkinter event object (unused).
        """
        selected_item = self.selection()
        self.selection_remove(selected_item)
        if selected_item:
            orb_ind = int(selected_item[0])
            self.highlight_orbital(orb_ind)
            self.selection_screen.current_mo_ind = orb_ind
            self.selection_screen.plot_orbital(orb_ind)
            self.selection_screen.update_nav_button_states()
