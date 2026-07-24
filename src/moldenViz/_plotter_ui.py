"""Private Tk and Qt user-interface helpers for :mod:`moldenViz.plotter`."""

from __future__ import annotations

import logging
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import TYPE_CHECKING, Any, cast

import matplotlib.colors as mcolors
import numpy as np
from qtpy.QtWidgets import QAction, QMenu  # pyright: ignore[reportPrivateImportUsage]
from shiboken6 import isValid

from ._config_module import Config
from .tabulator import GridType, Tabulator

if TYPE_CHECKING:
    from .models import MolecularOrbital
    from .plotter import Plotter


logger = logging.getLogger(__name__)

_MO_COLOR_SCHEMES = ['bwr', 'RdBu', 'seismic', 'coolwarm', 'PiYG']


def _mo_color_scheme_options(current_config: Config) -> tuple[list[str], str]:
    """Return dropdown options and the active MO color scheme.

    Parameters
    ----------
    current_config : Config
        Configuration containing the saved molecular-orbital colors.

    Returns
    -------
    tuple[list[str], str]
        Available dropdown values and the value that should be selected.
    """
    options = [*_MO_COLOR_SCHEMES, 'custom']
    if current_config.mo.custom_colors:
        return options, 'custom'
    if current_config.mo.color_scheme not in _MO_COLOR_SCHEMES:
        return [current_config.mo.color_scheme, *options], current_config.mo.color_scheme
    return options, current_config.mo.color_scheme


def _plotter_config() -> Config:
    """Return the plotter module's active configuration object.

    Returns
    -------
    Config
        Configuration singleton used by the active plotter module.
    """
    from .plotter import config  # ruff:ignore[import-outside-top-level]

    return config


class _ConfigProxy:
    """Forward UI configuration access to the plotter's active config."""

    def __getattr__(self, name: str) -> object:
        return getattr(_plotter_config(), name)

    def __setattr__(self, name: str, value: object) -> None:
        setattr(_plotter_config(), name, value)


config = cast(Config, _ConfigProxy())


class _PlotterUI:
    """Mixin implementing menus, dialogs, and settings for Plotter."""

    if TYPE_CHECKING:
        _CARTESIAN_GRID_SETTINGS_WINDOW_SIZE: str
        _SPHERICAL_GRID_SETTINGS_WINDOW_SIZE: str
        _molecule: Any
        _molecule_actors: list[Any]
        _only_molecule: bool
        _orb_actor: Any | None
        _pv_plotter: Any
        _selection_screen: Any | None
        tabulator: Tabulator
        _tk_root: Any

        def _clear_all(self) -> None: ...
        def _ensure_gtos_ready(self) -> bool: ...
        def are_atoms_visible(self) -> bool: ...
        def are_bonds_visible(self) -> bool: ...
        def _custom_cmap_from_colors(self, colors: list[str]) -> Any: ...
        def _load_molecule(self, config: Config) -> None: ...
        def plot_orbital(self, orb_ind: int) -> None: ...
        def toggle_atoms(self) -> None: ...
        def toggle_bonds(self) -> None: ...
        def toggle_molecule(self) -> None: ...
        def _update_mesh(self, i_points: Any, j_points: Any, k_points: Any, grid_type: GridType) -> None: ...

    def _override_clear_all_button(self) -> None:
        """Override the default "Clear All" action in the PyVista plotter's View menu."""
        view_menu = None
        main_menu = self._pv_plotter.main_menu
        if main_menu is None:
            raise RuntimeError('PyVista plotter does not have a main menu.')

        for action in main_menu.actions():
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
        settings_menu = QMenu('Settings', self._pv_plotter.app_window)

        # Add Settings submenu items
        if not self._only_molecule:
            grid_settings_action = QAction('Grid Settings', self._pv_plotter.app_window)
            grid_settings_action.triggered.connect(self._grid_settings_screen)
            settings_menu.addAction(grid_settings_action)

            mo_settings_action = QAction('MO Settings', self._pv_plotter.app_window)
            mo_settings_action.triggered.connect(self._mo_settings_screen)
            settings_menu.addAction(mo_settings_action)

        molecule_settings_action = QAction('Molecule Settings', self._pv_plotter.app_window)
        molecule_settings_action.triggered.connect(self._molecule_settings_screen)
        settings_menu.addAction(molecule_settings_action)

        color_settings_action = QAction('Color Settings', self._pv_plotter.app_window)
        color_settings_action.triggered.connect(self._color_settings_screen)
        settings_menu.addAction(color_settings_action)

        settings_menu.addSeparator()

        save_settings_action = QAction('Save Settings', self._pv_plotter.app_window)
        save_settings_action.triggered.connect(self._save_settings)
        settings_menu.addAction(save_settings_action)

        # Create Export menu with dropdown
        export_menu = QMenu('Export', self._pv_plotter.app_window)

        # Add Export submenu items
        if not self._only_molecule:
            export_data_action = QAction('Data', self._pv_plotter.app_window)
            export_data_action.triggered.connect(self._export_orbitals_dialog)
            export_menu.addAction(export_data_action)

        export_image_action = QAction('Image', self._pv_plotter.app_window)
        export_image_action.triggered.connect(self._export_image_dialog)
        export_menu.addAction(export_image_action)

        # Add menus to main menu bar
        main_menu = self._pv_plotter.main_menu
        if main_menu is None:
            raise RuntimeError('PyVista plotter does not have a main menu.')
        main_menu.addMenu(settings_menu)
        main_menu.addMenu(export_menu)

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
        assert self._selection_screen is not None

        file_format = format_var.get()
        scope = scope_var.get()
        logger.info('Export requested: format=%s, scope=%s.', file_format, scope)

        # Validate selection
        if scope == 'current' and self._selection_screen.current_mo_ind < 0:
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
            f'orbitals_all{ext}' if scope == 'all' else f'orbital_{self._selection_screen.current_mo_ind}{ext}'
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
            mo_index = self._selection_screen.current_mo_ind if scope == 'current' else None
            self.tabulator.export(file_path, mo_index=mo_index)
            messagebox.showinfo('Export Successful', f'Orbital(s) exported successfully to:\n{file_path}')
            logger.info('Export completed successfully to %s.', file_path)
            export_window.destroy()
        except (RuntimeError, ValueError) as e:
            logger.exception('Export failed during orbital export.')
            messagebox.showerror('Export Failed', f'Failed to export orbital(s):\n\n{e!s}')

    def _export_orbitals_dialog(self) -> None:
        """Open a dialog to configure and export molecular orbitals."""
        assert self._selection_screen is not None
        if not self._ensure_gtos_ready():
            return

        export_window = tk.Toplevel(self._tk_root)
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
            self._selection_screen.current_mo_ind + 1 if self._selection_screen.current_mo_ind >= 0 else 'None'
        )
        current_orb_radio = ttk.Radiobutton(
            main_frame,
            text=f'Current orbital (#{orbital_display})',
            variable=scope_var,
            value='current',
        )
        current_orb_radio.grid(row=4, column=0, columnspan=2, sticky=tk.W, padx=20)
        if self._selection_screen.current_mo_ind < 0:
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
        def _on_close() -> None:
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
        ttk.Button(button_frame, text='Cancel', command=_on_close).pack(side=tk.LEFT, padx=5)
        export_window.protocol('WM_DELETE_WINDOW', _on_close)

    def _export_image_dialog(self) -> None:
        """Open a dialog to export the current visualization as an image."""
        export_window = tk.Toplevel(self._tk_root)
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
            self._save_image(file_path, file_format, transparent)
            messagebox.showinfo('Export Successful', f'Image exported successfully to:\n{file_path}')
            export_window.destroy()
        except (RuntimeError, OSError, ValueError) as e:
            messagebox.showerror('Export Failed', f'Failed to export image:\n\n{e!s}')

    def _save_image(self, file_path: str, file_format: str, transparent: bool) -> None:
        """Save an image using the exporter for its file format.

        Parameters
        ----------
        file_path : str
            Destination path for the exported image.
        file_format : str
            Selected image format.
        transparent : bool
            Whether PNG output should use a transparent background.
        """
        if file_format in {'svg', 'pdf'}:
            self._pv_plotter.save_graphic(file_path)
        else:
            self._pv_plotter.screenshot(
                file_path,
                transparent_background=transparent if file_format == 'png' else False,
            )

    def _settings_parent(self) -> tk.Misc:
        """Return the appropriate parent widget for settings dialogs.

        Returns
        -------
        tk.Misc
            The parent widget for settings dialogs.
        """
        parent = self._selection_screen if self._selection_screen is not None else self._tk_root
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
        if self._selection_screen:
            return self._selection_screen.current_mo_ind
        return -1

    def _grid_settings_screen(self) -> None:
        """Open the grid settings window."""
        parent = self._settings_parent()
        self.grid_settings_window = tk.Toplevel(parent)
        self.grid_settings_window.title('Grid Settings')

        settings_frame = ttk.Frame(self.grid_settings_window)
        settings_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Grid parameters
        ttk.Label(settings_frame, text='MO Grid parameters').grid(row=0, column=0, padx=5, pady=5, columnspan=5)

        self.grid_type_radio_var = tk.StringVar()
        self.grid_type_radio_var.set(self.tabulator.grid_type.value)

        ttk.Label(settings_frame, text='Spherical grid:').grid(row=1, column=0, padx=5, pady=5)
        sph_grid_type_button = ttk.Radiobutton(
            settings_frame,
            variable=self.grid_type_radio_var,
            value=GridType.SPHERICAL.value,
            command=self._place_grid_params_frame,
        )

        ttk.Label(settings_frame, text='Cartesian grid:').grid(row=1, column=2, padx=5, pady=5)
        cart_grid_type_button = ttk.Radiobutton(
            settings_frame,
            variable=self.grid_type_radio_var,
            value=GridType.CARTESIAN.value,
            command=self._place_grid_params_frame,
        )

        sph_grid_type_button.grid(row=1, column=1, padx=5, pady=5)
        cart_grid_type_button.grid(row=1, column=3, padx=5, pady=5)

        self.sph_grid_params_frame = self._sph_grid_params_frame_widgets(settings_frame)
        self.cart_grid_params_frame = self._cart_grid_params_frame_widgets(settings_frame)

        self._place_grid_params_frame()

        # Reset button
        reset_button = ttk.Button(settings_frame, text='Reset', command=self._reset_grid_settings)
        reset_button.grid(row=8, column=0, padx=5, pady=5, columnspan=5)

        # Apply settings button
        apply_button = ttk.Button(settings_frame, text='Apply', command=self._apply_grid_settings)
        apply_button.grid(row=9, column=0, padx=5, pady=5, columnspan=5)

    def _mo_settings_screen(self) -> None:
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
        self.contour_entry.insert(0, str(self._contour))
        self.contour_entry.grid(row=1, column=0, padx=5, pady=5, sticky='ew')

        # Bind to apply changes on Enter key or focus out
        self.contour_entry.bind('<Return>', lambda _e: self._apply_mo_contour())
        self.contour_entry.bind('<FocusOut>', lambda _e: self._apply_mo_contour())

        # Opacity
        opacity_label = ttk.Label(settings_frame)
        opacity_label.grid(row=2, column=0, padx=5, pady=5, sticky='w')
        self.opacity_scale = ttk.Scale(
            settings_frame,
            length=200,
            command=self._on_opacity_change,
        )
        self.opacity_scale.set(self._opacity)
        self.opacity_scale.grid(row=3, column=0, padx=5, pady=5, sticky='ew')
        # Initialize label
        opacity_label.config(text=f'Molecular Orbital Opacity: {self._opacity:.2f}')

        # Configure grid column weight for proper resizing
        settings_frame.columnconfigure(0, weight=1)

        # Reset button
        reset_button = ttk.Button(settings_frame, text='Reset', command=self._reset_mo_settings)
        reset_button.grid(row=4, column=0, padx=5, pady=5, sticky='ew')

    def _on_opacity_change(self, val: str) -> None:
        """Handle opacity slider changes and apply immediately."""
        opacity = round(float(val), 2)
        self._opacity = opacity
        if self._orb_actor:
            self._orb_actor.GetProperty().SetOpacity(opacity)

        # Update label
        for widget in self.mo_settings_window.winfo_children():
            for child in widget.winfo_children():
                if isinstance(child, ttk.Label) and 'Molecular Orbital Opacity:' in child.cget('text'):
                    child.config(text=f'Molecular Orbital Opacity: {opacity:.2f}')
        logger.info('Set molecular orbital opacity to %.2f.', opacity)

    def _apply_mo_contour(self) -> None:
        """Apply contour changes immediately."""
        try:
            self._contour = float(self.contour_entry.get().strip())
            logger.info('Set molecular orbital contour to %.2f.', self._contour)
            # Replot the current orbital with the new contour
            idx = self._get_current_mo_index()
            if idx >= 0:
                self.plot_orbital(idx)
        except ValueError:
            pass  # Ignore invalid input

    def _update_settings_button_states(self) -> None:
        """Update the state of the settings buttons based on current plotter state."""
        if hasattr(self, 'show_atoms_var'):
            self.show_atoms_var.set(self.are_atoms_visible())
        if hasattr(self, 'show_bonds_var'):
            self.show_bonds_var.set(self.are_bonds_visible())

        config.molecule.atom.show = self.are_atoms_visible()
        config.molecule.bond.show = self.are_bonds_visible()

    def _molecule_settings_screen(self) -> None:
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
            command=self._on_molecule_opacity_change,
        )
        self.molecule_opacity_scale.set(self._molecule_opacity)
        self.molecule_opacity_scale.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky='ew')
        # Initialize label
        molecule_opacity_label.config(text=f'Molecule Opacity: {self._molecule_opacity:.2f}')

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
        reset_button = ttk.Button(settings_frame, text='Reset', command=self._reset_molecule_settings)
        reset_button.grid(row=13, column=0, columnspan=2, padx=5, pady=5, sticky='ew')

        # Apply settings button
        apply_button = ttk.Button(settings_frame, text='Apply', command=self._apply_molecule_settings)
        apply_button.grid(row=14, column=0, columnspan=2, padx=5, pady=5, sticky='ew')

    def _on_molecule_opacity_change(self, val: str) -> None:
        """Handle molecule opacity slider changes and apply immediately."""
        opacity = round(float(val), 2)
        self._molecule_opacity = opacity
        for actor in self._molecule_actors:
            actor.GetProperty().SetOpacity(opacity)
        # Update label
        for widget in self.molecule_settings_window.winfo_children():
            for child in widget.winfo_children():
                if isinstance(child, ttk.Label) and 'Molecule Opacity:' in child.cget('text'):
                    child.config(text=f'Molecule Opacity: {opacity:.2f}')
        logger.info('Set molecule opacity to %.2f.', opacity)

    def _apply_background_color(self) -> None:
        """Apply background color changes immediately."""
        try:
            color = self.background_color_entry.get().strip()
            self._set_background_color(color)
        except (ValueError, RuntimeError) as e:
            messagebox.showerror('Error', f'Failed to set background color: {e!s}')

    def _set_background_color(self, color: str) -> None:
        """Set the background color when the value is valid."""
        if mcolors.is_color_like(color):
            self._pv_plotter.set_background(color)
            logger.info('Set background color to %s.', color)
        else:
            messagebox.showerror('Invalid Input', f'"{color}" is not a valid color.')

    def _on_mo_color_scheme_change(self, _event: tk.Event) -> None:
        """Handle MO color scheme dropdown change to show/hide custom color entries."""
        if self.mo_color_scheme_var.get() == 'custom':
            for widget in self.mo_custom_color_widgets:
                widget.grid()
        else:
            for widget in self.mo_custom_color_widgets:
                widget.grid_remove()

    def _on_bond_color_type_change(self) -> None:
        """Handle bond color type change to show/hide bond color entry."""
        if self.bond_color_type_var.get() == 'uniform':
            self.bond_color_label.grid()
            self.bond_color_entry.grid()
        else:
            self.bond_color_label.grid_remove()
            self.bond_color_entry.grid_remove()

        self._apply_bond_color_settings()

    def _color_settings_screen(self) -> None:
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
        self.background_color_entry.bind('<Return>', lambda _e: self._apply_background_color())
        self.background_color_entry.bind('<FocusOut>', lambda _e: self._apply_background_color())

        # Separator
        ttk.Separator(settings_frame, orient='horizontal').grid(row=2, column=0, columnspan=2, sticky='ew', pady=10)

        if not self._only_molecule:
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
            color_schemes, default_scheme = _mo_color_scheme_options(config)

            self.mo_color_scheme_var = tk.StringVar(value=default_scheme)
            self.mo_color_scheme_dropdown = ttk.Combobox(
                settings_frame,
                state='readonly',
                textvariable=self.mo_color_scheme_var,
                values=color_schemes,
            )
            self.mo_color_scheme_dropdown.grid(row=4, column=1, padx=5, pady=5, sticky='w')
            self.mo_color_scheme_dropdown.bind('<<ComboboxSelected>>', self._on_mo_color_scheme_change)
            self.mo_color_scheme_dropdown.bind('<<ComboboxSelected>>', lambda _e: self._apply_mo_color_settings())

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
            command=self._on_bond_color_type_change,
        ).pack(side=tk.LEFT)
        ttk.Radiobutton(
            bond_color_frame,
            text='Split',
            variable=self.bond_color_type_var,
            value='split',
            command=self._on_bond_color_type_change,
        ).pack(side=tk.LEFT)

        # Bond color (for uniform type only)
        self.bond_color_label = ttk.Label(settings_frame, text='Bond Color:')
        self.bond_color_label.grid(row=10, column=0, padx=5, pady=5, sticky='w')
        self.bond_color_entry = ttk.Entry(settings_frame, width=15)
        self.bond_color_entry.insert(0, str(config.molecule.bond.color))
        self.bond_color_entry.grid(row=10, column=1, padx=5, pady=5, sticky='w')
        self.bond_color_entry.bind('<Return>', lambda _e: self._apply_bond_color_settings())

        # Hide bond color entry if split is selected
        if self.bond_color_type_var.get() == 'split':
            self.bond_color_label.grid_remove()
            self.bond_color_entry.grid_remove()

        # Configure grid column weights for proper resizing
        settings_frame.columnconfigure(0, weight=1)
        settings_frame.columnconfigure(1, weight=1)

        # Reset button
        reset_button = ttk.Button(settings_frame, text='Reset', command=self._reset_color_settings)
        reset_button.grid(row=11, column=0, columnspan=2, padx=5, pady=5, sticky='ew')

        # Apply settings button
        apply_button = ttk.Button(settings_frame, text='Apply', command=self._apply_color_settings)
        apply_button.grid(row=13, column=0, columnspan=2, padx=5, pady=5, sticky='ew')

    def _place_grid_params_frame(self) -> None:
        """Render the parameter frame that matches the selected grid type."""
        if self.grid_type_radio_var.get() == GridType.SPHERICAL.value:
            self.grid_settings_window.geometry(self._SPHERICAL_GRID_SETTINGS_WINDOW_SIZE)
            self.cart_grid_params_frame.grid_forget()
            self.sph_grid_params_frame.grid(row=2, column=0, padx=5, pady=5, rowspan=6, columnspan=4)
            self._sph_grid_params_frame_setup()
        else:
            self.grid_settings_window.geometry(self._CARTESIAN_GRID_SETTINGS_WINDOW_SIZE)
            self.sph_grid_params_frame.grid_forget()
            self.cart_grid_params_frame.grid(row=2, column=0, padx=5, pady=5, rowspan=6, columnspan=4)
            self._cart_grid_params_frame_setup()

    def _sph_grid_params_frame_widgets(self, master: ttk.Frame) -> ttk.Frame:
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

    def _cart_grid_params_frame_widgets(self, master: ttk.Frame) -> ttk.Frame:
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

    def _sph_grid_params_frame_setup(self) -> None:
        """Populate the spherical grid widgets with defaults or existing values."""
        self.radius_entry.delete(0, tk.END)
        self.radius_points_entry.delete(0, tk.END)
        self.theta_points_entry.delete(0, tk.END)
        self.phi_points_entry.delete(0, tk.END)

        # Previous grid was cartesian, so use default values
        if self.tabulator.grid_type == GridType.CARTESIAN:
            self.radius_entry.insert(
                0,
                str(max(config.grid.max_radius_multiplier * self._molecule.max_radius, config.grid.min_radius)),
            )
            self.radius_points_entry.insert(0, str(config.grid.spherical.num_r_points))
            self.theta_points_entry.insert(0, str(config.grid.spherical.num_theta_points))
            self.phi_points_entry.insert(0, str(config.grid.spherical.num_phi_points))
            return

        num_r, num_theta, num_phi = self.tabulator.grid_dimensions

        # The last point of the grid for sure has the largest r
        r, _, _ = Tabulator.cartesian_to_spherical(*self.tabulator.grid[-1, :])

        self.radius_entry.insert(0, str(r))
        self.radius_points_entry.insert(0, str(num_r))
        self.theta_points_entry.insert(0, str(num_theta))
        self.phi_points_entry.insert(0, str(num_phi))

    def _cart_grid_params_frame_setup(self) -> None:
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
        if self.tabulator.grid_type == GridType.SPHERICAL:
            r = max(config.grid.max_radius_multiplier * self._molecule.max_radius, config.grid.min_radius)

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

        x_num, y_num, z_num = self.tabulator.grid_dimensions
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

    def _reset_grid_settings(self) -> None:
        """Restore grid settings widgets back to configuration defaults."""
        self.grid_type_radio_var.set(config.grid.default_type)
        self._place_grid_params_frame()

        self.radius_entry.delete(0, tk.END)
        self.radius_entry.insert(
            0,
            str(max(config.grid.max_radius_multiplier * self._molecule.max_radius, config.grid.min_radius)),
        )

        self.radius_points_entry.delete(0, tk.END)
        self.radius_points_entry.insert(0, str(config.grid.spherical.num_r_points))

        self.theta_points_entry.delete(0, tk.END)
        self.theta_points_entry.insert(0, str(config.grid.spherical.num_theta_points))

        self.phi_points_entry.delete(0, tk.END)
        self.phi_points_entry.insert(0, str(config.grid.spherical.num_phi_points))

    def _reset_mo_settings(self) -> None:
        """Restore MO settings widgets back to configuration defaults."""
        self.contour_entry.delete(0, tk.END)
        self.contour_entry.insert(0, str(config.mo.contour))

        self.opacity_scale.set(config.mo.opacity)

        self._apply_mo_contour()  # Reapply contour with new value

    def _reset_molecule_settings(self) -> None:
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

        self._apply_molecule_settings()  # Reapply molecule settings with new values

    def _reset_color_settings(self) -> None:
        """Restore color settings widgets back to configuration defaults."""
        config = Config()  # Reload config to discard unsaved changes

        self.background_color_entry.delete(0, tk.END)
        self.background_color_entry.insert(0, str(config.background_color))

        # Reset MO color scheme dropdown
        color_schemes, default_scheme = _mo_color_scheme_options(config)
        self.mo_color_scheme_dropdown['values'] = color_schemes
        self.mo_color_scheme_var.set(default_scheme)

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

        self._apply_background_color()  # Reapply background color with new value
        self._apply_color_settings()  # Reapply MO and bond color settings with new values

    def _apply_grid_settings(self) -> None:
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
            xx, yy, zz = Tabulator.spherical_to_cartesian(rr, tt, pp)

            new_grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
            if not np.array_equal(new_grid, self.tabulator.grid):
                logger.info(
                    'Applying spherical grid: radius=%.3f (r=%d theta=%d phi=%d points).',
                    radius,
                    num_r_points,
                    num_theta_points,
                    num_phi_points,
                )
                self._update_mesh(r, theta, phi, GridType.SPHERICAL)

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
                logger.info(
                    'Applying cartesian grid: x=[%.2f, %.2f] (%d pts), y=[%.2f, %.2f] (%d pts), '
                    'z=[%.2f, %.2f] (%d pts).',
                    x_min,
                    x_max,
                    x_num,
                    y_min,
                    y_max,
                    y_num,
                    z_min,
                    z_max,
                    z_num,
                )
                self._update_mesh(x, y, z, GridType.CARTESIAN)

    def _apply_molecule_settings(self) -> None:
        """Validate UI inputs and apply the chosen molecule rendering parameters."""
        changes: list[str] = []

        try:
            bond_max_length = float(self.bond_max_length_entry.get())
            if bond_max_length != config.molecule.bond.max_length:
                config.molecule.bond.max_length = bond_max_length
                changes.append(f'bond_max_length={bond_max_length:.2f}')
        except ValueError:
            messagebox.showerror('Invalid Input', 'Bond Max Length must be a valid number.')
            return

        try:
            bond_radius = float(self.bond_radius_entry.get())
            if bond_radius != config.molecule.bond.radius:
                config.molecule.bond.radius = bond_radius
                changes.append(f'bond_radius={bond_radius:.2f}')
        except ValueError:
            messagebox.showerror('Invalid Input', 'Bond Radius must be a valid number.')
            return

        if changes:
            logger.info('Applying molecule settings: %s.', ', '.join(changes))
            self._load_molecule(config)
        else:
            logger.debug('Molecule settings unchanged; skipping reload.')

    def _apply_color_settings(self) -> None:
        """Apply both MO and bond color settings."""
        self._apply_mo_color_settings()
        self._apply_custom_mo_color_settings()
        self._apply_bond_color_settings()

    def _apply_mo_color_settings(self) -> None:
        """Validate UI inputs and apply the chosen MO color settings."""
        self._on_mo_color_scheme_change(tk.Event())  # Update visibility of custom color entries

        mo_color_scheme = self.mo_color_scheme_var.get().strip()
        if mo_color_scheme == 'custom':
            return

        if mo_color_scheme != config.mo.color_scheme:
            self._cmap = mo_color_scheme
            config.mo.color_scheme = mo_color_scheme
            config.mo.custom_colors = None
            logger.info('Applied MO color scheme: %s.', mo_color_scheme)

            idx = self._get_current_mo_index()
            if idx >= 0:
                self.plot_orbital(idx)

    def _apply_custom_mo_color_settings(self) -> None:
        """Validate UI inputs and apply the chosen MO color settings."""
        if self.mo_color_scheme_var.get().strip() != 'custom':
            return

        custom_colors = [self.mo_negative_color_entry.get().strip(), self.mo_positive_color_entry.get().strip()]

        if any(not mcolors.is_color_like(c) for c in custom_colors):
            messagebox.showerror('Invalid Input', 'One or more custom colors are not valid.')
            return

        if custom_colors != config.mo.custom_colors:
            config.mo.custom_colors = custom_colors
            self._cmap = self._custom_cmap_from_colors(custom_colors)
            logger.info('Applied custom MO colors: negative=%s, positive=%s.', custom_colors[0], custom_colors[1])

            idx = self._get_current_mo_index()
            if idx >= 0:
                self.plot_orbital(idx)

    def _apply_bond_color_settings(self) -> None:
        """Validate UI inputs and apply the chosen color settings."""
        changes: list[str] = []

        bond_color_type = self.bond_color_type_var.get().strip()
        if bond_color_type != config.molecule.bond.color_type:
            config.molecule.bond.color_type = bond_color_type
            changes.append(f'color_type={bond_color_type}')

        bond_color = self.bond_color_entry.get().strip()
        if self.bond_color_type_var.get() == 'uniform' and bond_color != config.molecule.bond.color:
            config.molecule.bond.color = bond_color
            changes.append(f'color={bond_color}')

        if changes:
            logger.info('Applying bond color settings: %s.', ', '.join(changes))
            self._load_molecule(config)
        else:
            logger.debug('Bond color settings unchanged; skipping reload.')

    @staticmethod
    def _save_settings() -> None:
        """Save current configuration to the user's custom config file."""
        try:
            config._save_current_config()  # ruff:ignore[private-member-access]
            messagebox.showinfo('Settings Saved', 'Configuration saved successfully to ~/.config/moldenViz/config.toml')
        except (OSError, ValueError) as e:
            messagebox.showerror('Save Error', f'Failed to save configuration: {e!s}')


class _OrbitalSelectionScreen(tk.Toplevel):
    """Modal dialog that lets users browse and configure molecular orbitals."""

    _SPHERICAL_GRID_SETTINGS_WINDOW_SIZE = '400x350'
    _CARTESIAN_GRID_SETTINGS_WINDOW_SIZE = '650x400'

    def __init__(self, plotter: Plotter) -> None:
        """Create the orbital selection dialog for a plotter instance.

        Parameters
        ----------
        plotter : Plotter
            Active plotter that supplies molecular orbital data.
        """
        super().__init__(plotter._tk_root)  # ruff:ignore[private-member-access]
        self.title('Orbitals')
        self.geometry('350x500')

        self._protocols()

        self.plotter = plotter
        self.current_mo_ind = -1  # Start with no orbital shown
        self._loading = False
        self._loading_label = ttk.Label(self, text='Tabulating orbitals...', anchor=tk.W)

        # Initialize export window attributes
        self._export_window = None
        self._export_current_orb_radio = None
        self._export_all_orb_radio = None

        nav_frame = ttk.Frame(self)
        nav_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.prev_button = ttk.Button(nav_frame, text='<< Previous', command=self._prev_plot)
        self.prev_button.pack(side=tk.LEFT, padx=5, pady=10)

        self.next_button = ttk.Button(nav_frame, text='Next >>', command=self._next_plot)
        self.next_button.pack(side=tk.RIGHT, padx=5, pady=10)

        self._update_nav_button_states()  # Update buttons for initial state

        self.orb_tv = _OrbitalsTreeview(self)
        self.orb_tv._populate_tree(  # ruff:ignore[private-member-access]
            self.plotter.tabulator.molecular_orbitals,
        )
        self.orb_tv.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    def _on_close(self) -> None:
        """Close the selection dialog and release GUI resources."""
        self.plotter._on_screen = False  # ruff:ignore[private-member-access]
        self.plotter._cancel_gto_future()  # ruff:ignore[private-member-access]
        self.plotter._pv_plotter.close()  # ruff:ignore[private-member-access]
        self.destroy()
        if self.plotter._tk_root and self.plotter._no_prev_tk_root:  # ruff:ignore[private-member-access]
            self.plotter._tk_root.quit()  # ruff:ignore[private-member-access]
            self.plotter._tk_root.destroy()  # ruff:ignore[private-member-access]

    def _protocols(self) -> None:
        """Attach standard close shortcuts to the dialog window."""
        self.protocol('WM_DELETE_WINDOW', self._on_close)
        self.bind('<Command-q>', lambda _event: self._on_close())
        self.bind('<Command-w>', lambda _event: self._on_close())
        self.bind('<Control-q>', lambda _event: self._on_close())
        self.bind('<Control-w>', lambda _event: self._on_close())

    def _next_plot(self) -> None:
        """Advance to the next molecular orbital."""
        if self._loading:
            return
        max_index = len(self.plotter.tabulator.molecular_orbitals) - 1
        if max_index < 0:
            return
        current = self.current_mo_ind
        new_index = 0 if current < 0 else min(current + 1, max_index)
        self.plotter.plot_orbital(new_index)
        if self.current_mo_ind >= 0:
            self.orb_tv._highlight_orbital(self.current_mo_ind)  # ruff:ignore[private-member-access]

    def _prev_plot(self) -> None:
        """Return to the previous molecular orbital."""
        if self._loading:
            return
        if self.current_mo_ind <= 0:
            return
        new_index = self.current_mo_ind - 1
        self.plotter.plot_orbital(new_index)
        self.orb_tv._highlight_orbital(self.current_mo_ind)  # ruff:ignore[private-member-access]

    def _set_loading_state(self, loading: bool, message: str = 'Tabulating orbitals...') -> None:
        """Toggle the inline loading label shown under the navigation buttons."""
        if loading:
            self._loading_label.config(text=message)
        if loading == self._loading:
            return
        self._loading = loading
        if loading:
            self._loading_label.pack(before=self.orb_tv, fill=tk.X, padx=10, pady=(0, 5))
            self.prev_button.config(state=tk.DISABLED)
            self.next_button.config(state=tk.DISABLED)
            self.orb_tv.configure(selectmode='none')
        else:
            self._loading_label.pack_forget()
            self.orb_tv.configure(selectmode='browse')
            self._update_nav_button_states()

    def _on_gtos_ready(self) -> None:
        """Handle the plotter callback when GTOs become available."""
        self._set_loading_state(False)

    def _update_nav_button_states(self) -> None:
        """Synchronize navigation button state with the current orbital index."""
        if self._loading:
            self.prev_button.config(state=tk.DISABLED)
            self.next_button.config(state=tk.DISABLED)
            return
        total = len(self.plotter.tabulator.molecular_orbitals)
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

    def _plot_orbital(self, orb_ind: int) -> None:
        """Render the selected orbital isosurface in the PyVista plotter.

        Parameters
        ----------
        orb_ind : int
            Index of the orbital to display; ``-1`` clears the current mesh.
        """
        if self._loading:
            return
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

        super().__init__(selection_screen, columns=columns, show='headings', height=20, selectmode='browse')

        for col, w in zip(columns, widths, strict=False):
            self.heading(col, text=col)
            self.column(col, width=w)

        self._selection_screen = selection_screen

        self.current_mo_ind = -1  # Start with no orbital shown

        # Configure tag
        self.tag_configure('highlight', background='lightblue')

        self.bind('<<TreeviewSelect>>', self._on_select)

    def _highlight_orbital(self, orb_ind: int) -> None:
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

    def _erase(self) -> None:
        """Remove all orbital entries from the tree view."""
        for item in self.get_children():
            self.delete(item)

    def _populate_tree(self, mos: list[MolecularOrbital]) -> None:
        """Populate the tree view with molecular orbital metadata.

        Parameters
        ----------
        mos : list[MolecularOrbital]
            Orbitals sourced from the parser.
        """
        self._erase()

        # Counts the number of MOs with a given symmetry
        mo_syms = list({mo.sym for mo in mos})
        mo_sym_count: dict[str, int] = dict.fromkeys(mo_syms, 0)
        for ind, mo in enumerate(mos):
            mo_sym_count[mo.sym] += 1
            self.insert('', 'end', iid=ind, values=(ind + 1, f'{mo.sym}.{mo_sym_count[mo.sym]}', mo.occ, mo.energy))

    def _on_select(self, _event: tk.Event) -> None:
        """Handle user selection events raised by the tree view.

        Parameters
        ----------
        _event : tk.Event
            Tkinter event object (unused).
        """
        if self._selection_screen._loading:  # ruff:ignore[private-member-access]
            return
        selected_item = self.selection()
        self.selection_remove(selected_item)
        if selected_item:
            orb_ind = int(selected_item[0])
            self._highlight_orbital(orb_ind)
            self._selection_screen.current_mo_ind = orb_ind
            self._selection_screen._plot_orbital(orb_ind)  # ruff:ignore[private-member-access]
            self._selection_screen._update_nav_button_states()  # ruff:ignore[private-member-access]
