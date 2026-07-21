API-001 | Change: Keep `Parser` as a supported package-root export. | Approve: [x] | Disapprove: [ ] | Comment:
API-002 | Change: Keep `Tabulator` as a supported package-root export. | Approve: [x] | Disapprove: [ ] | Comment:
API-003 | Change: Keep `Plotter` as a supported, lazily loaded package-root export. | Approve: [x] | Disapprove: [ ] | Comment:
API-004 | Change: Export `__version__` from the package root. | Approve: [x] | Disapprove: [ ] | Comment:
API-005 | Change: Promote parser `_Atom` to the public name `Atom` and export it from the package root. | Approve: [x] | Disapprove: [ ] | Comment:
API-006 | Change: Promote parser `_MolecularOrbital` to the public name `MolecularOrbital` and export it from the package root. | Approve: [x] | Disapprove: [ ] | Comment:
API-007 | Change: Promote parser `_Shell` to the public name `Shell` and export it from the package root. | Approve: [x] | Disapprove: [ ] | Comment:
API-008 | Change: Promote parser `_GTO` to the public name `GaussianPrimitive` and export it from the package root. | Approve: [x] | Disapprove: [ ] | Comment:
API-009 | Change: Promote `AtomType` to a supported public type and export it from the package root. | Approve: [x] | Disapprove: [ ] | Comment:
API-010 | Change: Export `GridType` from the package root. | Approve: [x] | Disapprove: [ ] | Comment:
API-011 | Change: Add a lightweight public `moldenViz.models` module for `Atom`, `AtomType`, `MolecularOrbital`, `Shell`, and `GaussianPrimitive`. | Approve: [x] | Disapprove: [ ] | Comment:
API-012 | Change: Preserve lazy imports so importing lightweight public models does not import PyVista, Qt, or Tk GUI code. | Approve: [x] | Disapprove: [ ] | Comment:
API-013 | Change: Define the supported API as package-root exports, explicitly documented members, and names in module `__all__` declarations. | Approve: [x] | Disapprove: [ ] | Comment:
API-014 | Change: Add an explicit `__all__` declaration to `moldenViz.parser`. | Approve: [x] | Disapprove: [ ] | Comment:
API-015 | Change: Add an explicit `__all__` declaration to `moldenViz.tabulator`. | Approve: [x] | Disapprove: [ ] | Comment:
API-016 | Change: Add `__all__ = ["Plotter"]` to `moldenViz.plotter`. | Approve: [x] | Disapprove: [ ] | Comment: Approval mark was corrected from the proposal text.
API-017 | Change: Add an explicit `__all__` declaration to `moldenViz.cli`. | Approve: [x] | Disapprove: [ ] | Comment:
API-018 | Change: Stop treating imported dependencies such as `Path`, `Any`, `NDArray`, NumPy, SciPy helpers, executors, plotting classes, and loggers as public module exports. | Approve: [x] | Disapprove: [ ] | Comment:
API-019 | Change: Replace blanket Sphinx `:members:` API publication with explicit supported-member lists. | Approve: [x] | Disapprove: [ ] | Comment:
API-020 | Change: Add API-contract tests that assert the exact package-root `__all__`. | Approve: [x] | Disapprove: [ ] | Comment:
API-021 | Change: Add API-contract tests for supported public import paths and exclusion of private names. | Approve: [x] | Disapprove: [ ] | Comment:
API-022 | Change: Add API-contract tests proving lightweight imports do not load GUI libraries or create configuration directories. | Approve: [x] | Disapprove: [ ] | Comment:
API-023 | Change: Keep `Parser.__init__(source, only_molecule=False)` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-024 | Change: Tighten the `Parser` source annotation so it no longer includes `Any`. | Approve: [x] | Disapprove: [ ] | Comment:
API-025 | Change: Keep `Parser.atoms` public and annotate it with the public `Atom` type. | Approve: [x] | Disapprove: [ ] | Comment:
API-026 | Change: Keep `Parser.shells` public and annotate it with the public `Shell` type. | Approve: [x] | Disapprove: [ ] | Comment:
API-027 | Change: Keep `Parser.mos` public and annotate it with the public `MolecularOrbital` type. | Approve: [x] | Disapprove: [ ] | Comment:
API-028 | Change: Keep `Parser.mo_coeffs` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-029 | Change: Rename `Parser.check_molden_format()` to `_check_molden_format()`. | Approve: [x] | Disapprove: [ ] | Comment:
API-030 | Change: Rename `Parser.divide_molden_lines()` to `_divide_molden_lines()`. | Approve: [x] | Disapprove: [ ] | Comment:
API-031 | Change: Rename `Parser.get_atoms()` to `_parse_atoms()`. | Approve: [x] | Disapprove: [ ] | Comment:
API-032 | Change: Rename `Parser.get_shells()` to `_parse_shells()`. | Approve: [x] | Disapprove: [ ] | Comment:
API-033 | Change: Rename `Parser.get_mos()` to `_parse_mos()`. | Approve: [x] | Disapprove: [ ] | Comment:
API-034 | Change: If unsorted molecular orbitals are needed publicly, expose that choice as a `Parser` constructor option instead of requiring a second call to `get_mos(sort=False)`. | Approve: [ ] | Disapprove: [ ] | Comment: I'm not sure on this, will need to expand on it later. Feel free to create an issue or discussion board on github on this.
API-035 | Change: Keep `Parser._gto_order()` private. | Approve: [x] | Disapprove: [ ] | Comment:
API-036 | Change: Rename public-looking `Parser.molden_lines` storage to `_molden_lines`. | Approve: [x] | Disapprove: [ ] | Comment:
API-037 | Change: Remove `ANGSTROM_TO_BOHR` from the `Parser` class and keep the conversion constant private. | Approve: [ ] | Disapprove: [x] | Comment: Should stay public (look at API-038)
API-038 | Change: If the unit-conversion constant is intended for users, expose it at module level as `BOHR_PER_ANGSTROM` instead of as a `Parser` class attribute. | Approve: [x] | Disapprove: [ ] | Comment:
API-039 | Change: Make `Parser(..., only_molecule=True)` initialize `shells`, `mos`, and `mo_coeffs` to documented empty or optional values rather than leaving the attributes absent. | Approve: [x] | Disapprove: [ ] | Comment:
API-040 | Change: Document and stabilize the public fields of `Atom`, `MolecularOrbital`, `Shell`, and `GaussianPrimitive`. | Approve: [x] | Disapprove: [ ] | Comment:
API-041 | Change: Rename shell computation caches such as normalization arrays with leading underscores before making `Shell` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-042 | Change: Replace public-looking `array_like_type` with the private alias `_MOIndices`. | Approve: [x] | Disapprove: [ ] | Comment:
API-043 | Change: Keep `_grid_creation_with_only_molecule_error()` private. | Approve: [x] | Disapprove: [ ] | Comment:
API-044 | Change: Keep `Tabulator.grid` as a public read-only property. | Approve: [x] | Disapprove: [ ] | Comment:
API-045 | Change: Remove the public `Tabulator.grid` setter and deleter. | Approve: [x] | Disapprove: [ ] | Comment: Clarified that changes should only be made through `set_grid()`.
API-046 | Change: Add a deliberate public `Tabulator.set_grid(...)` method only if arbitrary user-supplied grids are meant to be supported. | Approve: [x] | Disapprove: [ ] | Comment:
API-047 | Change: Keep `Tabulator.gtos` as a public read-only property. | Approve: [x] | Disapprove: [ ] | Comment:
API-048 | Change: Remove the public `Tabulator.gtos` deleter. | Approve: [x] | Disapprove: [ ] | Comment:
API-049 | Change: Add `Tabulator.clear_gtos()` only if manual cache invalidation is meant to be supported. | Approve: [ ] | Disapprove: [ ] | Comment: sno sure on this, will have to expand later.
API-050 | Change: Add a public read-only `Tabulator.grid_type: GridType` property. | Approve: [x] | Disapprove: [ ] | Comment:
API-051 | Change: Add a public read-only `Tabulator.grid_dimensions: tuple[int, int, int]` property. | Approve: [x] | Disapprove: [ ] | Comment: Approval mark was corrected from the proposal text.
API-052 | Change: Replace public `Tabulator.original_axes` with a read-only `grid_axes` property. | Approve: [x] | Disapprove: [ ] | Comment:
API-053 | Change: Keep `Tabulator._parser` private. | Approve: [x] | Disapprove: [ ] | Comment:
API-054 | Change: Stop `Plotter` and other internal code from directly accessing `Tabulator._grid_type`, `_grid_dimensions`, and `_axis_spacing()`. | Approve: [x] | Disapprove: [ ] | Comment:
API-055 | Change: Keep `Tabulator.cartesian_grid()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-056 | Change: Keep `Tabulator.spherical_grid()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-057 | Change: Keep `Tabulator.tabulate_gtos()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-058 | Change: Keep `Tabulator.tabulate_mos()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-059 | Change: Keep `Tabulator.export()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-060 | Change: Keep `Tabulator.export_vtk()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-061 | Change: Keep `Tabulator.export_cube()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-062 | Change: Keep `_axis_spacing()`, `_spherical_to_cartesian()`, `_cartesian_to_spherical()`, `_set_grid()`, `_tabulate_atom()`, `_tabulate_xlms()`, and `_check_bounds()` private. | Approve: [x] | Disapprove: [ ] | Comment:
API-063 | Change: Keep `Plotter.wait_for_gtos()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-064 | Change: Keep `Plotter.plot_orbital()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-065 | Change: Keep `Plotter.toggle_molecule()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-066 | Change: Keep `Plotter.toggle_atoms()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-067 | Change: Keep `Plotter.toggle_bonds()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-068 | Change: Keep `Plotter.is_molecule_visible()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-069 | Change: Keep `Plotter.are_atoms_visible()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-070 | Change: Keep `Plotter.are_bonds_visible()` public. | Approve: [x] | Disapprove: [ ] | Comment:
API-071 | Change: Keep `Plotter.tabulator` as a public attribute. | Approve: [x] | Disapprove: [ ] | Comment:
API-072 | Change: Rename `Plotter.custom_cmap_from_colors()` to `_custom_cmap_from_colors()`. | Approve: [x] | Disapprove: [ ] | Comment:
API-073 | Change: Rename `Plotter.load_molecule()` to `_load_molecule()`. | Approve: [x] | Disapprove: [ ] | Comment:
API-074 | Change: Rename `Plotter.update_mesh()` to `_update_mesh()`. | Approve: [x] | Disapprove: [ ] | Comment:
API-075 | Change: Keep Plotter background-tabulation methods, close-signal wiring, clearing, and mesh construction private. | Approve: [x] | Disapprove: [ ] | Comment:
API-076 | Change: Make Plotter window-size constants private. | Approve: [x] | Disapprove: [ ] | Comment:
API-077 | Change: Treat Plotter Tk/Qt widgets, PyVista objects, meshes, actors, selection windows, and callback state as private implementation attributes. | Approve: [x] | Disapprove: [ ] | Comment:
API-078 | Change: Rename public-looking `plotter_ui.py` to `_plotter_ui.py`. | Approve: [x] | Disapprove: [ ] | Comment:
API-079 | Change: Make all inherited export-dialog and settings-dialog methods from `_PlotterUI` private. | Approve: [x] | Disapprove: [ ] | Comment:
API-080 | Change: Make all inherited UI event and apply handlers from `_PlotterUI` private. | Approve: [x] | Disapprove: [ ] | Comment:
API-081 | Change: Make all inherited widget setup and reset helpers from `_PlotterUI` private. | Approve: [x] | Disapprove: [ ] | Comment:
API-082 | Change: Keep `_OrbitalSelectionScreen` and `_OrbitalsTreeview`, including their navigation and event methods, private. | Approve: [x] | Disapprove: [ ] | Comment:
API-083 | Change: Move `AtomType` out of `_config_module.py` so importing the public type does not load the configuration service. | Approve: [x] | Disapprove: [ ] | Comment:
API-084 | Change: Keep `SphericalGridConfig`, `CartesianGridConfig`, `GridConfig`, `MOConfig`, `AtomDisplayConfig`, `BondConfig`, `MoleculeConfig`, and `MainConfig` private until a public API accepts them. | Approve: [x] | Disapprove: [ ] | Comment:
API-085 | Change: Keep the `Config` loader class private. | Approve: [x] | Disapprove: [ ] | Comment:
API-086 | Change: Rename `Config.dict_to_namedspace()`, `merge_configs()`, `recursive_merge()`, `load_atom_types()`, `load_default_config()`, `load_custom_config()`, and `save_current_config()` as private implementation methods. | Approve: [x] | Disapprove: [ ] | Comment:
API-087 | Change: Move creation of `~/.config/moldenViz` from module import time to the configuration-save operation. | Approve: [x] | Disapprove: [ ] | Comment:
API-088 | Change: Keep rendered `_plotting_objects.Atom`, `Bond`, and `Molecule` private and distinct from the public parsed `Atom` model. | Approve: [x] | Disapprove: [ ] | Comment:
API-089 | Change: Keep rendered-object methods such as `remove_extra_bonds()`, `trim_ends()`, `get_atoms()`, and `add_meshes()` private. | Approve: [x] | Disapprove: [ ] | Comment:
API-090 | Change: Keep the public CLI command name `moldenViz` unchanged. | Approve: [x] | Disapprove: [ ] | Comment:
API-091 | Change: Rename CLI `ColorFormatter` to `_ColorFormatter`. | Approve: [x] | Disapprove: [ ] | Comment:
API-092 | Change: Keep `cli.main()` as console-entry infrastructure but do not document it as a general Python API. | Approve: [x] | Disapprove: [ ] | Comment:
API-093 | Change: Add `--only-molecule` as the canonical long CLI option. | Approve: [x] | Disapprove: [ ] | Comment:
API-094 | Change: Retain `--only_molecule` as a deprecated compatibility alias for at least one release. | Approve: [ ] | Disapprove: [x] | Comment: it should not be depricated
API-095 | Change: Keep the CLI flags `-m`, `--example`, `--version`, `--verbose`, `--debug`, and `--quiet` behavior unchanged. | Approve: [x] | Disapprove: [ ] | Comment:
API-096 | Change: Keep `acrolein`, `benzene`, `co`, `co2`, `furan`, `h2o`, `o2`, `prismane`, and `pyridine` public, but remove public `all_examples`. | Approve: [x] | Disapprove: [ ] | Comment: Clarified that only the individual examples should remain public.
API-097 | Change: Rename `molden_files_folder` to `_molden_files_folder`. | Approve: [x] | Disapprove: [ ] | Comment:
API-098 | Change: Rename `examples/get_example_files.py` to `examples/_get_example_files.py`. | Approve: [x] | Disapprove: [ ] | Comment:
API-099 | Change: Keep `examples._read_file()` private. | Approve: [x] | Disapprove: [ ] | Comment:
API-100 | Change: Introduce public names and properties additively before removing or renaming existing names. | Approve: [x] | Disapprove: [ ] | Comment:
API-101 | Change: Keep warning-emitting compatibility wrappers for currently documented parser methods for one release. | Approve: [ ] | Disapprove: [x] | Comment:
API-102 | Change: Keep compatibility wrappers for renamed public-looking Plotter methods for one release if release history shows external use. | Approve: [ ] | Disapprove: [x] | Comment:
API-103 | Change: Keep compatibility import shims for renamed public-looking modules until the next major release. | Approve: [ ] | Disapprove: [x] | Comment: The next release will be major, so breaking changes can be made now.
API-104 | Change: Remove deprecated wrappers and module shims in the next major release. | Approve: [x] | Disapprove: [ ] | Comment:
API-105 | Change: Publish an API stability policy explaining that underscored names and modules are unsupported. | Approve: [x] | Disapprove: [ ] | Comment:
