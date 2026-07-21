# moldenViz API surface review

## Scope

This review covers the complete user-visible surface of the package at commit
`83bf131`: root-package exports, importable modules, project-defined classes and
types, class methods, documented Python behavior, the `moldenViz` command, and
the bundled examples. It is a recommendation report only; no implementation API
has been changed.

For this review, a name is considered:

- **supported public API** when it is exported from `moldenViz`, explicitly
  documented for users, or deliberately accepted as a compatibility contract;
- **provisionally public** when Python naming makes it public and the generated
  docs publish it, even if it appears to be an implementation detail;
- **private** when its name or containing module starts with `_` and it is not
  deliberately re-exported;
- **accidentally exposed** when it is reachable from a public module only because
  that module imports it without defining `__all__`.

Python cannot prevent a user from importing an internal name. The goal is a
clear supported contract, not runtime access control.

## Executive conclusion

The root API is currently sensible but too small:

```python
from moldenViz import Parser, Plotter, Tabulator
```

Those three classes should remain public. The main correction is to expose the
domain types that users receive from them while hiding parsing, rendering, and
GUI implementation steps.

The recommended root API is:

```python
from moldenViz import (
    Atom,
    AtomType,
    GaussianPrimitive,
    GridType,
    MolecularOrbital,
    Parser,
    Plotter,
    Shell,
    Tabulator,
    __version__,
)
```

In particular:

- Promote the parser result types `_Atom`, `_MolecularOrbital`, `_Shell`, and
  `_GTO` to supported public models. `Parser.atoms`, `Parser.mos`, and
  `Parser.shells` already expose instances of these types, so calling the types
  private does not provide real encapsulation.
- Promote `AtomType`, as anticipated, but move it out of `_config_module.py`.
  Importing that module currently loads plotting/configuration dependencies and
  creates `~/.config/moldenViz`; a lightweight public model import should not
  have that side effect.
- Make `Parser.divide_molden_lines()` and the other constructor-only parsing
  stages private.
- Keep the useful computational `Tabulator` API public, and add public read-only
  grid metadata instead of requiring users and `Plotter` to read `_grid_type`
  and `_grid_dimensions`.
- Keep a small, coherent set of programmatic `Plotter` controls public. Make
  rendering construction, colormap construction, grid rebuilding, and all UI
  callback methods private.
- Add explicit `__all__` declarations to public modules. Without them, imported
  names such as `Path`, `Any`, `NDArray`, `Enum`, `ThreadPoolExecutor`, NumPy,
  SciPy helpers, loggers, and plotting classes appear in star imports and remain
  reachable as public-looking module attributes.
- Replace blanket Sphinx `:members:` publication with explicit members or
  `:exclude-members:`. Documentation should define the supported contract rather
  than infer it from the absence of a leading underscore.

## Current package-level contract

`src/moldenViz/__init__.py` defines `__all__` as `Parser`, `Plotter`, and
`Tabulator`. `Plotter` is loaded lazily, which avoids importing the heavy GUI
stack for parser/tabulator-only workflows. That design should stay.

| Name | Current status | Recommendation |
| --- | --- | --- |
| `Parser` | Explicit root export and documented | Keep public |
| `Tabulator` | Explicit root export and documented | Keep public |
| `Plotter` | Explicit lazy root export and documented | Keep public and lazy |
| `GridType` | Public in `moldenViz.tabulator`, absent at root | Promote to root |
| Parser result models | Private names returned by public attributes | Promote to root and a public models module |
| `AtomType` | Public-looking class inside a private module | Promote to root and a public models module |
| `__version__` | Available only through `moldenViz.__about__` | Export at root |
| Other configuration models | In a private module and not consumed by a public API | Keep private for now |
| Rendered `Atom`, `Bond`, `Molecule` | In private `_plotting_objects` module | Keep private |

The public submodules should be limited to `parser`, `tabulator`, `plotter`,
`models` (new), and `examples`. `cli` remains an entry-point implementation.
Rename `plotter_ui.py` to `_plotter_ui.py`, and preferably
`examples/get_example_files.py` to `examples/_get_example_files.py`, in the next
major release or with a compatibility shim.

## Project-defined symbol inventory

This section accounts for every project-defined type, class, function, and type
alias in `src/moldenViz`.

### `moldenViz.parser`

| Symbol | Current status | Recommended status | Reason |
| --- | --- | --- | --- |
| `_Atom` | Private type returned by `Parser.atoms` | Public `Atom` | Users directly inspect its label, atomic number, position, and shells |
| `_MolecularOrbital` | Private type returned by `Parser.mos` | Public `MolecularOrbital` | It is the supported representation of MO metadata |
| `_GTO` | Private type reachable through shells | Public `GaussianPrimitive` | Reachable public data; the proposed name states that it is one primitive |
| `_Shell` | Private type returned by `Parser.shells` | Public `Shell` | Reachable public data and part of the parsed model |
| `Parser` | Public | Public | Primary parsing entry point |

`Parser` member review:

| Member | Recommendation | Notes |
| --- | --- | --- |
| `__init__(source, only_molecule=False)` | Keep public | This is the intended operation; tighten the annotation instead of accepting `Any` |
| `atoms`, `shells`, `mos`, `mo_coeffs` | Keep public | These are the useful parsed results; give them public model annotations |
| `check_molden_format()` | Rename `_check_molden_format()` | Called during construction; users should receive its exception, not orchestrate it |
| `divide_molden_lines()` | Rename `_divide_molden_lines()` | Section-index discovery is a parsing implementation detail |
| `get_atoms()` | Rename `_parse_atoms()` | It depends on private section-index state and mutates the construction sequence |
| `get_shells()` | Rename `_parse_shells()` | It depends on atoms and section-index state |
| `get_mos(sort=True)` | Rename `_parse_mos()` | Reparsing is not a clean public operation; expose sorting as a constructor option if needed |
| `_gto_order()` | Keep private | Basis-order conversion helper |
| `molden_lines` | Rename `_molden_lines` | Normalized source storage is implementation state |
| `ANGSTROM_TO_BOHR` | Make private, or expose as `BOHR_PER_ANGSTROM` at module level | The current class placement is unrelated to parser identity |

The public model fields should be documented and stable. Computational caches
such as shell normalization arrays should become underscored fields so promoting
`Shell` does not accidentally freeze every implementation detail.

`only_molecule=True` currently returns from `Parser.__init__` before `shells`,
`mos`, and `mo_coeffs` exist. A public class should have a stable object shape.
Initialize these attributes to empty results, or annotate and set them to `None`,
and document the chosen behavior.

### `moldenViz.tabulator`

| Symbol | Current status | Recommended status | Reason |
| --- | --- | --- | --- |
| `array_like_type` | Public-looking lowercase type alias | Private `_MOIndices`, or public `MOIndices` only if users need it | It exists only to shorten one annotation and should not be a runtime API |
| `_grid_creation_with_only_molecule_error()` | Private | Private | Internal error factory |
| `GridType` | Public in submodule | Public at submodule and root | Users need a supported way to inspect grid kind |
| `Tabulator` | Public at submodule and root | Public | Primary numerical API |

Keep these `Tabulator` members public:

- `grid` and `gtos` as read-only data properties;
- `cartesian_grid()` and `spherical_grid()`;
- `tabulate_gtos()` and `tabulate_mos()`;
- `export()`, `export_vtk()`, and `export_cube()`.

The format-specific exporters were deliberately documented and listed as a
public v1.11 feature, so making them private now would be an API regression.

Keep the existing `_axis_spacing()`, `_spherical_to_cartesian()`,
`_cartesian_to_spherical()`, `_set_grid()`, `_tabulate_atom()`,
`_tabulate_xlms()`, and `_check_bounds()` private.

Change the following state exposure:

| Current member | Recommendation |
| --- | --- |
| Writable/deletable `grid` property | Prefer a read-only property; add a deliberate `set_grid(...)` only if arbitrary unstructured grids are supported |
| Deletable `gtos` property | Prefer a read-only property; add `clear_gtos()` only if cache invalidation is a supported operation |
| Public `original_axes` attribute | Replace with read-only `grid_axes` or make it private |
| Private `_grid_type` | Add read-only public `grid_type: GridType` |
| Private `_grid_dimensions` | Add read-only public `grid_dimensions: tuple[int, int, int]` |
| Private `_parser` | Keep private; add narrow `atoms`, `orbitals`, or `parser` properties only when a real user workflow requires them |

`Plotter` currently accesses `Tabulator._grid_type`, `_axis_spacing`, and a
nonexistent `gto_data` name while validating a supplied tabulator. A public
grid-metadata contract would remove this cross-class reliance on internals and
make third-party use less ambiguous.

### `moldenViz.plotter`

| Symbol | Current status | Recommended status |
| --- | --- | --- |
| `_describe_source()` | Private | Keep private |
| `Plotter` | Public at submodule and root | Keep public and lazy at root |

Keep this coherent programmatic control surface public:

- `wait_for_gtos()`;
- `plot_orbital()`;
- `toggle_molecule()`, `toggle_atoms()`, and `toggle_bonds()`;
- `is_molecule_visible()`, `are_atoms_visible()`, and
  `are_bonds_visible()`.

Make these currently public-looking methods private:

| Current method | Proposed name | Reason |
| --- | --- | --- |
| `custom_cmap_from_colors()` | `_custom_cmap_from_colors()` | Rendering helper, not an operation on a plotter |
| `load_molecule()` | `_load_molecule()` | Accepts the internal `Config` type and rebuilds internal actors |
| `update_mesh()` | `_update_mesh()` | Grid-settings callback and renderer rebuild step |

The existing `_schedule_gto_tabulation()`, `_on_gtos_ready()`,
`_apply_gtos_ready()`, `_ensure_gtos_ready()`, `_cancel_gto_future()`,
`_connect_pv_plotter_close_signal()`, `_clear_all()`, and `_create_mo_mesh()` are
correctly private.

The two window-size class constants and GUI object attributes (`tk_root`,
`pv_plotter`, selection windows, meshes, actors, and callback state) should be
documented as implementation state. `tabulator` is the one useful public
attribute because callers may supply and reuse it.

### `moldenViz.plotter_ui`

The module name is public-looking, but all four classes are correctly named as
private: `_ConfigProxy`, `_PlotterUI`, `_OrbitalSelectionScreen`, and
`_OrbitalsTreeview`. `_plotter_config()` is also correctly private. Rename the
module itself to `_plotter_ui.py` so its filename agrees with its role.

All non-underscored methods on these private classes are implementation details,
even though they are callable through `Plotter` inheritance. This includes:

- dialog and settings handlers: `export_orbitals_dialog`,
  `export_image_dialog`, `grid_settings_screen`, `mo_settings_screen`,
  `molecule_settings_screen`, `color_settings_screen`, and `save_settings`;
- event/apply handlers: `on_opacity_change`, `on_molecule_opacity_change`,
  `on_mo_color_scheme_change`, `on_bond_color_type_change`,
  `apply_mo_contour`, `apply_background_color`, `apply_grid_settings`,
  `apply_molecule_settings`, `apply_color_settings`,
  `apply_mo_color_settings`, `apply_custom_mo_color_settings`, and
  `apply_bond_color_settings`;
- widget setup/reset helpers: `update_settings_button_states`,
  `place_grid_params_frame`, `sph_grid_params_frame_widgets`,
  `cart_grid_params_frame_widgets`, `sph_grid_params_frame_setup`,
  `cart_grid_params_frame_setup`, `reset_grid_settings`, `reset_mo_settings`,
  `reset_molecule_settings`, and `reset_color_settings`;
- selection-window/tree handlers: `on_close`, `protocols`, `next_plot`,
  `prev_plot`, `set_loading_state`, `on_gtos_ready`,
  `update_nav_button_states`, `plot_orbital`, `highlight_orbital`, `erase`,
  `populate_tree`, and `on_select`.

These methods should either receive leading underscores or be excluded
explicitly from public documentation. There is no need to promise each Tk/Qt
callback as an independent user API.

### `moldenViz._config_module`

The leading underscore makes this entire module private today. Its
project-defined classes are `AtomType`, `SphericalGridConfig`,
`CartesianGridConfig`, `GridConfig`, `MOConfig`, `AtomDisplayConfig`,
`BondConfig`, `MoleculeConfig`, `MainConfig`, and `Config`.

Promote only `AtomType` now. Move it to a lightweight public `models.py` module
and re-export it at the package root. This directly supports atom inspection and
customization without declaring the whole current configuration loader stable.

Keep the remaining Pydantic configuration models and `Config` private until a
public constructor or method actually accepts them. Promoting schema classes
without a supported way to pass them into `Plotter` would create a misleading
API. If programmatic configuration is added later, promote the full cohesive
model tree together rather than selected nested classes.

All current `Config` methods are loader/UI implementation details:
`dict_to_namedspace()`, `merge_configs()`, `recursive_merge()`,
`load_atom_types()`, `load_default_config()`, `load_custom_config()`, and
`save_current_config()`. Their names should become private if `Config` remains.
The Pydantic validator methods are framework hooks and are private API even
though their names do not start with `_`.

The module creates the user config directory during import. Move directory
creation to the operation that writes a configuration. Public type imports
should be read-only and side-effect free.

### `moldenViz._plotting_objects`

The module and its `Atom`, `Bond`, and `Molecule` classes should remain private.
They represent PyVista meshes and inferred bonds, not parsed chemical data. This
is distinct from the recommended public parser `Atom` model. Their methods
`remove_extra_bonds()`, `trim_ends()`, `get_atoms()`, and `add_meshes()` are also
rendering implementation details; `_trim_atom_from_bond()` is already named
correctly.

Add `__all__ = ["Plotter"]` to `moldenViz.plotter` so importing `Molecule`,
`Config`, `GridType`, `Tabulator`, PyVista classes, or executor classes into that
module does not imply they are supported from that path.

### `moldenViz.cli`

| Symbol | Current status | Recommendation |
| --- | --- | --- |
| `_resolve_plotter()` | Private | Keep private |
| `ColorFormatter` | Public-looking class | Rename `_ColorFormatter` |
| `main()` | Console-script entry function | Keep as entry infrastructure; do not advertise as general Python API |

The public command name `moldenViz` and the current positional file/example,
molecule-only, version, and logging behaviors should remain stable.

There is one documentation/behavior mismatch: the code accepts
`--only_molecule`, while the CLI guide documents the conventional
`--only-molecule`. Add `--only-molecule` as the canonical spelling, retain
`--only_molecule` as a deprecated alias for at least one release, and keep `-m`.
The other flags (`--example`, `--version`, `--verbose`, `--debug`, and `--quiet`)
should stay unchanged.

### `moldenViz.examples`

The explicit `__all__` is good. Keep the documented `acrolein`, `benzene`, `co`,
`co2`, `furan`, `h2o`, `o2`, `prismane`, `pyridine`, and `all_examples` names
public.

`_read_file()` is correctly private. `molden_files_folder` is public-looking
only in the implementation module and should become `_molden_files_folder`.
The implementation module can become `_get_example_files.py`; users should
import from `moldenViz.examples`, as the documentation already demonstrates.

### Package metadata and magic hooks

`moldenViz.__init__.__getattr__()` is an internal language hook and should remain
undocumented. Extend it as needed to preserve lazy imports for heavy public
objects. `__version__` should be deliberately re-exported from the package root;
`__about__.py` itself does not need to be public.

## Recommended supported API after migration

### Root imports

```python
from moldenViz import (
    Atom,
    AtomType,
    GaussianPrimitive,
    GridType,
    MolecularOrbital,
    Parser,
    Plotter,
    Shell,
    Tabulator,
    __version__,
)
```

### `Parser`

```text
Parser(source, only_molecule=False)
  .atoms
  .shells
  .mos
  .mo_coeffs
```

Parsing phases and normalized source lines remain private.

### `Tabulator`

```text
Tabulator(source, only_molecule=False)
  .grid
  .gtos
  .grid_type
  .grid_dimensions
  .grid_axes
  .cartesian_grid(...)
  .spherical_grid(...)
  .tabulate_gtos()
  .tabulate_mos(...)
  .export(...)
  .export_vtk(...)
  .export_cube(...)
```

### `Plotter`

```text
Plotter(source, only_molecule=False, tabulator=None, tk_root=None)
  .tabulator
  .wait_for_gtos(...)
  .plot_orbital(...)
  .toggle_molecule()
  .toggle_atoms()
  .toggle_bonds()
  .is_molecule_visible()
  .are_atoms_visible()
  .are_bonds_visible()
```

GUI widgets, callbacks, mesh construction, actor management, and config loading
remain private.

## Compatibility and implementation order

Because the project identifies itself as production/stable and Sphinx currently
publishes all non-underscored methods on the main classes, treat removals as
deprecations even when they were not intended to be public.

### Phase 1: additive and non-breaking

1. Add explicit `__all__` declarations to `parser`, `tabulator`, `plotter`,
   `cli`, and the new public models module.
2. Add public model names and keep private aliases temporarily for internal and
   downstream compatibility.
3. Re-export the public models, `GridType`, and `__version__` from the root;
   preserve lazy loading for GUI-dependent names.
4. Add read-only `Tabulator.grid_type`, `grid_dimensions`, and `grid_axes`.
5. Add the canonical `--only-molecule` spelling while retaining the underscore
   spelling as an alias.
6. Replace blanket Sphinx member discovery with an explicit public-member list.
7. Add API-contract tests for root `__all__`, public import paths, private-name
   exclusion, lightweight imports, and CLI aliases.

### Phase 2: deprecate misleading names

1. Move parser construction to underscored methods. Leave warning-emitting
   wrappers for the currently documented public spellings for one release.
2. Make `custom_cmap_from_colors`, `load_molecule`, and `update_mesh` private,
   with deprecation wrappers if release history shows external use.
3. Rename the UI and example implementation modules with compatibility shims.
4. Move `AtomType` and other shared models out of the side-effecting config
   loader.
5. Stop internal code from accessing `Tabulator` private members.

### Phase 3: next major release

Remove deprecated wrappers and module shims. Publish a concise API stability
policy stating that root exports, documented members, and explicit module
`__all__` entries are supported; underscored names and modules are not.

## Tests that should guard the boundary

- `from moldenViz import ...` succeeds for every supported root name without
  importing PyVista/Qt unless `Plotter` is requested.
- `moldenViz.__all__` exactly matches the supported root contract.
- Each public module's `__all__` contains only project-owned supported names.
- `Parser.atoms`, `shells`, and `mos` contain the documented public model types.
- `Parser(..., only_molecule=True)` has a documented, stable set of attributes.
- `Tabulator.grid_type`, `grid_dimensions`, and `grid_axes` reflect both grid
  constructors and reset behavior.
- Sphinx builds list only supported methods.
- `moldenViz -m`, `--only-molecule`, and the temporary `--only_molecule` alias
  behave identically.
- Importing public lightweight models does not create configuration directories
  or import GUI libraries.

## Final recommendation

Keep `Parser`, `Tabulator`, `Plotter`, the computational tabulation/export
methods, the small programmatic plot controls, examples, and CLI behaviors as
the stable core. Promote the data users actually receive—especially `Atom` and
`AtomType`—and `GridType`. Privatize orchestration steps such as
`divide_molden_lines`, parser phase methods, mesh rebuilding, configuration
loading, and UI callbacks. Use explicit exports and explicit documentation so
future refactors do not accidentally expand the compatibility contract.
