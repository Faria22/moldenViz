# Code Health Findings

Notes gathered while reviewing the current implementation; use these as
starting points for future improvements.

## Coordinate Conversion Robustness
- File: `src/moldenViz/tabulator.py:75`
- Issue: `_cartesian_to_spherical` computes `theta = np.arccos(z / r)` without
  guarding the origin. When the grid includes `(0, 0, 0)` (default symmetric
  grids do) `r == 0`, producing `z / r` NaNs that propagate into the spherical
  harmonics tabulation.
- Potential fix: branch on `r == 0` and set `theta`/`phi` explicitly before
  applying `arccos`.
- Bonus: clip `z / r` to `[-1, 1]` to avoid NaNs caused by floating-point
  round-off in near-origin points.

## Plotter Event Loop Coupling
- File: `src/moldenViz/plotter.py:1`
- Issue: `Plotter.__init__` always runs `tk_root.mainloop()` when it created the
  root internally, which makes constructing a plotter block the caller’s
  thread. This prevents headless/testing usage where the caller wants to manage
  the event loop manually.
- Potential fix: move mainloop control out to the CLI/UI entry point and expose
  an explicit method (e.g., `run()`) so library consumers can opt-in to the GUI
  loop.

## Import-Time Side Effects
- File: `src/moldenViz/plotter.py:42`
- Issue: Instantiating `Config()` at module import forces configuration IO and
  GUI defaults to load in every context (including tests that just import
  `plotter`). This makes the module heavy to import and complicates dependency
  injection.
- Potential fix: lazily load the config inside `Plotter.__init__`, or accept a
  `Config` instance via constructor parameters so the module stays import-light.

## Layering Between Plotter and Tabulator
- Files: `src/moldenViz/plotter.py:247`, `src/moldenViz/plotter.py:1543`
- Issue: `Plotter` pokes into `Tabulator._parser` to grab atoms and orbital
  metadata, coupling two layers via a private attribute.
- Potential fix: expose dedicated accessors on `Tabulator` (e.g.,
  `tabulator.atoms`, `tabulator.mo_metadata(idx)`) so the plotting layer no
  longer depends on parser internals.

## Export Responsibilities in Tabulator
- File: `src/moldenViz/tabulator.py:494-579`
- Issue: `Tabulator` mixes grid tabulation logic with VTK and cube export code,
  which tightens coupling between computation and serialization.
- Potential fix: move exporters to helper modules or inject them as strategies
  so the core remains focused on generating data.

## Public API Definition for Tests
- File: `tests/_src_imports.py:29`
- Issue: Tests need to re-export private symbols (`_cartesian_to_spherical`,
  `_GTO`, etc.), signalling that the package lacks an official API surface for
  these helpers.
- Potential fix: promote the required utilities into a documented public module,
  allowing tests (and users) to import them without reaching into private names.
