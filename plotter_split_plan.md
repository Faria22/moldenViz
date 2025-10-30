# Plotter Module Split Plan

This document outlines the proposed extraction of the Tk/Qt user‑interface helpers
out of `src/moldenViz/plotter.py` into a dedicated companion module.

## Goals
- Reduce `plotter.py` size (1,900+ lines) and improve cohesion.
- Isolate the Tk/Qt UI classes so the remaining plotter core focuses on PyVista
  scene control, tabulation, and molecule rendering.
- Simplify headless testing by letting fixtures mock or skip the UI module.

## Proposed Module Layout

### `src/moldenViz/plotter.py`
- Keep `Plotter` and helpers directly related to molecule/orbital rendering.
- Provide a narrow interface the UI layer can call:
  - `plot_orbital`, `export_*`, grid setters, visibility toggles, etc.
- Host shared enums/constants or factories as needed (e.g., `PlotterMenus` factory).

### `src/moldenViz/plotter_ui.py`
- New module housing Tk widget classes and dialogs:
  - `_OrbitalSelectionScreen`
  - `_OrbitalsTreeview`
  - Menu/dialog builders (`grid_settings_screen`, `mo_settings_screen`, etc.).
- Accept a `Plotter` (or protocol) instance to access drawing/export operations.
- Keep Tk/Qt imports local to this module to avoid dragging UI deps into the core.

## Migration Steps
1. **Create `plotter_ui.py`:**
   - Move `_OrbitalSelectionScreen`, `_OrbitalsTreeview`, and related functions.
   - Replace direct attribute access with well-defined calls to the `Plotter` API.
2. **Update `plotter.py`:**
   - Import the new UI helpers.
   - Remove now duplicated Tk/Qt import statements.
   - Expose any additional methods the UI requires (consider a lightweight protocol).
3. **Adjust tests:**
   - Point UI-focused tests to the new module.
   - Add coverage for the delegated interfaces if missing.
4. **Manual validation:**
   - Run `hatch run all`.
   - Smoke-test the interactive plotter to ensure menus/dialogs still wire up.

## Notes & Risks
- Ensure the new module does not cause circular imports; avoid referencing
  `Plotter` at import time inside `plotter_ui.py` by using type-checking guards.
- Keep configuration objects (`Config`) and logging in one place to prevent divergence.
- Document the split in the developer docs (optional) if onboarding contributors rely on the old structure.
