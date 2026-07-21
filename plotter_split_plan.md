# Plotter Module Split Plan

Status: completed.

This document records the extraction of the Tk/Qt user-interface helpers
out of `src/moldenViz/plotter.py` into a dedicated companion module.

## Goals
- Reduce `plotter.py` size (1,900+ lines) and improve cohesion.
- Isolate the Tk/Qt UI classes so the remaining plotter core focuses on PyVista
  scene control, tabulation, and molecule rendering.
- Simplify headless testing by letting fixtures mock or skip the UI module.

## Module Layout

### `src/moldenViz/plotter.py`
- Keep `Plotter` and helpers directly related to molecule/orbital rendering.
- Provides the rendering interface the UI layer calls:
  - `plot_orbital`, `export_*`, grid setters, visibility toggles, etc.
- Host shared enums/constants or factories as needed (e.g., `PlotterMenus` factory).

### `src/moldenViz/plotter_ui.py`
- Houses the `_PlotterUI` mixin, Tk widget classes, and dialogs:
  - `_OrbitalSelectionScreen`
  - `_OrbitalsTreeview`
  - Menu/dialog builders (`grid_settings_screen`, `mo_settings_screen`, etc.).
- Uses a mixin so UI callbacks can access the active `Plotter` without circular runtime imports.
- Centralizes UI-specific Tk/Qt imports and configuration access in the companion module.

## Completed Migration
1. **Created `plotter_ui.py`:**
   - Moved `_OrbitalSelectionScreen`, `_OrbitalsTreeview`, and related functions.
   - Moved menus, dialogs, settings controls, and export UI into `_PlotterUI`.
2. **Updated `plotter.py`:**
   - Imports and inherits the new UI mixin.
   - Retains rendering, tabulation, lifecycle, and molecule visibility behavior.
3. **Adjusted tests:**
   - Pointed UI-focused patches to the new module.
   - Added a regression test that fixes the module boundary.
4. **Validated the split:**
   - `hatch run all` covers lint, types, and the mocked menu/dialog wiring.

## Notes & Risks
- `Plotter` is referenced only for type checking; runtime configuration access is lazy to avoid circular imports.
- Private UI classes remain re-exported from `plotter.py` for compatibility.
