# Screenshot Placeholders Reference

This document lists all screenshot placeholder locations added to the documentation. These placeholders indicate where actual screenshots should be inserted in the future to demonstrate moldenViz features.

## Getting Started (`getting-started.rst`)

### CLI Preview
- **File**: `_static/placeholder-cli.png` ✓ (existing)
  - **Location**: Quick CLI Preview section
  - **Description**: CLI window showing benzene with orbital controls
  - **Shows**: Full CLI interface with molecule visualization

### Python Preview
- **File**: `_static/placeholder-python.png` ✓ (existing)
  - **Location**: Quick Python Preview section
  - **Description**: Python plotter rendering an isosurface next to molecule controls
  - **Shows**: Interactive Python API plotter window

## CLI Guide (`cli-guide.rst`)

### Basic Rendering
- **File**: `_static/placeholder-full-viewer.png`
  - **Location**: Basic Rendering section
  - **Description**: Full moldenViz viewer showing molecule with orbital isosurface
  - **Shows**: Complete viewer interface with molecule structure and orbital visualization
  
- **File**: `_static/placeholder-molecule-only.png`
  - **Location**: Basic Rendering section (molecule-only mode)
  - **Description**: Molecule-only view without orbital surfaces
  - **Shows**: Viewer with just molecular structure (no orbitals)

### Bundled Examples
- **File**: `_static/placeholder-example-co.png`
  - **Location**: Bundled Examples section
  - **Description**: Carbon monoxide molecule with orbital visualization
  - **Shows**: CO molecule example with orbital isosurface

- **File**: `_static/placeholder-examples-grid.png`
  - **Location**: Bundled Examples section
  - **Description**: Grid showing various example molecules available in moldenViz
  - **Shows**: Collage of all built-in example molecules (co, o2, co2, h2o, benzene, prismane, pyridine, furan, acrolein)

## Python API Guide (`python-api.rst`)

### Plotting Molecules
- **File**: `_static/placeholder-plotter-comparison.png`
  - **Location**: Plotting Molecules section
  - **Description**: Side-by-side comparison of full orbital visualization vs molecule-only view
  - **Shows**: Split view demonstrating `only_molecule=True` vs `only_molecule=False`

- **File**: `_static/placeholder-controls-panel.png`
  - **Location**: Interactive Controls subsection
  - **Description**: Detailed view of the orbital controls panel showing sliders and buttons
  - **Shows**: Close-up of Tkinter control panel with all widgets visible

### Exporting from GUI
- **File**: `_static/placeholder-export-menu.png`
  - **Location**: Exporting from the GUI subsection
  - **Description**: Export menu in the plotter window showing Data and Image options
  - **Shows**: PyVista menu bar with Export menu expanded

- **File**: `_static/placeholder-data-export-dialog.png`
  - **Location**: Exporting Data subsection
  - **Description**: Data export dialog showing format and scope selection options
  - **Shows**: Export dialog with VTK/Cube format selection and current/all orbital options

- **File**: `_static/placeholder-image-export-dialog.png`
  - **Location**: Exporting Images subsection
  - **Description**: Image export dialog showing format options and transparent background checkbox
  - **Shows**: Image export dialog with PNG/JPEG/SVG/PDF options and transparency setting

- **File**: `_static/placeholder-exported-comparison.png`
  - **Location**: After image export instructions
  - **Description**: Comparison of different export formats (PNG, SVG, PDF) showing the same molecular visualization
  - **Shows**: Three side-by-side images demonstrating output from different export formats

### Tabulating Orbitals
- **File**: `_static/placeholder-grid-types.png`
  - **Location**: Tabulating Orbitals section
  - **Description**: Visualization comparing spherical and cartesian grid layouts
  - **Shows**: Diagram or visualization showing spherical vs cartesian grid point distribution

- **File**: `_static/placeholder-orbital-series.png`
  - **Location**: After tabulation examples
  - **Description**: Series of molecular orbitals showing HOMO-LUMO gap region
  - **Shows**: Multiple orbital visualizations arranged in sequence (e.g., HOMO-2, HOMO-1, HOMO, LUMO, LUMO+1)

### Exporting Volumetric Data
- **File**: `_static/placeholder-export-formats.png`
  - **Location**: Export Format Comparison subsection
  - **Description**: Comparison of VTK and cube file structures and visualization in different software
  - **Shows**: Split screen with ParaView (VTK) and another viewer (cube), possibly with file structure diagrams

- **File**: `_static/placeholder-batch-export.png`
  - **Location**: Batch Export Workflow subsection
  - **Description**: Workflow diagram showing grid creation, visualization, and export steps
  - **Shows**: Flowchart or diagram illustrating the complete workflow from parsing to export

## Configuration Reference (`configuration.rst`)

### Visualization Settings
- **File**: `_static/placeholder-background-comparison.png`
  - **Location**: Visualization Settings section
  - **Description**: Side-by-side comparison of white vs black background rendering
  - **Shows**: Same molecule rendered with different background colors

### Bond Settings
- **File**: `_static/placeholder-bond-styles.png`
  - **Location**: Bond Settings section
  - **Description**: Examples showing uniform vs split bond coloring and different radii
  - **Shows**: Multiple renderings demonstrating different bond visualization options

### Grid Settings
- **File**: `_static/placeholder-grid-settings-dialog.png`
  - **Location**: Grid Settings section
  - **Description**: Grid settings dialog showing spherical and cartesian grid configuration options
  - **Shows**: Dialog window with grid type selection and resolution controls

- **File**: `_static/placeholder-grid-resolution.png`
  - **Location**: After grid resolution examples
  - **Description**: Comparison of low vs high resolution grid effects on orbital visualization
  - **Shows**: Side-by-side orbitals showing smooth (high-res) vs blocky (low-res) isosurfaces

### Molecular Orbital Settings
- **File**: `_static/placeholder-colormap-comparison.png`
  - **Location**: Molecular Orbital Settings section
  - **Description**: Comparison of different colormaps applied to the same molecular orbital
  - **Shows**: Multiple versions of same orbital with different color schemes (bwr, RdBu, seismic, custom)

## Summary

### Existing Placeholders (2)
1. `placeholder-cli.png` - CLI interface screenshot
2. `placeholder-python.png` - Python API plotter screenshot

### New Placeholders Added (17)
1. `placeholder-full-viewer.png` - Complete viewer interface
2. `placeholder-molecule-only.png` - Molecule-only mode
3. `placeholder-example-co.png` - CO molecule example
4. `placeholder-examples-grid.png` - All example molecules grid
5. `placeholder-plotter-comparison.png` - Full vs molecule-only comparison
6. `placeholder-controls-panel.png` - Control panel close-up
7. `placeholder-export-menu.png` - Export menu
8. `placeholder-data-export-dialog.png` - Data export dialog
9. `placeholder-image-export-dialog.png` - Image export dialog
10. `placeholder-exported-comparison.png` - Export format comparison
11. `placeholder-grid-types.png` - Grid type visualization
12. `placeholder-orbital-series.png` - Orbital series
13. `placeholder-export-formats.png` - VTK vs cube comparison
14. `placeholder-batch-export.png` - Export workflow diagram
15. `placeholder-background-comparison.png` - Background color comparison
16. `placeholder-bond-styles.png` - Bond style examples
17. `placeholder-grid-settings-dialog.png` - Grid settings dialog
18. `placeholder-grid-resolution.png` - Resolution comparison
19. `placeholder-colormap-comparison.png` - Colormap comparison

**Total Placeholders: 19** (2 existing + 17 new)

## Guidelines for Creating Screenshots

When creating actual screenshots to replace these placeholders:

1. **Resolution**: Use consistent resolution (suggest 1920x1080 or scaled down to 1280x720)
2. **Format**: PNG format with transparent backgrounds where appropriate
3. **Consistency**: Use the same example molecule (benzene recommended) across related screenshots
4. **Annotations**: Consider adding arrows or labels for complex UI elements
5. **Quality**: Ensure good contrast and readability of text elements
6. **Size**: Optimize file sizes while maintaining quality (use PNG compression)

## Priority Order for Screenshot Creation

### High Priority (Core Features)
1. `placeholder-cli.png` and `placeholder-python.png` (update existing)
2. `placeholder-full-viewer.png` - Main interface
3. `placeholder-controls-panel.png` - Core UI element
4. `placeholder-export-menu.png` - Key feature
5. `placeholder-data-export-dialog.png` - Important export workflow
6. `placeholder-image-export-dialog.png` - Important export workflow

### Medium Priority (Feature Demonstrations)
1. `placeholder-plotter-comparison.png` - Shows flexibility
2. `placeholder-grid-types.png` - Shows technical capability
3. `placeholder-orbital-series.png` - Shows typical use case
4. `placeholder-colormap-comparison.png` - Visual customization
5. `placeholder-background-comparison.png` - Visual customization

### Lower Priority (Additional Context)
1. `placeholder-molecule-only.png` - Secondary mode
2. `placeholder-example-co.png` - One of many examples
3. `placeholder-examples-grid.png` - Overview
4. `placeholder-exported-comparison.png` - Format comparison
5. `placeholder-export-formats.png` - Technical details
6. `placeholder-batch-export.png` - Workflow diagram
7. `placeholder-bond-styles.png` - Configuration example
8. `placeholder-grid-settings-dialog.png` - Configuration dialog
9. `placeholder-grid-resolution.png` - Technical comparison
