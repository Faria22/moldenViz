Python API
==========

Use the Python API when you need to integrate ``moldenViz`` into scripts or notebooks, or when you need precise control over parsing and tabulating molecular orbitals.

Parsing Molden Files
--------------------

Read a Molden file and access its atoms and orbitals:

.. code-block:: python

   from moldenViz import Parser

   parser = Parser(filename='molden.inp')

   atoms = parser.atoms
   mos = parser.mos

To parse data already held in memory, pass the complete file as one string.
``filename`` and ``content`` are mutually exclusive:

.. code-block:: python

   from pathlib import Path

   content = Path('molden.inp').read_text(encoding='utf-8')
   parser = Parser(content=content)

``Tabulator`` and ``Plotter`` accept the same explicit input keywords. When a
configured tabulator is supplied to ``Plotter``, omit both input keywords and
use ``Plotter(tabulator=tab)``.

The returned objects use public parser model types, so they can also be
imported for annotations or construction:

.. code-block:: python

   from moldenViz import Atom, GaussianPrimitive, MolecularOrbital, Shell

Use ``moldenViz.parser.BOHR_PER_ANGSTROM`` when converting Angstrom coordinates
to the Bohr units used by parsed atom positions.

Skip molecular orbital parsing when you only need the structure:

.. code-block:: python

   parser = Parser(filename='molden.inp', only_molecule=True)

Choose whether molecular orbitals are sorted by energy or retain their order
in the source file:

.. code-block:: python

   energy_ordered = Parser(filename='molden.inp', mo_order='energy')
   file_ordered = Parser(filename='molden.inp', mo_order='file')

Example Molecules
-----------------

The ``moldenViz.examples`` module bundles molecules you can use without providing your own Molden file:

.. code-block:: python

   from moldenViz import Plotter
   from moldenViz.examples import co

   Plotter(content=co)

Plotting Molecules
------------------

The ``Plotter`` class renders atoms, bonds, and (optionally) orbital isosurfaces:

.. code-block:: python

   from moldenViz import Plotter

   # Plot molecule with orbitals
   Plotter(filename='molden.inp')

   # Plot only the molecular structure
   Plotter(filename='molden.inp', only_molecule=True)

``Plotter`` is the free-floating convenience API. It creates and runs a
``QApplication`` when no Qt application exists. Inside an existing PySide6
application it reuses that application, shows its window, and returns
immediately; retain the returned window for as long as it should remain open.
When ``Plotter`` owns the application, it paints a loading placeholder before
initializing the VTK-backed viewer so the standalone window appears promptly.

Embedding the Viewer
~~~~~~~~~~~~~~~~~~~~

Use the plain ``QWidget`` API when another PySide6 application owns the layout
and event loop:

.. code-block:: python

   from PySide6.QtWidgets import QApplication, QVBoxLayout, QWidget
   from moldenViz.qt import OrbitalViewer

   app = QApplication.instance() or QApplication([])
   page = QWidget()
   layout = QVBoxLayout(page)

   viewer = OrbitalViewer(parent=page)
   layout.addWidget(viewer)
   viewer.set_input(filename='molden.inp')

   page.show()
   app.exec()

Constructing ``OrbitalViewer`` never creates an application, opens a window,
or starts an event loop. Call ``viewer.close()`` before a page containing it is
discarded so the VTK render window is released promptly.

Host-provided Controls
~~~~~~~~~~~~~~~~~~~~~~

Hide the built-in panel when the host application supplies its own controls.
The viewer remains a single widget, but its selection, grid, appearance, and
export operations do not depend on the panel being visible:

.. code-block:: python

   viewer = OrbitalViewer(filename='molden.inp', parent=page, show_controls=False)

   # Populate a host-owned orbital picker.
   for index, orbital in enumerate(viewer.molecular_orbitals):
       host_picker.addItem(f'{index + 1}: {orbital.sym}')

   host_picker.currentIndexChanged.connect(viewer.show_orbital)
   viewer.update_appearance(contour=0.05, mo_opacity=0.8)
   viewer.set_spherical_grid(
       radius=5.0,
       radial_points=100,
       theta_points=60,
       phi_points=120,
   )

``current_orbital_index`` reports the current selection, with ``-1`` meaning
that no orbital is displayed. ``orbital_changed``, ``loading_changed``, and
``input_ready`` let host controls follow viewer state. Appearance updates are
partial: omitted values retain their current per-viewer setting. Cartesian
grids can be configured with ``set_cartesian_grid``; callers with pre-built
NumPy axes can use ``update_grid`` directly.

The built-in panel can be restored later with
``viewer.set_controls_visible(True)``. It remains synchronized with changes
made through the public viewer API.

For a free-floating viewer inside an application whose event loop is already
running, use ``window = Plotter(filename='molden.inp')``. This is the PySide6 replacement
for the former ``tk_root`` integration argument.

Embedded viewers do not open file or message dialogs. Export operations accept
explicit paths, errors are delivered through ``error_occurred``, and export
buttons emit ``export_requested`` so the host can present its own interface:

.. code-block:: python

   viewer.error_occurred.connect(show_host_error)
   viewer.export_requested.connect(handle_export_request)
   viewer.export_image('scene.png', file_format='png', transparent=True)

When the built-in controls are visible, the host must connect
``export_requested`` for those buttons to perform an action. If no receiver is
connected, the viewer logs a warning and does not open a dialog itself. Hosts
that provide their own controls can call ``export_data`` and ``export_image``
directly with explicit destinations.

Viewer configuration is per instance. Nested constructor overrides take
precedence over the defaults and the user's TOML file without changing another
viewer:

.. code-block:: python

   viewer = OrbitalViewer(
       filename='molden.inp',
       parent=page,
       config={'background_color': '#202124', 'mo': {'opacity': 0.8}},
   )

The visualization-specific ``AtomType`` model is available from the package
root when you need to describe atom display properties:

.. code-block:: python

   from moldenViz import AtomType

Interactive Controls
~~~~~~~~~~~~~~~~~~~~

The plotter window provides several interactive controls:

* **Orbital Selection**: Navigate through molecular orbitals using the control panel
* **Contour Adjustment**: Enter the isosurface contour and apply the change
* **Opacity Control**: Adjust orbital and molecule transparency in increments of 0.1
* **Color Control**: Choose molecular-orbital and background presets or enter custom colors
* **Bond Styling**: Select split or uniform bond colors and adjust radius in increments of 0.05
* **Grid Settings**: Change grid resolution and type (spherical/cartesian)
* **Export Options**: Access data and image export through the menu bar

When ``Plotter`` creates its own grid, it tabulates Gaussian-type orbitals in the background so the molecule window can appear first. Orbital controls become usable when that work finishes; failures are reported in the GUI instead of leaving the viewer silently unavailable.

Exporting from the GUI
~~~~~~~~~~~~~~~~~~~~~~

When using the ``Plotter`` GUI, you can export data or images from the Qt
control panel or menu bar:

1. Open the plotter and select the **Export** tab in the control panel.
2. Configure the data or image options described below.
3. Click the corresponding export button and choose a destination in the Qt
   file dialog.

The top-level **Export** menu provides shortcuts using the options currently
selected in the control panel.

**Exporting Data (Molecular Orbitals)**

1. Choose the data format:

   - **VTK (.vtk)**: Exports one orbital or all orbitals as point-data arrays on a structured grid
   - **Gaussian Cube (.cube)**: Exports a single orbital (cube format does not support multiple orbitals)

2. Select orbital scope:

   - **Current orbital**: Exports the currently displayed orbital
   - **All orbitals**: Exports all molecular orbitals (VTK format only)

3. Click **Export orbital data** and choose the save location.

The export uses the current grid configuration from the plotter, so adjust grid settings before exporting if needed.

**Exporting Images**

1. Choose the image format:

   - **PNG (.png)**: Raster format with optional transparent background
   - **JPEG (.jpg)**: Raster format (no transparency support)
   - **SVG (.svg)**: Vector format for scalable graphics
   - **PDF (.pdf)**: Vector format for publication-quality output

2. For PNG format, optionally enable **Transparent PNG background**.
3. Click **Export image** and choose the save location.

Image exports capture the current view exactly as displayed in the PyVista window, including all visible actors (molecule, orbitals, etc.).

Tabulating Orbitals
-------------------

Use ``Tabulator`` to build grids and evaluate molecular orbitals:

.. code-block:: python

   from moldenViz import Tabulator
   import numpy as np

   tab = Tabulator(filename='molden.inp')

   tab.spherical_grid(
       r=np.linspace(0, 5, 20),
       theta=np.linspace(0, np.pi, 20),
       phi=np.linspace(0, 2 * np.pi, 40)
   )

   tab.cartesian_grid(
       x=np.linspace(-2, 2, 20),
       y=np.linspace(-2, 2, 20),
       z=np.linspace(-2, 2, 20)
   )

Inspect structured-grid metadata through the read-only ``grid_type``,
``grid_dimensions``, and ``grid_axes`` properties. To supply an arbitrary point
grid, use ``set_grid``; direct assignment to ``grid`` is not supported:

.. code-block:: python

   tab.set_grid(np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]))

Tabulate specific molecular orbitals or ranges:

.. code-block:: python

   # Single orbital
   mo_data = tab.tabulate_mos(0)

   # Multiple orbitals by list
   mo_data = tab.tabulate_mos([0, 1, 4])

   # Range of orbitals
   mo_data = tab.tabulate_mos(range(1, 10, 2))

   # All orbitals
   mos_data = tab.tabulate_mos()

Advanced Workflows
------------------

Supply a pre-configured ``Tabulator`` to ``Plotter`` for re-use or fine-grained control over grid resolution:

.. code-block:: python

   from moldenViz import Tabulator, Plotter
   import numpy as np

   tab = Tabulator(filename='molden.inp')
   tab.cartesian_grid(
       x=np.linspace(-3, 3, 30),
       y=np.linspace(-3, 3, 30),
       z=np.linspace(-3, 3, 30)
   )

   Plotter(tabulator=tab)

When orbital rendering is enabled, pass a ``Tabulator`` only after it has a spherical or cartesian grid and cached GTO values. The grid helpers tabulate GTOs by default; if you call them with ``tabulate_gtos=False``, call ``tab.tabulate_gtos()`` before constructing ``Plotter``. Molecule-only viewers do not require cached GTOs.

The cartesian grid keeps spacing uniform—ideal for Gaussian cube exports—while the spherical grid matches the viewer defaults and keeps memory usage low for visual inspection. Pick the smallest grid that contains your molecule; doubling every axis multiplies memory use by eight.

GTO Concurrency
---------------

``Tabulator`` reuses one process-wide thread pool for GTO work instead of
constructing a pool for every grid. By default, tabulation uses no more than the
number of atoms, available CPUs, and four workers through 125,000 points. Larger
grids run sequentially because each concurrent atom adds large temporary arrays
and measured peak memory rises sharply at that scale:

.. code-block:: python

   # Default: at most four workers.
   tab = Tabulator(filename='molden.inp')

   # Deterministic sequential execution.
   sequential_tab = Tabulator(filename='molden.inp', max_workers=1)

   # A lower per-tabulation concurrency limit.
   two_worker_tab = Tabulator(filename='molden.inp', max_workers=2)

The configured ceiling is available as ``tab.max_workers``. Supplying an explicit
value overrides the large-grid sequential policy; values above four are still
clamped to the process-wide ceiling, and CPU and atom counts can lower it further.
Larger worker counts can reduce runtime on smaller grids, but each active atom
needs its own temporary coordinate, exponential, and solid-harmonic arrays. Use
``max_workers=1`` when minimizing peak memory matters more than throughput.
Multiple ``Tabulator`` instances share the same four-worker executor, so
simultaneous viewer and export work cannot create nested per-call pools or exceed
the documented process-wide concurrency bound.

.. _exporting-from-python:

Exporting Volumetric Data (v1.1+)
---------------------------------

You can export orbitals without opening the GUI. Create a grid, tabulate orbitals, and call the export method:

.. code-block:: python

   from moldenViz import Tabulator
   import numpy as np

   tab = Tabulator(filename='molecule.molden')
   tab.cartesian_grid(
       x=np.linspace(-8, 8, 120),
       y=np.linspace(-8, 8, 120),
       z=np.linspace(-8, 8, 120),
   )

   from pathlib import Path

   export_dir = Path('exports')
   export_dir.mkdir(exist_ok=True)

   # Export orbitals 15 and 16 to separate VTK and cube files
   for mo_index in (15, 16):
       tab.export(export_dir / f'orbital_{mo_index}.vtk', mo_index=mo_index)
       tab.export(export_dir / f'orbital_{mo_index}.cube', mo_index=mo_index)

The format-specific ``export_vtk`` and ``export_cube`` methods are also public when you need to call a writer directly. ``export`` is usually simpler because it selects the writer from the destination suffix.

**Export Format Comparison**

The table below compares VTK and Gaussian cube export formats:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Feature
     - VTK Format
     - Gaussian Cube
   * - Multiple orbitals
     - ✓ (one point-data array per orbital)
     - ✗ (single only)
   * - Grid type
     - Spherical or Cartesian
     - Cartesian only
   * - Software support
     - ParaView, VisIt
     - Most quantum chemistry viewers
   * - File size
     - Compact (binary available)
     - Larger (text format)

**Batch Export Workflow**

To reuse tabulation results in a notebook without re-computation:

.. code-block:: python

   tab = Tabulator(filename='molecule.molden')
   tab.spherical_grid(
       r=np.linspace(0, 10, 90),
       theta=np.linspace(0, np.pi, 60),
       phi=np.linspace(0, 2 * np.pi, 120),
   )

   # Keep tabulator to reuse precomputed GTOs
   Plotter(tabulator=tab)

   # Later, export the same grid to VTK
   tab.export('exports/spherical_0.vtk', mo_index=0)

Inspecting Parsed Data
----------------------

Loop over atoms, shells, and orbitals for deeper analysis:

.. code-block:: python

   from moldenViz import Parser

   parser = Parser(filename='molden.inp')

   for atom in parser.atoms:
       print(f"Atom: {atom.label}, Position: {atom.position}")
       for shell in atom.shells:
           print(f"Shell l={shell.l}, GTOs={len(shell.gtos)}")

   for i, mo in enumerate(parser.mos):
       print(f"MO {i}: Energy = {mo.energy}, Symmetry = {mo.sym}")
