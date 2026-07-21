Python API
==========

Use the Python API when you need to integrate ``moldenViz`` into scripts or notebooks, or when you need precise control over parsing and tabulating molecular orbitals.

Parsing Molden Files
--------------------

Read a Molden file and access its atoms and orbitals:

.. code-block:: python

   from moldenViz import Parser

   parser = Parser('molden.inp')

   atoms = parser.atoms
   mos = parser.mos

The returned objects use public model types, so they can also be imported for
annotations or construction:

.. code-block:: python

   from moldenViz import Atom, AtomType, GaussianPrimitive, MolecularOrbital, Shell

Use ``moldenViz.parser.BOHR_PER_ANGSTROM`` when converting Angstrom coordinates
to the Bohr units used by parsed atom positions.

Skip molecular orbital parsing when you only need the structure:

.. code-block:: python

   parser = Parser('molden.inp', only_molecule=True)

Choose whether molecular orbitals are sorted by energy or retain their order
in the source file:

.. code-block:: python

   energy_ordered = Parser('molden.inp', mo_order='energy')
   file_ordered = Parser('molden.inp', mo_order='file')

Example Molecules
-----------------

The ``moldenViz.examples`` module bundles molecules you can use without providing your own Molden file:

.. code-block:: python

   from moldenViz import Plotter
   from moldenViz.examples import co

   Plotter(co)

Plotting Molecules
------------------

The ``Plotter`` class renders atoms, bonds, and (optionally) orbital isosurfaces:

.. code-block:: python

   from moldenViz import Plotter

   # Plot molecule with orbitals
   Plotter('molden.inp')

   # Plot only the molecular structure
   Plotter('molden.inp', only_molecule=True)

Interactive Controls
~~~~~~~~~~~~~~~~~~~~

The plotter window provides several interactive controls:

* **Orbital Selection**: Navigate through molecular orbitals using the control panel
* **Contour Adjustment**: Modify isosurface contour levels in real-time
* **Opacity Control**: Adjust transparency of orbital surfaces and molecule
* **Grid Settings**: Change grid resolution and type (spherical/cartesian)
* **Export Options**: Access data and image export through the menu bar

When ``Plotter`` creates its own grid, it tabulates Gaussian-type orbitals in the background so the molecule window can appear first. Orbital controls become usable when that work finishes; failures are reported in the GUI instead of leaving the viewer silently unavailable.

Exporting from the GUI
~~~~~~~~~~~~~~~~~~~~~~

When using the ``Plotter`` GUI, you can export data or images from the PyVista plotter window:

1. Open the plotter (the Orbitals window appears automatically when plotting with orbitals enabled)
2. Click the **Export** menu in the PyVista plotter menubar (next to File, View, and Tools)
3. Choose the export type:

   - **Data**: Export molecular orbital data
   - **Image**: Export the current visualization as an image

**Exporting Data (Molecular Orbitals)**

1. Select **Data** from the Export menu
2. Choose your export format:

   - **VTK (.vtk)**: Exports one orbital or all orbitals as point-data arrays on a structured grid
   - **Gaussian Cube (.cube)**: Exports a single orbital (cube format does not support multiple orbitals)

3. Select orbital scope:

   - **Current orbital**: Exports the currently displayed orbital
   - **All orbitals**: Exports all molecular orbitals (VTK format only)

4. Click **Export** and choose the save location

The export uses the current grid configuration from the plotter, so adjust grid settings before exporting if needed.

**Exporting Images**

1. Select **Image** from the Export menu
2. Choose your image format:

   - **PNG (.png)**: Raster format with optional transparent background
   - **JPEG (.jpg)**: Raster format (no transparency support)
   - **SVG (.svg)**: Vector format for scalable graphics
   - **PDF (.pdf)**: Vector format for publication-quality output

3. For PNG format, optionally enable **Transparent background** to remove the background color
4. Click **Export** and choose the save location

Image exports capture the current view exactly as displayed in the PyVista window, including all visible actors (molecule, orbitals, etc.).

Tabulating Orbitals
-------------------

Use ``Tabulator`` to build grids and evaluate molecular orbitals:

.. code-block:: python

   from moldenViz import Tabulator
   import numpy as np

   tab = Tabulator('molden.inp')

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

   tab = Tabulator('molden.inp')
   tab.cartesian_grid(
       x=np.linspace(-3, 3, 30),
       y=np.linspace(-3, 3, 30),
       z=np.linspace(-3, 3, 30)
   )

   Plotter('molden.inp', tabulator=tab)

The cartesian grid keeps spacing uniform—ideal for Gaussian cube exports—while the spherical grid matches the viewer defaults and keeps memory usage low for visual inspection. Pick the smallest grid that contains your molecule; doubling every axis multiplies memory use by eight.

.. _exporting-from-python:

Exporting Volumetric Data (v1.1+)
---------------------------------

You can export orbitals without opening the GUI. Create a grid, tabulate orbitals, and call the export method:

.. code-block:: python

   from moldenViz import Tabulator
   import numpy as np

   tab = Tabulator('molecule.molden')
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

   tab = Tabulator('molecule.molden')
   tab.spherical_grid(
       r=np.linspace(0, 10, 90),
       theta=np.linspace(0, np.pi, 60),
       phi=np.linspace(0, 2 * np.pi, 120),
   )

   # Keep tabulator to reuse precomputed GTOs
   Plotter('molecule.molden', tabulator=tab)

   # Later, export the same grid to VTK
   tab.export('exports/spherical_0.vtk', mo_index=0)

Inspecting Parsed Data
----------------------

Loop over atoms, shells, and orbitals for deeper analysis:

.. code-block:: python

   from moldenViz import Parser

   parser = Parser('molden.inp')

   for atom in parser.atoms:
       print(f"Atom: {atom.label}, Position: {atom.position}")
       for shell in atom.shells:
           print(f"Shell l={shell.l}, GTOs={len(shell.gtos)}")

   for i, mo in enumerate(parser.mos):
       print(f"MO {i}: Energy = {mo.energy}, Symmetry = {mo.sym}")
