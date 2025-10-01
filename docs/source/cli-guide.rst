Command Line
============

The ``moldenViz`` CLI offers a fast way to inspect Molden files and explore bundled examples. Run ``moldenViz -h`` for a full option list; the tables and recipes below summarise each capability.

Primary Options
----------------

.. list-table::
   :header-rows: 1
   :widths: 25 55 20

   * - Flag / Argument
     - Description
     - Default
   * - ``FILE``
     - Positional argument; path to a Molden file to visualise.
     - ``None`` (required unless ``--example`` is used)
   * - ``-e``, ``--example``
     - Load a bundled example molecule (``co``, ``o2``, ``co2``, ``h2o``, ``benzene``, ``prismane``, ``pyridine``, ``furan``, ``acrolein``).
     - None
   * - ``-m``, ``--only-molecule``
     - Skip orbital surfaces and render only the nuclear framework.
     - ``False``
   * - ``--export-vtk PATH``
     - Export selected orbitals as VTK structured grids alongside launching the viewer (v1.1+).
     - Disabled
   * - ``--export-cube PATH``
     - Write selected orbitals to Gaussian cube files (v1.1+).
     - Disabled
   * - ``--orbitals INDICES``
     - Comma-separated list or range expression for orbitals to export (defaults to HOMO/LUMO pair when omitted).
     - Auto-selected
   * - ``--grid cartesian|spherical``
     - Choose the export grid type. The CLI defaults to spherical (cartesian when ``--export-cube`` is used).
     - ``spherical``
   * - ``--resolution Nx,Ny,Nz``
     - Override the grid resolution used for export (interpreted as ``r,theta,phi`` counts for spherical grids).
     - Configuration defaults
   * - ``-h``, ``--help``
     - Display inline help and exit.
     - ---

Basic Rendering
----------------

Render molecular orbitals from a Molden file:

.. code-block:: bash

   moldenViz path/to/file.molden

Show only the molecular structure (skip orbital surfaces):

.. code-block:: bash

   moldenViz path/to/file.molden -m

Bundled Examples
----------------

Use one of the shipped example molecules when you just want to explore the plotting experience:

.. code-block:: bash

   moldenViz -e co

Tips
----

- Combine ``-e`` with ``-m`` to explore only the geometry of an example.
- Use configuration overrides (:doc:`see details <configuration>`) to change colors, bond lengths, or orbital contours before launching the CLI.
- Pass ``--resolution`` to down-sample large grids before exporting to keep file sizes manageable.

Volumetric Export (v1.1+)
-------------------------

Exporting new volumetric data does not block the interactive window; exports run first, then the GUI opens.

.. code-block:: bash

   # Export orbitals 15, 16, and 17 to both VTK and cube files
   moldenViz molecule.molden --orbitals 15,16,17 \
       --export-vtk exports/molecule_{index}.vtk \
       --export-cube exports/molecule_{index}.cube

VTK exports require PyVista (installed with the default extras). Cube files mandate a cartesian grid; the CLI switches automatically when ``--export-cube`` is present. Use braces like ``{index}`` or ``{mo}`` in file names to generate one file per orbital.

You can mix example molecules with export flags as well:

.. code-block:: bash

   moldenViz -e benzene --export-vtk exports/benzene_{index}.vtk --grid cartesian --resolution 80,80,80

Common recipes
--------------

Use the following patterns as building blocks:

.. code-block:: bash

   # Export only the structure (no orbitals) while keeping the GUI open
   moldenViz my.molden -m

.. code-block:: bash

   # Use a spherical grid with 120 radial points and auto-selected orbitals
   moldenViz my.molden --export-cube exports/my_{index}.cube --resolution 120,64,128

.. code-block:: bash

   # Use your ~/.config/moldenViz/config.toml overrides and export to VTK
   moldenViz my.molden --export-vtk exports/my_{index}.vtk

For additional Python-based export options, see :ref:`exporting-from-python` in the :doc:`Python API guide <python-api>`.
