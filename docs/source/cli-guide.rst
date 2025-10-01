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
