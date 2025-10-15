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

This launches the interactive viewer with both the molecule structure and orbital visualization:

.. image:: _static/placeholder-full-viewer.png
   :alt: Full moldenViz viewer showing molecule with orbital isosurface
   :align: center
   :class: screenshot-placeholder

The interactive window allows you to:

* Select and view different molecular orbitals
* Adjust visualization parameters (contour, opacity, colors)
* Switch between spherical and cartesian grids
* Export orbital data or visualization images

Show only the molecular structure (skip orbital surfaces):

.. code-block:: bash

   moldenViz path/to/file.molden -m

.. image:: _static/placeholder-molecule-only.png
   :alt: Molecule-only view without orbital surfaces
   :align: center
   :class: screenshot-placeholder

Bundled Examples
----------------

Use one of the shipped example molecules when you just want to explore the plotting experience:

.. code-block:: bash

   moldenViz -e co

.. image:: _static/placeholder-example-co.png
   :alt: Carbon monoxide molecule with orbital visualization
   :align: center
   :class: screenshot-placeholder

Available examples include: ``co``, ``o2``, ``co2``, ``h2o``, ``benzene``, ``prismane``, ``pyridine``, ``furan``, ``acrolein``

.. image:: _static/placeholder-examples-grid.png
   :alt: Grid showing various example molecules available in moldenViz
   :align: center
   :class: screenshot-placeholder

Tips
----

- Combine ``-e`` with ``-m`` to explore only the geometry of an example.
- Use configuration overrides (:doc:`see details <configuration>`) to change colors, bond lengths, or orbital contours before launching the CLI.
Common recipes
--------------

Use the following patterns as building blocks:

.. code-block:: bash

   # Only plot the molecule inside my.molden
   moldenViz my.molden -m

.. code-block:: bash

   # Plot only the molecule of one of the examples
   moldenViz -e co -m

For additional Python-based export options, see :ref:`exporting-from-python` in the :doc:`Python API guide <python-api>`.
