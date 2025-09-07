API Reference
=============

This section provides detailed documentation for all modules, classes, and functions in moldenViz.

Parser Module
-------------

The parser module provides functionality to read and parse molden files.

.. py:class:: Parser(source, only_molecule=False)

   Parser for molden files.

   :param source: The path to the molden file, or the lines from the file.
   :type source: str | list[str]
   :param only_molecule: Only parse the atoms and skip molecular orbitals. Default is ``False``.
   :type only_molecule: bool, optional

   **Attributes:**

   * ``atoms`` (list): A list of Atom objects containing the label, atomic number, and position for each atom.
   * ``mos`` (list): A list of molecular orbital objects (only available when ``only_molecule=False``).

   **Example Usage:**

   .. code-block:: python

      from moldenViz import Parser
      
      # Parse from file
      parser = Parser('molden.inp')
      
      # Access atoms and molecular orbitals
      atoms = parser.atoms
      mos = parser.mos
      
      # Parse only molecule structure (skip MOs)
      parser = Parser('molden.inp', only_molecule=True)

Plotter Module
--------------

The plotter module provides 3D visualization capabilities for molecules and molecular orbitals.

.. py:class:: Plotter(source, only_molecule=False, tabulator=None)

   Create interactive 3D plots of molecules and molecular orbitals.

   :param source: The path to the molden file, lines from the file, or a Parser object.
   :type source: str | list[str] | Parser
   :param only_molecule: Only plot the molecule structure without orbitals. Default is ``False``.
   :type only_molecule: bool, optional
   :param tabulator: Pre-configured tabulator object for custom grids.
   :type tabulator: Tabulator, optional

   **Example Usage:**

   .. code-block:: python

      from moldenViz import Plotter
      
      # Plot molecule with orbitals
      Plotter('molden.inp')
      
      # Plot only the molecule structure
      Plotter('molden.inp', only_molecule=True)
      
      # Use example molecules
      from moldenViz.examples import benzene
      Plotter(benzene)

Tabulator Module
----------------

The tabulator module provides functionality to create grids and tabulate molecular orbitals on those grids.

.. py:class:: Tabulator(source, only_molecule=False)

   Create grids and tabulate Gaussian-type orbitals (GTOs) and molecular orbitals.

   :param source: The path to the molden file, lines from the file, or a Parser object.
   :type source: str | list[str] | Parser
   :param only_molecule: Only load molecule structure. Cannot create grids when True.
   :type only_molecule: bool, optional

   **Attributes:**

   * `grid` (numpy.ndarray): The generated grid points.
   * `gtos` (numpy.ndarray): Tabulated Gaussian-type orbital data on the grid.

   **Methods:**

   .. py:method:: spherical_grid(r, theta, phi)

      Create a spherical coordinate grid.

      :param r: Radial distances.
      :param theta: Polar angles.
      :param phi: Azimuthal angles.

   .. py:method:: cartesian_grid(x, y, z)

      Create a Cartesian coordinate grid.

      :param x: X coordinates.
      :param y: Y coordinates.
      :param z: Z coordinates.

   .. py:method:: tabulate_mos(indices=None)

      Tabulate molecular orbitals on the current grid.

      :param indices: Specific orbital indices to tabulate. If None, tabulates all.
      :type indices: int | list[int] | range, optional
      :returns: Tabulated molecular orbital data.
      :rtype: numpy.ndarray

   **Example Usage:**

   .. code-block:: python

      from moldenViz import Tabulator
      import numpy as np
      
      # Create tabulator
      tab = Tabulator('molden.inp')
      
      # Create spherical grid
      tab.spherical_grid(
          r=np.linspace(0, 5, 20),
          theta=np.linspace(0, np.pi, 20),
          phi=np.linspace(0, 2 * np.pi, 40)
      )
      
      # Tabulate molecular orbitals
      mo_data = tab.tabulate_mos([0, 1, 2])  # First three orbitals

Examples Module
---------------

The examples module provides pre-defined molecular structures for testing and demonstration.

**Available Examples:**

* `co`
* `o2`
* `co2`
* `h2o`
* `benzene`
* `prismane`
* `pyridine`
* `furan`
* `acrolein`

**Example Usage:**

.. code-block:: python

   from moldenViz.examples import benzene, h2o
   from moldenViz import Plotter
   
   # Use example molecules
   Plotter(benzene)
   Plotter(h2o, only_molecule=True)
