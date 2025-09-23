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

Skip molecular orbital parsing when you only need the structure:

.. code-block:: python

   parser = Parser('molden.inp', only_molecule=True)

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
