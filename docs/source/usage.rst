Usage Guide
===========

This page provides comprehensive usage examples for the ``moldenViz`` library.

Installation
------------

Install moldenViz from PyPI:

.. code-block:: bash

   pip install moldenViz

Prerequisites
~~~~~~~~~~~~~

If you want to use the plotter, make sure Python has access to ``tkinter``:

.. code-block:: bash

   python3 -m tkinter

If Python doesn't have access to ``tkinter``, you can install it:

**macOS:**

.. code-block:: bash

   brew install python-tk

**Ubuntu:**

.. code-block:: bash

   sudo apt-get install python-tk

Command Line Interface
----------------------

Basic Usage
~~~~~~~~~~~

Run moldenViz with a molden file:

.. code-block:: bash

   moldenViz molden.inp

Use example molecules:

.. code-block:: bash

   moldenViz -e co

Plot only the molecule structure (without orbitals):

.. code-block:: bash

   moldenViz molden.inp -m

Get help:

.. code-block:: bash

   moldenViz -h

Available Examples
~~~~~~~~~~~~~~~~~~

The following example molecules are available:

- ``co``
- ``o2`` 
- ``co2``
- ``h2o``
- ``benzene``
- ``prismane``
- ``pyridine``
- ``furan``
- ``acrolein``

Python API
----------

Basic Parsing
~~~~~~~~~~~~~

Parse a molden file and access molecular data:

.. code-block:: python

   from moldenViz import Parser
   
   # Parse from file
   parser = Parser('molden.inp')
   
   # Access atoms and molecular orbitals
   atoms = parser.atoms
   mos = parser.mos
   
   # Parse only molecule structure (skip MOs)
   parser = Parser('molden.inp', only_molecule=True)

Using Examples in Python
~~~~~~~~~~~~~~~~~~~~~~~~~

Work with example molecules:

.. code-block:: python

   from moldenViz import Plotter
   from moldenViz.examples import co
   
   # Plot example molecule
   Plotter(co)

Plotting Molecules
~~~~~~~~~~~~~~~~~~

Basic plotting:

.. code-block:: python

   from moldenViz import Plotter
   
   # Plot molecule with orbitals
   Plotter('molden.inp')
   
   # Plot only the molecule structure
   Plotter('molden.inp', only_molecule=True)

Tabulating Orbitals
~~~~~~~~~~~~~~~~~~~

Create grids and tabulate molecular orbitals:

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
   
   # Create cartesian grid
   tab.cartesian_grid(
       x=np.linspace(-2, 2, 20),
       y=np.linspace(-2, 2, 20),
       z=np.linspace(-2, 2, 20)
   )
   
   # Check grid and GTO data
   print(tab.grid.shape)
   print(tab.gtos.shape)

Tabulating Molecular Orbitals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Tabulate specific molecular orbitals:

.. code-block:: python

   # Tabulate single orbital
   mo_data = tab.tabulate_mos(0)
   
   # Tabulate multiple orbitals by list
   mo_data = tab.tabulate_mos([0, 1, 4])
   
   # Tabulate range of orbitals
   mo_data = tab.tabulate_mos(range(1, 10, 2))
   
   # Tabulate all orbitals
   mos_data = tab.tabulate_mos()

Advanced Usage
--------------

Custom Tabulator with Plotter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use a pre-configured tabulator with the plotter:

.. code-block:: python

   from moldenViz import Tabulator, Plotter
   import numpy as np
   
   # Create custom tabulator
   tab = Tabulator('molden.inp')
   tab.cartesian_grid(
       x=np.linspace(-3, 3, 30),
       y=np.linspace(-3, 3, 30), 
       z=np.linspace(-3, 3, 30)
   )
   
   # Use with plotter
   Plotter('molden.inp', tabulator=tab)

Working with Parser Data
~~~~~~~~~~~~~~~~~~~~~~~~

Access detailed molecular structure information:

.. code-block:: python

   from moldenViz import Parser
   
   parser = Parser('molden.inp')
   
   # Access atoms
   for atom in parser.atoms:
       print(f"Atom: {atom.label}, Position: {atom.position}")
       
   # Access molecular orbitals
   for i, mo in enumerate(parser.mos):
       print(f"MO {i}: Energy = {mo.energy}, Symmetry = {mo.sym}")
       
   # Access shells and basis functions
   for atom in parser.atoms:
       for shell in atom.shells:
           print(f"Shell l={shell.l}, GTOs={len(shell.gtos)}")

Error Handling
--------------

Common issues and solutions:

**File not found:**

.. code-block:: python

   try:
       parser = Parser('nonexistent.inp')
   except FileNotFoundError:
       print("Molden file not found")

**Invalid molden format:**

.. code-block:: python

   try:
       parser = Parser('invalid.inp')
   except ValueError as e:
       print(f"Invalid molden file: {e}")

**Grid creation with only_molecule:**

.. code-block:: python

   tab = Tabulator('molden.inp', only_molecule=True)
   try:
       tab.cartesian_grid(x, y, z)  # This will fail
   except RuntimeError:
       print("Cannot create grids when only_molecule=True")
