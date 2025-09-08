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

- `co`
- `o2` 
- `co2`
- `h2o`
- `benzene`
- `prismane`
- `pyridine`
- `furan`
- `acrolein`

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

Configuration
-------------

moldenViz supports customization through configuration files. The configuration system uses a hierarchical approach where custom settings override defaults.

Configuration Files
~~~~~~~~~~~~~~~~~~~

**Custom Configuration**: User-specific settings at ``~/.config/moldenViz/config.toml``

Available Configuration Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bond Settings
^^^^^^^^^^^^^

Control bond appearance and behavior:

.. code-block:: toml

   [molecule.bond]
   show = true                    # Show/hide bonds
   max_length = 4.0              # Maximum bond length to display
   color_type = 'uniform'        # 'uniform' or 'split' coloring
   color = 'grey'                # Bond color (hex code or color name)
   radius = 0.15                 # Bond radius/thickness

**Bond Color Types:**

- ``uniform``: All bonds use the same color
- ``split``: Bonds are colored based on the atoms they connect

**Example - Change bond color and thickness:**

.. code-block:: toml

   [molecule.bond]
   radius = 0.25
   color_type = 'split'

^^^^^^^^^^^^^

Control grid generation for molecular orbital calculations:

.. code-block:: toml

   [grid]
   min_radius = 5                # Minimum radius for grid generation
   max_radius_multiplier = 2     # Grid extends to max_radius_multiplier * molecular_size

   [grid.spherical]
   num_r_points = 100           # Radial grid points
   num_theta_points = 60        # Theta grid points  
   num_phi_points = 120         # Phi grid points

   [grid.cartesian]
   num_x_points = 100           # X-axis grid points
   num_y_points = 100           # Y-axis grid points
   num_z_points = 100           # Z-axis grid points

**Example - Increase grid resolution:**

.. code-block:: toml

   [grid]
   min_radius = 3
   
   [grid.cartesian]
   num_x_points = 150
   num_y_points = 150
   num_z_points = 150

Molecular Orbital Settings
^^^^^^^^^^^^^^^^^^^^^^^^^^

Control MO visualization:

.. code-block:: toml

   [mo]
   contour = 0.1                # Isosurface contour value
   opacity = 1.0                # MO opacity (0.0 to 1.0)

Atom Settings
^^^^^^^^^^^^^

Control atom display:

.. code-block:: toml

   [molecule.atom]
   show = true                  # Show/hide atoms

   [molecule]
   opacity = 1.0                # Overall molecule opacity

**Example - Semi-transparent molecule:**

.. code-block:: toml

   [molecule]
   opacity = 0.7

   [mo]
   opacity = 0.8

Customizing Atom Types
~~~~~~~~~~~~~~~~~~~~~~

You can customize the appearance of specific atom types by their atomic number:

.. code-block:: toml

   [Atom.1]        # Hydrogen (atomic number 1)
   color = "FF0000"           # Red (hex without #)
   radius = 0.3

   [Atom.6]        # Carbon (atomic number 6)  
   color = "00FF00"           # Green
   radius = 0.8

   [Atom.8]        # Oxygen (atomic number 8)
   color = "0000FF"           # Blue
   radius = 0.6

**Available atom properties:**

- ``name``: Atom symbol (e.g., 'H', 'C', 'O')
- ``color``: Hex color code without # (e.g., 'FF0000' for red)
- ``radius``: Atom display radius (positive float)
- ``max_num_bonds``: Maximum bonds the atom can form

**Example - Colorful atom scheme:**

.. code-block:: toml

   [Atom.1]    # Hydrogen - bright white
   color = "FFFFFF"
   radius = 0.25

   [Atom.6]    # Carbon - dark gray  
   color = "404040"
   radius = 0.7

   [Atom.7]    # Nitrogen - blue
   color = "0080FF"
   radius = 0.65

   [Atom.8]    # Oxygen - red
   color = "FF4040"
   radius = 0.6

Complete Configuration Example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here's a complete custom configuration file that demonstrates various customizations:

.. code-block:: toml

   # ~/.config/moldenViz/config.toml
   
   # Enable smooth shading
   smooth_shading = true
   
   # Grid settings for higher resolution
   [grid]
   min_radius = 3
   max_radius_multiplier = 2.5
   
   [grid.spherical]
   num_r_points = 120
   num_theta_points = 80
   num_phi_points = 160
   
   [grid.cartesian]
   num_x_points = 120
   num_y_points = 120
   num_z_points = 120
   
   # Molecular orbital settings
   [mo]
   contour = 0.05
   opacity = 0.8
   
   # Molecule appearance
   [molecule]
   opacity = 0.9
   
   [molecule.atom]
   show = true
   
   [molecule.bond]
   show = true
   max_length = 3.5
   color_type = 'split'
   radius = 0.12
   
   # Custom atom colors and sizes
   [Atom.1]    # Hydrogen

   radius = 0.3
   
   [Atom.6]    # Carbon
   color = "303030"
   radius = 0.75
   
   [Atom.7]    # Nitrogen  
   color = "0060FF"
   radius = 0.65
   
   [Atom.8]    # Oxygen
   color = "FF3030"
   radius = 0.6

Using Custom Configurations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you create a custom configuration file, moldenViz will automatically load and apply your settings:

.. code-block:: python

   from moldenViz import Plotter
   from moldenViz.examples import benzene
   
   # Uses your custom configuration automatically
   Plotter(benzene)

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

**Invalid configuration:**

.. code-block:: python

   # If your ~/.config/moldenViz/config.toml has invalid settings
   try:
       from moldenViz import Plotter
       Plotter('molden.inp')
   except ValueError as e:
       print(f"Configuration error: {e}")
       # Fix the configuration file and try again
