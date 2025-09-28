Configuration Reference
=======================

``moldenViz`` looks for a TOML configuration file at ``~/.config/moldenViz/config.toml``. Any values you set there override the built-in defaults.

Default Configuration
---------------------

The full set of built-in defaults is bundled with the package.  Use it as a
reference when you only want to override a handful of options:

.. literalinclude:: ../../src/moldenViz/default_configs/config.toml
   :language: toml
   :caption: Default configuration shipped with moldenViz

Configuration Files
-------------------

Create the directory (if needed) and add a config file:

.. code-block:: bash

   mkdir -p ~/.config/moldenViz
   $EDITOR ~/.config/moldenViz/config.toml

Bond Settings
-------------

Tune how bonds are rendered:

.. code-block:: toml

   [molecule.bond]
   show = true          # Show or hide bonds entirely
   max_length = 4.0     # Maximum bond length to display
   color_type = 'uniform'  # 'uniform' or 'split'
   color = 'grey'       # Hex string or colour name
   radius = 0.15        # Cylinder radius

Switch to per-atom colouring and thicker cylinders:

.. code-block:: toml

   [molecule.bond]
   color_type = 'split'
   radius = 0.25

Grid Settings
-------------

Control the grid used when tabulating molecular orbitals:

.. code-block:: toml

   [grid]
   min_radius = 5
   max_radius_multiplier = 2

   [grid.spherical]
   num_r_points = 100
   num_theta_points = 60
   num_phi_points = 120

   [grid.cartesian]
   num_x_points = 100
   num_y_points = 100
   num_z_points = 100

Increase the resolution in the cartesian grid only:

.. code-block:: toml

   [grid.cartesian]
   num_x_points = 150
   num_y_points = 150
   num_z_points = 150

Molecular Orbital Settings
--------------------------

Adjust the appearance of orbital isosurfaces:

.. code-block:: toml

   [mo]
   contour = 0.1
   opacity = 1.0

Atom Settings
-------------

Drive molecule-wide and per-element rendering tweaks:

.. code-block:: toml

   [molecule]
   opacity = 0.9

   [molecule.atom]
   show = true

Per-atom overrides use the ``Atom.<atomic number>`` table name:

.. code-block:: toml

   [Atom.1]  # Hydrogen
   color = "FFFFFF"
   radius = 0.25

   [Atom.6]  # Carbon
   color = "404040"
   radius = 0.7

   [Atom.8]  # Oxygen
   color = "FF4040"
   radius = 0.6

Control the maximum bonding for a specific element (for example, iron):

.. code-block:: toml

   [Atom.26]
   max_num_bonds = 6

Worked Examples
---------------

The snippets below combine the most common overrides so you can copy them into ``~/.config/moldenViz/config.toml`` and adjust as needed.

Reduce bond clutter on dense organic rings while highlighting a metal centre:

.. code-block:: toml

   [molecule]
   opacity = 0.85

   [molecule.bond]
   max_length = 2.0
   radius = 0.08
   color_type = 'split'

   [Atom.6]  # Carbon
   color = "2A2A2A"
   radius = 0.7

   [Atom.8]  # Oxygen
   color = "FF4040"
   radius = 0.6

   [Atom.26]  # Iron centre
   color = "FFB347"
   radius = 1.1
   max_num_bonds = 4

Tighten radius settings for a very large biomolecule to keep rendering performant:

.. code-block:: toml

   [grid.cartesian]
   num_x_points = 80
   num_y_points = 80
   num_z_points = 80

   [mo]
   contour = 0.07
   opacity = 0.6

   [Atom.1]
   radius = 0.25

   [Atom.6]
   radius = 0.65
   max_num_bonds = 3

   [molecule.bond]
   radius = 0.05
   show = true

For troubleshooting tips tied to these overrides, see :doc:`Troubleshooting <troubleshooting>`, especially the sections on configuration validation and export errors.

Complete Example
----------------

Combine several tweaks in one file:

.. code-block:: toml

   smooth_shading = true

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

   [mo]
   contour = 0.05
   opacity = 0.8

   [molecule]
   opacity = 0.9

   [molecule.bond]
   show = true
   max_length = 3.5
   color_type = 'split'
   radius = 0.12

   [Atom.1]
   color = "FFFFFF"
   radius = 0.3

   [Atom.6]
   color = "303030"
   radius = 0.75

   [Atom.7]
   color = "0060FF"
   radius = 0.65

   [Atom.8]
   color = "FF3030"
   radius = 0.6

Applying Custom Settings
------------------------

The CLI and API automatically load your config file:

.. code-block:: python

   from moldenViz import Plotter
   from moldenViz.examples import benzene

   Plotter(benzene)

No extra arguments are required—the overrides take effect as soon as the file exists.
