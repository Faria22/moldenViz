API Reference
=============

Only the entries listed on this page are part of the supported Python API.
Underscored names and modules are implementation details; see
:doc:`API stability <api-stability>`.

Package exports
---------------

``moldenViz`` exports ``Atom``, ``AtomType``, ``GaussianPrimitive``,
``GridType``, ``MolecularOrbital``, ``Parser``, ``Plotter``, ``Shell``,
``Tabulator``, and ``__version__``.

The public ``moldenViz.qt`` module exports ``OrbitalViewer`` for embedding and
``ViewerConfig`` for validated per-instance configuration.

The public ``moldenViz.testing`` module exports ``without_rendering`` and
``NullInteractor`` for Qt host smoke tests that do not exercise VTK rendering.

Parser data models
------------------

.. autoclass:: moldenViz.models.Atom
   :members:

.. autoclass:: moldenViz.models.GaussianPrimitive

.. autoclass:: moldenViz.models.MolecularOrbital
   :members:

.. autoclass:: moldenViz.models.Shell

.. data:: moldenViz.parser.BOHR_PER_ANGSTROM

   Number of Bohr in one Angstrom.

Visualization configuration
---------------------------

``AtomType`` is a visualization-specific model. Its supported import path is
``moldenViz.AtomType``; it is not part of the parser-focused
``moldenViz.models`` namespace.

.. autoclass:: moldenViz.AtomType
   :members:

Parser
------

.. autoclass:: moldenViz.parser.Parser
   :show-inheritance:

Tabulator
---------

.. autoclass:: moldenViz.tabulator.Tabulator
   :members: grid, gtos, has_gtos, atoms, molecular_orbitals, grid_type, grid_dimensions, grid_axes, set_grid, set_gtos, clear_gtos, spherical_to_cartesian, cartesian_to_spherical, cartesian_grid, spherical_grid, tabulate_gtos, tabulate_mos, export, export_vtk, export_cube
   :member-order: bysource
   :show-inheritance:

.. autoclass:: moldenViz.tabulator.GridType
   :members:

Plotter
-------

.. autoclass:: moldenViz.plotter.Plotter
   :members: wait_for_gtos, plot_orbital, toggle_molecule, toggle_atoms, toggle_bonds, is_molecule_visible, are_atoms_visible, are_bonds_visible
   :member-order: bysource
   :show-inheritance:

OrbitalViewer
-------------

.. autoclass:: moldenViz.qt.OrbitalViewer
   :members: config, gtos_ready, axes_visible, controls_visible, current_orbital_index, has_export_handler, molecular_orbitals, set_input, show_orbital, set_axes_visible, set_controls_visible, update_grid, set_spherical_grid, set_cartesian_grid, apply_appearance, update_appearance, set_background_color, export_data, export_image, save_settings, wait_for_gtos, close
   :member-order: bysource
   :show-inheritance:

.. autoclass:: moldenViz.qt.ViewerConfig
   :members:

Testing Qt hosts
----------------

.. autofunction:: moldenViz.testing.without_rendering

.. autoclass:: moldenViz.testing.NullInteractor
