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

Data models
-----------

.. autoclass:: moldenViz.models.Atom
   :members:

.. autoclass:: moldenViz.models.AtomType
   :members:

.. autoclass:: moldenViz.models.GaussianPrimitive

.. autoclass:: moldenViz.models.MolecularOrbital
   :members:

.. autoclass:: moldenViz.models.Shell

.. data:: moldenViz.parser.BOHR_PER_ANGSTROM

   Number of Bohr in one Angstrom.

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
