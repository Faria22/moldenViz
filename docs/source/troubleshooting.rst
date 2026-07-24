Troubleshooting
===============

Use this page to diagnose common issues when running ``moldenViz``.

Parser Exceptions
-----------------

- ``ValueError: Mixed spherical and Cartesian basis declarations are not supported.`` – the file contains both spherical ``5D``/``7F``/``9G`` and Cartesian ``6D``/``10F``/``15G`` markers. Re-export it with one consistent basis representation.
- ``ValueError: No supported ... basis declaration found.`` – add the standard Molden spherical or Cartesian basis markers emitted by your quantum-chemistry producer.
- ``ValueError: Invalid shell label`` – occurs when the ``[GTO]`` section contains unexpected angular momentum labels. Confirm the file adheres to the `Molden format specification <https://www.theochem.ru.nl/molden/molden_format.html>`_.
- ``RuntimeError: Grid creation is not allowed when `only_molecule` is set to `True`.`` – raised if you request grids or exports while skipping orbitals. Re-run without ``--only-molecule``.

File Not Found
--------------

.. code-block:: python

   from moldenViz import Parser
   
   try:
       Parser('missing-file.inp')
   except FileNotFoundError:
       print('Molden file not found')

Invalid Molden Format
---------------------

.. code-block:: python

   from moldenViz import Parser
   
   try:
       Parser('invalid.inp')
   except ValueError as exc:
       print(f'Invalid molden file: {exc}')

Grids With ``only_molecule=True``
---------------------------------

``Tabulator`` cannot build grids when you skip molecular orbital data:

.. code-block:: python

   from moldenViz import Tabulator
   
   tab = Tabulator('molden.inp', only_molecule=True)
   
   try:
       tab.cartesian_grid(x, y, z)
   except RuntimeError:
       print('Cannot create grids when only_molecule=True')

Configuration Errors
--------------------

Invalid entries in ``~/.config/moldenViz/config.toml`` raise ``ValueError`` the next time you import plotting classes:

.. code-block:: python

   from moldenViz import Plotter
   
   try:
       Plotter('molden.inp')
   except ValueError as exc:
       print(f'Configuration error: {exc}')
       print('Review your TOML configuration and try again')

Export Errors
-------------

- ``ValueError: Orbital selection out of bounds`` – ensure the indices passed to ``--orbitals`` or ``Tabulator.export_*`` fall within the available range reported by ``Parser.mos``.
- ``RuntimeError: Tabulator grid is undefined`` – create a grid (cartesian or spherical) before calling an export method.
- ``RuntimeError: Cube export requires a cartesian grid`` – Gaussian cube files expect a rectilinear grid; re-run the export with ``--grid cartesian`` or ``--export-cube`` only.

If a problem persists, run ``moldenViz -h`` to confirm the CLI supports the options you are using and check the :doc:`Configuration Reference <configuration>` for grid defaults.
