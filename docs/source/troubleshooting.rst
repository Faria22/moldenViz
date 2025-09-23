Troubleshooting
===============

Use this page to diagnose common issues when running ``moldenViz``.

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
