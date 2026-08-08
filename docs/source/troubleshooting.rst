Troubleshooting
===============

Use this page to diagnose common issues when running ``moldenViz``.

Parser Exceptions
-----------------

- ``ValueError: Unsupported basis functions`` – the parser only accepts spherical Gaussian functions. Re-export your Molden file with spherical GTOs enabled.
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

Invalid entries in ``~/.config/moldenViz/config.toml`` raise ``ValueError`` the
next time you construct a viewer:

.. code-block:: python

   from moldenViz import Plotter
   
   try:
       Plotter('molden.inp')
   except ValueError as exc:
       print(f'Configuration error: {exc}')
       print('Review your TOML configuration and try again')

Repeated Slow Startup
---------------------

PyVista uses Matplotlib for colors and text. The first launch may take longer
while Matplotlib creates its font cache, but later launches should reuse it. If
every launch prints ``Matplotlib is building the font cache`` or a fontconfig
error, configure a persistent writable cache directory before starting
``moldenViz``:

.. code-block:: bash

   export MPLCONFIGDIR=/path/to/a/writable/matplotlib-cache

Use a directory owned by the current user and keep it between launches. A
temporary directory avoids the write error but forces the cache to be rebuilt
when that directory is removed.

Qt Application Errors
---------------------

``RuntimeError: OrbitalViewer requires an existing QApplication`` means the
embeddable widget was constructed before its host initialized PySide6. Create
``QApplication`` first, then add the viewer to a layout. Do not call
``app.exec()`` from the viewer or from a button handler; the host application
owns that loop.

If an embedded page is rebuilt, call ``viewer.close()`` during teardown. This
cancels result delivery and explicitly releases the VTK render window before
Qt destroys the surrounding layout.

Headless Qt Tests
-----------------

``pyvistaqt.QtInteractor`` requires a working OpenGL context. Constructing a
real ``OrbitalViewer`` with ``QT_QPA_PLATFORM=offscreen`` can terminate the
process inside VTK before Python can raise an exception. moldenViz detects this
unsupported combination and raises ``RuntimeError`` first.

For page construction, layout, and lifecycle smoke tests, use the supported
non-rendering context:

.. code-block:: python

   from moldenViz.testing import without_rendering

   def test_orbital_page(qapp):
       with without_rendering():
           page = build_orbital_page()
           assert page.viewer is not None
           page.close()

The contained ``NullInteractor`` records basic scene operations but does not
create a VTK render window. Rendering and image-output tests still require a
platform plugin with a real OpenGL context, such as a virtual display.

Export Errors
-------------

- ``ValueError: Orbital selection out of bounds`` – ensure the indices passed to ``--orbitals`` or ``Tabulator.export_*`` fall within the available range reported by ``Parser.mos``.
- ``RuntimeError: Tabulator grid is undefined`` – create a grid (cartesian or spherical) before calling an export method.
- ``RuntimeError: Cube export requires a cartesian grid`` – Gaussian cube files expect a rectilinear grid; re-run the export with ``--grid cartesian`` or ``--export-cube`` only.

If a problem persists, run ``moldenViz -h`` to confirm the CLI supports the options you are using and check the :doc:`Configuration Reference <configuration>` for grid defaults.
