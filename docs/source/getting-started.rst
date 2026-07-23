Getting Started
===============

``moldenViz`` lets you inspect Molden outputs either from the command line or directly in Python. This page covers the essentials so you can install the package and render your first molecule quickly.

What is a Molden file?
----------------------

A Molden file stores the atomic structure and molecular orbital information produced by many quantum-chemistry codes (GAUSSIAN, ORCA, Molpro, and others). Each file contains ``[Atoms]`` coordinates, contracted Gaussian shells in ``[GTO]``, and orbital coefficients under ``[MO]``. ``moldenViz`` currently reads Molden files that use **spherical Gaussian functions**; Cartesian variants (sometimes labelled 5D/9G) raise an error during parsing. For the full format specification, refer to the `official Molden format description <https://www.theochem.ru.nl/molden/molden_format.html>`_.

Installation
------------

Install the core package from PyPI for parsing and tabulation:

.. code-block:: bash

   pip install moldenViz

The core package does not install Qt, PyVista, or other visualization
dependencies. To use the interactive viewer or the ``moldenViz`` CLI, install
the GUI extra:

.. code-block:: bash

   pip install 'moldenViz[gui]'

The GUI extra uses PySide6 as the supported Qt binding. The root package keeps
``Atom``, ``Parser``, ``Tabulator``, and the other model types as intentional
eager imports; importing them loads NumPy and Pydantic. ``Plotter`` remains a
lazy import, so core-only workflows do not load or require the GUI stack.

Prerequisites
-------------

The interactive plotting window relies on ``tkinter``. Verify that your Python install can import it:

.. code-block:: bash

   python3 -m tkinter

If that check fails, install ``tkinter`` manually:

* macOS:

  .. code-block:: bash

     brew install python-tk

* Ubuntu:

  .. code-block:: bash

     sudo apt-get install python3-tk

Quick CLI Preview
-----------------

Render a provided example molecule directly from the command line:

.. code-block:: bash

   moldenViz -e benzene

This will open an interactive 3D visualization window showing the benzene molecule with its molecular orbitals:

The viewer provides controls to:

* Navigate through different molecular orbitals
* Adjust contour levels and opacity
* Modify grid resolution
* Export data and images

Use your own Molden file instead:

.. code-block:: bash

   moldenViz my-calculation.molden

Get comfortable with the logging toggles early:

.. code-block:: bash

   moldenViz -v my-calculation.molden
   moldenViz -q my-calculation.molden

These flags help you audit which steps the CLI performs or keep output quiet during batch runs. See :doc:`Command Line Guide <cli-guide>` for more combinations.

Quick Python Preview
--------------------

Create plots programmatically by importing the high-level API:

.. code-block:: python

   from moldenViz import Plotter
   from moldenViz.examples import benzene

   Plotter(benzene)

This opens an interactive plotter window with full orbital visualization and control panel:

The Python API provides the same interactive capabilities as the CLI, with additional programmatic control over grid settings, tabulation, and export workflows.

Next Steps
----------

Once you can run the basics, move on to the topic-focused guides:

* :doc:`Command Line Guide <cli-guide>`
* :doc:`Python API Guide <python-api>`
* :doc:`Configuration Reference <configuration>`
* :doc:`Troubleshooting <troubleshooting>`
