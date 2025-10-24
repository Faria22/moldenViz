Getting Started
===============

``moldenViz`` lets you inspect Molden outputs either from the command line or directly in Python. This page covers the essentials so you can install the package and render your first molecule quickly.

What is a Molden file?
----------------------

A Molden file stores the atomic structure and molecular orbital information produced by many quantum-chemistry codes (GAUSSIAN, ORCA, Molpro, and others). Each file contains ``[Atoms]`` coordinates, contracted Gaussian shells in ``[GTO]``, and orbital coefficients under ``[MO]``. ``moldenViz`` currently reads Molden files that use **spherical Gaussian functions**; Cartesian variants (sometimes labelled 5D/9G) raise an error during parsing. For the full format specification, refer to the `official Molden format description <https://www.theochem.ru.nl/molden/molden_format.html>`_.

Installation
------------

Install the package from PyPI:

.. code-block:: bash

   pip install moldenViz

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

.. image:: _static/placeholder-cli.png
   :alt: Screenshot of the moldenViz CLI window showing benzene with orbital controls
   :align: center
   :class: screenshot-placeholder

Use your own Molden file instead:

.. code-block:: bash

   moldenViz my-calculation.molden

Get comfortable with the logging toggles early:

.. code-block:: bash

   moldenViz --version
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

.. image:: _static/placeholder-python.png
   :alt: Screenshot of the moldenViz Python plotter rendering an isosurface next to molecule controls
   :align: center
   :class: screenshot-placeholder

Next Steps
----------

Once you can run the basics, move on to the topic-focused guides:

* :doc:`Command Line Guide <cli-guide>`
* :doc:`Python API Guide <python-api>`
* :doc:`Configuration Reference <configuration>`
* :doc:`Troubleshooting <troubleshooting>`
