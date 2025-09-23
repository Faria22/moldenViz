moldenViz Documentation
=======================

.. image:: https://img.shields.io/pypi/v/moldenviz.svg
   :target: https://pypi.org/project/moldenviz
   :alt: PyPI Version

.. image:: https://img.shields.io/pypi/pyversions/moldenviz.svg
   :target: https://pypi.org/project/moldenviz
   :alt: Python Versions

.. raw:: html

   <div style="margin-top: 20px;"></div>

``moldenViz`` is a Python package for parsing Molden files, tabulating molecular orbitals, and visualizing them through intuitive 3D plots. It provides both a command-line interface and a Python API for working with quantum chemistry calculation results.

Key Features
------------

* **Molden File Parsing**: Read and parse molden files to extract molecular structures and orbital data
* **3D Visualization**: Interactive plotting of molecules and molecular orbitals using PyVista
* **Grid Tabulation**: Create custom grids and tabulate Gaussian-type orbitals (GTOs)
* **CLI and Python API**: Use from command line or integrate into Python scripts
* **Example Molecules**: Built-in example molecules for testing and demonstration

Quick Start
-----------

Install moldenViz:

.. code-block:: bash

   pip install moldenViz

Plot a molecule from the command line:

.. code-block:: bash

   moldenViz -e benzene

Or use the Python API:

.. code-block:: python

   from moldenViz import Plotter
   from moldenViz.examples import benzene

   Plotter(benzene)

For a guided walk-through, continue with :doc:`Getting Started <getting-started>`.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   getting-started
   cli-guide
   python-api
   configuration
   troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: API Documentation

   api

.. toctree::
   :maxdepth: 1
   :caption: Additional Information

   GitHub Repository <https://github.com/Faria22/moldenViz>
   PyPI Package <https://pypi.org/project/moldenviz>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
