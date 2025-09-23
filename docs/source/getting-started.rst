Getting Started
===============

``moldenViz`` lets you inspect Molden outputs either from the command line or directly in Python. This page covers the essentials so you can install the package and render your first molecule quickly.

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

Use your own Molden file instead:

.. code-block:: bash

   moldenViz my-calculation.molden

Quick Python Preview
--------------------

Create plots programmatically by importing the high-level API:

.. code-block:: python

   from moldenViz import Plotter
   from moldenViz.examples import benzene
   
   Plotter(benzene)

Next Steps
----------

Once you can run the basics, move on to the topic-focused guides:

* :doc:`Command Line Guide <cli-guide>`
* :doc:`Python API Guide <python-api>`
* :doc:`Configuration Reference <configuration>`
* :doc:`Troubleshooting <troubleshooting>`
