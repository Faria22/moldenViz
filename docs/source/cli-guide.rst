Command Line
============

The ``moldenViz`` CLI offers a fast way to inspect Molden files and bundled examples. Run ``moldenViz -h`` for a full option list, and use the sections below for the most common workflows.

Basic Usage
-----------

Render molecular orbitals from a Molden file:

.. code-block:: bash

   moldenViz path/to/file.molden

Show only the molecular structure (skip orbital surfaces):

.. code-block:: bash

   moldenViz path/to/file.molden -m

List the available commands and flags:

.. code-block:: bash

   moldenViz -h

Bundled Examples
----------------

Use one of the shipped example molecules when you just want to explore the plotting experience:

.. code-block:: bash

   moldenViz -e co

Available examples:

- ``co``
- ``o2``
- ``co2``
- ``h2o``
- ``benzene``
- ``prismane``
- ``pyridine``
- ``furan``
- ``acrolein``

Tips
----

- Combine ``-e`` with ``-m`` to explore only the geometry of an example.
- Run with ``-v`` to increase logging output when debugging an issue.
- Use configuration overrides (:doc:`see details <configuration>`) to change colours, bond lengths, or orbital contours before launching the CLI.
