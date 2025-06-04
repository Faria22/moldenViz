Usage
=====

This page provides basic usage examples for the ``moldenViz`` library.

Prerequisites
-------------

Ensure you have a Molden file available (e.g., ``molden.inp``).

Basic Parsing
-------------

To get started, you can parse a Molden file using the ``Parser`` class. This will extract atomic information, Gaussian Type Orbitals (GTOs), and Molecular Orbitals (MOs).

.. code-block:: python

   from moldenViz import Parser

   # Replace 'path/to/your/sample_molden.inp' with the actual file path
   molden_file_path = 'path/to/your/sample_molden.inp'

   try:
       parser = Parser(filename=molden_file_path)

       # Access parsed data
       print(f"Parsed {len(parser.atoms)} atoms.")
       for atom in parser.atoms:
           print(f"Atom: {atom.label}, Position: {atom.position}, Num GTOs: {len(atom.gtos)}") # [cite: faria22/moldenviz/
