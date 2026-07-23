API Stability
=============

``moldenViz`` distinguishes supported public API from implementation details.
This policy applies beginning with version 2.0.

Supported API
-------------

A Python name is supported when it is exported from the package root, listed
in a public module's ``__all__``, or explicitly included in the API reference.
The documented ``moldenViz`` command and its flags are also supported.

The individual datasets exported by ``moldenViz.examples`` are supported.
The internal example registry used by the command-line interface is not public.

The ``moldenViz.models`` namespace is reserved for parser and core result
types. Visualization configuration types are not exported from that module;
``AtomType`` is supported through ``moldenViz.AtomType`` only.

Private API
-----------

Names and modules beginning with an underscore are implementation details.
They may change or be removed without a deprecation period. Objects imported
inside a public module are not supported merely because Python permits direct
attribute access; consult that module's ``__all__`` instead.

Compatibility
-------------

Supported API changes follow semantic versioning. Breaking changes are reserved
for major releases. Deprecation periods may be provided for supported names in
minor releases, but private names do not receive compatibility guarantees.
