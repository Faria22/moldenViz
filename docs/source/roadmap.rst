Roadmap
=======

This page summarises the major themes planned for upcoming releases. GitHub
milestones are the source of truth for issue-level scope and progress.

Next Release
------------

`v2.0 <https://github.com/Faria22/moldenViz/milestone/1>`_ brings together the
substantial improvements completed since v1.11:

- Faster, more memory-conscious GTO and MO tabulation through matrix
  multiplication, bounded work, primitive-exponential reuse, and direct
  solid-harmonic evaluation in Cartesian coordinates.
- Responsive and safer background grid and orbital computation, with stronger
  lifecycle and concurrency coverage.
- Clearer public APIs and more focused plotter implementation boundaries.
- A lightweight core installation and fewer required runtime dependencies.
- ``uv``-based packaging and development workflows, dedicated ASV benchmarks,
  expanded performance and peak-memory coverage, and Python 3.14 support.
- Recent parser, configuration, plotting, documentation, and CI fixes.

The remaining release work is tracked in `issue #118
<https://github.com/Faria22/moldenViz/issues/118>`_.

Future
------

`v3.0 <https://github.com/Faria22/moldenViz/milestone/2>`_ is planned to add
reading and tabulation of Molden files that contain Cartesian Gaussian basis
functions. This work is tracked in `issue #117
<https://github.com/Faria22/moldenViz/issues/117>`_.

This is separate from v2.0's solid-harmonic performance work: evaluating the
currently supported spherical basis on Cartesian coordinates does not add
support for Cartesian-basis coefficient layouts in Molden files.

Shipped
-------

- **v1.11** – Public format-specific tabulator exporters, lazy plotter imports, and hardened spherical-coordinate conversion.
- **v1.10** – Background orbital tabulation for a more responsive plotter, broader plotter regression coverage, and tabulation improvements.
- **v1.9** – CLI version and logging controls, coloured progress messages, and automated release workflows.
- **v1.8** – Image export (PNG, JPEG, SVG, PDF), enhanced GUI export dialogs, improved export workflows.
- **v1.1** – Volumetric export to VTK/Gaussian cube, expanded documentation, CLI option reference.
- **v1.0** – Configuration system overhaul, revised docs, PyVista plotter improvements.
- **v0.3** – CLI for interactive viewing of Molden files and bundled examples.
- **v0.2** – PyVista-based 3D plotter with Tkinter controls.
- **v0.1** – Initial parser and Python API.

How to Influence the Roadmap
----------------------------

- Review the `GitHub milestones
  <https://github.com/Faria22/moldenViz/milestones>`_ and open an issue
  describing your workflow.
- Vote on existing issues so we can prioritise the most impactful requests.
- Contribute prototype implementations or documentation to accelerate features.
