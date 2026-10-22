Roadmap
=======

GitHub `milestones <https://github.com/Faria22/moldenViz/milestones>`_ are the
source of truth for issue-level scope and progress. The latest published
release is `v2.3.1
<https://github.com/Faria22/moldenViz/releases/tag/v2.3.1>`_.

Maintenance Plans
-----------------

The next two patch releases preserve the supported public API, CLI, and Python
support floor. No release dates are committed.

- `v2.3.2 <https://github.com/Faria22/moldenViz/milestone/5>`_ comes first:
  refresh the roadmap, fix PyPI metadata and README rendering, clean up CI and
  release-milestone automation, compress documentation screenshots, report test
  coverage in CI, and define a deliberate Ruff update policy.
- `v2.3.3 <https://github.com/Faria22/moldenViz/milestone/6>`_ follows that
  coverage and CI baseline: split the Qt viewer into focused modules and create
  the GTO thread pool lazily. These internal refactors must preserve supported
  imports, viewer lifecycle behaviour, tabulation results, and performance;
  they must not delay v2.3.2.

New public APIs and other feature work require separate planning for a future
minor release, outside these maintenance patches.

Shipped
-------

- `v2.3.1 <https://github.com/Faria22/moldenViz/releases/tag/v2.3.1>`_ – Viewer
  screenshots and control-panel improvements, including repainting native
  controls after display-scale changes.
- `v2.3.0 <https://github.com/Faria22/moldenViz/releases/tag/v2.3.0>`_ – Axes
  visibility toggle.
- `v2.2.0 <https://github.com/Faria22/moldenViz/releases/tag/v2.2.0>`_ – Atom-only
  molecule parsing.
- `v2.1.0 <https://github.com/Faria22/moldenViz/releases/tag/v2.1.0>`_ – Adaptive
  contour grids.
- `v2.0.0 <https://github.com/Faria22/moldenViz/releases/tag/v2.0.0>`_ – Embeddable
  PySide6 viewer replacing the Tk viewer; NumPy-only core installation with a
  ``[gui]`` extra; matrix-multiplication MO tabulation and chunked GTO tabulation;
  generated solid-harmonic fast paths through g shells; safer asynchronous grid
  and orbital computation; public API and semantic-versioning policy; ``uv``
  packaging and development workflows, ASV benchmarks, and Python 3.14 support.
- **v1.11** – Public format-specific tabulator exporters, lazy plotter imports, and hardened spherical-coordinate conversion.
- **v1.10** – Background orbital tabulation for a more responsive plotter, broader plotter regression coverage, and tabulation improvements.
- **v1.9** – CLI version and logging controls, coloured progress messages, and automated release workflows.
- **v1.8** – Image export (PNG, JPEG, SVG, PDF), enhanced GUI export dialogs, improved export workflows.
- **v1.1** – Volumetric export to VTK/Gaussian cube, expanded documentation, CLI option reference.
- **v1.0** – Configuration system overhaul, revised docs, PyVista plotter improvements.
- **v0.3** – CLI for interactive viewing of Molden files and bundled examples.
- **v0.2** – PyVista-based 3D plotter with Tkinter controls.
- **v0.1** – Initial parser and Python API.

Not Planned
-----------

The following proposals were closed as not planned and are outside the current
roadmap:

- Cartesian-basis Molden reading and tabulation (`#117
  <https://github.com/Faria22/moldenViz/issues/117>`_).
- h and higher angular-momentum shells (`#119
  <https://github.com/Faria22/moldenViz/issues/119>`_).
- Rust acceleration for orbital tabulation (`#104
  <https://github.com/Faria22/moldenViz/issues/104>`_).
- Periodic systems and CP2K crystal-orbital support (`#67
  <https://github.com/Faria22/moldenViz/issues/67>`_).

The shipped solid-harmonic optimizations evaluate the supported spherical basis
on Cartesian coordinates; they do not add Cartesian-basis coefficient layouts
or support for shells above g.

How to Influence the Roadmap
----------------------------

- Review the `GitHub milestones
  <https://github.com/Faria22/moldenViz/milestones>`_ and open an issue
  describing your workflow.
- Vote on existing issues so we can prioritise the most impactful requests.
- Contribute prototype implementations or documentation to accelerate features.
