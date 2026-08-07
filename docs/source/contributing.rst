Contributing
============

We welcome bug reports, feature ideas, documentation updates, and code contributions. Use the steps below to get set up and build confidence before opening a pull request.

Report Issues
-------------

- File issues on GitHub with a clear title, the version you are using, and (when relevant) a Molden file snippet that reproduces the problem.
- Tag documentation-only requests so they can be picked up quickly.

Set Up a Development Environment
--------------------------------

Install uv, just 1.54 or newer, and Python 3.10 or newer. Then:

.. code-block:: bash

   git clone https://github.com/Faria22/moldenViz.git
   cd moldenViz
   just sync

Run Tests and Linters
---------------------

Always run the test suite and static checks before pushing:

.. code-block:: bash

   just all

The combined ``all`` recipe runs linting, type checks, tests, and the
documentation build in sequence. To invoke an individual stage, use:

.. code-block:: bash

   just test
   just lint
   just typecheck
   just docs

Linting auto-fixes are available with ``just lint --fix``. All Python commands
run through uv with the committed lockfile.

Build Documentation
-------------------

Docs live under ``docs/`` and use Sphinx:

.. code-block:: bash

   just docs

Open ``docs/build/html/index.html`` in a browser to preview changes.

Plotter Architecture
--------------------

Keep changes to the interactive plotter within these boundaries:

- ``qt.py`` owns the embeddable ``OrbitalViewer``, Qt controls, queued result
  delivery, per-viewer state, and explicit render-window teardown. It must not
  create ``QApplication``, start an event loop, show itself, or open dialogs.
  Viewer state and operations must remain usable when the built-in control
  panel is hidden; the panel mirrors the viewer rather than owning selection.
- ``plotter.py`` owns the free-floating ``QMainWindow``, application/event-loop
  ownership, native file dialogs, message boxes, and menu shortcuts.
- ``_plotter_rendering.py`` owns PyVista scene and orbital rendering.
  Keep PyVista, PyVistaQt, Matplotlib, and visualization configuration imports
  behind viewer construction so importing ``Plotter`` does not initialize the
  rendering backend.
- ``_plotter_jobs.py`` owns background-job state and must remain independent
  of Qt and PyVista so it can be tested without a GUI.
- ``tabulator.py`` exposes the parsed data and computation operations needed
  by those layers; plotter modules must not access private ``Tabulator``
  fields.

The underscored modules are internal implementation details. Keep
``moldenViz.Plotter`` and ``moldenViz.plotter.Plotter`` as the supported
standalone entry points and ``moldenViz.qt.OrbitalViewer`` as the supported
embedding entry point. Background workers may compute numerical arrays, but
all Qt and VTK updates must be delivered on the owning GUI thread.

The package root uses module ``__getattr__`` to preserve its public re-exports
without importing NumPy or GUI dependencies for ``import moldenViz``. Add new
root exports to the lazy import map and retain them in ``__all__``.

Pull Request Checklist
----------------------

- Reference the issue number and describe behaviour changes.
- Mention new CLI flags or configuration keys in the docs and changelog.
- Attach screenshots for UI tweaks.
- Confirm ``just all`` succeeds.

Prepare a Release
-----------------

``pyproject.toml`` is the single source of truth for the package version. To
preview a patch release without changing files, commits, tags, or remotes:

.. code-block:: bash

   just release patch -d

Replace ``patch`` with another bump supported by ``uv version`` when needed.
The long form of the flag is ``--dry-run``. Running the recipe without the
dry-run flag requires a clean working tree. It
updates the version and lockfile, validates the project, builds the
distributions, commits the version bump, creates the corresponding ``v`` tag,
and atomically pushes the commit and tag to ``origin``.

Roadmap and Community
---------------------

See :doc:`Roadmap <roadmap>` for upcoming work. Share renderings, notebooks, or integration ideas in GitHub discussions—community examples help us improve defaults and documentation.
