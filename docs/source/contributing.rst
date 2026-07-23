Contributing
============

We welcome bug reports, feature ideas, documentation updates, and code contributions. Use the steps below to get set up and build confidence before opening a pull request.

Report Issues
-------------

- File issues on GitHub with a clear title, the version you are using, and (when relevant) a Molden file snippet that reproduces the problem.
- Tag documentation-only requests so they can be picked up quickly.

Set Up a Development Environment
--------------------------------

.. code-block:: bash

   git clone https://github.com/Faria22/moldenViz.git
   cd moldenViz
   pip install -e '.[dev,gui]'

Run Tests and Linters
---------------------

Always run the test suite and static checks before pushing:

.. code-block:: bash

   hatch run all

The combined ``all`` script runs tests, lint, and type checks in sequence. To invoke an individual stage, use:

.. code-block:: bash

   hatch run test
   hatch run lint
   hatch run typecheck

Linting auto-fixes are available with ``hatch run lint --fix``.

Build Documentation
-------------------

Docs live under ``docs/`` and use Sphinx:

.. code-block:: bash

   make -C docs clean html

Open ``docs/build/html/index.html`` in a browser to preview changes.

Plotter Architecture
--------------------

Keep changes to the interactive plotter within these boundaries:

- ``plotter.py`` coordinates construction and connects the other layers.
- ``_plotter_ui.py`` owns Tk and Qt menus, widgets, and dialogs.
- ``_plotter_rendering.py`` owns PyVista scene and orbital rendering.
- ``_plotter_jobs.py`` owns background-job state and must remain independent
  of Tk, Qt, and PyVista so it can be tested without a GUI.
- ``tabulator.py`` exposes the parsed data and computation operations needed
  by those layers; plotter modules must not access private ``Tabulator``
  fields.

The underscored modules are internal implementation details. Keep
``moldenViz.Plotter`` and ``moldenViz.plotter.Plotter`` as the supported
public entry points, and pass dependencies through narrow methods rather than
introducing imports between the UI and rendering modules.

Pull Request Checklist
----------------------

- Reference the issue number and describe behaviour changes.
- Mention new CLI flags or configuration keys in the docs and changelog.
- Attach screenshots for UI tweaks.
- Confirm ``pytest``, ``ruff check``, ``basedpyright``, and ``make -C docs html`` all succeed.

Roadmap and Community
---------------------

See :doc:`Roadmap <roadmap>` for upcoming work. Share renderings, notebooks, or integration ideas in GitHub discussions—community examples help us improve defaults and documentation.
