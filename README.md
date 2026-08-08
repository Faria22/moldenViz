# moldenViz

[![PyPI - Version](https://img.shields.io/pypi/v/moldenviz.svg)](https://pypi.org/project/moldenviz)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/moldenviz.svg)](https://pypi.org/project/moldenviz)
[![Documentation Status](https://readthedocs.org/projects/moldenviz/badge/?version=latest)](https://moldenviz.readthedocs.io/en/latest/?badge=latest)

-----

## Installation

Install the core parser and tabulator:

```console
pip install moldenViz
```

Or add it to a project with uv:

```console
uv add moldenViz
```

Install the interactive viewer and CLI dependencies with the GUI extra:

```console
pip install 'moldenViz[gui]'
```

With uv:

```console
uv add 'moldenViz[gui]'
```

Or install the CLI as a standalone uv tool:

```console
uv tool install 'moldenViz[gui]'
```

This makes the ``moldenViz`` command available on your ``PATH``.

The GUI extra uses PySide6 as its only GUI toolkit. The standalone CLI creates
and runs the Qt application for you; applications that embed the viewer keep
ownership of their existing Qt event loop.

## Quick start

- Launch the viewer with an example molecule:

  ```console
  moldenViz -e benzene
  ```

- Review the [CLI guide](docs/source/cli-guide.rst) for version checks, verbosity toggles, and other flags you can pass to
  ``moldenViz``.

- Use the Python API for scripted workflows:

  ```python
  from moldenViz import Parser, Tabulator

  parser = Parser(filename='my.molden')
  tabulator = Tabulator(filename='my.molden')
  ```

- With the GUI extra installed, launch a viewer from Python:

  ```python
  from moldenViz import Plotter

  Plotter(filename='my.molden')
  ```

- Embed the same controls and renderer in an existing PySide6 application:

  ```python
  from moldenViz.qt import OrbitalViewer

  viewer = OrbitalViewer(parent=page)
  page.layout().addWidget(viewer)
  viewer.set_input(filename='my.molden')
  ```

  A host-owned dashboard can hide the built-in controls and drive the public
  viewer API directly:

  ```python
  viewer = OrbitalViewer(filename='my.molden', parent=page, show_controls=False)
  viewer.show_orbital(0)
  viewer.update_appearance(contour=0.05, mo_opacity=0.8)
  viewer.set_controls_visible(True)  # Restore the built-in panel when needed.
  ```

Full CLI usage, configuration examples, and API walkthroughs live in the docs.

## Documentation

Latest docs: https://moldenviz.readthedocs.io/en/latest/

## Roadmap

GitHub milestones are the source of truth for release scope and progress:

- [v2.0](https://github.com/Faria22/moldenViz/milestone/1) is the next release, bringing together the recent performance, responsiveness, reliability, API, packaging, and developer-tooling improvements.
- [v3.0](https://github.com/Faria22/moldenViz/milestone/2) is planned to add reading and tabulation of Cartesian-basis Molden files.

See the [documentation roadmap](https://moldenviz.readthedocs.io/en/latest/roadmap.html) for a high-level summary.

## Contributing

Guidelines for reporting issues, running tests, and building docs are in the [Contributing guide](https://moldenviz.readthedocs.io/en/latest/contributing.html).
