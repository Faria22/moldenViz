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

Install the interactive viewer and CLI dependencies with the GUI extra:

```console
pip install 'moldenViz[gui]'
```

The GUI extra uses PySide6 as its supported Qt binding. ``moldenViz`` also uses
``tkinter``, which Python distributions commonly provide separately. If
``python3 -m tkinter`` fails, install the tkinter package provided by your
operating system (``brew install python-tk`` on macOS,
``sudo apt-get install python3-tk`` on Ubuntu).

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

  parser = Parser('my.molden')
  tabulator = Tabulator(parser)
  ```

- With the GUI extra installed, launch a viewer from Python:

  ```python
  from moldenViz import Plotter

  Plotter('my.molden')
  ```

Full CLI usage, configuration examples, and API walkthroughs live in the docs.

## Documentation

Latest docs: https://moldenviz.readthedocs.io/en/latest/

## Roadmap

Major milestones and planned features are tracked in the [Roadmap](https://moldenviz.readthedocs.io/en/latest/roadmap.html). Highlights:

- ✅ v1.11 – Public tabulator exporters, lazy plotter imports, and numerical hardening.
- ✅ v1.10 – Responsive background orbital tabulation and expanded plotter coverage.
- ✅ v1.9 – CLI version and logging controls with coloured progress messages.
- ✅ v1.8 – Image export (PNG, JPEG, SVG, PDF) and enhanced GUI export dialogs.
- ✅ v1.1 – VTK/cube export, expanded CLI reference, richer docs.
- ▶️ v2.0 – Cartesian basis support

## Contributing

Guidelines for reporting issues, running tests, and building docs are in the [Contributing guide](https://moldenviz.readthedocs.io/en/latest/contributing.html).
