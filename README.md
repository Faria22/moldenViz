# moldenViz

[![PyPI - Version](https://img.shields.io/pypi/v/moldenviz.svg)](https://pypi.org/project/moldenviz)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/moldenviz.svg)](https://pypi.org/project/moldenviz)
[![Documentation Status](https://readthedocs.org/projects/moldenviz/badge/?version=latest)](https://moldenviz.readthedocs.io/en/latest/?badge=latest)

-----

## Installation

```console
pip install moldenViz
```

``moldenViz`` uses ``tkinter`` for its GUI. If ``python3 -m tkinter`` fails, install the tkinter package provided by your operating system (``brew install python-tk`` on macOS, ``sudo apt-get install python3-tk`` on Ubuntu).

## Quick start

- Launch the viewer with an example molecule:

  ```console
  moldenViz -e benzene
  ```

- Use the Python API for scripted workflows:

  ```python
  from moldenViz import Plotter
  Plotter('my.molden')
  ```

Full CLI usage, configuration examples, and API walkthroughs live in the docs.

## Documentation

Latest docs: https://moldenviz.readthedocs.io/en/latest/

## Roadmap

Major milestones and planned features are tracked in the [Roadmap](https://moldenviz.readthedocs.io/en/latest/roadmap.html). Highlights:

- ✅ v1.1 – VTK/cube export, expanded CLI reference, richer docs.
- ▶️ v2.0 – cartesian basis support

## Contributing

Guidelines for reporting issues, running tests, and building docs are in the [Contributing guide](https://moldenviz.readthedocs.io/en/latest/contributing.html).
