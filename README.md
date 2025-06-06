# moldenViz

[![PyPI - Version](https://img.shields.io/pypi/v/moldenviz.svg)](https://pypi.org/project/moldenviz)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/moldenviz.svg)](https://pypi.org/project/moldenviz)
[![Documentation Status](https://readthedocs.org/projects/moldenviz/badge/?version=latest)](https://moldenviz.readthedocs.io/en/latest/?badge=latest)

-----

## Installation

```console
pip install moldenViz
```

**Note:** If you want to use the plotter, make sure python has access to `tkinter`
```console
python3 -m tkinter
```

If python doesn't have access to `tkinter`, then you can install it with
### macOS
```console
brew install python-tk
```

### Ubuntu
```console
sudo apt-get install python-tk
```

## Usage
If you don't have a molden file (then why are you using this package?) you can get some examples by doing
```python
from moldenViz.examples import co

Plotter(co)
```
the total list of examples is:
- co
- o2
- co2
- h2o
- benzene
- prismane
- pyridine
- furan
- acrolein

In the next examples, I'll be using `'molden.inp'`, but you can replace it with you molden-file path, or one of the examples.

### Plotter
You can use the plotter to plot just the molecule
```python
from moldenViz import Plotter

Plotter('molden.inp', only_molecule=True)
```
or to plot the molecular orbitals
```python
from moldenViz import Plotter

Plotter('molden.inp')
```

### Tabulator
You can use `moldenViz` to tabulate the GTOs and molecular orbitals
```python
from moldenViz import Tabulator
import numpy as np

tab = Tabulator('molden.inp')

# Using a spherical grid
tab.spherical_grid(
    r = np.linspace(0, 5, 20),
    theta = np.linspace(0, np.pi, 20)
    phi = np.linspace(0, 2 * np.pi, 40)
)

# Or a cartesian grid
tab.cartesian_grid(
    x = np.linspace(-2, 2, 20)
    y = np.linspace(-2, 2, 20)
    z = np.linspace(-2, 2, 20)
)

print(tab.grid.shape)
print(tab.gtos_data.shape)
```

And to tabulate a molecular orbital
```python
mo_data = tab.tabulate_mos(0)
```
or a list
```python
mo_data = tab.tabulate_mos([0,1,4])
```
or a range
```python
mo_data = tab.tabulate_mos(range(1, 10, 2))
```
or all of them
```python
mos_data = tab.tabulate_mos()
```

## Documentation
You can find the documentation [here](https://moldenviz.readthedocs.io/en/latest/).

## Roadmap
- [x] 0.1: Basic commands inside python.
  - [x] 0.1.5: Documentation support
- [x] 0.2: Plotter.
- [ ] 0.3: CLI options.
