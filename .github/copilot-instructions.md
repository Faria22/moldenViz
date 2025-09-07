# moldenViz
moldenViz is a Python package for visualizing and analyzing Molden files containing molecular orbital data. It provides CLI tools and Python APIs for parsing Molden files, tabulating Gaussian-type orbitals (GTOs), and creating 3D molecular visualizations.

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Working Effectively
- Bootstrap, build, and test the repository:
  - `sudo apt-get update && sudo apt-get install -y python3-tk` -- installs tkinter (required for plotting functionality)
  - `python3 -m pip install -e .` -- takes 1-2 minutes to install all dependencies. NEVER CANCEL. Set timeout to 180+ seconds.
  - `python3 -m pytest tests/ -v` -- runs 43 tests in ~2 seconds. NEVER CANCEL. Set timeout to 30+ seconds.
- Linting and formatting:
  - `python3 -m pip install --user ruff`
  - `ruff check .` -- lint checking, takes <1 second
  - `ruff format .` -- code formatting, takes <1 second  
  - Always run `ruff check .` and `ruff format .` before committing changes
- Documentation:
  - `python3 -m pip install --user -r docs/requirements.txt` -- takes 1-2 minutes. NEVER CANCEL. Set timeout to 180+ seconds.
  - `cd docs && make html` -- builds documentation in ~4 seconds. NEVER CANCEL. Set timeout to 60+ seconds.

## Validation
- ALWAYS manually validate any new code changes via Python API testing when modifying core functionality.
- ALWAYS run through at least one complete end-to-end scenario after making changes:
  ```python
  from moldenViz import Parser, Tabulator
  from moldenViz.examples import co
  import numpy as np
  
  # Test parsing
  parser = Parser(co)
  print(f"Parsed {len(parser.atoms)} atoms, {len(parser.mos)} MOs")
  
  # Test tabulation
  tab = Tabulator(co)
  tab.cartesian_grid(x=np.linspace(-2,2,5), y=np.linspace(-2,2,5), z=np.linspace(-2,2,5))
  mo_data = tab.tabulate_mos(0)
  print(f"MO tabulation successful: {mo_data.shape}")
  ```
- You can build and test the core functionality (Parser, Tabulator), however GUI plotting requires special setup in headless environments.
- CLI functionality works for help and validation: `moldenViz --help`
- Always run `pytest tests/` and ensure all 43 tests pass before considering changes complete.

## GUI and Plotting Limitations
- **IMPORTANT**: GUI plotting (`moldenViz -e co` or `Plotter()` class) requires a graphical environment or virtual display setup.
- In headless environments, GUI plotting will fail with Qt/EGL errors - this is expected.
- For testing changes to plotting code, use virtual display: `xvfb-run -a moldenViz -e co -m`
- Core functionality (Parser, Tabulator) works perfectly in all environments.

## Common Tasks
The following are outputs from frequently run commands. Reference them instead of viewing, searching, or running bash commands to save time.

### Repository structure
```
ls -la /home/runner/work/moldenViz/moldenViz/
.git/
.gitignore
.readthedocs.yaml       # Read the Docs configuration
.ruff.toml             # Ruff linting configuration
LICENSE
README.md
docs/                  # Sphinx documentation
pyproject.toml         # Main project configuration
src/moldenViz/         # Main package source
tests/                 # Test suite
```

### Key files
- `pyproject.toml` - Main project configuration, dependencies, build system (hatchling)
- `src/moldenViz/__init__.py` - Package entry point, exports Parser, Plotter, Tabulator
- `src/moldenViz/_cli.py` - CLI entry point for `moldenViz` command
- `src/moldenViz/parser.py` - Core Molden file parsing functionality
- `src/moldenViz/plotter.py` - 3D visualization using PyVista/PySide6
- `src/moldenViz/tabulator.py` - GTO tabulation and grid generation
- `src/moldenViz/examples/` - Built-in example molecules (co, h2o, benzene, etc.)
- `.ruff.toml` - Linting configuration (strict settings, uses numpy docstring style)

### Package configuration (pyproject.toml)
```toml
[project]
name = "moldenViz"
requires-python = ">=3.8"
dependencies = ["numpy", "pyvista", "pyvistaqt", "scipy", "PySide6", "toml"]

[project.scripts]
moldenViz = "moldenViz._cli:main"
```

### CLI help output
```
$ moldenViz --help
usage: moldenViz [-h] [-m] [-e molecule] [file]

positional arguments:
  file                  Optional molden file path

options:
  -h, --help            show this help message and exit
  -m, --only_molecule   Only plots the molecule
  -e molecule, --example molecule
                        Load example molecule. Options are: co, o2, co2, h2o, benzene, prismane, pyridine, furan, acrolein
```

### Test execution
```
$ pytest tests/ -v
================================================= test session starts ==================================================
platform linux -- Python 3.12.3, pytest-7.4.4, pluggy-1.4.0 -- /usr/bin/python3
cachedir: .pytest_cache
rootdir: /home/runner/work/moldenViz/moldenViz
collected 43 items

tests/test_parser.py::test_section_indices_order PASSED                                                          [  2%]
tests/test_parser.py::test_gaussian_normalization_positive PASSED                                                [  4%]
tests/test_parser.py::test_atomic_orbital_permutation PASSED                                                     [  6%]
tests/test_parser.py::test_atom_labels PASSED                                                                    [  9%]
tests/test_parser.py::test_basis_and_mo_dimensions PASSED                                                        [ 11%]
tests/test_parser.py::test_mo_energies_are_sorted PASSED                                                         [ 13%]
tests/test_parser.py::test_parser_invalid_input_type[None] PASSED                                                [ 16%]
[... continues for all 43 tests ...]
================================================== 43 passed in 1.69s ==================================================
```

### Ruff linting
- Running `ruff check .` typically shows minor warnings about missing docstrings and test files accessing private members
- Running `ruff format .` should show "X files already formatted" if code is properly formatted
- The project uses strict linting with numpy docstring conventions

## Development Tips
- Always work in development mode: `pip install -e .`
- Example molecules are available: co, o2, co2, h2o, benzene, prismane, pyridine, furan, acrolein
- Core classes: `Parser` (file parsing), `Tabulator` (GTO calculations), `Plotter` (visualization)
- When making changes to parsing: always test with `Parser(moldenViz.examples.co)`
- When making changes to tabulation: always test with `Tabulator(moldenViz.examples.co)` and grid setup
- When making changes to plotting: test core functionality first, GUI plotting second
- The package supports both file paths and example names for input sources
- Grid types supported: cartesian_grid() and spherical_grid()
- MO tabulation supports single indices, lists, ranges, or None (all MOs)

## Expected Timing
- Package installation: 1-2 minutes
- Test suite: 2 seconds for all 43 tests
- Linting: <1 second
- Documentation build: 4 seconds
- Core API validation: <1 second per operation

## Architecture
- **Parser**: Reads Molden files, extracts atoms, basis functions, molecular orbitals
- **Tabulator**: Creates 3D grids, tabulates GTOs and MOs at grid points
- **Plotter**: Visualizes molecules and molecular orbitals using PyVista/Qt
- **CLI**: Provides command-line interface with example molecules and options
- Built with modern Python packaging (pyproject.toml, hatchling)
- Uses scientific Python stack (NumPy, SciPy) for computational core
- Uses PyVista/PySide6 for 3D visualization