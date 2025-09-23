# Repository Guidelines

## Project Structure & Module Organization
- `src/moldenViz/`: core package housing the parser, tabulator, plotting widgets, configuration loader, and CLI entry point.
- `tests/`: pytest suites that mirror the source modules (`test_parser.py`, etc.). Add new cases by matching the target module name.
- `docs/`: Sphinx project (`make html`) for authoring user docs; generated output lands in `docs/build/html/`.
- `dist/` and `Library/`: build artifacts—leave untouched unless packaging.

## Build, Test, and Development Commands
- Install dependencies: `pip install -e .[dev]` (adds pytest, ruff, basedpyright, mypy).
- Run tests: `pytest` (append `--maxfail=1 -k name` while debugging). Use `pytest --cov=moldenViz` for coverage.
- Lint & format: `ruff check src tests` and `ruff format src tests`; static types via `basedpyright src tests`.
- Docs: `make -C docs html`; run `make -C docs clean` first when autosummary signatures change.

## Coding Style & Naming Conventions
- Follow PEP 8 with 4-space indents; modules/functions `snake_case`, classes `PascalCase`.
- Docstrings follow NumPy style (`Parameters`, `Returns`, `Raises`); keep literal parameter names (e.g., `theta`, `phi`).
- Keep imports explicit and ordered (`from moldenViz.parser import Parser`). Avoid circular references by centralising shared helpers in `_plotting_objects.py`.
- Run `ruff format` before commits; it sorts imports and enforces consistent spacing.

## Testing Guidelines
- Place tests under `tests/` using the pattern `test_<feature>.py` with functions `test_<behavior>`.
- Focus on parser edge cases, grid/grid-type conversions, and CLI flows. Mock PyVista-heavy paths when feasible.
- Use `pytest --cov=moldenViz --cov-report=term-missing` before submission and ensure new features include regression tests.

## Commit & Pull Request Guidelines
- Commits stay small and descriptive in present tense (`Add grid validation`) as shown in `git log`.
- PRs should link issues, outline behaviour changes, enumerate tests run, and note doc updates if user-facing behaviour shifts.
- Confirm `pytest`, `ruff check`, `basedpyright`, and `make -C docs html` succeed locally before requesting review.
