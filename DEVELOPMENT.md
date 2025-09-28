# Development Guide

This document collects the essentials for contributing to `moldenViz`—from environment setup
through validation tasks—so you can get productive quickly and stay consistent with the rest of
the project.

## Environment Setup
- Ensure Python 3.8 or newer is available (`pyproject.toml` sets `requires-python = ">=3.8"`).
- Create and activate a virtual environment (`python -m venv .venv` followed by `source .venv/bin/activate` on Unix-like systems or `.venv\\Scripts\\activate` on Windows).
- Install development dependencies in editable mode: `pip install -e .[dev]`.
- If you plan to use the interactive plotter, confirm that `tkinter` is installed
  (`python -m tkinter` should open a test window).

## Project Layout
- `src/moldenViz/`: Core library modules (parser, tabulator, plotting widgets, configuration loader, CLI).
- `tests/`: Pytest suites that mirror the source modules (add new tests via `test_<module>.py`).
- `docs/`: Sphinx documentation project; generated HTML lives in `docs/build/html/`.
- `dist/` & `Library/`: Build artifacts—leave untouched unless you are packaging a release.

## Daily Development Tasks
- **Formatting & linting**: `ruff format src tests` then `ruff check src tests`.
- **Static typing**: `basedpyright src tests`.
- **Unit tests**: `pytest` (optionally `--maxfail=1 -k <pattern>` while iterating).
- **Coverage**: `pytest --cov=moldenViz --cov-report=term-missing` before opening a PR.
- **Docs build**: `make -C docs clean` then `make -C docs html` when signatures or API docs change.

## Coding Conventions
- Follow PEP 8 with 4-space indentation; prefer explicit imports (`from moldenViz.parser import Parser`).
- Name modules/functions with `snake_case`, classes with `PascalCase`, and keep shared plotting helpers in `_plotting_objects.py` to avoid circular imports.
- Write NumPy-style docstrings (`Parameters`, `Returns`, `Raises`) and use literal parameter names in documentation.
- Default to ASCII for new files; introduce non-ASCII only when already in use and justified.
- Keep code comments sparse and focused on intent or non-obvious context.

## Contribution Workflow
- Start feature work or fixes on a dedicated branch.
- Keep commits scoped and present tense (e.g., `Add grid validation`).
- Run `ruff format`, `ruff check`, `basedpyright`, `pytest`, and (if relevant) `make -C docs html` before pushing.
- When preparing a PR, note behaviour changes, link related issues, list test commands executed, and mention documentation updates if user-facing behaviour changed.

## Troubleshooting Tips
- Use `pytest --maxfail=1` to stop on the first failure while debugging.
- If PyVista-related tests are slow or flaky, isolate them with `-k` filters or mock heavy rendering paths.
- Rerun `make -C docs clean` before rebuilding docs when autosummary or API signatures change to avoid stale pages.
