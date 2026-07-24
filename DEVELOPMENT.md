# Development Guide

This document collects the essentials for contributing to `moldenViz`—from environment setup
through validation tasks—so you can get productive quickly and stay consistent with the rest of
the project.

## Environment Setup
- Install uv, just 1.54 or newer, and Python 3.10 or newer
  (`pyproject.toml` sets `requires-python = ">=3.10"`).
- Create the locked development environment with `just sync`.
- If you plan to use the interactive plotter, confirm that `tkinter` is installed
  (`uv run python -m tkinter` should open a test window).

## Project Layout
- `src/moldenViz/`: Core library modules (parser, tabulator, plotting widgets, configuration loader, CLI).
- `tests/`: Pytest suites that mirror the source modules (add new tests via `test_<module>.py`).
- `benchmarks/`: ASV performance and peak-memory suites.
- `docs/`: Sphinx documentation project; generated HTML lives in `docs/build/html/`.
- `dist/` & `Library/`: Build artifacts—leave untouched unless you are packaging a release.

## Daily Development Tasks
- **Formatting & linting**: `just format` then `just lint`.
- **Static typing**: `just typecheck`.
- **Unit tests**: `just test` (optionally append `--maxfail=1 -k <pattern>` while iterating).
- **Coverage**: `just cov` before opening a PR.
- **Benchmark smoke test**: `just bench-smoke`.
- **Current-commit benchmarks**: `just bench`.
- **Docs build**: `just docs`.
- **All required checks**: `just all`.
- **Build wheel and source distributions**: `just build`.

See the [benchmark guide](https://moldenviz.readthedocs.io/en/latest/benchmarks.html)
for the full matrix, revision comparisons, result storage, and CI policy.

## Generated Numerical Kernels

The direct NumPy real solid-harmonic kernels through `l=4` are generated from
the finite polynomial definition with development-only SymPy. Regenerate the
committed runtime module from the repository root:

```console
uv run --with sympy==1.14.0 python scripts/generate_solid_harmonics.py
```

Confirm that a checkout is reproducible without rewriting the module:

```console
uv run --with sympy==1.14.0 python scripts/generate_solid_harmonics.py --check
```

The runtime package continues to depend only on NumPy. The general
finite-polynomial implementation remains in `Tabulator` as the correctness
oracle and fallback above `l=4`.

## Coding Conventions
- Follow PEP 8 with 4-space indentation; prefer explicit imports (`from moldenViz.parser import Parser`).
- Name modules/functions with `snake_case`, classes with `PascalCase`, and keep shared plotting helpers in `_plotting_objects.py` to avoid circular imports.
- Write NumPy-style docstrings (`Parameters`, `Returns`, `Raises`) and use literal parameter names in documentation.
- Default to ASCII for new files; introduce non-ASCII only when already in use and justified.
- Keep code comments sparse and focused on intent or non-obvious context.

## Contribution Workflow
- Start feature work or fixes on a dedicated branch.
- Keep commits scoped and present tense (e.g., `Add grid validation`).
- Run `just all` before pushing.
- When preparing a PR, note behaviour changes, link related issues, list test commands executed, and mention documentation updates if user-facing behaviour changed.

## Troubleshooting Tips
- Use `just test --maxfail=1` to stop on the first failure while debugging.
- If PyVista-related tests are slow or flaky, isolate them with `-k` filters or mock heavy rendering paths.
- Rerun `just docs` when autosummary or API signatures change to avoid stale pages.

## Release Workflow

The version in `pyproject.toml` is the single source of truth. Preview a patch
release without changing files or Git state:

```console
just release patch -d
```

Replace `patch` with any bump supported by `uv version`, such as `minor` or
`major`. The long form of the flag is `--dry-run`. Without the dry-run flag,
the recipe bumps the version with uv, updates
`uv.lock`, runs all checks and builds the distributions, commits the two version
files, creates a `v<version>` tag, and atomically pushes the commit and tag to
`origin`. Start from a clean working tree.
