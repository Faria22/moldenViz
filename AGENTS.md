# Repository Guidelines

- `docs/`: Sphinx project (`make html`) for authoring user docs; generated output lands in `docs/build/html/`.
- After every change, run `hatch run all`; if anything fails, fix the code before continuing.
- Run tests: `pytest` (append `--maxfail=1 -k name` while debugging). Use `pytest --cov=moldenViz` for coverage.
- Lint & format: `ruff check src tests` and `ruff format src tests`; static types via `basedpyright src tests`.
- Docs: `make -C docs html`; run `make -C docs clean` first when autosummary signatures change.
- Follow PEP 8 with 4-space indents; modules/functions `snake_case`, classes `PascalCase`.
- Docstrings follow NumPy style (`Parameters`, `Returns`, `Raises`); keep literal parameter names (e.g., `theta`, `phi`).
- Run `ruff format` before commits; it sorts imports and enforces consistent spacing.
- Confirm `pytest`, `ruff check`, `basedpyright`, and `make -C docs html` succeed locally before requesting review.
