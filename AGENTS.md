# Repository Guidelines

- `docs/`: Sphinx project (`make html`) for authoring user docs; generated output lands in `docs/build/html/`.
- After every change, run `just all`; if anything fails, fix the code before continuing.
- Run tests: `just test` (append `--maxfail=1 -k name` while debugging). Use `just cov` for coverage.
- Lint & format: `just lint` and `just format`; static types via `just typecheck`.
- Docs: `just docs`; it cleans the Sphinx output before rebuilding.
- Follow PEP 8 with 4-space indents; modules/functions `snake_case`, classes `PascalCase`.
- Docstrings follow NumPy style (`Parameters`, `Returns`, `Raises`); keep literal parameter names (e.g., `theta`, `phi`).
- Run `just format` before commits; Ruff sorts imports and enforces consistent spacing.
- Confirm `just all` succeeds locally before requesting review.
