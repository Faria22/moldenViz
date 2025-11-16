# Code Health Findings

Notes gathered while reviewing the current implementation; use these as
starting points for future improvements.

## Coordinate Conversion Robustness
- File: `src/moldenViz/tabulator.py:75`
- Issue: `_cartesian_to_spherical` computes `theta = np.arccos(z / r)` without
  guarding the origin. When the grid includes `(0, 0, 0)` (default symmetric
  grids do) `r == 0`, producing `z / r` NaNs that propagate into the spherical
  harmonics tabulation.
- Potential fix: branch on `r == 0` and set `theta`/`phi` explicitly before
  applying `arccos`.
- Bonus: clip `z / r` to `[-1, 1]` to avoid NaNs caused by floating-point
  round-off in near-origin points.
