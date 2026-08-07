"""Adaptive Cartesian grid construction for molecular-orbital contours."""

from __future__ import annotations

import ast
from typing import TYPE_CHECKING

import numpy as np
import pyvista as pv

if TYPE_CHECKING:
    from collections.abc import Iterable

    from numpy.typing import NDArray

AdaptiveScale = float | tuple[float, float, float]
_AXIS_COUNT = 3


def parse_scale(value: str) -> AdaptiveScale:
    """Parse a scalar or three-axis adaptive scale from settings text.

    Returns
    -------
    float or tuple of float
        Parsed and validated scale.
    """
    try:
        parsed = ast.literal_eval(value.strip())
    except (SyntaxError, ValueError) as exc:
        raise ValueError('Adaptive scale must be a number or three-number tuple.') from exc
    if isinstance(parsed, int | float) and not isinstance(parsed, bool):
        result: AdaptiveScale = float(parsed)
    elif isinstance(parsed, list | tuple) and len(parsed) == _AXIS_COUNT:
        result = (float(parsed[0]), float(parsed[1]), float(parsed[2]))
    else:
        raise ValueError('Adaptive scale must be a number or three-number tuple.')
    normalize_scale(result)
    return result


def format_scale(scale: AdaptiveScale) -> str:
    """Format an adaptive scale for the settings entry.

    Returns
    -------
    str
        Scalar or tuple representation suitable for parsing again.
    """
    if isinstance(scale, int | float):
        return str(float(scale))
    return f'({scale[0]}, {scale[1]}, {scale[2]})'


def normalize_scale(scale: AdaptiveScale) -> tuple[float, float, float]:
    """Return a validated three-axis adaptive-grid scale.

    Parameters
    ----------
    scale : float or tuple of float
        Fine-grid resolution multiplier. Fine spacing is approximately the
        coarse spacing divided by this value.

    Returns
    -------
    tuple of float
        Scale for the x, y, and z axes.

    Raises
    ------
    ValueError
        If a scale is non-finite or less than one.
    """
    values = (scale, scale, scale) if isinstance(scale, int | float) else scale
    normalized = tuple(float(value) for value in values)
    if len(normalized) != _AXIS_COUNT or any(not np.isfinite(value) or value < 1.0 for value in normalized):
        raise ValueError('Adaptive scale must contain one or three finite values greater than or equal to 1.')
    return normalized


def crossed_cell_ids(
    coarse_grid: pv.StructuredGrid,
    mo_values: NDArray[np.floating],
    contour: float,
) -> NDArray[np.int64]:
    """Find coarse cells crossed by any available molecular orbital.

    PyVista first generates both signed contour surfaces for every MO. The
    centers of the generated contour faces are then mapped back to cells in the
    source structured grid.

    Parameters
    ----------
    coarse_grid : pyvista.StructuredGrid
        Coarse Cartesian grid.
    mo_values : NDArray[np.floating]
        Molecular-orbital values shaped ``(n_points, n_mos)``.
    contour : float
        Positive contour magnitude.

    Returns
    -------
    NDArray[np.int64]
        Sorted unique source-cell identifiers.
    """
    if contour <= 0:
        raise ValueError('Contour must be greater than zero.')
    values = np.asarray(mo_values)
    if values.ndim == 1:
        values = values[:, None]
    if values.shape[0] != coarse_grid.n_points:
        raise ValueError('MO values must have one row for every coarse-grid point.')

    crossed: set[int] = set()
    work_grid = coarse_grid.copy(deep=True)
    for mo_index in range(values.shape[1]):
        work_grid.point_data['orbital'] = values[:, mo_index]
        surface = work_grid.contour([-contour, contour], scalars='orbital')
        if surface.n_cells == 0:
            continue
        containing = work_grid.find_containing_cell(np.asarray(surface.cell_centers().points))
        crossed.update(int(cell_id) for cell_id in np.asarray(containing) if cell_id >= 0)
    return np.asarray(sorted(crossed), dtype=np.int64)


def _axis_fractions(scale: float) -> NDArray[np.float64]:
    fractions = np.arange(0.0, 1.0, 1.0 / scale, dtype=float)
    if fractions.size == 0 or not np.isclose(fractions[0], 0.0):
        fractions = np.insert(fractions, 0, 0.0)
    if not np.isclose(fractions[-1], 1.0):
        fractions = np.append(fractions, 1.0)
    else:
        fractions[-1] = 1.0
    return fractions


def refined_grid(
    coarse_grid: pv.StructuredGrid,
    cell_ids: Iterable[int],
    scale: AdaptiveScale,
) -> pv.UnstructuredGrid:
    """Subdivide selected coarse cells into a conforming hexahedral grid.

    Parameters
    ----------
    coarse_grid : pyvista.StructuredGrid
        Coarse Cartesian source grid.
    cell_ids : iterable of int
        Source cells to refine.
    scale : float or tuple of float
        Resolution multiplier for one or all three axes.

    Returns
    -------
    pyvista.UnstructuredGrid
        Refined cells with shared boundary points deduplicated.
    """
    scales = normalize_scale(scale)
    fractions = tuple(_axis_fractions(value) for value in scales)
    point_ids: dict[tuple[float, float, float], int] = {}
    points: list[tuple[float, float, float]] = []
    hexes: list[list[int]] = []

    def point_id(point: tuple[float, float, float]) -> int:
        existing = point_ids.get(point)
        if existing is not None:
            return existing
        new_id = len(points)
        point_ids[point] = new_id
        points.append(point)
        return new_id

    for cell_id in sorted({int(value) for value in cell_ids}):
        bounds = coarse_grid.get_cell(cell_id).bounds
        axes = tuple(
            low + (high - low) * axis_fractions
            for low, high, axis_fractions in zip(
                (bounds.x_min, bounds.y_min, bounds.z_min),
                (bounds.x_max, bounds.y_max, bounds.z_max),
                fractions,
                strict=True,
            )
        )
        for ix in range(len(axes[0]) - 1):
            for iy in range(len(axes[1]) - 1):
                for iz in range(len(axes[2]) - 1):
                    corners = (
                        (axes[0][ix], axes[1][iy], axes[2][iz]),
                        (axes[0][ix + 1], axes[1][iy], axes[2][iz]),
                        (axes[0][ix + 1], axes[1][iy + 1], axes[2][iz]),
                        (axes[0][ix], axes[1][iy + 1], axes[2][iz]),
                        (axes[0][ix], axes[1][iy], axes[2][iz + 1]),
                        (axes[0][ix + 1], axes[1][iy], axes[2][iz + 1]),
                        (axes[0][ix + 1], axes[1][iy + 1], axes[2][iz + 1]),
                        (axes[0][ix], axes[1][iy + 1], axes[2][iz + 1]),
                    )
                    hexes.append(
                        [point_id((float(corner[0]), float(corner[1]), float(corner[2]))) for corner in corners],
                    )

    if not hexes:
        return pv.UnstructuredGrid()
    cells = np.column_stack((np.full(len(hexes), 8, dtype=np.int64), np.asarray(hexes, dtype=np.int64))).ravel()
    cell_types = np.full(len(hexes), pv.CellType.HEXAHEDRON, dtype=np.uint8)
    return pv.UnstructuredGrid(cells, cell_types, np.asarray(points))
