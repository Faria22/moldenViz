"""Adaptive Cartesian grid construction for molecular-orbital contours."""

from __future__ import annotations

import ast
from importlib import import_module
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterable

    import pyvista as pv
    from numpy.typing import NDArray

AdaptiveScale = float | tuple[float, float, float]
_AXIS_COUNT = 3
_HEXAHEDRON_POINT_COUNT = 8
_VOXEL_TO_HEXAHEDRON = (0, 1, 3, 2, 4, 5, 7, 6)


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

    Raises
    ------
    ValueError
        If the source is not an axis-aligned Cartesian grid or a cell ID is
        outside the source grid.
    """
    pv = import_module('pyvista')
    scales = normalize_scale(scale)
    selected_cell_ids = np.asarray(sorted({int(value) for value in cell_ids}), dtype=np.int64)
    if selected_cell_ids.size == 0:
        return pv.UnstructuredGrid()
    if selected_cell_ids[0] < 0 or selected_cell_ids[-1] >= coarse_grid.n_cells:
        raise ValueError('Cell IDs must identify cells in the coarse grid.')

    fractions = tuple(_axis_fractions(value) for value in scales)
    coarse_points = np.asarray(coarse_grid.points)
    coarse_axes = tuple(np.unique(coarse_points[:, axis]).astype(np.float64) for axis in range(_AXIS_COUNT))
    if np.prod([len(axis) for axis in coarse_axes]) != coarse_grid.n_points:
        raise ValueError('Adaptive refinement requires an axis-aligned Cartesian grid.')

    fine_axes: list[NDArray[np.float64]] = []
    subdivisions: list[int] = []
    for axis, axis_fractions in zip(coarse_axes, fractions, strict=True):
        subdivisions.append(len(axis_fractions) - 1)
        interval_points = axis[:-1, None] + np.diff(axis)[:, None] * axis_fractions[:-1]
        fine_axes.append(np.concatenate((interval_points.ravel(), axis[-1:])))

    coarse_centers = np.asarray(coarse_grid.cell_centers().points)[selected_cell_ids]
    coarse_indices = tuple(
        np.searchsorted(axis, coarse_centers[:, dimension], side='right') - 1
        for dimension, axis in enumerate(coarse_axes)
    )
    fine_indices = tuple(
        indices[:, None] * count + np.arange(count, dtype=np.int64)[None, :]
        for indices, count in zip(coarse_indices, subdivisions, strict=True)
    )
    fine_x = fine_indices[0][:, :, None, None]
    fine_y = fine_indices[1][:, None, :, None]
    fine_z = fine_indices[2][:, None, None, :]
    num_x_cells = len(fine_axes[0]) - 1
    num_y_cells = len(fine_axes[1]) - 1
    selected_fine_cells = (fine_x + num_x_cells * (fine_y + num_y_cells * fine_z)).ravel()

    full_grid = pv.RectilinearGrid(*fine_axes)
    result = full_grid.extract_cells(selected_fine_cells)

    # ``extract_cells`` may materialize float32 points even when the rectilinear
    # axes are float64. Reconstruct them from VTK's source IDs so tabulation sees
    # the same coordinates as direct subdivision.
    original_point_ids = np.asarray(result.point_data['vtkOriginalPointIds'])
    num_x_points = len(fine_axes[0])
    num_y_points = len(fine_axes[1])
    point_x = original_point_ids % num_x_points
    point_y = (original_point_ids // num_x_points) % num_y_points
    point_z = original_point_ids // (num_x_points * num_y_points)
    points = np.column_stack((fine_axes[0][point_x], fine_axes[1][point_y], fine_axes[2][point_z]))

    voxel_cells = result.cells.reshape(-1, _HEXAHEDRON_POINT_COUNT + 1)
    hexahedron_points = voxel_cells[:, 1:][:, _VOXEL_TO_HEXAHEDRON]
    hexahedron_cells = np.column_stack(
        (
            np.full(result.n_cells, _HEXAHEDRON_POINT_COUNT, dtype=np.int64),
            hexahedron_points,
        ),
    ).ravel()
    cell_types = np.full(result.n_cells, pv.CellType.HEXAHEDRON, dtype=np.uint8)
    return pv.UnstructuredGrid(hexahedron_cells, cell_types, points)
