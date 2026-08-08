"""Tests for cell-local adaptive Cartesian grids."""
# ruff:file-ignore[import-private-name, undocumented-public-function, magic-value-comparison]

import numpy as np
import pytest
import pyvista as pv

from moldenViz._adaptive_grid import (
    crossed_cell_ids,
    normalize_scale,
    parse_scale,
    refined_grid,
)


def structured_grid(edge_points: int = 3) -> pv.StructuredGrid:
    """Return a small Cartesian structured grid.

    Returns
    -------
    pyvista.StructuredGrid
        Uniform test grid.
    """
    axis = np.linspace(-1.0, 1.0, edge_points)
    x, y, z = np.meshgrid(axis, axis, axis, indexing='ij')
    return pv.StructuredGrid(x, y, z)


def test_crossed_cells_are_unioned_across_signed_contours_and_mos() -> None:
    grid = structured_grid()
    x, y, _ = grid.points.T
    mo_values = np.column_stack((x, y))

    actual = crossed_cell_ids(grid, mo_values, 0.5)

    assert actual.size == grid.n_cells


def test_crossed_cells_uses_pyvista_contour_mapping(monkeypatch: pytest.MonkeyPatch) -> None:
    grid = structured_grid()
    calls: list[np.ndarray] = []
    original = pv.StructuredGrid.find_containing_cell

    def record_find(self: pv.StructuredGrid, points: np.ndarray) -> np.ndarray:
        calls.append(points)
        return np.asarray(original(self, points))

    monkeypatch.setattr(pv.StructuredGrid, 'find_containing_cell', record_find)
    crossed_cell_ids(grid, grid.points[:, 0], 0.5)

    assert calls
    assert calls[0].shape[1] == 3


def test_crossed_cells_handles_empty_contours() -> None:
    grid = structured_grid()
    actual = crossed_cell_ids(grid, np.zeros(grid.n_points), 2.0)
    assert actual.size == 0


@pytest.mark.parametrize(
    ('value', 'expected'),
    [
        ('1.5', (1.5, 1.5, 1.5)),
        ('(1.5, 2, 3.25)', (1.5, 2.0, 3.25)),
    ],
)
def test_parse_and_normalize_scale(value: str, expected: tuple[float, float, float]) -> None:
    assert normalize_scale(parse_scale(value)) == expected


@pytest.mark.parametrize('value', ['0.99', '(1, 0.5, 2)', '(1, 2)', 'not-a-scale'])
def test_parse_scale_rejects_invalid_values(value: str) -> None:
    with pytest.raises(ValueError, match=r'scale|Scale'):
        parse_scale(value)


def test_fractional_refinement_is_conforming_across_adjacent_cells() -> None:
    grid = structured_grid()
    fine = refined_grid(grid, [0, 1], (1.5, 2.0, 1.0))

    assert fine.n_cells == 8
    assert fine.n_points < fine.n_cells * 8
    assert set(fine.celltypes) == {pv.CellType.HEXAHEDRON}
    assert np.isclose(fine.bounds.x_min, -1.0)
    assert np.isclose(fine.bounds.x_max, 1.0)


@pytest.mark.parametrize(
    ('scale', 'expected_cells', 'expected_points'),
    [
        (1.5, 16, 45),
        ((2.0, 3.0, 2.5), 36, 80),
    ],
)
def test_vectorized_refinement_preserves_fractional_and_anisotropic_geometry(
    scale: float | tuple[float, float, float],
    expected_cells: int,
    expected_points: int,
) -> None:
    grid = structured_grid()
    fine = refined_grid(grid, [0, 1], scale)

    assert fine.n_cells == expected_cells
    assert fine.n_points == expected_points
    assert np.asarray(fine.points).dtype == np.float64
    assert set(fine.celltypes) == {pv.CellType.HEXAHEDRON}


def test_refinement_deduplicates_repeated_cell_ids() -> None:
    grid = structured_grid()
    expected = refined_grid(grid, [0, 1], 2.0)
    actual = refined_grid(grid, [1, 0, 1, 0], 2.0)

    assert actual.n_cells == expected.n_cells
    assert actual.n_points == expected.n_points
    np.testing.assert_allclose(
        np.unique(actual.points, axis=0),
        np.unique(expected.points, axis=0),
    )


def test_refined_cells_map_back_to_selected_coarse_cells() -> None:
    grid = structured_grid()
    selected = np.asarray([0, 3, 6])
    fine = refined_grid(grid, selected, (2.0, 3.0, 2.5))

    containing = grid.find_containing_cell(np.asarray(fine.cell_centers().points))
    actual_ids, counts = np.unique(containing, return_counts=True)

    np.testing.assert_array_equal(actual_ids, selected)
    np.testing.assert_array_equal(counts, np.full(selected.size, 18))


@pytest.mark.parametrize('cell_id', [-1, 8])
def test_refinement_rejects_invalid_cell_ids(cell_id: int) -> None:
    with pytest.raises(ValueError, match='Cell IDs'):
        refined_grid(structured_grid(), [cell_id], 2.0)


def test_scale_one_keeps_one_cell_per_selected_cell() -> None:
    grid = structured_grid()
    fine = refined_grid(grid, [0, 1], 1.0)
    assert fine.n_cells == 2


def test_refining_no_cells_returns_empty_grid() -> None:
    fine = refined_grid(structured_grid(), [], 5.0)
    assert fine.n_cells == 0
    assert fine.n_points == 0
