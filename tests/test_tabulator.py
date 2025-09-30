"""Tests for the Tabulator class and related functions."""

from pathlib import Path

import numpy as np
import pytest

from tests._src_imports import (
    Tabulator,
    _cartesian_to_spherical,
    _spherical_to_cartesian,
    array_like_type,
)

MOLDEN_PATH = Path(__file__).with_name('sample_molden.inp')


def test_spherical_cartesian_roundtrip() -> None:
    """Test roundtrip conversion between spherical and Cartesian coordinates."""
    rng = np.random.default_rng(seed=42)
    r_vals = rng.uniform(0.1, 10.0, size=100)
    theta_vals = rng.uniform(0.0, np.pi, size=100)
    phi_vals = rng.uniform(-np.pi, np.pi, size=100)

    x, y, z = _spherical_to_cartesian(r_vals, theta_vals, phi_vals)
    r2, theta2, phi2 = _cartesian_to_spherical(x, y, z)

    np.testing.assert_allclose(r_vals, r2, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(theta_vals, theta2, rtol=1e-12, atol=1e-12)
    assert np.allclose(phi_vals, phi2)


def test_tabulate_gtos_requires_grid() -> None:
    """Test that tabulate_gtos raises RuntimeError if grid is not set."""
    tab = Tabulator(str(MOLDEN_PATH))
    with pytest.raises(RuntimeError):
        tab.tabulate_gtos()


def test_cartesian_grid_shape() -> None:
    """Test that the Cartesian grid is created with the correct shape."""
    tab = Tabulator(str(MOLDEN_PATH))
    x, y, z = np.linspace(-1, 1, 3), np.linspace(-1, 1, 4), np.linspace(-1, 1, 2)
    tab.cartesian_grid(x, y, z, tabulate_gtos=False)
    assert tab.grid is not None
    assert tab.grid.shape == (len(x) * len(y) * len(z), 3)


def test_spherical_grid_shape() -> None:
    """Test that the spherical grid is created with the correct shape."""
    tab = Tabulator(str(MOLDEN_PATH))
    r, theta, phi = np.r_[1.0, 2.0], np.r_[0.0, np.pi / 2, np.pi], np.r_[-np.pi, 0.0, np.pi / 2, np.pi]
    tab.spherical_grid(r, theta, phi, tabulate_gtos=False)
    assert tab.grid is not None
    assert tab.grid.shape == (len(r) * len(theta) * len(phi), 3)


@pytest.mark.parametrize(('lmax', 'num_points'), [(0, 10), (3, 25), (5, 50)])
def test_tabulate_xlms_shape(lmax: int, num_points: int) -> None:
    """Test that _tabulate_xlms returns an array of the correct shape."""
    theta = np.linspace(0.0, np.pi, num_points, dtype=float)
    phi = np.linspace(-np.pi, np.pi, num_points, dtype=float)

    xlms = Tabulator._tabulate_xlms(theta, phi, lmax)  # noqa: SLF001

    assert xlms.shape == (lmax + 1, 2 * lmax + 1, num_points)


@pytest.mark.parametrize('mo_inds', [None, 0, [0], [0, 1, 2], [0, 1, 2, 3, 4], range(1, 10)])
def test_tabulate_mos(mo_inds: int | array_like_type | None) -> None:
    """Test that tabulate_mos returns an array of the correct shape."""
    tab = Tabulator(str(MOLDEN_PATH))
    tab.cartesian_grid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))
    mo_data = tab.tabulate_mos(mo_inds)

    assert mo_data is not None

    if mo_inds is None:
        assert mo_data.shape == (125, 177)
    elif isinstance(mo_inds, int):
        assert mo_data.shape == (125,)
    else:
        assert mo_data.shape == (125, len(mo_inds))


@pytest.mark.parametrize('mo_inds', [-1, range(0), range(-1, 1), [0, -1], [1, 2, 3, -1], [0, 178]])
def test_invalid_mo_inds(mo_inds: int | array_like_type | None) -> None:
    """Test that tabulate_mos raises ValueError for invalid mo_inds."""
    tab = Tabulator(str(MOLDEN_PATH))
    tab.cartesian_grid(np.linspace(-1, 1, 5), np.linspace(-1, 1, 5), np.linspace(-1, 1, 5))

    with pytest.raises(ValueError, match=r'Provided mo_ind.* Please provide valid .*'):
        tab.tabulate_mos(mo_inds)


def test_export_cube_creates_file(tmp_path: Path) -> None:
    """Ensure exporting a cube file writes the expected artifact."""
    tab = Tabulator(str(MOLDEN_PATH))
    axis = np.linspace(-1.0, 1.0, 2)
    tab.cartesian_grid(axis, axis, axis)

    cube_file_path = tmp_path / 'orbital.cube'
    tab.export(cube_file_path, mo_index=0)

    cube_path = cube_file_path
    assert cube_path.exists()

    contents = cube_path.read_text(encoding='ascii').splitlines()
    assert 'Molecular orbital 0' in contents[1]
    header_tokens = contents[2].split()
    assert int(header_tokens[0]) == len(tab._parser.atoms)  # noqa: SLF001


def test_export_cube_requires_cartesian_grid(tmp_path: Path) -> None:
    """Cube export should fail when the grid is spherical."""
    tab = Tabulator(str(MOLDEN_PATH))
    r, theta, phi = np.r_[1.0, 2.0], np.r_[0.0, np.pi / 2], np.r_[-np.pi, 0.0]
    tab.spherical_grid(r, theta, phi)

    with pytest.raises(RuntimeError, match='Cube exports are only supported'):
        tab.export(tmp_path / 'orbital.cube', mo_index=0)


def test_export_cube_requires_mo_index(tmp_path: Path) -> None:
    """Cube export must receive an orbital index."""
    tab = Tabulator(str(MOLDEN_PATH))
    axis = np.linspace(-0.5, 0.5, 2)
    tab.cartesian_grid(axis, axis, axis)

    with pytest.raises(ValueError, match='Cube exports require'):
        tab.export(tmp_path / 'orbital.cube')


def test_export_vtk_writes_multiblock(tmp_path: Path) -> None:
    """VTK export should emit a multiblock file with molecule and atom data."""
    pv = pytest.importorskip('pyvista')

    tab = Tabulator(str(MOLDEN_PATH))
    axis = np.linspace(-0.5, 0.5, 2)
    tab.cartesian_grid(axis, axis, axis)

    vtk_path = tmp_path / 'dataset.vtk'
    tab.export(vtk_path)

    assert vtk_path.exists()

    mos_data = pv.read(vtk_path)
    assert isinstance(mos_data, pv.StructuredGrid)

    assert 'mo_0' in mos_data.point_data
