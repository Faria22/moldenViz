"""Tests for the Tabulator class and related functions."""

from __future__ import annotations

import warnings
from math import factorial
from pathlib import Path

import numpy as np
import pytest

import moldenViz.tabulator as tabulator_module
from moldenViz import Atom, GaussianPrimitive, Shell
from moldenViz.tabulator import GridType, Tabulator

MOLDEN_PATH = Path(__file__).with_name('sample_molden.inp')


def _tabulate_xlms(theta: np.ndarray, phi: np.ndarray, lmax: int) -> np.ndarray:
    """Return normalized real spherical harmonics for test comparisons.

    Returns
    -------
    np.ndarray
        Real spherical harmonics indexed by ``[l, m, point]``.
    """
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    plms = np.zeros((lmax + 1, lmax + 1, theta.size), dtype=float)

    for m in range(lmax + 1):
        if m == 0:
            plms[0, 0] = 1.0
        else:
            plms[m, m] = (2 * m - 1) * sin_theta * plms[m - 1, m - 1]

        if m < lmax:
            plms[m + 1, m] = (2 * m + 1) * cos_theta * plms[m, m]

        for l in range(m + 2, lmax + 1):
            plms[l, m] = ((2 * l - 1) * cos_theta * plms[l - 1, m] - (l + m - 1) * plms[l - 2, m]) / (l - m)

    xlms = np.zeros((lmax + 1, 2 * lmax + 1, theta.size), dtype=float)
    for l in range(lmax + 1):
        for m in range(l + 1):
            normalization = np.sqrt(
                (2 * l + 1) * factorial(l - m) / (4 * np.pi * factorial(l + m)),
            )
            if m == 0:
                xlms[l, 0] = normalization * plms[l, 0]
            else:
                scale = np.sqrt(2) * normalization * plms[l, m]
                xlms[l, -m] = scale * np.sin(m * phi)
                xlms[l, m] = scale * np.cos(m * phi)

    return xlms


def test_spherical_cartesian_roundtrip() -> None:
    """Test roundtrip conversion between spherical and Cartesian coordinates."""
    rng = np.random.default_rng(seed=42)
    r_vals = rng.uniform(0, 10.0, size=100)
    theta_vals = rng.uniform(0.0, np.pi, size=100)
    phi_vals = rng.uniform(-np.pi, np.pi, size=100)

    x, y, z = Tabulator.spherical_to_cartesian(r_vals, theta_vals, phi_vals)
    r2, theta2, phi2 = Tabulator.cartesian_to_spherical(x, y, z)

    np.testing.assert_allclose(r_vals, r2, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(theta_vals, theta2, rtol=1e-12, atol=1e-12)
    assert np.allclose(phi_vals, phi2)


def test_cartesian_to_spherical_handles_zero_radius() -> None:
    """Zero-radius points should not emit warnings or NaNs."""
    x = np.array([0.0])
    y = np.array([0.0])
    z = np.array([0.0])

    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        r, theta, phi = Tabulator.cartesian_to_spherical(x, y, z)

    assert np.isfinite(r).all()
    assert np.isfinite(theta).all()
    assert np.isfinite(phi).all()
    np.testing.assert_allclose(r, 0.0)
    np.testing.assert_allclose(theta, 0.0)


def test_tabulate_gtos_requires_grid() -> None:
    """Test that tabulate_gtos raises RuntimeError if grid is not set."""
    tab = Tabulator(str(MOLDEN_PATH))
    with pytest.raises(RuntimeError):
        tab.tabulate_gtos()


def test_tabulate_gtos_cached_values_cover_all_coeffs() -> None:
    """Ensure tabulate_gtos populates every MO coefficient on the grid."""
    tab = Tabulator(str(MOLDEN_PATH))
    axis = np.linspace(-1.0, 1.0, 4)
    tab.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    gto_data = tab.tabulate_gtos()

    expected_points = axis.size**3
    expected_coeffs = tab._parser.mo_coeffs.shape[1]  # ruff:ignore[private-member-access]
    assert gto_data.shape == (expected_points, expected_coeffs)
    assert np.all(np.isfinite(gto_data))
    assert tab.gtos is gto_data  # Cached array should match the return value
    assert tab.has_gtos


def test_compute_gtos_uses_explicit_grid_without_updating_cache() -> None:
    """Explicit-grid computation should not read or update live grid state."""
    tab = Tabulator(str(MOLDEN_PATH))
    live_axis = np.linspace(-1.0, 1.0, 2)
    tab.cartesian_grid(live_axis, live_axis, live_axis, tabulate_gtos=False)
    live_grid = tab.grid.copy()
    snapshot = np.array([[0.0, 0.0, 0.0], [0.25, -0.5, 0.75]])

    gtos = tab.compute_gtos(snapshot)

    assert gtos.shape == (snapshot.shape[0], tab._parser.mo_coeffs.shape[1])  # ruff:ignore[private-member-access]
    np.testing.assert_array_equal(tab.grid, live_grid)
    assert not tab.has_gtos


@pytest.mark.parametrize('point_chunk_size', [1, 17, 10_000])
def test_compute_gtos_chunks_match_full_grid(point_chunk_size: int) -> None:
    """Point chunks should preserve GTO values, shape, and basis ordering."""
    tab = Tabulator(str(MOLDEN_PATH))
    axis = np.linspace(-1.0, 1.0, 5)
    tab.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    expected = tab.compute_gtos(tab.grid, point_chunk_size=None)
    actual = tab.compute_gtos(tab.grid, point_chunk_size=point_chunk_size)

    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize('point_chunk_size', [0, -1, True, 1.5])
def test_compute_gtos_rejects_invalid_point_chunk_size(point_chunk_size: object) -> None:
    """Chunk sizes must be positive integers or None."""
    tab = Tabulator(str(MOLDEN_PATH))
    axis = np.linspace(-1.0, 1.0, 2)
    tab.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    with pytest.raises(ValueError, match='positive integer or None'):
        tab.compute_gtos(tab.grid, point_chunk_size=point_chunk_size)  # type: ignore[arg-type]


def test_compute_gtos_default_bounds_worker_point_slices(monkeypatch: pytest.MonkeyPatch) -> None:
    """The default policy should keep each worker task at or below 32,768 points."""
    tab = Tabulator(str(MOLDEN_PATH))
    grid = np.zeros((32_769, 3))
    chunk_lengths: list[int] = []
    block_shapes: list[tuple[int, int]] = []

    def fake_tabulate_atom(
        chunk: np.ndarray,
        _atom: Atom,
        atom_block: np.ndarray,
    ) -> None:
        chunk_lengths.append(chunk.shape[0])
        block_shapes.append(atom_block.shape)
        atom_block[:] = 0.0

    monkeypatch.setattr('moldenViz.tabulator.os.cpu_count', lambda: 1)
    monkeypatch.setattr(tab, '_tabulate_atom', fake_tabulate_atom)

    actual = tab.compute_gtos(grid)

    num_atoms = len(tab._parser.atoms)  # ruff:ignore[private-member-access]
    expected_block_shapes = [
        (chunk_length, sum(2 * shell.l + 1 for shell in atom.shells))
        for atom in tab._parser.atoms  # ruff:ignore[private-member-access]
        for chunk_length in (32_768, 1)
    ]
    assert actual.shape == (grid.shape[0], tab._parser.mo_coeffs.shape[1])  # ruff:ignore[private-member-access]
    assert chunk_lengths == [32_768, 1] * num_atoms
    assert block_shapes == expected_block_shapes


def test_gto_worker_policy_is_bounded(monkeypatch: pytest.MonkeyPatch) -> None:
    """Default and explicit worker counts should respect CPU and atom limits."""
    default_worker_count = 4
    explicit_worker_count = 2
    monkeypatch.setattr(tabulator_module.os, 'cpu_count', lambda: 64)

    default_tabulator = Tabulator(str(MOLDEN_PATH))
    explicit_tabulator = Tabulator(str(MOLDEN_PATH), max_workers=explicit_worker_count)
    capped_tabulator = Tabulator(str(MOLDEN_PATH), max_workers=64)

    assert default_tabulator.max_workers == default_worker_count
    assert explicit_tabulator.max_workers == explicit_worker_count
    assert capped_tabulator.max_workers == default_worker_count

    monkeypatch.setattr(tabulator_module.os, 'cpu_count', lambda: 1)
    assert default_tabulator.max_workers == 1
    assert explicit_tabulator.max_workers == 1
    assert capped_tabulator.max_workers == 1


def test_default_gto_workers_switch_to_sequential_for_large_grids() -> None:
    """Default concurrency should avoid costly large-grid memory amplification."""
    largest_parallel_grid = 125_000
    tabulator = Tabulator(str(MOLDEN_PATH))
    explicit_tabulator = Tabulator(str(MOLDEN_PATH), max_workers=4)

    assert tabulator._workers_for_grid(largest_parallel_grid) == tabulator.max_workers  # ruff:ignore[private-member-access]
    assert tabulator._workers_for_grid(largest_parallel_grid + 1) == 1  # ruff:ignore[private-member-access]
    assert (
        explicit_tabulator._workers_for_grid(largest_parallel_grid + 1)  # ruff:ignore[private-member-access]
        == explicit_tabulator.max_workers
    )


def test_gto_worker_policy_rejects_invalid_limits() -> None:
    """Worker limits must be positive integers."""
    with pytest.raises(TypeError, match='positive integer'):
        Tabulator(str(MOLDEN_PATH), max_workers=True)
    with pytest.raises(ValueError, match='at least 1'):
        Tabulator(str(MOLDEN_PATH), max_workers=0)


@pytest.mark.parametrize(
    'molden_path',
    [
        MOLDEN_PATH,
        Path(__file__).parents[1] / 'src/moldenViz/examples/molden_files/pyridine.inp',
    ],
)
def test_parallel_and_sequential_gtos_are_equivalent(molden_path: Path) -> None:
    """Bounded parallel work should preserve sequential numerical results."""
    axis = np.linspace(-2.0, 2.0, 5)
    sequential = Tabulator(str(molden_path), max_workers=1)
    parallel = Tabulator(str(molden_path), max_workers=4)
    sequential.cartesian_grid(axis, axis, axis, tabulate_gtos=False)
    parallel.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    sequential_gtos = sequential.tabulate_gtos()
    parallel_gtos = parallel.tabulate_gtos()

    np.testing.assert_array_equal(parallel_gtos, sequential_gtos)


def test_gto_calls_reuse_process_executor(monkeypatch: pytest.MonkeyPatch) -> None:
    """Parallel calls should not construct a fresh executor."""

    def fail_executor_creation(*_args: object, **_kwargs: object) -> None:
        raise AssertionError('GTO tabulation created a new executor.')

    monkeypatch.setattr(tabulator_module, 'ThreadPoolExecutor', fail_executor_creation)
    tabulator = Tabulator(str(MOLDEN_PATH), max_workers=2)
    axis = np.linspace(-1.0, 1.0, 3)
    tabulator.cartesian_grid(axis, axis, axis, tabulate_gtos=False)

    tabulator.tabulate_gtos()
    tabulator.tabulate_gtos()


def test_tabulate_atom_reuses_exponentials_for_compatible_shells(monkeypatch: pytest.MonkeyPatch) -> None:
    """Compatible shells should share exponentials while retaining their prefactors."""

    def normalized_shell(l: int, exponents: list[float], coefficients: list[float]) -> Shell:
        shell = Shell(l, [GaussianPrimitive(exp, coeff) for exp, coeff in zip(exponents, coefficients, strict=True)])
        shell._normalize()  # ruff:ignore[private-member-access]
        return shell

    s_shell = normalized_shell(0, [0.5, 1.5], [0.25, 0.75])
    p_shell = normalized_shell(1, [0.5, 1.5], [0.8, 0.2])
    d_shell = normalized_shell(2, [2.0], [1.0])
    atom = Atom('X', 0, np.zeros(3), [s_shell, p_shell, d_shell])

    tab = Tabulator(str(MOLDEN_PATH))
    grid = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, -0.25, 0.75],
            [-1.0, 0.5, 0.25],
            [1.5, 1.0, -0.5],
        ],
    )
    tab.set_grid(grid)
    actual = np.empty((grid.shape[0], 9))

    original_exp = np.exp
    exponential_shapes: list[tuple[int, ...]] = []

    def tracked_exp(values: np.ndarray) -> np.ndarray:
        exponential_shapes.append(values.shape)
        return original_exp(values)

    monkeypatch.setattr(np, 'exp', tracked_exp)
    tab._tabulate_atom(grid, atom, actual)  # ruff:ignore[private-member-access]

    r_sq = np.einsum('ij,ij->i', grid, grid)
    solid_harmonics = Tabulator._tabulate_real_solid_harmonics(grid, 2)  # ruff:ignore[private-member-access]
    expected = np.empty_like(actual)
    block_cursor = 0
    for shell in atom.shells:
        num_m = 2 * shell.l + 1
        contraction = shell._prefactor @ original_exp(  # ruff:ignore[private-member-access]
            -shell._gto_exps[:, None] * r_sq[None, :],  # ruff:ignore[private-member-access]
        )
        expected[:, block_cursor : block_cursor + num_m] = (
            contraction[:, None] * solid_harmonics[shell.l, np.arange(-shell.l, shell.l + 1), ...].T
        )
        block_cursor += num_m

    assert exponential_shapes == [(2, grid.shape[0]), (1, grid.shape[0])]
    assert not np.array_equal(s_shell._prefactor, p_shell._prefactor)  # ruff:ignore[private-member-access]
    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=1e-14)


def test_clear_gtos_releases_cache_and_reports_missing_data() -> None:
    """Manual cache eviction should retain the grid and expose a clear state."""
    tab = Tabulator(str(MOLDEN_PATH))
    axis = np.linspace(-1.0, 1.0, 2)
    tab.cartesian_grid(axis, axis, axis)
    grid = tab.grid

    tab.clear_gtos()

    assert not tab.has_gtos
    assert tab._gtos is None  # ruff:ignore[private-member-access]
    assert tab.grid is grid
    with pytest.raises(RuntimeError, match=r'Call tabulate_gtos\(\) first'):
        _ = tab.gtos

    tab.clear_gtos()


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


def test_set_grid_is_the_explicit_arbitrary_grid_mutator() -> None:
    """Arbitrary grids should reset structured metadata and cached GTOs."""
    tab = Tabulator(str(MOLDEN_PATH))
    axis = np.linspace(-1.0, 1.0, 2)
    tab.cartesian_grid(axis, axis, axis)

    new_grid = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    tab.set_grid(new_grid)

    assert tab.grid is new_grid
    assert tab.grid_type is GridType.UNKNOWN
    assert tab.grid_dimensions == (0, 0, 0)
    assert tab.grid_axes is None
    assert not tab.has_gtos


def test_set_grid_exits_early_when_only_molecule_is_enabled() -> None:
    """Molecule-only tabulators should reject grid creation before validation."""
    tab = Tabulator(str(MOLDEN_PATH), only_molecule=True)

    with pytest.raises(RuntimeError, match='Grid creation is not allowed'):
        tab.set_grid(None)


def test_private_set_grid_rejects_unknown_grid_type() -> None:
    """Structured grids require a known coordinate system."""
    tab = Tabulator(str(MOLDEN_PATH))
    axis = np.linspace(-1.0, 1.0, 2)

    with pytest.raises(ValueError, match='Grid type cannot be unknown'):
        tab._set_grid(axis, axis, axis, GridType.UNKNOWN)  # ruff:ignore[private-member-access]


def test_grid_property_is_read_only() -> None:
    """Grid replacement should go through ``set_grid`` rather than assignment."""
    tab = Tabulator(str(MOLDEN_PATH))

    with pytest.raises(AttributeError):
        tab.grid = np.zeros((1, 3))  # type: ignore[misc]


@pytest.mark.parametrize(('lmax', 'num_points'), [(0, 10), (3, 25), (5, 50)])
def test_tabulate_xlms_shape(lmax: int, num_points: int) -> None:
    """Test that _tabulate_xlms returns an array of the correct shape."""
    theta = np.linspace(0.0, np.pi, num_points, dtype=float)
    phi = np.linspace(-np.pi, np.pi, num_points, dtype=float)

    xlms = _tabulate_xlms(theta, phi, lmax)

    assert xlms.shape == (lmax + 1, 2 * lmax + 1, num_points)


@pytest.mark.parametrize('l', range(5))
def test_cartesian_solid_harmonics_match_spherical_kernel(l: int) -> None:
    """Cartesian polynomials should preserve every supported l, m value."""
    rng = np.random.default_rng(seed=8300 + l)
    points = rng.uniform(-4.0, 4.0, size=(500, 3))
    r, theta, phi = Tabulator.cartesian_to_spherical(*points.T)

    actual = Tabulator._tabulate_real_solid_harmonics(points, l)  # ruff:ignore[private-member-access]
    spherical = _tabulate_xlms(theta, phi, l)

    for m in range(-l, l + 1):
        expected = r**l * spherical[l, m]
        np.testing.assert_allclose(actual[l, m], expected, rtol=1e-11, atol=1e-11)


def test_cartesian_solid_harmonics_handle_origin_and_axes() -> None:
    """Every supported solid harmonic should be finite at Cartesian edge cases."""
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ],
    )

    actual = Tabulator._tabulate_real_solid_harmonics(points, 4)  # ruff:ignore[private-member-access]
    r, theta, phi = Tabulator.cartesian_to_spherical(*points.T)
    spherical = _tabulate_xlms(theta, phi, 4)

    assert np.all(np.isfinite(actual))
    for l in range(5):
        for m in range(-l, l + 1):
            np.testing.assert_allclose(actual[l, m], r**l * spherical[l, m], rtol=1e-10, atol=5e-8)
    np.testing.assert_allclose(actual[0, 0], 1 / np.sqrt(4 * np.pi))
    np.testing.assert_allclose(actual[1:, :, 0], 0.0, atol=0.0)
    np.testing.assert_allclose(actual[1, 1], np.sqrt(3 / (4 * np.pi)) * points[:, 0])
    np.testing.assert_allclose(actual[1, -1], np.sqrt(3 / (4 * np.pi)) * points[:, 1])
    np.testing.assert_allclose(actual[1, 0], np.sqrt(3 / (4 * np.pi)) * points[:, 2])


@pytest.mark.parametrize(
    'molden_path',
    [
        MOLDEN_PATH,
        Path(__file__).parents[1] / 'src/moldenViz/examples/molden_files/pyridine.inp',
    ],
)
def test_cartesian_gto_tabulation_matches_spherical_implementation(
    monkeypatch: pytest.MonkeyPatch,
    molden_path: Path,
) -> None:
    """Representative Molden inputs should match the previous GTO kernel."""
    tab = Tabulator(str(molden_path))
    axis = np.linspace(-2.0, 2.0, 7)
    tab.cartesian_grid(axis, axis, axis, tabulate_gtos=False)
    actual = tab.tabulate_gtos()

    def spherical_solid_harmonics(
        centered_grid: np.ndarray,
        lmax: int,
    ) -> np.ndarray:
        r, theta, phi = Tabulator.cartesian_to_spherical(*centered_grid.T)
        xlms = _tabulate_xlms(theta, phi, lmax)
        return xlms * np.stack([r**l for l in range(lmax + 1)])[:, None, :]

    monkeypatch.setattr(
        Tabulator,
        '_tabulate_real_solid_harmonics',
        staticmethod(spherical_solid_harmonics),
    )
    expected = tab.tabulate_gtos()

    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=2e-8)


@pytest.mark.parametrize('mo_inds', [None, 0, [0], [0, 1, 2], [0, 1, 2, 3, 4], range(1, 10)])
def test_tabulate_mos(mo_inds: int | list[int] | range | None) -> None:
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


@pytest.mark.parametrize('mo_inds', [0, [0, 1, 2], None])
def test_tabulate_mos_matches_sum_reduction(mo_inds: int | list[int] | None) -> None:
    """Matrix contractions should match the previous sum reduction."""
    tab = Tabulator(str(MOLDEN_PATH))
    axis = np.linspace(-1.0, 1.0, 5)
    tab.cartesian_grid(axis, axis, axis)

    mo_data = tab.tabulate_mos(mo_inds)
    if isinstance(mo_inds, int):
        expected = np.sum(tab.gtos * tab._parser.mo_coeffs[mo_inds][None, :], axis=1)  # ruff:ignore[private-member-access]
    else:
        indices = list(range(len(tab._parser.mos))) if mo_inds is None else mo_inds  # ruff:ignore[private-member-access]
        mo_coeffs = tab._parser.mo_coeffs[indices]  # ruff:ignore[private-member-access]
        expected = np.sum(tab.gtos[:, None, :] * mo_coeffs[None, ...], axis=2)

    np.testing.assert_allclose(mo_data, expected, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize('mo_inds', [-1, range(0), range(-1, 1), [0, -1], [1, 2, 3, -1], [0, 178]])
def test_invalid_mo_inds(mo_inds: int | list[int] | range | None) -> None:
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
    assert int(header_tokens[0]) == len(tab._parser.atoms)  # ruff:ignore[private-member-access]


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
