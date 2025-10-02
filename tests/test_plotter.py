"""Unit tests for the Plotter class."""

import re
from collections.abc import Iterator
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pytest

from tests._src_imports import GridType, Plotter, Tabulator, plotter_module, tabulator_module

# Use the same sample molden file as other tests
MOLDEN_PATH = Path(__file__).with_name('sample_molden.inp')


class MockTabulator(Tabulator):
    """Mock tabulator for testing error conditions."""

    def __init__(
        self,
        has_grid: bool = True,
        has_gto_data: bool = True,
        grid_type: GridType = GridType.SPHERICAL,
        original_axes: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None,
    ) -> None:
        # Always create attributes, then optionally remove them
        self._grid = np.array([[0, 0, 0], [1, 1, 1]])
        self.grid = self._grid
        self.gto_data = np.array([[1, 2, 3], [4, 5, 6]])
        self._grid_type = grid_type
        self._grid_dimensions = (2, 2, 2)
        self.original_axes: tuple[np.ndarray, np.ndarray, np.ndarray] | None

        # Set original_axes if provided, otherwise use a uniform default
        if original_axes is not None:
            self.original_axes = original_axes
        else:
            # Default uniform axes
            self.original_axes = (
                np.array([0.0, 1.0, 2.0]),
                np.array([0.0, 1.0, 2.0]),
                np.array([0.0, 1.0, 2.0]),
            )

        # Remove attributes if requested
        if not has_grid:
            delattr(self, 'grid')
        if not has_gto_data:
            delattr(self, 'gto_data')

        # Mock parser with atoms for Molecule creation
        self._parser = Mock()
        self._parser.atoms = []


@pytest.fixture(autouse=True)
def patch_plotter_dependencies() -> Iterator[None]:
    """Patch GUI-heavy dependencies used by Plotter.

    Yields
    ------
    None
        Allows tests to run with patched GUI components.
    """
    patcher_tk = patch.object(plotter_module, 'tk')
    patcher_pv = patch.object(plotter_module, 'pv')
    patcher_background_plotter = patch.object(plotter_module, 'BackgroundPlotter')
    patcher_molecule = patch.object(plotter_module, 'Molecule')

    mock_tk = patcher_tk.start()
    mock_pv = patcher_pv.start()
    mock_background_plotter = patcher_background_plotter.start()
    mock_molecule = patcher_molecule.start()

    mock_tk.return_value = Mock()
    mock_pv.return_value = Mock()

    mock_plotter_instance = Mock()
    mock_plotter_instance.show_axes.return_value = None
    mock_background_plotter.return_value = mock_plotter_instance

    mock_molecule_instance = Mock()
    mock_molecule_instance.add_meshes.return_value = []
    mock_molecule_instance.max_radius = 5.0
    mock_molecule.return_value = mock_molecule_instance

    try:
        yield
    finally:
        patcher_tk.stop()
        patcher_pv.stop()
        patcher_background_plotter.stop()
        patcher_molecule.stop()


def test_tabulator_missing_grid_attribute_raises_error() -> None:
    """Test that Plotter requires tabulators with a grid attribute."""
    mock_tabulator = MockTabulator(has_grid=False, has_gto_data=True)

    with pytest.raises(ValueError, match=re.escape('Tabulator does not have grid attribute.')):
        Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)


def test_tabulator_missing_gto_data_with_only_molecule_false_raises_error() -> None:
    """Test that Plotter enforces presence of GTO data when needed."""
    mock_tabulator = MockTabulator(has_grid=True, has_gto_data=False)

    with pytest.raises(ValueError, match=re.escape('Tabulator does not have tabulated GTOs.')):
        Plotter(str(MOLDEN_PATH), only_molecule=False, tabulator=mock_tabulator)


def test_tabulator_missing_gto_data_with_only_molecule_true_succeeds() -> None:
    """Test that Plotter allows missing GTO data when only molecules are rendered."""
    mock_tabulator = MockTabulator(has_grid=True, has_gto_data=False)

    plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)

    assert plotter is not None


def test_tabulator_unknown_grid_type_raises_error() -> None:
    """Test that unsupported grid types raise a ValueError."""
    mock_tabulator = MockTabulator(has_grid=True, has_gto_data=True, grid_type=GridType.UNKNOWN)

    with pytest.raises(
        ValueError,
        match=re.escape('The plotter only supports spherical and cartesian grids.'),
    ):
        Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)


def test_tabulator_spherical_grid_type_succeeds() -> None:
    """Test that spherical tabulator grids are accepted."""
    mock_tabulator = MockTabulator(has_grid=True, has_gto_data=True, grid_type=GridType.SPHERICAL)

    plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)

    assert plotter.tabulator == mock_tabulator


def test_tabulator_cartesian_grid_type_succeeds() -> None:
    """Test that cartesian tabulator grids are accepted."""
    mock_tabulator = MockTabulator(has_grid=True, has_gto_data=True, grid_type=GridType.CARTESIAN)

    plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)

    assert plotter.tabulator == mock_tabulator


def test_valid_tabulator_with_all_attributes_succeeds() -> None:
    """Test that a fully populated tabulator can be used directly."""
    mock_tabulator = MockTabulator(has_grid=True, has_gto_data=True, grid_type=GridType.SPHERICAL)

    plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)

    assert plotter.tabulator == mock_tabulator
    assert plotter.tabulator is mock_tabulator


def test_none_tabulator_creates_default_tabulator() -> None:
    """Test that Plotter instantiates a Tabulator when none is provided."""
    with patch.object(plotter_module, 'Tabulator') as mock_tabulator_class:
        mock_tabulator_instance = Mock()
        mock_tabulator_instance._parser.atoms = []  # noqa: SLF001
        mock_tabulator_instance._grid_type = GridType.SPHERICAL  # noqa: SLF001
        mock_tabulator_instance.grid = np.array([[0, 0, 0]])
        mock_tabulator_instance.grid_dimensions = (1, 1, 1)
        mock_tabulator_class.return_value = mock_tabulator_instance

        plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=None)

    mock_tabulator_class.assert_called_once_with(str(MOLDEN_PATH), only_molecule=True)
    assert plotter.tabulator == mock_tabulator_instance


def test_real_tabulator_instance_validation() -> None:
    """Test validation when using a real Tabulator lacking a grid."""
    with patch.object(tabulator_module, 'Parser') as mock_parser_class:
        mock_parser = Mock()
        mock_parser.atoms = []
        mock_parser_class.return_value = mock_parser

        real_tabulator = Tabulator(str(MOLDEN_PATH), only_molecule=True)

    assert real_tabulator._grid_type == GridType.UNKNOWN  # noqa: SLF001
    assert not hasattr(real_tabulator, 'grid')

    with pytest.raises(ValueError, match=re.escape('Tabulator does not have grid attribute.')):
        Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=real_tabulator)


def test_non_uniform_grid_x_axis_raises_error() -> None:
    """Test that non-uniform x-axis spacing is rejected."""
    non_uniform_x = np.array([0.0, 1.0, 2.5])
    uniform_y = np.array([0.0, 1.0, 2.0])
    uniform_z = np.array([0.0, 1.0, 2.0])

    mock_tabulator = MockTabulator(
        has_grid=True,
        has_gto_data=True,
        grid_type=GridType.CARTESIAN,
        original_axes=(non_uniform_x, uniform_y, uniform_z),
    )

    with pytest.raises(ValueError, match=re.escape('x-axis must be evenly spaced.')):
        Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)


def test_non_uniform_grid_y_axis_raises_error() -> None:
    """Test that non-uniform y-axis spacing is rejected."""
    uniform_x = np.array([0.0, 1.0, 2.0])
    non_uniform_y = np.array([0.0, 1.0, 3.0])
    uniform_z = np.array([0.0, 1.0, 2.0])

    mock_tabulator = MockTabulator(
        has_grid=True,
        has_gto_data=True,
        grid_type=GridType.CARTESIAN,
        original_axes=(uniform_x, non_uniform_y, uniform_z),
    )

    with pytest.raises(ValueError, match=re.escape('y-axis must be evenly spaced.')):
        Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)


def test_non_uniform_grid_z_axis_raises_error() -> None:
    """Test that non-uniform z-axis spacing is rejected."""
    uniform_x = np.array([0.0, 1.0, 2.0])
    uniform_y = np.array([0.0, 1.0, 2.0])
    non_uniform_z = np.array([0.0, 0.5, 2.0])

    mock_tabulator = MockTabulator(
        has_grid=True,
        has_gto_data=True,
        grid_type=GridType.CARTESIAN,
        original_axes=(uniform_x, uniform_y, non_uniform_z),
    )

    with pytest.raises(ValueError, match=re.escape('z-axis must be evenly spaced.')):
        Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)


def test_uniform_grid_succeeds() -> None:
    """Test that uniform axes pass Plotter validation."""
    uniform_x = np.array([0.0, 1.0, 2.0, 3.0])
    uniform_y = np.array([0.0, 0.5, 1.0, 1.5])
    uniform_z = np.array([-1.0, 0.0, 1.0])

    mock_tabulator = MockTabulator(
        has_grid=True,
        has_gto_data=True,
        grid_type=GridType.CARTESIAN,
        original_axes=(uniform_x, uniform_y, uniform_z),
    )

    plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)

    assert plotter.tabulator == mock_tabulator


def test_no_original_axes_succeeds() -> None:
    """Test that missing original axes does not raise an error."""
    mock_tabulator = MockTabulator(
        has_grid=True,
        has_gto_data=True,
        grid_type=GridType.CARTESIAN,
        original_axes=None,
    )
    mock_tabulator.original_axes = None

    plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)

    assert plotter is not None


@pytest.fixture
def plotter_with_mock_tabulator() -> tuple[Plotter, MockTabulator]:
    """Return a plotter paired with a mock tabulator.

    Returns
    -------
    tuple[Plotter, MockTabulator]
        The plotter under test and its backing tabulator mock.
    """
    mock_tabulator = MockTabulator(has_grid=True, has_gto_data=True, grid_type=GridType.SPHERICAL)
    mock_tabulator.cartesian_grid = Mock()
    mock_tabulator.spherical_grid = Mock()

    plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)

    return plotter, mock_tabulator


def test_update_mesh_cartesian_grid(plotter_with_mock_tabulator: tuple[Plotter, MockTabulator]) -> None:
    """Test that cartesian grids trigger the cartesian update hook."""
    plotter, mock_tabulator = plotter_with_mock_tabulator
    x_points = np.array([0, 1, 2])
    y_points = np.array([0, 1])
    z_points = np.array([0, 1, 2, 3])

    plotter.update_mesh(x_points, y_points, z_points, GridType.CARTESIAN)

    mock_tabulator.cartesian_grid.assert_called_once_with(x_points, y_points, z_points)
    mock_tabulator.spherical_grid.assert_not_called()


def test_update_mesh_spherical_grid(plotter_with_mock_tabulator: tuple[Plotter, MockTabulator]) -> None:
    """Test that spherical grids trigger the spherical update hook."""
    plotter, mock_tabulator = plotter_with_mock_tabulator
    r_points = np.array([0, 1, 2])
    theta_points = np.array([0, np.pi / 2, np.pi])
    phi_points = np.array([0, np.pi, 2 * np.pi])

    plotter.update_mesh(r_points, theta_points, phi_points, GridType.SPHERICAL)

    mock_tabulator.spherical_grid.assert_called_once_with(r_points, theta_points, phi_points)
    mock_tabulator.cartesian_grid.assert_not_called()


def test_update_mesh_unknown_grid_type_raises_error(
    plotter_with_mock_tabulator: tuple[Plotter, MockTabulator],
) -> None:
    """Test that update_mesh rejects unsupported grid types."""
    plotter, mock_tabulator = plotter_with_mock_tabulator
    points = np.array([0, 1, 2])

    with pytest.raises(
        ValueError,
        match=re.escape('The plotter only supports spherical and cartesian grids.'),
    ):
        plotter.update_mesh(points, points, points, GridType.UNKNOWN)

    mock_tabulator.cartesian_grid.assert_not_called()
    mock_tabulator.spherical_grid.assert_not_called()
