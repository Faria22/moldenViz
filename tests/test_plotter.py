# ruff: noqa: PT027, PT009, SLF001
# type: ignore[reportAttributeAccessIssue]
"""Unit tests for the Plotter class."""

import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np

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
    ) -> None:
        # Always create attributes, then optionally remove them
        self._grid = np.array([[0, 0, 0], [1, 1, 1]])
        self.grid = self._grid
        self.gto_data = np.array([[1, 2, 3], [4, 5, 6]])
        self._grid_type = grid_type
        self._grid_dimensions = (2, 2, 2)

        # Remove attributes if requested
        if not has_grid:
            delattr(self, 'grid')
        if not has_gto_data:
            delattr(self, 'gto_data')

        # Mock parser with atoms for Molecule creation
        self._parser = Mock()
        self._parser.atoms = []


class TestPlotterTabulatorValidation(unittest.TestCase):
    """Test validation of tabulator parameter in Plotter.__init__."""

    def setUp(self) -> None:
        """Set up mocks to avoid GUI components during testing."""
        # Mock all GUI-related components
        self.patcher_tk = patch.object(plotter_module, 'tk')
        self.patcher_pv = patch.object(plotter_module, 'pv')
        self.patcher_background_plotter = patch.object(plotter_module, 'BackgroundPlotter')
        self.patcher_molecule = patch.object(plotter_module, 'Molecule')

        self.mock_tk = self.patcher_tk.start()
        self.mock_pv = self.patcher_pv.start()
        self.mock_background_plotter = self.patcher_background_plotter.start()
        self.mock_molecule = self.patcher_molecule.start()

        # Configure mocks
        mock_plotter_instance = Mock()
        mock_plotter_instance.show_axes.return_value = None
        self.mock_background_plotter.return_value = mock_plotter_instance

        mock_molecule_instance = Mock()
        mock_molecule_instance.add_meshes.return_value = []
        mock_molecule_instance.max_radius = 5.0
        self.mock_molecule.return_value = mock_molecule_instance

    def tearDown(self) -> None:
        """Clean up patches."""
        self.patcher_tk.stop()
        self.patcher_pv.stop()
        self.patcher_background_plotter.stop()
        self.patcher_molecule.stop()

    def test_tabulator_missing_grid_attribute_raises_error(self) -> None:
        """Test that ValueError is raised when tabulator lacks grid attribute."""
        mock_tabulator = MockTabulator(has_grid=False, has_gto_data=True)

        with self.assertRaises(ValueError) as context:
            Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)

        self.assertEqual(str(context.exception), 'Tabulator does not have grid attribute.')

    def test_tabulator_missing_gto_data_with_only_molecule_false_raises_error(self) -> None:
        """Test that ValueError is raised when tabulator lacks gto_data and only_molecule=False."""
        mock_tabulator = MockTabulator(has_grid=True, has_gto_data=False)

        with self.assertRaises(ValueError) as context:
            Plotter(str(MOLDEN_PATH), only_molecule=False, tabulator=mock_tabulator)

        self.assertEqual(str(context.exception), 'Tabulator does not have tabulated GTOs.')

    def test_tabulator_missing_gto_data_with_only_molecule_true_succeeds(self) -> None:
        """Test that missing gto_data is allowed when only_molecule=True."""
        mock_tabulator = MockTabulator(has_grid=True, has_gto_data=False)

        # This should not raise an exception
        try:
            plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)
            # If we get here, the test passed
            self.assertIsNotNone(plotter)
        except ValueError:
            self.fail('Plotter raised ValueError when only_molecule=True and gto_data is missing')

    def test_tabulator_unknown_grid_type_raises_error(self) -> None:
        """Test that ValueError is raised when tabulator has UNKNOWN grid type."""
        mock_tabulator = MockTabulator(has_grid=True, has_gto_data=True, grid_type=GridType.UNKNOWN)

        with self.assertRaises(ValueError) as context:
            Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)

        self.assertEqual(str(context.exception), 'The plotter only supports spherical and cartesian grids.')

    def test_tabulator_spherical_grid_type_succeeds(self) -> None:
        """Test that spherical grid type is accepted."""
        mock_tabulator = MockTabulator(has_grid=True, has_gto_data=True, grid_type=GridType.SPHERICAL)

        # This should not raise an exception
        try:
            plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)
            self.assertIsNotNone(plotter)
            self.assertEqual(plotter.tabulator, mock_tabulator)
        except ValueError:
            self.fail('Plotter raised ValueError with valid spherical grid type')

    def test_tabulator_cartesian_grid_type_succeeds(self) -> None:
        """Test that cartesian grid type is accepted."""
        mock_tabulator = MockTabulator(has_grid=True, has_gto_data=True, grid_type=GridType.CARTESIAN)

        # This should not raise an exception
        try:
            plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)
            self.assertIsNotNone(plotter)
            self.assertEqual(plotter.tabulator, mock_tabulator)
        except ValueError:
            self.fail('Plotter raised ValueError with valid cartesian grid type')

    def test_valid_tabulator_with_all_attributes_succeeds(self) -> None:
        """Test that a fully valid tabulator is accepted."""
        mock_tabulator = MockTabulator(has_grid=True, has_gto_data=True, grid_type=GridType.SPHERICAL)

        try:
            plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)
            self.assertIsNotNone(plotter)
            self.assertEqual(plotter.tabulator, mock_tabulator)
            # Verify that the provided tabulator is used, not a new one created
            self.assertIs(plotter.tabulator, mock_tabulator)
        except ValueError:
            self.fail('Plotter raised ValueError with fully valid tabulator')

    @patch.object(plotter_module, 'Tabulator')
    def test_none_tabulator_creates_default_tabulator(self, mock_tabulator_class: MockTabulator) -> None:
        """Test that passing None for tabulator creates a default Tabulator."""
        # Create a mock instance for the default tabulator
        mock_tabulator_instance = Mock()
        mock_tabulator_instance._parser.atoms = []
        mock_tabulator_instance._grid_type = GridType.SPHERICAL
        mock_tabulator_instance.grid = np.array([[0, 0, 0]])
        mock_tabulator_instance.grid_dimensions = (1, 1, 1)
        mock_tabulator_class.return_value = mock_tabulator_instance

        plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=None)

        # Verify that a new Tabulator was created
        mock_tabulator_class.assert_called_once_with(str(MOLDEN_PATH), only_molecule=True)
        self.assertEqual(plotter.tabulator, mock_tabulator_instance)

    def test_real_tabulator_instance_validation(self) -> None:
        """Test validation with a real Tabulator instance that has no grid set."""
        # This test uses a real Tabulator instance to verify the validation works
        # with actual objects, not just mocks
        with patch.object(tabulator_module, 'Parser') as mock_parser_class:
            # Mock the parser to avoid needing a real molden file
            mock_parser = Mock()
            mock_parser.atoms = []
            mock_parser_class.return_value = mock_parser

            # Create a real Tabulator instance but don't set up a grid
            real_tabulator = Tabulator(str(MOLDEN_PATH), only_molecule=True)

            # This tabulator should have _grid_type as UNKNOWN and no grid attribute
            self.assertEqual(real_tabulator._grid_type, GridType.UNKNOWN)
            self.assertFalse(hasattr(real_tabulator, 'grid'))

            # Test that it raises the appropriate error
            with self.assertRaises(ValueError) as context:
                Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=real_tabulator)

            self.assertEqual(str(context.exception), 'Tabulator does not have grid attribute.')


class TestPlotterUpdateMesh(unittest.TestCase):
    """Test the update_mesh method of Plotter."""

    def setUp(self) -> None:
        """Set up mocks and a plotter instance."""
        # Mock all GUI-related components
        self.patcher_tk = patch.object(plotter_module, 'tk')
        self.patcher_pv = patch.object(plotter_module, 'pv')
        self.patcher_background_plotter = patch.object(plotter_module, 'BackgroundPlotter')
        self.patcher_molecule = patch.object(plotter_module, 'Molecule')

        self.mock_tk = self.patcher_tk.start()
        self.mock_pv = self.patcher_pv.start()
        self.mock_background_plotter = self.patcher_background_plotter.start()
        self.mock_molecule = self.patcher_molecule.start()

        # Configure mocks
        mock_plotter_instance = Mock()
        mock_plotter_instance.show_axes.return_value = None
        self.mock_background_plotter.return_value = mock_plotter_instance

        mock_molecule_instance = Mock()
        mock_molecule_instance.add_meshes.return_value = []
        mock_molecule_instance.max_radius = 5.0
        self.mock_molecule.return_value = mock_molecule_instance

        # Create a plotter with a mock tabulator
        self.mock_tabulator = MockTabulator(has_grid=True, has_gto_data=True, grid_type=GridType.SPHERICAL)
        self.mock_tabulator.cartesian_grid = Mock()
        self.mock_tabulator.spherical_grid = Mock()

        self.plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=self.mock_tabulator)

    def tearDown(self) -> None:
        """Clean up patches."""
        self.patcher_tk.stop()
        self.patcher_pv.stop()
        self.patcher_background_plotter.stop()
        self.patcher_molecule.stop()

    def test_update_mesh_cartesian_grid(self) -> None:
        """Test update_mesh with cartesian grid type."""
        x_points = np.array([0, 1, 2])
        y_points = np.array([0, 1])
        z_points = np.array([0, 1, 2, 3])

        self.plotter.update_mesh(x_points, y_points, z_points, GridType.CARTESIAN)

        # Verify that cartesian_grid was called with correct parameters
        self.mock_tabulator.cartesian_grid.assert_called_once_with(x_points, y_points, z_points)
        self.mock_tabulator.spherical_grid.assert_not_called()

    def test_update_mesh_spherical_grid(self) -> None:
        """Test update_mesh with spherical grid type."""
        r_points = np.array([0, 1, 2])
        theta_points = np.array([0, np.pi / 2, np.pi])
        phi_points = np.array([0, np.pi, 2 * np.pi])

        self.plotter.update_mesh(r_points, theta_points, phi_points, GridType.SPHERICAL)

        # Verify that spherical_grid was called with correct parameters
        self.mock_tabulator.spherical_grid.assert_called_once_with(r_points, theta_points, phi_points)
        self.mock_tabulator.cartesian_grid.assert_not_called()

    def test_update_mesh_unknown_grid_type_raises_error(self) -> None:
        """Test that update_mesh raises ValueError for UNKNOWN grid type."""
        points = np.array([0, 1, 2])

        with self.assertRaises(ValueError) as context:
            self.plotter.update_mesh(points, points, points, GridType.UNKNOWN)

        self.assertEqual(str(context.exception), 'The plotter only supports spherical and cartesian grids.')

        # Verify that neither grid method was called
        self.mock_tabulator.cartesian_grid.assert_not_called()
        self.mock_tabulator.spherical_grid.assert_not_called()


if __name__ == '__main__':
    unittest.main()
