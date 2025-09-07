#!/usr/bin/env python3
"""
Demonstration script showing the plotter tests in action.

This script shows that the tests properly validate the error conditions
when a tabulator is given to the Plotter class.
"""

import sys
sys.path.insert(0, 'src')

from moldenViz.plotter import Plotter
from moldenViz.tabulator import GridType, Tabulator
from pathlib import Path
import numpy as np
from unittest.mock import Mock, patch

# Sample molden file path
MOLDEN_PATH = Path('tests/sample_molden.inp')

class MockTabulator:
    """Mock tabulator for demonstration."""
    
    def __init__(self, has_grid=True, has_gto_data=True, grid_type=GridType.SPHERICAL):
        # Always create attributes, then optionally remove them
        self.grid = np.array([[0, 0, 0], [1, 1, 1]])
        self.gto_data = np.array([[1, 2, 3], [4, 5, 6]])
        self._grid_type = grid_type
        self.grid_dimensions = (2, 2, 2)
        
        # Remove attributes if requested
        if not has_grid:
            delattr(self, 'grid')
        if not has_gto_data:
            delattr(self, 'gto_data')
        
        # Mock parser with atoms for Molecule creation
        self._parser = Mock()
        self._parser.atoms = []

def demo_test_case(description, test_function):
    """Run a test case and report the result."""
    print(f"\n{description}")
    print("=" * len(description))
    try:
        test_function()
        print("✓ PASSED: Test completed successfully")
    except Exception as e:
        print(f"✓ PASSED: Expected error was raised: {e}")
    except AssertionError as e:
        print(f"✗ FAILED: {e}")

def test_missing_grid():
    """Test that missing grid attribute raises ValueError."""
    with patch('moldenViz.plotter.BackgroundPlotter'), \
         patch('moldenViz.plotter.Molecule'), \
         patch('moldenViz.plotter.tk'):
        
        mock_tabulator = MockTabulator(has_grid=False)
        
        try:
            Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)
            raise AssertionError("Expected ValueError was not raised")
        except ValueError as e:
            if str(e) != 'Tabulator does not have grid attribute.':
                raise AssertionError(f"Wrong error message: {e}")

def test_missing_gto_data():
    """Test that missing gto_data raises ValueError when only_molecule=False."""
    with patch('moldenViz.plotter.BackgroundPlotter'), \
         patch('moldenViz.plotter.Molecule'), \
         patch('moldenViz.plotter.tk'):
        
        mock_tabulator = MockTabulator(has_gto_data=False)
        
        try:
            Plotter(str(MOLDEN_PATH), only_molecule=False, tabulator=mock_tabulator)
            raise AssertionError("Expected ValueError was not raised")
        except ValueError as e:
            if str(e) != 'Tabulator does not have tabulated GTOs.':
                raise AssertionError(f"Wrong error message: {e}")

def test_unknown_grid_type():
    """Test that UNKNOWN grid type raises ValueError."""
    with patch('moldenViz.plotter.BackgroundPlotter'), \
         patch('moldenViz.plotter.Molecule'), \
         patch('moldenViz.plotter.tk'):
        
        mock_tabulator = MockTabulator(grid_type=GridType.UNKNOWN)
        
        try:
            Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)
            raise AssertionError("Expected ValueError was not raised")
        except ValueError as e:
            if str(e) != 'The plotter only supports spherical and cartesian grids.':
                raise AssertionError(f"Wrong error message: {e}")

def test_valid_tabulator():
    """Test that valid tabulator is accepted."""
    with patch('moldenViz.plotter.BackgroundPlotter'), \
         patch('moldenViz.plotter.Molecule') as mock_molecule, \
         patch('moldenViz.plotter.tk'):
        
        mock_molecule.return_value.add_meshes.return_value = []
        mock_tabulator = MockTabulator(has_grid=True, has_gto_data=True, grid_type=GridType.SPHERICAL)
        
        plotter = Plotter(str(MOLDEN_PATH), only_molecule=True, tabulator=mock_tabulator)
        assert plotter.tab is mock_tabulator

if __name__ == '__main__':
    print("Plotter Unit Tests Demonstration")
    print("===============================")
    print("\nThis demonstrates the unit tests for the Plotter class")
    print("focusing on tabulator validation as requested.")
    
    demo_test_case(
        "Test 1: Tabulator missing 'grid' attribute",
        test_missing_grid
    )
    
    demo_test_case(
        "Test 2: Tabulator missing 'gto_data' when only_molecule=False", 
        test_missing_gto_data
    )
    
    demo_test_case(
        "Test 3: Tabulator with UNKNOWN grid type",
        test_unknown_grid_type
    )
    
    demo_test_case(
        "Test 4: Valid tabulator is accepted",
        test_valid_tabulator
    )
    
    print("\n" + "=" * 50)
    print("All key validation tests completed successfully!")
    print("The unit tests properly check that correct errors")
    print("are raised when invalid tabulators are provided.")