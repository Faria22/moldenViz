"""Tests for Tabulator grid creation using config parameters."""

from pathlib import Path

import pytest

from tests._src_imports import Tabulator

MOLDEN_PATH = Path(__file__).with_name('sample_molden.inp')


def test_cartesian_grid_with_no_params() -> None:
    """Test that cartesian_grid works with no parameters using config defaults."""
    tab = Tabulator(str(MOLDEN_PATH))
    
    # Should not raise an error
    tab.cartesian_grid()
    
    # Grid should be created
    assert tab.grid is not None
    assert tab.grid.shape[1] == 3  # Should have 3 columns (x, y, z)
    assert tab.grid.shape[0] > 0  # Should have at least one point
    
    # GTOs should be tabulated by default
    assert hasattr(tab, '_gtos')


def test_spherical_grid_with_no_params() -> None:
    """Test that spherical_grid works with no parameters using config defaults."""
    tab = Tabulator(str(MOLDEN_PATH))
    
    # Should not raise an error
    tab.spherical_grid()
    
    # Grid should be created
    assert tab.grid is not None
    assert tab.grid.shape[1] == 3  # Should have 3 columns (x, y, z)
    assert tab.grid.shape[0] > 0  # Should have at least one point
    
    # GTOs should be tabulated by default
    assert hasattr(tab, '_gtos')


def test_cartesian_grid_no_gto_tabulation() -> None:
    """Test that cartesian_grid with no params respects tabulate_gtos=False."""
    tab = Tabulator(str(MOLDEN_PATH))
    
    tab.cartesian_grid(tabulate_gtos=False)
    
    # Grid should be created
    assert tab.grid is not None
    
    # GTOs should not be tabulated
    assert not hasattr(tab, '_gtos')


def test_spherical_grid_no_gto_tabulation() -> None:
    """Test that spherical_grid with no params respects tabulate_gtos=False."""
    tab = Tabulator(str(MOLDEN_PATH))
    
    tab.spherical_grid(tabulate_gtos=False)
    
    # Grid should be created
    assert tab.grid is not None
    
    # GTOs should not be tabulated
    assert not hasattr(tab, '_gtos')


def test_cartesian_grid_with_only_molecule() -> None:
    """Test that grid creation with config defaults fails when only_molecule is True."""
    tab = Tabulator(str(MOLDEN_PATH), only_molecule=True)
    
    with pytest.raises(RuntimeError, match='Grid creation is not allowed'):
        tab.cartesian_grid()


def test_spherical_grid_with_only_molecule() -> None:
    """Test that grid creation with config defaults fails when only_molecule is True."""
    tab = Tabulator(str(MOLDEN_PATH), only_molecule=True)
    
    with pytest.raises(RuntimeError, match='Grid creation is not allowed'):
        tab.spherical_grid()
