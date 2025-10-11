"""Unit tests for the configuration save functionality."""

from pathlib import Path

import pytest
import toml

from tests._src_imports import config_module

BACKGROUND_COLOR = 'black'
MO_CONTOUR = 0.2
MO_OPACITY = 0.5
GRID_MIN_RADIUS = 10
BOND_MAX_LENGTH = 5.0


def test_save_current_config_creates_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that save_current_config creates the config file."""
    # Temporarily change the CUSTOM_CONFIG_PATH to a temp directory
    test_config_dir = tmp_path / '.config' / 'moldenViz'
    test_config_dir.mkdir(parents=True, exist_ok=True)
    test_config_path = test_config_dir / 'config.toml'
    monkeypatch.setattr(config_module, 'CUSTOM_CONFIG_PATH', test_config_path)

    # Create a config instance
    config = config_module.Config()

    # Save the config
    config.save_current_config()

    # Verify the file was created
    assert test_config_path.exists()


def test_save_current_config_preserves_values(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that save_current_config correctly saves all configuration values."""
    # Temporarily change the CUSTOM_CONFIG_PATH to a temp directory
    test_config_dir = tmp_path / '.config' / 'moldenViz'
    test_config_dir.mkdir(parents=True, exist_ok=True)
    test_config_path = test_config_dir / 'config.toml'
    monkeypatch.setattr(config_module, 'CUSTOM_CONFIG_PATH', test_config_path)

    # Create a config instance
    config = config_module.Config()

    # Modify some values
    config.config.background_color = BACKGROUND_COLOR
    config.config.mo.contour = MO_CONTOUR
    config.config.mo.opacity = MO_OPACITY
    config.config.grid.min_radius = GRID_MIN_RADIUS
    config.config.molecule.bond.max_length = BOND_MAX_LENGTH

    # Save the config
    config.save_current_config()

    # Load and verify the saved config
    with test_config_path.open('r') as f:
        saved_config = toml.load(f)

    assert saved_config['background_color'] == BACKGROUND_COLOR
    assert saved_config['MO']['contour'] == MO_CONTOUR
    assert saved_config['MO']['opacity'] == MO_OPACITY
    assert saved_config['grid']['min_radius'] == GRID_MIN_RADIUS
    assert saved_config['molecule']['bond']['max_length'] == BOND_MAX_LENGTH


def test_save_current_config_with_custom_colors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that save_current_config correctly saves custom MO colors."""
    # Temporarily change the CUSTOM_CONFIG_PATH to a temp directory
    test_config_dir = tmp_path / '.config' / 'moldenViz'
    test_config_dir.mkdir(parents=True, exist_ok=True)
    test_config_path = test_config_dir / 'config.toml'
    monkeypatch.setattr(config_module, 'CUSTOM_CONFIG_PATH', test_config_path)

    # Create a config instance
    config = config_module.Config()

    # Set custom colors
    config.config.mo.custom_colors = ['blue', 'red']

    # Save the config
    config.save_current_config()

    # Load and verify the saved config
    with test_config_path.open('r') as f:
        saved_config = toml.load(f)

    assert saved_config['MO']['custom_colors'] == ['blue', 'red']


def test_save_current_config_without_custom_colors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that save_current_config doesn't include custom_colors when None."""
    # Temporarily change the CUSTOM_CONFIG_PATH to a temp directory
    test_config_dir = tmp_path / '.config' / 'moldenViz'
    test_config_dir.mkdir(parents=True, exist_ok=True)
    test_config_path = test_config_dir / 'config.toml'
    monkeypatch.setattr(config_module, 'CUSTOM_CONFIG_PATH', test_config_path)

    # Create a config instance
    config = config_module.Config()

    # Ensure custom_colors is None
    config.config.mo.custom_colors = None

    # Save the config
    config.save_current_config()

    # Load and verify the saved config
    with test_config_path.open('r') as f:
        saved_config = toml.load(f)

    assert 'custom_colors' not in saved_config['MO']
