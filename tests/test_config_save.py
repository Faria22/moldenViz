"""Unit tests for the configuration save functionality."""

from pathlib import Path

import pytest
import toml

from tests._src_imports import config_module


def test_save_current_config_creates_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that save_current_config creates the config file."""
    # Temporarily change the custom_configs_dir to a temp directory
    test_config_dir = tmp_path / '.config' / 'moldenViz'
    test_config_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(config_module, 'custom_configs_dir', test_config_dir)

    # Create a config instance
    config = config_module.Config()

    # Save the config
    config.save_current_config()

    # Verify the file was created
    config_file = test_config_dir / 'config.toml'
    assert config_file.exists()


def test_save_current_config_preserves_values(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that save_current_config correctly saves all configuration values."""
    # Temporarily change the custom_configs_dir to a temp directory
    test_config_dir = tmp_path / '.config' / 'moldenViz'
    test_config_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(config_module, 'custom_configs_dir', test_config_dir)

    # Create a config instance
    config = config_module.Config()

    # Modify some values
    config.config.background_color = 'black'
    config.config.mo.contour = 0.2
    config.config.mo.opacity = 0.5
    config.config.grid.min_radius = 10
    config.config.molecule.bond.max_length = 5.0

    # Save the config
    config.save_current_config()

    # Load and verify the saved config
    config_file = test_config_dir / 'config.toml'
    with config_file.open('r') as f:
        saved_config = toml.load(f)

    assert saved_config['background_color'] == 'black'
    assert saved_config['MO']['contour'] == 0.2
    assert saved_config['MO']['opacity'] == 0.5
    assert saved_config['grid']['min_radius'] == 10
    assert saved_config['molecule']['bond']['max_length'] == 5.0


def test_save_current_config_with_custom_colors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that save_current_config correctly saves custom MO colors."""
    # Temporarily change the custom_configs_dir to a temp directory
    test_config_dir = tmp_path / '.config' / 'moldenViz'
    test_config_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(config_module, 'custom_configs_dir', test_config_dir)

    # Create a config instance
    config = config_module.Config()

    # Set custom colors
    config.config.mo.custom_colors = ['blue', 'red']

    # Save the config
    config.save_current_config()

    # Load and verify the saved config
    config_file = test_config_dir / 'config.toml'
    with config_file.open('r') as f:
        saved_config = toml.load(f)

    assert saved_config['MO']['custom_colors'] == ['blue', 'red']


def test_save_current_config_without_custom_colors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that save_current_config doesn't include custom_colors when None."""
    # Temporarily change the custom_configs_dir to a temp directory
    test_config_dir = tmp_path / '.config' / 'moldenViz'
    test_config_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.setattr(config_module, 'custom_configs_dir', test_config_dir)

    # Create a config instance
    config = config_module.Config()

    # Ensure custom_colors is None
    config.config.mo.custom_colors = None

    # Save the config
    config.save_current_config()

    # Load and verify the saved config
    config_file = test_config_dir / 'config.toml'
    with config_file.open('r') as f:
        saved_config = toml.load(f)

    assert 'custom_colors' not in saved_config['MO']
