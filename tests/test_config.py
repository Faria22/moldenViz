"""Unit tests for the configuration module."""

import pytest
from pydantic import ValidationError

from tests._src_imports import config_module


def test_default_color_scheme() -> None:
    """Test that default color scheme is 'bwr'."""
    mo_config = config_module.MOConfig()
    assert mo_config.color_scheme == 'bwr'
    assert mo_config.custom_colors is None


def test_valid_color_scheme() -> None:
    """Test that valid matplotlib colormap names are accepted."""
    valid_cmaps = ['bwr', 'RdBu', 'seismic', 'coolwarm', 'RdYlBu', 'viridis']
    for cmap in valid_cmaps:
        mo_config = config_module.MOConfig(color_scheme=cmap)
        assert mo_config.color_scheme == cmap


def test_invalid_color_scheme_raises_error() -> None:
    """Test that invalid colormap name raises ValidationError."""
    message = 'Color scheme must be a valid matplotlib colormap'
    with pytest.raises(ValidationError, match=message):
        config_module.MOConfig(color_scheme='invalid_colormap_name')


def test_custom_colors_with_two_valid_colors() -> None:
    """Test that custom_colors accepts two valid colors."""
    mo_config = config_module.MOConfig(custom_colors=['blue', 'red'])
    assert mo_config.custom_colors == ['blue', 'red']


def test_custom_colors_with_hex_colors() -> None:
    """Test that custom_colors accepts hex color codes."""
    mo_config = config_module.MOConfig(custom_colors=['#0000FF', '#FF0000'])
    assert mo_config.custom_colors == ['#0000FF', '#FF0000']


def test_custom_colors_with_mixed_formats() -> None:
    """Test that custom_colors accepts mixed color formats."""
    mo_config = config_module.MOConfig(custom_colors=['blue', '#FF0000'])
    assert mo_config.custom_colors == ['blue', '#FF0000']


def test_custom_colors_none_is_allowed() -> None:
    """Test that None is allowed for custom_colors."""
    mo_config = config_module.MOConfig(custom_colors=None)
    assert mo_config.custom_colors is None


def test_custom_colors_with_invalid_color_raises_error() -> None:
    """Test that invalid color in custom_colors raises ValidationError."""
    message = 'Custom color must be a valid matplotlib color'
    with pytest.raises(ValidationError, match=message):
        config_module.MOConfig(custom_colors=['blue', 'not_a_color'])


def test_custom_colors_with_one_color_raises_error() -> None:
    """Test that custom_colors with only one color raises ValidationError."""
    with pytest.raises(ValidationError, match='at least 2 items'):
        config_module.MOConfig(custom_colors=['blue'])


def test_custom_colors_with_three_colors_raises_error() -> None:
    """Test that custom_colors with three colors raises ValidationError."""
    with pytest.raises(ValidationError, match='at most 2 items'):
        config_module.MOConfig(custom_colors=['blue', 'red', 'green'])


def test_custom_colors_with_empty_list_raises_error() -> None:
    """Test that empty list for custom_colors raises ValidationError."""
    with pytest.raises(ValidationError, match='at least 2 items'):
        config_module.MOConfig(custom_colors=[])
