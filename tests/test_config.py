"""Unit tests for the configuration module."""

import unittest

from pydantic import ValidationError

from tests._src_imports import config_module


class TestMOConfig(unittest.TestCase):
    """Test MOConfig colormap configuration options."""

    def test_default_color_scheme(self) -> None:
        """Test that default color scheme is 'bwr'."""
        mo_config = config_module.MOConfig()
        self.assertEqual(mo_config.color_scheme, 'bwr')
        self.assertIsNone(mo_config.custom_colors)

    def test_valid_color_scheme(self) -> None:
        """Test that valid matplotlib colormap names are accepted."""
        valid_cmaps = ['bwr', 'RdBu', 'seismic', 'coolwarm', 'RdYlBu', 'viridis']
        for cmap in valid_cmaps:
            mo_config = config_module.MOConfig(color_scheme=cmap)
            self.assertEqual(mo_config.color_scheme, cmap)

    def test_invalid_color_scheme_raises_error(self) -> None:
        """Test that invalid colormap name raises ValidationError."""
        with self.assertRaises(ValidationError) as context:
            config_module.MOConfig(color_scheme='invalid_colormap_name')
        
        self.assertIn('Color scheme must be a valid matplotlib colormap', str(context.exception))

    def test_custom_colors_with_two_valid_colors(self) -> None:
        """Test that custom_colors accepts two valid colors."""
        mo_config = config_module.MOConfig(custom_colors=['blue', 'red'])
        self.assertEqual(mo_config.custom_colors, ['blue', 'red'])

    def test_custom_colors_with_hex_colors(self) -> None:
        """Test that custom_colors accepts hex color codes."""
        mo_config = config_module.MOConfig(custom_colors=['#0000FF', '#FF0000'])
        self.assertEqual(mo_config.custom_colors, ['#0000FF', '#FF0000'])

    def test_custom_colors_with_mixed_formats(self) -> None:
        """Test that custom_colors accepts mixed color formats."""
        mo_config = config_module.MOConfig(custom_colors=['blue', '#FF0000'])
        self.assertEqual(mo_config.custom_colors, ['blue', '#FF0000'])

    def test_custom_colors_none_is_allowed(self) -> None:
        """Test that None is allowed for custom_colors."""
        mo_config = config_module.MOConfig(custom_colors=None)
        self.assertIsNone(mo_config.custom_colors)

    def test_custom_colors_with_invalid_color_raises_error(self) -> None:
        """Test that invalid color in custom_colors raises ValidationError."""
        with self.assertRaises(ValidationError) as context:
            config_module.MOConfig(custom_colors=['blue', 'not_a_color'])
        
        self.assertIn('Custom color must be a valid matplotlib color', str(context.exception))

    def test_custom_colors_with_one_color_raises_error(self) -> None:
        """Test that custom_colors with only one color raises ValidationError."""
        with self.assertRaises(ValidationError) as context:
            config_module.MOConfig(custom_colors=['blue'])
        
        self.assertIn('at least 2 items', str(context.exception))

    def test_custom_colors_with_three_colors_raises_error(self) -> None:
        """Test that custom_colors with three colors raises ValidationError."""
        with self.assertRaises(ValidationError) as context:
            config_module.MOConfig(custom_colors=['blue', 'red', 'green'])
        
        self.assertIn('at most 2 items', str(context.exception))

    def test_custom_colors_with_empty_list_raises_error(self) -> None:
        """Test that empty list for custom_colors raises ValidationError."""
        with self.assertRaises(ValidationError) as context:
            config_module.MOConfig(custom_colors=[])
        
        self.assertIn('at least 2 items', str(context.exception))


if __name__ == '__main__':
    unittest.main()
