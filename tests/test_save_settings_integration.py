"""Integration test to verify Save Settings button functionality in plotter."""

from typing import Any
from unittest.mock import MagicMock, patch

from tests._src_imports import plotter_module


def test_save_settings_method_exists_in_orbital_selection_screen() -> None:
    """Test that save_settings method exists in _OrbitalSelectionScreen."""
    # Access the private class through the module
    orbital_selection_screen = plotter_module._OrbitalSelectionScreen  # noqa: SLF001

    # Check that the method exists
    assert hasattr(orbital_selection_screen, 'save_settings')

    # Verify it's a callable method
    assert callable(orbital_selection_screen.save_settings)


@patch('moldenViz.plotter.messagebox')
@patch('moldenViz.plotter.config')
def test_save_settings_success(mock_config: Any, mock_messagebox: Any) -> None:
    """Test that save_settings calls config.save_current_config and shows success message."""
    # Access the private class through the module
    orbital_selection_screen_class = plotter_module._OrbitalSelectionScreen  # noqa: SLF001

    # Set up mocks
    mock_config.save_current_config = MagicMock()

    # Create a minimal instance (without full initialization)
    # We only need to test the save_settings method
    selection_screen = orbital_selection_screen_class.__new__(orbital_selection_screen_class)

    # Call the method
    selection_screen.save_settings()

    # Verify that save_current_config was called
    mock_config.save_current_config.assert_called_once()

    # Verify that success message was shown
    mock_messagebox.showinfo.assert_called_once()
    args = mock_messagebox.showinfo.call_args[0]
    assert 'Settings Saved' in args[0]
    assert 'Configuration saved successfully' in args[1]


@patch('moldenViz.plotter.messagebox')
@patch('moldenViz.plotter.config')
def test_save_settings_handles_oserror(mock_config: Any, mock_messagebox: Any) -> None:
    """Test that save_settings handles OSError gracefully."""
    # Access the private class through the module
    orbital_selection_screen_class = plotter_module._OrbitalSelectionScreen  # noqa: SLF001

    # Set up mock to raise OSError
    mock_config.save_current_config = MagicMock(side_effect=OSError('Permission denied'))

    # Create a minimal instance
    selection_screen = orbital_selection_screen_class.__new__(orbital_selection_screen_class)

    # Call the method
    selection_screen.save_settings()

    # Verify that error message was shown
    mock_messagebox.showerror.assert_called_once()
    args = mock_messagebox.showerror.call_args[0]
    assert 'Save Error' in args[0]
    assert 'Failed to save configuration' in args[1]
    assert 'Permission denied' in args[1]


@patch('moldenViz.plotter.messagebox')
@patch('moldenViz.plotter.config')
def test_save_settings_handles_valueerror(mock_config: Any, mock_messagebox: Any) -> None:
    """Test that save_settings handles ValueError gracefully."""
    # Access the private class through the module
    orbital_selection_screen_class = plotter_module._OrbitalSelectionScreen  # noqa: SLF001

    # Set up mock to raise ValueError
    mock_config.save_current_config = MagicMock(side_effect=ValueError('Invalid config'))

    # Create a minimal instance
    selection_screen = orbital_selection_screen_class.__new__(orbital_selection_screen_class)

    # Call the method
    selection_screen.save_settings()

    # Verify that error message was shown
    mock_messagebox.showerror.assert_called_once()
    args = mock_messagebox.showerror.call_args[0]
    assert 'Save Error' in args[0]
    assert 'Failed to save configuration' in args[1]
    assert 'Invalid config' in args[1]
