"""Integration test to verify Save Settings button functionality in plotter."""

from typing import Any
from unittest.mock import MagicMock, patch

import pytest


class _ModuleProxy:
    def __init__(self) -> None:
        self._module: Any | None = None

    def _load(self) -> Any:
        if self._module is None:
            self._module = pytest.importorskip('moldenViz.plotter')
        return self._module

    def __getattr__(self, item: str) -> Any:
        return getattr(self._load(), item)


plotter_module = _ModuleProxy()


def test_save_settings_method_exists_in_plotter() -> None:
    """Test that save_settings method exists on Plotter."""
    plotter_class = plotter_module.Plotter

    assert hasattr(plotter_class, 'save_settings')
    assert callable(plotter_class.save_settings)


@patch('moldenViz.plotter.messagebox')
@patch('moldenViz.plotter.config')
def test_save_settings_success(mock_config: Any, mock_messagebox: Any) -> None:
    """Test that save_settings calls config.save_current_config and shows success message."""
    # Set up mocks
    mock_config.save_current_config = MagicMock()

    # Call the method
    plotter_module.Plotter.save_settings()

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
    # Set up mock to raise OSError
    mock_config.save_current_config = MagicMock(side_effect=OSError('Permission denied'))

    # Call the method
    plotter_module.Plotter.save_settings()

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
    # Set up mock to raise ValueError
    mock_config.save_current_config = MagicMock(side_effect=ValueError('Invalid config'))

    # Call the method
    plotter_module.Plotter.save_settings()

    # Verify that error message was shown
    mock_messagebox.showerror.assert_called_once()
    args = mock_messagebox.showerror.call_args[0]
    assert 'Save Error' in args[0]
    assert 'Failed to save configuration' in args[1]
    assert 'Invalid config' in args[1]
