"""Tests for the ColoredFormatter."""

import logging

from moldenViz._cli import ColoredFormatter  # noqa: PLC2701


def test_colored_formatter_has_correct_colors() -> None:
    """Verify that ColoredFormatter has the expected color codes."""
    assert ColoredFormatter.COLORS['ERROR'] == '\033[38;5;196m'
    assert ColoredFormatter.COLORS['WARNING'] == '\033[38;5;208m'
    assert ColoredFormatter.COLORS['INFO'] == '\033[38;5;34m'
    assert ColoredFormatter.COLORS['DEBUG'] == '\033[38;5;27m'
    assert ColoredFormatter.COLORS['RESET'] == '\033[0m'


def test_colored_formatter_applies_color_to_debug() -> None:
    """Verify that DEBUG messages get blue color."""
    formatter = ColoredFormatter('%(levelname)s: %(message)s')
    record = logging.LogRecord(
        name='test',
        level=logging.DEBUG,
        pathname='',
        lineno=0,
        msg='test message',
        args=(),
        exc_info=None,
    )
    formatted = formatter.format(record)
    assert '\033[38;5;27m' in formatted  # Blue color for DEBUG
    assert '\033[0m' in formatted  # Reset color
    assert 'DEBUG' in formatted
    assert 'test message' in formatted


def test_colored_formatter_applies_color_to_info() -> None:
    """Verify that INFO messages get green color."""
    formatter = ColoredFormatter('%(levelname)s: %(message)s')
    record = logging.LogRecord(
        name='test',
        level=logging.INFO,
        pathname='',
        lineno=0,
        msg='test message',
        args=(),
        exc_info=None,
    )
    formatted = formatter.format(record)
    assert '\033[38;5;34m' in formatted  # Green color for INFO
    assert '\033[0m' in formatted  # Reset color
    assert 'INFO' in formatted
    assert 'test message' in formatted


def test_colored_formatter_applies_color_to_warning() -> None:
    """Verify that WARNING messages get orange color."""
    formatter = ColoredFormatter('%(levelname)s: %(message)s')
    record = logging.LogRecord(
        name='test',
        level=logging.WARNING,
        pathname='',
        lineno=0,
        msg='test message',
        args=(),
        exc_info=None,
    )
    formatted = formatter.format(record)
    assert '\033[38;5;208m' in formatted  # Orange color for WARNING
    assert '\033[0m' in formatted  # Reset color
    assert 'WARNING' in formatted
    assert 'test message' in formatted


def test_colored_formatter_applies_color_to_error() -> None:
    """Verify that ERROR messages get red color."""
    formatter = ColoredFormatter('%(levelname)s: %(message)s')
    record = logging.LogRecord(
        name='test',
        level=logging.ERROR,
        pathname='',
        lineno=0,
        msg='test message',
        args=(),
        exc_info=None,
    )
    formatted = formatter.format(record)
    assert '\033[38;5;196m' in formatted  # Red color for ERROR
    assert '\033[0m' in formatted  # Reset color
    assert 'ERROR' in formatted
    assert 'test message' in formatted


def test_colored_formatter_restores_original_levelname() -> None:
    """Verify that the formatter restores the original levelname after formatting."""
    formatter = ColoredFormatter('%(levelname)s: %(message)s')
    record = logging.LogRecord(
        name='test',
        level=logging.INFO,
        pathname='',
        lineno=0,
        msg='test message',
        args=(),
        exc_info=None,
    )
    original_levelname = record.levelname
    formatter.format(record)
    assert record.levelname == original_levelname
