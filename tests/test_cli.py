"""Tests for the command-line interface."""

from __future__ import annotations

import importlib
import logging
import sys
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pytest

from moldenViz import __about__

if TYPE_CHECKING:
    from collections.abc import Generator
    from types import ModuleType


@pytest.fixture(autouse=True)
def reset_root_logger() -> Generator[None, None, None]:
    """Ensure each test starts with a clean logging configuration.

    Yields
    ------
    None
        Allows the test body to run with a pristine logger before cleanup.
    """
    root = logging.getLogger()
    for handler in root.handlers[:]:
        root.removeHandler(handler)
    root.setLevel(logging.WARNING)
    yield
    for handler in root.handlers[:]:
        root.removeHandler(handler)
    root.setLevel(logging.WARNING)


def _reload_cli(monkeypatch: pytest.MonkeyPatch, plotter: Any | None = None) -> ModuleType:
    """Return a freshly reloaded CLI module with an optional Plotter patch.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Fixture used to patch the CLI module during the test.
    plotter : Any | None, optional
        Replacement callable for `Plotter` to capture arguments in tests.

    Returns
    -------
    ModuleType
        Reloaded CLI module ready for invocation.
    """
    cli = importlib.import_module('moldenViz._cli')
    cli = importlib.reload(cli)
    if plotter is not None:
        monkeypatch.setattr(cli, '_resolve_plotter', lambda: plotter)
    return cli


def test_cli_version(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    """Verify ``--version`` prints the package version and exits successfully."""
    cli = _reload_cli(monkeypatch)
    monkeypatch.setattr(sys, 'argv', ['moldenViz', '--version'])

    with pytest.raises(SystemExit) as exc:
        cli.main()

    assert exc.value.code == 0
    captured = capsys.readouterr().out
    assert __about__.__version__ in captured


@pytest.mark.parametrize(
    ('flags', 'expected_level'),
    [
        ([], logging.WARNING),
        (['-v'], logging.INFO),
        (['-d'], logging.DEBUG),
        (['-q'], logging.ERROR),
    ],
)
def test_cli_logging_levels(
    monkeypatch: pytest.MonkeyPatch,
    flags: list[str],
    expected_level: int,
) -> None:
    """Ensure verbosity flags set the expected logging level."""
    calls: dict[str, Any] = {}

    def fake_plotter(source: Any, *, only_molecule: bool = False) -> None:
        calls['source'] = source
        calls['only_molecule'] = only_molecule

    cli = _reload_cli(monkeypatch, fake_plotter)
    sample_file = Path(__file__).with_name('sample_molden.inp')
    monkeypatch.setattr(sys, 'argv', ['moldenViz', *flags, str(sample_file)])

    cli.main()

    assert calls['source'] == str(sample_file)
    assert calls['only_molecule'] is False
    assert logging.getLogger().getEffectiveLevel() == expected_level


def test_cli_example_dispatch(monkeypatch: pytest.MonkeyPatch) -> None:
    """Confirm the ``--example`` flag dispatches bundled Molden data."""
    calls: dict[str, Any] = {}

    def fake_plotter(source: Any, *, only_molecule: bool = False) -> None:
        calls['source'] = source
        calls['only_molecule'] = only_molecule

    cli = _reload_cli(monkeypatch, fake_plotter)
    monkeypatch.setattr(sys, 'argv', ['moldenViz', '--example', 'co'])

    cli.main()

    assert isinstance(calls['source'], list)
    assert len(calls['source']) > 0
    assert calls['only_molecule'] is False
