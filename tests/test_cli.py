"""Tests for the command line interface."""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Iterator

import pytest

from moldenViz.__about__ import __version__
from moldenViz._cli import main
from moldenViz.examples.get_example_files import all_examples

SAMPLE_FILE = str(Path(__file__).with_name('sample_molden.inp'))


@pytest.fixture
def plotter_spy(monkeypatch: pytest.MonkeyPatch) -> list[tuple[object, bool]]:
    """Collect invocations of :class:`moldenViz.plotter.Plotter`."""

    calls: list[tuple[object, bool]] = []

    def _fake_plotter(path: object, *, only_molecule: bool = False) -> None:
        calls.append((path, only_molecule))

    monkeypatch.setattr('moldenViz._cli.Plotter', _fake_plotter)
    return calls


@pytest.fixture
def fresh_logging() -> Iterator[None]:
    """Ensure ``logging.basicConfig`` reconfigures the root logger."""

    original_handlers = logging.root.handlers[:]
    original_level = logging.getLogger().level
    logging.root.handlers = []
    yield
    logging.root.handlers = original_handlers
    logging.getLogger().setLevel(original_level)


def test_cli_version_flag(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    plotter_spy: list[tuple[object, bool]],
    fresh_logging: Iterator[None],
) -> None:
    """``--version`` prints the package version and exits successfully."""

    monkeypatch.setattr(sys, 'argv', ['moldenViz', '--version'])

    with pytest.raises(SystemExit) as exc:
        main()

    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out
    assert plotter_spy == []


def test_cli_invokes_plotter_with_file(
    monkeypatch: pytest.MonkeyPatch, plotter_spy: list[tuple[object, bool]], fresh_logging: Iterator[None]
) -> None:
    """Providing a file path dispatches to :class:`Plotter`."""

    monkeypatch.setattr(sys, 'argv', ['moldenViz', SAMPLE_FILE])

    main()

    assert plotter_spy == [(SAMPLE_FILE, False)]
    assert logging.getLogger().getEffectiveLevel() == logging.WARNING


def test_cli_invokes_plotter_with_example(
    monkeypatch: pytest.MonkeyPatch, plotter_spy: list[tuple[object, bool]], fresh_logging: Iterator[None]
) -> None:
    """``--example`` loads the requested bundled molecule."""

    monkeypatch.setattr(sys, 'argv', ['moldenViz', '-e', 'benzene'])

    main()

    assert plotter_spy == [(all_examples['benzene'], False)]


@pytest.mark.parametrize(
    ('argv', 'expected_level'),
    [
        (['moldenViz', '-v', SAMPLE_FILE], logging.INFO),
        (['moldenViz', '-d', SAMPLE_FILE], logging.DEBUG),
        (['moldenViz', '-q', SAMPLE_FILE], logging.ERROR),
    ],
)
def test_cli_sets_requested_logging_levels(
    monkeypatch: pytest.MonkeyPatch,
    plotter_spy: list[tuple[object, bool]],
    fresh_logging: Iterator[None],
    argv: list[str],
    expected_level: int,
) -> None:
    """Verbosity flags adjust the root logger level."""

    monkeypatch.setattr(sys, 'argv', argv)

    main()

    assert plotter_spy == [(SAMPLE_FILE, False)]
    assert logging.getLogger().getEffectiveLevel() == expected_level
