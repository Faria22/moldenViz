"""Tests for GUI-independent Plotter job orchestration."""
# ruff:file-ignore[import-private-name, magic-value-comparison, no-self-use, undocumented-public-function, undocumented-public-method]

from __future__ import annotations

from concurrent.futures import Executor, Future
from typing import Any

import pytest

from moldenViz._plotter_jobs import BackgroundJob


class InlineExecutor(Executor):
    """Executor that completes submitted work deterministically."""

    def submit(self, fn: Any, /, *args: Any, **kwargs: Any) -> Future[Any]:
        future: Future[Any] = Future()
        try:
            future.set_result(fn(*args, **kwargs))
        except Exception as exc:  # ruff:ignore[blind-except]
            future.set_exception(exc)
        return future


def test_background_job_delivers_success_through_supplied_dispatcher() -> None:
    callbacks: list[Any] = []
    successes: list[tuple[int, float]] = []
    errors: list[Exception] = []
    job = BackgroundJob[int](InlineExecutor(), callbacks.append)

    assert job.start(
        lambda: 42,
        on_success=lambda result, elapsed: successes.append((result, elapsed)),
        on_error=errors.append,
    )
    assert job.pending
    assert successes == []

    callbacks.pop()()

    assert successes[0][0] == 42
    assert successes[0][1] >= 0.0
    assert errors == []
    assert not job.pending


def test_background_job_wait_delivers_result_only_once() -> None:
    callbacks: list[Any] = []
    successes: list[int] = []
    job = BackgroundJob[int](InlineExecutor(), callbacks.append)
    job.start(
        lambda: 7,
        on_success=lambda result, _elapsed: successes.append(result),
        on_error=lambda exc: pytest.fail(f'unexpected error: {exc}'),
    )

    assert job.wait() == 7
    callbacks.pop()()

    assert successes == [7]
    assert not job.pending


def test_background_job_ignores_cancelled_generation() -> None:
    callbacks: list[Any] = []
    successes: list[int] = []
    job = BackgroundJob[int](InlineExecutor(), callbacks.append)
    job.start(
        lambda: 1,
        on_success=lambda result, _elapsed: successes.append(result),
        on_error=lambda exc: pytest.fail(f'unexpected error: {exc}'),
    )

    job.cancel()
    callbacks.pop()()

    assert successes == []
    assert not job.pending


def test_background_job_reports_worker_exception() -> None:
    callbacks: list[Any] = []
    errors: list[Exception] = []
    job = BackgroundJob[int](InlineExecutor(), callbacks.append)

    def fail() -> int:
        raise ValueError('bad grid')

    job.start(
        fail,
        on_success=lambda _result, _elapsed: pytest.fail('success callback ran'),
        on_error=errors.append,
    )
    callbacks.pop()()

    assert len(errors) == 1
    assert isinstance(errors[0], ValueError)
    assert str(errors[0]) == 'bad grid'
