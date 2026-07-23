"""Background job orchestration used by the interactive plotter."""

from __future__ import annotations

from collections.abc import Callable
from threading import Lock
from time import perf_counter
from typing import TYPE_CHECKING, Generic, TypeVar

if TYPE_CHECKING:
    from concurrent.futures import Executor, Future

_ResultT = TypeVar('_ResultT')
_Dispatch = Callable[[Callable[[], None]], None]


class BackgroundJob(Generic[_ResultT]):
    """Run one replaceable background job and deliver its result.

    The controller deliberately has no dependency on Tk, Qt, or PyVista. The
    caller supplies a dispatcher that transfers completion work to its owning
    event loop.

    Parameters
    ----------
    executor : Executor
        Executor used to run the background callable.
    dispatch : Callable[[Callable[[], None]], None]
        Function that schedules a no-argument completion callback.
    """

    def __init__(self, executor: Executor, dispatch: _Dispatch) -> None:
        self._executor = executor
        self._dispatch = dispatch
        self._lock = Lock()
        self._generation = 0
        self._future: Future[_ResultT] | None = None
        self._started_at: float | None = None
        self._on_success: Callable[[_ResultT, float], None] | None = None
        self._on_error: Callable[[Exception], None] | None = None

    @property
    def future(self) -> Future[_ResultT] | None:
        """The active future, if a job is pending."""
        with self._lock:
            return self._future

    @property
    def pending(self) -> bool:
        """Whether a job is currently active."""
        return self.future is not None

    def start(
        self,
        work: Callable[[], _ResultT],
        *,
        on_success: Callable[[_ResultT, float], None],
        on_error: Callable[[Exception], None],
    ) -> bool:
        """Start ``work`` unless another job is already active.

        Returns
        -------
        bool
            ``True`` when a new job was submitted, otherwise ``False``.
        """
        with self._lock:
            if self._future is not None:
                return False
            self._generation += 1
            generation = self._generation
            self._started_at = perf_counter()
            self._on_success = on_success
            self._on_error = on_error
            future = self._executor.submit(work)
            self._future = future

        future.add_done_callback(
            lambda completed, generation=generation: self._dispatch(
                lambda: self._finish(generation, completed),
            ),
        )
        return True

    def wait(self, timeout: float | None = None) -> _ResultT:
        """Wait for the active job and deliver its result exactly once.

        Parameters
        ----------
        timeout : float | None, optional
            Maximum number of seconds to wait.

        Returns
        -------
        _ResultT
            Result returned by the background callable.

        Raises
        ------
        RuntimeError
            If no job is active.
        """
        with self._lock:
            future = self._future
            generation = self._generation
        if future is None:
            raise RuntimeError('Background job has not been scheduled.')

        result = future.result(timeout=timeout)
        self._finish(generation, future)
        return result

    def cancel(self) -> None:
        """Invalidate and cancel the active job when possible."""
        with self._lock:
            future = self._future
            self._generation += 1
            self._clear_locked()
        if future is not None and not future.done():
            future.cancel()

    def _finish(self, generation: int, future: Future[_ResultT]) -> None:
        with self._lock:
            if generation != self._generation or future is not self._future:
                return

            started_at = self._started_at
            on_success = self._on_success
            on_error = self._on_error
            self._clear_locked()

        if future.cancelled():
            return
        try:
            result = future.result()
        except Exception as exc:  # ruff:ignore[blind-except]
            if on_error is not None:
                on_error(exc)
            return

        elapsed = perf_counter() - started_at if started_at is not None else 0.0
        if on_success is not None:
            on_success(result, elapsed)

    def _clear_locked(self) -> None:
        self._future = None
        self._started_at = None
        self._on_success = None
        self._on_error = None
