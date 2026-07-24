"""Tests covering asynchronous GTO tabulation behavior in Plotter."""
# ruff:file-ignore[private-member-access]

from __future__ import annotations

import threading
import time
from collections import UserDict
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any, cast

import numpy as np

from moldenViz import plotter as plotter_module

if TYPE_CHECKING:
    from collections.abc import Callable

    import pytest

_REAL_SELECTION_SCREEN = plotter_module._OrbitalSelectionScreen


class DummyMesh(UserDict):
    """Stub StructuredGrid used in place of the PyVista mesh."""

    def contour(self, *_args: object, **_kwargs: object) -> DummyMesh:
        """Return the same mesh instance for fluent-style chaining.

        Returns
        -------
        DummyMesh
            Self so method chains continue to work in tests.
        """
        return self


class FakeSignal:
    """Record callbacks connected to a Qt-style signal."""

    def __init__(self) -> None:
        self.callbacks: list[Callable[[], None]] = []

    def connect(self, callback: Callable[[], None]) -> None:
        """Record a callback for later deterministic delivery."""
        self.callbacks.append(callback)

    def emit(self) -> None:
        """Deliver the signal to every connected callback."""
        for callback in self.callbacks:
            callback()


class FakeBackgroundPlotter:
    """Lightweight BackgroundPlotter stand-in for headless tests."""

    def __init__(self, *_args: object, **_kwargs: object) -> None:
        self.app_window = SimpleNamespace(signal_close=FakeSignal())
        self.main_menu = SimpleNamespace(actions=list, addMenu=lambda *_args, **_kwargs: None)
        self._last_mesh_args: tuple[tuple[object, ...], dict[str, object]] | None = None
        self.closed = False

    def set_background(self, *_args: object, **_kwargs: object) -> None:  # pragma: no cover - trivial
        """Ignore background updates for headless runs."""

    def show_axes(self) -> None:  # pragma: no cover - trivial
        """Pretend to show axes to match the PyVista API."""

    def add_mesh(self, *_args: object, **_kwargs: object) -> object:  # pragma: no cover - trivial
        """Record mesh creation requests and return a dummy actor.

        Returns
        -------
        object
            Stand-in actor handle resembling PyVista's return value.
        """
        self._last_mesh_args = (_args, _kwargs)
        return object()

    def remove_actor(self, *_args: object, **_kwargs: object) -> None:  # pragma: no cover - trivial
        """Ignore mesh removal events."""

    def close(self) -> None:  # pragma: no cover - trivial
        """Simulate closing the plotter window."""
        self.closed = True

    def update(self) -> None:  # pragma: no cover - trivial
        """Expose the update hook expected by Plotter."""


class FakeTk:
    """Minimal Tk root substitute with an explicitly pumped event loop."""

    def __init__(self) -> None:
        self.owner_thread_id = threading.get_ident()
        self.tk_call_thread_ids: list[int] = []
        self.callbacks: dict[str, tuple[object, tuple[object, ...]]] = {}
        self._next_callback_id = 0
        self.quit_calls = 0
        self.destroy_calls = 0

    def after(self, _delay: int, callback: object, *args: object) -> str:
        """Queue callbacks for explicit execution by the owning test thread.

        Returns
        -------
        str
            Identifier that can be passed to :meth:`after_cancel`.
        """
        self.tk_call_thread_ids.append(threading.get_ident())
        assert self.tk_call_thread_ids[-1] == self.owner_thread_id
        assert callable(callback)
        callback_id = f'after-{self._next_callback_id}'
        self._next_callback_id += 1
        self.callbacks[callback_id] = (callback, args)
        return callback_id

    def after_cancel(self, callback_id: str) -> None:
        """Remove a queued callback."""
        self.tk_call_thread_ids.append(threading.get_ident())
        assert self.tk_call_thread_ids[-1] == self.owner_thread_id
        self.callbacks.pop(callback_id, None)

    def run_one_callback(self) -> bool:
        """Run one queued callback as the owning Tk thread.

        Returns
        -------
        bool
            Whether a callback was available and executed.
        """
        if not self.callbacks:
            return False
        callback_id = next(iter(self.callbacks))
        callback, args = self.callbacks.pop(callback_id)
        assert callable(callback)
        callback(*args)
        return True

    def mainloop(self) -> None:  # pragma: no cover - trivial
        """Expose Tk's mainloop hook for compatibility."""

    def quit(self) -> None:  # pragma: no cover - trivial
        """Satisfy the Tk API contract for quitting."""
        self.quit_calls += 1

    def destroy(self) -> None:  # pragma: no cover - trivial
        """Simulate tearing down the Tk root."""
        self.destroy_calls += 1


class FakeSelectionScreen:
    """Stub orbital selection dialog used to capture loading indicator states."""

    def __init__(self, plotter: plotter_module.Plotter) -> None:
        self.plotter = plotter
        self.current_mo_ind = -1
        self.loading_states: list[bool] = []
        self.messages: list[str] = []
        self.destroyed = False

    def _set_loading_state(self, loading: bool, message: str = 'Tabulating orbitals...') -> None:
        """Mirror the dialog's loading indicator for assertions."""
        self.loading_states.append(loading)
        self.messages.append(message)

    def _on_gtos_ready(self) -> None:
        """Transition the dialog back to idle once GTOs arrive."""
        self._set_loading_state(False)

    def _update_nav_button_states(self) -> None:  # pragma: no cover - trivial
        """Stub navigation button updates."""

    def winfo_exists(self) -> bool:  # pragma: no cover - trivial
        """Report that the dialog window is still alive.

        Returns
        -------
        bool
            True while the fake selection dialog is attached to a plotter.
        """
        return not self.destroyed

    def destroy(self) -> None:  # pragma: no cover - trivial
        """Implement the Tk destroy hook."""
        self.destroyed = True


def _configure_plotter_env(monkeypatch: pytest.MonkeyPatch) -> None:
    """Patch heavyweight UI helpers so Plotter can be instantiated in tests."""

    def fake_load_molecule(self: plotter_module.Plotter, _config: object) -> None:
        dummy = SimpleNamespace(max_radius=1.0, atoms=[])
        self._molecule = cast(Any, dummy)
        self._molecule_actors = []
        self._atom_actors = []
        self._bond_actors = []

    monkeypatch.setattr(plotter_module, 'BackgroundPlotter', FakeBackgroundPlotter)
    monkeypatch.setattr(plotter_module, '_OrbitalSelectionScreen', FakeSelectionScreen)
    monkeypatch.setattr(plotter_module.Plotter, '_add_orbital_menus_to_pv_plotter', lambda _plotter: None)
    monkeypatch.setattr(plotter_module.Plotter, '_connect_pv_plotter_close_signal', lambda _plotter: None)
    monkeypatch.setattr(plotter_module.Plotter, '_override_clear_all_button', lambda _plotter: None)
    monkeypatch.setattr(plotter_module.Plotter, '_load_molecule', fake_load_molecule)
    monkeypatch.setattr(plotter_module.Plotter, '_create_mo_mesh', lambda _plotter: DummyMesh())


def _sample_molden() -> str:
    return str(Path(__file__).with_name('sample_molden.inp'))


def _fake_tk_root() -> Any:
    """Return a FakeTk instance typed loosely for Plotter construction.

    Returns
    -------
    Any
        A FakeTk handle treated permissively for type checking.
    """
    return cast(Any, FakeTk())


def _pump_until(root: FakeTk, condition: Callable[[], bool], timeout: float = 1.0) -> None:
    """Pump fake Tk callbacks until ``condition`` becomes true."""
    deadline = time.monotonic() + timeout
    while not condition():
        assert time.monotonic() < deadline
        root.run_one_callback()
        time.sleep(0.001)


def test_plotter_defers_gto_tabulation(monkeypatch: pytest.MonkeyPatch) -> None:
    """Default Plotter grids should tabulate GTOs in the background."""
    _configure_plotter_env(monkeypatch)
    start_event = threading.Event()
    finish_event = threading.Event()

    def fake_compute_gtos(_tabulator: plotter_module.Tabulator, grid: np.ndarray) -> np.ndarray:
        start_event.set()
        finish_event.wait()
        return np.ones((grid.shape[0], 1))

    original_cartesian = plotter_module.Tabulator.cartesian_grid
    cartesian_args: dict[str, bool] = {}

    def fake_cartesian(
        self: plotter_module.Tabulator,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        tabulate_gtos: bool = True,
    ) -> None:
        cartesian_args['tabulate_gtos'] = tabulate_gtos
        original_cartesian(self, x, y, z, tabulate_gtos)

    monkeypatch.setattr(plotter_module.Tabulator, 'cartesian_grid', fake_cartesian)
    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', fake_compute_gtos)
    monkeypatch.setattr(plotter_module.config.grid, 'default_type', 'cartesian', raising=False)

    plotter: plotter_module.Plotter | None = None
    try:
        root = _fake_tk_root()
        plotter = plotter_module.Plotter(_sample_molden(), tk_root=root)

        assert cartesian_args['tabulate_gtos'] is False
        assert start_event.wait(timeout=1.0)
        assert plotter._gto_future is not None
        assert not plotter._gto_future.done()
        assert plotter._selection_screen is not None
        assert isinstance(plotter._selection_screen, FakeSelectionScreen)
        assert plotter._selection_screen.loading_states == [True]
        assert plotter._selection_screen.messages == ['Tabulating orbitals...']
    finally:
        finish_event.set()
        if plotter is not None and plotter._gto_future is not None and not plotter._gtos_ready:
            plotter.wait_for_gtos()

    assert plotter is not None
    assert plotter._gtos_ready is True
    assert plotter._selection_screen is not None
    assert isinstance(plotter._selection_screen, FakeSelectionScreen)
    assert plotter._selection_screen.loading_states[-1] is False
    assert set(root.tk_call_thread_ids) == {root.owner_thread_id}


def test_wait_for_gtos_populates_data(monkeypatch: pytest.MonkeyPatch) -> None:
    """wait_for_gtos should block until background work finishes and expose data."""
    _configure_plotter_env(monkeypatch)

    def fake_compute_gtos(_tabulator: plotter_module.Tabulator, grid: np.ndarray) -> np.ndarray:
        return np.full((grid.shape[0], 1), 7.0)

    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', fake_compute_gtos)

    plotter = plotter_module.Plotter(_sample_molden(), tk_root=_fake_tk_root())
    plotter.wait_for_gtos()

    assert plotter._gtos_ready is True
    assert plotter._selection_screen is not None
    assert isinstance(plotter._selection_screen, FakeSelectionScreen)
    assert plotter._selection_screen.loading_states == [True, False]
    np.testing.assert_array_equal(
        plotter.tabulator.gtos,
        np.full((plotter.tabulator.grid.shape[0], 1), 7.0),
    )
    assert plotter._gto_future is None


def test_replacing_grid_discards_running_generation(monkeypatch: pytest.MonkeyPatch) -> None:
    """Only the newest grid snapshot may update the Tabulator and UI."""
    _configure_plotter_env(monkeypatch)
    first_started = threading.Event()
    release_first = threading.Event()
    second_started = threading.Event()
    release_second = threading.Event()

    def controlled_compute(_tabulator: plotter_module.Tabulator, grid: np.ndarray) -> np.ndarray:
        if not first_started.is_set():
            first_started.set()
            release_first.wait()
            return np.full((grid.shape[0], 1), 1.0)
        second_started.set()
        release_second.wait()
        return np.full((grid.shape[0], 1), 2.0)

    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', controlled_compute)
    plotter = plotter_module.Plotter(_sample_molden(), tk_root=_fake_tk_root())
    assert first_started.wait(timeout=1.0)
    selection_screen = plotter._selection_screen
    assert isinstance(selection_screen, FakeSelectionScreen)
    original_mesh = plotter._orb_mesh

    new_axis = np.linspace(-2.0, 2.0, 3)
    plotter._update_mesh(new_axis, new_axis, new_axis, plotter_module.GridType.CARTESIAN)
    release_first.set()
    assert second_started.wait(timeout=1.0)

    assert not plotter.tabulator.has_gtos
    assert plotter._orb_mesh is original_mesh
    assert selection_screen.loading_states == [True, True]

    release_second.set()
    plotter.wait_for_gtos()

    np.testing.assert_array_equal(
        plotter.tabulator.gtos,
        np.full((new_axis.size**3, 1), 2.0),
    )
    assert plotter._orb_mesh is not original_mesh
    assert selection_screen.loading_states == [True, True, False]


def test_rapid_grid_replacements_only_apply_latest(monkeypatch: pytest.MonkeyPatch) -> None:
    """Several running generations must leave only the newest grid installed."""
    _configure_plotter_env(monkeypatch)
    started = [threading.Event() for _ in range(3)]
    releases = [threading.Event() for _ in range(3)]
    call_lock = threading.Lock()
    call_count = 0

    def controlled_compute(_tabulator: plotter_module.Tabulator, grid: np.ndarray) -> np.ndarray:
        nonlocal call_count
        with call_lock:
            call_index = call_count
            call_count += 1
        started[call_index].set()
        releases[call_index].wait()
        return np.full((grid.shape[0], 1), call_index + 1.0)

    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', controlled_compute)
    plotter = plotter_module.Plotter(_sample_molden(), tk_root=_fake_tk_root())

    try:
        assert started[0].wait(timeout=1.0)
        middle_axis = np.linspace(-2.0, 2.0, 3)
        plotter._update_mesh(
            middle_axis,
            middle_axis,
            middle_axis,
            plotter_module.GridType.CARTESIAN,
        )
        releases[0].set()
        assert started[1].wait(timeout=1.0)

        latest_axis = np.linspace(-3.0, 3.0, 4)
        plotter._update_mesh(
            latest_axis,
            latest_axis,
            latest_axis,
            plotter_module.GridType.CARTESIAN,
        )
        releases[1].set()
        assert started[2].wait(timeout=1.0)

        assert not plotter.tabulator.has_gtos
        releases[2].set()
        plotter.wait_for_gtos()
    finally:
        for release in releases:
            release.set()

    np.testing.assert_array_equal(
        plotter.tabulator.gtos,
        np.full((latest_axis.size**3, 1), 3.0),
    )
    latest_axes = plotter.tabulator.grid_axes
    assert latest_axes is not None
    np.testing.assert_array_equal(latest_axes[0], latest_axis)
    assert isinstance(plotter._selection_screen, FakeSelectionScreen)
    assert plotter._selection_screen.loading_states == [True, True, True, False]


def test_simultaneous_plotters_keep_results_isolated(monkeypatch: pytest.MonkeyPatch) -> None:
    """Independent Plotters must receive their own serialized worker results."""
    _configure_plotter_env(monkeypatch)
    started = [threading.Event(), threading.Event()]
    releases = [threading.Event(), threading.Event()]
    call_lock = threading.Lock()
    call_count = 0

    def controlled_compute(_tabulator: plotter_module.Tabulator, grid: np.ndarray) -> np.ndarray:
        nonlocal call_count
        with call_lock:
            call_index = call_count
            call_count += 1
        started[call_index].set()
        releases[call_index].wait()
        return np.full((grid.shape[0], 1), call_index + 4.0)

    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', controlled_compute)
    first = plotter_module.Plotter(_sample_molden(), tk_root=_fake_tk_root())
    second = plotter_module.Plotter(_sample_molden(), tk_root=_fake_tk_root())

    try:
        assert started[0].wait(timeout=1.0)
        assert not started[1].is_set()
        releases[0].set()
        assert started[1].wait(timeout=1.0)
        releases[1].set()
        first.wait_for_gtos()
        second.wait_for_gtos()
    finally:
        for release in releases:
            release.set()

    np.testing.assert_array_equal(
        first.tabulator.gtos,
        np.full((first.tabulator.grid.shape[0], 1), 4.0),
    )
    np.testing.assert_array_equal(
        second.tabulator.gtos,
        np.full((second.tabulator.grid.shape[0], 1), 5.0),
    )


def test_closing_plotter_discards_running_generation(monkeypatch: pytest.MonkeyPatch) -> None:
    """A result finishing after close must not update Tabulator or UI state."""
    _configure_plotter_env(monkeypatch)
    started = threading.Event()
    release = threading.Event()
    callback_finished = threading.Event()

    def controlled_compute(_tabulator: plotter_module.Tabulator, grid: np.ndarray) -> np.ndarray:
        started.set()
        release.wait()
        callback_finished.set()
        return np.ones((grid.shape[0], 1))

    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', controlled_compute)
    plotter = plotter_module.Plotter(_sample_molden(), tk_root=_fake_tk_root())
    assert started.wait(timeout=1.0)
    selection_screen = plotter._selection_screen
    assert isinstance(selection_screen, FakeSelectionScreen)
    original_mesh = plotter._orb_mesh

    plotter._on_screen = False
    plotter._cancel_gto_future()
    release.set()
    assert callback_finished.wait(timeout=1.0)

    assert not plotter.tabulator.has_gtos
    assert plotter._gtos_ready is False
    assert plotter._orb_mesh is original_mesh
    assert selection_screen.loading_states == [True]


def test_selection_window_close_during_work_preserves_custom_root(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Selection-window close must discard work without owning a supplied root."""
    _configure_plotter_env(monkeypatch)
    started = threading.Event()
    release = threading.Event()
    finished = threading.Event()

    def controlled_compute(_tabulator: plotter_module.Tabulator, grid: np.ndarray) -> np.ndarray:
        started.set()
        release.wait()
        finished.set()
        return np.ones((grid.shape[0], 1))

    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', controlled_compute)
    root = cast(FakeTk, _fake_tk_root())
    plotter = plotter_module.Plotter(_sample_molden(), tk_root=cast(Any, root))

    try:
        assert started.wait(timeout=1.0)
        selection_screen = cast(FakeSelectionScreen, plotter._selection_screen)
        _REAL_SELECTION_SCREEN._on_close(cast(Any, selection_screen))
    finally:
        release.set()
    assert finished.wait(timeout=1.0)

    assert not plotter.tabulator.has_gtos
    assert selection_screen.destroyed
    assert cast(FakeBackgroundPlotter, plotter._pv_plotter).closed
    assert root.quit_calls == 0
    assert root.destroy_calls == 0


def test_plotter_window_close_during_work_preserves_custom_root(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """PyVista close must discard work without quitting a supplied Tk root."""
    _configure_plotter_env(monkeypatch)
    started = threading.Event()
    release = threading.Event()
    finished = threading.Event()

    def controlled_compute(_tabulator: plotter_module.Tabulator, grid: np.ndarray) -> np.ndarray:
        started.set()
        release.wait()
        finished.set()
        return np.ones((grid.shape[0], 1))

    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', controlled_compute)
    root = cast(FakeTk, _fake_tk_root())
    plotter = plotter_module.Plotter(_sample_molden(), tk_root=cast(Any, root))
    plotter_module._PlotterRendering._connect_pv_plotter_close_signal(plotter)

    try:
        assert started.wait(timeout=1.0)
        signal = cast(FakeSignal, plotter._pv_plotter.app_window.signal_close)
        signal.emit()
    finally:
        release.set()
    assert finished.wait(timeout=1.0)

    assert not plotter.tabulator.has_gtos
    assert isinstance(plotter._selection_screen, FakeSelectionScreen)
    assert plotter._selection_screen.destroyed
    assert root.quit_calls == 0
    assert root.destroy_calls == 0


def test_gto_success_is_delivered_by_tk_thread(monkeypatch: pytest.MonkeyPatch) -> None:
    """A worker should publish data without making any Tk calls."""
    _configure_plotter_env(monkeypatch)
    worker_thread_ids: list[int] = []

    def fake_compute_gtos(_tabulator: plotter_module.Tabulator, grid: np.ndarray) -> np.ndarray:
        worker_thread_ids.append(threading.get_ident())
        return np.ones((grid.shape[0], 1))

    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', fake_compute_gtos)
    root = cast(FakeTk, _fake_tk_root())
    plotter = plotter_module.Plotter(_sample_molden(), tk_root=cast(Any, root))

    _pump_until(root, lambda: plotter._gtos_ready)

    assert worker_thread_ids
    assert worker_thread_ids[0] != root.owner_thread_id
    assert set(root.tk_call_thread_ids) == {root.owner_thread_id}


def test_gto_success_is_delivered_by_real_tcl_event_loop(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A headless Tcl event loop should deliver completion on its owner thread."""
    _configure_plotter_env(monkeypatch)
    delivery_thread_ids: list[int] = []
    owner_thread_id = threading.get_ident()
    original_apply = plotter_module.Plotter._apply_gtos_ready

    def fake_compute_gtos(_tabulator: plotter_module.Tabulator, grid: np.ndarray) -> np.ndarray:
        return np.ones((grid.shape[0], 1))

    def record_apply(
        self: plotter_module.Plotter,
        result: plotter_module._GTOResult,
        elapsed: float,
    ) -> None:
        delivery_thread_ids.append(threading.get_ident())
        original_apply(self, result, elapsed)

    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', fake_compute_gtos)
    monkeypatch.setattr(plotter_module.Plotter, '_apply_gtos_ready', record_apply)
    root = cast(Any, plotter_module.tk.Tcl())
    plotter = plotter_module.Plotter(_sample_molden(), tk_root=root)

    try:
        deadline = time.monotonic() + 1.0
        while not plotter._gtos_ready:
            assert time.monotonic() < deadline
            root.dooneevent(0)
    finally:
        plotter._on_screen = False
        plotter._cancel_gto_future()

    assert delivery_thread_ids == [owner_thread_id]


def test_gto_failure_is_delivered_by_tk_thread(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A background failure should show its message from the Tk thread."""
    _configure_plotter_env(monkeypatch)
    shown_on_threads: list[int] = []

    def fail_compute(_tabulator: plotter_module.Tabulator, _grid: np.ndarray) -> np.ndarray:
        raise ValueError('broken grid')

    def fake_showerror(_title: str, _message: str) -> None:
        shown_on_threads.append(threading.get_ident())

    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', fail_compute)
    monkeypatch.setattr(plotter_module.messagebox, 'showerror', fake_showerror)
    root = cast(FakeTk, _fake_tk_root())
    plotter_module.Plotter(_sample_molden(), tk_root=cast(Any, root))

    _pump_until(root, lambda: bool(shown_on_threads))

    assert shown_on_threads == [root.owner_thread_id]
    assert set(root.tk_call_thread_ids) == {root.owner_thread_id}


def test_failed_grid_replacement_preserves_previous_state(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A failed replacement should restore the previously usable grid."""
    _configure_plotter_env(monkeypatch)
    fail_replacement = False

    def controlled_compute(_tabulator: plotter_module.Tabulator, grid: np.ndarray) -> np.ndarray:
        if fail_replacement:
            raise ValueError('broken replacement')
        return np.ones((grid.shape[0], 1))

    shown_errors: list[str] = []
    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', controlled_compute)
    monkeypatch.setattr(
        plotter_module.messagebox,
        'showerror',
        lambda _title, message: shown_errors.append(message),
    )
    root = cast(FakeTk, _fake_tk_root())
    plotter = plotter_module.Plotter(_sample_molden(), tk_root=cast(Any, root))
    _pump_until(root, lambda: plotter._gtos_ready)

    previous_grid = plotter.tabulator.grid
    previous_gtos = plotter.tabulator.gtos
    previous_mesh = plotter._orb_mesh
    fail_replacement = True
    new_axis = np.linspace(-2.0, 2.0, 3)

    plotter._update_mesh(new_axis, new_axis, new_axis, plotter_module.GridType.CARTESIAN)
    _pump_until(root, lambda: bool(shown_errors))

    assert plotter._gtos_ready
    assert plotter.tabulator.grid is previous_grid
    assert plotter.tabulator.gtos is previous_gtos
    assert plotter._orb_mesh is previous_mesh
    assert isinstance(plotter._selection_screen, FakeSelectionScreen)
    assert plotter._selection_screen.loading_states[-2:] == [True, False]
    assert 'broken replacement' in shown_errors[0]


def test_close_stops_gto_delivery(monkeypatch: pytest.MonkeyPatch) -> None:
    """Closing the UI should cancel polling and discard queued delivery."""
    _configure_plotter_env(monkeypatch)
    finish_event = threading.Event()

    def fake_compute_gtos(_tabulator: plotter_module.Tabulator, grid: np.ndarray) -> np.ndarray:
        finish_event.wait()
        return np.ones((grid.shape[0], 1))

    monkeypatch.setattr(plotter_module.Tabulator, 'compute_gtos', fake_compute_gtos)
    root = cast(FakeTk, _fake_tk_root())
    plotter = plotter_module.Plotter(_sample_molden(), tk_root=cast(Any, root))

    plotter._on_screen = False
    plotter._cancel_gto_future()
    finish_event.set()
    time.sleep(0.01)

    assert not root.callbacks
    assert plotter._gtos_ready is False
