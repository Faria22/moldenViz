"""Tests covering asynchronous GTO tabulation behavior in Plotter."""

from __future__ import annotations

import threading
from collections import UserDict
from pathlib import Path
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any, cast

import numpy as np

from moldenViz import plotter as plotter_module

if TYPE_CHECKING:
    import pytest


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


class FakeBackgroundPlotter:
    """Lightweight BackgroundPlotter stand-in for headless tests."""

    def __init__(self, *_args: object, **_kwargs: object) -> None:
        self.app_window = SimpleNamespace(signal_close=SimpleNamespace(connect=lambda _cb: None))
        self.main_menu = SimpleNamespace(actions=list, addMenu=lambda *_args, **_kwargs: None)
        self._last_mesh_args: tuple[tuple[object, ...], dict[str, object]] | None = None

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

    def update(self) -> None:  # pragma: no cover - trivial
        """Expose the update hook expected by Plotter."""


class FakeTk:
    """Minimal Tk root substitute used to execute callbacks immediately."""

    def __init__(self) -> None:
        self.idle_callbacks: list[tuple[object, tuple[object, ...]]] = []

    def after_idle(self, callback: object, *args: object) -> None:
        """Invoke callbacks immediately while tracking invocations for tests."""
        assert callable(callback)
        self.idle_callbacks.append((callback, args))
        callback(*args)

    def mainloop(self) -> None:  # pragma: no cover - trivial
        """Expose Tk's mainloop hook for compatibility."""

    def quit(self) -> None:  # pragma: no cover - trivial
        """Satisfy the Tk API contract for quitting."""

    def destroy(self) -> None:  # pragma: no cover - trivial
        """Simulate tearing down the Tk root."""


class FakeSelectionScreen:
    """Stub orbital selection dialog used to capture loading indicator states."""

    def __init__(self, plotter: plotter_module.Plotter) -> None:
        self.plotter = plotter
        self.current_mo_ind = -1
        self.loading_states: list[bool] = []
        self.messages: list[str] = []

    def set_loading_state(self, loading: bool, message: str = 'Tabulating orbitals...') -> None:
        """Mirror the dialog's loading indicator for assertions."""
        self.loading_states.append(loading)
        self.messages.append(message)

    def on_gtos_ready(self) -> None:
        """Transition the dialog back to idle once GTOs arrive."""
        self.set_loading_state(False)

    def update_nav_button_states(self) -> None:  # pragma: no cover - trivial
        """Stub navigation button updates."""

    def winfo_exists(self) -> bool:  # pragma: no cover - trivial
        """Report that the dialog window is still alive.

        Returns
        -------
        bool
            True while the fake selection dialog is attached to a plotter.
        """
        return bool(self.plotter)

    def destroy(self) -> None:  # pragma: no cover - trivial
        """Implement the Tk destroy hook."""


def _configure_plotter_env(monkeypatch: pytest.MonkeyPatch) -> None:
    """Patch heavyweight UI helpers so Plotter can be instantiated in tests."""

    def fake_load_molecule(self: plotter_module.Plotter, _config: object) -> None:
        dummy = SimpleNamespace(max_radius=1.0, atoms=[])
        self.molecule = cast(Any, dummy)
        self.molecule_actors = []
        self.atom_actors = []
        self.bond_actors = []

    monkeypatch.setattr(plotter_module, 'BackgroundPlotter', FakeBackgroundPlotter)
    monkeypatch.setattr(plotter_module, '_OrbitalSelectionScreen', FakeSelectionScreen)
    monkeypatch.setattr(plotter_module.Plotter, '_add_orbital_menus_to_pv_plotter', lambda _plotter: None)
    monkeypatch.setattr(plotter_module.Plotter, '_connect_pv_plotter_close_signal', lambda _plotter: None)
    monkeypatch.setattr(plotter_module.Plotter, '_override_clear_all_button', lambda _plotter: None)
    monkeypatch.setattr(plotter_module.Plotter, 'load_molecule', fake_load_molecule)
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


def test_plotter_defers_gto_tabulation(monkeypatch: pytest.MonkeyPatch) -> None:
    """Default Plotter grids should tabulate GTOs in the background."""
    _configure_plotter_env(monkeypatch)
    start_event = threading.Event()
    finish_event = threading.Event()

    def fake_tabulate_gtos(_tabulator: plotter_module.Tabulator) -> np.ndarray:
        start_event.set()
        finish_event.wait()
        return np.ones((1, 1))

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
    monkeypatch.setattr(plotter_module.Tabulator, 'tabulate_gtos', fake_tabulate_gtos)

    plotter: plotter_module.Plotter | None = None
    try:
        plotter = plotter_module.Plotter(_sample_molden(), tk_root=_fake_tk_root())

        assert cartesian_args['tabulate_gtos'] is False
        assert start_event.wait(timeout=1.0)
        assert plotter._gto_future is not None  # noqa: SLF001
        assert not plotter._gto_future.done()  # noqa: SLF001
        assert plotter.selection_screen is not None
        assert isinstance(plotter.selection_screen, FakeSelectionScreen)
        assert plotter.selection_screen.loading_states == [True]
        assert plotter.selection_screen.messages == ['Tabulating orbitals...']
    finally:
        finish_event.set()
        if plotter is not None and plotter._gto_future is not None and not plotter._gtos_ready:  # noqa: SLF001
            plotter.wait_for_gtos()

    assert plotter is not None
    assert plotter._gtos_ready is True  # noqa: SLF001
    assert plotter.selection_screen is not None
    assert isinstance(plotter.selection_screen, FakeSelectionScreen)
    assert plotter.selection_screen.loading_states[-1] is False


def test_wait_for_gtos_populates_data(monkeypatch: pytest.MonkeyPatch) -> None:
    """wait_for_gtos should block until background work finishes and expose data."""
    _configure_plotter_env(monkeypatch)

    def fake_tabulate_gtos(_tabulator: plotter_module.Tabulator) -> np.ndarray:
        return np.full((2, 1), 7.0)

    monkeypatch.setattr(plotter_module.Tabulator, 'tabulate_gtos', fake_tabulate_gtos)

    plotter = plotter_module.Plotter(_sample_molden(), tk_root=_fake_tk_root())
    plotter.wait_for_gtos()

    assert plotter._gtos_ready is True  # noqa: SLF001
    assert plotter.selection_screen is not None
    assert isinstance(plotter.selection_screen, FakeSelectionScreen)
    assert plotter.selection_screen.loading_states == [True, False]
    np.testing.assert_array_equal(plotter.tabulator.gtos, np.full((2, 1), 7.0))
    assert plotter._gto_future is None  # noqa: SLF001
