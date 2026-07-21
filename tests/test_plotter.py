"""Tests covering Plotter behaviours called out by plotter_coverage_gaps.md."""
# ruff: noqa: D101, D102, D103

from __future__ import annotations

from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterator

import numpy as np
import pytest
from matplotlib import colors as mcolors

from moldenViz import Tabulator

plotter_module = pytest.importorskip('moldenViz.plotter')


class _GridTypeProxy:
    def __getattr__(self, item: str) -> Any:
        return getattr(plotter_module.GridType, item)


GridType = _GridTypeProxy()


class DummySignal:
    def __init__(self) -> None:
        self.callbacks: list = []

    def connect(self, callback: object) -> None:
        self.callbacks.append(callback)

    def disconnect(self) -> bool:
        if self.callbacks:
            self.callbacks.pop()
            return True
        return False


class DummyActor:
    def __init__(self) -> None:
        self.visible = True
        self.opacity = 1.0

    def SetVisibility(self, value: bool) -> None:  # noqa: N802
        self.visible = bool(value)

    def GetVisibility(self) -> bool:  # noqa: N802
        return self.visible

    def GetProperty(self) -> SimpleNamespace:  # noqa: N802
        return SimpleNamespace(SetOpacity=lambda val: setattr(self, 'opacity', val))


class DummyMenuBar:
    def __init__(self) -> None:
        self.menus: list = []

    def addMenu(self, menu: object) -> None:  # noqa: N802 - Qt naming
        self.menus.append(menu)

    @staticmethod
    def actions() -> list:
        return []


class DummyMenuAction:
    def __init__(self, text: str, menu: DummyMenuWithActions | None = None) -> None:
        self._text = text
        self._menu = menu
        self.triggered = DummySignal()

    def text(self) -> str:
        return self._text

    def menu(self) -> DummyMenuWithActions | None:
        return self._menu


class DummyMenuWithActions:
    def __init__(self, title: str) -> None:
        self.title = title
        self._actions: list[DummyMenuAction] = []

    def addAction(self, action: DummyMenuAction) -> None:  # noqa: N802 - match Qt API
        self._actions.append(action)

    def addSeparator(self) -> None:  # noqa: N802 - match Qt API
        self._actions.append(DummyMenuAction('---'))

    def actions(self) -> list[DummyMenuAction]:
        return self._actions


class DummyBackgroundPlotter:
    def __init__(self, *_args: object, **_kwargs: object) -> None:
        self.background: str | None = None
        self.removed_actors: list = []
        self.added_meshes: list = []
        self.saved_graphic: str | None = None
        self.screenshot_calls: list[tuple[str, bool]] = []
        self.update_count = 0
        self.app_window = SimpleNamespace(signal_close=DummySignal())
        self.main_menu = DummyMenuBar()

    def set_background(self, color: str) -> None:
        self.background = color

    def show_axes(self) -> None:  # pragma: no cover - noop stub
        pass

    def remove_actor(self, actor: object) -> None:
        self.removed_actors.append(actor)

    def add_mesh(self, mesh: object, **kwargs: object) -> DummyActor:
        actor = DummyActor()
        self.added_meshes.append((mesh, kwargs))
        return actor

    def update(self) -> None:  # pragma: no cover - noop stub
        self.update_count += 1

    def save_graphic(self, path: str) -> None:
        self.saved_graphic = path

    def screenshot(self, path: str, transparent_background: bool = False) -> None:
        self.screenshot_calls.append((path, transparent_background))
        self.screenshot_args = (path, transparent_background)


class MenuAwareMainMenu:
    def __init__(self) -> None:
        self.menus: list[Any] = []
        self.view_menu = DummyMenuWithActions('View')
        self.clear_action = DummyMenuAction('Clear All')
        self.view_menu.addAction(self.clear_action)
        self._view_action = DummyMenuAction('View', menu=self.view_menu)

    def addMenu(self, menu: Any) -> None:  # noqa: N802 - match Qt API
        self.menus.append(menu)

    def actions(self) -> list[DummyMenuAction]:
        return [self._view_action]


class MenuAwareBackgroundPlotter(DummyBackgroundPlotter):
    def __init__(self, *_args: object, **_kwargs: object) -> None:
        super().__init__(*_args, **_kwargs)
        self.main_menu = MenuAwareMainMenu()


class FakeQAction:
    def __init__(self, text: str, _parent: Any) -> None:
        self._text = text
        self.triggered = DummySignal()

    def text(self) -> str:
        return self._text


class FakeQMenu:
    def __init__(self, title: str, _parent: Any) -> None:
        self.title = title
        self._actions: list[FakeQAction] = []
        self.separators = 0

    def addAction(self, action: FakeQAction) -> None:  # noqa: N802 - match Qt API
        self._actions.append(action)

    def addSeparator(self) -> None:  # noqa: N802 - match Qt API
        self.separators += 1

    def actions(self) -> list[FakeQAction]:
        return self._actions


class DummyMolecule:
    def __init__(self, atoms: list, _config: object) -> None:
        self.atoms = atoms
        self.max_radius = 1.0

    def add_meshes(  # noqa: PLR6301
        self,
        _plotter: DummyBackgroundPlotter,
        opacity: float,
    ) -> tuple[list[DummyActor], list[DummyActor], list[DummyActor]]:
        actors = [DummyActor()]
        atoms = [DummyActor()]
        bonds = [DummyActor()]
        for actor in actors:
            actor.opacity = opacity
        return actors, atoms, bonds


class DummySelectionScreen:
    def __init__(self, plotter: Any) -> None:
        self.plotter = plotter
        self.current_mo_ind = -1
        self.destroyed = False
        self._loading = False
        self.loading_events: list[tuple[bool, str]] = []
        self.last_loading_message = 'Tabulating orbitals...'

    def update_nav_button_states(self) -> None:  # pragma: no cover - noop stub
        pass

    def set_loading_state(self, loading: bool, message: str = 'Tabulating orbitals...') -> None:
        self._loading = loading
        self.last_loading_message = message
        self.loading_events.append((loading, message))

    def on_gtos_ready(self) -> None:
        self.set_loading_state(False)

    def destroy(self) -> None:  # pragma: no cover - noop stub
        self.destroyed = True

    def winfo_exists(self) -> bool:
        return not self.destroyed


class DummyStructuredGrid:
    def __init__(self) -> None:
        self.points: np.ndarray | None = None
        self.dimensions: tuple[int, int, int] | None = None
        self.arrays: dict[str, np.ndarray] = {}

    def __setitem__(self, key: str, value: np.ndarray) -> None:
        """Store array on the fake grid."""
        self.arrays[key] = value

    def contour(self, values: list[float]) -> dict:
        return {'levels': tuple(values), 'points': self.points}


class DummyTk:
    def __init__(self) -> None:
        self.withdrawn = False
        self.mainloop_calls = 0
        self.quit_calls = 0

    def withdraw(self) -> None:
        self.withdrawn = True

    def mainloop(self) -> None:
        self.mainloop_calls += 1

    def quit(self) -> None:
        self.quit_calls += 1


class DummyWindow:
    def __init__(self) -> None:
        self.destroyed = False

    def destroy(self) -> None:
        self.destroyed = True


class DummyVar:
    def __init__(self, value: Any = '') -> None:
        self._value = value

    def get(self) -> Any:
        return self._value

    def set(self, value: Any) -> None:
        self._value = value

    def trace_add(self, mode: str, callback: Any) -> None:
        self._last_trace = (mode, callback)


class DummyEntry:
    def __init__(self, value: str = '') -> None:
        self.value = value

    def get(self) -> str:
        return self.value

    def insert(self, _index: int, text: str) -> None:
        self.value = text

    def delete(self, _start: int, _end: int | None = None) -> None:
        self.value = ''


class DummyLabelWidget:
    def __init__(self, text: str) -> None:
        self.text = text

    def cget(self, _key: str) -> str:
        return self.text

    def config(self, **kwargs: Any) -> None:
        if 'text' in kwargs:
            self.text = kwargs['text']


class DummyContainer:
    def __init__(self, children: list[Any]) -> None:
        self._children = children

    def winfo_children(self) -> list[Any]:
        return self._children


class SimpleWidget:
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self.args = args
        self.kwargs = kwargs
        self.children: list[Any] = []
        self.value = ''

    def grid(self, *_args: Any, **_kwargs: Any) -> SimpleWidget:
        return self

    def pack(self, *_args: Any, **_kwargs: Any) -> SimpleWidget:
        return self

    def columnconfigure(self, *args: Any, **kwargs: Any) -> None:
        self._columnconfigure = (args, kwargs)

    def rowconfigure(self, *args: Any, **kwargs: Any) -> None:
        self._rowconfigure = (args, kwargs)

    def config(self, **kwargs: Any) -> SimpleWidget:
        self.kwargs.update(kwargs)
        return self

    def bind(self, *_args: Any, **_kwargs: Any) -> SimpleWidget:
        return self

    def grid_remove(self) -> SimpleWidget:
        return self

    def grid_forget(self) -> SimpleWidget:
        return self

    def winfo_children(self) -> list[Any]:
        return self.children

    def set(self, value: Any) -> None:
        self.value = value


class SimpleEntry(SimpleWidget):
    def insert(self, _index: int, text: str) -> None:
        self.value = text

    def delete(self, *_args: Any) -> None:
        self.value = ''

    def get(self) -> str:
        return self.value


class SimpleCombobox(SimpleWidget):
    def __init__(self, *args: Any, textvariable: DummyVar | None = None, values: Any = None, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.textvariable = textvariable
        self.values = values or []
        self._items: dict[str, Any] = {}

    def current(self, index: int) -> None:
        if self.textvariable is not None and 0 <= index < len(self.values):
            self.textvariable.set(self.values[index])

    def __setitem__(self, key: str, value: Any) -> None:
        """Store Tk-style option value."""
        self._items[key] = value

    def __getitem__(self, key: str) -> Any:
        """Return Tk-style option value.

        Returns
        -------
        Any
            Stored value for the requested option.
        """
        return self._items.get(key)


class SimpleToplevel(SimpleWidget):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self.protocol_handlers: dict[str, Any] = {}
        self.destroyed = False
        self.parent_args = args

    def title(self, text: str) -> None:
        self._title = text

    def geometry(self, size: str) -> None:
        self._geometry = size

    def destroy(self) -> None:
        self.destroyed = True

    def protocol(self, name: str, handler: Any) -> None:
        self.protocol_handlers[name] = handler

    def winfo_exists(self) -> bool:
        return not self.destroyed


class SimpleTreeview(SimpleWidget):
    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._headings: dict[str, str] = {}
        self._columns: dict[str, int] = {}
        self._items: dict[Any, dict[str, Any]] = {}
        self._selection: tuple[str, ...] = ()
        self._bindings: dict[str, Any] = {}
        self._tags: dict[str, Any] = {}
        self.seen: Any = None

    def heading(self, column: str, text: str) -> None:
        self._headings[column] = text

    def column(self, column: str, width: int) -> None:
        self._columns[column] = width

    def tag_configure(self, tag: str, **options: Any) -> None:
        self._tags[tag] = options

    def bind(self, event: str, handler: Any) -> SimpleWidget:
        self._bindings[event] = handler
        return super().bind(event, handler)

    def insert(self, parent: str, index: str, iid: Any, values: Any = None) -> Any:
        self._items[iid] = {'parent': parent, 'index': index, 'values': values or (), 'tags': ()}
        return iid

    def item(self, iid: Any, **kwargs: Any) -> dict[str, Any]:
        item = self._items.setdefault(iid, {'values': (), 'tags': ()})
        if 'tags' in kwargs:
            item['tags'] = kwargs['tags']
        return item

    def see(self, iid: Any) -> None:
        self.seen = iid

    def get_children(self) -> list[Any]:
        return list(self._items.keys())

    def delete(self, iid: Any) -> None:
        self._items.pop(iid, None)

    def selection(self) -> tuple[str, ...]:
        return self._selection

    def selection_remove(self, items: Any) -> None:
        removals = {str(item) for item in items}
        self._selection = tuple(entry for entry in self._selection if entry not in removals)

    def set_selection(self, iid: Any) -> None:
        self._selection = (str(iid),)


class FakeTabulator:
    def __init__(self, _source: Any = None, only_molecule: bool = False, **_kwargs: Any) -> None:
        grid_enum = plotter_module.GridType
        self._grid_type = grid_enum.SPHERICAL
        self._grid_dimensions = (1, 1, 1)
        self.grid = np.zeros((1, 3))
        self.original_axes = None
        self.gto_data = object()
        self._gtos = np.zeros((1, 1))
        self._parser = SimpleNamespace(
            atoms=[SimpleNamespace(symbol='H', coords=(0.0, 0.0, 0.0))],
            mos=[SimpleNamespace(sym='s', spin='alpha', occ=2.0, energy=-0.5)],
        )
        self.export_calls: list[tuple[str, int | None]] = []
        self.only_molecule = only_molecule

    def spherical_grid(
        self,
        r: np.ndarray,
        theta: np.ndarray,
        phi: np.ndarray,
        tabulate_gtos: bool = True,
    ) -> None:
        rr, tt, pp = np.meshgrid(r, theta, phi, indexing='ij')
        xx, yy, zz = Tabulator._spherical_to_cartesian(rr, tt, pp)  # noqa: SLF001
        self.grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        self._grid_type = GridType.SPHERICAL
        self._grid_dimensions = (len(r), len(theta), len(phi))
        self.original_axes = (r, theta, phi)
        if tabulate_gtos:
            self._gtos = np.zeros((self.grid.shape[0], 1))

    def cartesian_grid(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        tabulate_gtos: bool = True,
    ) -> None:
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        self.grid = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        self._grid_type = GridType.CARTESIAN
        self._grid_dimensions = (len(x), len(y), len(z))
        self.original_axes = (x, y, z)
        if tabulate_gtos:
            self._gtos = np.zeros((self.grid.shape[0], 1))

    def tabulate_mos(self, _orb_ind: int) -> np.ndarray:
        return np.arange(self.grid.shape[0], dtype=float)

    def export(self, path: str, mo_index: int | None) -> None:
        self.export_calls.append((path, mo_index))

    def tabulate_gtos(self) -> np.ndarray:
        gtos = np.ones((self.grid.shape[0], 1))
        self._gtos = gtos
        return gtos


def seed_tabulator_with_cartesian_grid(tabulator: FakeTabulator) -> None:
    points = np.linspace(0.0, 1.0, 2)
    tabulator.cartesian_grid(points, points, points)


class RecordingTabulator(FakeTabulator):
    def __init__(self, source: Any = None, only_molecule: bool = False, **kwargs: Any) -> None:
        super().__init__(source, only_molecule=only_molecule, **kwargs)
        self.spherical_calls: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
        self.cartesian_calls: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []

    def spherical_grid(
        self,
        r: np.ndarray,
        theta: np.ndarray,
        phi: np.ndarray,
        tabulate_gtos: bool = True,
    ) -> None:
        self.spherical_calls.append((r, theta, phi))
        super().spherical_grid(r, theta, phi, tabulate_gtos=tabulate_gtos)

    def cartesian_grid(
        self,
        x: np.ndarray,
        y: np.ndarray,
        z: np.ndarray,
        tabulate_gtos: bool = True,
    ) -> None:
        self.cartesian_calls.append((x, y, z))
        super().cartesian_grid(x, y, z, tabulate_gtos=tabulate_gtos)


class RootRecorder:
    def __init__(self) -> None:
        self.quit_calls = 0
        self.destroy_calls = 0

    def quit(self) -> None:
        self.quit_calls += 1

    def destroy(self) -> None:
        self.destroy_calls += 1


class PVRecorder:
    def __init__(self) -> None:
        self.closed = False

    def close(self) -> None:
        self.closed = True


class SelectionPlotter:
    def __init__(self) -> None:
        self.tk_root = RootRecorder()
        self._no_prev_tk_root = True
        self.tabulator = FakeTabulator()
        self.tabulator._parser.mos = [  # noqa: SLF001
            SimpleNamespace(sym='s', spin='alpha', occ=2.0, energy=-0.5),
            SimpleNamespace(sym='p', spin='alpha', occ=1.0, energy=-0.1),
            SimpleNamespace(sym='d', spin='beta', occ=0.0, energy=0.2),
        ]
        self.pv_plotter = PVRecorder()
        self.on_screen = True
        self.selection_screen: Any | None = None
        self.plot_calls: list[int] = []

    def plot_orbital(self, idx: int) -> None:
        self.plot_calls.append(idx)
        if self.selection_screen is not None:
            self.selection_screen.current_mo_ind = idx

    def _cancel_gto_future(self) -> None:  # pragma: no cover - stubbed out
        pass


@pytest.fixture
def plotter_env(monkeypatch: pytest.MonkeyPatch) -> Any:
    """Provide helper factory plus global patches for Plotter instantiation.

    Returns
    -------
    Any
        Helper object that creates patched Plotter instances.
    """
    monkeypatch.setattr(plotter_module, 'config', plotter_module.Config())
    monkeypatch.setattr(plotter_module, 'BackgroundPlotter', DummyBackgroundPlotter)
    monkeypatch.setattr(plotter_module, 'Molecule', DummyMolecule)
    monkeypatch.setattr(plotter_module, '_OrbitalSelectionScreen', DummySelectionScreen)
    monkeypatch.setattr(plotter_module.pv, 'StructuredGrid', DummyStructuredGrid)
    monkeypatch.setattr(plotter_module.pv, 'pyvista_ndarray', lambda arr: arr)
    monkeypatch.setattr(plotter_module.Plotter, '_add_orbital_menus_to_pv_plotter', lambda _self: None)
    monkeypatch.setattr(plotter_module.Plotter, '_override_clear_all_button', lambda _self: None)

    class Env:
        @staticmethod
        def make_root() -> DummyTk:
            return DummyTk()

        @staticmethod
        def make_tabulator() -> FakeTabulator:
            return FakeTabulator()

        def make_plotter(
            self,
            *,
            tabulator: FakeTabulator | None = None,
            only_molecule: bool = False,
            root: DummyTk | None = None,
        ) -> Any:
            root = root or self.make_root()
            tabulator = tabulator or self.make_tabulator()
            return plotter_module.Plotter('dummy', tabulator=tabulator, only_molecule=only_molecule, tk_root=root)

    return Env()


def test_describe_source_reports_list_length() -> None:
    assert plotter_module._describe_source('sample.molden') == 'sample.molden'  # noqa: SLF001
    assert plotter_module._describe_source(['a', 'b', 'c']) == '3 molden lines'  # noqa: SLF001


def test_custom_cmap_from_colors_uses_endpoints() -> None:
    cmap = plotter_module.Plotter.custom_cmap_from_colors(['red', 'blue'])
    assert cmap.name == 'custom_mo'
    assert cmap(0) == pytest.approx(mcolors.to_rgba('red'))
    almost_one = np.nextafter(1.0, 0.0)
    assert cmap(almost_one) == pytest.approx(mcolors.to_rgba('blue'))


def test_settings_parent_prefers_selection_screen(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    assert plotter._settings_parent() is plotter.selection_screen  # noqa: SLF001


def test_settings_parent_requires_root(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.selection_screen = None
    plotter.tk_root = None
    with pytest.raises(RuntimeError):
        plotter._settings_parent()  # noqa: SLF001


def test_get_current_mo_index_tracks_selection(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    target_index = 4
    plotter.selection_screen.current_mo_ind = target_index
    assert plotter._get_current_mo_index() == target_index  # noqa: SLF001
    plotter.selection_screen = None
    missing_index = -1
    assert plotter._get_current_mo_index() == missing_index  # noqa: SLF001


def test_do_export_with_all_scope_calls_tabulator(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.selection_screen.current_mo_ind = 2
    export_window = DummyWindow()

    messagebox_calls: list[tuple[str, str]] = []

    monkeypatch.setattr(plotter_module.filedialog, 'asksaveasfilename', lambda **_kwargs: '/tmp/export.vtk')
    monkeypatch.setattr(
        plotter_module.messagebox,
        'showinfo',
        lambda title, msg: messagebox_calls.append((title, msg)),
    )
    monkeypatch.setattr(
        plotter_module.messagebox,
        'showerror',
        lambda *_args, **_kwargs: pytest.fail('showerror called'),
    )

    plotter._do_export(export_window, DummyVar('vtk'), DummyVar('all'))  # noqa: SLF001

    assert plotter.tabulator.export_calls == [('/tmp/export.vtk', None)]
    assert export_window.destroyed
    assert messagebox_calls


def test_do_export_requires_selected_orbital(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.selection_screen.current_mo_ind = -1

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(
        plotter_module.messagebox,
        'showerror',
        lambda title, msg: errors.append((title, msg)),
    )
    monkeypatch.setattr(
        plotter_module.filedialog,
        'asksaveasfilename',
        lambda **_kwargs: pytest.fail('Dialog should not open'),
    )

    plotter._do_export(DummyWindow(), DummyVar('vtk'), DummyVar('current'))  # noqa: SLF001
    assert errors
    assert 'No orbital' in errors[0][1]


def test_do_export_rejects_cube_all(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.selection_screen.current_mo_ind = 0

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(
        plotter_module.messagebox,
        'showerror',
        lambda title, msg: errors.append((title, msg)),
    )
    monkeypatch.setattr(
        plotter_module.filedialog,
        'asksaveasfilename',
        lambda **_kwargs: pytest.fail('Dialog should not open'),
    )

    plotter._do_export(DummyWindow(), DummyVar('cube'), DummyVar('all'))  # noqa: SLF001
    assert errors
    assert 'Cube format' in errors[0][1]


def test_do_export_uses_current_index(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.selection_screen.current_mo_ind = 3
    export_window = DummyWindow()

    monkeypatch.setattr(plotter_module.filedialog, 'asksaveasfilename', lambda **_kwargs: '/tmp/single.vtk')
    monkeypatch.setattr(plotter_module.messagebox, 'showinfo', lambda *_args, **_kwargs: None)
    monkeypatch.setattr(plotter_module.messagebox, 'showerror', lambda *_args, **_kwargs: None)

    plotter._do_export(export_window, DummyVar('vtk'), DummyVar('current'))  # noqa: SLF001
    assert plotter.tabulator.export_calls == [('/tmp/single.vtk', 3)]


def test_do_export_handles_errors(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.selection_screen.current_mo_ind = 0
    export_window = DummyWindow()

    def boom(*_args: Any, **_kwargs: Any) -> None:
        raise RuntimeError('fail')

    plotter.tabulator.export = boom  # type: ignore[assignment]
    monkeypatch.setattr(plotter_module.filedialog, 'asksaveasfilename', lambda **_kwargs: '/tmp/single.vtk')

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(plotter_module.messagebox, 'showerror', lambda title, msg: errors.append((title, msg)))
    monkeypatch.setattr(plotter_module.messagebox, 'showinfo', lambda *_args, **_kwargs: None)

    plotter._do_export(export_window, DummyVar('vtk'), DummyVar('current'))  # noqa: SLF001
    assert errors


def test_plotter_rejects_tabulator_without_grid(plotter_env: Any) -> None:
    bad_tab = SimpleNamespace(_grid_type=GridType.SPHERICAL, gto_data=object(), original_axes=None)
    with pytest.raises(ValueError, match='grid attribute'):
        plotter_env.make_plotter(tabulator=bad_tab)


def test_plotter_rejects_unknown_grid_type(plotter_env: Any) -> None:
    bad_tab = SimpleNamespace(
        grid=np.zeros((1, 3)),
        _grid_type=GridType.UNKNOWN,
        gto_data=object(),
        original_axes=None,
    )
    with pytest.raises(ValueError, match='only supports spherical and cartesian'):
        plotter_env.make_plotter(tabulator=bad_tab)


def test_plotter_requires_tabulated_gtos(plotter_env: Any) -> None:
    tabulator = plotter_env.make_tabulator()
    del tabulator.gto_data

    with pytest.raises(ValueError, match='tabulated GTOs'):
        plotter_env.make_plotter(tabulator=tabulator)


def test_export_orbitals_dialog_sets_attributes(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    install_fake_tk_widgets(monkeypatch)
    plotter = plotter_env.make_plotter()
    plotter.selection_screen.current_mo_ind = 1

    plotter.export_orbitals_dialog()

    assert isinstance(plotter._export_window, SimpleToplevel)  # noqa: SLF001
    assert plotter._export_current_orb_radio is not None  # noqa: SLF001
    assert plotter._export_all_orb_radio is not None  # noqa: SLF001


def test_export_image_dialog_builds_controls(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    install_fake_tk_widgets(monkeypatch)
    plotter = plotter_env.make_plotter()

    plotter.export_image_dialog()


@pytest.mark.usefixtures('plotter_env')
def test_plotter_creates_internal_tk_root(monkeypatch: pytest.MonkeyPatch) -> None:
    created: list[DummyTk] = []

    def fake_tk() -> DummyTk:
        root = DummyTk()
        created.append(root)
        return root

    monkeypatch.setattr(plotter_module, 'Tabulator', RecordingTabulator)
    monkeypatch.setattr(plotter_module.tk, 'Tk', fake_tk)

    plotter = plotter_module.Plotter('dummy', only_molecule=True)

    assert created
    assert created[0].withdrawn
    assert created[0].mainloop_calls == 1
    assert plotter.selection_screen is None


def test_plotter_generates_default_spherical_grid(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    monkeypatch.setattr(plotter_module, 'Tabulator', RecordingTabulator)
    plotter_module.config.grid.default_type = 'spherical'
    plotter_module.config.grid.spherical.num_r_points = 2
    plotter_module.config.grid.spherical.num_theta_points = 2
    plotter_module.config.grid.spherical.num_phi_points = 2
    plotter_module.config.grid.min_radius = 1
    plotter_module.config.grid.max_radius_multiplier = 1

    plotter = plotter_module.Plotter('dummy', tk_root=plotter_env.make_root())

    tabulator = plotter.tabulator
    assert isinstance(tabulator, RecordingTabulator)
    assert tabulator.spherical_calls
    r, theta, phi = tabulator.spherical_calls[0]
    expected_size = 2
    assert r.size == expected_size
    assert theta.size == expected_size
    assert phi.size == expected_size


def test_plotter_generates_default_cartesian_grid(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    monkeypatch.setattr(plotter_module, 'Tabulator', RecordingTabulator)
    plotter_module.config.grid.default_type = 'cartesian'
    plotter_module.config.grid.cartesian.num_x_points = 2
    plotter_module.config.grid.cartesian.num_y_points = 3
    plotter_module.config.grid.cartesian.num_z_points = 4
    plotter_module.config.grid.max_radius_multiplier = 1
    plotter_module.config.grid.min_radius = 1

    plotter = plotter_module.Plotter('dummy', tk_root=plotter_env.make_root())

    tabulator = plotter.tabulator
    assert isinstance(tabulator, RecordingTabulator)
    assert tabulator.cartesian_calls
    x, y, z = tabulator.cartesian_calls[0]
    expected_x = 2
    expected_y = 3
    expected_z = 4
    assert x.size == expected_x
    assert y.size == expected_y
    assert z.size == expected_z


def test_plotter_initializes_custom_colormap(plotter_env: Any) -> None:
    plotter_module.config.mo.custom_colors = ['navy', 'white']
    plotter = plotter_env.make_plotter()
    assert getattr(plotter.cmap, 'name', '') == 'custom_mo'


def test_plotter_builds_menus_and_overrides_clear(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(plotter_module, 'config', plotter_module.Config())
    monkeypatch.setattr(plotter_module, 'Tabulator', FakeTabulator)
    monkeypatch.setattr(plotter_module, 'BackgroundPlotter', MenuAwareBackgroundPlotter)
    monkeypatch.setattr(plotter_module, 'Molecule', DummyMolecule)
    monkeypatch.setattr(plotter_module, '_OrbitalSelectionScreen', DummySelectionScreen)
    monkeypatch.setattr(plotter_module.pv, 'StructuredGrid', DummyStructuredGrid)
    monkeypatch.setattr(plotter_module.pv, 'pyvista_ndarray', lambda arr: arr)
    monkeypatch.setattr(plotter_module, 'QMenu', FakeQMenu)
    monkeypatch.setattr(plotter_module, 'QAction', FakeQAction)
    monkeypatch.setattr(plotter_module, 'isValid', lambda _action: True)

    plotter = plotter_module.Plotter('dummy', tk_root=DummyTk())

    main_menu = plotter.pv_plotter.main_menu
    assert [menu.title for menu in main_menu.menus] == ['Settings', 'Export']
    settings_action_texts = [action.text() for action in main_menu.menus[0].actions()]
    assert {'Grid Settings', 'MO Settings', 'Molecule Settings'} <= set(settings_action_texts)
    assert main_menu.clear_action.triggered.callbacks[-1] == plotter._clear_all  # noqa: SLF001


def test_plotter_validates_axis_spacing(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    tabulator = plotter_env.make_tabulator()
    axes = (
        np.array([0.0, 0.5, 1.0]),
        np.array([0.0, 0.5, 1.0]),
        np.array([0.0, 0.5, 1.0]),
    )
    tabulator.original_axes = axes

    calls: list[str] = []
    monkeypatch.setattr(
        plotter_module.Tabulator,
        '_axis_spacing',
        lambda _axis, name: calls.append(name) or 0.5,
    )

    plotter_env.make_plotter(tabulator=tabulator)
    assert calls == ['x', 'y', 'z']


def test_plotter_only_molecule_skips_selection_screen(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter(only_molecule=True)
    assert not hasattr(plotter, 'orb_mesh')
    assert plotter.selection_screen is None


def test_apply_grid_settings_updates_spherical_grid(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.grid_type_radio_var = DummyVar(GridType.SPHERICAL.value)
    plotter.radius_entry = DummyEntry('2.0')
    plotter.radius_points_entry = DummyEntry('2')
    plotter.theta_points_entry = DummyEntry('2')
    plotter.phi_points_entry = DummyEntry('2')
    plotter.selection_screen.current_mo_ind = 0

    captured: dict[str, Any] = {}

    def fake_update(i_points: np.ndarray, j_points: np.ndarray, k_points: np.ndarray, grid_type: Any) -> None:
        captured['args'] = (i_points, j_points, k_points, grid_type)

    replotted: list[int] = []

    def remember(idx: int) -> None:
        replotted.append(idx)

    plotter.update_mesh = fake_update  # type: ignore[assignment]
    plotter.plot_orbital = remember  # type: ignore[assignment]

    plotter.apply_grid_settings()

    i_points, j_points, k_points, grid_type = captured['args']
    assert grid_type == GridType.SPHERICAL
    expected_points = 2
    assert i_points.shape[0] == expected_points
    assert j_points.shape[0] == expected_points
    assert k_points.shape[0] == expected_points
    assert replotted == [0]


def test_apply_grid_settings_cartesian_validation_shows_error(
    monkeypatch: pytest.MonkeyPatch,
    plotter_env: Any,
) -> None:
    plotter = plotter_env.make_plotter()
    plotter.grid_type_radio_var = DummyVar(GridType.CARTESIAN.value)
    plotter.x_min_entry = DummyEntry('-1.0')
    plotter.x_max_entry = DummyEntry('1.0')
    plotter.x_num_points_entry = DummyEntry('0')
    plotter.y_min_entry = DummyEntry('-1.0')
    plotter.y_max_entry = DummyEntry('1.0')
    plotter.y_num_points_entry = DummyEntry('5')
    plotter.z_min_entry = DummyEntry('-1.0')
    plotter.z_max_entry = DummyEntry('1.0')
    plotter.z_num_points_entry = DummyEntry('5')

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(
        plotter_module.messagebox,
        'showerror',
        lambda title, msg: errors.append((title, msg)),
    )

    plotter.apply_grid_settings()
    assert errors


def test_apply_grid_settings_rejects_nonpositive_radius(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.grid_type_radio_var = DummyVar(GridType.SPHERICAL.value)
    plotter.radius_entry = DummyEntry('0')
    plotter.radius_points_entry = DummyEntry('1')
    plotter.theta_points_entry = DummyEntry('1')
    plotter.phi_points_entry = DummyEntry('1')

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(plotter_module.messagebox, 'showerror', lambda title, msg: errors.append((title, msg)))

    plotter.apply_grid_settings()
    assert errors
    assert 'Radius' in errors[0][1]


def test_apply_grid_settings_rejects_nonpositive_points(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.grid_type_radio_var = DummyVar(GridType.SPHERICAL.value)
    plotter.radius_entry = DummyEntry('1.0')
    plotter.radius_points_entry = DummyEntry('0')
    plotter.theta_points_entry = DummyEntry('1')
    plotter.phi_points_entry = DummyEntry('1')

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(plotter_module.messagebox, 'showerror', lambda title, msg: errors.append((title, msg)))

    plotter.apply_grid_settings()
    assert errors
    assert 'Number of points' in errors[0][1]
    assert 'greater than zero' in errors[0][1]


def test_grid_settings_screen_creates_entries(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    install_fake_tk_widgets(monkeypatch)
    plotter = plotter_env.make_plotter()
    seed_tabulator_with_cartesian_grid(plotter.tabulator)

    plotter.grid_settings_screen()

    assert hasattr(plotter, 'radius_entry')
    assert hasattr(plotter, 'x_min_entry')


def test_mo_settings_screen_initializes_controls(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    install_fake_tk_widgets(monkeypatch)
    plotter = plotter_env.make_plotter()

    plotter.mo_settings_screen()

    assert hasattr(plotter, 'contour_entry')
    assert hasattr(plotter, 'opacity_scale')


def test_molecule_settings_screen_initializes_entries(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    install_fake_tk_widgets(monkeypatch)
    plotter = plotter_env.make_plotter()

    plotter.molecule_settings_screen()

    assert hasattr(plotter, 'bond_max_length_entry')
    assert hasattr(plotter, 'bond_radius_entry')


def test_color_settings_screen_initializes_entries(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    install_fake_tk_widgets(monkeypatch)
    plotter = plotter_env.make_plotter()

    plotter.color_settings_screen()

    assert isinstance(plotter.background_color_entry, SimpleEntry)


def test_color_settings_screen_populates_custom_scheme(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    install_fake_tk_widgets(monkeypatch)
    plotter_module.config.mo.color_scheme = 'bespoke'
    plotter_module.config.mo.custom_colors = ['black', 'white']
    plotter = plotter_env.make_plotter()

    plotter.color_settings_screen()

    assert plotter.mo_color_scheme_var.get() == 'bespoke'
    assert plotter.mo_negative_color_entry.get() == 'black'
    assert plotter.mo_positive_color_entry.get() == 'white'


def test_on_mo_color_scheme_change_toggles_widgets(plotter_env: Any) -> None:
    class Tracker(SimpleWidget):
        def __init__(self) -> None:
            super().__init__()
            self.visible = False

        def grid(self, *_args: Any, **_kwargs: Any) -> SimpleWidget:
            self.visible = True
            return super().grid(*_args, **_kwargs)

        def grid_remove(self) -> SimpleWidget:
            self.visible = False
            return self

    plotter = plotter_env.make_plotter()
    plotter.mo_custom_color_widgets = [Tracker()]
    plotter.mo_color_scheme_var = DummyVar('custom')

    plotter.on_mo_color_scheme_change(SimpleNamespace())
    assert plotter.mo_custom_color_widgets[0].visible

    plotter.mo_color_scheme_var.set('coolwarm')
    plotter.on_mo_color_scheme_change(SimpleNamespace())
    assert not plotter.mo_custom_color_widgets[0].visible


def test_on_bond_color_type_change_updates_visibility(plotter_env: Any) -> None:
    class Tracker(SimpleWidget):
        def __init__(self) -> None:
            super().__init__()
            self.visible = False

        def grid(self, *_args: Any, **_kwargs: Any) -> SimpleWidget:
            self.visible = True
            return self

        def grid_remove(self) -> SimpleWidget:
            self.visible = False
            return self

    plotter = plotter_env.make_plotter()
    plotter.bond_color_type_var = DummyVar('uniform')
    plotter.bond_color_label = Tracker()
    plotter.bond_color_entry = Tracker()

    calls: list[str] = []
    plotter.apply_bond_color_settings = lambda: calls.append('run')  # type: ignore[assignment]

    plotter.on_bond_color_type_change()
    assert plotter.bond_color_label.visible
    plotter.bond_color_type_var.set('gradient')
    plotter.on_bond_color_type_change()
    assert not plotter.bond_color_label.visible
    assert calls


def test_reset_grid_settings_restores_defaults(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    install_fake_tk_widgets(monkeypatch)
    plotter = plotter_env.make_plotter()
    plotter.grid_settings_screen()
    plotter.radius_entry.insert(0, '9.0')

    plotter.reset_grid_settings()
    expected = max(
        plotter_module.config.grid.max_radius_multiplier * plotter.molecule.max_radius,
        plotter_module.config.grid.min_radius,
    )
    assert float(plotter.radius_entry.get()) == pytest.approx(expected)


def test_reset_mo_settings_restores_inputs(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    install_fake_tk_widgets(monkeypatch)
    plotter = plotter_env.make_plotter()
    plotter.mo_settings_screen()
    plotter.contour_entry.insert(0, '9.5')
    plotter.opacity_scale.value = 0.1

    plotter.reset_mo_settings()
    assert float(plotter.contour_entry.get()) == pytest.approx(plotter_module.config.mo.contour)
    assert plotter.opacity_scale.value == pytest.approx(plotter_module.config.mo.opacity)


def test_reset_molecule_settings_restores_config(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    install_fake_tk_widgets(monkeypatch)
    plotter = plotter_env.make_plotter()
    plotter.molecule_settings_screen()
    plotter.atom_actors[0].SetVisibility(False)
    plotter.bond_actors[0].SetVisibility(False)
    plotter.bond_max_length_entry.insert(0, '7.7')

    plotter.reset_molecule_settings()
    assert plotter.are_atoms_visible()
    assert plotter.are_bonds_visible()
    assert float(plotter.bond_max_length_entry.get()) == pytest.approx(
        plotter_module.config.molecule.bond.max_length,
    )


def test_reset_color_settings_restores_entries(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    install_fake_tk_widgets(monkeypatch)
    plotter_module.config.molecule.bond.color_type = 'uniform'
    plotter_module.config.molecule.bond.color = 'yellow'
    plotter = plotter_env.make_plotter()
    plotter.color_settings_screen()
    plotter.background_color_entry.insert(0, 'black')
    plotter.mo_color_scheme_var.set('custom')
    plotter.mo_negative_color_entry.insert(0, 'purple')
    plotter.bond_color_type_var.set('split')

    plotter.reset_color_settings()
    assert plotter.background_color_entry.get() == str(plotter_module.config.background_color)
    assert plotter.bond_color_type_var.get() == plotter_module.config.molecule.bond.color_type
    assert plotter.bond_color_entry.get() == str(plotter_module.config.molecule.bond.color)


def test_apply_mo_contour_replots_current(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.contour_entry = DummyEntry('0.25')
    plotter.selection_screen.current_mo_ind = 1

    replotted: list[int] = []

    def remember(idx: int) -> None:
        replotted.append(idx)

    plotter.plot_orbital = remember  # type: ignore[assignment]

    plotter.apply_mo_contour()
    assert plotter.contour == pytest.approx(0.25)
    assert replotted == [1]


def test_apply_custom_mo_color_settings_updates_scheme(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.mo_color_scheme_var = DummyVar('custom')
    plotter.mo_negative_color_entry = DummyEntry('navy')
    plotter.mo_positive_color_entry = DummyEntry('gold')
    plotter.selection_screen.current_mo_ind = 0

    replotted: list[int] = []

    def remember(idx: int) -> None:
        replotted.append(idx)

    plotter.plot_orbital = remember  # type: ignore[assignment]

    plotter.apply_custom_mo_color_settings()

    assert plotter_module.config.mo.custom_colors == ['navy', 'gold']
    assert replotted == [0]


def test_apply_custom_mo_color_settings_rejects_invalid(
    monkeypatch: pytest.MonkeyPatch,
    plotter_env: Any,
) -> None:
    plotter = plotter_env.make_plotter()
    plotter.mo_color_scheme_var = DummyVar('custom')
    plotter.mo_negative_color_entry = DummyEntry('navy')
    plotter.mo_positive_color_entry = DummyEntry('not-a-color')

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(
        plotter_module.messagebox,
        'showerror',
        lambda title, msg: errors.append((title, msg)),
    )

    plotter.apply_custom_mo_color_settings()
    assert errors
    assert 'custom colors' in errors[0][1]


def test_apply_mo_color_settings_switches_scheme(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.mo_color_scheme_var = DummyVar('viridis')
    plotter.selection_screen.current_mo_ind = 2
    plotter.on_mo_color_scheme_change = lambda *_args: None  # type: ignore[assignment]

    replotted: list[int] = []

    def remember(idx: int) -> None:
        replotted.append(idx)

    plotter.plot_orbital = remember  # type: ignore[assignment]

    plotter.apply_mo_color_settings()

    assert plotter_module.config.mo.color_scheme == 'viridis'
    assert plotter.cmap == 'viridis'
    assert replotted == [2]


def test_apply_color_settings_runs_all_handlers(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    calls: list[str] = []
    plotter.apply_mo_color_settings = lambda: calls.append('mo')  # type: ignore[assignment]
    plotter.apply_custom_mo_color_settings = lambda: calls.append('custom')  # type: ignore[assignment]
    plotter.apply_bond_color_settings = lambda: calls.append('bond')  # type: ignore[assignment]

    plotter.apply_color_settings()
    assert calls == ['mo', 'custom', 'bond']


def test_apply_mo_color_settings_returns_for_custom(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.mo_color_scheme_var = DummyVar('custom')
    plotter.on_mo_color_scheme_change = lambda *_args: None  # type: ignore[assignment]

    plotter.apply_mo_color_settings()
    assert plotter.cmap == plotter_module.config.mo.color_scheme


def test_apply_mo_color_settings_same_value_no_replot(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    current_scheme = plotter_module.config.mo.color_scheme
    plotter.mo_color_scheme_var = DummyVar(current_scheme)
    plotter.on_mo_color_scheme_change = lambda *_args: None  # type: ignore[assignment]

    calls: list[int] = []

    def remember_plot(idx: int) -> None:
        calls.append(idx)

    plotter.plot_orbital = remember_plot  # type: ignore[assignment]
    plotter.selection_screen.current_mo_ind = 1

    plotter.apply_mo_color_settings()
    assert not calls


def test_apply_mo_color_settings_skips_replot_without_selection(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.mo_color_scheme_var = DummyVar('viridis')
    plotter.selection_screen.current_mo_ind = -1
    plotter.on_mo_color_scheme_change = lambda *_args: None  # type: ignore[assignment]

    calls: list[int] = []

    def remember_plot(idx: int) -> None:
        calls.append(idx)

    plotter.plot_orbital = remember_plot  # type: ignore[assignment]

    plotter.apply_mo_color_settings()
    assert not calls


def test_apply_custom_mo_color_settings_non_custom_returns(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.mo_color_scheme_var = DummyVar('coolwarm')

    plotter.apply_custom_mo_color_settings()
    assert plotter.cmap == plotter_module.config.mo.color_scheme


def test_apply_background_color_updates_plotter(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.background_color_entry = DummyEntry('navy')
    plotter.background_color_var = DummyVar('navy')

    plotter.apply_background_color()

    assert plotter.pv_plotter.background == 'navy'


def test_apply_background_color_rejects_invalid(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.background_color_entry = DummyEntry('not_a_color')
    plotter.background_color_var = DummyVar('not_a_color')

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(plotter_module.messagebox, 'showerror', lambda title, msg: errors.append((title, msg)))

    plotter.apply_background_color()
    assert errors


def test_apply_background_color_handles_exception(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.background_color_entry = DummyEntry('navy')

    def boom(_color: str) -> None:
        raise ValueError('explode')

    plotter.pv_plotter.set_background = boom  # type: ignore[assignment]

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(plotter_module.messagebox, 'showerror', lambda title, msg: errors.append((title, msg)))

    plotter.apply_background_color()
    assert errors


def test_apply_bond_color_settings_reloads_molecule(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.bond_color_type_var = DummyVar('split')
    plotter.bond_color_entry = DummyEntry('red')

    reloads: list[Any] = []

    def remember(cfg: Any) -> None:
        reloads.append(cfg)

    plotter.load_molecule = remember  # type: ignore[assignment]

    plotter.apply_bond_color_settings()

    assert plotter_module.config.molecule.bond.color_type == 'split'
    assert reloads


def test_apply_bond_color_settings_no_reload_when_unchanged(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    cfg = plotter_module.config.molecule.bond
    plotter.bond_color_type_var = DummyVar(cfg.color_type)
    plotter.bond_color_entry = DummyEntry(cfg.color)

    reloads: list[Any] = []

    def remember(cfg: Any) -> None:
        reloads.append(cfg)

    plotter.load_molecule = remember  # type: ignore[assignment]

    plotter.apply_bond_color_settings()
    assert not reloads


def test_apply_bond_color_settings_updates_uniform_color(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.bond_color_type_var = DummyVar('uniform')
    plotter.bond_color_entry = DummyEntry('cyan')

    reloads: list[Any] = []

    def remember(cfg: Any) -> None:
        reloads.append(cfg)

    plotter.load_molecule = remember  # type: ignore[assignment]

    plotter.apply_bond_color_settings()
    assert plotter_module.config.molecule.bond.color == 'cyan'
    assert reloads


def test_on_opacity_change_updates_actor(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    monkeypatch.setattr(plotter_module.ttk, 'Label', DummyLabelWidget)
    plotter = plotter_env.make_plotter()
    plotter.orb_actor = DummyActor()
    label = DummyLabelWidget('Molecular Orbital Opacity: 1.00')
    container = DummyContainer([label])
    plotter.mo_settings_window = DummyContainer([container])

    plotter.on_opacity_change('0.33')

    assert plotter.opacity == pytest.approx(0.33)
    assert plotter.orb_actor.opacity == pytest.approx(0.33)
    assert '0.33' in label.text


def test_on_molecule_opacity_change_updates_actors(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    monkeypatch.setattr(plotter_module.ttk, 'Label', DummyLabelWidget)
    plotter = plotter_env.make_plotter()
    label = DummyLabelWidget('Molecule Opacity: 1.00')
    plotter.molecule_settings_window = DummyContainer([DummyContainer([label])])

    plotter.on_molecule_opacity_change('0.45')

    assert plotter.molecule_opacity == pytest.approx(0.45)
    assert all(actor.opacity == pytest.approx(0.45) for actor in plotter.molecule_actors)
    assert '0.45' in label.text


def test_clear_all_hides_actors_and_resets_selection(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.selection_screen.current_mo_ind = 3
    plotter.orb_actor = DummyActor()

    plotter._clear_all()  # noqa: SLF001

    assert plotter.orb_actor is None
    assert plotter.selection_screen.current_mo_ind == -1
    assert all(not actor.GetVisibility() for actor in plotter.molecule_actors)


def test_toggle_atoms_flips_visibility(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    initial_visibility = plotter.atom_actors[0].GetVisibility()

    plotter.toggle_atoms()

    assert plotter.atom_actors[0].GetVisibility() is (not initial_visibility)


def test_apply_grid_settings_updates_cartesian(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.grid_type_radio_var = DummyVar(GridType.CARTESIAN.value)
    plotter.x_min_entry = DummyEntry('-1.0')
    plotter.x_max_entry = DummyEntry('1.0')
    plotter.x_num_points_entry = DummyEntry('2')
    plotter.y_min_entry = DummyEntry('-2.0')
    plotter.y_max_entry = DummyEntry('2.0')
    plotter.y_num_points_entry = DummyEntry('3')
    plotter.z_min_entry = DummyEntry('-3.0')
    plotter.z_max_entry = DummyEntry('3.0')
    plotter.z_num_points_entry = DummyEntry('4')
    plotter.selection_screen.current_mo_ind = 0

    captured: dict[str, Any] = {}

    def fake_update(i_points: np.ndarray, j_points: np.ndarray, k_points: np.ndarray, grid_type: Any) -> None:
        captured['args'] = (i_points, j_points, k_points, grid_type)

    replotted: list[int] = []

    def remember(idx: int) -> None:
        replotted.append(idx)

    plotter.update_mesh = fake_update  # type: ignore[assignment]
    plotter.plot_orbital = remember  # type: ignore[assignment]

    plotter.apply_grid_settings()

    i_points, j_points, k_points, grid_type = captured['args']
    assert grid_type == GridType.CARTESIAN
    assert np.isclose(i_points[[0, -1]], [-1.0, 1.0]).all()
    expected_j_points = 3
    expected_k_points = 4
    assert j_points.size == expected_j_points
    assert k_points.size == expected_k_points
    assert replotted == [0]


def test_apply_molecule_settings_updates_config(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.bond_max_length_entry = DummyEntry('5.5')
    plotter.bond_radius_entry = DummyEntry('0.22')

    reloads: list[Any] = []

    def remember(cfg: Any) -> None:
        reloads.append(cfg)

    plotter.load_molecule = remember  # type: ignore[assignment]

    plotter.apply_molecule_settings()

    new_max = 5.5
    new_radius = 0.22
    assert pytest.approx(plotter_module.config.molecule.bond.max_length) == new_max
    assert pytest.approx(plotter_module.config.molecule.bond.radius) == new_radius
    assert reloads


def test_apply_molecule_settings_validates_numbers(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.bond_max_length_entry = DummyEntry('bad-number')
    plotter.bond_radius_entry = DummyEntry('0.22')

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(
        plotter_module.messagebox,
        'showerror',
        lambda title, msg: errors.append((title, msg)),
    )

    plotter.apply_molecule_settings()
    assert errors
    assert 'Bond Max Length' in errors[0][1]


def test_apply_molecule_settings_rejects_invalid_radius(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.bond_max_length_entry = DummyEntry('1.0')
    plotter.bond_radius_entry = DummyEntry('not-number')

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(plotter_module.messagebox, 'showerror', lambda title, msg: errors.append((title, msg)))

    plotter.apply_molecule_settings()
    assert errors
    assert 'Bond Radius' in errors[0][1]


def test_apply_molecule_settings_no_changes_skip_reload(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    cfg = plotter_module.config
    plotter.bond_max_length_entry = DummyEntry(str(cfg.molecule.bond.max_length))
    plotter.bond_radius_entry = DummyEntry(str(cfg.molecule.bond.radius))

    reloads: list[Any] = []

    def remember(cfg: Any) -> None:
        reloads.append(cfg)

    plotter.load_molecule = remember  # type: ignore[assignment]

    plotter.apply_molecule_settings()
    assert not reloads


def test_load_molecule_removes_existing_actors(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    original_actor = plotter.molecule_actors[0]

    plotter.load_molecule(plotter_module.config)
    assert original_actor in plotter.pv_plotter.removed_actors


def test_apply_mo_contour_invalid_input(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.contour_entry = DummyEntry('not-a-number')
    plotter.selection_screen.current_mo_ind = 0

    replotted: list[int] = []

    def remember(idx: int) -> None:
        replotted.append(idx)

    plotter.plot_orbital = remember  # type: ignore[assignment]

    plotter.apply_mo_contour()
    assert replotted == []
    assert plotter.contour == pytest.approx(plotter_module.config.mo.contour)


def test_do_image_export_vector(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    window = DummyWindow()
    format_var = DummyVar('svg')
    transparent_var = DummyVar(False)

    infos: list[tuple[str, str]] = []
    monkeypatch.setattr(plotter_module.filedialog, 'asksaveasfilename', lambda **_kwargs: '/tmp/export.svg')
    monkeypatch.setattr(plotter_module.messagebox, 'showinfo', lambda title, msg: infos.append((title, msg)))
    monkeypatch.setattr(
        plotter_module.messagebox,
        'showerror',
        lambda *_args, **_kwargs: pytest.fail('Unexpected error'),
    )

    plotter._do_image_export(window, format_var, transparent_var)  # noqa: SLF001

    assert plotter.pv_plotter.saved_graphic == '/tmp/export.svg'
    assert window.destroyed
    assert infos


def test_do_image_export_png_respects_transparency(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    window = DummyWindow()
    format_var = DummyVar('png')
    transparent_var = DummyVar(True)

    monkeypatch.setattr(plotter_module.filedialog, 'asksaveasfilename', lambda **_kwargs: '/tmp/export.png')
    monkeypatch.setattr(plotter_module.messagebox, 'showinfo', lambda *_args, **_kwargs: None)
    monkeypatch.setattr(plotter_module.messagebox, 'showerror', lambda *_args, **_kwargs: None)

    plotter._do_image_export(window, format_var, transparent_var)  # noqa: SLF001

    assert plotter.pv_plotter.screenshot_calls == [('/tmp/export.png', True)]


def test_do_image_export_handles_failures(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    window = DummyWindow()
    format_var = DummyVar('png')
    transparent_var = DummyVar(False)

    monkeypatch.setattr(plotter_module.filedialog, 'asksaveasfilename', lambda **_kwargs: '/tmp/export.png')

    def boom(*_args: Any, **_kwargs: Any) -> None:
        raise RuntimeError('boom')

    plotter.pv_plotter.screenshot = boom  # type: ignore[assignment]

    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(plotter_module.messagebox, 'showerror', lambda title, msg: errors.append((title, msg)))
    monkeypatch.setattr(plotter_module.messagebox, 'showinfo', lambda *_args, **_kwargs: None)

    plotter._do_image_export(window, format_var, transparent_var)  # noqa: SLF001
    assert errors
    assert 'Failed to export image' in errors[0][1]


def test_do_image_export_cancel(monkeypatch: pytest.MonkeyPatch, plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    window = DummyWindow()
    format_var = DummyVar('png')
    transparent_var = DummyVar(False)

    monkeypatch.setattr(plotter_module.filedialog, 'asksaveasfilename', lambda **_kwargs: '')

    plotter._do_image_export(window, format_var, transparent_var)  # noqa: SLF001
    assert not plotter.pv_plotter.screenshot_calls


def test_save_settings_reports_success(monkeypatch: pytest.MonkeyPatch) -> None:
    saves: list[str] = []
    monkeypatch.setattr(plotter_module.config, 'save_current_config', lambda: saves.append('saved'))

    infos: list[tuple[str, str]] = []
    monkeypatch.setattr(plotter_module.messagebox, 'showinfo', lambda title, msg: infos.append((title, msg)))

    plotter_module.Plotter.save_settings()
    assert saves == ['saved']
    assert infos


def test_save_settings_reports_error(monkeypatch: pytest.MonkeyPatch) -> None:
    def boom() -> None:
        raise OSError('disk full')

    monkeypatch.setattr(plotter_module.config, 'save_current_config', boom)
    errors: list[tuple[str, str]] = []
    monkeypatch.setattr(plotter_module.messagebox, 'showerror', lambda title, msg: errors.append((title, msg)))

    plotter_module.Plotter.save_settings()
    assert errors


def test_plot_orbital_creates_actor(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.selection_screen.current_mo_ind = -1

    plotter.plot_orbital(0)

    assert plotter.selection_screen.current_mo_ind == 0
    assert isinstance(plotter.orb_actor, DummyActor)
    assert 'orbital' in plotter.orb_mesh.arrays
    assert plotter.pv_plotter.added_meshes


def test_plot_orbital_minus_one_clears_scene(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.plot_orbital(0)

    plotter.plot_orbital(-1)

    assert plotter.orb_actor is None
    assert plotter.selection_screen.current_mo_ind == -1
    assert plotter.pv_plotter.removed_actors


def test_toggle_bonds_triggers_update(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    initial_visibility = plotter.bond_actors[0].GetVisibility()

    plotter.toggle_bonds()

    assert plotter.bond_actors[0].GetVisibility() is (not initial_visibility)
    assert plotter.pv_plotter.update_count == 1


def test_update_mesh_rebuilds_structured_grid(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    x = np.linspace(-1, 1, 2)
    y = np.linspace(-1, 1, 2)
    z = np.linspace(-1, 1, 2)

    plotter.update_mesh(x, y, z, GridType.CARTESIAN)

    expected_points = x.size * y.size * z.size
    assert plotter.tabulator.grid.shape[0] == expected_points
    assert plotter.orb_mesh.points.shape[0] == expected_points


def test_update_mesh_rejects_unknown_grid_type(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    points = np.linspace(0.0, 1.0, 2)

    with pytest.raises(ValueError, match='only supports spherical'):
        plotter.update_mesh(points, points, points, GridType.UNKNOWN)


def test_update_mesh_handles_spherical_grid(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    r = np.linspace(0.0, 1.0, 2)

    plotter.update_mesh(r, r, r, GridType.SPHERICAL)
    assert plotter.tabulator._grid_type == GridType.SPHERICAL  # noqa: SLF001


def test_toggle_molecule_balances_atom_and_bond_visibility(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.atom_actors[0].SetVisibility(False)
    plotter.bond_actors[0].SetVisibility(True)
    plotter.show_atoms_var = DummyVar(True)
    plotter.show_bonds_var = DummyVar(True)

    plotter.toggle_molecule()

    assert plotter.are_atoms_visible()
    assert plotter.are_bonds_visible()
    assert plotter_module.config.molecule.atom.show is plotter.are_atoms_visible()
    assert plotter_module.config.molecule.bond.show is plotter.are_bonds_visible()


def test_toggle_molecule_toggles_scene_when_states_match(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    visible_before = plotter.molecule_actors[0].GetVisibility()

    plotter.toggle_molecule()

    assert plotter.molecule_actors[0].GetVisibility() is (not visible_before)
    assert plotter.pv_plotter.update_count == 1


def test_visibility_helpers_reflect_actor_state(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    assert plotter.is_molecule_visible()
    assert plotter.are_atoms_visible()
    assert plotter.are_bonds_visible()

    plotter.molecule_actors[0].SetVisibility(False)
    plotter.atom_actors[0].SetVisibility(False)
    plotter.bond_actors[0].SetVisibility(False)

    assert not plotter.is_molecule_visible()
    assert not plotter.are_atoms_visible()
    assert not plotter.are_bonds_visible()

    plotter.molecule_actors = []
    plotter.atom_actors = []
    plotter.bond_actors = []
    assert plotter.is_molecule_visible() is False
    assert plotter.are_atoms_visible() is False
    assert plotter.are_bonds_visible() is False


def test_update_settings_button_states(plotter_env: Any) -> None:
    plotter = plotter_env.make_plotter()
    plotter.show_atoms_var = DummyVar(False)
    plotter.show_bonds_var = DummyVar(False)
    plotter.atom_actors[0].SetVisibility(False)
    plotter.bond_actors[0].SetVisibility(True)

    plotter.update_settings_button_states()

    assert plotter.show_atoms_var.get() is False
    assert plotter.show_bonds_var.get() is True
    assert plotter_module.config.molecule.atom.show is False
    assert plotter_module.config.molecule.bond.show is True


@pytest.mark.usefixtures('plotter_env')
def test_pv_plotter_close_signal_closes_windows(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(plotter_module, 'Tabulator', RecordingTabulator)
    monkeypatch.setattr(plotter_module.tk, 'Tk', DummyTk)
    plotter = plotter_module.Plotter('dummy')
    callbacks = plotter.pv_plotter.app_window.signal_close.callbacks
    assert callbacks

    callbacks[0]()

    assert not plotter.on_screen
    assert plotter.selection_screen.destroyed
    assert plotter.tk_root.quit_calls == 1


def install_fake_tk_widgets(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(plotter_module.tk, 'Toplevel', SimpleToplevel)
    monkeypatch.setattr(plotter_module.tk, 'StringVar', DummyVar)
    monkeypatch.setattr(plotter_module.tk, 'BooleanVar', DummyVar)

    widget_names = [
        'Frame',
        'Label',
        'Radiobutton',
        'Button',
        'Scale',
        'Checkbutton',
        'Separator',
    ]
    for name in widget_names:
        monkeypatch.setattr(plotter_module.ttk, name, SimpleWidget)
    monkeypatch.setattr(plotter_module.ttk, 'Entry', SimpleEntry)
    monkeypatch.setattr(plotter_module.ttk, 'Combobox', SimpleCombobox)


@pytest.fixture
def selection_screen_ui(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    original_screen_bases = plotter_module._OrbitalSelectionScreen.__bases__  # noqa: SLF001
    original_tree_bases = plotter_module._OrbitalsTreeview.__bases__  # noqa: SLF001
    plotter_module._OrbitalSelectionScreen.__bases__ = (SimpleToplevel,)  # noqa: SLF001
    plotter_module._OrbitalsTreeview.__bases__ = (SimpleTreeview,)  # noqa: SLF001
    monkeypatch.setattr(plotter_module.ttk, 'Frame', SimpleWidget)
    monkeypatch.setattr(plotter_module.ttk, 'Button', SimpleWidget)
    monkeypatch.setattr(plotter_module.ttk, 'Label', SimpleWidget)
    yield
    plotter_module._OrbitalSelectionScreen.__bases__ = original_screen_bases  # noqa: SLF001
    plotter_module._OrbitalsTreeview.__bases__ = original_tree_bases  # noqa: SLF001


def test_orbital_selection_screen_navigation_and_close(selection_screen_ui: None) -> None:
    _ = selection_screen_ui
    plotter = SelectionPlotter()
    screen = plotter_module._OrbitalSelectionScreen(plotter)  # noqa: SLF001
    plotter.selection_screen = screen

    screen._export_current_orb_radio = SimpleWidget()  # noqa: SLF001
    screen.update_nav_button_states()
    assert screen.prev_button.kwargs['state'] == plotter_module.tk.DISABLED
    assert screen.next_button.kwargs['state'] == plotter_module.tk.NORMAL
    assert 'None' in screen._export_current_orb_radio.kwargs['text']  # noqa: SLF001
    assert screen._export_current_orb_radio.kwargs['state'] == plotter_module.tk.DISABLED  # noqa: SLF001

    screen.current_mo_ind = 0
    screen.update_nav_button_states()
    assert screen.next_button.kwargs['state'] == plotter_module.tk.NORMAL
    assert '#1' in screen._export_current_orb_radio.kwargs['text']  # noqa: SLF001

    screen.current_mo_ind = -1
    screen.next_plot()
    screen.next_plot()
    assert plotter.plot_calls[:2] == [0, 1]
    assert screen.orb_tv.current_mo_ind == 1

    screen.prev_plot()
    assert plotter.plot_calls[-1] == 0
    screen.prev_plot()  # Should be no-op when at the start
    assert plotter.plot_calls[-1] == 0

    plotter.tabulator._parser.mos = []  # type: ignore[attr-defined]  # noqa: SLF001
    before = plotter.plot_calls.copy()
    screen.current_mo_ind = -1
    screen.next_plot()
    assert plotter.plot_calls == before

    screen.plot_orbital(0)
    assert screen.current_mo_ind == 0

    screen.on_close()
    assert not plotter.on_screen
    assert plotter.pv_plotter.closed
    assert screen.destroyed
    assert plotter.tk_root.quit_calls == 1
    assert plotter.tk_root.destroy_calls == 1


def test_orbitals_treeview_populate_and_select(selection_screen_ui: None) -> None:
    _ = selection_screen_ui
    plot_calls: list[int] = []
    updates: list[int] = []

    def record_plot(idx: int) -> None:
        plot_calls.append(idx)

    def record_update() -> None:
        updates.append(1)

    selection_screen = SimpleNamespace(
        current_mo_ind=-1,
        plot_orbital=record_plot,
        update_nav_button_states=record_update,
        _loading=False,
    )

    tree = plotter_module._OrbitalsTreeview(selection_screen)  # noqa: SLF001
    mos = [
        SimpleNamespace(sym='s', occ=2.0, energy=-0.5),
        SimpleNamespace(sym='s', occ=1.0, energy=-0.3),
        SimpleNamespace(sym='p', occ=1.0, energy=-0.1),
    ]
    tree.populate_tree(mos)
    expected_children = 3
    assert len(tree.get_children()) == expected_children
    assert tree._items[1]['values'][1].startswith('s.')  # noqa: SLF001

    tree.current_mo_ind = 1
    tree.highlight_orbital(2)
    assert tree._items[1]['tags'] == ('!hightlight',)  # noqa: SLF001
    assert tree._items[2]['tags'] == ('highlight',)  # noqa: SLF001
    expected_seen = 2
    assert tree.seen == expected_seen

    tree.set_selection(1)
    tree.on_select(SimpleNamespace())
    assert selection_screen.current_mo_ind == 1
    assert plot_calls == [1]
    assert updates
    assert tree.current_mo_ind == 1
    assert tree.selection() == ()

    tree.erase()
    assert tree.get_children() == []
