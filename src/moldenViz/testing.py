"""Supported helpers for testing Qt hosts without a rendering context."""

from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
from typing import TYPE_CHECKING

from PySide6.QtWidgets import QWidget

if TYPE_CHECKING:
    from collections.abc import Iterator

__all__ = ['NullInteractor', 'without_rendering']


class _NullActor:
    """Minimal VTK actor stand-in used by :class:`NullInteractor`."""

    def __init__(self) -> None:
        self.visible = True
        self.opacity = 1.0

    def GetVisibility(self) -> bool:  # ruff: ignore[invalid-function-name]
        """Return the recorded visibility state.

        Returns
        -------
        bool
            Recorded actor visibility.
        """
        return self.visible

    def SetVisibility(self, visible: bool) -> None:  # ruff: ignore[invalid-function-name]
        """Record a visibility change."""
        self.visible = visible

    def GetProperty(self) -> _NullActor:  # ruff: ignore[invalid-function-name]
        """Return the object that records actor properties.

        Returns
        -------
        _NullActor
            This actor stand-in.
        """
        return self

    def SetOpacity(self, opacity: float) -> None:  # ruff: ignore[invalid-function-name]
        """Record an opacity change."""
        self.opacity = opacity


class NullInteractor(QWidget):
    """Non-rendering ``QtInteractor`` replacement for UI smoke tests.

    This widget records scene operations but creates no VTK render window or
    OpenGL context. It is intended only for layout, lifecycle, and controller
    tests performed inside :func:`without_rendering`.
    """

    def __init__(self, parent: QWidget | None = None, **_kwargs: object) -> None:
        super().__init__(parent)
        self.background = ''
        self.closed_count = 0
        self.actors: list[_NullActor] = []
        self.saved_graphic: Path | None = None
        self.saved_screenshot: tuple[Path, bool] | None = None
        self.reset_count = 0

    def set_background(self, color: str) -> None:
        """Record the requested background color."""
        self.background = color

    def show_axes(self) -> None:
        """Accept an axes request without rendering it."""

    def add_mesh(self, _mesh: object, **_kwargs: object) -> _NullActor:
        """Record a mesh addition and return a minimal actor stand-in.

        Returns
        -------
        _NullActor
            Actor recording later visibility and opacity changes.
        """
        actor = _NullActor()
        self.actors.append(actor)
        return actor

    def remove_actor(self, actor: _NullActor) -> None:
        """Remove a previously recorded actor."""
        if actor in self.actors:
            self.actors.remove(actor)

    def update(self) -> None:
        """Accept a render update without drawing."""

    def close(self) -> None:  # type: ignore[override]
        """Record explicit render-window cleanup."""
        self.closed_count += 1

    def save_graphic(self, path: str | Path) -> None:
        """Record a requested vector export destination."""
        self.saved_graphic = Path(path)

    def screenshot(self, path: str | Path, *, transparent_background: bool) -> None:
        """Record a requested raster export destination."""
        self.saved_screenshot = Path(path), transparent_background

    def reset_camera(self) -> None:
        """Record a camera-reset request."""
        self.reset_count += 1


_INTERACTOR_OVERRIDE: ContextVar[type[QWidget] | None] = ContextVar('moldenviz_interactor_override', default=None)


def _get_interactor_override() -> type[QWidget] | None:
    """Return the active test interactor override, if any.

    Returns
    -------
    type[QWidget] | None
        Active override class or ``None``.
    """
    return _INTERACTOR_OVERRIDE.get()


@contextmanager
def without_rendering() -> Iterator[None]:
    """Use :class:`NullInteractor` for viewers constructed in this context.

    This allows Qt host pages to be smoke-tested under the ``offscreen`` Qt
    platform without asking VTK to create an unavailable OpenGL context.

    Yields
    ------
    None
        Control while the non-rendering interactor override is active.
    """
    token = _INTERACTOR_OVERRIDE.set(NullInteractor)
    try:
        yield
    finally:
        _INTERACTOR_OVERRIDE.reset(token)
