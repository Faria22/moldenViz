"""Qt delivery regressions for asynchronous viewer work."""
# ruff:file-ignore[import-private-name, undocumented-public-function]

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from threading import get_ident

from PySide6.QtCore import QCoreApplication
from PySide6.QtWidgets import QApplication

from moldenViz.qt import _CompletionDispatcher


def test_completion_dispatcher_delivers_on_qt_thread() -> None:
    application = QApplication.instance() or QApplication([])
    dispatcher = _CompletionDispatcher(application)
    owner_thread = get_ident()
    delivered: list[int] = []

    with ThreadPoolExecutor(max_workers=1) as executor:
        executor.submit(lambda: dispatcher.dispatch(lambda: delivered.append(get_ident()))).result()
    while not delivered:
        QCoreApplication.processEvents()

    assert delivered == [owner_thread]
