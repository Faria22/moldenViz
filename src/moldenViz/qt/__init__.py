"""Qt-native embeddable molecular-orbital viewer.

Implementation lives in focused viewer, panels, workers, and export modules.
"""

from importlib import import_module
from typing import TYPE_CHECKING, Any

from .viewer import OrbitalViewer

if TYPE_CHECKING:
    from moldenViz._config_module import MainConfig

    ViewerConfig = MainConfig

__all__ = ['OrbitalViewer', 'ViewerConfig']


def __getattr__(name: str) -> Any:
    """Load the public configuration model only when it is requested.

    Returns
    -------
    Any
        Requested module attribute.
    """
    if name == 'ViewerConfig':
        viewer_config = import_module('moldenViz._config_module').MainConfig
        globals()[name] = viewer_config
        return viewer_config
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')
