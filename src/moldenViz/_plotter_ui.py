"""Compatibility aliases for the former private Tk user-interface module."""

from .qt import OrbitalControlPanel

_OrbitalSelectionScreen = OrbitalControlPanel
_PlotterUI = OrbitalControlPanel

__all__: list[str] = []
