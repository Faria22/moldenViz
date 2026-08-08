"""Integration coverage for explicit per-viewer setting persistence."""
# ruff:file-ignore[import-private-name, undocumented-public-function]

from __future__ import annotations

from typing import TYPE_CHECKING

from moldenViz import _config_module as config_module

if TYPE_CHECKING:
    from pathlib import Path

    import pytest


def test_config_overrides_are_persisted_only_when_saved(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    custom_path = tmp_path / 'config.toml'
    monkeypatch.setattr(config_module, 'CUSTOM_CONFIG_PATH', custom_path)
    config = config_module.Config({'background_color': '#123456'})

    assert not custom_path.exists()
    config._save_current_config()  # ruff: ignore[private-member-access]

    saved = config_module.toml.load(custom_path)
    assert saved['background_color'] == '#123456'
