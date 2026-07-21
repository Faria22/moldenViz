"""Tests for molecular plotting objects."""

# ruff:file-ignore[import-private-name]

from typing import cast

import numpy as np
import pytest
import pyvista as pv

from moldenViz._config_module import Config
from moldenViz._plotting_objects import Atom, Bond


def _axis_extents(mesh: pv.PolyData, origin: np.ndarray, axis: np.ndarray) -> tuple[float, float]:
    projections = (mesh.points - origin) @ axis
    return cast(float, projections.min()), cast(float, projections.max())


def test_trim_ends_shortens_issue_58_geometry_without_boolean_operation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The reported 2.76 Angstrom N-C bond is shortened analytically."""
    atom_a = Atom(7, np.array([0.0, 0.0, 0.0]))
    atom_b = Atom(6, np.array([2.76, 0.0, 0.0]))
    bond = Bond(atom_a, atom_b)

    def fail_boolean(*_args: object, **_kwargs: object) -> None:
        pytest.fail('Bond trimming must not invoke a VTK boolean operation')

    monkeypatch.setattr(pv.PolyData, 'boolean_difference', fail_boolean)

    bond.trim_ends()

    assert isinstance(bond.mesh, pv.PolyData)
    assert _axis_extents(bond.mesh, atom_a.center, np.array([1.0, 0.0, 0.0])) == pytest.approx((0.4, 2.36))


def test_trim_ends_preserves_split_bond_colors_and_surface_endpoints() -> None:
    """Split bonds retain both colours and meet at the adjusted midpoint."""
    config = Config()
    config.molecule.bond.color_type = 'split'
    atom_a = Atom(1, np.array([0.0, 0.0, 0.0]))
    atom_b = Atom(9, np.array([0.0, 0.0, 2.0]))
    bond = Bond(atom_a, atom_b, config)

    bond.trim_ends()

    assert isinstance(bond.mesh, list)
    assert bond.color == [atom_a.atom_type.color, atom_b.atom_type.color]
    split = (2.0 + atom_a.atom_type.radius - atom_b.atom_type.radius) / 2
    axis = np.array([0.0, 0.0, 1.0])
    assert _axis_extents(bond.mesh[0], atom_a.center, axis) == pytest.approx((atom_a.atom_type.radius, split))
    assert _axis_extents(bond.mesh[1], atom_a.center, axis) == pytest.approx(
        (split, 2.0 - atom_b.atom_type.radius),
    )


def test_trim_ends_discards_bond_contained_by_overlapping_atoms(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Atoms whose radii overlap do not produce a reversed bond mesh."""
    atom_a = Atom(6, np.array([0.0, 0.0, 0.0]))
    atom_b = Atom(7, np.array([0.5, 0.0, 0.0]))
    bond = Bond(atom_a, atom_b)

    bond.trim_ends()

    assert bond.mesh is None
    assert 'Bond is entirely contained' in caplog.text
