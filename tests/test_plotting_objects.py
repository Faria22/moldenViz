"""Tests for molecular plotting objects."""

# ruff:file-ignore[import-private-name, private-member-access]

from typing import cast

import numpy as np
import pytest
import pyvista as pv

from moldenViz._config_module import Config
from moldenViz._plotting_objects import Atom, Bond, Molecule
from moldenViz.models import Atom as ParsedAtom


def _axis_extents(mesh: pv.PolyData, origin: np.ndarray, axis: np.ndarray) -> tuple[float, float]:
    projections = (mesh.points - origin) @ axis
    return cast(float, projections.min()), cast(float, projections.max())


def _assert_cap_on_atom_surface(
    mesh: pv.PolyData,
    atom: Atom,
    axis: np.ndarray,
    end: str,
) -> None:
    """Assert that one cylinder rim lies on an atom's spherical surface."""
    projections = (mesh.points - atom.center) @ axis
    cap_projection = projections.min() if end == 'start' else projections.max()
    cap_points = mesh.points[np.isclose(projections, cap_projection)]
    distances = np.linalg.norm(cap_points - atom.center, axis=1)
    np.testing.assert_allclose(distances, atom.atom_type.radius, rtol=1e-6, atol=1e-6)


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

    bond._trim_ends()

    assert isinstance(bond.mesh, pv.PolyData)
    axis = np.array([1.0, 0.0, 0.0])
    trim_distance = np.sqrt(atom_a.atom_type.radius**2 - bond.radius**2)
    assert _axis_extents(bond.mesh, atom_a.center, axis) == pytest.approx(
        (trim_distance, bond.length - trim_distance),
    )
    _assert_cap_on_atom_surface(bond.mesh, atom_a, axis, 'start')
    _assert_cap_on_atom_surface(bond.mesh, atom_b, -axis, 'start')


def test_trim_ends_preserves_split_bond_colors_and_surface_endpoints() -> None:
    """Split bonds retain both colours and meet at the adjusted midpoint."""
    config = Config()
    config.molecule.bond.color_type = 'split'
    atom_a = Atom(1, np.array([0.0, 0.0, 0.0]))
    atom_b = Atom(9, np.array([0.0, 0.0, 2.0]))
    bond = Bond(atom_a, atom_b, config)

    bond._trim_ends()

    assert isinstance(bond.mesh, list)
    assert bond.color == [atom_a.atom_type.color, atom_b.atom_type.color]
    trim_a = np.sqrt(atom_a.atom_type.radius**2 - bond.radius**2)
    trim_b = np.sqrt(atom_b.atom_type.radius**2 - bond.radius**2)
    split = (2.0 + trim_a - trim_b) / 2
    axis = np.array([0.0, 0.0, 1.0])
    assert _axis_extents(bond.mesh[0], atom_a.center, axis) == pytest.approx((trim_a, split))
    assert _axis_extents(bond.mesh[1], atom_a.center, axis) == pytest.approx(
        (split, 2.0 - trim_b),
    )
    _assert_cap_on_atom_surface(bond.mesh[0], atom_a, axis, 'start')
    _assert_cap_on_atom_surface(bond.mesh[1], atom_b, -axis, 'start')


def test_trim_ends_discards_bond_contained_by_overlapping_atoms(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Atoms whose radii overlap do not produce a reversed bond mesh."""
    atom_a = Atom(6, np.array([0.0, 0.0, 0.0]))
    atom_b = Atom(7, np.array([0.5, 0.0, 0.0]))
    bond = Bond(atom_a, atom_b)

    bond._trim_ends()

    assert bond.mesh is None
    assert 'Bond is entirely contained' in caplog.text


def test_zero_length_bond_is_discarded_before_mesh_construction(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Coincident atoms do not trigger zero-vector PyVista operations."""
    center = np.array([1.0, 2.0, 3.0])

    bond = Bond(Atom(6, center), Atom(7, center))

    assert bond.mesh is None
    assert 'Cannot render zero-length bond' in caplog.text


def test_molecule_infers_bonds_from_numpy_pairwise_distances() -> None:
    """Only atom pairs within the configured distance should be bonded."""
    atoms = [
        ParsedAtom('H', 1, np.array([0.0, 0.0, 0.0]), []),
        ParsedAtom('H', 1, np.array([1.0, 0.0, 0.0]), []),
        ParsedAtom('H', 1, np.array([10.0, 0.0, 0.0]), []),
    ]

    molecule = Molecule(atoms)

    assert len(molecule.atoms[0].bonds) == 1
    assert len(molecule.atoms[1].bonds) == 1
    assert molecule.atoms[0].bonds[0] is molecule.atoms[1].bonds[0]
    assert molecule.atoms[2].bonds == []
