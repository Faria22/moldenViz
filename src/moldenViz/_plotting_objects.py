"""Utility objects used by the Plotter to build molecular meshes."""

from __future__ import annotations

import logging
from enum import Enum
from typing import TYPE_CHECKING, cast

import numpy as np
import pyvista as pv

from ._config_module import Config
from .models import Atom as ParsedAtom
from .models import AtomType

if TYPE_CHECKING:
    from numpy.typing import NDArray

logger = logging.getLogger(__name__)

config = Config()

ATOM_TYPES = config.atom_types

# Default atom type for invalid atomic numbers
ATOM_X = AtomType(name='X', color='000000', radius=1.0, max_num_bonds=0)


class Atom:
    """Represents an atom in 3D space for visualization purposes.

    Parameters
    ----------
    atomic_number : int
        The atomic number of the element.
    center : NDArray[np.floating]
        The 3D coordinates of the atom center.
    """

    def __init__(
        self,
        atomic_number: int,
        center: NDArray[np.floating],
    ) -> None:
        """Initialize an atom for visualization.

        Parameters
        ----------
        atomic_number : int
            Atomic number that determines colour, radius, and bond limits.
        center : NDArray[np.floating]
            Cartesian coordinates of the atom centre in Angstroms.
        """
        self.atom_type = ATOM_TYPES.get(atomic_number, ATOM_X)
        if self.atom_type is ATOM_X:
            logger.warning(
                "Invalid atomic number: %d. Atom type could not be determined. Using atom 'X' instead.",
                atomic_number,
            )

        self.center = np.array(center)
        self.mesh = pv.Sphere(center=center, radius=self.atom_type.radius)
        self.bonds: list[Bond] = []

    def _remove_extra_bonds(self) -> None:
        """Clip bonds so the atom respects its configured maximum.

        Notes
        -----
        Bonds remain attached to both atoms, but the meshes are cleared for any
        discarded bonds so they are not rendered by PyVista.
        """
        if len(self.bonds) <= self.atom_type.max_num_bonds:
            return

        self.bonds.sort(key=lambda x: x.length)

        for bond in self.bonds[self.atom_type.max_num_bonds :]:
            bond.mesh = None


class Bond:
    """Represents a chemical bond between two atoms for visualization.

    Parameters
    ----------
    atom_a : Atom
        The first atom in the bond.
    atom_b : Atom
        The second atom in the bond.
    """

    class ColorType(Enum):
        """Enumeration for bond color types."""

        UNIFORM = 'uniform'
        SPLIT = 'split'

    def __init__(self, atom_a: Atom, atom_b: Atom, config: Config = config) -> None:
        """Initialize a bond between two atoms for visualization.

        Parameters
        ----------
        atom_a : Atom
            First atom participating in the bond.
        atom_b : Atom
            Second atom participating in the bond.
        """
        bond_vec = atom_a.center - atom_b.center
        center = (atom_a.center + atom_b.center) / 2

        length = cast(float, np.linalg.norm(bond_vec))
        self.length = length
        self.radius = config.molecule.bond.radius
        self.color_type = self.ColorType(config.molecule.bond.color_type.lower())
        self.mesh: pv.PolyData | list[pv.PolyData] | None
        self.atom_a = atom_a
        self.atom_b = atom_b
        self.plotted = False

        if self.color_type is self.ColorType.UNIFORM:
            self.color = config.molecule.bond.color
        else:
            self.color = [atom_a.atom_type.color, atom_b.atom_type.color]

        if length <= np.finfo(float).eps:
            self.mesh = None
            logger.warning(
                'Cannot render zero-length bond between atoms %s and %s.',
                atom_a.atom_type.name,
                atom_b.atom_type.name,
            )
            return

        if self.color_type is self.ColorType.UNIFORM:
            self.mesh = pv.Cylinder(
                radius=self.radius,
                center=center,
                height=length,
                direction=bond_vec,
            )
        elif self.color_type is self.ColorType.SPLIT:
            atom_radii_adjustement = bond_vec * (atom_b.atom_type.radius - atom_a.atom_type.radius) / length

            center_a = (atom_a.center + center + atom_radii_adjustement / 2) / 2
            center_b = (atom_b.center + center + atom_radii_adjustement / 2) / 2

            atom_radii_adjustement_length = cast(float, np.linalg.norm(atom_radii_adjustement))
            sign = 1 if atom_b.atom_type.radius <= atom_a.atom_type.radius else -1

            mesh_a = pv.Cylinder(
                radius=self.radius,
                center=center_a,
                height=(length + sign * atom_radii_adjustement_length) / 2,
                direction=bond_vec,
            )

            mesh_b = pv.Cylinder(
                radius=self.radius,
                center=center_b,
                height=(length - sign * atom_radii_adjustement_length) / 2,
                direction=bond_vec,
            )

            self.mesh = [mesh_a, mesh_b]

    def _cylinder_between(
        self,
        point_a: NDArray[np.floating],
        point_b: NDArray[np.floating],
    ) -> pv.PolyData | None:
        """Create a cylinder spanning two points.

        Parameters
        ----------
        point_a : NDArray[np.floating]
            Centre of the first cylinder cap.
        point_b : NDArray[np.floating]
            Centre of the second cylinder cap.

        Returns
        -------
        pv.PolyData or None
            Cylinder between the points, or ``None`` when it has no length.
        """
        direction = point_b - point_a
        height = cast(float, np.linalg.norm(direction))
        if height <= np.finfo(float).eps:
            return None

        return pv.Cylinder(
            radius=self.radius,
            center=(point_a + point_b) / 2,
            height=height,
            direction=direction,
        )

    def _trim_distance(self, atom: Atom) -> float:
        """Return the axial distance where the cylinder wall meets an atom.

        The bond surface is a cylinder at ``self.radius`` from its axis.
        Intersecting that cylinder with an atom sphere gives a right triangle
        whose hypotenuse is the atom radius.

        Parameters
        ----------
        atom : Atom
            Atom whose spherical surface bounds the bond.

        Returns
        -------
        float
            Distance from the atom centre along the bond axis.
        """
        squared_distance = atom.atom_type.radius**2 - self.radius**2
        return float(np.sqrt(max(squared_distance, 0.0)))

    def _trim_ends(self) -> None:
        """Shorten the bond analytically so it ends at each atom surface.

        Rebuilding cylinders avoids VTK boolean subtraction, which can abort
        the process for otherwise valid atom and bond geometries.
        """
        if self.mesh is None:
            return

        axis = (self.atom_b.center - self.atom_a.center) / self.length
        trim_a = self._trim_distance(self.atom_a)
        trim_b = self._trim_distance(self.atom_b)
        end_a = self.atom_a.center + axis * trim_a
        end_b = self.atom_b.center - axis * trim_b

        if np.dot(end_b - end_a, axis) <= np.finfo(float).eps:
            self.mesh = None
        elif self.color_type is self.ColorType.SPLIT:
            split = self.atom_a.center + axis * (self.length + trim_a - trim_b) / 2
            mesh_a = self._cylinder_between(end_a, split)
            mesh_b = self._cylinder_between(split, end_b)
            self.mesh = [mesh_a, mesh_b] if mesh_a is not None and mesh_b is not None else None
        else:
            self.mesh = self._cylinder_between(end_a, end_b)

        if self.mesh is None:
            logger.warning(
                'Bond is entirely contained by atoms %s and %s.',
                self.atom_a.atom_type.name,
                self.atom_b.atom_type.name,
            )


class Molecule:
    """Composite object storing rendered atoms and inferred bonds."""

    def __init__(self, atoms: list[ParsedAtom], config: Config = config) -> None:
        """Initialize a molecule from parsed atom data.

        Parameters
        ----------
        atoms : list[ParsedAtom]
            Parsed atoms emitted by :class:`moldenViz.parser.Parser`.
        """
        self.config = config

        # Max radius is used later for plotting
        self.max_radius = 0

        self._get_atoms(atoms)

    def _get_atoms(self, atoms: list[ParsedAtom]) -> None:
        """Convert parsed atoms to visualization atoms and create bonds.

        Parameters
        ----------
        atoms : list[ParsedAtom]
            List of parsed atom objects.
        """
        atomic_numbers = [atom.atomic_number for atom in atoms]
        atom_centers = np.asarray([atom.position for atom in atoms], dtype=float)
        self.atoms = list(map(Atom, atomic_numbers, atom_centers))
        self.max_radius = np.max(np.linalg.norm(atom_centers, axis=1))

        atom_a_indices, atom_b_indices = np.triu_indices(len(atom_centers), k=1)
        pairwise_distances = np.linalg.norm(
            atom_centers[atom_a_indices] - atom_centers[atom_b_indices],
            axis=1,
        )
        within_bond_length = pairwise_distances < self.config.molecule.bond.max_length
        bond_indices = zip(
            atom_a_indices[within_bond_length],
            atom_b_indices[within_bond_length],
            strict=False,
        )

        if self.config.molecule.bond.show:
            for atom_a_ind, atom_b_ind in bond_indices:
                bond = Bond(self.atoms[atom_a_ind], self.atoms[atom_b_ind], self.config)
                self.atoms[atom_a_ind].bonds.append(bond)
                self.atoms[atom_b_ind].bonds.append(bond)

            for atom in self.atoms:
                atom._remove_extra_bonds()  # ruff:ignore[private-member-access]

    def _add_meshes(self, plotter: pv.Plotter, opacity: float = config.molecule.opacity) -> tuple[list[pv.Actor], ...]:
        """Add all molecule meshes (atoms and bonds) to the PyVista plotter.

        Parameters
        ----------
        plotter : pv.Plotter
            The PyVista plotter to add meshes to.
        opacity : float, optional
            The opacity level for the molecule meshes. Default from config.

        Returns
        -------
        tuple[list[pv.Actor], ...]
            A list containing all added actors, a list for the atom actors, and one for the bond actors.
        """
        atom_actors = []
        bond_actors = []
        for atom in self.atoms:
            if self.config.molecule.atom.show:
                atom_actors.append(
                    plotter.add_mesh(
                        atom.mesh,
                        color=atom.atom_type.color,
                        smooth_shading=self.config.smooth_shading,
                        opacity=opacity,
                    ),
                )

            for bond in atom.bonds:
                if bond.plotted or bond.mesh is None:
                    continue

                bond._trim_ends()  # ruff:ignore[private-member-access]
                if bond.mesh is None:
                    continue

                if isinstance(bond.mesh, list):
                    for mesh, color in zip(bond.mesh, bond.color, strict=False):
                        bond_actors.append(plotter.add_mesh(mesh, color=color, opacity=opacity))
                else:
                    if not isinstance(bond.color, str):
                        raise TypeError('Bond color should be a string for uniform color type.')
                    bond_actors.append(plotter.add_mesh(bond.mesh, color=bond.color, opacity=opacity))
                bond.plotted = True

        return atom_actors + bond_actors, atom_actors, bond_actors
