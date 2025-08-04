import logging
from dataclasses import dataclass
from enum import Enum
from typing import cast

import numpy as np
import pyvista as pv
from numpy.typing import NDArray
from scipy.spatial.distance import pdist, squareform

from ._config_module import Config
from .parser import _Atom

logger = logging.getLogger(__name__)


@dataclass
class AtomType:
    name: str
    color: str
    radius: float
    max_num_bonds: int


config = Config()


ATOM_TYPES = config.atom_types

# Default atom type for invalid atomic numbers
ATOM_X = AtomType('X', 'black', 1, 0)


class Atom:
    def __init__(
        self,
        atomic_number: int,
        center: NDArray[np.floating],
    ) -> None:
        self.atom_type = ATOM_TYPES.get(atomic_number, ATOM_X)
        if self.atom_type is ATOM_X:
            logger.warning(
                "Invalid atomic number: %d. Atom type could not be determined. Using atom 'X' instead.",
                atomic_number,
            )

        self.center = np.array(center)
        self.mesh = pv.Sphere(center=center, radius=self.atom_type.radius)
        self.bonds: list[Bond] = []

    def remove_extra_bonds(self) -> None:
        """Remove the longest bonds if there are more bonds than `max_bonds`."""
        if len(self.bonds) <= self.atom_type.max_num_bonds:
            return

        self.bonds.sort(key=lambda x: x.length)

        for bond in self.bonds[self.atom_type.max_num_bonds :]:
            bond.mesh = None


class Bond:
    class ColorType(Enum):
        UNIFORM = 'uniform'
        SPLIT = 'split'

    def __init__(self, atom_a: Atom, atom_b: Atom) -> None:
        bond_vec = atom_a.center - atom_b.center
        center = (atom_a.center + atom_b.center) / 2

        length = cast(float, np.linalg.norm(bond_vec))
        self.length = length

        if config.molecule.bond.color_type.lower() == self.ColorType.UNIFORM.value:
            self.mesh = pv.Cylinder(
                radius=config.molecule.bond.radius,
                center=center,
                height=length,
                direction=bond_vec,
            )
            self.color = config.molecule.bond.color
        elif config.molecule.bond.color_type.lower() == self.ColorType.SPLIT.value:
            center = atom_a.center + atom_b.center
            center_a = center / 4
            center_b = center * 3 / 4

            mesh_a = pv.Cylinder(
                radius=config.molecule.bond.radius,
                center=center_a,
                height=length / 2,
                direction=bond_vec,
            )

            mesh_b = pv.Cylinder(
                radius=config.molecule.bond.radius,
                center=center_b,
                height=length / 2,
                direction=bond_vec,
            )

            self.mesh = [mesh_a, mesh_b]
            self.color = [atom_a.atom_type.color, atom_b.atom_type.color]
        else:
            raise ValueError(
                f'Invalid bond color type: {config.molecule.bond.color_type}. '
                f'Expected one of {[color_type.value for color_type in self.ColorType]}.',
            )

        self.atom_a = atom_a
        self.atom_b = atom_b

        self.plotted = False

    def trim_ends(self) -> None:
        """Remove the ends of the bond that are going into the atoms."""
        if self.mesh is None:
            return

        warning = False
        if isinstance(self.mesh, list):
            self.mesh = [mesh.triangulate() - atom.mesh for mesh, atom in zip(self.mesh, [self.atom_a, self.atom_b])]
            if any(mesh.n_points == 0 for mesh in self.mesh):
                warning = True
        else:
            self.mesh = self.mesh.triangulate() - self.atom_a.mesh - self.atom_b.mesh
            if self.mesh.n_points == 0:
                warning = True

        if warning:
            logger.warning(
                'Error: Bond mesh is empty between atoms %s and %s.',
                self.atom_a.atom_type.name,
                self.atom_b.atom_type.name,
            )
            self.mesh = None


class Molecule:
    def __init__(self, atoms: list[_Atom]) -> None:
        # Max radius is used later for plotting
        self.max_radius = 0

        self.get_atoms(atoms)

    def get_atoms(self, atoms: list[_Atom]) -> None:
        atomic_numbers = [atom.atomic_number for atom in atoms]
        atom_centers = [atom.position for atom in atoms]
        self.atoms = list(map(Atom, atomic_numbers, atom_centers))
        self.max_radius = np.max(np.linalg.norm(atom_centers, axis=1))

        distances = squareform(pdist(atom_centers))  # Compute pairwise distances
        mask = np.triu(np.ones_like(distances, dtype=bool), k=1)  # Ensure boolean mask
        indices = np.where((distances < config.molecule.max_bond_lenght) & mask)  # Apply mask

        for atom_a_ind, atom_b_ind in zip(indices[0], indices[1]):
            bond = Bond(self.atoms[atom_a_ind], self.atoms[atom_b_ind])
            self.atoms[atom_a_ind].bonds.append(bond)
            self.atoms[atom_b_ind].bonds.append(bond)

        for atom in self.atoms:
            atom.remove_extra_bonds()

    def add_meshes(self, plotter: pv.Plotter, opacity: float = config.molecule.opacity) -> list[pv.Actor]:
        actors = []
        for atom in self.atoms:
            actors.append(plotter.add_mesh(atom.mesh, color=atom.atom_type.color, smooth_shading=config.smooth_shading))
            for bond in atom.bonds:
                if bond.plotted or bond.mesh is None:
                    continue

                bond.trim_ends()
                if isinstance(bond.mesh, list):
                    for mesh, color in zip(bond.mesh, bond.color):
                        actors.append(plotter.add_mesh(mesh, color=color, opacity=opacity))
                else:
                    actors.append(plotter.add_mesh(bond.mesh, color=bond.color, opacity=opacity))
                bond.plotted = True

        return actors
