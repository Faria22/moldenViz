import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import numpy as np
import pyvista as pv
from _config_module import load_atom_types
from numpy.typing import NDArray
from scipy.spatial.distance import pdist, squareform

from .parser import _Atom

logger = logging.getLogger(__name__)


@dataclass
class AtomType:
    name: str
    color: str
    radius: float
    max_bonds: int


# Atom type based on charge
# Colors are based on CPK color scheme
# https://sciencenotes.org/molecule-atom-colors-cpk-colors/
with Path('atoms.json').open('r') as f:
    atom_types = json.load(f)

ATOM_TYPES = {int(k): AtomType(**v) for k, v in atom_types.items()}
ATOM_TYPES = load_atom_types()

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
        if len(self.bonds) <= self.atom_type.max_bonds:
            return

        self.bonds.sort(key=lambda x: x.length)

        for bond in self.bonds[self.atom_type.max_bonds :]:
            bond.mesh = None


class Bond:
    def __init__(self, atom_a: Atom, atom_b: Atom, radius: float = 0.15) -> None:
        center = (atom_a.center + atom_b.center) / 2

        bond_vec = atom_a.center - atom_b.center
        length = cast(float, np.linalg.norm(bond_vec))

        self.length = length

        self.mesh = pv.Cylinder(radius=radius, center=center, height=length, direction=bond_vec)
        self.atom_a = atom_a
        self.atom_b = atom_b

        self.color = 'grey'
        self.plotted = False

    def trim_ends(self) -> None:
        """Remove the ends of the bond that are going into the atoms."""
        if self.mesh is None:
            return

        self.mesh = self.mesh.triangulate() - self.atom_a.mesh - self.atom_b.mesh

        if self.mesh.n_points == 0:
            logger.warning(
                'Error: Bond mesh is empty between atoms %s and %s.',
                self.atom_a.atom_type.name,
                self.atom_b.atom_type.name,
            )
            self.mesh = None


class Molecule:
    def __init__(self, atoms: list[_Atom], max_bond_length: float = 4) -> None:
        # Max_bond_length helps the program skip over any bonds that should not exist
        self.max_bond_length = max_bond_length

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
        indices = np.where((distances < self.max_bond_length) & mask)  # Apply mask

        for atom_a_ind, atom_b_ind in zip(indices[0], indices[1]):
            bond = Bond(self.atoms[atom_a_ind], self.atoms[atom_b_ind])
            self.atoms[atom_a_ind].bonds.append(bond)
            self.atoms[atom_b_ind].bonds.append(bond)

        for atom in self.atoms:
            atom.remove_extra_bonds()

    def add_meshes(self, plotter: pv.Plotter, opacity: float = 1) -> list[pv.Actor]:
        actors = []
        for atom in self.atoms:
            actors.append(plotter.add_mesh(atom.mesh, color=atom.atom_type.color, smooth_shading=True))
            for bond in atom.bonds:
                if bond.plotted or bond.mesh is None:
                    continue

                bond.trim_ends()
                actors.append(plotter.add_mesh(bond.mesh, color=bond.color, opacity=opacity))
                bond.plotted = True

        return actors
