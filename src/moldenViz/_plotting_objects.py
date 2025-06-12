import logging
from dataclasses import dataclass
from typing import cast

import numpy as np
import pyvista as pv
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
ATOM_TYPES = {
    1: AtomType('H', 'white', 0.2, 1),
    2: AtomType('He', 'cyan', 0.286, 0),
    3: AtomType('Li', 'lightyellow', 0.34, 1),
    4: AtomType('Be', 'lightgreen', 0.589, 2),
    5: AtomType('B', 'brown', 0.415, 3),
    6: AtomType('C', 'dimgrey', 0.4, 4),
    7: AtomType('N', 'blue', 0.4, 3),
    8: AtomType('O', 'red', 0.4, 2),
    9: AtomType('F', 'purple', 0.32, 1),
    10: AtomType('Ne', 'cyan', 0.423, 0),
    11: AtomType('Na', 'lightblue', 0.485, 1),
    12: AtomType('Mg', 'darkgreen', 0.55, 2),
    13: AtomType('Al', 'lightgrey', 0.675, 3),
    14: AtomType('Si', 'darkgrey', 0.6, 4),
    15: AtomType('P', 'orange', 0.525, 3),
    16: AtomType('S', 'yellow', 0.51, 2),
    17: AtomType('Cl', 'green', 0.495, 1),
    18: AtomType('Ar', 'cyan', 0.508, 0),
    19: AtomType('K', 'navy', 0.665, 1),
    20: AtomType('Ca', 'snow', 0.495, 2),
    21: AtomType('Sc', 'lightblue', 0.72, 3),
    22: AtomType('Ti', 'lightblue', 0.735, 6),
    23: AtomType('V', 'lightblue', 0.665, 6),
    24: AtomType('Cr', 'lightblue', 0.675, 6),
    25: AtomType('Mn', 'lightblue', 0.675, 6),
    26: AtomType('Fe', 'darkorange', 0.67, 6),
    27: AtomType('Co', 'darkorange', 0.615, 6),
    28: AtomType('Ni', 'darkorange', 0.75, 6),
    29: AtomType('Cu', 'darkorange', 0.76, 6),
    30: AtomType('Zn', 'darkorange', 0.725, 6),
    31: AtomType('Ga', 'darkorange', 0.61, 3),
    32: AtomType('Ge', 'orange', 0.585, 4),
    33: AtomType('As', 'orange', 0.605, 3),
    34: AtomType('Se', 'orange', 0.61, 2),
    35: AtomType('Br', 'orange', 0.605, 1),
    36: AtomType('Kr', 'cyan', 0.524, 0),
    37: AtomType('Rb', 'orange', 0.735, 1),
    38: AtomType('Sr', 'orange', 0.56, 2),
    39: AtomType('Y', 'orange', 0.89, 3),
    40: AtomType('Zr', 'darkorange', 0.78, 6),
    41: AtomType('Nb', 'darkorange', 0.74, 6),
    42: AtomType('Mo', 'darkorange', 0.735, 6),
    43: AtomType('Tc', 'darkorange', 0.675, 6),
    44: AtomType('Ru', 'darkorange', 0.7, 6),
    45: AtomType('Rh', 'darkorange', 0.725, 6),
    46: AtomType('Pd', 'darkorange', 0.75, 6),
    47: AtomType('Ag', 'darkorange', 0.795, 6),
    48: AtomType('Cd', 'darkorange', 0.845, 6),
    49: AtomType('In', 'orange', 0.815, 3),
    50: AtomType('Sn', 'orange', 0.73, 4),
    51: AtomType('Sb', 'orange', 0.73, 3),
    52: AtomType('Te', 'orange', 0.735, 2),
    53: AtomType('I', 'orange', 0.7, 1),
    54: AtomType('Xe', 'cyan', 0.577, 0),
    55: AtomType('Cs', 'orange', 0.835, 1),
    56: AtomType('Ba', 'orange', 0.67, 2),
    57: AtomType('La', 'orange', 0.935, 3),
    58: AtomType('Ce', 'darkorange', 0.915, 6),
    59: AtomType('Pr', 'darkorange', 0.91, 6),
    60: AtomType('Nd', 'darkorange', 0.905, 6),
    61: AtomType('Pm', 'darkorange', 0.9, 6),
    62: AtomType('Sm', 'darkorange', 0.9, 6),
    63: AtomType('Eu', 'darkorange', 0.995, 6),
    64: AtomType('Gd', 'darkorange', 0.895, 6),
    65: AtomType('Tb', 'darkorange', 0.88, 6),
    66: AtomType('Dy', 'darkorange', 0.875, 6),
    67: AtomType('Ho', 'darkorange', 0.87, 6),
    68: AtomType('Er', 'darkorange', 0.865, 6),
    69: AtomType('Tm', 'darkorange', 0.86, 6),
    70: AtomType('Yb', 'darkorange', 0.97, 6),
    71: AtomType('Lu', 'darkorange', 0.86, 6),
    72: AtomType('Hf', 'darkorange', 0.785, 6),
    73: AtomType('Ta', 'darkorange', 0.715, 6),
    74: AtomType('W', 'darkorange', 0.685, 6),
    75: AtomType('Re', 'darkorange', 0.675, 6),
    76: AtomType('Os', 'darkorange', 0.685, 6),
    77: AtomType('Ir', 'darkorange', 0.66, 6),
    78: AtomType('Pt', 'darkorange', 0.75, 6),
    79: AtomType('Au', 'darkorange', 0.75, 6),
    80: AtomType('Hg', 'darkorange', 0.85, 6),
    81: AtomType('Tl', 'orange', 0.775, 3),
    82: AtomType('Pb', 'orange', 0.77, 4),
    83: AtomType('Bi', 'orange', 0.77, 3),
    84: AtomType('Po', 'orange', 0.84, 2),
    85: AtomType('At', 'orange', 1, 1),
    86: AtomType('Rn', 'cyan', 1, 0),
    87: AtomType('Fr', 'darkorange', 1, 1),
    88: AtomType('Ra', 'darkorange', 0.95, 2),
    89: AtomType('Ac', 'darkorange', 0.94, 3),
    90: AtomType('Th', 'darkorange', 0.895, 6),
    91: AtomType('Pa', 'darkorange', 0.805, 6),
    92: AtomType('U', 'darkorange', 0.79, 6),
    93: AtomType('Np', 'darkorange', 0.775, 6),
    94: AtomType('Pu', 'darkorange', 1, 6),
    95: AtomType('Am', 'darkorange', 0.755, 6),
    96: AtomType('Cm', 'darkorange', 1, 6),
    97: AtomType('Bk', 'darkorange', 1, 6),
    98: AtomType('Cf', 'darkorange', 0.765, 6),
    99: AtomType('Es', 'darkorange', 0.76, 6),
    100: AtomType('Fm', 'darkorange', 0.755, 6),
}

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
            actors.append(plotter.add_mesh(atom.mesh, color=atom.atom_type.color))
            for bond in atom.bonds:
                if bond.plotted or bond.mesh is None:
                    continue

                bond.trim_ends()
                actors.append(plotter.add_mesh(bond.mesh, color=bond.color, opacity=opacity))
                bond.plotted = True

        return actors
