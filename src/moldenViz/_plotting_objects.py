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
# Colors are based on CPK color scheme
# https://sciencenotes.org/molecule-atom-colors-cpk-colors/
ATOM_TYPES = {
    1: AtomType('H', 'FFFFFF', 0.2, 1),
    2: AtomType('He', 'D9FFFF', 0.286, 0),
    3: AtomType('Li', 'CC80FF', 0.34, 1),
    4: AtomType('Be', 'C2FF00', 0.589, 2),
    5: AtomType('B', 'FFB5B5', 0.415, 3),
    6: AtomType('C', '909090', 0.4, 4),
    7: AtomType('N', '3050F8', 0.4, 3),
    8: AtomType('O', 'FF0D0D', 0.4, 2),
    9: AtomType('F', '90E050', 0.32, 1),
    10: AtomType('Ne', 'B3E3F5', 0.423, 0),
    11: AtomType('Na', 'AB5CF2', 0.485, 1),
    12: AtomType('Mg', '8AFF00', 0.55, 2),
    13: AtomType('Al', 'BFA6A6', 0.675, 3),
    14: AtomType('Si', 'F0C8A0', 0.6, 4),
    15: AtomType('P', 'FF8000', 0.525, 5),
    16: AtomType('S', 'FFFF30', 0.51, 2),
    17: AtomType('Cl', '1FF01F', 0.495, 1),
    18: AtomType('Ar', '80D1E3', 0.508, 0),
    19: AtomType('K', '8F40D4', 0.665, 1),
    20: AtomType('Ca', '3DFF00', 0.495, 2),
    21: AtomType('Sc', 'E6E6E6', 0.72, 3),
    22: AtomType('Ti', 'BFC2C7', 0.735, 6),
    23: AtomType('V', 'A6A6AB', 0.665, 6),
    24: AtomType('Cr', '8A99C7', 0.675, 6),
    25: AtomType('Mn', '9C7AC7', 0.675, 6),
    26: AtomType('Fe', 'E06633', 0.67, 6),
    27: AtomType('Co', 'F090A0', 0.615, 6),
    28: AtomType('Ni', '50D050', 0.75, 6),
    29: AtomType('Cu', 'C88033', 0.76, 6),
    30: AtomType('Zn', '7D80B0', 0.725, 6),
    31: AtomType('Ga', 'C28F8F', 0.61, 3),
    32: AtomType('Ge', '668F8F', 0.585, 4),
    33: AtomType('As', 'BD80E3', 0.605, 3),
    34: AtomType('Se', 'FFA100', 0.61, 2),
    35: AtomType('Br', 'A62929', 0.605, 1),
    36: AtomType('Kr', '5CB8D1', 0.524, 0),
    37: AtomType('Rb', '702EB0', 0.735, 1),
    38: AtomType('Sr', '00FF00', 0.56, 2),
    39: AtomType('Y', '94FFFF', 0.89, 3),
    40: AtomType('Zr', '94E0E0', 0.78, 6),
    41: AtomType('Nb', '73C2C9', 0.74, 6),
    42: AtomType('Mo', '54B5B5', 0.735, 6),
    43: AtomType('Tc', '3B9E9E', 0.675, 6),
    44: AtomType('Ru', '248F8F', 0.7, 6),
    45: AtomType('Rh', '0A7D8C', 0.725, 6),
    46: AtomType('Pd', '006985', 0.75, 6),
    47: AtomType('Ag', 'C0C0C0', 0.795, 6),
    48: AtomType('Cd', 'FFD98F', 0.845, 6),
    49: AtomType('In', 'A67573', 0.815, 3),
    50: AtomType('Sn', '668080', 0.73, 4),
    51: AtomType('Sb', '9E63B5', 0.73, 3),
    52: AtomType('Te', 'D47A00', 0.735, 2),
    53: AtomType('I', '940094', 0.7, 1),
    54: AtomType('Xe', '429EB0', 0.577, 0),
    55: AtomType('Cs', '57178F', 0.835, 1),
    56: AtomType('Ba', '00C900', 0.67, 2),
    57: AtomType('La', '70D4FF', 0.935, 3),
    58: AtomType('Ce', 'FFFFC7', 0.915, 6),
    59: AtomType('Pr', 'D9FFC7', 0.91, 6),
    60: AtomType('Nd', 'C7FFC7', 0.905, 6),
    61: AtomType('Pm', 'A3FFC7', 0.9, 6),
    62: AtomType('Sm', '8FFFC7', 0.9, 6),
    63: AtomType('Eu', '61FFC7', 0.995, 6),
    64: AtomType('Gd', '45FFC7', 0.895, 6),
    65: AtomType('Tb', '30FFC7', 0.88, 6),
    66: AtomType('Dy', '1FFFC7', 0.875, 6),
    67: AtomType('Ho', '00FF9C', 0.87, 6),
    68: AtomType('Er', '00E675', 0.865, 6),
    69: AtomType('Tm', '00D452', 0.86, 6),
    70: AtomType('Yb', '00BF38', 0.97, 6),
    71: AtomType('Lu', '00AB24', 0.86, 6),
    72: AtomType('Hf', '4DC2FF', 0.785, 6),
    73: AtomType('Ta', '4DA6FF', 0.715, 6),
    74: AtomType('W', '2194D6', 0.685, 6),
    75: AtomType('Re', '267DAB', 0.675, 6),
    76: AtomType('Os', '266696', 0.685, 6),
    77: AtomType('Ir', '175487', 0.66, 6),
    78: AtomType('Pt', 'D0D0E0', 0.75, 6),
    79: AtomType('Au', 'FFD123', 0.75, 6),
    80: AtomType('Hg', 'B8B8D0', 0.85, 6),
    81: AtomType('Tl', 'A6544D', 0.775, 3),
    82: AtomType('Pb', '575961', 0.77, 4),
    83: AtomType('Bi', '9E4FB5', 0.77, 3),
    84: AtomType('Po', 'AB5C00', 0.84, 2),
    85: AtomType('At', '754F45', 1, 1),
    86: AtomType('Rn', '428296', 1, 0),
    87: AtomType('Fr', '420066', 1, 1),
    88: AtomType('Ra', '007D00', 0.95, 2),
    89: AtomType('Ac', '70ABFA', 0.94, 3),
    90: AtomType('Th', '00BAFF', 0.895, 6),
    91: AtomType('Pa', '00A1FF', 0.805, 6),
    92: AtomType('U', '008FFF', 0.79, 6),
    93: AtomType('Np', '0080FF', 0.775, 6),
    94: AtomType('Pu', '006BFF', 1, 6),
    95: AtomType('Am', '545CF2', 0.755, 6),
    96: AtomType('Cm', '785CE3', 1, 6),
    97: AtomType('Bk', '8A4FE3', 1, 6),
    98: AtomType('Cf', 'A136D4', 0.765, 6),
    99: AtomType('Es', 'B31FD4', 0.76, 6),
    100: AtomType('Fm', 'B31FBA', 0.755, 6),
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
            actors.append(plotter.add_mesh(atom.mesh, color=atom.atom_type.color, smooth_shading=True))
            for bond in atom.bonds:
                if bond.plotted or bond.mesh is None:
                    continue

                bond.trim_ends()
                actors.append(plotter.add_mesh(bond.mesh, color=bond.color, opacity=opacity))
                bond.plotted = True

        return actors
