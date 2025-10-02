"""Utility objects used by the Plotter to build molecular meshes."""

import logging
from enum import Enum
from typing import cast

import numpy as np
import pyvista as pv
from numpy.typing import NDArray
from scipy.spatial.distance import pdist, squareform

from ._config_module import AtomType, Config
from .parser import _Atom

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

    def remove_extra_bonds(self) -> None:
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

    def __init__(self, atom_a: Atom, atom_b: Atom) -> None:
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
        self.mesh: pv.PolyData | list[pv.PolyData] | None

        if config.molecule.bond.color_type.lower() == self.ColorType.UNIFORM.value:
            self.mesh = pv.Cylinder(
                radius=config.molecule.bond.radius,
                center=center,
                height=length,
                direction=bond_vec,
            )
            self.color = config.molecule.bond.color
        elif config.molecule.bond.color_type.lower() == self.ColorType.SPLIT.value:
            atom_radii_adjustement = bond_vec * (atom_b.atom_type.radius - atom_a.atom_type.radius) / length

            center_a = (atom_a.center + center + atom_radii_adjustement / 2) / 2
            center_b = (atom_b.center + center + atom_radii_adjustement / 2) / 2

            atom_radii_adjustement_length = cast(float, np.linalg.norm(atom_radii_adjustement))
            sign = 1 if atom_b.atom_type.radius <= atom_a.atom_type.radius else -1

            mesh_a = pv.Cylinder(
                radius=config.molecule.bond.radius,
                center=center_a,
                height=(length + sign * atom_radii_adjustement_length) / 2,
                direction=bond_vec,
            )

            mesh_b = pv.Cylinder(
                radius=config.molecule.bond.radius,
                center=center_b,
                height=(length - sign * atom_radii_adjustement_length) / 2,
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

    @staticmethod
    def _trim_atom_from_bond(bond_mesh: pv.PolyData, atom_mesh: pv.PolyData) -> pv.PolyData:
        """Trim the bond mesh to remove parts that intrude into the atom mesh.

        Parameters
        ----------
        bond_mesh : pv.PolyData
            The mesh representing the bond.
        atom_mesh : pv.PolyData
            The mesh representing the atom.

        Returns
        -------
        pv.PolyData
            The trimmed bond mesh.
        """
        mesh = cast(pv.PolyData, bond_mesh.triangulate()) - atom_mesh
        return cast(pv.PolyData, mesh)

    def trim_ends(self) -> None:
        """Trim bond geometry so it does not intrude into atom spheres."""
        if self.mesh is None:
            return

        warning = False
        if isinstance(self.mesh, list):
            self.mesh = [
                self._trim_atom_from_bond(mesh, atom.mesh) for mesh, atom in zip(self.mesh, (self.atom_a, self.atom_b))
            ]
            if any(mesh.n_points == 0 for mesh in self.mesh):
                warning = True
        else:
            self.mesh = self._trim_atom_from_bond(self.mesh, self.atom_a.mesh)
            self.mesh = self._trim_atom_from_bond(self.mesh, self.atom_b.mesh)
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
    """Composite object storing rendered atoms and inferred bonds."""

    def __init__(self, atoms: list[_Atom]) -> None:
        """Initialize a molecule from parsed atom data.

        Parameters
        ----------
        atoms : list[_Atom]
            Parsed atoms emitted by :class:`moldenViz.parser.Parser`.
        """
        # Max radius is used later for plotting
        self.max_radius = 0

        self.get_atoms(atoms)

    def get_atoms(self, atoms: list[_Atom]) -> None:
        """Convert parsed atoms to visualization atoms and create bonds.

        Parameters
        ----------
        atoms : list[_Atom]
            List of parsed atom objects.
        """
        atomic_numbers = [atom.atomic_number for atom in atoms]
        atom_centers = [atom.position for atom in atoms]
        self.atoms = list(map(Atom, atomic_numbers, atom_centers))
        self.max_radius = np.max(np.linalg.norm(atom_centers, axis=1))

        distances = squareform(pdist(atom_centers))  # Compute pairwise distances
        mask = np.triu(np.ones_like(distances, dtype=bool), k=1)  # Ensure boolean mask
        indices = np.where((distances < config.molecule.bond.max_length) & mask)  # Apply mask

        for atom_a_ind, atom_b_ind in zip(indices[0], indices[1]):
            bond = Bond(self.atoms[atom_a_ind], self.atoms[atom_b_ind])
            self.atoms[atom_a_ind].bonds.append(bond)
            self.atoms[atom_b_ind].bonds.append(bond)

        for atom in self.atoms:
            atom.remove_extra_bonds()

    def add_meshes(self, plotter: pv.Plotter, opacity: float = config.molecule.opacity) -> list[pv.Actor]:
        """Add all molecule meshes (atoms and bonds) to the PyVista plotter.

        Parameters
        ----------
        plotter : pv.Plotter
            The PyVista plotter to add meshes to.
        opacity : float, optional
            The opacity level for the molecule meshes. Default from config.

        Returns
        -------
        list[pv.Actor]
            List of PyVista actors that were added to the plotter.
        """
        actors = []
        for atom in self.atoms:
            if config.molecule.atom.show:
                actors.append(
                    plotter.add_mesh(
                        atom.mesh,
                        color=atom.atom_type.color,
                        smooth_shading=config.smooth_shading,
                        opacity=opacity,
                    ),
                )
            for bond in atom.bonds:
                if bond.plotted or bond.mesh is None or not config.molecule.bond.show:
                    continue

                bond.trim_ends()
                if bond.mesh is None:
                    continue

                if isinstance(bond.mesh, list):
                    for mesh, color in zip(bond.mesh, bond.color):
                        actors.append(plotter.add_mesh(mesh, color=color, opacity=opacity))
                else:
                    if not isinstance(bond.color, str):
                        raise TypeError('Bond color should be a string for uniform color type.')
                    actors.append(plotter.add_mesh(bond.mesh, color=bond.color, opacity=opacity))
                bond.plotted = True

        return actors
