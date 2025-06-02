"""Read and parse a molden file."""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from scipy.special import gamma

logger = logging.getLogger(__name__)


@dataclass
class _Atom:
    label: str
    atomic_number: int
    position: NDArray[np.floating]
    gtos: list['_Gto']


@dataclass
class _MolecularOrbital:
    symmetry: str
    energy: float
    spin: str
    occ: int
    coeffs: NDArray[np.floating]


class _GaussianPrimitive:
    def __init__(self, exp: float, coeff: float) -> None:
        self.exp = exp
        self.coeff = coeff

        self.norm = 0.0

    def normalize(self, l: int) -> None:
        # See (Jiyun Kuang and C D Lin 1997 J. Phys. B: At. Mol. Opt. Phys. 30 2529)
        # page 2532 for the normalization factor
        self.norm = np.sqrt(2 * (2 * self.exp) ** (l + 1.5) / gamma(l + 1.5))


class _Gto:
    def __init__(self, l: int, prims: list[_GaussianPrimitive]) -> None:
        self.l = l
        self.prims = prims

        self.norm = 0.0

    def normalize(self) -> None:
        # See (Jiyun Kuang and C D Lin 1997 J. Phys. B: At. Mol. Opt. Phys. 30 2529)
        # equation 18 and 20 for the normalization factor
        for prim in self.prims:
            prim.normalize(self.l)

        overlap = 0.0
        for i_prim in self.prims:
            for j_prim in self.prims:
                overlap += (
                    i_prim.coeff
                    * j_prim.coeff
                    * (2 * np.sqrt(i_prim.exp * j_prim.exp) / (i_prim.exp + j_prim.exp)) ** (self.l + 1.5)
                )

        self.norm = 1 / np.sqrt(overlap)


class Parser:
    """Parser for molden files."""

    def __init__(
        self,
        filename: Optional[str] = None,
        molden_lines: Optional[list[str]] = None,
    ) -> None:
        """Initialize the Parser with either a filename or molden lines.

        Args:
        ----
            filename (Optional[str]): The path to the molden file.
            molden_lines (Optional[list[str]]): A list of lines from a molden file.

        Raises:
        ------
            ValueError: If both 'filename' and 'molden_lines' are provided,
                        or if neither is provided.

        """
        if filename and molden_lines is not None:
            raise ValueError("Provide either 'filename' or 'molden_lines', not both.")

        if filename:
            with Path(filename).open('r') as file:
                self.molden_lines = file.readlines()
        elif molden_lines is not None:
            if not molden_lines:
                raise ValueError("'molden_lines' was provided but is empty.")

            self.molden_lines = molden_lines
        else:
            raise ValueError("Must provide either 'filename' or 'molden_lines'.")

        self.molden_lines = [
            line.strip() for line in self.molden_lines
        ]  # Remove leading/trailing whitespace and newlines

        self.check_molden_format()

        self.atom_ind, self.gto_ind, self.mo_ind = self.divide_molden_lines()

        self.atoms = self.get_atoms()
        self.gtos = self.get_gtos()
        self.mos = self.get_mos()

        self._sort_mos()

    def check_molden_format(self) -> None:
        """Check if the provided molden lines conform to the expected format.

        Raises
        ------
            ValueError: If the molden lines do not contain the required sections
                        or if they are in an unsupported format.

        """
        logger.info('Checking molden format...')
        if not self.molden_lines:
            raise ValueError('The provided molden lines are empty.')

        if not any('[Atoms]' in line for line in self.molden_lines):
            raise ValueError("No '[Atoms]' section found in the molden file.")

        if not any('[GTO]' in line for line in self.molden_lines):
            raise ValueError("No '[GTO]' section found in the molden file.")

        if not any('[MO]' in line for line in self.molden_lines):
            raise ValueError("No '[MO]' section found in the molden file.")

        if not any(orbs in line for orbs in ['5D', '9G'] for line in self.molden_lines):
            raise ValueError('Cartesian orbitals functions are not currently supported.')

        logger.info('Molden format check passed.')

    def divide_molden_lines(self) -> tuple[int, int, int]:
        """Divide the molden lines into sections for atoms, GTOs, and MOs.

        Returns
        -------
            tuple[int, int, int]: Indices of the '[Atoms]', '[GTO]', and '[MO]' lines.

        Raises
        ------
            ValueError: If the molden lines do not contain the required sections.

        """
        logger.info('Dividing molden lines into sections...')
        if '[Atoms] AU' in self.molden_lines:
            atom_ind = self.molden_lines.index('[Atoms] AU')
        elif '[Atoms] Angs' in self.molden_lines:
            atom_ind = self.molden_lines.index('[Atoms] Angs')
        else:
            raise ValueError("No '[Atoms] (AU/Angs)' section found in the molden file.")

        gto_ind = self.molden_lines.index('[GTO]')

        mo_ind = self.molden_lines.index('[MO]')

        logger.info('Finished dividing molden lines.')
        return atom_ind, gto_ind, mo_ind

    def get_atoms(self) -> list[_Atom]:
        """Parse the atoms from the molden file.

        Returns
        -------
            list[Atom]: A list of Atom objects containing the label, atomic number,
            and position for each atom.

        """
        logger.info('Parsing atoms...')
        atoms = []
        for line in self.molden_lines[self.atom_ind + 1 : self.gto_ind]:
            label, _, atomic_number, *coords = line.split()

            position = np.array([float(coord) for coord in coords], dtype=float)

            atoms.append(_Atom(label, int(atomic_number), position, []))

        logger.info('Parsed %s atoms.', len(atoms))
        return atoms

    def get_gtos(self) -> list[_Gto]:
        """Parse the Gaussian-type orbitals (GTOs) from the molden file.

        Returns
        -------
            list[_Gto]: A list of Gto objects containing the atom, angular
            momentum quantum number (l), and Gaussian primitives for each gto.

        Raises
        ------
            ValueError: If the shell label is not supported or if the GTOs are not
            formatted correctly in the molden file.

        """
        logger.info('Parsing GTO lines...')

        shell_lables = ['s', 'p', 'd', 'f', 'g']

        lines = iter(self.molden_lines[self.gto_ind + 1 : self.mo_ind])

        gtos = []
        for atom in self.atoms:
            logger.debug('Parsing GTOs for atom: %s', atom.label)
            _ = next(lines)  # Skip atom index

            # Read shells until a blank line
            while True:
                line = next(lines)
                if not line:
                    break

                shell_label, number_of_primitives, _ = line.split()
                if shell_label not in shell_lables:
                    raise ValueError(f"Shell label '{shell_label}' is currently not supported.")

                prims = []
                for _ in range(int(number_of_primitives)):
                    exp, coeff = next(lines).split()
                    prims.append(_GaussianPrimitive(float(exp), float(coeff)))

                gto = _Gto(shell_lables.index(shell_label), prims)
                gto.normalize()

                atom.gtos.append(gto)
                gtos.append(gto)

        logger.info('Parsed %s GTOs.', len(gtos))
        return gtos

    def get_mos(self) -> list[_MolecularOrbital]:
        """Parse the molecular orbitals (MOs) from the molden file.

        Returns
        -------
            list[MolecularOrbital]: A list of MolecularOrbital objects containing
            the symmetry, energy, and coefficients for each MO.

        """
        logger.info('Parsing MO coefficients...')

        num_atomic_orbs = sum(2 * gto.l + 1 for gto in self.gtos)

        order = self._atomic_orbs_order()

        lines = self.molden_lines[self.mo_ind + 1 :]
        total_num_mos = sum('Sym=' in line for line in lines)
        lines = iter(lines)

        mos = []
        for _ in range(total_num_mos):
            _, sym = next(lines).split()

            energy_line = next(lines)
            energy = float(energy_line.split()[1])

            _, spin = next(lines).split()

            occ_line = next(lines)
            occ = int(float(occ_line.split()[1]))

            coeffs = []
            for _ in range(num_atomic_orbs):
                _, coeff = next(lines).split()
                coeffs.append(coeff)

            mo = _MolecularOrbital(
                symmetry=sym,
                energy=energy,
                spin=spin,
                occ=occ,
                coeffs=np.array(coeffs, dtype=float)[order],
            )

            mos.append(mo)

        logger.info('Parsed MO coefficients.')
        return mos

    def _atomic_orbs_order(self) -> list[int]:
        """Return the order of the atomic orbitals in the molden file.

        Molden defines the order of the orbitals as 0, 1, -1, 2, -2, ...
        We want it to be -l, -l + 1, ..., l - 1, l.

        Note: For l = 1, the order is 1, -1, 0, which is different from the
        general pattern. This is handled separately.

        Returns
        -------
            list[int]: The order of the atomic orbitals.

        """
        order = []
        ind = 0
        for gto in self.gtos:
            l = gto.l
            if l == 1:
                order.extend([ind + 1, ind + 2, ind])
            else:
                order.extend([ind + i for i in range(2 * l, -1, -2)])
                order.extend([ind + i for i in range(1, 2 * l, 2)])
            ind += 2 * l + 1

        return order

    def _sort_mos(self) -> None:
        """Sort the MOs by energy."""
        logger.info('Sorting MOs by energy...')
        self.mos.sort(key=lambda mo: mo.energy)
        logger.info('MOs sorted by energy.')
