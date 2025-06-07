"""Read and parse a molden file."""

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.special import gamma

logger = logging.getLogger(__name__)


@dataclass
class _Atom:
    label: str
    atomic_number: int
    position: NDArray[np.floating]
    shells: list['_Shell']


@dataclass
class _MolecularOrbital:
    sym: str
    energy: float
    spin: str
    occ: int
    coeffs: NDArray[np.floating]


class _GTO:
    def __init__(self, exp: float, coeff: float) -> None:
        self.exp = exp
        self.coeff = coeff

        self.norm = 0.0

    def normalize(self, l: int) -> None:
        # See (Jiyun Kuang and C D Lin 1997 J. Phys. B: At. Mol. Opt. Phys. 30 2529)
        # page 2532 for the normalization factor
        self.norm = np.sqrt(2 * (2 * self.exp) ** (l + 1.5) / gamma(l + 1.5))


class _Shell:
    def __init__(self, l: int, gtos: list[_GTO]) -> None:
        self.l = l
        self.gtos = gtos

        self.norm = 0.0

    def normalize(self) -> None:
        # See (Jiyun Kuang and C D Lin 1997 J. Phys. B: At. Mol. Opt. Phys. 30 2529)
        # equation 18 and 20 for the normalization factor
        for gto in self.gtos:
            gto.normalize(self.l)

        overlap = 0.0
        for i_gto in self.gtos:
            for j_gto in self.gtos:
                overlap += (
                    i_gto.coeff
                    * j_gto.coeff
                    * (2 * np.sqrt(i_gto.exp * j_gto.exp) / (i_gto.exp + j_gto.exp)) ** (self.l + 1.5)
                )

        self.norm = 1 / np.sqrt(overlap)


class Parser:
    """Parser for molden files.

    Args
    ----
        source: str | list[str]
            The path to the molden file, or the lines from the file.

        only_molecule: bool, optional
            Only parse the atoms and skip molecular orbitals.
            Default is `False`.
    """

    ANGSTROM_TO_BOHR = 1.8897259886

    def __init__(
        self,
        source: str | list[str],
        only_molecule: bool = False,
    ) -> None:
        """Initialize the Parser with either a filename or molden lines."""
        if isinstance(source, str):
            with Path(source).open('r') as file:
                self.molden_lines = file.readlines()
        elif isinstance(source, list):
            self.molden_lines = source

        # Remove leading/trailing whitespace and newline characters
        self.molden_lines = [line.strip() for line in self.molden_lines]

        self.check_molden_format()

        self.atom_ind, self.gto_ind, self.mo_ind = self.divide_molden_lines()

        self.atoms = self.get_atoms()

        if only_molecule:
            return

        self.shells = self.get_shells()
        self.mos = self.get_mos()

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
        angs = 'Angs' in self.molden_lines[self.atom_ind]

        atoms = []
        for line in self.molden_lines[self.atom_ind + 1 : self.gto_ind]:
            label, _, atomic_number, *coords = line.split()

            position = np.array([float(coord) for coord in coords], dtype=float)
            if angs:
                position *= self.ANGSTROM_TO_BOHR

            atoms.append(_Atom(label, int(atomic_number), position, []))

        logger.info('Parsed %s atoms.', len(atoms))
        return atoms

    def get_shells(self) -> list[_Shell]:
        """Parse the Gaussian-type orbitals (GTOs) from the molden file.

        Returns
        -------
            list[_Shell]: A list of `_Shell` objects containing the atom, angular
            momentum quantum number (l), and GTOs for each shell.

        Raises
        ------
            ValueError: If the shell label is not supported or if the GTOs are not
            formatted correctly in the molden file.

        """
        logger.info('Parsing GTO lines...')

        shell_lables = ['s', 'p', 'd', 'f', 'g']

        lines = iter(self.molden_lines[self.gto_ind + 1 : self.mo_ind])

        shells = []
        for atom in self.atoms:
            logger.debug('Parsing GTOs for atom: %s', atom.label)
            _ = next(lines)  # Skip atom index

            # Read shells until a blank line
            while True:
                line = next(lines)
                if not line:
                    break

                shell_label, num_gtos, _ = line.split()
                if shell_label not in shell_lables:
                    raise ValueError(f"Shell label '{shell_label}' is currently not supported.")

                gtos = []
                for _ in range(int(num_gtos)):
                    exp, coeff = next(lines).split()
                    gtos.append(_GTO(float(exp), float(coeff)))

                shell = _Shell(shell_lables.index(shell_label), gtos)
                shell.normalize()

                atom.shells.append(shell)
                shells.append(shell)

        logger.info('Parsed %s GTOs.', len(shells))
        return shells

    def get_mos(self, sort: bool = True) -> list[_MolecularOrbital]:
        """Parse the molecular orbitals (MOs) from the molden file.

        Args
        ----
            mo_inds: int, ArrayLike]
                Indices of the MOs to tabulate. If None, all MOs are tabulated.
            sort: bool
                If true (default), returns the mos sorted by energy. If false, returns the mos in the order
                given in the molden file.


        Returns
        -------
            list[MolecularOrbital]: A list of MolecularOrbital objects containing
            the symmetry, energy, and coefficients for each MO.

        """
        logger.info('Parsing MO coefficients...')

        num_total_gtos = sum(2 * gto.l + 1 for gto in self.shells)

        order = self._gto_order()

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
            for _ in range(num_total_gtos):
                _, coeff = next(lines).split()
                coeffs.append(coeff)

            mo = _MolecularOrbital(
                sym=sym,
                energy=energy,
                spin=spin,
                occ=occ,
                coeffs=np.array(coeffs, dtype=float)[order],
            )

            mos.append(mo)

        logger.info('Parsed MO coefficients.')

        if sort:
            mos = self.sort_mos(mos)

        return mos

    def _gto_order(self) -> list[int]:
        """Return the order of the GTOs in the molden file.

        Molden defines the order of the orbitals as m = 0, 1, -1, 2, -2, ...
        We want it to be m = -l, -l + 1, ..., l - 1, l.

        Note: For l = 1, the order is 1, -1, 0, which is different from the
        general pattern. This is handled separately.

        Returns
        -------
            list[int]: The order of the atomic orbitals.

        """
        order = []
        ind = 0
        for shell in self.shells:
            l = shell.l
            if l == 1:
                order.extend([ind + 1, ind + 2, ind])
            else:
                order.extend([ind + i for i in range(2 * l, -1, -2)])
                order.extend([ind + i for i in range(1, 2 * l, 2)])
            ind += 2 * l + 1

        return order

    @staticmethod
    def sort_mos(mos: list[_MolecularOrbital]) -> list[_MolecularOrbital]:
        """Sort a list of MOs by energy.

        Args
        ----
            mos: list[_MolecularOrbital]
                A list of `_MolecularOrbital` objects to be sorted.

        Returns
        -------
            list[_MolecularOrbital]: A new list containing the `_MolecularOrbital` objects sorted
            by their energy in ascending order.

        """
        logger.info('Sorting MOs by energy...')
        mos = sorted(mos, key=lambda mo: mo.energy)
        logger.info('MOs sorted by energy.')

        return mos
