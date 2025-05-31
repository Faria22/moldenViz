import logging
from dataclasses import dataclass

from pathlib import Path
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from scipy.special import gamma

logger = logging.getLogger(__name__)


@dataclass
class Atom:
    label: str
    atomic_number: int
    position: NDArray[np.float64]


class MolecularOrbital:
    def __init__(self, symmetry: str, energy: float, coeffs: NDArray[np.float64]) -> None:
        self.symmetry = symmetry
        self.energy = energy
        self.coeffs = coeffs


class GaussianPrimitive:
    def __init__(self, exp: float, coeff: float) -> None:
        self.exp = exp
        self.coeff = coeff

        self.norm = 0.0

    def normalize(self, l: int) -> None:
        # See (Jiyun Kuang and C D Lin 1997 J. Phys. B: At. Mol. Opt. Phys. 30 2529)
        # page 2532 for the normalization factor
        self.norm = np.sqrt(2 * (2 * l) ** (l + 1.5) / gamma(l + 1.5))


class GtoShell:
    def __init__(self, atom: 'Atom', l: int, prims: list[GaussianPrimitive]) -> None:
        self.atom = atom
        self.l = l
        self.prims = prims
        self.norm = 0

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
                    * (2 * np.sqrt(i_prim.exp * j_prim.exp) / i_prim.exp + j_prim.exp) ** (self.l + 1.5)
                )

        self.norm = 1 / overlap


class Parser:
    def __init__(
        self,
        filename: Optional[str] = None,
        molden_lines: Optional[list[str]] = None,
    ) -> None:
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
        self.gto_shells = self.get_gto_shells()
        self.mo_coeffs = self.get_mos()

    def check_molden_format(self) -> None:
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

    def get_atoms(self) -> list[Atom]:
        logger.info('Parsing atoms...')
        atoms = []
        for line in self.molden_lines[self.atom_ind + 1 : self.gto_ind]:
            label, _, atomic_number, *coords = line.split()

            position = np.array([float(coord) for coord in coords], dtype=np.float64)

            atoms.append(Atom(label, int(atomic_number), position))

        logger.info(f'Parsed {len(atoms)} atoms.')
        return atoms

    def get_gto_shells(self) -> list[GtoShell]:
        logger.info('Parsing GTO lines...')

        shell_lables = ['s', 'p', 'd', 'f', 'g']

        lines = iter(self.molden_lines[self.gto_ind + 1 : self.mo_ind])

        gto_shells = []
        for atom in self.get_atoms():
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
                    prims.append(GaussianPrimitive(float(exp), float(coeff)))

                gto_shell = GtoShell(atom, shell_lables.index(shell_label), prims)
                gto_shell.normalize()

                gto_shells.append(gto_shell)

        logger.info(f'Parsed {len(gto_shells)} GTO shells.')
        return gto_shells

    def get_mos(self, mo_list: Optional[list[int]] = None) -> list[MolecularOrbital]:
        """
        Parses the molecular orbitals (MOs) from the molden file.

        Args:
            mo_list (Optional[list[int]]): A list of MO indices to parse. If None,
            all MOs will be parsed.

        Returns:
            list[MolecularOrbital]: A list of MolecularOrbital objects containing
            the symmetry, energy, and coefficients for each MO.
        """
        logger.info('Parsing MO coefficients...')

        num_atomic_orbs = sum(2 * shell.l + 1 for shell in self.gto_shells)

        if mo_list is not None:
            if not mo_list:
                raise ValueError("The provided 'mo_list' is empty.")

        order = self.atomic_orbs_order()

        lines = self.molden_lines[self.mo_ind + 1 :]
        total_num_mos = sum('Sym=' in line for line in lines)
        lines = iter(lines)

        mos = []
        for mo_ind in range(total_num_mos):
            if mo_list:
                if mo_ind not in mo_list:
                    continue

            logger.debug('Parsing MO %d', mo_ind + 1)
            _, sym = next(lines).split()

            energy_line = next(lines)
            energy = float(energy_line.split()[1])

            # Skip the next two lines which are not needed
            for _ in range(2):
                next(lines)

            coeffs = []
            for _ in range(num_atomic_orbs):
                _, coeff = next(lines).split()
                coeffs.append(float(coeff))

            mo = MolecularOrbital(
                symmetry=sym,
                energy=energy,
                coeffs=np.array(coeffs, dtype=np.float64)[order],
            )

            mos.append(mo)

        logger.info('Parsed MO coefficients.')
        return mos

    def atomic_orbs_order(self) -> None:
        """
        Return the order of the atomic orbitals in the molden file.

        Molden defines the order of the orbitals as 0, 1, -1, 2, -2, ...
        We want it to be -l, -l + 1, ..., l - 1, l.

        Note: For l = 1, the order is 1, -1, 0, which is different from the
        general pattern. This is handled separately.
        """

        order = []
        ind = 0
        for shell in self.gto_shells:
            l = shell.l
            if l == 1:
                order.extend([ind + 1, ind + 2, ind])
            else:
                order.extend([ind + i for i in range(2 * l, -1, -2)])
                order.extend([ind + i for i in range(1, 2 * l, 2)])
            ind += 2 * l + 1
