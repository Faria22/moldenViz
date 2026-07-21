"""Read and parse a molden file."""

import logging
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from .models import Atom, GaussianPrimitive, MolecularOrbital, Shell

__all__ = ['BOHR_PER_ANGSTROM', 'Atom', 'GaussianPrimitive', 'MolecularOrbital', 'Parser', 'Shell']

BOHR_PER_ANGSTROM = 1.8897259886

logger = logging.getLogger(__name__)


class Parser:
    """Parser for molden files.

    Parameters
    ----------
    source : str | list[str]
        The path to the molden file, or the lines from the file.
    only_molecule : bool, optional
        Only parse the atoms and skip molecular orbitals.
        Default is `False`.
    mo_order : {'energy', 'file'}, optional
        Order molecular orbitals by ascending energy or preserve their order
        in the Molden file. Default is ``'energy'``.

    Attributes
    ----------
    atoms : list[Atom]
        A list of Atom objects containing the label, atomic number,
        and position for each atom.
    shells : list[Shell]
        A list of `Shell` objects containing the atom, angular
        momentum quantum number (l), and GTOs for each shell.
    mos : list[MolecularOrbital]
        A list of MolecularOrbital objects containing the symmetry,
        energy, and coefficients for each MO.
    mo_coeffs : NDArray[np.floating]
        A 2D array containing all molecular orbital coefficients, where
        each row represents the coefficients for one molecular orbital.

    Raises
    ------
    TypeError
        If the source is not a valid molden file path, or molden file lines.
    ValueError
        If ``mo_order`` is not ``'energy'`` or ``'file'``.
    """

    def __init__(
        self,
        source: str | list[str],
        only_molecule: bool = False,
        mo_order: Literal['energy', 'file'] = 'energy',
    ) -> None:
        """Initialize the Parser with either a filename or molden lines."""
        if mo_order not in {'energy', 'file'}:
            raise ValueError("'mo_order' must be either 'energy' or 'file'.")

        if isinstance(source, str):
            with Path(source).open('r') as file:
                self._molden_lines = file.readlines()
        elif isinstance(source, list):
            self._molden_lines = source
        else:
            raise TypeError('Source must be a filename (str) or list of lines (list[str]).')

        # Remove leading/trailing whitespace and newline characters
        self._molden_lines = [line.strip() for line in self._molden_lines]

        self._check_molden_format()

        self._atom_ind, self._gto_ind, self._mo_ind = self._divide_molden_lines()

        self.atoms = self._parse_atoms()
        self.shells: list[Shell] = []
        self.mos: list[MolecularOrbital] = []
        self.mo_coeffs: NDArray[np.floating] = np.empty((0, 0), dtype=float)

        if only_molecule:
            return

        self.shells = self._parse_shells()
        self.mos, self.mo_coeffs = self._parse_mos(sort=mo_order == 'energy')

    def _check_molden_format(self) -> None:
        """Check if the provided molden lines conform to the expected format.

        Raises
        ------
        ValueError
            If the molden lines do not contain the required sections
            or if they are in an unsupported format.

        """
        logger.info('Checking molden format...')
        if not self._molden_lines:
            raise ValueError('The provided molden lines are empty.')

        if not any('[Atoms]' in line for line in self._molden_lines):
            raise ValueError("No '[Atoms]' section found in the molden file.")

        if not any('[GTO]' in line for line in self._molden_lines):
            raise ValueError("No '[GTO]' section found in the molden file.")

        if not any('[MO]' in line for line in self._molden_lines):
            raise ValueError("No '[MO]' section found in the molden file.")

        if not any(orbs in line for orbs in ['5D', '9G'] for line in self._molden_lines):
            raise ValueError('Cartesian orbitals functions are not currently supported.')

        logger.info('Molden format check passed.')

    def _divide_molden_lines(self) -> tuple[int, int, int]:
        """Divide the molden lines into sections for atoms, GTOs, and MOs.

        Returns
        -------
        tuple[int, int, int]
            Indices of the '[Atoms]', '[GTO]', and '[MO]' lines.

        Raises
        ------
        ValueError
            If the molden lines do not contain the required sections.

        """
        logger.info('Dividing molden lines into sections...')
        if '[Atoms] AU' in self._molden_lines:
            atom_ind = self._molden_lines.index('[Atoms] AU')
        elif '[Atoms] Angs' in self._molden_lines:
            atom_ind = self._molden_lines.index('[Atoms] Angs')
        else:
            raise ValueError('No (AU/Angs) in [Atoms] section found in the molden file.')

        gto_ind = self._molden_lines.index('[GTO]')

        mo_ind = self._molden_lines.index('[MO]')

        logger.info('Finished dividing molden lines.')
        return atom_ind, gto_ind, mo_ind

    def _parse_atoms(self) -> list[Atom]:
        """Parse the atoms from the molden file.

        Returns
        -------
        list[Atom]
            A list of Atom objects containing the label, atomic number,
            and position for each atom.

        """
        logger.info('Parsing atoms...')
        angs = 'Angs' in self._molden_lines[self._atom_ind]

        atoms = []
        for line in self._molden_lines[self._atom_ind + 1 : self._gto_ind]:
            label, _, atomic_number, *coords = line.split()

            position = np.array([float(coord) for coord in coords], dtype=float)
            if angs:
                position *= BOHR_PER_ANGSTROM

            atoms.append(Atom(label, int(atomic_number), position, []))

        logger.info('Parsed %s atoms.', len(atoms))
        return atoms

    def _parse_shells(self) -> list[Shell]:
        """Parse the Gaussian-type orbitals (GTOs) from the molden file.

        Returns
        -------
        list[Shell]
            A list of `Shell` objects containing the atom, angular
            momentum quantum number (l), and GTOs for each shell.

        Raises
        ------
        ValueError
            If the shell label is not supported or if the GTOs are not
            formatted correctly in the molden file.

        """
        logger.info('Parsing GTO lines...')

        shell_labels = ['s', 'p', 'd', 'f', 'g']

        lines = iter(self._molden_lines[self._gto_ind + 1 : self._mo_ind])

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
                if shell_label not in shell_labels:
                    raise ValueError(f"Shell label '{shell_label}' is currently not supported.")

                gtos = []
                for _ in range(int(num_gtos)):
                    exp, coeff = next(lines).split()
                    gtos.append(GaussianPrimitive(float(exp), float(coeff)))

                shell = Shell(shell_labels.index(shell_label), gtos)
                shell._normalize()  # ruff:ignore[private-member-access]

                atom.shells.append(shell)
                shells.append(shell)

        logger.info('Parsed %s GTOs.', len(shells))
        return shells

    def _parse_mos(self, sort: bool = True) -> tuple[list[MolecularOrbital], NDArray[np.floating]]:
        """Parse the molecular orbitals (MOs) from the molden file.

        Parameters
        ----------
        sort : bool, optional
            If true (default), returns the MOs sorted by energy. If false,
            returns the MOs in the order given in the molden file.

        Returns
        -------
        tuple[list[MolecularOrbital], NDArray[np.floating]]
            Two-item tuple: the first element contains the parsed molecular
            orbitals (symmetry, energy, spin, occupation), and the second is a
            2D NumPy array of orbital coefficients shaped
            ``(num_mos, num_basis_functions)``.
        """
        logger.info('Parsing MO coefficients...')

        num_total_gtos = sum(2 * gto.l + 1 for gto in self.shells)

        order = self._gto_order()

        lines = self._molden_lines[self._mo_ind + 1 :]
        total_num_mos = sum('Sym=' in line for line in lines)

        lines = iter(lines)

        mos = []
        mo_coeffs = np.empty((total_num_mos, num_total_gtos), dtype=float)

        for mo_ind in range(total_num_mos):
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

            # Store coefficients in shared array
            mo_coeffs[mo_ind] = np.array(coeffs, dtype=float)[order]

            mo = MolecularOrbital(
                sym=sym,
                energy=energy,
                spin=spin,
                occ=occ,
            )

            mos.append(mo)

        logger.info('Parsed MO coefficients.')

        if sort:
            # Sort MOs and reorder mo_coeffs to match
            sorted_indices = sorted(range(len(mos)), key=lambda i: mos[i].energy)
            mos = [mos[i] for i in sorted_indices]
            mo_coeffs = mo_coeffs[sorted_indices]

        return mos, mo_coeffs

    def _gto_order(self) -> list[int]:
        """Return the order of the GTOs in the molden file.

        Molden defines the order of the orbitals as m = 0, 1, -1, 2, -2, ...
        We want it to be m = -l, -l + 1, ..., l - 1, l.

        Note: For l = 1, the order is 1, -1, 0, which is different from the
        general pattern. This is handled separately.

        Returns
        -------
        list[int]
            The order of the atomic orbitals.

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
