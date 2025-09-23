"""Read and parse a molden file."""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.special import gamma

logger = logging.getLogger(__name__)


@dataclass
class _Atom:
    """Represents an atom with its properties and associated shells.

    Parameters
    ----------
    label : str
        The atomic label (e.g., 'C', 'O', 'H').
    atomic_number : int
        The atomic number of the element.
    position : NDArray[np.floating]
        The 3D position coordinates of the atom.
    shells : list[_Shell]
        List of electron shells associated with this atom.
    """

    label: str
    atomic_number: int
    position: NDArray[np.floating]
    shells: list['_Shell']


@dataclass
class _MolecularOrbital:
    """Represents a molecular orbital with its properties.

    Parameters
    ----------
    sym : str
        The symmetry label of the molecular orbital.
    energy : float
        The energy of the molecular orbital.
    spin : str
        The spin state of the molecular orbital ('Alpha' or 'Beta').
    occ : int
        The occupation number of the molecular orbital.
    """

    sym: str
    energy: float
    spin: str
    occ: int


class _GTO:
    """Represents a Gaussian-type orbital (GTO) primitive.

    Parameters
    ----------
    exp : float
        The exponent of the Gaussian function.
    coeff : float
        The coefficient of the Gaussian function.
    """

    def __init__(self, exp: float, coeff: float) -> None:
        """Initialize a GTO primitive with exponent and coefficient."""
        self.exp = exp
        self.coeff = coeff

        self.norm = 0.0

    def normalize(self, l: int) -> None:
        """Normalize the GTO primitive.

        Parameters
        ----------
        l : int
            The angular momentum quantum number.

        Notes
        -----
        Uses the normalization factor from Jiyun Kuang and C D Lin 1997
        J. Phys. B: At. Mol. Opt. Phys. 30 2529, page 2532.
        """
        # See (Jiyun Kuang and C D Lin 1997 J. Phys. B: At. Mol. Opt. Phys. 30 2529)
        # page 2532 for the normalization factor
        self.norm = np.sqrt(2 * (2 * self.exp) ** (l + 1.5) / gamma(l + 1.5))


class _Shell:
    """Represents an electron shell containing multiple GTO primitives.

    Parameters
    ----------
    l : int
        The angular momentum quantum number of the shell.
    gtos : list[_GTO]
        List of GTO primitives that compose this shell.
    """

    def __init__(self, l: int, gtos: list[_GTO]) -> None:
        """Initialize a shell with angular momentum and GTO primitives."""
        self.l = l
        self.gtos = gtos

        self.norm = 0.0

    def normalize(self) -> None:
        """Normalize the shell.

        Notes
        -----
        Uses the normalization factor from Jiyun Kuang and C D Lin 1997
        J. Phys. B: At. Mol. Opt. Phys. 30 2529, equations 18 and 20.
        """
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

    Parameters
    ----------
    source : str | list[str]
        The path to the molden file, or the lines from the file.
    only_molecule : bool, optional
        Only parse the atoms and skip molecular orbitals.
        Default is `False`.

    Attributes
    ----------
    atoms : list[_Atom]
        A list of Atom objects containing the label, atomic number,
        and position for each atom.
    shells : list[_Shell]
        A list of `_Shell` objects containing the atom, angular
        momentum quantum number (l), and GTOs for each shell.
    mos : list[_MolecularOrbital]
        A list of MolecularOrbital objects containing the symmetry,
        energy, and coefficients for each MO.
    mo_coeffs : NDArray[np.floating]
        A 2D array containing all molecular orbital coefficients, where
        each row represents the coefficients for one molecular orbital.

    Raises
    ------
    TypeError
        If the source is not a valid molden file path, or molden file lines.
    """

    ANGSTROM_TO_BOHR = 1.8897259886

    def __init__(
        self,
        source: str | list[str] | Any,
        only_molecule: bool = False,
    ) -> None:
        """Initialize the Parser with either a filename or molden lines."""
        if isinstance(source, str):
            with Path(source).open('r') as file:
                self.molden_lines = file.readlines()
        elif isinstance(source, list):
            self.molden_lines = source
        else:
            raise TypeError('Source must be a filename (str) or list of lines (list[str]).')

        # Remove leading/trailing whitespace and newline characters
        self.molden_lines = [line.strip() for line in self.molden_lines]

        self.check_molden_format()

        self._atom_ind, self._gto_ind, self._mo_ind = self.divide_molden_lines()

        self.atoms = self.get_atoms()

        if only_molecule:
            return

        self.shells = self.get_shells()
        self.mos, self.mo_coeffs = self.get_mos()

    def check_molden_format(self) -> None:
        """Check if the provided molden lines conform to the expected format.

        Raises
        ------
        ValueError
            If the molden lines do not contain the required sections
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
        tuple[int, int, int]
            Indices of the '[Atoms]', '[GTO]', and '[MO]' lines.

        Raises
        ------
        ValueError
            If the molden lines do not contain the required sections.

        """
        logger.info('Dividing molden lines into sections...')
        if '[Atoms] AU' in self.molden_lines:
            atom_ind = self.molden_lines.index('[Atoms] AU')
        elif '[Atoms] Angs' in self.molden_lines:
            atom_ind = self.molden_lines.index('[Atoms] Angs')
        else:
            raise ValueError('No (AU/Angs) in [Atoms] section found in the molden file.')

        gto_ind = self.molden_lines.index('[GTO]')

        mo_ind = self.molden_lines.index('[MO]')

        logger.info('Finished dividing molden lines.')
        return atom_ind, gto_ind, mo_ind

    def get_atoms(self) -> list[_Atom]:
        """Parse the atoms from the molden file.

        Returns
        -------
        list[_Atom]
            A list of Atom objects containing the label, atomic number,
            and position for each atom.

        """
        logger.info('Parsing atoms...')
        angs = 'Angs' in self.molden_lines[self._atom_ind]

        atoms = []
        for line in self.molden_lines[self._atom_ind + 1 : self._gto_ind]:
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
        list[_Shell]
            A list of `_Shell` objects containing the atom, angular
            momentum quantum number (l), and GTOs for each shell.

        Raises
        ------
        ValueError
            If the shell label is not supported or if the GTOs are not
            formatted correctly in the molden file.

        """
        logger.info('Parsing GTO lines...')

        shell_labels = ['s', 'p', 'd', 'f', 'g']

        lines = iter(self.molden_lines[self._gto_ind + 1 : self._mo_ind])

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
                    gtos.append(_GTO(float(exp), float(coeff)))

                shell = _Shell(shell_labels.index(shell_label), gtos)
                shell.normalize()

                atom.shells.append(shell)
                shells.append(shell)

        logger.info('Parsed %s GTOs.', len(shells))
        return shells

    def get_mos(self, sort: bool = True) -> tuple[list[_MolecularOrbital], NDArray[np.floating]]:
        """Parse the molecular orbitals (MOs) from the molden file.

        Parameters
        ----------
        sort : bool, optional
            If true (default), returns the MOs sorted by energy. If false,
            returns the MOs in the order given in the molden file.

        Returns
        -------
        tuple[list[_MolecularOrbital], NDArray[np.floating]]
            Two-item tuple: the first element contains the parsed molecular
            orbitals (symmetry, energy, spin, occupation), and the second is a
            2D NumPy array of orbital coefficients shaped
            ``(num_mos, num_basis_functions)``.
        """
        logger.info('Parsing MO coefficients...')

        num_total_gtos = sum(2 * gto.l + 1 for gto in self.shells)

        order = self._gto_order()

        lines = self.molden_lines[self._mo_ind + 1 :]
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

            mo = _MolecularOrbital(
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
