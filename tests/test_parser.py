"""Unit tests for the Molden file parser."""

from pathlib import Path
from typing import Any

import numpy as np
import pytest

from tests._src_imports import _GTO, Parser, _Shell

# ----------------------------------------------------------------------
# utilities
# ----------------------------------------------------------------------
MOLDEN_PATH = Path(__file__).with_name('sample_molden.inp')


@pytest.fixture(scope='session')
def parser_obj() -> Parser:
    """
    Parser built once per test session from the reference Molden file.

    Returns
    -------
        Parser object
    """
    return Parser(str(MOLDEN_PATH))


# ----------------------------------------------------------------------
# basic structural sanity
# ----------------------------------------------------------------------
def test_section_indices_order(parser_obj: Parser) -> None:
    """Check if section indices are in the correct order."""
    assert parser_obj._atom_ind < parser_obj._gto_ind < parser_obj._mo_ind  # noqa: SLF001


def test_gaussian_normalization_positive() -> None:
    """Check if Gaussian type orbitals and shells are normalized correctly."""
    gto = _GTO(0.8, 0.5)
    gto.normalize(l=2)
    shell = _Shell(2, [gto])
    shell.normalize()
    assert gto.norm > 0.0
    assert shell.norm > 0.0


def test_atomic_orbital_permutation(parser_obj: Parser) -> None:
    """Check if the permutation of atomic orbitals is a valid one."""
    order = parser_obj._gto_order()  # noqa: SLF001
    assert sorted(order) == list(range(len(order)))


def test_atom_labels(parser_obj: Parser) -> None:
    """Check if atom labels are loaded correctly."""
    labels = [atm.label for atm in parser_obj.atoms]
    assert labels == ['Br', 'C_a', 'C_b', 'C_c', 'C_d', 'H']


def test_basis_and_mo_dimensions(parser_obj: Parser) -> None:
    """Check number of MOs and GTOs against known values."""
    num_mos = 177
    assert len(parser_obj.mos) == num_mos

    num_gtos = sum(2 * shell.l + 1 for shell in parser_obj.shells)
    # Check that mo_coeffs has the right shape
    assert parser_obj.mo_coeffs.shape == (num_mos, num_gtos)


def test_mo_energies_are_sorted(parser_obj: Parser) -> None:
    """Molecular orbital energies must be sorted in ascending order."""
    energies = np.asarray([mo.energy for mo in parser_obj.mos])
    assert np.all(np.diff(energies) >= 0.0)


@pytest.mark.parametrize('source', [None, 1, 1.0, {}, set()])
def test_parser_invalid_input_type(source: Any) -> None:
    """
    Parser must raise TypeError if input is not str or list of str.

    Raises
    ------
        TypeError
    """
    with pytest.raises(TypeError):
        Parser(source)


# ----------------------------------------------------------------------
# reproducibility checks
# ----------------------------------------------------------------------
def test_file_vs_lines_consistency(tmp_path: Path) -> None:
    """Parsing via filename or via pre-read lines must give identical results."""
    lines = MOLDEN_PATH.read_text().splitlines(True)

    p_from_lines = Parser(lines)

    tmp_file = tmp_path / 'copy.molden'
    tmp_file.write_text(''.join(lines))
    p_from_file = Parser(str(tmp_file))

    # Quick invariants - if these match, deeper structures are identical
    assert [a.atomic_number for a in p_from_lines.atoms] == [a.atomic_number for a in p_from_file.atoms]
    assert [mo.energy for mo in p_from_lines.mos] == [mo.energy for mo in p_from_file.mos]
