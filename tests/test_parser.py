# ruff: noqa
from pathlib import Path

import numpy as np
import pytest

from moldenViz.parser import Parser, _GTO, _Shell

# ----------------------------------------------------------------------
# utilities
# ----------------------------------------------------------------------
MOLDEN_PATH = Path(__file__).with_name('sample_molden.inp')


@pytest.fixture(scope='session')
def parser_obj() -> Parser:
    """Parser built once per test session from the reference Molden file."""
    return Parser(str(MOLDEN_PATH))


# ----------------------------------------------------------------------
# basic structural sanity
# ----------------------------------------------------------------------
def test_section_indices_order(parser_obj: Parser) -> None:
    assert parser_obj.atom_ind < parser_obj.gto_ind < parser_obj.mo_ind


def test_gaussian_normalization_positive() -> None:
    gto = _GTO(0.8, 0.5)
    gto.normalize(l=2)
    shell = _Shell(2, [gto])
    shell.normalize()
    assert gto.norm > 0.0 and shell.norm > 0.0


def test_atomic_orbital_permutation(parser_obj: Parser) -> None:
    order = parser_obj._gto_order()
    assert sorted(order) == list(range(len(order)))


def test_atom_labels(parser_obj: Parser) -> None:
    labels = [atm.label for atm in parser_obj.atoms]
    assert labels == ['Br', 'C_a', 'C_b', 'C_c', 'C_d', 'H']


def test_basis_and_mo_dimensions(parser_obj: Parser) -> None:
    n_basis = sum(2 * shell.l + 1 for shell in parser_obj.shells)

    assert len(parser_obj.mos) == 177

    # every MO coefficient vector must have that length
    for mo in parser_obj.mos:
        assert len(mo.coeffs) == n_basis


def test_mo_energies_are_sorted(parser_obj: Parser) -> None:
    energies = np.asarray([mo.energy for mo in parser_obj.mos])
    assert np.all(np.diff(energies) >= 0.0)


# ----------------------------------------------------------------------
# reproducibility checks
# ----------------------------------------------------------------------
def test_file_vs_lines_consistency(tmp_path) -> None:
    """Parsing via filename or via pre-read lines must give identical results."""
    lines = MOLDEN_PATH.read_text().splitlines(True)

    p_from_lines = Parser(lines)

    tmp_file = tmp_path / 'copy.molden'
    tmp_file.write_text(''.join(lines))
    p_from_file = Parser(str(tmp_file))

    # Quick invariants - if these match, deeper structures are identical
    assert [a.atomic_number for a in p_from_lines.atoms] == [a.atomic_number for a in p_from_file.atoms]
    assert [mo.energy for mo in p_from_lines.mos] == [mo.energy for mo in p_from_file.mos]
