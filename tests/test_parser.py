"""Unit tests for the Molden file parser."""

from pathlib import Path
from typing import Any

import numpy as np
import pytest

import moldenViz.parser as parser_module
from moldenViz import GaussianPrimitive, Shell

Parser = parser_module.Parser

# ----------------------------------------------------------------------
# utilities
# ----------------------------------------------------------------------
MOLDEN_PATH = Path(__file__).with_name('sample_molden.inp')
PYSCF_SPHERICAL_PATH = Path(__file__).parent / 'fixtures/pyscf/co-cc-pvqz-spherical.molden'
CO_ATOM_COUNT = 2
PYSCF_LMAX = 4


@pytest.fixture(scope='session')
def parser_obj() -> Parser:
    """
    Parser built once per test session from the reference Molden file.

    Returns
    -------
        Parser object
    """
    return Parser(filename=MOLDEN_PATH)


# ----------------------------------------------------------------------
# basic structural sanity
# ----------------------------------------------------------------------
def test_section_indices_order(parser_obj: Parser) -> None:
    """Check if section indices are in the correct order."""
    assert parser_obj._atom_ind < parser_obj._gto_ind < parser_obj._mo_ind  # ruff:ignore[private-member-access]


def test_gaussian_normalization_positive() -> None:
    """Check if Gaussian type orbitals and shells are normalized correctly."""
    gto = GaussianPrimitive(0.8, 0.5)
    gto._normalize(l=2)  # ruff:ignore[private-member-access]
    shell = Shell(2, [gto])
    shell._normalize()  # ruff:ignore[private-member-access]
    assert gto._norm > 0.0  # ruff:ignore[private-member-access]
    assert shell._norm > 0.0  # ruff:ignore[private-member-access]


def test_atomic_orbital_permutation(parser_obj: Parser) -> None:
    """Check if the permutation of atomic orbitals is a valid one."""
    order = parser_obj._gto_order()  # ruff:ignore[private-member-access]
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


def test_mo_order_can_preserve_file_order() -> None:
    """The public constructor option should control molecular-orbital order."""
    lines = MOLDEN_PATH.read_text().splitlines(True)
    mo_start = next(index for index, line in enumerate(lines) if line.strip() == '[MO]')
    prefix = lines[: mo_start + 1]
    mo_lines = lines[mo_start + 1 :]
    block_starts = [index for index, line in enumerate(mo_lines) if line.strip().startswith('Sym=')]
    block_length = block_starts[1] - block_starts[0]
    blocks = [mo_lines[start : start + block_length] for start in block_starts]
    reversed_lines = prefix + [line for block in reversed(blocks) for line in block]

    original_file_ordered = Parser(content=''.join(lines), mo_order='file')
    reversed_content = ''.join(reversed_lines)
    file_ordered = Parser(content=reversed_content, mo_order='file')
    energy_ordered = Parser(content=reversed_content, mo_order='energy')

    file_energies = [mo.energy for mo in file_ordered.mos]
    assert file_energies == list(reversed([mo.energy for mo in original_file_ordered.mos]))
    assert [mo.energy for mo in energy_ordered.mos] == sorted(file_energies)


def test_invalid_mo_order_is_rejected() -> None:
    """Only the documented molecular-orbital order values are accepted."""
    with pytest.raises(ValueError, match='mo_order'):
        Parser(filename=MOLDEN_PATH, mo_order='input')  # type: ignore[arg-type]


@pytest.mark.parametrize('filename', [1, 1.0, {}, set()])
def test_parser_rejects_invalid_filename_type(filename: Any) -> None:
    """Filename input must be a string or path-like object."""
    with pytest.raises(TypeError, match='filename'):
        Parser(filename=filename)


@pytest.mark.parametrize('content', [1, 1.0, {}, set(), ['line']])
def test_parser_rejects_invalid_content_type(content: Any) -> None:
    """Content input must be one complete string, not pre-split lines."""
    with pytest.raises(TypeError, match='content'):
        Parser(content=content)


def test_parser_requires_exactly_one_input() -> None:
    """Filename and content inputs must be mutually exclusive and required."""
    with pytest.raises(ValueError, match='Exactly one'):
        Parser()
    with pytest.raises(ValueError, match='Exactly one'):
        Parser(filename=MOLDEN_PATH, content=MOLDEN_PATH.read_text(encoding='utf-8'))


def test_parser_rejects_removed_positional_input() -> None:
    """The removed positional source argument must not remain callable."""
    with pytest.raises(TypeError):
        Parser(str(MOLDEN_PATH))  # type: ignore[misc]


# ----------------------------------------------------------------------
# reproducibility checks
# ----------------------------------------------------------------------
def test_filename_vs_content_consistency(tmp_path: Path) -> None:
    """Parsing via filename or complete content must give identical results."""
    content = MOLDEN_PATH.read_text(encoding='utf-8')
    p_from_content = Parser(content=content)

    tmp_file = tmp_path / 'copy.molden'
    tmp_file.write_text(content, encoding='utf-8')
    p_from_file = Parser(filename=tmp_file)

    assert [a.atomic_number for a in p_from_content.atoms] == [a.atomic_number for a in p_from_file.atoms]
    np.testing.assert_allclose(
        [atom.position for atom in p_from_content.atoms],
        [atom.position for atom in p_from_file.atoms],
    )
    assert [shell.l for shell in p_from_content.shells] == [shell.l for shell in p_from_file.shells]
    assert p_from_content.mos == p_from_file.mos
    np.testing.assert_allclose(p_from_content.mo_coeffs, p_from_file.mo_coeffs)


def test_empty_content_uses_existing_format_validation() -> None:
    """Empty content should fail as an empty Molden document."""
    with pytest.raises(ValueError, match='empty'):
        Parser(content='')


def test_only_molecule_has_stable_result_attributes() -> None:
    """Molecule-only parsing should expose empty orbital result containers."""
    parser = Parser(filename=MOLDEN_PATH, only_molecule=True)

    assert parser.shells == []
    assert parser.mos == []
    assert parser.mo_coeffs.shape == (0, 0)


def test_only_molecule_accepts_atoms_without_orbital_sections() -> None:
    """Molecule-only parsing should require only a valid atom section."""
    parser = Parser(
        content="""[Molden Format]
[Atoms] Angs
H 1 1 1.0 2.0 3.0
""",
        only_molecule=True,
    )

    assert len(parser.atoms) == 1
    np.testing.assert_allclose(parser.atoms[0].position, np.array([1.0, 2.0, 3.0]) * parser_module.BOHR_PER_ANGSTROM)
    assert parser.shells == []
    assert parser.mos == []
    assert parser.mo_coeffs.shape == (0, 0)


def test_only_molecule_ignores_invalid_orbital_sections() -> None:
    """Molecule-only parsing should not validate or parse trailing orbital data."""
    parser = Parser(
        content="""[Molden Format]
[Atoms] AU
H 1 1 0.0 0.0 0.0
[GTO]
invalid basis data
[MO]
invalid molecular orbital data
""",
        only_molecule=True,
    )

    assert len(parser.atoms) == 1
    assert parser.shells == []
    assert parser.mos == []
    assert parser.mo_coeffs.shape == (0, 0)


def test_full_parser_still_requires_orbital_sections() -> None:
    """Full parsing should retain its existing section requirements."""
    content = """[Molden Format]
[Atoms] AU
H 1 1 0.0 0.0 0.0
"""

    with pytest.raises(ValueError, match=r"No '\[GTO\]' section"):
        Parser(content=content)


def test_pyscf_spherical_fixture() -> None:
    """PySCF's parenthesized atom units and lowercase basis tags should parse."""
    parser = Parser(filename=PYSCF_SPHERICAL_PATH)

    assert len(parser.atoms) == CO_ATOM_COUNT
    assert max(shell.l for shell in parser.shells) == PYSCF_LMAX
    assert parser.mo_coeffs.shape == (12, 110)
