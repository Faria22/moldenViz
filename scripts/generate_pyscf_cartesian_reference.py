"""Generate trusted PySCF Cartesian AO values for the issue 117 tests.

Run from the repository root with::

    uv run --with pyscf==2.14.0 python scripts/generate_pyscf_cartesian_reference.py

PySCF evaluates Cartesian AOs in its internal component order and without
per-component normalization. This script applies the same normalization and
Molden ordering used by ``pyscf.tools.molden.orbital_coeff``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pyscf
from pyscf.tools import molden

DEFAULT_SOURCE = Path('tests/fixtures/pyscf/co-cc-pvqz-cartesian.molden')
DEFAULT_OUTPUT = Path('tests/fixtures/pyscf/co-cc-pvqz-cartesian-reference.json')
REFERENCE_POINTS = np.array(
    [
        [0.0, 0.0, 0.0],
        [0.25, -0.5, 1.0],
        [-1.5, 0.75, -0.25],
        [2.0, -1.0, 0.5],
        [1e-10, 0.0, 1.0],
        [0.0, -1e-10, -1.0],
    ],
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    """Generate the normalized AO reference in Molden component order."""
    args = _parse_args()
    molecule, _, mo_coefficients, *_ = molden.load(str(args.source))
    if not molecule.cart:
        raise ValueError(f'{args.source} does not declare a Cartesian basis.')

    pyscf_values = molecule.eval_gto('GTOval_cart', REFERENCE_POINTS)
    overlap_norms = np.sqrt(molecule.intor('int1e_ovlp').diagonal())
    normalized_values = pyscf_values / overlap_norms
    molden_order = np.asarray(molden.order_ao_index(molecule))
    molden_values = normalized_values[:, molden_order]
    mo_values = pyscf_values @ mo_coefficients

    payload = {
        'ao_values_molden_order': molden_values.tolist(),
        'mo_values': mo_values.tolist(),
        'points_bohr': REFERENCE_POINTS.tolist(),
        'pyscf_version': pyscf.__version__,
        'source': str(args.source),
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n')
    print(f'wrote {args.output}')  # ruff: ignore[print]


if __name__ == '__main__':
    main()
