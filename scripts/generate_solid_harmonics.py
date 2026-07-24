"""Generate direct NumPy real solid-harmonic kernels through ``l=4``.

Run from the repository root with::

    uv run --with sympy==1.14.0 python scripts/generate_solid_harmonics.py

SymPy is a development-only code-generation dependency. The generated runtime
module depends only on NumPy.
"""

from __future__ import annotations

import argparse
from math import factorial
from pathlib import Path
from typing import TYPE_CHECKING

import sympy as sp

if TYPE_CHECKING:
    from collections.abc import Iterable

DEFAULT_OUTPUT = Path('src/moldenViz/_generated_solid_harmonics.py')
LMAX = 4


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--output',
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f'Generated module path (default: {DEFAULT_OUTPUT}).',
    )
    parser.add_argument(
        '--check',
        action='store_true',
        help='Fail instead of writing when the committed module is stale.',
    )
    return parser.parse_args()


def _solid_harmonic(
    x: sp.Symbol,
    y: sp.Symbol,
    z: sp.Symbol,
    l: int,
    m: int,
) -> sp.Expr:
    """Build one normalized real solid-harmonic polynomial.

    Returns
    -------
    sp.Expr
        Expanded symbolic polynomial.
    """
    absolute_m = abs(m)
    xy_polynomial = sp.expand((x + sp.I * y) ** absolute_m)
    xy_factor = sp.re(xy_polynomial) if m >= 0 else sp.im(xy_polynomial)
    r_sq = x**2 + y**2 + z**2

    legendre_factor = 0
    for k in range((l - absolute_m) // 2 + 1):
        coefficient = (
            (-1) ** k
            * factorial(2 * l - 2 * k)
            / (2**l * factorial(k) * factorial(l - k) * factorial(l - 2 * k - absolute_m))
        )
        legendre_factor += coefficient * z ** (l - absolute_m - 2 * k) * r_sq**k

    normalization = sp.sqrt(
        sp.Rational((2 * l + 1) * factorial(l - absolute_m), 4 * factorial(l + absolute_m)) / sp.pi,
    )
    if absolute_m:
        normalization *= sp.sqrt(2)
    return sp.expand(normalization * xy_factor * legendre_factor)


def _floating_polynomial(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, sp.Symbol, sp.Symbol],
) -> sp.Expr:
    """Replace exact symbolic coefficients with deterministic binary floats.

    Returns
    -------
    sp.Expr
        Polynomial with floating-point coefficients.
    """
    polynomial = sp.Poly(expression, *variables)
    result = 0
    for powers, coefficient in polynomial.terms():
        term = sp.Float(coefficient.evalf(17), 17)
        for variable, power in zip(variables, powers, strict=True):
            term *= variable**power
        result += term
    return sp.expand(result)


def _python(expression: sp.Expr) -> str:
    return sp.pycode(expression, standard='python3')


def _function_lines(lmax: int, expressions: Iterable[tuple[int, int, sp.Expr]]) -> list[str]:
    indexed_expressions = list(expressions)
    replacements, reduced = sp.cse(
        [expression for _, _, expression in indexed_expressions],
        symbols=sp.numbered_symbols('cse'),
        optimizations='basic',
        order='canonical',
    )
    replacement_symbols = {symbol for symbol, _ in replacements}
    replacement_map = dict(replacements)

    def required_replacements(expression: sp.Expr) -> set[sp.Symbol]:
        required = expression.free_symbols & replacement_symbols
        pending = list(required)
        while pending:
            symbol = pending.pop()
            dependencies = replacement_map[symbol].free_symbols & replacement_symbols
            new_dependencies = dependencies - required
            required.update(new_dependencies)
            pending.extend(new_dependencies)
        return required

    operations: list[tuple[str, str | sp.Symbol, sp.Expr]] = []
    emitted: set[sp.Symbol] = set()
    for (l, m, _), expression in zip(indexed_expressions, reduced, strict=True):
        required = required_replacements(expression)
        for symbol, replacement in replacements:
            if symbol in required and symbol not in emitted:
                operations.append(('replacement', symbol, replacement))
                emitted.add(symbol)
        operations.append(('output', f'values[{l}, {m}, :]', expression))

    last_uses: dict[sp.Symbol, int] = {}
    for operation_index, (_, _, expression) in enumerate(operations):
        for symbol in expression.free_symbols & replacement_symbols:
            last_uses[symbol] = operation_index

    lines = [
        f'def _tabulate_l{lmax}(x: FloatArray, y: FloatArray, z: FloatArray) -> FloatArray:',
        f'    values = np.zeros(({lmax + 1}, {2 * lmax + 1}, x.size), dtype=float)',
    ]
    for operation_index, (_, target, expression) in enumerate(operations):
        lines.append(f'    {target} = {_python(expression)}')
        expired = sorted(str(symbol) for symbol, last_use in last_uses.items() if last_use == operation_index)
        if expired:
            lines.append(f'    del {", ".join(expired)}')
    lines.extend(['    return values', '', ''])
    return lines


def generate_module() -> str:
    """Return the deterministic generated module contents.

    Returns
    -------
    str
        Complete Python source for the generated runtime module.
    """
    x, y, z = sp.symbols('x y z', real=True)
    variables = (x, y, z)
    expressions = {
        (l, m): _floating_polynomial(_solid_harmonic(x, y, z, l, m), variables)
        for l in range(LMAX + 1)
        for m in range(-l, l + 1)
    }

    lines = [
        '"""Generated direct real solid-harmonic kernels. Do not edit by hand."""',
        '',
        (
            '# ruff: file-ignore[unused-function-argument, missing-whitespace-around-arithmetic-operator, '
            'line-too-long, docstring-missing-returns]'
        ),
        '# fmt: off',
        '',
        'import numpy as np',
        'from numpy.typing import NDArray',
        '',
        'FloatArray = NDArray[np.floating]',
        '',
        '',
    ]
    for lmax in range(LMAX + 1):
        selected = [(l, m, expressions[l, m]) for l in range(lmax + 1) for m in range(-l, l + 1)]
        lines.extend(_function_lines(lmax, selected))

    lines.extend(
        [
            '_GENERATED_KERNELS = (_tabulate_l0, _tabulate_l1, _tabulate_l2, _tabulate_l3, _tabulate_l4)',
            '',
            '',
            'def tabulate_solid_harmonics(',
            '    x: FloatArray,',
            '    y: FloatArray,',
            '    z: FloatArray,',
            '    lmax: int,',
            ') -> FloatArray:',
            '    """Evaluate generated real solid harmonics through ``lmax``."""',
            '    return _GENERATED_KERNELS[lmax](x, y, z)',
            '',
        ],
    )
    return '\n'.join(lines)


def main() -> None:
    """Write the generated runtime module."""
    args = _parse_args()
    generated = generate_module()
    if args.check:
        if not args.output.exists() or args.output.read_text() != generated:
            raise SystemExit(f'{args.output} is stale; rerun this script without --check.')
        print(f'{args.output} is up to date')  # ruff: ignore[print]
        return

    args.output.write_text(generated)
    print(f'wrote {args.output}')  # ruff: ignore[print]


if __name__ == '__main__':
    main()
