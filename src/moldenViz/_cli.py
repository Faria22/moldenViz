"""CLI entrance point."""

import argparse

from .examples.get_example_files import all_examples
from .plotter import Plotter


def main() -> None:
    parser = argparse.ArgumentParser(prog='moldenViz')
    source = parser.add_mutually_exclusive_group(required=True)

    source.add_argument('file', nargs='?', default=None, help='Optional molden file path', type=str)
    parser.add_argument('-m', '--only_molecule', action='store_true', help='Only plots the molecule')
    source.add_argument(
        '-e',
        '--example',
        nargs='?',
        type=str,
        metavar='molecule',
        choices=all_examples.keys(),
        help='Load example %(metavar)s. Options are: %(choices)s',
    )

    args = parser.parse_args()

    Plotter(
        args.file or all_examples[args.example],
        only_molecule=args.only_molecule,
    )


if __name__ == '__main__':
    main()
