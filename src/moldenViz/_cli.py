"""CLI entrance point."""

import argparse
import logging

from moldenViz.__about__ import __version__

from .examples.get_example_files import all_examples
from .plotter import Plotter


def main() -> None:
    """Entry point for the moldenViz command-line interface.

    Parses command line arguments and launches the plotter with the specified
    molden file or example molecule. Supports options to plot only the molecule
    structure without molecular orbitals.
    """
    parser = argparse.ArgumentParser(prog='moldenViz')
    source = parser.add_mutually_exclusive_group(required=True)

    source.add_argument('file', nargs='?', default=None, help='Optional molden file path', type=str)
    parser.add_argument('-m', '--only_molecule', action='store_true', help='Only plots the molecule')
    source.add_argument(
        '-e',
        '--example',
        type=str,
        metavar='molecule',
        choices=all_examples.keys(),
        help='Load example %(metavar)s. Options are: %(choices)s',
    )
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')

    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument('-v', '--verbose', action='store_true', help='Increase logging verbosity to INFO')
    verbosity.add_argument('-d', '--debug', action='store_true', help='Enable debug logging')
    verbosity.add_argument('-q', '--quiet', action='store_true', help='Reduce logging output to errors only')

    args = parser.parse_args()

    if args.debug:
        level = logging.DEBUG
    elif args.verbose:
        level = logging.INFO
    elif args.quiet:
        level = logging.ERROR
    else:
        level = logging.WARNING

    logging.basicConfig(level=level, format='%(levelname)s %(name)s: %(message)s', force=True)

    Plotter(
        args.file or all_examples[args.example],
        only_molecule=args.only_molecule,
    )


if __name__ == '__main__':
    main()
