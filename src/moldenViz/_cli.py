"""CLI entrance point."""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Sequence

from .__about__ import __version__
from .examples.get_example_files import all_examples
from .plotter import Plotter


def _configure_logging(level: int) -> None:
    """Configure the root logger with the requested level."""

    logging.basicConfig(level=level)
    logging.getLogger().setLevel(level)


def main(argv: Sequence[str] | None = None) -> None:
    """Entry point for the moldenViz command-line interface.

    Parses command line arguments and launches the plotter with the specified
    molden file or example molecule. Supports options to plot only the molecule
    structure without molecular orbitals and to configure logging verbosity.
    """

    parser = argparse.ArgumentParser(prog='moldenViz')
    parser.set_defaults(log_level=None)
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument(
        '-v',
        '--verbose',
        dest='log_level',
        action='store_const',
        const=logging.INFO,
        help='Enable info level logging output.',
    )
    verbosity.add_argument(
        '-d',
        '--debug',
        dest='log_level',
        action='store_const',
        const=logging.DEBUG,
        help='Enable debug level logging output.',
    )
    verbosity.add_argument(
        '-q',
        '--quiet',
        dest='log_level',
        action='store_const',
        const=logging.ERROR,
        help='Limit logging output to errors only.',
    )

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

    args = parser.parse_args(sys.argv[1:] if argv is None else argv)

    log_level = args.log_level if args.log_level is not None else logging.WARNING
    _configure_logging(log_level)

    Plotter(
        args.file or all_examples[args.example],
        only_molecule=args.only_molecule,
    )


if __name__ == '__main__':
    main()
