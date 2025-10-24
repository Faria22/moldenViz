"""CLI entrance point."""

import argparse
import logging
from typing import Any

from .__about__ import __version__
from .examples.get_example_files import all_examples
from .plotter import Plotter

logger = logging.getLogger(__name__)


def main() -> None:
    """Entry point for the moldenViz command-line interface.

    Parses command line arguments and launches the plotter with the specified
    molden file or example molecule. Supports options to plot only the molecule
    structure without molecular orbitals.
    """
    parser = argparse.ArgumentParser(prog='moldenViz')
    parser.add_argument('-V', '--version', action='version', version=f'%(prog)s {__version__}')
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

    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument('-v', '--verbose', action='store_true', help='Enable INFO level logging output')
    verbosity.add_argument('-d', '--debug', action='store_true', help='Enable DEBUG level logging output')
    verbosity.add_argument('-q', '--quiet', action='store_true', help='Reduce logging output to errors only')

    args = parser.parse_args()

    level = logging.WARNING
    if args.debug:
        level = logging.DEBUG
    elif args.verbose:
        level = logging.INFO
    elif args.quiet:
        level = logging.ERROR

    logging.basicConfig(
        level=level,
        format='%(levelname)s %(name)s: %(message)s',
        force=True,
    )

    logger.debug('Parsed CLI arguments: %s', _redact_args(vars(args)))

    source_path = args.file or all_examples[args.example]
    source_label = args.file if args.file else f'example {args.example}'
    logger.info('Launching plotter for %s', source_label)

    Plotter(
        source_path,
        only_molecule=args.only_molecule,
    )


def _redact_args(arguments: dict[str, Any]) -> dict[str, Any]:
    """Return a shallow copy of CLI arguments with filesystem paths redacted."""

    redacted = arguments.copy()
    if redacted.get('file'):
        redacted['file'] = '<redacted path>'
    return redacted


if __name__ == '__main__':
    main()
