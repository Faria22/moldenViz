"""CLI entrance point."""

import argparse
import logging

from .__about__ import __version__
from .examples.get_example_files import all_examples
from .plotter import Plotter

logger = logging.getLogger(__name__)


class ColoredFormatter(logging.Formatter):
    """Logging formatter that adds color codes to log levels."""

    COLORS = {
        'ERROR': '\033[38;5;196m',  # Bright Red
        'WARNING': '\033[38;5;208m',  # Bright Orange
        'INFO': '\033[38;5;34m',  # Bright Green
        'DEBUG': '\033[38;5;27m',  # Bright Blue
        'RESET': '\033[0m',  # Reset to default color
    }

    def format(self, record: logging.LogRecord) -> str:
        """Format the log record with colors.

        Parameters
        ----------
        record : logging.LogRecord
            The log record to format.

        Returns
        -------
        str
            The formatted log message with color codes.
        """
        # Get color for this level, default to no color
        color = self.COLORS.get(record.levelname, '')
        reset = self.COLORS['RESET'] if color else ''

        # Store original levelname
        original_levelname = record.levelname

        # Add color to levelname if we have a color for this level
        if color:
            record.levelname = f'{color}{record.levelname}{reset}'

        # Format the message
        result = super().format(record)

        # Restore original levelname
        record.levelname = original_levelname

        return result


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

    verbosity_group = parser.add_argument_group('verbosity')
    verbosity = verbosity_group.add_mutually_exclusive_group()
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

    # Configure logging with colored formatter
    handler = logging.StreamHandler()
    handler.setFormatter(ColoredFormatter('%(levelname)s %(name)s: %(message)s'))
    logging.basicConfig(level=level, handlers=[handler], force=True)

    logger.debug('Parsed CLI arguments: %s', vars(args))

    source_path = args.file or all_examples[args.example]
    source_label = args.file or f'example {args.example}'
    logger.info('Launching plotter for %s', source_label)

    Plotter(
        source_path,
        only_molecule=args.only_molecule,
    )


if __name__ == '__main__':
    main()
