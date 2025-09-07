"""Get example files from folder and make them available to the package.

This module provides access to pre-loaded molecular structures for demonstration
and testing purposes. All examples are stored as lists of lines from molden files.
"""

from pathlib import Path


def _read_file(path: Path) -> list[str]:
    """Read a molden file and return its contents as a list of lines.

    Parameters
    ----------
    path : Path
        Path to the molden file to read.

    Returns
    -------
    list[str]
        List of lines from the molden file.
    """
    with path.open('r') as f:
        return f.readlines()


molden_files_folder = Path(__file__).parent / 'molden_files'

# Load example molecular structures
co = _read_file(molden_files_folder / 'co.inp')
o2 = _read_file(molden_files_folder / 'o2.inp')
co2 = _read_file(molden_files_folder / 'co2.inp')
h2o = _read_file(molden_files_folder / 'h2o.inp')
benzene = _read_file(molden_files_folder / 'benzene.inp')
prismane = _read_file(molden_files_folder / 'prismane.inp')
pyridine = _read_file(molden_files_folder / 'pyridine.inp')
furan = _read_file(molden_files_folder / 'furan.inp')
acrolein = _read_file(molden_files_folder / 'acrolein.inp')

#: Dictionary mapping example names to their molden file contents
all_examples = {
    'co': co,
    'o2': o2,
    'co2': co2,
    'h2o': h2o,
    'benzene': benzene,
    'prismane': prismane,
    'pyridine': pyridine,
    'furan': furan,
    'acrolein': acrolein,
}
