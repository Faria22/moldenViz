"""Get example files from folder and make them available to the package."""

from pathlib import Path


def _read_file(path: Path) -> list[str]:
    with path.open('r') as f:
        return f.readlines()


molden_files_folder = Path(__file__).parent / 'molden_files'

co = _read_file(molden_files_folder / 'co.inp')
o2 = _read_file(molden_files_folder / 'o2.inp')
co2 = _read_file(molden_files_folder / 'co2.inp')
h2o = _read_file(molden_files_folder / 'h2o.inp')
benzene = _read_file(molden_files_folder / 'benzene.inp')
prismane = _read_file(molden_files_folder / 'prismane.inp')
pyridine = _read_file(molden_files_folder / 'pyridine.inp')
furan = _read_file(molden_files_folder / 'furan.inp')
acrolein = _read_file(molden_files_folder / 'acrolein.inp')

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
