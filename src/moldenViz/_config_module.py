import json
from pathlib import Path

import toml

from ._plotting_objects import AtomType

default_configs_dir = Path(__file__).parent / 'default_configs'

custom_configs_dir = Path().home() / '.config/moldenViz'
custom_configs_dir.mkdir(parents=True, exist_ok=True)


def load_atom_types() -> dict[int, AtomType]:
    """Load atom types from the JSON file.

    Atom type based on charge
    Colors are based on CPK color scheme
    https://sciencenotes.org/molecule-atom-colors-cpk-colors/

    Returns
    -------
        dict[int, AtomType]: A dictionary mapping atomic numbers to AtomType objects.

    """
    with (default_configs_dir / 'atom_types.json').open('r') as f:
        atom_types = json.load(f)

    atom_types = {int(k): AtomType(**v) for k, v in atom_types.items()}

    # TODO: Add custom atom types from the custom config if they exist
    return atom_types


def load_default_config() -> dict:
    """Load default configuration from the TOML file.

    Returns
    -------
        dict: The custom configuration dictionary.

    """
    default_config_path = default_configs_dir / 'config.toml'
    if not default_config_path.exists():
        raise FileNotFoundError(f'Default configuration file not found at {default_config_path}. ')

    with default_config_path.open('r') as f:
        return toml.load(f)


def load_custom_congig() -> dict:
    """Load custom configuration from the TOML file.

    Returns
    -------
        dict: The custom configuration dictionary.

    """
    custom_config_path = custom_configs_dir / 'config.toml'
    if not custom_config_path.exists():
        return {}

    with custom_config_path.open('r') as f:
        return toml.load(f)
