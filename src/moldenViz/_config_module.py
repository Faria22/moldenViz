import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import toml

from ._plotting_objects import AtomType

default_configs_dir = Path(__file__).parent / 'default_configs'

custom_configs_dir = Path().home() / '.config/moldenViz'
custom_configs_dir.mkdir(parents=True, exist_ok=True)


class Config:
    """Configuration class to manage default and custom configurations."""

    def __init__(self) -> None:
        default_config = self.load_default_config()
        custom_config = self.load_custom_config()

        atoms_custom_config = custom_config.pop('Atom', {})

        self.config = self.dict_to_namedspace(default_config | custom_config)

        self.atom_types = self.load_atom_types(atoms_custom_config)

    @staticmethod
    def dict_to_namedspace(d: dict) -> SimpleNamespace:
        """Convert a dictionary to a SimpleNamespace for attribute-style access.

        Parameters
        ----------
            d: dict: The dictionary to convert.

        Returns
        -------
            SimpleNamespace: A SimpleNamespace object with attributes corresponding to the dictionary keys.
        """
        return SimpleNamespace(**{k: Config.dict_to_namedspace(v) if isinstance(v, dict) else v for k, v in d.items()})

    def __getattr__(self, item: str) -> Any:
        """Get an attribute from the configuration.

        Returns
        -------
            object: The value of the requested configuration item.
        """
        if not hasattr(self.config, item):
            raise AttributeError(f"No attribute '{item}' found in the configurations.")

        return self.config.item

    @staticmethod
    def load_atom_types(atoms_custom_config: dict) -> dict[int, AtomType]:
        """Load default atom types from the JSON file and custom atom types from the custom config.

        Atom type based on atomic number
        Colors are based on CPK color scheme https://sciencenotes.org/molecule-atom-colors-cpk-colors/
        Radius come from default molden values
        Max number of bonds is set to author's best guess. Please, if you have better values, let me know!

        Returns
        -------
            dict[int, AtomType]: A dictionary mapping atomic numbers to AtomType objects.

        """
        with (default_configs_dir / 'atom_types.json').open('r') as f:
            atom_types = json.load(f)

        atom_types = {int(k): AtomType(**v) for k, v in atom_types.items()}

        for atomic_number_str, atom_properties in atoms_custom_config.items():
            if atomic_number_str == 'show':
                continue  # Skip the 'show' key, which is not an atomic number

            try:
                atomic_number = int(atomic_number_str)
            except ValueError:
                raise ValueError(f'Invalid atomic number in custom configuration: {atomic_number_str}') from ValueError
            if not atom_types.get(atomic_number):
                raise ValueError('Invalid atomic number in custom configuration: %d', atomic_number)

            for prop, value in atom_properties.items():
                if hasattr(atom_types[atomic_number], prop):
                    setattr(atom_types[atomic_number], prop, value)
                else:
                    raise ValueError(f'Invalid property "{prop}" for atom in custom configuration.')

        return atom_types

    @staticmethod
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

    @staticmethod
    def load_custom_config() -> dict:
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
