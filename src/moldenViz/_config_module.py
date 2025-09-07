import json
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import toml

default_configs_dir = Path(__file__).parent / 'default_configs'

custom_configs_dir = Path().home() / '.config/moldenViz'
custom_configs_dir.mkdir(parents=True, exist_ok=True)


@dataclass
class AtomType:
    """Represents the properties of an atom type for visualization.

    Parameters
    ----------
    name : str
        The name/symbol of the atom type (e.g., 'C', 'H', 'O').
    color : str
        The color to use for visualizing this atom type.
    radius : float
        The radius for displaying this atom type.
    max_num_bonds : int
        The maximum number of bonds this atom type can form.
    """

    name: str
    color: str
    radius: float
    max_num_bonds: int


class Config:
    """Configuration class to manage default and custom configurations."""

    def __init__(self) -> None:
        default_config = self.load_default_config()
        custom_config = self.load_custom_config()

        atoms_custom_config = custom_config.pop('Atom', {})

        self.config = self.merge_configs(default_config, custom_config)

        self.atom_types = self.load_atom_types(atoms_custom_config)

    @staticmethod
    def dict_to_namedspace(d: dict) -> SimpleNamespace:
        """Convert a dictionary to a SimpleNamespace for attribute-style access.

        Parameters
        ----------
        d : dict
            The dictionary to convert.

        Returns
        -------
        SimpleNamespace
            A SimpleNamespace object with attributes corresponding to the dictionary keys.
        """
        return SimpleNamespace(**{k: Config.dict_to_namedspace(v) if isinstance(v, dict) else v for k, v in d.items()})

    @staticmethod
    def merge_configs(default_config: dict, custom_config: dict) -> SimpleNamespace:
        """Merge multiple configuration dictionaries into a single SimpleNamespace.

        Parameters
        ----------
        default_config : dict
            The default configuration dictionary.
        custom_config : dict
            The custom configuration dictionary to merge with defaults.

        Returns
        -------
        SimpleNamespace
            A SimpleNamespace object with attributes corresponding to the merged configuration items.
        """
        return Config.dict_to_namedspace(Config.recursive_merge(default_config, custom_config))

    @staticmethod
    def recursive_merge(default: dict, custom: dict) -> dict:
        """Recursively merge two dictionaries.

        Parameters
        ----------
        default : dict
            The default dictionary.
        custom : dict
            The custom dictionary to merge with default.

        Returns
        -------
        dict
            The merged dictionary.
        """
        merged = default.copy()
        for k, v in custom.items():
            if isinstance(v, dict) and isinstance(default.get(k), dict):
                merged[k] = Config.recursive_merge(default[k], v)
            else:
                merged[k] = v
        return merged

    def __getattr__(self, item: str) -> Any:
        """Get an attribute from the configuration.

        Parameters
        ----------
        item : str
            The name of the configuration attribute to retrieve.

        Returns
        -------
        Any
            The value of the requested configuration item.

        Raises
        ------
        AttributeError
            If the requested attribute is not found in the configuration.
        """
        if not hasattr(self.config, item):
            raise AttributeError(f"No attribute '{item}' found in the configurations.")

        return getattr(self.config, item)

    @staticmethod
    def load_atom_types(atoms_custom_config: dict) -> dict[int, AtomType]:
        """Load default atom types from the JSON file and custom atom types from the custom config.

        Atom type based on atomic number
        Colors are based on CPK color scheme https://sciencenotes.org/molecule-atom-colors-cpk-colors/
        Radius come from default molden values
        Max number of bonds is set to author's best guess. Please, if you have better values, let me know!

        Parameters
        ----------
        atoms_custom_config : dict
            Custom configuration for atom types.

        Returns
        -------
        dict[int, AtomType]
            A dictionary mapping atomic numbers to AtomType objects.
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
        dict
            The default configuration dictionary.

        Raises
        ------
        FileNotFoundError
            If the default configuration file is not found.
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
        dict
            The custom configuration dictionary. Empty dict if file doesn't exist.
        """
        custom_config_path = custom_configs_dir / 'config.toml'
        if not custom_config_path.exists():
            return {}

        with custom_config_path.open('r') as f:
            return toml.load(f)
