# type: ignore[reportArgumentType]
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Literal

import matplotlib.colors as mcolors
import toml
from pydantic import BaseModel, Field, field_validator

default_configs_dir = Path(__file__).parent / 'default_configs'

custom_configs_dir = Path().home() / '.config/moldenViz'
custom_configs_dir.mkdir(parents=True, exist_ok=True)


class AtomType(BaseModel):
    """Represents the properties of an atom type for visualization.

    Parameters
    ----------
    name : str
        The name/symbol of the atom type (e.g., 'C', 'H', 'O').
    color : str
        The color to use for visualizing this atom type (hex color without #).
    radius : float
        The radius for displaying this atom type (must be positive).
    max_num_bonds : int
        The maximum number of bonds this atom type can form (non-negative).
    """

    name: str = Field(..., min_length=1, max_length=3, description='Atom symbol')
    color: str = Field(..., pattern=r'^[0-9A-Fa-f]{6}$', description='Hex color code without #')
    radius: float = Field(..., gt=0, description='Atom radius (must be positive)')
    max_num_bonds: int = Field(..., ge=0, description='Maximum number of bonds')


class SphericalGridConfig(BaseModel):
    """Configuration for spherical grid parameters."""

    num_r_points: int = Field(100, gt=0, description='Number of radial points (1-1000)')
    num_theta_points: int = Field(60, gt=0, description='Number of theta points (1-1000)')
    num_phi_points: int = Field(120, gt=0, description='Number of phi points (1-1000)')


class CartesianGridConfig(BaseModel):
    """Configuration for cartesian grid parameters."""

    num_x_points: int = Field(100, gt=0, description='Number of x points (1-1000)')
    num_y_points: int = Field(100, gt=0, description='Number of y points (1-1000)')
    num_z_points: int = Field(100, gt=0, description='Number of z points (1-1000)')


class GridConfig(BaseModel):
    """Configuration for grid generation."""

    min_radius: int = Field(5, gt=0, description='Minimum radius (1-100)')
    max_radius_multiplier: int = Field(2, gt=0, description='Max radius multiplier (1-10)')
    spherical: SphericalGridConfig = Field(default_factory=SphericalGridConfig)
    cartesian: CartesianGridConfig = Field(default_factory=CartesianGridConfig)


class MOConfig(BaseModel):
    """Configuration for molecular orbital visualization."""

    contour: float = Field(0.1, gt=0, description='Contour value')
    opacity: float = Field(1.0, ge=0, le=1, description='Opacity (0-1)')


class AtomDisplayConfig(BaseModel):
    """Configuration for atom display."""

    show: bool = Field(True, description='Whether to show atoms')


class BondConfig(BaseModel):
    """Configuration for bond display."""

    show: bool = Field(True, description='Whether to show bonds')
    max_length: float = Field(4.0, gt=0, description='Maximum bond length')
    color_type: Literal['uniform', 'split'] = Field('uniform', description='Bond color type')
    color: str = Field('grey', description='Bond color (hex code or common name)')
    radius: float = Field(0.15, gt=0, description='Bond radius')

    @field_validator('color')
    @classmethod
    def validate_color(cls, v: str) -> str:
        """Validate color using matplotlib.colors.is_color_like.

        Returns
        -------
        str
            The validated color string.
        """
        if mcolors.is_color_like(v):
            return v
        raise ValueError(f'Color must be a valid matplotlib color. Got: {v}')


class MoleculeConfig(BaseModel):
    """Configuration for molecule display."""

    opacity: float = Field(1.0, ge=0, le=1, description='Molecule opacity (0-1)')
    atom: AtomDisplayConfig = Field(default_factory=AtomDisplayConfig)
    bond: BondConfig = Field(default_factory=BondConfig)


class MainConfig(BaseModel):
    """Main configuration model for moldenViz."""

    smooth_shading: bool = Field(True, description='Enable smooth shading')
    grid: GridConfig = Field(default_factory=GridConfig)
    mo: MOConfig = Field(default_factory=MOConfig)
    molecule: MoleculeConfig = Field(default_factory=MoleculeConfig)

    class ConfigDict:
        populate_by_name = True


class Config:
    """Configuration class to manage default and custom configurations."""

    def __init__(self) -> None:
        default_config = self.load_default_config()
        custom_config = self.load_custom_config()

        atoms_custom_config = custom_config.pop('Atom', {})

        # Validate and merge configuration using pydantic
        merged_config_dict = self.recursive_merge(default_config, custom_config)

        # Validate the merged configuration with pydantic
        try:
            self._pydantic_config = MainConfig(**merged_config_dict)
        except Exception as e:
            raise ValueError(f'Invalid configuration: {e}') from e

        # Convert to SimpleNamespace for backward compatibility
        self.config = self.dict_to_namedspace(self._pydantic_config.model_dump(by_alias=True))

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
            atom_types_data = json.load(f)

        # Validate and create AtomType objects using pydantic
        atom_types = {}
        for k, v in atom_types_data.items():
            try:
                atom_types[int(k)] = AtomType(**v)
            except Exception as e:  # noqa: PERF203
                raise ValueError(f'Invalid atom type data for atomic number {k}: {e}') from e

        for atomic_number_str, atom_properties in atoms_custom_config.items():
            if atomic_number_str == 'show':
                continue  # Skip the 'show' key, which is not an atomic number

            try:
                atomic_number = int(atomic_number_str)
            except ValueError:
                raise ValueError(f'Invalid atomic number in custom configuration: {atomic_number_str}') from ValueError

            if atomic_number not in atom_types:
                raise ValueError(f'Invalid atomic number in custom configuration: {atomic_number}')

            # Update the atom type with custom properties
            current_atom = atom_types[atomic_number]
            updated_data = current_atom.dict()

            for prop, value in atom_properties.items():
                if prop in updated_data:
                    updated_data[prop] = value
                else:
                    raise ValueError(f'Invalid property "{prop}" for atom in custom configuration.')

            # Validate the updated atom type
            try:
                atom_types[atomic_number] = AtomType(**updated_data)
            except Exception as e:
                raise ValueError(f'Invalid custom atom type for atomic number {atomic_number}: {e}') from e

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
