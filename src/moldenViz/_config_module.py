# type: ignore[reportArgumentType]
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Literal

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import toml
from pydantic import BaseModel, Field, field_validator

# Global config directory paths
DEFAULT_CONFIGS_DIR = Path(__file__).parent / 'default_configs'
CUSTOM_CONFIGS_DIR = Path().home() / '.config/moldenViz'
CUSTOM_CONFIGS_DIR.mkdir(parents=True, exist_ok=True)

# Global config file paths
DEFAULT_CONFIG_PATH = DEFAULT_CONFIGS_DIR / 'config.toml'
CUSTOM_CONFIG_PATH = CUSTOM_CONFIGS_DIR / 'config.toml'
ATOM_TYPES_PATH = DEFAULT_CONFIGS_DIR / 'atom_types.json'

# Maintain backwards compatibility
default_configs_dir = DEFAULT_CONFIGS_DIR
custom_configs_dir = CUSTOM_CONFIGS_DIR


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

    num_r_points: int = Field(100, gt=0, description='Number of radial points')
    num_theta_points: int = Field(60, gt=0, description='Number of theta points')
    num_phi_points: int = Field(120, gt=0, description='Number of phi points')


class CartesianGridConfig(BaseModel):
    """Configuration for cartesian grid parameters."""

    num_x_points: int = Field(100, gt=0, description='Number of x points')
    num_y_points: int = Field(100, gt=0, description='Number of y points')
    num_z_points: int = Field(100, gt=0, description='Number of z points')


class GridConfig(BaseModel):
    """Configuration for grid generation."""

    min_radius: int = Field(5, gt=0, description='Minimum radius')
    max_radius_multiplier: int = Field(2, gt=0, description='Max radius multiplier')
    default_type: Literal['spherical', 'cartesian'] = Field('spherical', description='Default grid type')
    spherical: SphericalGridConfig = Field(default_factory=SphericalGridConfig)
    cartesian: CartesianGridConfig = Field(default_factory=CartesianGridConfig)


class MOConfig(BaseModel):
    """Configuration for molecular orbital visualization."""

    contour: float = Field(0.1, gt=0, description='Contour value')
    opacity: float = Field(1.0, ge=0, le=1, description='Opacity (0-1)')
    color_scheme: str = Field(
        'bwr',
        description='Colormap for MO visualization (diverging colormaps recommended)',
    )
    custom_colors: list[str] | None = Field(
        None,
        min_length=2,
        max_length=2,
        description='Custom two colors for MO visualization [negative, positive]',
    )

    @field_validator('color_scheme')
    @classmethod
    def validate_color_scheme(cls, v: str) -> str:
        """Validate color scheme is a valid matplotlib colormap.

        Parameters
        ----------
        v : str
            The color scheme name to validate.

        Returns
        -------
        str
            The validated color scheme string.

        Raises
        ------
        ValueError
            If the color scheme is not a valid matplotlib colormap.
        """
        try:
            plt.get_cmap(v)
        except ValueError as e:
            raise ValueError(f'Color scheme must be a valid matplotlib colormap. Got: {v}') from e
        else:
            return v

    @field_validator('custom_colors')
    @classmethod
    def validate_custom_colors(cls, v: list[str] | None) -> list[str] | None:
        """Validate custom colors using matplotlib.colors.is_color_like.

        Parameters
        ----------
        v : list[str] | None
            The list of colors to validate.

        Returns
        -------
        list[str] | None
            The validated list of colors or None.

        Raises
        ------
        ValueError
            If any color in the list is not a valid matplotlib color.
        """
        if v is None:
            return v

        for color in v:
            if not mcolors.is_color_like(color):
                raise ValueError(f'Custom color must be a valid matplotlib color. Got: {color}')
        return v


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
    background_color: str = Field('white', description='Background color for 3D visualization')
    grid: GridConfig = Field(default_factory=GridConfig)
    mo: MOConfig = Field(default_factory=MOConfig)
    molecule: MoleculeConfig = Field(default_factory=MoleculeConfig)

    @field_validator('background_color')
    @classmethod
    def validate_background_color(cls, v: str) -> str:
        """Validate background color using matplotlib.colors.is_color_like.

        Parameters
        ----------
        v : str
            The background color to validate.

        Returns
        -------
        str
            The validated background color string.

        Raises
        ------
        ValueError
            If the background color is not a valid matplotlib color.
        """
        if mcolors.is_color_like(v):
            return v
        raise ValueError(f'Background color must be a valid matplotlib color. Got: {v}')

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
        with ATOM_TYPES_PATH.open('r') as f:
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
        if not DEFAULT_CONFIG_PATH.exists():
            raise FileNotFoundError(f'Default configuration file not found at {DEFAULT_CONFIG_PATH}. ')

        with DEFAULT_CONFIG_PATH.open('r') as f:
            return toml.load(f)

    @staticmethod
    def load_custom_config() -> dict:
        """Load custom configuration from the TOML file.

        Returns
        -------
        dict
            The custom configuration dictionary. Empty dict if file doesn't exist.
        """
        if not CUSTOM_CONFIG_PATH.exists():
            return {}

        with CUSTOM_CONFIG_PATH.open('r') as f:
            return toml.load(f)

    def save_current_config(self) -> None:
        """Save the current configuration to the custom config file.

        This method writes the current configuration values to ~/.config/moldenViz/config.toml,
        preserving the TOML structure.
        """
        # Build the configuration dict from the current values
        config_dict = {
            'smooth_shading': self.config.smooth_shading,
            'background_color': self.config.background_color,
            'grid': {
                'min_radius': self.config.grid.min_radius,
                'max_radius_multiplier': self.config.grid.max_radius_multiplier,
                'spherical': {
                    'num_r_points': self.config.grid.spherical.num_r_points,
                    'num_theta_points': self.config.grid.spherical.num_theta_points,
                    'num_phi_points': self.config.grid.spherical.num_phi_points,
                },
                'cartesian': {
                    'num_x_points': self.config.grid.cartesian.num_x_points,
                    'num_y_points': self.config.grid.cartesian.num_y_points,
                    'num_z_points': self.config.grid.cartesian.num_z_points,
                },
            },
            'MO': {
                'contour': self.config.mo.contour,
                'opacity': self.config.mo.opacity,
                'color_scheme': self.config.mo.color_scheme,
            },
            'molecule': {
                'opacity': self.config.molecule.opacity,
                'atom': {
                    'show': self.config.molecule.atom.show,
                },
                'bond': {
                    'show': self.config.molecule.bond.show,
                    'max_length': self.config.molecule.bond.max_length,
                    'color_type': self.config.molecule.bond.color_type,
                    'color': self.config.molecule.bond.color,
                    'radius': self.config.molecule.bond.radius,
                },
            },
        }

        # Add custom_colors if set
        if self.config.mo.custom_colors is not None:
            config_dict['MO']['custom_colors'] = self.config.mo.custom_colors

        # Write to file
        with CUSTOM_CONFIG_PATH.open('w') as f:
            toml.dump(config_dict, f)
