"""Configuration module with validation using pydantic-style validation patterns."""

import json
import re
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Union

import toml

default_configs_dir = Path(__file__).parent / 'default_configs'

custom_configs_dir = Path().home() / '.config/moldenViz'
custom_configs_dir.mkdir(parents=True, exist_ok=True)


class ConfigurationError(Exception):
    """Custom exception for configuration validation errors."""


class ConfigValidator:
    """Configuration validator providing type checking and validation."""

    @staticmethod
    def validate_type(value: Any, expected_type: type, field_name: str) -> Any:
        """Validate that a value is of the expected type.
        
        Returns
        -------
        Any
            The validated value.
        """
        if not isinstance(value, expected_type):
            raise ConfigurationError(
                f"Configuration field '{field_name}' must be of type {expected_type.__name__}, "
                f"got {type(value).__name__}",
            )
        return value

    @staticmethod
    def validate_positive_number(value: float, field_name: str) -> Union[int, float]:
        """Validate that a number is positive.
        
        Returns
        -------
        Union[int, float]
            The validated positive number.
        """
        ConfigValidator.validate_type(value, (int, float), field_name)
        if value <= 0:
            raise ConfigurationError(f"Configuration field '{field_name}' must be positive, got {value}")
        return value

    @staticmethod
    def validate_range(value: float, min_val: float, max_val: float, field_name: str) -> Union[int, float]:
        """Validate that a value is within a specified range.
        
        Returns
        -------
        Union[int, float]
            The validated value within range.
        """
        ConfigValidator.validate_type(value, (int, float), field_name)
        if not (min_val <= value <= max_val):
            raise ConfigurationError(
                f"Configuration field '{field_name}' must be between {min_val} and {max_val}, got {value}",
            )
        return value

    @staticmethod
    def validate_color(value: str, field_name: str) -> str:
        """Validate that a string is a valid hex color.
        
        Returns
        -------
        str
            The validated hex color (normalized to uppercase).
        """
        ConfigValidator.validate_type(value, str, field_name)
        # Remove # if present and check if it's a valid hex color
        color = value.lstrip('#')
        if not re.match(r'^[0-9A-Fa-f]{6}$', color):
            raise ConfigurationError(
                f"Configuration field '{field_name}' must be a valid hex color (6 hex digits), got '{value}'",
            )
        return color.upper()  # Normalize to uppercase

    @staticmethod
    def validate_choices(value: Any, choices: List[Any], field_name: str) -> Any:
        """Validate that a value is one of the allowed choices.
        
        Returns
        -------
        Any
            The validated choice value.
        """
        if value not in choices:
            raise ConfigurationError(
                f"Configuration field '{field_name}' must be one of {choices}, got '{value}'",
            )
        return value

    @staticmethod
    def validate_opacity(value: float, field_name: str) -> Union[int, float]:
        """Validate opacity value (0.0 to 1.0).
        
        Returns
        -------
        Union[int, float]
            The validated opacity value.
        """
        return ConfigValidator.validate_range(value, 0.0, 1.0, field_name)


@dataclass
class AtomType:
    """Represents the properties of an atom type for visualization.

    Parameters
    ----------
    name : str
        The name/symbol of the atom type (e.g., 'C', 'H', 'O').
    color : str
        The color to use for visualizing this atom type (hex color code).
    radius : float
        The radius for displaying this atom type (must be positive).
    max_num_bonds : int
        The maximum number of bonds this atom type can form (0-10).
    """

    name: str
    color: str
    radius: float
    max_num_bonds: int

    def __post_init__(self) -> None:
        """Validate atom type properties after initialization."""
        try:
            # Validate name
            self.name = ConfigValidator.validate_type(self.name, str, 'name').strip()
            if not self.name:
                raise ConfigurationError('Atom name cannot be empty')

            # Validate color
            self.color = ConfigValidator.validate_color(self.color, 'color')

            # Validate radius
            self.radius = ConfigValidator.validate_positive_number(self.radius, 'radius')

            # Validate max_num_bonds
            self.max_num_bonds = ConfigValidator.validate_type(self.max_num_bonds, int, 'max_num_bonds')
            self.max_num_bonds = ConfigValidator.validate_range(self.max_num_bonds, 0, 10, 'max_num_bonds')

        except (ValueError, TypeError) as e:
            raise ConfigurationError(f'Invalid atom type configuration: {e}') from e

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'AtomType':
        """Create AtomType from dictionary with validation.
        
        Returns
        -------
        AtomType
            The validated AtomType instance.
        """
        required_fields = ['name', 'color', 'radius', 'max_num_bonds']
        missing_fields = [field for field in required_fields if field not in data]
        if missing_fields:
            raise ConfigurationError(f'Missing required atom type fields: {missing_fields}')

        return cls(
            name=data['name'],
            color=data['color'],
            radius=data['radius'],
            max_num_bonds=data['max_num_bonds'],
        )


class Config:
    """Configuration class to manage default and custom configurations with validation."""

    def __init__(self) -> None:
        try:
            default_config = self.load_default_config()
            custom_config = self.load_custom_config()

            # Validate configurations before processing
            default_config = self._validate_config(default_config)
            custom_config = self._validate_config(custom_config)

            atoms_custom_config = custom_config.pop('Atom', {})

            self.config = self.merge_configs(default_config, custom_config)

            self.atom_types = self.load_atom_types(atoms_custom_config)
        except Exception as e:
            raise ConfigurationError(f'Failed to initialize configuration: {e}') from e

    def _validate_config(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Validate a configuration dictionary.

        Returns
        -------
        Dict[str, Any]
            Validated configuration dictionary.
        """
        if not config:  # Empty config is valid (for custom configs)
            return config

        validated_config = {}

        # Validate top-level boolean settings
        if 'smooth_shading' in config:
            validated_config['smooth_shading'] = ConfigValidator.validate_type(
                config['smooth_shading'], bool, 'smooth_shading',
            )

        # Validate grid configuration
        if 'grid' in config:
            validated_config['grid'] = self._validate_grid_config(config['grid'])

        # Validate MO configuration
        if 'MO' in config:
            validated_config['MO'] = self._validate_mo_config(config['MO'])

        # Validate molecule configuration
        if 'molecule' in config:
            validated_config['molecule'] = self._validate_molecule_config(config['molecule'])

        # Pass through Atom configuration for separate processing
        if 'Atom' in config:
            validated_config['Atom'] = config['Atom']

        return validated_config

    def _validate_grid_config(self, grid_config: Dict[str, Any]) -> Dict[str, Any]:
        """Validate grid configuration section.
        
        Returns
        -------
        Dict[str, Any]
            Validated grid configuration.
        """
        validated = {}

        if 'min_radius' in grid_config:
            validated['min_radius'] = ConfigValidator.validate_positive_number(
                grid_config['min_radius'], 'grid.min_radius',
            )

        if 'max_radius_multiplier' in grid_config:
            validated['max_radius_multiplier'] = ConfigValidator.validate_positive_number(
                grid_config['max_radius_multiplier'], 'grid.max_radius_multiplier',
            )

        if 'spherical' in grid_config:
            validated['spherical'] = self._validate_spherical_config(grid_config['spherical'])

        if 'cartesian' in grid_config:
            validated['cartesian'] = self._validate_cartesian_config(grid_config['cartesian'])

        return validated

    def _validate_spherical_config(self, spherical_config: Dict[str, Any]) -> Dict[str, Any]:
        """Validate spherical grid configuration.
        
        Returns
        -------
        Dict[str, Any]
            Validated spherical configuration.
        """
        validated = {}

        grid_point_fields = ['num_r_points', 'num_theta_points', 'num_phi_points']
        for field in grid_point_fields:
            if field in spherical_config:
                validated[field] = ConfigValidator.validate_range(
                    spherical_config[field], 10, 1000, f'grid.spherical.{field}',
                )

        return validated

    def _validate_cartesian_config(self, cartesian_config: Dict[str, Any]) -> Dict[str, Any]:
        """Validate cartesian grid configuration.
        
        Returns
        -------
        Dict[str, Any]
            Validated cartesian configuration.
        """
        validated = {}

        grid_point_fields = ['num_x_points', 'num_y_points', 'num_z_points']
        for field in grid_point_fields:
            if field in cartesian_config:
                validated[field] = ConfigValidator.validate_range(
                    cartesian_config[field], 10, 1000, f'grid.cartesian.{field}',
                )

        return validated

    def _validate_mo_config(self, mo_config: Dict[str, Any]) -> Dict[str, Any]:
        """Validate molecular orbital configuration.
        
        Returns
        -------
        Dict[str, Any]
            Validated MO configuration.
        """
        validated = {}

        if 'contour' in mo_config:
            validated['contour'] = ConfigValidator.validate_positive_number(
                mo_config['contour'], 'MO.contour',
            )

        if 'opacity' in mo_config:
            validated['opacity'] = ConfigValidator.validate_opacity(
                mo_config['opacity'], 'MO.opacity',
            )

        return validated

    def _validate_molecule_config(self, molecule_config: Dict[str, Any]) -> Dict[str, Any]:
        """Validate molecule configuration.
        
        Returns
        -------
        Dict[str, Any]
            Validated molecule configuration.
        """
        validated = {}

        if 'opacity' in molecule_config:
            validated['opacity'] = ConfigValidator.validate_opacity(
                molecule_config['opacity'], 'molecule.opacity',
            )

        if 'atom' in molecule_config:
            validated['atom'] = self._validate_atom_config(molecule_config['atom'])

        if 'bond' in molecule_config:
            validated['bond'] = self._validate_bond_config(molecule_config['bond'])

        return validated

    def _validate_atom_config(self, atom_config: Dict[str, Any]) -> Dict[str, Any]:
        """Validate atom display configuration.
        
        Returns
        -------
        Dict[str, Any]
            Validated atom configuration.
        """
        validated = {}

        if 'show' in atom_config:
            validated['show'] = ConfigValidator.validate_type(
                atom_config['show'], bool, 'molecule.atom.show',
            )

        return validated

    def _validate_bond_config(self, bond_config: Dict[str, Any]) -> Dict[str, Any]:
        """Validate bond configuration.
        
        Returns
        -------
        Dict[str, Any]
            Validated bond configuration.
        """
        validated = {}

        if 'show' in bond_config:
            validated['show'] = ConfigValidator.validate_type(
                bond_config['show'], bool, 'molecule.bond.show',
            )

        if 'max_length' in bond_config:
            validated['max_length'] = ConfigValidator.validate_positive_number(
                bond_config['max_length'], 'molecule.bond.max_length',
            )

        if 'color_type' in bond_config:
            validated['color_type'] = ConfigValidator.validate_choices(
                bond_config['color_type'], ['uniform', 'element'], 'molecule.bond.color_type',
            )

        if 'color' in bond_config:
            validated['color'] = ConfigValidator.validate_type(
                bond_config['color'], str, 'molecule.bond.color',
            )

        if 'radius' in bond_config:
            validated['radius'] = ConfigValidator.validate_positive_number(
                bond_config['radius'], 'molecule.bond.radius',
            )

        return validated

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
        try:
            with (default_configs_dir / 'atom_types.json').open('r') as f:
                atom_types_data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            raise ConfigurationError(f'Failed to load atom types from JSON: {e}') from e

        # Create AtomType objects with validation
        atom_types = {}
        for atomic_number_str, atom_data in atom_types_data.items():
            atomic_number = int(atomic_number_str)
            try:
                atom_types[atomic_number] = AtomType.from_dict(atom_data)
            except (ValueError, ConfigurationError) as e:
                raise ConfigurationError(
                    f'Invalid atom type data for atomic number {atomic_number_str}: {e}',
                ) from e

        # Apply custom configurations with validation
        for atomic_number_str, atom_properties in atoms_custom_config.items():
            if atomic_number_str == 'show':
                continue  # Skip the 'show' key, which is not an atomic number

            try:
                atomic_number = int(atomic_number_str)
            except ValueError as e:
                raise ConfigurationError(
                    f'Invalid atomic number in custom configuration: {atomic_number_str}',
                ) from e

            if atomic_number not in atom_types:
                raise ConfigurationError(
                    f'Invalid atomic number in custom configuration: {atomic_number}. '
                    f'Atomic number must be between 1 and {max(atom_types.keys())}',
                )

            # Apply custom properties with validation
            for prop, value in atom_properties.items():
                if not hasattr(atom_types[atomic_number], prop):
                    raise ConfigurationError(
                        f'Invalid property "{prop}" for atom in custom configuration. '
                        f'Valid properties are: name, color, radius, max_num_bonds',
                    )

                # Create a new AtomType with the updated property to trigger validation
                current_atom = atom_types[atomic_number]
                atom_dict = {
                    'name': current_atom.name,
                    'color': current_atom.color,
                    'radius': current_atom.radius,
                    'max_num_bonds': current_atom.max_num_bonds,
                }
                atom_dict[prop] = value

                try:
                    atom_types[atomic_number] = AtomType.from_dict(atom_dict)
                except ConfigurationError as e:
                    raise ConfigurationError(
                        f'Invalid value for property "{prop}" of atom {atomic_number}: {e}',
                    ) from e

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
        ConfigurationError
            If the default configuration file is not found or contains invalid data.
        """
        default_config_path = default_configs_dir / 'config.toml'
        if not default_config_path.exists():
            raise ConfigurationError(f'Default configuration file not found at {default_config_path}')

        try:
            with default_config_path.open('r') as f:
                return toml.load(f)
        except (toml.TomlDecodeError, OSError) as e:
            raise ConfigurationError(f'Failed to load default configuration: {e}') from e

    @staticmethod
    def load_custom_config() -> dict:
        """Load custom configuration from the TOML file.

        Returns
        -------
        dict
            The custom configuration dictionary. Empty dict if file doesn't exist.

        Raises
        ------
        ConfigurationError
            If the custom configuration file exists but contains invalid data.
        """
        custom_config_path = custom_configs_dir / 'config.toml'
        if not custom_config_path.exists():
            return {}

        try:
            with custom_config_path.open('r') as f:
                return toml.load(f)
        except (toml.TomlDecodeError, OSError) as e:
            raise ConfigurationError(f'Failed to load custom configuration from {custom_config_path}: {e}') from e
