"""Tests for the configuration module validation functionality."""

try:
    import pytest
    HAS_PYTEST = True
except ImportError:
    HAS_PYTEST = False
    
    # Simple pytest.raises replacement
    class MockRaises:
        def __init__(self, exception_type, match=None):
            self.exception_type = exception_type
            self.match = match
        
        def __enter__(self):
            return self
        
        def __exit__(self, exc_type, exc_val, exc_tb):
            if exc_type is None:
                raise AssertionError(f'Expected {self.exception_type.__name__} but no exception was raised')
            if not issubclass(exc_type, self.exception_type):
                return False  # Let the exception propagate
            if self.match and self.match not in str(exc_val):
                raise AssertionError(f'Expected exception message to contain "{self.match}" but got "{exc_val}"')
            return True  # Suppress the expected exception
    
    class MockPytest:
        raises = MockRaises
    
    pytest = MockPytest()

from moldenViz._config_module import Config, AtomType, ConfigurationError, ConfigValidator


class TestConfigValidator:
    """Test the ConfigValidator class methods."""
    
    def test_validate_type_valid(self):
        """Test type validation with valid inputs."""
        assert ConfigValidator.validate_type("test", str, "field") == "test"
        assert ConfigValidator.validate_type(42, int, "field") == 42
        assert ConfigValidator.validate_type(3.14, float, "field") == 3.14
    
    def test_validate_type_invalid(self):
        """Test type validation with invalid inputs."""
        with pytest.raises(ConfigurationError, match="must be of type str"):
            ConfigValidator.validate_type(42, str, "field")
    
    def test_validate_positive_number_valid(self):
        """Test positive number validation with valid inputs."""
        assert ConfigValidator.validate_positive_number(5, "field") == 5
        assert ConfigValidator.validate_positive_number(3.14, "field") == 3.14
    
    def test_validate_positive_number_invalid(self):
        """Test positive number validation with invalid inputs."""
        with pytest.raises(ConfigurationError, match="must be positive"):
            ConfigValidator.validate_positive_number(-5, "field")
        with pytest.raises(ConfigurationError, match="must be positive"):
            ConfigValidator.validate_positive_number(0, "field")
    
    def test_validate_range_valid(self):
        """Test range validation with valid inputs."""
        assert ConfigValidator.validate_range(5, 0, 10, "field") == 5
        assert ConfigValidator.validate_range(0, 0, 10, "field") == 0
        assert ConfigValidator.validate_range(10, 0, 10, "field") == 10
    
    def test_validate_range_invalid(self):
        """Test range validation with invalid inputs."""
        with pytest.raises(ConfigurationError, match="must be between"):
            ConfigValidator.validate_range(-1, 0, 10, "field")
        with pytest.raises(ConfigurationError, match="must be between"):
            ConfigValidator.validate_range(11, 0, 10, "field")
    
    def test_validate_color_valid(self):
        """Test color validation with valid inputs."""
        assert ConfigValidator.validate_color("FF0000", "field") == "FF0000"
        assert ConfigValidator.validate_color("00ff00", "field") == "00FF00"  # Normalized to uppercase
        assert ConfigValidator.validate_color("#123ABC", "field") == "123ABC"  # Hash removed
    
    def test_validate_color_invalid(self):
        """Test color validation with invalid inputs."""
        with pytest.raises(ConfigurationError, match="valid hex color"):
            ConfigValidator.validate_color("invalid", "field")
        with pytest.raises(ConfigurationError, match="valid hex color"):
            ConfigValidator.validate_color("FF00", "field")  # Too short
        with pytest.raises(ConfigurationError, match="valid hex color"):
            ConfigValidator.validate_color("GG0000", "field")  # Invalid hex
    
    def test_validate_choices_valid(self):
        """Test choice validation with valid inputs."""
        assert ConfigValidator.validate_choices("a", ["a", "b", "c"], "field") == "a"
        assert ConfigValidator.validate_choices(1, [1, 2, 3], "field") == 1
    
    def test_validate_choices_invalid(self):
        """Test choice validation with invalid inputs."""
        with pytest.raises(ConfigurationError, match="must be one of"):
            ConfigValidator.validate_choices("d", ["a", "b", "c"], "field")
    
    def test_validate_opacity_valid(self):
        """Test opacity validation with valid inputs."""
        assert ConfigValidator.validate_opacity(0.0, "field") == 0.0
        assert ConfigValidator.validate_opacity(0.5, "field") == 0.5
        assert ConfigValidator.validate_opacity(1.0, "field") == 1.0
    
    def test_validate_opacity_invalid(self):
        """Test opacity validation with invalid inputs."""
        with pytest.raises(ConfigurationError, match="must be between"):
            ConfigValidator.validate_opacity(-0.1, "field")
        with pytest.raises(ConfigurationError, match="must be between"):
            ConfigValidator.validate_opacity(1.1, "field")


class TestAtomType:
    """Test the AtomType class validation."""
    
    def test_valid_atom_type(self):
        """Test creation of valid atom types."""
        atom = AtomType(name="C", color="909090", radius=0.4, max_num_bonds=4)
        assert atom.name == "C"
        assert atom.color == "909090"
        assert atom.radius == 0.4
        assert atom.max_num_bonds == 4
    
    def test_from_dict_valid(self):
        """Test creation from valid dictionary."""
        data = {"name": "H", "color": "FFFFFF", "radius": 0.2, "max_num_bonds": 1}
        atom = AtomType.from_dict(data)
        assert atom.name == "H"
        assert atom.color == "FFFFFF"
    
    def test_from_dict_missing_fields(self):
        """Test creation from dictionary with missing fields."""
        data = {"name": "H", "color": "FFFFFF", "radius": 0.2}  # Missing max_num_bonds
        with pytest.raises(ConfigurationError, match="Missing required"):
            AtomType.from_dict(data)
    
    def test_invalid_name(self):
        """Test atom type with invalid name."""
        with pytest.raises(ConfigurationError, match="empty"):
            AtomType(name="", color="FF0000", radius=0.4, max_num_bonds=4)
    
    def test_invalid_color(self):
        """Test atom type with invalid color."""
        with pytest.raises(ConfigurationError, match="hex color"):
            AtomType(name="C", color="invalid", radius=0.4, max_num_bonds=4)
    
    def test_invalid_radius(self):
        """Test atom type with invalid radius."""
        with pytest.raises(ConfigurationError, match="positive"):
            AtomType(name="C", color="FF0000", radius=-0.1, max_num_bonds=4)
    
    def test_invalid_max_bonds(self):
        """Test atom type with invalid max bonds."""
        with pytest.raises(ConfigurationError, match="between 0 and 10"):
            AtomType(name="C", color="FF0000", radius=0.4, max_num_bonds=15)


class TestConfig:
    """Test the Config class validation."""
    
    def test_config_loading(self):
        """Test basic configuration loading."""
        config = Config()
        assert hasattr(config, 'grid')
        assert hasattr(config, 'MO')
        assert hasattr(config, 'molecule')
        assert len(config.atom_types) > 0
    
    def test_grid_config_validation(self):
        """Test grid configuration validation."""
        config = Config()
        
        # Valid grid config
        valid_grid = {"min_radius": 10, "max_radius_multiplier": 3}
        result = config._validate_grid_config(valid_grid)
        assert result["min_radius"] == 10
        assert result["max_radius_multiplier"] == 3
        
        # Invalid grid config
        with pytest.raises(ConfigurationError):
            config._validate_grid_config({"min_radius": -5})
    
    def test_mo_config_validation(self):
        """Test MO configuration validation."""
        config = Config()
        
        # Valid MO config
        valid_mo = {"contour": 0.1, "opacity": 0.8}
        result = config._validate_mo_config(valid_mo)
        assert result["contour"] == 0.1
        assert result["opacity"] == 0.8
        
        # Invalid MO config
        with pytest.raises(ConfigurationError):
            config._validate_mo_config({"opacity": 1.5})  # opacity > 1
    
    def test_molecule_config_validation(self):
        """Test molecule configuration validation."""
        config = Config()
        
        # Valid molecule config
        valid_molecule = {
            "opacity": 0.9,
            "atom": {"show": True},
            "bond": {"show": True, "max_length": 2.0, "color_type": "uniform"}
        }
        result = config._validate_molecule_config(valid_molecule)
        assert result["opacity"] == 0.9
        assert result["atom"]["show"] is True
        assert result["bond"]["color_type"] == "uniform"
        
        # Invalid molecule config
        with pytest.raises(ConfigurationError):
            config._validate_molecule_config({
                "bond": {"color_type": "invalid_type"}
            })


if __name__ == "__main__":
    # Simple test runner for when pytest is not available
    import sys
    
    def run_tests():
        """Simple test runner."""
        test_classes = [TestConfigValidator, TestAtomType, TestConfig]
        total_tests = 0
        passed_tests = 0
        
        for test_class in test_classes:
            instance = test_class()
            methods = [method for method in dir(instance) if method.startswith('test_')]
            
            for method_name in methods:
                total_tests += 1
                try:
                    method = getattr(instance, method_name)
                    method()
                    print(f"✓ {test_class.__name__}.{method_name}")
                    passed_tests += 1
                except Exception as e:
                    print(f"✗ {test_class.__name__}.{method_name}: {e}")
        
        print(f"\nTests: {passed_tests}/{total_tests} passed")
        return passed_tests == total_tests
    
    success = run_tests()
    sys.exit(0 if success else 1)