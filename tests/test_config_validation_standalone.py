#!/usr/bin/env python3
"""Standalone tests for configuration validation without pytest dependency."""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

# Import directly from the module file to avoid tkinter issues
import json
import re
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Union
import toml

# Set __file__ for the module
module_file = Path(__file__).parent.parent / 'src' / 'moldenViz' / '_config_module.py'

# Execute the module content to get the classes
with open(module_file, 'r') as f:
    module_content = f.read()

# Replace __file__ in the module content
module_content = module_content.replace(
    "default_configs_dir = Path(__file__).parent / 'default_configs'",
    f"default_configs_dir = Path('{module_file.parent}') / 'default_configs'"
)

exec(module_content)


def run_validation_tests():
    """Run comprehensive validation tests."""
    print("Running configuration validation tests...")
    
    total_tests = 0
    passed_tests = 0
    
    def test(description, test_func):
        nonlocal total_tests, passed_tests
        total_tests += 1
        try:
            test_func()
            print(f"✓ {description}")
            passed_tests += 1
        except Exception as e:
            print(f"✗ {description}: {e}")
    
    # Test ConfigValidator
    test("ConfigValidator type validation", lambda: (
        ConfigValidator.validate_type("test", str, "field") == "test" and
        ConfigValidator.validate_type(42, int, "field") == 42
    ))
    
    def test_invalid_type():
        try:
            ConfigValidator.validate_type(42, str, "field")
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError as e:
            if "must be of type str" not in str(e):
                raise AssertionError(f"Unexpected error message: {e}")
    
    test("ConfigValidator invalid type rejection", test_invalid_type)
    
    test("ConfigValidator positive number validation", lambda: (
        ConfigValidator.validate_positive_number(5, "field") == 5 and
        ConfigValidator.validate_positive_number(3.14, "field") == 3.14
    ))
    
    def test_negative_number():
        try:
            ConfigValidator.validate_positive_number(-5, "field")
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError as e:
            if "must be positive" not in str(e):
                raise AssertionError(f"Unexpected error message: {e}")
    
    test("ConfigValidator negative number rejection", test_negative_number)
    
    test("ConfigValidator range validation", lambda: (
        ConfigValidator.validate_range(5, 0, 10, "field") == 5 and
        ConfigValidator.validate_range(0, 0, 10, "field") == 0 and
        ConfigValidator.validate_range(10, 0, 10, "field") == 10
    ))
    
    def test_out_of_range():
        try:
            ConfigValidator.validate_range(-1, 0, 10, "field")
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError as e:
            if "must be between" not in str(e):
                raise AssertionError(f"Unexpected error message: {e}")
    
    test("ConfigValidator out of range rejection", test_out_of_range)
    
    test("ConfigValidator color validation", lambda: (
        ConfigValidator.validate_color("FF0000", "field") == "FF0000" and
        ConfigValidator.validate_color("00ff00", "field") == "00FF00" and  # Normalized
        ConfigValidator.validate_color("#123ABC", "field") == "123ABC"  # Hash removed
    ))
    
    def test_invalid_color():
        try:
            ConfigValidator.validate_color("invalid", "field")
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError as e:
            if "valid hex color" not in str(e):
                raise AssertionError(f"Unexpected error message: {e}")
    
    test("ConfigValidator invalid color rejection", test_invalid_color)
    
    test("ConfigValidator choices validation", lambda: (
        ConfigValidator.validate_choices("a", ["a", "b", "c"], "field") == "a" and
        ConfigValidator.validate_choices(1, [1, 2, 3], "field") == 1
    ))
    
    def test_invalid_choice():
        try:
            ConfigValidator.validate_choices("d", ["a", "b", "c"], "field")
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError as e:
            if "must be one of" not in str(e):
                raise AssertionError(f"Unexpected error message: {e}")
    
    test("ConfigValidator invalid choice rejection", test_invalid_choice)
    
    test("ConfigValidator opacity validation", lambda: (
        ConfigValidator.validate_opacity(0.0, "field") == 0.0 and
        ConfigValidator.validate_opacity(0.5, "field") == 0.5 and
        ConfigValidator.validate_opacity(1.0, "field") == 1.0
    ))
    
    def test_invalid_opacity():
        try:
            ConfigValidator.validate_opacity(1.5, "field")
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError as e:
            if "must be between" not in str(e):
                raise AssertionError(f"Unexpected error message: {e}")
    
    test("ConfigValidator invalid opacity rejection", test_invalid_opacity)
    
    # Test AtomType
    def test_valid_atom():
        atom = AtomType(name="C", color="909090", radius=0.4, max_num_bonds=4)
        return (atom.name == "C" and atom.color == "909090" and 
                atom.radius == 0.4 and atom.max_num_bonds == 4)
    
    test("AtomType valid creation", test_valid_atom)
    
    def test_atom_from_dict():
        data = {"name": "H", "color": "FFFFFF", "radius": 0.2, "max_num_bonds": 1}
        atom = AtomType.from_dict(data)
        return atom.name == "H" and atom.color == "FFFFFF"
    
    test("AtomType from_dict valid", test_atom_from_dict)
    
    def test_atom_missing_fields():
        try:
            data = {"name": "H", "color": "FFFFFF", "radius": 0.2}  # Missing max_num_bonds
            AtomType.from_dict(data)
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError as e:
            if "Missing required" not in str(e):
                raise AssertionError(f"Unexpected error message: {e}")
    
    test("AtomType missing fields rejection", test_atom_missing_fields)
    
    def test_atom_invalid_name():
        try:
            AtomType(name="", color="FF0000", radius=0.4, max_num_bonds=4)
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError as e:
            if "empty" not in str(e):
                raise AssertionError(f"Unexpected error message: {e}")
    
    test("AtomType invalid name rejection", test_atom_invalid_name)
    
    def test_atom_invalid_color():
        try:
            AtomType(name="C", color="invalid", radius=0.4, max_num_bonds=4)
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError as e:
            if "hex color" not in str(e):
                raise AssertionError(f"Unexpected error message: {e}")
    
    test("AtomType invalid color rejection", test_atom_invalid_color)
    
    def test_atom_invalid_radius():
        try:
            AtomType(name="C", color="FF0000", radius=-0.1, max_num_bonds=4)
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError as e:
            if "positive" not in str(e):
                raise AssertionError(f"Unexpected error message: {e}")
    
    test("AtomType invalid radius rejection", test_atom_invalid_radius)
    
    def test_atom_invalid_bonds():
        try:
            AtomType(name="C", color="FF0000", radius=0.4, max_num_bonds=15)
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError as e:
            if "between 0 and 10" not in str(e):
                raise AssertionError(f"Unexpected error message: {e}")
    
    test("AtomType invalid bonds rejection", test_atom_invalid_bonds)
    
    # Test Config
    def test_config_loading():
        config = Config()
        return (hasattr(config, 'grid') and hasattr(config, 'MO') and 
                hasattr(config, 'molecule') and len(config.atom_types) > 0)
    
    test("Config basic loading", test_config_loading)
    
    def test_grid_validation():
        config = Config()
        valid_grid = {"min_radius": 10, "max_radius_multiplier": 3}
        result = config._validate_grid_config(valid_grid)
        return result["min_radius"] == 10 and result["max_radius_multiplier"] == 3
    
    test("Config grid validation", test_grid_validation)
    
    def test_invalid_grid():
        config = Config()
        try:
            config._validate_grid_config({"min_radius": -5})
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError:
            pass  # Expected
    
    test("Config invalid grid rejection", test_invalid_grid)
    
    def test_mo_validation():
        config = Config()
        valid_mo = {"contour": 0.1, "opacity": 0.8}
        result = config._validate_mo_config(valid_mo)
        return result["contour"] == 0.1 and result["opacity"] == 0.8
    
    test("Config MO validation", test_mo_validation)
    
    def test_invalid_mo():
        config = Config()
        try:
            config._validate_mo_config({"opacity": 1.5})  # opacity > 1
            raise AssertionError("Should have raised ConfigurationError")
        except ConfigurationError:
            pass  # Expected
    
    test("Config invalid MO rejection", test_invalid_mo)
    
    print(f"\nTest Results: {passed_tests}/{total_tests} tests passed")
    
    if passed_tests == total_tests:
        print("✓ All configuration validation tests passed!")
        return True
    else:
        print("✗ Some tests failed")
        return False


if __name__ == "__main__":
    success = run_validation_tests()
    sys.exit(0 if success else 1)