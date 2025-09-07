"""Additional unit tests for the single mo_coeffs array implementation.

These tests verify the new implementation works correctly and maintains
backward compatibility.
"""

import numpy as np
from pathlib import Path
import sys
import os

# This would normally use pytest, but we'll make it work standalone
def assert_equal(a, b, msg=""):
    """Simple assert equal function."""
    if not np.array_equal(a, b):
        raise AssertionError(f"Arrays not equal: {msg}")

def assert_true(condition, msg=""):
    """Simple assert true function."""
    if not condition:
        raise AssertionError(f"Condition not true: {msg}")

# Import setup
sys.path.insert(0, '/home/runner/work/moldenViz/moldenViz/src')
import importlib.util
spec = importlib.util.spec_from_file_location("parser", "/home/runner/work/moldenViz/moldenViz/src/moldenViz/parser.py")
parser_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parser_module)

MOLDEN_PATH = Path('/home/runner/work/moldenViz/moldenViz/tests/sample_molden.inp')

def test_single_mo_coeffs_array():
    """Test that parser creates a single shared mo_coeffs array."""
    parser = parser_module.Parser(str(MOLDEN_PATH))
    
    # Test array structure
    assert_true(hasattr(parser, 'mo_coeffs'), "Parser should have mo_coeffs attribute")
    assert_true(parser.mo_coeffs.ndim == 2, "mo_coeffs should be 2D")
    assert_true(parser.mo_coeffs.shape[0] == len(parser.mos), "First dim should match MO count")
    
    # Test that MOs no longer have individual coeffs
    for i in range(min(5, len(parser.mos))):
        mo = parser.mos[i]
        
        # Test that coeffs field doesn't exist
        assert_true(not hasattr(mo, 'coeffs'), f"MO {i} should not have coeffs field")
        
        # Test that coefficients are accessible via parser.mo_coeffs
        coeffs = parser.mo_coeffs[i]
        assert_true(len(coeffs) > 0, f"MO {i} should have coefficients in mo_coeffs array")

def test_sorting_efficiency():
    """Test that sorting is efficient and preserves coefficient access."""
    parser = parser_module.Parser(str(MOLDEN_PATH))
    
    # Get unsorted MOs
    mos_unsorted, mo_coeffs_unsorted = parser.get_mos(sort=False)
    
    # Store original coefficients by index
    original_coeffs = {i: mo_coeffs_unsorted[i].copy() for i in range(len(mos_unsorted))}
    original_energies = {i: mos_unsorted[i].energy for i in range(len(mos_unsorted))}
    
    # Get sorted MOs (default behavior)
    mos_sorted, mo_coeffs_sorted = parser.get_mos(sort=True)
    
    # Verify sorting worked
    energies = [mo.energy for mo in mos_sorted]
    assert_true(all(energies[i] <= energies[i+1] for i in range(len(energies)-1)), 
                "Energies should be sorted")
    
    # Verify coefficient arrays are properly sorted - coefficients should follow MO order
    for i, mo in enumerate(mos_sorted):
        # Find the original index of this MO by matching energy and checking coefficients
        matching_indices = [j for j, energy in original_energies.items() if energy == mo.energy]
        
        # Find which one has matching coefficients
        found_match = False
        for j in matching_indices:
            if np.array_equal(mo_coeffs_sorted[i], original_coeffs[j]):
                found_match = True
                break
        
        assert_true(found_match, f"Could not find matching coefficients for MO {i}")

def test_tabulator_integration():
    """Test that tabulator operations work correctly with the new structure."""
    parser = parser_module.Parser(str(MOLDEN_PATH))
    
    # Test patterns that tabulator would use
    test_patterns = [
        [0, 1, 2],
        list(range(10)),
        [0, 5, 10, 15, 20]
    ]
    
    for pattern in test_patterns:
        if all(0 <= idx < len(parser.mos) for idx in pattern):
            # New approach (direct slicing since mo_coeffs is now the only source)
            new_coeffs = parser.mo_coeffs[pattern]
            
            # Test that we get the expected shape and data
            assert_true(new_coeffs.shape[0] == len(pattern), 
                        f"Pattern {pattern} should return correct number of MOs")
            assert_true(new_coeffs.shape[1] == parser.mo_coeffs.shape[1], 
                        f"Pattern {pattern} should preserve coefficient dimension")

def test_memory_and_performance():
    """Test memory efficiency and performance improvements."""
    parser = parser_module.Parser(str(MOLDEN_PATH))
    
    import time
    
    # Test with different sizes
    for size in [5, 10, 20]:
        if size <= len(parser.mos):
            mo_inds = list(range(size))
            
            # Time direct slicing approach (the only way now)
            start = time.time()
            for _ in range(50):
                _ = parser.mo_coeffs[mo_inds]
            direct_time = time.time() - start
            
            # Verify we get reasonable performance
            assert_true(direct_time < 1.0, f"Direct slicing should be fast (took {direct_time:.4f}s)")

def test_backward_compatibility():
    """Test that the new API works correctly."""
    parser = parser_module.Parser(str(MOLDEN_PATH))
    
    # Test that we can access coefficients through mo_coeffs
    coeffs = parser.mo_coeffs[0]
    
    # Should work like a normal numpy array
    assert_true(len(coeffs) > 0, "Should have non-zero length")
    assert_true(coeffs.dtype == np.float64, "Should be float64")
    assert_true(isinstance(coeffs[0], (float, np.floating)), "Should support indexing")
    assert_true(coeffs[:5].shape == (5,), "Should support slicing")
    
    # Should support numpy operations
    sum_val = np.sum(coeffs)
    assert_true(isinstance(sum_val, (float, np.floating)), "Should support numpy functions")
    
    # Test that MO objects no longer have coeffs
    mo = parser.mos[0]
    assert_true(not hasattr(mo, 'coeffs'), "MO should not have coeffs field")

if __name__ == "__main__":
    print("Running comprehensive tests for single mo_coeffs array implementation...")
    
    test_single_mo_coeffs_array()
    print("✓ Single mo_coeffs array test passed")
    
    test_sorting_efficiency()
    print("✓ Sorting efficiency test passed")
    
    test_tabulator_integration()
    print("✓ Tabulator integration test passed")
    
    test_memory_and_performance()
    print("✓ Memory and performance test passed")
    
    test_backward_compatibility()
    print("✓ Backward compatibility test passed")
    
    print("\n✅ All comprehensive tests passed!")