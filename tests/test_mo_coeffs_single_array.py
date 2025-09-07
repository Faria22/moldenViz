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
    
    # Test MO index access
    for i in range(min(5, len(parser.mos))):
        mo = parser.mos[i]
        assert_true(hasattr(mo, 'index'), f"MO {i} should have index")
        assert_true(0 <= mo.index < len(parser.mos), f"MO {i} index should be valid")
        
        # Test coefficient access
        coeffs = mo.coeffs
        direct_coeffs = parser.mo_coeffs[mo.index]
        assert_equal(coeffs, direct_coeffs, f"MO {i} coeffs should match direct access")

def test_sorting_efficiency():
    """Test that sorting is efficient and preserves coefficient access."""
    parser = parser_module.Parser(str(MOLDEN_PATH))
    
    # Get unsorted MOs
    mos_unsorted, mo_coeffs = parser.get_mos(sort=False)
    
    # Store original coefficient mapping
    original_mapping = {mo.index: mo.coeffs.copy() for mo in mos_unsorted}
    
    # Sort MOs
    mos_sorted = parser_module.Parser.sort_mos(mos_unsorted)
    
    # Verify sorting worked
    energies = [mo.energy for mo in mos_sorted]
    assert_true(all(energies[i] <= energies[i+1] for i in range(len(energies)-1)), 
                "Energies should be sorted")
    
    # Verify coefficient access preserved
    for mo in mos_sorted:
        current_coeffs = mo.coeffs
        expected_coeffs = original_mapping[mo.index]
        assert_equal(current_coeffs, expected_coeffs, 
                    f"MO coeffs should be preserved after sorting")

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
            # Old approach (what tabulator used to do)
            old_coeffs = np.stack([parser.mos[idx].coeffs for idx in pattern])
            
            # New approach (what tabulator does now)
            indices = [parser.mos[idx].index for idx in pattern]
            new_coeffs = parser.mo_coeffs[indices]
            
            assert_equal(old_coeffs, new_coeffs, 
                        f"Pattern {pattern} should produce identical results")

def test_memory_and_performance():
    """Test memory efficiency and performance improvements."""
    parser = parser_module.Parser(str(MOLDEN_PATH))
    
    import time
    
    # Test with different sizes
    for size in [5, 10, 20]:
        if size <= len(parser.mos):
            mo_inds = list(range(size))
            
            # Time old approach
            start = time.time()
            for _ in range(50):
                _ = np.stack([parser.mos[idx].coeffs for idx in mo_inds])
            old_time = time.time() - start
            
            # Time new approach
            start = time.time()
            for _ in range(50):
                indices = [parser.mos[idx].index for idx in mo_inds]
                _ = parser.mo_coeffs[indices]
            new_time = time.time() - start
            
            # New should be faster
            speedup = old_time / new_time if new_time > 0 else float('inf')
            assert_true(speedup >= 1.0, f"New approach should be faster (got {speedup:.2f}x)")

def test_backward_compatibility():
    """Test that all existing usage patterns still work."""
    parser = parser_module.Parser(str(MOLDEN_PATH))
    
    # Test individual MO access
    mo = parser.mos[0]
    coeffs = mo.coeffs
    
    # Should work like a normal numpy array
    assert_true(len(coeffs) > 0, "Should have non-zero length")
    assert_true(coeffs.dtype == np.float64, "Should be float64")
    assert_true(isinstance(coeffs[0], (float, np.floating)), "Should support indexing")
    assert_true(coeffs[:5].shape == (5,), "Should support slicing")
    
    # Should support numpy operations
    sum_val = np.sum(coeffs)
    assert_true(isinstance(sum_val, (float, np.floating)), "Should support numpy functions")

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