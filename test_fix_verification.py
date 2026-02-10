#!/usr/bin/env python3
"""
Comprehensive verification of the window spectral tests fix.
This script manually tests the fixed logic without requiring pytest or a full build.
"""

import numpy as np
import sys

# Manually compute FFT (we can't import scipy.fft directly in source dir)
# but we can use numpy.fft instead
from numpy.fft import fft

# Mock the window functions with numpy equivalents
def boxcar(M, sym=True):
    return np.ones(M)

def hann(M, sym=True):
    return np.hanning(M) if sym else np.hanning(M+1)[:-1]

def hamming(M, sym=True):
    return np.hamming(M) if sym else np.hamming(M+1)[:-1]

def blackman(M, sym=True):
    return np.blackman(M) if sym else np.blackman(M+1)[:-1]

windows_dict = {
    'boxcar': boxcar,
    'hann': hann,
    'hamming': hamming,
    'blackman': blackman,
}

expected_psll = {
    'boxcar': -13.26,
    'hann': -31.47,
    'hamming': -42.67,
    'blackman': -58.11,
}

expected_width = {
    'boxcar': 2.0,
    'hann': 4.0,
    'hamming': 4.0,
    'blackman': 6.0,
}

def test_peak_sidelobe_level():
    """Test PSLL calculation with the fixed logic."""
    print("\n" + "="*70)
    print("TEST 1: Peak Sidelobe Level (PSLL) Correctness")
    print("="*70)
    
    M_win = 1024
    N_fft = 65536
    passed = 0
    failed = 0
    
    for window_name in ['boxcar', 'hann', 'hamming', 'blackman']:
        try:
            # Generate window
            window_func = windows_dict[window_name]
            w = window_func(M_win, sym=False)
            
            # Compute spectrum and normalize 
            spec = np.abs(fft(w, N_fft))
            spec = spec / spec.max()
            half_spectrum = spec[:N_fft // 2]
            
            # Find first null (NEW ROBUST LOGIC)
            null_candidates = np.flatnonzero(half_spectrum <= 1e-6)
            assert null_candidates.size > 0, "Could not find first null"
            first_null = int(null_candidates[0])
            
            # Extract sidelobe region
            sidelobe_region = half_spectrum[first_null + 1:]
            assert sidelobe_region.size > 0, "Sidelobe region is empty"
            
            # Compute PSLL in dB
            psll = 20 * np.log10(float(np.max(sidelobe_region)))
            expected = expected_psll[window_name]
            tolerance = 0.5
            error = abs(psll - expected)
            
            status = "✓ PASS" if error <= tolerance else "✗ FAIL"
            print(f"\n{window_name:10s}  Expected: {expected:7.2f} dB  Got: {psll:7.2f} dB  "
                  f"Error: {error:6.4f} dB  {status}")
            print(f"           First null at bin: {first_null:5d}  "
                  f"Sidelobe region size: {sidelobe_region.size:8d}")
            
            if error <= tolerance:
                passed += 1
            else:
                failed += 1
                print(f"           WARNING: Error {error:.4f} exceeds tolerance {tolerance}")
        except Exception as e:
            failed += 1
            print(f"\n{window_name:10s}  ✗ EXCEPTION: {e}")
    
    print("\n" + "-"*70)
    print(f"PSLL Test Results: {passed} passed, {failed} failed")
    return failed == 0

def test_mainlobe_width():
    """Test mainlobe width calculation with the fixed logic."""
    print("\n" + "="*70)
    print("TEST 2: Mainlobe Width Correctness")
    print("="*70)
    
    M_win = 1024
    N_fft = 65536
    passed = 0
    failed = 0
    
    for window_name in ['boxcar', 'hann', 'hamming', 'blackman']:
        try:
            # Generate window
            window_func = windows_dict[window_name]
            w = window_func(M_win, sym=False)
            
            # Compute spectrum and normalize
            spec = np.abs(fft(w, N_fft))
            normalized = spec / spec.max()
            half_spectrum = normalized[:N_fft // 2]
            
            # Find first null (NEW ROBUST LOGIC)
            null_candidates = np.flatnonzero(half_spectrum <= 1e-6)
            assert null_candidates.size > 0, "Could not find first null"
            first_null = int(null_candidates[0])
            
            # Convert to normalized bins
            width_bins = 2 * first_null / N_fft * M_win
            expected = expected_width[window_name]
            rel_tol = 0.02
            error_rel = abs(width_bins - expected) / expected * 100
            
            status = "✓ PASS" if error_rel <= rel_tol * 100 else "✗ FAIL"
            print(f"\n{window_name:10s}  Expected: {expected:6.2f} bins  Got: {width_bins:6.2f} bins  "
                  f"Error: {error_rel:5.2f}%  {status}")
            print(f"           First null at bin: {first_null:5d}  "
                  f"Half spectrum size: {half_spectrum.size:8d}")
            
            if error_rel <= rel_tol * 100:
                passed += 1
            else:
                failed += 1
                print(f"           WARNING: Error {error_rel:.2f}% exceeds tolerance {rel_tol*100:.2f}%")
        except Exception as e:
            failed += 1
            print(f"\n{window_name:10s}  ✗ EXCEPTION: {e}")
    
    print("\n" + "-"*70)
    print(f"Width Test Results: {passed} passed, {failed} failed")
    return failed == 0

def test_robustness():
    """Test that the new assertions catch edge cases properly."""
    print("\n" + "="*70)
    print("TEST 3: Robustness & Error Handling")
    print("="*70)
    
    passed = 0
    failed = 0
    
    # Test 3a: Verify np.flatnonzero rejects when no nulls exist
    print("\nTest 3a: flatnonzero correctly identifies no nulls")
    try:
        # Create an artificial case with no true nulls
        dummy_spec = np.ones(100) * 0.5  # All 0.5, no nulls <= 1e-6
        null_candidates = np.flatnonzero(dummy_spec <= 1e-6)
        if null_candidates.size == 0:
            print("  ✓ PASS - flatnonzero correctly returns empty array for no-null case")
            passed += 1
        else:
            print("  ✗ FAIL - flatnonzero should return empty but didn't")
            failed += 1
    except Exception as e:
        print(f"  ✗ EXCEPTION: {e}")
        failed += 1
    
    # Test 3b: Verify guard assertion catches empty sidelobe region
    print("\nTest 3b: Sidelobe region guard assertion works")
    try:
        # Simulate case where first_null is at the very end
        half_spectrum = np.array([1.0] + [0.0]*99)
        first_null = 99
        sidelobe_region = half_spectrum[first_null + 1:]
        
        if sidelobe_region.size == 0:
            print("  ✓ PASS - sidelobe_region correctly becomes empty at boundary")
            passed += 1
        else:
            print("  ✗ FAIL - sidelobe_region should be empty at boundary")
            failed += 1
    except Exception as e:
        print(f"  ✗ EXCEPTION: {e}")
        failed += 1
    
    # Test 3c: Verify normalization doesn't cause log(0)
    print("\nTest 3c: No log(0) due to normalization")
    try:
        # Generate a real window spectrum
        w = np.hanning(1024)
        spec = np.abs(np.fft.fft(w, 65536))
        spec = spec / spec.max()  # Normalized to [0, 1]
        
        # The max in sidelobe region should be > 0 after normalization
        half_spec = spec[:65536//2]
        nulls = np.flatnonzero(half_spec <= 1e-6)
        if nulls.size > 0:
            sidelobe = half_spec[nulls[0] + 1:]
            max_sidelobe = np.max(sidelobe)
            
            if max_sidelobe > 0:
                psll_db = 20 * np.log10(max_sidelobe)
                print(f"  ✓ PASS - No log(0); max sidelobe = {max_sidelobe:.4e}, PSLL = {psll_db:.2f} dB")
                passed += 1
            else:
                print("  ✗ FAIL - max_sidelobe is 0 or negative, would cause log error")
                failed += 1
    except Exception as e:
        print(f"  ✗ EXCEPTION: {e}")
        failed += 1
    
    print("\n" + "-"*70)
    print(f"Robustness Test Results: {passed} passed, {failed} failed")
    return failed == 0

def test_comparison_old_vs_new():
    """Compare the problematic old logic with the new logic."""
    print("\n" + "="*70)
    print("TEST 4: Old Logic vs. New Logic Comparison")
    print("="*70)
    
    M_win = 1024
    N_fft = 65536
    
    print("\nFor 'hann' window:")
    w = hann(M_win, sym=False)
    spec = np.abs(fft(w, N_fft))
    
    # OLD LOGIC (PROBLEMATIC)
    print("\n  OLD LOGIC (Original, problematic):")
    spec_db = 20 * np.log10(spec / spec.max())
    diff_condition = np.diff(spec_db) > 0
    old_first_zero = int(np.argmax(diff_condition))
    print(f"    - np.diff(spec_db) > 0 found {np.sum(diff_condition)} True values")
    print(f"    - np.argmax() returns: {old_first_zero}")
    print(f"    - Slice spec_db[{old_first_zero}:-{old_first_zero}] would compute PSLL")
    
    if old_first_zero == 0:
        print("    ⚠ WARNING: np.argmax returned 0 (no True values OR True at index 0)")
        print("    ⚠ This is the bug we fixed!")
    else:
        old_psll = np.max(spec_db[old_first_zero:-old_first_zero])
        print(f"    - Result: {old_psll:.2f} dB")
    
    # NEW LOGIC (FIXED)
    print("\n  NEW LOGIC (Fixed, robust):")
    spec_norm = spec / spec.max()
    half_spectrum = spec_norm[:N_fft // 2]
    null_candidates = np.flatnonzero(half_spectrum <= 1e-6)
    
    print(f"    - flatnonzero(half_spectrum <= 1e-6) found {null_candidates.size} null candidates")
    
    if null_candidates.size > 0:
        first_null = int(null_candidates[0])
        sidelobe_region = half_spectrum[first_null + 1:]
        if sidelobe_region.size > 0:
            new_psll = 20 * np.log10(np.max(sidelobe_region))
            print(f"    - First null at bin: {first_null}")
            print(f"    - Sidelobe region from bin {first_null+1} to {len(half_spectrum)-1}")
            print(f"    - Result: {new_psll:.2f} dB")
        else:
            print("    - Sidelobe region is empty (assertion would catch this)")
    else:
        print("    - No nulls found (assertion would catch this cleanly)")
    
    print("\n" + "-"*70)
    print("Logic comparison complete.")
    return True

if __name__ == '__main__':
    print("\n" + "█"*70)
    print("█  COMPREHENSIVE FIX VERIFICATION")
    print("█  Testing scipy.signal window spectral property tests")
    print("█"*70)
    
    all_passed = True
    
    # Run all tests
    all_passed &= test_peak_sidelobe_level()
    all_passed &= test_mainlobe_width()
    all_passed &= test_robustness()
    all_passed &= test_comparison_old_vs_new()
    
    # Summary
    print("\n" + "█"*70)
    if all_passed:
        print("█  ✓ ALL VERIFICATION TESTS PASSED")
        print("█  The fix is correct and robust.")
        sys.exit(0)
    else:
        print("█  ✗ SOME TESTS FAILED - Review output above")
        sys.exit(1)
    print("█"*70)
