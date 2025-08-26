"""Test that trapezoid.fit now converges properly with default starting values."""

import numpy as np
import scipy.stats
import pytest


def test_trapezoid_fit_convergence():
    """Test that trapezoid.fit converges without starting values."""
    
    # Generate test data from a trapezoidal distribution
    true_c, true_d = 0.3, 0.7
    true_dist = scipy.stats.trapezoid(true_c, true_d, 0, 1)
    rvs = true_dist.rvs(1000, random_state=42)
    
    # Test 1: Fit without starting values (should now work)
    fitted_args = scipy.stats.trapezoid.fit(rvs, floc=0, scale=1)
    fitted_c, fitted_d = fitted_args[0], fitted_args[1]
    
    # Should not converge to triangular distribution (c=d=1)
    assert not (abs(fitted_c - 1.0) < 0.01 and abs(fitted_d - 1.0) < 0.01), \
        "fit() still converging to triangular distribution (c=d=1)"
    
    # Should converge to reasonable values
    assert 0 <= fitted_c <= 1, f"fitted c={fitted_c} not in [0,1]"
    assert 0 <= fitted_d <= 1, f"fitted d={fitted_d} not in [0,1]"
    assert fitted_d >= fitted_c, f"fitted d={fitted_d} < c={fitted_c}"


def test_trapezoid_fit_with_explicit_values():
    """Test that trapezoid.fit still works with explicit starting values."""
    
    # Generate test data
    true_dist = scipy.stats.trapezoid(0.4, 0.6, 0, 1)
    rvs = true_dist.rvs(500, random_state=42)
    
    # Test with explicit starting values
    fitted_args = scipy.stats.trapezoid.fit(rvs, 0.25, 0.75, floc=0, scale=1)
    fitted_c, fitted_d = fitted_args[0], fitted_args[1]
    
    # Should work as before
    assert 0 <= fitted_c <= 1
    assert 0 <= fitted_d <= 1
    assert fitted_d >= fitted_c


def test_trapezoid_fit_methods():
    """Test that both MLE and MM methods work."""
    
    # Generate test data
    true_dist = scipy.stats.trapezoid(0.35, 0.65, 0, 1)
    rvs = true_dist.rvs(800, random_state=42)
    
    # Test MLE method
    mle_args = scipy.stats.trapezoid.fit(rvs, floc=0, scale=1, method='MLE')
    mle_c, mle_d = mle_args[0], mle_args[1]
    assert 0 <= mle_c <= 1 and 0 <= mle_d <= 1 and mle_d >= mle_c
    
    # Test MM method
    mm_args = scipy.stats.trapezoid.fit(rvs, floc=0, scale=1, method='MM')
    mm_c, mm_d = mm_args[0], mm_args[1]
    assert 0 <= mm_c <= 1 and 0 <= mm_d <= 1 and mm_d >= mm_c


def test_trapezoid_fit_edge_cases():
    """Test edge cases for trapezoid.fit."""
    
    # Edge case: c=0, d=1 (uniform distribution)
    edge_dist = scipy.stats.trapezoid(0, 1, 0, 1)
    edge_rvs = edge_dist.rvs(300, random_state=42)
    
    edge_fitted = scipy.stats.trapezoid.fit(edge_rvs, floc=0, scale=1)
    edge_c, edge_d = edge_fitted[0], edge_fitted[1]
    
    # Should handle edge case properly
    assert 0 <= edge_c <= 1 and 0 <= edge_d <= 1 and edge_d >= edge_c


if __name__ == "__main__":
    # Run tests
    test_trapezoid_fit_convergence()
    test_trapezoid_fit_with_explicit_values()
    test_trapezoid_fit_methods()
    test_trapezoid_fit_edge_cases()
    print("All tests passed! The trapezoid.fit fix is working correctly.")
