import numpy as np
import pytest
from numpy.testing import (
    assert_allclose, assert_array_equal, assert_, assert_equal
)

from scipy.stats._grenander import grenander


def _integral_of_hist(knots, slopes):
    widths = np.diff(knots)
    return np.sum(widths * slopes)

def test_input_validation():
    rng = np.random.default_rng(0)

    # empty
    with pytest.raises(ValueError, match="no observations"):
        grenander([])

    x = rng.exponential(size=10)

    # support_min strict inequality
    with pytest.raises(ValueError, match="strictly greater"):
        grenander(np.r_[0.0, x], support_min=0.0)

    with pytest.raises(ValueError, match="strictly greater"):
        grenander(x, support_min=x.min())

    # duplicates forbidden
    x2 = np.array([0.1, 0.2, 0.2, 0.3])
    with pytest.raises(ValueError, match="duplicate"):
        grenander(x2)

    # if assume_sorted=True, duplicates should still be caught
    x3 = np.array([0.1, 0.2, 0.2, 0.3])
    with pytest.raises(ValueError, match="duplicate"):
        grenander(x3, assume_sorted=True)


def test_assume_sorted_equivalence_and_fitted_ordering():
    rng = np.random.default_rng(123)
    x = rng.exponential(size=200)

    g_unsorted = grenander(x, assume_sorted=False)

    xs = np.sort(x)
    g_sorted = grenander(xs, assume_sorted=True)

    # "function" parts should agree (cdf/pdf), regardless of data ordering
    grid = np.linspace(0, xs.max() * 1.1, 200)
    assert_allclose(g_unsorted.cdf(grid), g_sorted.cdf(grid), 
                    rtol=0, atol=1e-14)
    assert_allclose(g_unsorted.pdf(grid), g_sorted.pdf(grid), 
                    rtol=0, atol=1e-14)

    # fitted is returned in original order for assume_sorted=False;
    # it should match the sorted version when permuted by argsort.
    I = np.argsort(x)
    assert_allclose(g_unsorted.fitted[I], g_sorted.fitted, rtol=0, atol=0)


def test_basic_shapes_and_monotonicity_invariants():
    rng = np.random.default_rng(1)
    x = rng.exponential(size=300)
    g = grenander(x)

    n = x.size

    # shapes
    assert_equal(g.fitted.shape, (n,))
    assert_(g.knots.ndim == 1)
    assert_(g.heights.ndim == 1)
    assert_(g.slopes.ndim == 1)
    assert_equal(g.knots.size, g.slopes.size + 1)
    assert_equal(g.heights.size, g.knots.size)

    # knots increasing, heights increasing, slopes nonincreasing
    assert_(np.all(np.diff(g.knots) > 0))
    assert_(np.all(np.diff(g.heights) > 0))
    assert_(np.all(np.diff(g.slopes) < 0))

    # first knot is support_min; last height is 1
    assert_allclose(g.knots[0], g.support_min, rtol=0, atol=0)
    assert_allclose(g.heights[-1], 1.0, rtol=0, atol=0)

    # integral is 1 (up to float error)
    assert_allclose(_integral_of_hist(g.knots, g.slopes), 1.0, 
                    rtol=1e-13, atol=1e-13)


def test_pdf_cdf_boundary_behavior():
    rng = np.random.default_rng(2)
    x = rng.exponential(size=100)
    g = grenander(x, support_min=0.0)

    # below support_min: cdf=0, pdf=0
    assert_allclose(g.cdf([-1.0, -0.5, 0.0]), [0.0, 0.0, 0.0])
    assert_allclose(g.pdf([-1.0, -0.5]), [0.0, 0.0])

    # above last knot: cdf=1, pdf=0
    big = g.knots[-1] + 10.0
    assert_allclose(g.cdf(big), 1.0)
    assert_allclose(g.pdf(big), 0.0)


def test_cdf_is_consistent_with_histogram_representation():
    rng = np.random.default_rng(4)
    x = rng.exponential(size=250)
    g = grenander(x, support_min=0.0)

    # On any point inside an interval [k_j, k_{j+1}),
    # cdf(x) = height[j] + slope[j]*(x - k_j)
    # We'll test midpoints of each bin.
    mid = 0.5 * (g.knots[:-1] + g.knots[1:])
    expected = g.heights[:-1] + g.slopes * (mid - g.knots[:-1])
    assert_allclose(g.cdf(mid), expected, rtol=0, atol=2e-14)

    # CDF right at knots should match heights (linear interpolation)
    assert_allclose(g.cdf(g.knots), g.heights, rtol=0, atol=0)


def test_pdf_cdf_compatibility_at_boundary():
    """Test that PDF and CDF are consistent at the support boundary."""
    rng = np.random.default_rng(1234567891)
    x = rng.exponential(size=100)
    g = grenander(x)

    k = g.knots[-1]  # Upper boundary of support
    eps = 1e-10

    # Test points: just before, at, and just after boundary
    z_before = k - eps
    z_at = k
    z_after = k + eps

    # PDF should be non-zero at and just before boundary
    assert g.pdf(z_before) > 0, "PDF should be positive just before boundary"
    assert g.pdf(z_at) > 0, "PDF should be positive at boundary"

    # PDF should be zero after boundary
    assert g.pdf(z_after) == 0, "PDF should be zero after boundary"

    # CDF should reach 1.0 at the boundary
    assert_allclose(g.cdf(z_at), 1.0, rtol=1e-10)

    # CDF should stay at 1.0 after boundary
    assert_allclose(g.cdf(z_after), 1.0, rtol=1e-10)

    # PDF and CDF should be consistent: where PDF=0, CDF should be flat
    # i.e., derivative of CDF should be zero where PDF is zero
    cdf_before = g.cdf(z_before)
    cdf_at = g.cdf(z_at)
    cdf_after = g.cdf(z_after)

    # CDF should increase from before to at (since PDF > 0)
    assert cdf_at > cdf_before, "CDF should increase where PDF > 0"

    # CDF should be flat from at to after (since PDF = 0)
    assert_allclose(cdf_at, cdf_after, rtol=1e-10,
                    atol=0, err_msg="CDF should be constant where PDF = 0")