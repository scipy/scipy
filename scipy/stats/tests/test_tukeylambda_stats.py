
import numpy as np
from numpy.testing import assert_allclose, assert_equal

from scipy.stats._tukeylambda_stats import tukeylambda_variance, tukeylambda_kurtosis


def test_tukeylambda_stats_known_exact():
    """Compare results with some known exact formulas."""
    # Some exact values of the Tukey Lambda variance and kurtosis:
    # lambda   var      kurtosis
    #   0     pi**2/3     6/5     (logistic distribution)
    #  0.5    4 - pi    (5/3 - pi/2)/(pi/4 - 1)**2 - 3
    #   1      1/3       -6/5     (uniform distribution on (-1,1))
    #   2      1/12      -6/5     (uniform distribution on (-1/2, 1/2))

    # lambda = 0
    var = tukeylambda_variance(0)
    assert_allclose(var, np.pi**2 / 3, atol=1e-9)
    kurt = tukeylambda_kurtosis(0)
    assert_allclose(kurt, 1.2, atol=1e-9)

    # lambda = 0.5
    var = tukeylambda_variance(0.5)
    assert_allclose(var, 4 - np.pi, rtol=1e-8)
    kurt = tukeylambda_kurtosis(0.5)
    desired = (5./3 - np.pi/2) / (np.pi/4 - 1)**2 - 3
    assert_allclose(kurt, desired, atol=1e-9)

    # lambda = 1
    var = tukeylambda_variance(1)
    assert_allclose(var, 1.0 / 3, rtol=1e-8)
    kurt = tukeylambda_kurtosis(1)
    assert_allclose(kurt, -1.2, rtol=1e-8)    

    # lambda = 2
    var = tukeylambda_variance(2)
    assert_allclose(var, 1.0 / 12, rtol=1e-8)
    kurt = tukeylambda_kurtosis(2)
    assert_allclose(kurt, -1.2, rtol=1e-8)


def test_tukeylambda_stats_mpmath():
    """Compare results with some values that were computed using mpmath."""

    # For -0.00161 < lam < 0.00161, the absolute error of the variance should
    # be less than 1e-9.  
    #
    # For -0.065 < lam < 0.065, the absolute error of the kurtosis should
    # be less than 1e-10.
    #
    # Outside these intervals, we get whatever accuracy the regular formulas
    # provide.  We'll use rtol=1e-8 for these.

    r8 = dict(rtol=1e-8)
    a9 = dict(atol=1e-9, rtol=0)
    a10 = dict(atol=1e-10, rtol=0)
    data = [
        # lambda        variance                    kurtosis
        [-0.1,     4.78050217874253547,  r8,  3.78559520346454510,   r8],
        [-0.0649,  4.16428023599895777,  r8,  2.52019675947435718,  a10],
        [-0.001,   3.30128380368804765,  a9,  1.21452460083542988,  a10],
        [ 0.001,   3.27850775595073868,  a9,  1.18560634779287585,  a10],
        [ 0.03125, 2.95927803254613125,  r8,  0.804487555161819980, a10],
        [ 0.0649,  2.65282386754100378,  r8,  0.476834119532774540, a10],
        [ 1.2,     0.242153920578588392, r8, -1.23428047169049726,   r8],
        [ 10.0,  9.52375797577035934e-4, r8,  2.37810697355144933,   r8],
    ]

    for lam, var_expected, vartol, kurt_expected, kurttol in data:
        var = tukeylambda_variance(lam)
        assert_allclose(var, var_expected, **vartol)
        kurt = tukeylambda_kurtosis(lam)
        assert_allclose(kurt, kurt_expected, **kurttol)


def test_tukeylambda_stats_invalid():
    """Test values of lambda outside the domains of the functions."""
    lam = [-1.0, -0.5]
    var = tukeylambda_variance(lam)
    assert_equal(var, np.array([np.nan, np.inf]))

    lam = [-1.0, -0.25]
    kurt = tukeylambda_kurtosis(lam)
    assert_equal(kurt, np.array([np.nan, np.inf]))
