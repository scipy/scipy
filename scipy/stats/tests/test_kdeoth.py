


from scipy import stats
import numpy as np
from numpy.testing import assert_almost_equal, assert_

def test_kde_1d():
    #some basic tests comparing to normal distribution
    np.random.seed(8765678)
    n_basesample = 500
    xn = np.random.randn(n_basesample)
    xnmean = xn.mean()
    xnstd = xn.std(ddof=1)

    # get kde for original sample
    gkde = stats.gaussian_kde(xn)
    
    # evaluate the density funtion for the kde for some points
    xs = np.linspace(-7,7,501)
    kdepdf = gkde.evaluate(xs)
    normpdf = stats.norm.pdf(xs, loc=xnmean, scale=xnstd)
    intervall = xs[1] - xs[0]
    
    assert_(np.sum((kdepdf - normpdf)**2)*intervall < 0.01)
    prob1 = gkde.integrate_box_1d(xnmean, np.inf)
    prob2 = gkde.integrate_box_1d(-np.inf, xnmean)
    assert_almost_equal(prob1, 0.5, decimal=1)
    assert_almost_equal(prob2, 0.5, decimal=1)
    assert_almost_equal(gkde.integrate_box(xnmean, np.inf), prob1, decimal=13)
    assert_almost_equal(gkde.integrate_box(-np.inf, xnmean), prob2, decimal=13)
    
    assert_almost_equal(gkde.integrate_kde(gkde),
                        (kdepdf**2).sum()*intervall, decimal=2)
    assert_almost_equal(gkde.integrate_gaussian(xnmean, xnstd**2),
                        (kdepdf*normpdf).sum()*intervall, decimal=2) 
