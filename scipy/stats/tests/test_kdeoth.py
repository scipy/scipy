from scipy import stats
import numpy as np
from numpy.testing import assert_almost_equal, assert_, \
    assert_array_almost_equal, assert_array_almost_equal_nulp


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


def test_kde_bandwidth_method():
    def scotts_factor(kde_obj):
        """Same as default, just check that it works."""
        return np.power(kde_obj.n, -1./(kde_obj.d+4))

    np.random.seed(8765678)
    n_basesample = 50
    xn = np.random.randn(n_basesample)

    gkde = stats.gaussian_kde(xn)
    gkde2 = stats.gaussian_kde(xn, bw_method=scotts_factor)

    xs = np.linspace(-7,7,51)
    kdepdf = gkde.evaluate(xs)
    kdepdf2 = gkde2.evaluate(xs)
    assert_almost_equal(gkde(xs), gkde2(xs))


# Subclasses that should stay working (extracted from various sources).
# Unfortunately the earlier design of gaussian_kde made it necessary for users
# to create these kinds of subclasses, or call _compute_covariance() directly.

class _kde_subclass1(stats.gaussian_kde):
    def __init__(self, dataset):
        self.dataset = np.atleast_2d(dataset)
        self.d, self.n = self.dataset.shape
        self.covariance_factor = self.scotts_factor
        self._compute_covariance()

class _kde_subclass2(stats.gaussian_kde):
    def __init__(self, dataset):
        self.covariance_factor = self.scotts_factor
        super(_kde_subclass2, self).__init__(dataset)

class _kde_subclass3(stats.gaussian_kde):
    def __init__(self, dataset, covariance):
        self.covariance = covariance
        stats.gaussian_kde.__init__(self, dataset)

    def _compute_covariance(self):
        self.inv_cov = np.linalg.inv(self.covariance)
        self._norm_factor = np.sqrt(np.linalg.det(2*np.pi * self.covariance)) \
                                   * self.n

class _kde_subclass4(stats.gaussian_kde):
    def __init__(self, dataset):
        self.covariance_factor = self.silverman_factor
        stats.gaussian_kde.__init__(self, dataset)


def test_gaussian_kde_subclassing():
    x1 = np.array([-7, -5, 1, 4, 5], dtype=np.float)
    xs = np.linspace(-10, 10, num=50)

    # gaussian_kde itself
    kde = stats.gaussian_kde(x1)
    ys = kde(xs)

    # subclass 1
    kde1 = _kde_subclass1(x1)
    y1 = kde1(xs)
    assert_array_almost_equal_nulp(ys, y1, nulp=10)

    # subclass 2
    kde2 = _kde_subclass2(x1)
    y2 = kde2(xs)
    assert_array_almost_equal_nulp(ys, y2, nulp=10)

    # subclass 3
    kde3 = _kde_subclass3(x1, kde.covariance)
    y3 = kde3(xs)
    assert_array_almost_equal_nulp(ys, y3, nulp=10)

    # subclass 4, decimal=2 due to scott vs. silverman
    kde4 = _kde_subclass4(x1)
    y4 = kde4(xs)
    assert_array_almost_equal(ys, ys, decimal=2)

    # Not a subclass, but check for use of _compute_covariance()
    kde5 = kde
    kde5.covariance_factor = lambda: kde.factor
    kde._compute_covariance()
    assert_array_almost_equal_nulp(ys, y1, nulp=10)

