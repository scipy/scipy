"""
Test functions for multivariate normal distributions.

"""
from __future__ import division, print_function, absolute_import

from numpy.testing import run_module_suite, assert_allclose, assert_equal

import numpy
import numpy as np

import scipy.linalg
import scipy.stats._multivariate
from scipy.stats import multivariate_normal
from scipy.stats import norm

from scipy.stats._multivariate import _psd_pinv_decomposed_log_pdet

from scipy.integrate import romb


def test_scalar_values():
    # When evaluated on scalar data, the pdf should return a scalar
    x, mean, cov = 1.5, 1.7, 2.5
    pdf = multivariate_normal.pdf(x, mean, cov)
    assert(pdf.ndim == 0)

    # When evaluated on a single vector, the pdf should return a scalar
    x = np.random.randn(5)
    mean = np.random.randn(5)
    cov = np.abs(np.random.randn(5))  # Diagonal values for cov. matrix
    pdf = multivariate_normal.pdf(x, mean, cov)
    assert(pdf.ndim == 0)


def test_logpdf():
    # Check that the log of the pdf is in fact the logpdf
    x = np.random.randn(5)
    mean = np.random.randn(5)
    cov = np.abs(np.random.randn(5))
    d1 = multivariate_normal.logpdf(x, mean, cov)
    d2 = multivariate_normal.pdf(x, mean, cov)
    assert_allclose(d1, np.log(d2))


def test_large_pseudo_determinant():
    # Check that large pseudo-determinants are handled appropriately.

    # Construct a singular diagonal covariance matrix
    # whose pseudo determinant overflows double precision.
    large_total_log = 1000.0
    npos = 100
    nzero = 2
    large_entry = np.exp(large_total_log / npos)
    n = npos + nzero
    cov = np.zeros((n, n), dtype=float)
    np.fill_diagonal(cov, large_entry)
    cov[-nzero:, -nzero:] = 0

    # Check some determinants.
    assert_equal(scipy.linalg.det(cov), 0)
    assert_equal(scipy.linalg.det(cov[:npos, :npos]), np.inf)

    # np.linalg.slogdet is only available in numpy 1.6+
    # but scipy currently supports numpy 1.5.1.
    #assert_allclose(np.linalg.slogdet(cov[:npos, :npos]), (1, large_total_log))

    # Check the pseudo-determinant.
    U, log_pdet = scipy.stats._multivariate._psd_pinv_decomposed_log_pdet(cov)
    assert_allclose(log_pdet, large_total_log)


def test_broadcasting():
    n = 4

    # Construct a random covariance matrix.
    data = np.random.randn(n, n)
    cov = np.dot(data, data.T)
    mean = np.random.randn(n)

    # Construct an ndarray which can be interpreted as
    # a 2x3 array whose elements are random data vectors.
    X = np.random.randn(2, 3, n)

    # Check that multiple data points can be evaluated at once.
    for i in range(2):
        for j in range(3):
            actual = multivariate_normal.pdf(X[i, j], mean, cov)
            desired = multivariate_normal.pdf(X, mean, cov)[i, j]
            assert_allclose(actual, desired)


def test_normal_1D():
    # The probability density function for a 1D normal variable should
    # agree with the standard normal distribution in scipy.stats.distributions
    x = np.linspace(0, 2, 10)
    mean, cov = 1.2, 0.9
    scale=cov**0.5
    d1 = norm.pdf(x, mean, scale)
    d2 = multivariate_normal.pdf(x, mean, cov)
    assert_allclose(d1, d2)


def test_marginalization():
    # Integrating out one of the variables of a 2D Gaussian should
    # yield a 1D Gaussian
    mean = np.array([2.5, 3.5])
    cov = np.array([[.5, 0.2], [0.2, .6]])
    n = 2**8 + 1  # Number of samples
    delta = 6 / (n - 1)  # Grid spacing

    v = np.linspace(0, 6, n)
    xv, yv = np.meshgrid(v, v)
    pos = np.empty((n, n, 2))
    pos[:, :, 0] = xv
    pos[:, :, 1] = yv
    pdf = multivariate_normal.pdf(pos, mean, cov)

    # Marginalize over x and y axis
    margin_x = romb(pdf, delta, axis=0)
    margin_y = romb(pdf, delta, axis=1)

    # Compare with standard normal distribution
    gauss_x = norm.pdf(v, loc=mean[0], scale=cov[0, 0]**0.5)
    gauss_y = norm.pdf(v, loc=mean[1], scale=cov[1, 1]**0.5)
    assert_allclose(margin_x, gauss_x, rtol=1e-2, atol=1e-2)
    assert_allclose(margin_y, gauss_y, rtol=1e-2, atol=1e-2)


def test_frozen():
    # The frozen distribution should agree with the regular one
    x = np.random.randn(5)
    mean = np.random.randn(5)
    cov = np.abs(np.random.randn(5))
    norm_frozen = multivariate_normal(mean, cov)
    assert_allclose(norm_frozen.pdf(x), multivariate_normal.pdf(x, mean, cov))
    assert_allclose(norm_frozen.logpdf(x),
                    multivariate_normal.logpdf(x, mean, cov))


def test_pseudodet_pinv():
    # Make sure that pseudo-inverse and pseudo-det agree on cutoff

    # Assemble covariance matrix with large and small eigenvalues
    n = 7
    x = np.random.randn(n, n)
    cov = np.dot(x, x.T)
    s, u = scipy.linalg.eigh(cov)
    s = 0.5 * np.ones(n)
    s[0] = 1.0; s[-1] = 1e-7
    cov = np.dot(u, np.dot(np.diag(s), u.T))

    # Set cond so that the lowest eigenvalue is below the cutoff
    cond = 1e-5
    U, log_pdet = _psd_pinv_decomposed_log_pdet(cov, cond)
    pinv = np.dot(U, U.T)
    _, log_pdet_pinv = _psd_pinv_decomposed_log_pdet(pinv, cond)

    # Check that the log pseudo-determinant agrees with the sum
    # of the logs of all but the smallest eigenvalue
    assert_allclose(log_pdet, np.sum(np.log(s[:-1])))
    # Check that the pseudo-determinant of the pseudo-inverse
    # agrees with 1 / pseudo-determinant
    assert_allclose(-log_pdet, log_pdet_pinv)


def test_R_values():
    # Compare the multivariate pdf with some values precomputed
    # in R version 3.0.1 (2013-05-16) on Mac OS X 10.6.

    # The values below were generated by the following R-script:
    # > library(mnormt)
    # > x <- seq(0, 2, length=5)
    # > y <- 3*x - 2
    # > z <- x + cos(y)
    # > mu <- c(1, 3, 2)
    # > Sigma <- matrix(c(1,2,0,2,5,0.5,0,0.5,3), 3, 3)
    # > r_pdf <- dmnorm(cbind(x,y,z), mu, Sigma)
    r_pdf = np.array([0.0002214706, 0.0013819953, 0.0049138692,
                      0.0103803050, 0.0140250800])

    x = np.linspace(0, 2, 5)
    y = 3 * x - 2
    z = x + np.cos(y)
    r = np.array([x, y, z]).T

    mean = np.array([1, 3, 2], 'd')
    cov = np.array([[1, 2, 0], [2, 5, .5], [0, .5, 3]], 'd')

    pdf = multivariate_normal.pdf(r, mean, cov)
    assert_allclose(pdf, r_pdf, atol=1e-10)


if __name__ == "__main__":
    run_module_suite()
