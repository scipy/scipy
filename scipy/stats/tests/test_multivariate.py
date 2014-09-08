"""
Test functions for multivariate normal distributions.

"""
from __future__ import division, print_function, absolute_import

from numpy.testing import (
    assert_allclose,
    assert_almost_equal,
    assert_array_almost_equal,
    assert_equal,
    assert_raises,
    run_module_suite,
)

import numpy
import numpy as np

import scipy.linalg
from scipy.stats._multivariate import _PSD, _lnB
from scipy.stats import multivariate_normal
from scipy.stats import dirichlet, beta
from scipy.stats import norm

from scipy.integrate import romb


def test_input_shape():
    mu = np.arange(3)
    cov = np.identity(2)
    assert_raises(ValueError, multivariate_normal.pdf, (0, 1), mu, cov)
    assert_raises(ValueError, multivariate_normal.pdf, (0, 1, 2), mu, cov)


def test_scalar_values():
    np.random.seed(1234)

    # When evaluated on scalar data, the pdf should return a scalar
    x, mean, cov = 1.5, 1.7, 2.5
    pdf = multivariate_normal.pdf(x, mean, cov)
    assert_equal(pdf.ndim, 0)

    # When evaluated on a single vector, the pdf should return a scalar
    x = np.random.randn(5)
    mean = np.random.randn(5)
    cov = np.abs(np.random.randn(5))  # Diagonal values for cov. matrix
    pdf = multivariate_normal.pdf(x, mean, cov)
    assert_equal(pdf.ndim, 0)


def test_logpdf():
    # Check that the log of the pdf is in fact the logpdf
    np.random.seed(1234)
    x = np.random.randn(5)
    mean = np.random.randn(5)
    cov = np.abs(np.random.randn(5))
    d1 = multivariate_normal.logpdf(x, mean, cov)
    d2 = multivariate_normal.pdf(x, mean, cov)
    assert_allclose(d1, np.log(d2))


def test_rank():
    # Check that the rank is detected correctly.
    np.random.seed(1234)
    n = 4
    mean = np.random.randn(n)
    for expected_rank in range(1, n + 1):
        s = np.random.randn(n, expected_rank)
        cov = np.dot(s, s.T)
        distn = multivariate_normal(mean, cov, allow_singular=True)
        assert_equal(distn.cov_info.rank, expected_rank)


def _sample_orthonormal_matrix(n):
    M = np.random.randn(n, n)
    u, s, v = scipy.linalg.svd(M)
    return u


def test_degenerate_distributions():
    for n in range(1, 5):
        x = np.random.randn(n)
        for k in range(1, n + 1):
            # Sample a small covariance matrix.
            s = np.random.randn(k, k)
            cov_kk = np.dot(s, s.T)

            # Embed the small covariance matrix into a larger low rank matrix.
            cov_nn = np.zeros((n, n))
            cov_nn[:k, :k] = cov_kk

            # Define a rotation of the larger low rank matrix.
            u = _sample_orthonormal_matrix(n)
            cov_rr = np.dot(u, np.dot(cov_nn, u.T))
            y = np.dot(u, x)

            # Check some identities.
            distn_kk = multivariate_normal(np.zeros(k), cov_kk,
                                           allow_singular=True)
            distn_nn = multivariate_normal(np.zeros(n), cov_nn,
                                           allow_singular=True)
            distn_rr = multivariate_normal(np.zeros(n), cov_rr,
                                           allow_singular=True)
            assert_equal(distn_kk.cov_info.rank, k)
            assert_equal(distn_nn.cov_info.rank, k)
            assert_equal(distn_rr.cov_info.rank, k)
            pdf_kk = distn_kk.pdf(x[:k])
            pdf_nn = distn_nn.pdf(x)
            pdf_rr = distn_rr.pdf(y)
            assert_allclose(pdf_kk, pdf_nn)
            assert_allclose(pdf_kk, pdf_rr)
            logpdf_kk = distn_kk.logpdf(x[:k])
            logpdf_nn = distn_nn.logpdf(x)
            logpdf_rr = distn_rr.logpdf(y)
            assert_allclose(logpdf_kk, logpdf_nn)
            assert_allclose(logpdf_kk, logpdf_rr)


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
    # assert_allclose(np.linalg.slogdet(cov[:npos, :npos]),
    #                 (1, large_total_log))

    # Check the pseudo-determinant.
    psd = _PSD(cov)
    assert_allclose(psd.log_pdet, large_total_log)


def test_broadcasting():
    np.random.seed(1234)
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
    scale = cov**0.5
    d1 = norm.pdf(x, mean, scale)
    d2 = multivariate_normal.pdf(x, mean, cov)
    assert_allclose(d1, d2)


def test_marginalization():
    # Integrating out one of the variables of a 2D Gaussian should
    # yield a 1D Gaussian
    mean = np.array([2.5, 3.5])
    cov = np.array([[.5, 0.2], [0.2, .6]])
    n = 2 ** 8 + 1  # Number of samples
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
    gauss_x = norm.pdf(v, loc=mean[0], scale=cov[0, 0] ** 0.5)
    gauss_y = norm.pdf(v, loc=mean[1], scale=cov[1, 1] ** 0.5)
    assert_allclose(margin_x, gauss_x, rtol=1e-2, atol=1e-2)
    assert_allclose(margin_y, gauss_y, rtol=1e-2, atol=1e-2)


def test_frozen():
    # The frozen distribution should agree with the regular one
    np.random.seed(1234)
    x = np.random.randn(5)
    mean = np.random.randn(5)
    cov = np.abs(np.random.randn(5))
    norm_frozen = multivariate_normal(mean, cov)
    assert_allclose(norm_frozen.pdf(x), multivariate_normal.pdf(x, mean, cov))
    assert_allclose(norm_frozen.logpdf(x),
                    multivariate_normal.logpdf(x, mean, cov))


def test_pseudodet_pinv():
    # Make sure that pseudo-inverse and pseudo-det agree on cutoff

    # Assemble random covariance matrix with large and small eigenvalues
    np.random.seed(1234)
    n = 7
    x = np.random.randn(n, n)
    cov = np.dot(x, x.T)
    s, u = scipy.linalg.eigh(cov)
    s = 0.5 * np.ones(n)
    s[0] = 1.0
    s[-1] = 1e-7
    cov = np.dot(u, np.dot(np.diag(s), u.T))

    # Set cond so that the lowest eigenvalue is below the cutoff
    cond = 1e-5
    psd = _PSD(cov, cond=cond)
    psd_pinv = _PSD(psd.pinv, cond=cond)

    # Check that the log pseudo-determinant agrees with the sum
    # of the logs of all but the smallest eigenvalue
    assert_allclose(psd.log_pdet, np.sum(np.log(s[:-1])))
    # Check that the pseudo-determinant of the pseudo-inverse
    # agrees with 1 / pseudo-determinant
    assert_allclose(-psd.log_pdet, psd_pinv.log_pdet)


def test_exception_nonsquare_cov():
    cov = [[1, 2, 3], [4, 5, 6]]
    assert_raises(ValueError, _PSD, cov)


def test_exception_nonfinite_cov():
    cov_nan = [[1, 0], [0, np.nan]]
    assert_raises(ValueError, _PSD, cov_nan)
    cov_inf = [[1, 0], [0, np.inf]]
    assert_raises(ValueError, _PSD, cov_inf)


def test_exception_non_psd_cov():
    cov = [[1, 0], [0, -1]]
    assert_raises(ValueError, _PSD, cov)


def test_exception_singular_cov():
    np.random.seed(1234)
    x = np.random.randn(5)
    mean = np.random.randn(5)
    cov = np.ones((5, 5))
    e = np.linalg.LinAlgError
    assert_raises(e, multivariate_normal, mean, cov)
    assert_raises(e, multivariate_normal.pdf, x, mean, cov)
    assert_raises(e, multivariate_normal.logpdf, x, mean, cov)


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


def test_multivariate_normal_rvs_zero_covariance():
    mean = np.zeros(2)
    covariance = np.zeros((2, 2))
    model = multivariate_normal(mean, covariance, allow_singular=True)
    sample = model.rvs()
    assert_equal(sample, [0, 0])


def test_rvs_shape():
    # Check that rvs parses the mean and covariance correctly, and returns
    # an array of the right shape
    N = 300
    d = 4
    sample = multivariate_normal.rvs(mean=np.zeros(d), cov=1, size=N)
    assert_equal(sample.shape, (N, d))

    sample = multivariate_normal.rvs(mean=None,
                                     cov=np.array([[2, .1], [.1, 1]]),
                                     size=N)
    assert_equal(sample.shape, (N, 2))

    u = multivariate_normal(mean=0, cov=1)
    sample = u.rvs(N)
    assert_equal(sample.shape, (N, ))


def test_large_sample():
    # Generate large sample and compare sample mean and sample covariance
    # with mean and covariance matrix.

    np.random.seed(2846)

    n = 3
    mean = np.random.randn(n)
    M = np.random.randn(n, n)
    cov = np.dot(M, M.T)
    size = 5000

    sample = multivariate_normal.rvs(mean, cov, size)

    assert_allclose(numpy.cov(sample.T), cov, rtol=1e-1)
    assert_allclose(sample.mean(0), mean, rtol=1e-1)


def test_entropy():
    np.random.seed(2846)

    n = 3
    mean = np.random.randn(n)
    M = np.random.randn(n, n)
    cov = np.dot(M, M.T)

    rv = multivariate_normal(mean, cov)

    # Check that frozen distribution agrees with entropy function
    assert_almost_equal(rv.entropy(), multivariate_normal.entropy(mean, cov))
    # Compare entropy with manually computed expression involving
    # the sum of the logs of the eigenvalues of the covariance matrix
    eigs = np.linalg.eig(cov)[0]
    desired = 1 / 2 * (n * (np.log(2 * np.pi) + 1) + np.sum(np.log(eigs)))
    assert_almost_equal(desired, rv.entropy())


def test_lnB():
    alpha = np.array([1, 1, 1])
    desired = .5  # e^lnB = 1/2 for [1, 1, 1]

    assert_almost_equal(np.exp(_lnB(alpha)), desired)


def test_frozen_dirichlet():
    np.random.seed(2846)

    n = np.random.randint(1, 32)
    alpha = np.random.uniform(10e-10, 100, n)

    d = dirichlet(alpha)

    assert_equal(d.var(), dirichlet.var(alpha))
    assert_equal(d.mean(), dirichlet.mean(alpha))
    assert_equal(d.entropy(), dirichlet.entropy(alpha))
    num_tests = 10
    for i in range(num_tests):
        x = np.random.uniform(10e-10, 100, n)
        x /= np.sum(x)
        assert_equal(d.pdf(x[:-1]), dirichlet.pdf(x[:-1], alpha))
        assert_equal(d.logpdf(x[:-1]), dirichlet.logpdf(x[:-1], alpha))


def test_simple_values():
    alpha = np.array([1, 1])
    d = dirichlet(alpha)

    assert_almost_equal(d.mean(), 0.5)
    assert_almost_equal(d.var(), 1. / 12.)

    b = beta(1, 1)
    assert_almost_equal(d.mean(), b.mean())
    assert_almost_equal(d.var(), b.var())


def test_K_and_K_minus_1_calls_equal():
    # Test that calls with K and K-1 entries yield the same results.

    np.random.seed(2846)

    n = np.random.randint(1, 32)
    alpha = np.random.uniform(10e-10, 100, n)

    d = dirichlet(alpha)
    num_tests = 10
    for i in range(num_tests):
        x = np.random.uniform(10e-10, 100, n)
        x /= np.sum(x)
        assert_almost_equal(d.pdf(x[:-1]), d.pdf(x))


def test_multiple_entry_calls():
    # Test that calls with multiple x vectors as matrix work

    np.random.seed(2846)

    n = np.random.randint(1, 32)
    alpha = np.random.uniform(10e-10, 100, n)
    d = dirichlet(alpha)

    num_tests = 10
    num_multiple = 5
    xm = None
    for i in range(num_tests):
        for m in range(num_multiple):
            x = np.random.uniform(10e-10, 100, n)
            x /= np.sum(x)
            if xm is not None:
                xm = np.vstack((xm, x))
            else:
                xm = x
        rm = d.pdf(xm.T)
        rs = None
        for xs in xm:
            r = d.pdf(xs)
            if rs is not None:
                rs = np.append(rs, r)
            else:
                rs = r
        assert_array_almost_equal(rm, rs)


def test_2D_dirichlet_is_beta():
    np.random.seed(2846)

    alpha = np.random.uniform(10e-10, 100, 2)
    d = dirichlet(alpha)
    b = beta(alpha[0], alpha[1])

    num_tests = 10
    for i in range(num_tests):
        x = np.random.uniform(10e-10, 100, 2)
        x /= np.sum(x)
        assert_almost_equal(b.pdf(x), d.pdf([x]))

    assert_almost_equal(b.mean(), d.mean()[0])
    assert_almost_equal(b.var(), d.var()[0])


def test_dimensions_mismatch():
    # Regression test for GH #3493. Check that setting up a PDF with a mean of
    # length M and a covariance matrix of size (N, N), where M != N, raises a
    # ValueError with an informative error message.

    mu = np.array([0.0, 0.0])
    sigma = np.array([[1.0]])

    assert_raises(ValueError, multivariate_normal, mu, sigma)

    # A simple check that the right error message was passed along. Checking
    # that the entire message is there, word for word, would be somewhat
    # fragile, so we just check for the leading part.
    try:
        multivariate_normal(mu, sigma)
    except ValueError as e:
        msg = "Dimension mismatch"
        assert_equal(str(e)[:len(msg)], msg)


if __name__ == "__main__":
    run_module_suite()
