"""
Test functions for multivariate normal distributions.

"""
from __future__ import division, print_function, absolute_import

import pickle

from numpy.testing import (assert_allclose, assert_almost_equal,
                           assert_array_almost_equal, assert_equal,
                           assert_array_less, assert_raises, assert_,
                           run_module_suite, TestCase)

from test_continuous_basic import check_distribution_rvs

import numpy
import numpy as np

import scipy.linalg
from scipy.stats._multivariate import _PSD, _lnB
from scipy.stats import multivariate_normal
from scipy.stats import matrix_normal
from scipy.stats import special_ortho_group, ortho_group
from scipy.stats import random_correlation
from scipy.stats import dirichlet, beta
from scipy.stats import wishart, invwishart, chi2, invgamma
from scipy.stats import norm
from scipy.stats import ks_2samp

from scipy.integrate import romb

from common_tests import check_random_state_property


class TestMultivariateNormal(TestCase):
    def test_input_shape(self):
        mu = np.arange(3)
        cov = np.identity(2)
        assert_raises(ValueError, multivariate_normal.pdf, (0, 1), mu, cov)
        assert_raises(ValueError, multivariate_normal.pdf, (0, 1, 2), mu, cov)

    def test_scalar_values(self):
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

    def test_logpdf(self):
        # Check that the log of the pdf is in fact the logpdf
        np.random.seed(1234)
        x = np.random.randn(5)
        mean = np.random.randn(5)
        cov = np.abs(np.random.randn(5))
        d1 = multivariate_normal.logpdf(x, mean, cov)
        d2 = multivariate_normal.pdf(x, mean, cov)
        assert_allclose(d1, np.log(d2))

    def test_logpdf_default_values(self):
        # Check that the log of the pdf is in fact the logpdf
        # with default parameters Mean=None and cov = 1
        np.random.seed(1234)
        x = np.random.randn(5)
        d1 = multivariate_normal.logpdf(x)
        d2 = multivariate_normal.pdf(x)
        # check whether default values are being used
        d3 = multivariate_normal.logpdf(x, None, 1)
        d4 = multivariate_normal.pdf(x, None, 1)
        assert_allclose(d1, np.log(d2))
        assert_allclose(d3, np.log(d4))

    def test_rank(self):
        # Check that the rank is detected correctly.
        np.random.seed(1234)
        n = 4
        mean = np.random.randn(n)
        for expected_rank in range(1, n + 1):
            s = np.random.randn(n, expected_rank)
            cov = np.dot(s, s.T)
            distn = multivariate_normal(mean, cov, allow_singular=True)
            assert_equal(distn.cov_info.rank, expected_rank)

    def test_degenerate_distributions(self):

        def _sample_orthonormal_matrix(n):
            M = np.random.randn(n, n)
            u, s, v = scipy.linalg.svd(M)
            return u

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

    def test_large_pseudo_determinant(self):
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
        assert_allclose(np.linalg.slogdet(cov[:npos, :npos]),
                        (1, large_total_log))

        # Check the pseudo-determinant.
        psd = _PSD(cov)
        assert_allclose(psd.log_pdet, large_total_log)

    def test_broadcasting(self):
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

    def test_normal_1D(self):
        # The probability density function for a 1D normal variable should
        # agree with the standard normal distribution in scipy.stats.distributions
        x = np.linspace(0, 2, 10)
        mean, cov = 1.2, 0.9
        scale = cov**0.5
        d1 = norm.pdf(x, mean, scale)
        d2 = multivariate_normal.pdf(x, mean, cov)
        assert_allclose(d1, d2)

    def test_marginalization(self):
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

    def test_frozen(self):
        # The frozen distribution should agree with the regular one
        np.random.seed(1234)
        x = np.random.randn(5)
        mean = np.random.randn(5)
        cov = np.abs(np.random.randn(5))
        norm_frozen = multivariate_normal(mean, cov)
        assert_allclose(norm_frozen.pdf(x), multivariate_normal.pdf(x, mean, cov))
        assert_allclose(norm_frozen.logpdf(x),
                        multivariate_normal.logpdf(x, mean, cov))

    def test_pseudodet_pinv(self):
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

    def test_exception_nonsquare_cov(self):
        cov = [[1, 2, 3], [4, 5, 6]]
        assert_raises(ValueError, _PSD, cov)

    def test_exception_nonfinite_cov(self):
        cov_nan = [[1, 0], [0, np.nan]]
        assert_raises(ValueError, _PSD, cov_nan)
        cov_inf = [[1, 0], [0, np.inf]]
        assert_raises(ValueError, _PSD, cov_inf)

    def test_exception_non_psd_cov(self):
        cov = [[1, 0], [0, -1]]
        assert_raises(ValueError, _PSD, cov)

    def test_exception_singular_cov(self):
        np.random.seed(1234)
        x = np.random.randn(5)
        mean = np.random.randn(5)
        cov = np.ones((5, 5))
        e = np.linalg.LinAlgError
        assert_raises(e, multivariate_normal, mean, cov)
        assert_raises(e, multivariate_normal.pdf, x, mean, cov)
        assert_raises(e, multivariate_normal.logpdf, x, mean, cov)

    def test_R_values(self):
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

    def test_multivariate_normal_rvs_zero_covariance(self):
        mean = np.zeros(2)
        covariance = np.zeros((2, 2))
        model = multivariate_normal(mean, covariance, allow_singular=True)
        sample = model.rvs()
        assert_equal(sample, [0, 0])

    def test_rvs_shape(self):
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

    def test_large_sample(self):
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

    def test_entropy(self):
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

    def test_lnB(self):
        alpha = np.array([1, 1, 1])
        desired = .5  # e^lnB = 1/2 for [1, 1, 1]

        assert_almost_equal(np.exp(_lnB(alpha)), desired)

class TestMatrixNormal(TestCase):

    def test_bad_input(self):
        # Check that bad inputs raise errors
        num_rows = 4
        num_cols = 3
        M = 0.3 * np.ones((num_rows,num_cols))
        U = 0.5 * np.identity(num_rows) + 0.5 * np.ones((num_rows, num_rows))
        V = 0.7 * np.identity(num_cols) + 0.3 * np.ones((num_cols, num_cols))

        # Incorrect dimensions
        assert_raises(ValueError, matrix_normal, np.zeros((5,4,3)))
        assert_raises(ValueError, matrix_normal, M, np.zeros(10), V)
        assert_raises(ValueError, matrix_normal, M, U, np.zeros(10))
        assert_raises(ValueError, matrix_normal, M, U, U)
        assert_raises(ValueError, matrix_normal, M, V, V)
        assert_raises(ValueError, matrix_normal, M.T, U, V)

        # Singular covariance
        e = np.linalg.LinAlgError
        assert_raises(e, matrix_normal, M, U, np.ones((num_cols, num_cols)))
        assert_raises(e, matrix_normal, M, np.ones((num_rows, num_rows)), V)

    def test_default_inputs(self):
        # Check that default argument handling works
        num_rows = 4
        num_cols = 3
        M = 0.3 * np.ones((num_rows,num_cols))
        U = 0.5 * np.identity(num_rows) + 0.5 * np.ones((num_rows, num_rows))
        V = 0.7 * np.identity(num_cols) + 0.3 * np.ones((num_cols, num_cols))
        Z = np.zeros((num_rows, num_cols))
        Zr = np.zeros((num_rows, 1))
        Zc = np.zeros((1, num_cols))
        Ir = np.identity(num_rows)
        Ic = np.identity(num_cols)
        I1 = np.identity(1)

        assert_equal(matrix_normal.rvs(mean=M, rowcov=U, colcov=V).shape,
                     (num_rows, num_cols))
        assert_equal(matrix_normal.rvs(mean=M).shape,
                     (num_rows, num_cols))
        assert_equal(matrix_normal.rvs(rowcov=U).shape,
                     (num_rows, 1))
        assert_equal(matrix_normal.rvs(colcov=V).shape,
                     (1, num_cols))
        assert_equal(matrix_normal.rvs(mean=M, colcov=V).shape,
                     (num_rows, num_cols))
        assert_equal(matrix_normal.rvs(mean=M, rowcov=U).shape,
                     (num_rows, num_cols))
        assert_equal(matrix_normal.rvs(rowcov=U, colcov=V).shape,
                     (num_rows, num_cols))

        assert_equal(matrix_normal(mean=M).rowcov, Ir)
        assert_equal(matrix_normal(mean=M).colcov, Ic)
        assert_equal(matrix_normal(rowcov=U).mean, Zr)
        assert_equal(matrix_normal(rowcov=U).colcov, I1)
        assert_equal(matrix_normal(colcov=V).mean, Zc)
        assert_equal(matrix_normal(colcov=V).rowcov, I1)
        assert_equal(matrix_normal(mean=M, rowcov=U).colcov, Ic)
        assert_equal(matrix_normal(mean=M, colcov=V).rowcov, Ir)
        assert_equal(matrix_normal(rowcov=U, colcov=V).mean, Z)

    def test_covariance_expansion(self):
        # Check that covariance can be specified with scalar or vector
        num_rows = 4
        num_cols = 3
        M = 0.3 * np.ones((num_rows,num_cols))
        Uv = 0.2*np.ones(num_rows)
        Us = 0.2
        Vv = 0.1*np.ones(num_cols)
        Vs = 0.1

        Ir = np.identity(num_rows)
        Ic = np.identity(num_cols)

        assert_equal(matrix_normal(mean=M, rowcov=Uv, colcov=Vv).rowcov,
                     0.2*Ir)
        assert_equal(matrix_normal(mean=M, rowcov=Uv, colcov=Vv).colcov,
                     0.1*Ic)
        assert_equal(matrix_normal(mean=M, rowcov=Us, colcov=Vs).rowcov,
                     0.2*Ir)
        assert_equal(matrix_normal(mean=M, rowcov=Us, colcov=Vs).colcov,
                     0.1*Ic)

    def test_frozen_matrix_normal(self):
        for i in range(1,5):
            for j in range(1,5):
                M = 0.3 * np.ones((i,j))
                U = 0.5 * np.identity(i) + 0.5 * np.ones((i,i))
                V = 0.7 * np.identity(j) + 0.3 * np.ones((j,j))

                frozen = matrix_normal(mean=M, rowcov=U, colcov=V)

                rvs1 = frozen.rvs(random_state=1234)
                rvs2 = matrix_normal.rvs(mean=M, rowcov=U, colcov=V,
                                         random_state=1234)
                assert_equal(rvs1, rvs2)

                X = frozen.rvs(random_state=1234)

                pdf1 = frozen.pdf(X)
                pdf2 = matrix_normal.pdf(X, mean=M, rowcov=U, colcov=V)
                assert_equal(pdf1, pdf2)

                logpdf1 = frozen.logpdf(X)
                logpdf2 = matrix_normal.logpdf(X, mean=M, rowcov=U, colcov=V)
                assert_equal(logpdf1, logpdf2)

    def test_matches_multivariate(self):
        # Check that the pdfs match those obtained by vectorising and
        # treating as a multivariate normal.
        for i in range(1,5):
            for j in range(1,5):
                M = 0.3 * np.ones((i,j))
                U = 0.5 * np.identity(i) + 0.5 * np.ones((i,i))
                V = 0.7 * np.identity(j) + 0.3 * np.ones((j,j))

                frozen = matrix_normal(mean=M, rowcov=U, colcov=V)
                X = frozen.rvs(random_state=1234)
                pdf1 = frozen.pdf(X)
                logpdf1 = frozen.logpdf(X)

                vecX = X.T.flatten()
                vecM = M.T.flatten()
                cov = np.kron(V,U)
                pdf2 = multivariate_normal.pdf(vecX, mean=vecM, cov=cov)
                logpdf2 = multivariate_normal.logpdf(vecX, mean=vecM, cov=cov)

                assert_allclose(pdf1, pdf2, rtol=1E-10)
                assert_allclose(logpdf1, logpdf2, rtol=1E-10)

    def test_array_input(self):
        # Check array of inputs has the same output as the separate entries.
        num_rows = 4
        num_cols = 3
        M = 0.3 * np.ones((num_rows,num_cols))
        U = 0.5 * np.identity(num_rows) + 0.5 * np.ones((num_rows, num_rows))
        V = 0.7 * np.identity(num_cols) + 0.3 * np.ones((num_cols, num_cols))
        N = 10

        frozen = matrix_normal(mean=M, rowcov=U, colcov=V)
        X1 = frozen.rvs(size=N, random_state=1234)
        X2 = frozen.rvs(size=N, random_state=4321)
        X = np.concatenate((X1[np.newaxis,:,:,:],X2[np.newaxis,:,:,:]), axis=0)
        assert_equal(X.shape, (2, N, num_rows, num_cols))

        array_logpdf = frozen.logpdf(X)
        assert_equal(array_logpdf.shape, (2, N))
        for i in range(2):
            for j in range(N):
                separate_logpdf = matrix_normal.logpdf(X[i,j], mean=M,
                                                       rowcov=U, colcov=V)
                assert_allclose(separate_logpdf, array_logpdf[i,j], 1E-10)

    def test_moments(self):
        # Check that the sample moments match the parameters
        num_rows = 4
        num_cols = 3
        M = 0.3 * np.ones((num_rows,num_cols))
        U = 0.5 * np.identity(num_rows) + 0.5 * np.ones((num_rows, num_rows))
        V = 0.7 * np.identity(num_cols) + 0.3 * np.ones((num_cols, num_cols))
        N = 1000

        frozen = matrix_normal(mean=M, rowcov=U, colcov=V)
        X = frozen.rvs(size=N, random_state=1234)

        sample_mean = np.mean(X,axis=0)
        assert_allclose(sample_mean, M, atol=0.1)

        sample_colcov = np.cov(X.reshape(N*num_rows,num_cols).T)
        assert_allclose(sample_colcov, V, atol=0.1)

        sample_rowcov = np.cov(np.swapaxes(X,1,2).reshape(
                                                        N*num_cols,num_rows).T)
        assert_allclose(sample_rowcov, U, atol=0.1)

class TestDirichlet(TestCase):

    def test_frozen_dirichlet(self):
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

    def test_numpy_rvs_shape_compatibility(self):
        np.random.seed(2846)
        alpha = np.array([1.0, 2.0, 3.0])
        x = np.random.dirichlet(alpha, size=7)
        assert_equal(x.shape, (7, 3))
        assert_raises(ValueError, dirichlet.pdf, x, alpha)
        assert_raises(ValueError, dirichlet.logpdf, x, alpha)
        dirichlet.pdf(x.T, alpha)
        dirichlet.pdf(x.T[:-1], alpha)
        dirichlet.logpdf(x.T, alpha)
        dirichlet.logpdf(x.T[:-1], alpha)

    def test_alpha_with_zeros(self):
        np.random.seed(2846)
        alpha = [1.0, 0.0, 3.0]
        x = np.random.dirichlet(alpha, size=7).T
        assert_raises(ValueError, dirichlet.pdf, x, alpha)
        assert_raises(ValueError, dirichlet.logpdf, x, alpha)

    def test_alpha_with_negative_entries(self):
        np.random.seed(2846)
        alpha = [1.0, -2.0, 3.0]
        x = np.random.dirichlet(alpha, size=7).T
        assert_raises(ValueError, dirichlet.pdf, x, alpha)
        assert_raises(ValueError, dirichlet.logpdf, x, alpha)

    def test_data_with_zeros(self):
        alpha = np.array([1.0, 2.0, 3.0, 4.0])
        x = np.array([0.1, 0.0, 0.2, 0.7])
        assert_raises(ValueError, dirichlet.pdf, x, alpha)
        assert_raises(ValueError, dirichlet.logpdf, x, alpha)

    def test_data_with_negative_entries(self):
        alpha = np.array([1.0, 2.0, 3.0, 4.0])
        x = np.array([0.1, -0.1, 0.3, 0.7])
        assert_raises(ValueError, dirichlet.pdf, x, alpha)
        assert_raises(ValueError, dirichlet.logpdf, x, alpha)

    def test_data_with_too_large_entries(self):
        alpha = np.array([1.0, 2.0, 3.0, 4.0])
        x = np.array([0.1, 1.1, 0.3, 0.7])
        assert_raises(ValueError, dirichlet.pdf, x, alpha)
        assert_raises(ValueError, dirichlet.logpdf, x, alpha)

    def test_data_too_deep_c(self):
        alpha = np.array([1.0, 2.0, 3.0])
        x = np.ones((2, 7, 7)) / 14
        assert_raises(ValueError, dirichlet.pdf, x, alpha)
        assert_raises(ValueError, dirichlet.logpdf, x, alpha)

    def test_alpha_too_deep(self):
        alpha = np.array([[1.0, 2.0], [3.0, 4.0]])
        x = np.ones((2, 2, 7)) / 4
        assert_raises(ValueError, dirichlet.pdf, x, alpha)
        assert_raises(ValueError, dirichlet.logpdf, x, alpha)

    def test_alpha_correct_depth(self):
        alpha = np.array([1.0, 2.0, 3.0])
        x = np.ones((3, 7)) / 3
        dirichlet.pdf(x, alpha)
        dirichlet.logpdf(x, alpha)

    def test_non_simplex_data(self):
        alpha = np.array([1.0, 2.0, 3.0])
        x = np.ones((3, 7)) / 2
        assert_raises(ValueError, dirichlet.pdf, x, alpha)
        assert_raises(ValueError, dirichlet.logpdf, x, alpha)

    def test_data_vector_too_short(self):
        alpha = np.array([1.0, 2.0, 3.0, 4.0])
        x = np.ones((2, 7)) / 2
        assert_raises(ValueError, dirichlet.pdf, x, alpha)
        assert_raises(ValueError, dirichlet.logpdf, x, alpha)

    def test_data_vector_too_long(self):
        alpha = np.array([1.0, 2.0, 3.0, 4.0])
        x = np.ones((5, 7)) / 5
        assert_raises(ValueError, dirichlet.pdf, x, alpha)
        assert_raises(ValueError, dirichlet.logpdf, x, alpha)

    def test_simple_values(self):
        alpha = np.array([1, 1])
        d = dirichlet(alpha)

        assert_almost_equal(d.mean(), 0.5)
        assert_almost_equal(d.var(), 1. / 12.)

        b = beta(1, 1)
        assert_almost_equal(d.mean(), b.mean())
        assert_almost_equal(d.var(), b.var())

    def test_K_and_K_minus_1_calls_equal(self):
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

    def test_multiple_entry_calls(self):
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

    def test_2D_dirichlet_is_beta(self):
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


def test_multivariate_normal_dimensions_mismatch():
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


class TestWishart(TestCase):
    def test_scale_dimensions(self):
        # Test that we can call the Wishart with various scale dimensions

        # Test case: dim=1, scale=1
        true_scale = np.array(1, ndmin=2)
        scales = [
            1,                    # scalar
            [1],                  # iterable
            np.array(1),          # 0-dim
            np.r_[1],             # 1-dim
            np.array(1, ndmin=2)  # 2-dim
        ]
        for scale in scales:
            w = wishart(1, scale)
            assert_equal(w.scale, true_scale)
            assert_equal(w.scale.shape, true_scale.shape)

        # Test case: dim=2, scale=[[1,0]
        #                          [0,2]
        true_scale = np.array([[1,0],
                               [0,2]])
        scales = [
            [1,2],             # iterable
            np.r_[1,2],        # 1-dim
            np.array([[1,0],   # 2-dim
                      [0,2]])
        ]
        for scale in scales:
            w = wishart(2, scale)
            assert_equal(w.scale, true_scale)
            assert_equal(w.scale.shape, true_scale.shape)

        # We cannot call with a df < dim
        assert_raises(ValueError, wishart, 1, np.eye(2))

        # We cannot call with a 3-dimension array
        scale = np.array(1, ndmin=3)
        assert_raises(ValueError, wishart, 1, scale)

    def test_quantile_dimensions(self):
        # Test that we can call the Wishart rvs with various quantile dimensions

        # If dim == 1, consider x.shape = [1,1,1]
        X = [
            1,                      # scalar
            [1],                    # iterable
            np.array(1),            # 0-dim
            np.r_[1],               # 1-dim
            np.array(1, ndmin=2),   # 2-dim
            np.array([1], ndmin=3)  # 3-dim
        ]

        w = wishart(1,1)
        density = w.pdf(np.array(1, ndmin=3))
        for x in X:
            assert_equal(w.pdf(x), density)

        # If dim == 1, consider x.shape = [1,1,*]
        X = [
            [1,2,3],                     # iterable
            np.r_[1,2,3],                # 1-dim
            np.array([1,2,3], ndmin=3)   # 3-dim
        ]

        w = wishart(1,1)
        density = w.pdf(np.array([1,2,3], ndmin=3))
        for x in X:
            assert_equal(w.pdf(x), density)

        # If dim == 2, consider x.shape = [2,2,1]
        # where x[:,:,*] = np.eye(1)*2
        X = [
            2,                    # scalar
            [2,2],                # iterable
            np.array(2),          # 0-dim
            np.r_[2,2],           # 1-dim
            np.array([[2,0],
                      [0,2]]),    # 2-dim
            np.array([[2,0],
                      [0,2]])[:,:,np.newaxis]  # 3-dim
        ]

        w = wishart(2,np.eye(2))
        density = w.pdf(np.array([[2,0],
                                  [0,2]])[:,:,np.newaxis])
        for x in X:
            assert_equal(w.pdf(x), density)

    def test_frozen(self):
        # Test that the frozen and non-frozen Wishart gives the same answers

        # Construct an arbitrary positive definite scale matrix
        dim = 4
        scale = np.diag(np.arange(dim)+1)
        scale[np.tril_indices(dim, k=-1)] = np.arange(dim * (dim-1) // 2)
        scale = np.dot(scale.T, scale)

        # Construct a collection of positive definite matrices to test the PDF
        X = []
        for i in range(5):
            x = np.diag(np.arange(dim)+(i+1)**2)
            x[np.tril_indices(dim, k=-1)] = np.arange(dim * (dim-1) // 2)
            x = np.dot(x.T, x)
            X.append(x)
        X = np.array(X).T

        # Construct a 1D and 2D set of parameters
        parameters = [
            (10, 1, np.linspace(0.1, 10, 5)),  # 1D case
            (10, scale, X)
        ]

        for (df, scale, x) in parameters:
            w = wishart(df, scale)
            assert_equal(w.var(), wishart.var(df, scale))
            assert_equal(w.mean(), wishart.mean(df, scale))
            assert_equal(w.mode(), wishart.mode(df, scale))
            assert_equal(w.entropy(), wishart.entropy(df, scale))
            assert_equal(w.pdf(x), wishart.pdf(x, df, scale))

    def test_1D_is_chisquared(self):
        # The 1-dimensional Wishart with an identity scale matrix is just a
        # chi-squared distribution.
        # Test variance, mean, entropy, pdf
        # Kolgomorov-Smirnov test for rvs
        np.random.seed(482974)

        sn = 500
        dim = 1
        scale = np.eye(dim)

        df_range = np.arange(1, 10, 2, dtype=float)
        X = np.linspace(0.1,10,num=10)
        for df in df_range:
            w = wishart(df, scale)
            c = chi2(df)

            # Statistics
            assert_allclose(w.var(), c.var())
            assert_allclose(w.mean(), c.mean())
            assert_allclose(w.entropy(), c.entropy())

            # PDF
            assert_allclose(w.pdf(X), c.pdf(X))

            # rvs
            rvs = w.rvs(size=sn)
            args = (df,)
            alpha = 0.01
            check_distribution_rvs('chi2', args, alpha, rvs)

    def test_is_scaled_chisquared(self):
        # The 2-dimensional Wishart with an arbitrary scale matrix can be
        # transformed to a scaled chi-squared distribution.
        # For :math:`S \sim W_p(V,n)` and :math:`\lambda \in \mathbb{R}^p` we have
        # :math:`\lambda' S \lambda \sim \lambda' V \lambda \times \chi^2(n)`
        np.random.seed(482974)

        sn = 500
        df = 10
        dim = 4
        # Construct an arbitrary positive definite matrix
        scale = np.diag(np.arange(4)+1)
        scale[np.tril_indices(4, k=-1)] = np.arange(6)
        scale = np.dot(scale.T, scale)
        # Use :math:`\lambda = [1, \dots, 1]'`
        lamda = np.ones((dim,1))
        sigma_lamda = lamda.T.dot(scale).dot(lamda).squeeze()
        w = wishart(df, sigma_lamda)
        c = chi2(df, scale=sigma_lamda)

        # Statistics
        assert_allclose(w.var(), c.var())
        assert_allclose(w.mean(), c.mean())
        assert_allclose(w.entropy(), c.entropy())

        # PDF
        X = np.linspace(0.1,10,num=10)
        assert_allclose(w.pdf(X), c.pdf(X))

        # rvs
        rvs = w.rvs(size=sn)
        args = (df,0,sigma_lamda)
        alpha = 0.01
        check_distribution_rvs('chi2', args, alpha, rvs)


class TestInvwishart(TestCase):
    def test_frozen(self):
        # Test that the frozen and non-frozen inverse Wishart gives the same
        # answers

        # Construct an arbitrary positive definite scale matrix
        dim = 4
        scale = np.diag(np.arange(dim)+1)
        scale[np.tril_indices(dim, k=-1)] = np.arange(dim*(dim-1)/2)
        scale = np.dot(scale.T, scale)

        # Construct a collection of positive definite matrices to test the PDF
        X = []
        for i in range(5):
            x = np.diag(np.arange(dim)+(i+1)**2)
            x[np.tril_indices(dim, k=-1)] = np.arange(dim*(dim-1)/2)
            x = np.dot(x.T, x)
            X.append(x)
        X = np.array(X).T

        # Construct a 1D and 2D set of parameters
        parameters = [
            (10, 1, np.linspace(0.1, 10, 5)),  # 1D case
            (10, scale, X)
        ]

        for (df, scale, x) in parameters:
            iw = invwishart(df, scale)
            assert_equal(iw.var(), invwishart.var(df, scale))
            assert_equal(iw.mean(), invwishart.mean(df, scale))
            assert_equal(iw.mode(), invwishart.mode(df, scale))
            assert_allclose(iw.pdf(x), invwishart.pdf(x, df, scale))

    def test_1D_is_invgamma(self):
        # The 1-dimensional inverse Wishart with an identity scale matrix is
        # just an inverse gamma distribution.
        # Test variance, mean, pdf
        # Kolgomorov-Smirnov test for rvs
        np.random.seed(482974)

        sn = 500
        dim = 1
        scale = np.eye(dim)

        df_range = np.arange(5, 20, 2, dtype=float)
        X = np.linspace(0.1,10,num=10)
        for df in df_range:
            iw = invwishart(df, scale)
            ig = invgamma(df/2, scale=1./2)

            # Statistics
            assert_allclose(iw.var(), ig.var())
            assert_allclose(iw.mean(), ig.mean())

            # PDF
            assert_allclose(iw.pdf(X), ig.pdf(X))

            # rvs
            rvs = iw.rvs(size=sn)
            args = (df/2, 0, 1./2)
            alpha = 0.01
            check_distribution_rvs('invgamma', args, alpha, rvs)

    def test_wishart_invwishart_2D_rvs(self):
        dim = 3
        df = 10

        # Construct a simple non-diagonal positive definite matrix
        scale = np.eye(dim)
        scale[0,1] = 0.5
        scale[1,0] = 0.5

        # Construct frozen Wishart and inverse Wishart random variables
        w = wishart(df, scale)
        iw = invwishart(df, scale)

        # Get the generated random variables from a known seed
        np.random.seed(248042)
        w_rvs = wishart.rvs(df, scale)
        np.random.seed(248042)
        frozen_w_rvs = w.rvs()
        np.random.seed(248042)
        iw_rvs = invwishart.rvs(df, scale)
        np.random.seed(248042)
        frozen_iw_rvs = iw.rvs()

        # Manually calculate what it should be, based on the Bartlett (1933)
        # decomposition of a Wishart into D A A' D', where D is the Cholesky
        # factorization of the scale matrix and A is the lower triangular matrix
        # with the square root of chi^2 variates on the diagonal and N(0,1)
        # variates in the lower triangle.
        np.random.seed(248042)
        covariances = np.random.normal(size=3)
        variances = np.r_[
            np.random.chisquare(df),
            np.random.chisquare(df-1),
            np.random.chisquare(df-2),
        ]**0.5

        # Construct the lower-triangular A matrix
        A = np.diag(variances)
        A[np.tril_indices(dim, k=-1)] = covariances

        # Wishart random variate
        D = np.linalg.cholesky(scale)
        DA = D.dot(A)
        manual_w_rvs = np.dot(DA, DA.T)

        # inverse Wishart random variate
        # Supposing that the inverse wishart has scale matrix `scale`, then the
        # random variate is the inverse of a random variate drawn from a Wishart
        # distribution with scale matrix `inv_scale = np.linalg.inv(scale)`
        iD = np.linalg.cholesky(np.linalg.inv(scale))
        iDA = iD.dot(A)
        manual_iw_rvs = np.linalg.inv(np.dot(iDA, iDA.T))

        # Test for equality
        assert_allclose(w_rvs, manual_w_rvs)
        assert_allclose(frozen_w_rvs, manual_w_rvs)
        assert_allclose(iw_rvs, manual_iw_rvs)
        assert_allclose(frozen_iw_rvs, manual_iw_rvs)


class TestSpecialOrthoGroup(TestCase):
    def test_reproducibility(self):
        np.random.seed(514)
        x = special_ortho_group.rvs(3)
        expected = np.array([[0.99394515, -0.04527879, 0.10011432],
                             [-0.04821555, 0.63900322, 0.76769144],
                             [-0.09873351, -0.76787024, 0.63295101]])
        assert_array_almost_equal(x, expected)

        random_state = np.random.RandomState(seed=514)
        x = special_ortho_group.rvs(3, random_state=random_state)
        assert_array_almost_equal(x, expected)

    def test_invalid_dim(self):
        assert_raises(ValueError, special_ortho_group.rvs, None)
        assert_raises(ValueError, special_ortho_group.rvs, (2, 2))
        assert_raises(ValueError, special_ortho_group.rvs, 1)
        assert_raises(ValueError, special_ortho_group.rvs, 2.5)

    def test_frozen_matrix(self):
        dim = 7
        frozen = special_ortho_group(dim)

        rvs1 = frozen.rvs(random_state=1234)
        rvs2 = special_ortho_group.rvs(dim, random_state=1234)

        assert_equal(rvs1, rvs2)

    def test_det_and_ortho(self):
        xs = [special_ortho_group.rvs(dim)
              for dim in range(2,12)
              for i in range(3)]

        # Test that determinants are always +1
        dets = [np.linalg.det(x) for x in xs]
        assert_allclose(dets, [1.]*30, rtol=1e-13)

        # Test that these are orthogonal matrices
        for x in xs:
            assert_array_almost_equal(np.dot(x, x.T),
                                      np.eye(x.shape[0]))

    def test_haar(self):
        # Test that the distribution is constant under rotation
        # Every column should have the same distribution
        # Additionally, the distribution should be invariant under another rotation

        # Generate samples
        dim = 5
        samples = 1000  # Not too many, or the test takes too long
        ks_prob = 0.39  # ...so don't expect much precision
        np.random.seed(514)
        xs = special_ortho_group.rvs(dim, size=samples)

        # Dot a few rows (0, 1, 2) with unit vectors (0, 2, 4, 3),
        #   effectively picking off entries in the matrices of xs.
        #   These projections should all have the same disribution,
        #     establishing rotational invariance. We use the two-sided
        #     KS test to confirm this.
        #   We could instead test that angles between random vectors
        #     are uniformly distributed, but the below is sufficient.
        #   It is not feasible to consider all pairs, so pick a few.
        els = ((0,0), (0,2), (1,4), (2,3))
        #proj = {(er, ec): [x[er][ec] for x in xs] for er, ec in els}
        proj = dict(((er, ec), sorted([x[er][ec] for x in xs])) for er, ec in els)
        pairs = [(e0, e1) for e0 in els for e1 in els if e0 > e1]
        ks_tests = [ks_2samp(proj[p0], proj[p1])[1] for (p0, p1) in pairs]
        assert_array_less([ks_prob]*len(pairs), ks_tests)

class TestOrthoGroup(TestCase):
    def test_reproducibility(self):
        np.random.seed(514)
        x = ortho_group.rvs(3)
        x2 = ortho_group.rvs(3, random_state=514)
        # Note this matrix has det -1, distinguishing O(N) from SO(N)
        expected = np.array([[0.993945, -0.045279, 0.100114],
                             [-0.048216, -0.998469, 0.02711],
                             [-0.098734, 0.031773, 0.994607]])
        assert_array_almost_equal(x, expected)
        assert_array_almost_equal(x2, expected)
        assert_almost_equal(np.linalg.det(x), -1)

    def test_invalid_dim(self):
        assert_raises(ValueError, ortho_group.rvs, None)
        assert_raises(ValueError, ortho_group.rvs, (2, 2))
        assert_raises(ValueError, ortho_group.rvs, 1)
        assert_raises(ValueError, ortho_group.rvs, 2.5)

    def test_det_and_ortho(self):
        xs = [ortho_group.rvs(dim)
              for dim in range(2,12)
              for i in range(3)]

        # Test that determinants are always +1
        dets = [np.fabs(np.linalg.det(x)) for x in xs]
        assert_allclose(dets, [1.]*30, rtol=1e-13)

        # Test that these are orthogonal matrices
        for x in xs:
            assert_array_almost_equal(np.dot(x, x.T),
                                      np.eye(x.shape[0]))

    def test_haar(self):
        # Test that the distribution is constant under rotation
        # Every column should have the same distribution
        # Additionally, the distribution should be invariant under another rotation

        # Generate samples
        dim = 5
        samples = 1000  # Not too many, or the test takes too long
        ks_prob = 0.39  # ...so don't expect much precision
        np.random.seed(518)  # Note that the test is sensitive to seed too
        xs = ortho_group.rvs(dim, size=samples)

        # Dot a few rows (0, 1, 2) with unit vectors (0, 2, 4, 3),
        #   effectively picking off entries in the matrices of xs.
        #   These projections should all have the same disribution,
        #     establishing rotational invariance. We use the two-sided
        #     KS test to confirm this.
        #   We could instead test that angles between random vectors
        #     are uniformly distributed, but the below is sufficient.
        #   It is not feasible to consider all pairs, so pick a few.
        els = ((0,0), (0,2), (1,4), (2,3))
        #proj = {(er, ec): [x[er][ec] for x in xs] for er, ec in els}
        proj = dict(((er, ec), sorted([x[er][ec] for x in xs])) for er, ec in els)
        pairs = [(e0, e1) for e0 in els for e1 in els if e0 > e1]
        ks_tests = [ks_2samp(proj[p0], proj[p1])[1] for (p0, p1) in pairs]
        assert_array_less([ks_prob]*len(pairs), ks_tests)

class TestRandomCorrelation(TestCase):
    def test_reproducibility(self):
        np.random.seed(514)
        eigs = (.5, .8, 1.2, 1.5)
        x = random_correlation.rvs((.5, .8, 1.2, 1.5))
        x2 = random_correlation.rvs((.5, .8, 1.2, 1.5), random_state=514)
        expected = np.array([[1., -0.20387311, 0.18366501, -0.04953711],
                             [-0.20387311, 1., -0.24351129, 0.06703474],
                             [0.18366501, -0.24351129, 1., 0.38530195],
                             [-0.04953711, 0.06703474, 0.38530195, 1.]])
        assert_array_almost_equal(x, expected)
        assert_array_almost_equal(x2, expected)

    def test_invalid_eigs(self):
        assert_raises(ValueError, random_correlation.rvs, None)
        assert_raises(ValueError, random_correlation.rvs, 'test')
        assert_raises(ValueError, random_correlation.rvs, 2.5)
        assert_raises(ValueError, random_correlation.rvs, [2.5])
        assert_raises(ValueError, random_correlation.rvs, [[1,2],[3,4]])
        assert_raises(ValueError, random_correlation.rvs, [2.5, -.5])
        assert_raises(ValueError, random_correlation.rvs, [1, 2, .1])

    def test_definition(self):
        # Test the defintion of a correlation matrix in several dimensions:
        #
        # 1. Det is product of eigenvalues (and positive by construction
        #    in examples)
        # 2. 1's on diagonal
        # 3. Matrix is symmetric

        def norm(i, e):
            return i*e/sum(e)

        np.random.seed(123)

        eigs = [norm(i, np.random.uniform(size=i)) for i in range(2, 6)]
        eigs.append([4,0,0,0])

        ones = [[1.]*len(e) for e in eigs]
        xs = [random_correlation.rvs(e) for e in eigs]

        # Test that determinants are products of eigenvalues
        #   These are positive by construction
        # Could also test that the eigenvalues themselves are correct,
        #   but this seems sufficient.
        dets = [np.fabs(np.linalg.det(x)) for x in xs]
        dets_known = [np.prod(e) for e in eigs]
        assert_allclose(dets, dets_known, rtol=1e-13, atol=1e-13)

        # Test for 1's on the diagonal
        diags = [np.diag(x) for x in xs]
        for a, b in zip(diags, ones):
            assert_allclose(a, b, rtol=1e-13)

        # Correlation matrices are symmetric
        for x in xs:
            assert_allclose(x, x.T, rtol=1e-13)

    def test_to_corr(self):
        # Check some corner cases in to_corr

        # ajj == 1
        m = np.array([[0.1, 0], [0, 1]], dtype=float)
        m = random_correlation._to_corr(m)
        assert_allclose(m, np.array([[1, 0], [0, 0.1]]))

        # Floating point overflow; fails to compute the correct
        # rotation, but should still produce some valid rotation
        # rather than infs/nans
        with np.errstate(over='ignore'):
            g = np.array([[0, 1], [-1, 0]])

            m0 = np.array([[1e300, 0], [0, np.nextafter(1, 0)]], dtype=float)
            m = random_correlation._to_corr(m0.copy())
            assert_allclose(m, g.T.dot(m0).dot(g))

            m0 = np.array([[0.9, 1e300], [1e300, 1.1]], dtype=float)
            m = random_correlation._to_corr(m0.copy())
            assert_allclose(m, g.T.dot(m0).dot(g))

        # Zero discriminant; should set the first diag entry to 1
        m0 = np.array([[2, 1], [1, 2]], dtype=float)
        m = random_correlation._to_corr(m0.copy())
        assert_allclose(m[0,0], 1)

        # Slightly negative discriminant; should be approx correct still
        m0 = np.array([[2 + 1e-7, 1], [1, 2]], dtype=float)
        m = random_correlation._to_corr(m0.copy())
        assert_allclose(m[0,0], 1)


def check_pickling(distfn, args):
    # check that a distribution instance pickles and unpickles
    # pay special attention to the random_state property

    # save the random_state (restore later)
    rndm = distfn.random_state

    distfn.random_state = 1234
    distfn.rvs(*args, size=8)
    s = pickle.dumps(distfn)
    r0 = distfn.rvs(*args, size=8)

    unpickled = pickle.loads(s)
    r1 = unpickled.rvs(*args, size=8)
    assert_equal(r0, r1)

    # restore the random_state
    distfn.random_state = rndm


def test_random_state_property():
    scale = np.eye(3)
    scale[0, 1] = 0.5
    scale[1, 0] = 0.5
    dists = [
        [multivariate_normal, ()],
        [dirichlet, (np.array([1.]), )],
        [wishart, (10, scale)],
        [invwishart, (10, scale)]
    ]
    for distfn, args in dists:
        check_random_state_property(distfn, args)
        check_pickling(distfn, args)

if __name__ == "__main__":
    run_module_suite()
