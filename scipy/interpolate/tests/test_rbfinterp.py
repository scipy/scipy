import numpy as np
from numpy.testing import assert_, assert_allclose, assert_raises
from scipy.spatial import cKDTree
from scipy.interpolate.rbfinterp import (
    _NAME_TO_MIN_DEGREE, _NAME_TO_FUNC, _SCALE_INVARIANT,
    _vandermonde, _distance,
    RBFInterpolator, KNearestRBFInterpolator
    )


def _1d_test_function(x):
    # test function used in Wahba's "Spline Models for Observational Data
    # domain ~= (0, 3), range ~= (-1.0, 0.2)"
    x = x[:, 0]
    y =  4.26*(np.exp(-x) - 4*np.exp(-2*x) + 3*np.exp(-3*x))
    return y


def _2d_test_function(x):
    # Franke's test function
    # domain ~= (0, 1) X (0, 1), range ~= (0.0, 1.2)"
    x1, x2 = x[:, 0], x[:, 1]
    term1 = 0.75 * np.exp(-(9*x1-2)**2/4 - (9*x2-2)**2/4)
    term2 = 0.75 * np.exp(-(9*x1+1)**2/49 - (9*x2+1)/10)
    term3 = 0.5 * np.exp(-(9*x1-7)**2/4 - (9*x2-3)**2/4)
    term4 = -0.2 * np.exp(-(9*x1-4)**2 - (9*x2-7)**2)
    y = term1 + term2 + term3 + term4
    return y


def _is_conditionally_positive_definite(kernel, m):
    # Tests whether the kernel is conditionally positive definite of order m.
    # See chapter 8 of Fasshauer's "Meshfree Approximation Methods with
    # MATLAB".
    np.random.seed(0)
    nx = 10
    ntests = 100
    for ndim in [1, 2, 3, 4, 5]:
        for _ in range(ntests):
            x = np.random.normal(0.0, 1.0, (nx, ndim))
            # add a small value to the diagonals in case the kernel _should_ be
            # c.p.d but it is not due to numerical precision.
            A = kernel(_distance(x, x)) + 1e-12*np.eye(nx)
            P = _vandermonde(x, m - 1)
            Q, R  = np.linalg.qr(P, mode='complete')
            # Q2 forms a basis spanning the space where P.T.dot(x) = 0. Project
            # A onto this space, and then see if it is positive definite using
            # the Cholesky decomposition. If not, then the kernel is not c.p.d
            # of order m
            Q2 = Q[:, P.shape[1]:]
            B = Q2.T.dot(A).dot(Q2)
            try:
                np.linalg.cholesky(B)
            except np.linalg.LinAlgError:
                return False

    return True


def test_conditionally_positive_definite():
    # Test if each kernel in _NAME_TO_FUNC is conditionally positive definite
    # of order m, where m comes from _NAME_TO_MIN_DEGREE. This is a necessary
    # condition for the smoothed RBF interpolant to be well-posed in general
    for kernel_name, kernel_func in _NAME_TO_FUNC.items():
        m = _NAME_TO_MIN_DEGREE.get(kernel_name, -1) + 1
        assert_(_is_conditionally_positive_definite(kernel_func, m))


class _TestRBFInterpolator:
    def test_scale_invariance_1d(self):
        # Verify that the functions in _SCALE_INVARIANT are insensitive to the
        # shape parameter (when smoothing == 0) in 1-D
        np.random.seed(0)
        x = np.random.uniform(0.0, 3.0, (50, 1))
        y = _1d_test_function(x)
        xitp = np.random.uniform(0.0, 3.0, (50, 1))
        for kernel in _SCALE_INVARIANT:
            yitp1 = self.build(x, y, epsilon=1.0, kernel=kernel)(xitp)
            yitp2 = self.build(x, y, epsilon=2.0, kernel=kernel)(xitp)
            assert_allclose(yitp1, yitp2, atol=1e-8)

    def test_scale_invariance_2d(self):
        # Verify that the functions in _SCALE_INVARIANT are insensitive to the
        # shape parameter (when smoothing == 0) in 2-D
        np.random.seed(0)
        x = np.random.uniform(0.0, 1.0, (100, 2))
        y = _2d_test_function(x)
        xitp = np.random.uniform(0.0, 1.0, (100, 2))
        for kernel in _SCALE_INVARIANT:
            yitp1 = self.build(x, y, epsilon=1.0, kernel=kernel)(xitp)
            yitp2 = self.build(x, y, epsilon=2.0, kernel=kernel)(xitp)
            assert_allclose(yitp1, yitp2, atol=1e-8)

    def test_extreme_domains(self):
        # Make sure the interpolant remains numerically stable for very
        # large/small domains
        np.random.seed(0)
        scale = 1e50
        shift = 1e56

        x = np.random.uniform(0.0, 1.0, (100, 2))
        y = _2d_test_function(x)
        xitp = np.random.uniform(0.0, 1.0, (100, 2))

        for kernel in _NAME_TO_FUNC.keys():
            if kernel in _SCALE_INVARIANT:
                yitp1 = self.build(x, y, kernel=kernel)(xitp)
                yitp2 = self.build(
                    x*scale + shift, y,
                    kernel=kernel
                    )(xitp*scale + shift)
            else:
                yitp1 = self.build(x, y, epsilon=5.0, kernel=kernel)(xitp)
                yitp2 = self.build(
                    x*scale + shift, y,
                    epsilon=5.0/scale,
                    kernel=kernel
                    )(xitp*scale + shift)

            assert_allclose(yitp1, yitp2, atol=1e-8)

    def test_polynomial_reproduction(self):
        # If the observed data comes from a polynomial, then the interpolant
        # should be able to reproduce the polynomial exactly, provided that
        # `degree` is sufficiently high
        np.random.seed(0)
        ndim = 2
        degree = 3

        x = np.random.normal(0.0, 1.0, (50, ndim))
        xitp = np.random.uniform(0.0, 1.0, (50, ndim))

        P = _vandermonde(x, degree)
        Pitp = _vandermonde(xitp, degree)

        poly_coeffs = np.random.normal(0.0, 1.0, P.shape[1])

        y = P.dot(poly_coeffs)
        yitp1 = Pitp.dot(poly_coeffs)
        yitp2 = self.build(x, y, degree=degree)(xitp)

        assert_allclose(yitp1, yitp2, atol=1e-8)

    def test_chunking(self):
        # Make sure we get exactly the same results with and without chunking
        # the interpolation points
        np.random.seed(0)

        x = np.random.uniform(0.0, 3.0, (50, 1))
        y = _1d_test_function(x)
        xitp = np.random.uniform(0.0, 3.0, (50, 1))
        interp = self.build(x, y)
        yitp1 = interp(xitp, chunk_size=1)
        yitp2 = interp(xitp, chunk_size=30)
        yitp3 = interp(xitp, chunk_size=50)
        yitp4 = interp(xitp, chunk_size=None)
        assert_allclose(yitp1, yitp2)
        assert_allclose(yitp1, yitp3)
        assert_allclose(yitp1, yitp4)

    def test_vector_data(self):
        # Make sure interpolating a vector field is the same as interpolating
        # each component separately
        np.random.seed(0)

        x = np.random.uniform(0.0, 1.0, (100, 2))
        y = np.array([_2d_test_function(x),
                      _2d_test_function(x[:, ::-1])]).T
        xitp = np.random.uniform(0.0, 1.0, (100, 2))
        yitp1 = self.build(x, y)(xitp)
        yitp2 = self.build(x, y[:, 0])(xitp)
        yitp3 = self.build(x, y[:, 1])(xitp)
        assert_allclose(yitp1[:, 0], yitp2)
        assert_allclose(yitp1[:, 1], yitp3)

    def test_complex_data(self):
        # Interpolating complex input should be the same as interpolating the
        # real and complex components
        np.random.seed(0)

        x = np.random.uniform(0.0, 1.0, (100, 2))
        y = _2d_test_function(x) + 1j*_2d_test_function(x[:, ::-1])
        xitp = np.random.uniform(0.0, 1.0, (100, 2))
        yitp1 = self.build(x, y)(xitp)
        yitp2 = self.build(x, y.real)(xitp)
        yitp3 = self.build(x, y.imag)(xitp)
        assert_allclose(yitp1.real, yitp2)
        assert_allclose(yitp1.imag, yitp3)

    def test_interpolation_misfit_1d(self):
        # Make sure that each kernel, with its default `degree` and an
        # appropriate `epsilon`, does a good job at interpolation in 1d
        np.random.seed(0)

        x = np.random.uniform(0.0, 3.0, (50, 1))
        y = _1d_test_function(x)

        # use min/max of x so we dont extrapolate
        xitp = np.random.uniform(np.min(x), np.max(x), (50, 1))
        ytrue = _1d_test_function(xitp)

        for kernel in _NAME_TO_FUNC.keys():
            yitp = self.build(x, y, epsilon=0.6, kernel=kernel)(xitp)
            mse = np.mean((yitp - ytrue)**2)
            assert_(mse < 1.0e-4)

    def test_interpolation_misfit_2d(self):
        # Make sure that each kernel, with its default `degree` and an
        # appropriate `epsilon`, does a good job at interpolation in 2d
        np.random.seed(0)

        x = np.random.uniform(0.0, 1.0, (100, 2))
        y = _2d_test_function(x)

        xitp = np.random.uniform(0.0, 1.0, (100, 2))
        ytrue = _2d_test_function(xitp)

        for kernel in _NAME_TO_FUNC.keys():
            yitp = self.build(x, y, epsilon=5.0, kernel=kernel)(xitp)
            mse = np.mean((yitp - ytrue)**2)
            assert_(mse < 2.0e-4)

    def test_smoothing_misfit(self):
        # Make sure we can find a smoothing parameter for each kernel that
        # removes a sufficient amount of noise
        np.random.seed(0)

        noise = 0.2
        rmse_tol = 0.1
        smoothing_range = 10**np.linspace(-4, 1, 20)

        x = np.random.uniform(0.0, 3.0, (100, 1))
        y = _1d_test_function(x) + np.random.normal(0.0, noise, (100,))
        ytrue = _1d_test_function(x)
        for kernel in _NAME_TO_FUNC.keys():
            rmse_within_tol = False
            for smoothing in smoothing_range:
                ysmooth = self.build(
                    x, y,
                    epsilon=1.0,
                    smoothing=smoothing,
                    kernel=kernel)(x)
                rmse = np.sqrt(np.mean((ysmooth - ytrue)**2))
                if rmse < rmse_tol:
                    rmse_within_tol = True
                    break

            assert_(rmse_within_tol)

    def test_array_smoothing(self):
        # test using an array for the smoothing parameter to give less weight
        # to a known outlier
        np.random.seed(0)
        ndim = 1
        degree = 2

        x = np.random.normal(0.0, 1.0, (50, ndim))
        P = _vandermonde(x, degree)
        poly_coeffs = np.random.normal(0.0, 1.0, P.shape[1])
        y = P.dot(poly_coeffs)
        y_with_outlier = np.copy(y)
        y_with_outlier[10] += 1.0
        smoothing = np.zeros((50,))
        smoothing[10] = 1000.0
        yitp = self.build(x, y_with_outlier, smoothing=smoothing)(x)
        # should be able to reproduce the uncorrupted data almost exactly
        assert_allclose(yitp, y, atol=1e-5)

    def test_inconsistent_dimensions_error(self):
        # ValueError should be raised if the observations points and
        # interpolation points have a different number of dimensions
        np.random.seed(0)
        x = np.random.uniform(0.0, 1.0, (100, 2))
        y = _2d_test_function(x)
        xitp = np.random.uniform(0.0, 1.0, (100, 1))
        with assert_raises(ValueError):
            self.build(x, y)(xitp)


class TestRBFInterpolator(_TestRBFInterpolator):
    def build(self, *args, **kwargs):
        return RBFInterpolator(*args, **kwargs)

    def test_smoothing_limit_1d(self):
        # For large smoothing parameters, the interpolant should approach a
        # least squares fit of a polynomial with the specified degree
        np.random.seed(0)

        degree = 3
        smoothing = 1e8

        x = np.random.uniform(0.0, 3.0, (50, 1))
        y = _1d_test_function(x)

        xitp = np.random.uniform(0.0, 3.0, (50, 1))
        yitp1 = self.build(
            x, y,
            degree=degree,
            smoothing=smoothing
            )(xitp)

        P = _vandermonde(x, degree)
        Pitp = _vandermonde(xitp, degree)
        yitp2 = Pitp.dot(np.linalg.lstsq(P, y, rcond=None)[0])
        assert_allclose(yitp1, yitp2, atol=1e-8)

    def test_smoothing_limit_2d(self):
        # For large smoothing parameters, the interpolant should approach a
        # least squares fit of a polynomial with the specified degree
        np.random.seed(0)

        degree = 3
        smoothing = 1e8

        x = np.random.uniform(0.0, 1.0, (100, 2))
        y = _2d_test_function(x)

        xitp = np.random.uniform(0.0, 1.0, (100, 2))
        yitp1 = self.build(
            x, y,
            degree=degree,
            smoothing=smoothing
            )(xitp)

        P = _vandermonde(x, degree)
        Pitp = _vandermonde(xitp, degree)
        yitp2 = Pitp.dot(np.linalg.lstsq(P, y, rcond=None)[0])

        assert_allclose(yitp1, yitp2, atol=1e-8)


class TestKNearestRBFInterpolator20(_TestRBFInterpolator):
    # KNearestRBFInterpolator using 20 nearest neighbors
    def build(self, *args, **kwargs):
        return KNearestRBFInterpolator(*args, **kwargs, k=20)

    def test_equivalent_to_rbf_interpolator(self):
        # Make sure this is equivalent to using RBFInterpolator with the 20
        # nearest observations
        np.random.seed(0)
        x = np.random.uniform(0.0, 3.0, (50, 1))
        y = _1d_test_function(x)
        xitp = np.random.uniform(0.0, 3.0, (50, 1))
        yitp1 = self.build(x, y)(xitp)

        yitp2 = []
        tree = cKDTree(x)
        for xi in xitp:
            _, nbr = tree.query(xi, 20)
            yitp2.append(RBFInterpolator(x[nbr], y[nbr])(xi[None])[0])

        assert_allclose(yitp1, yitp2, atol=1e-8)


class TestKNearestRBFInterpolatorInf(TestRBFInterpolator):
    # KNearestRBFInterpolator using all points. This should behave exactly like
    # RBFInterpolator
    def build(self, *args, **kwargs):
        return KNearestRBFInterpolator(*args, **kwargs, k=np.inf)

    def test_equivalent_to_rbf_interpolator(self):
        np.random.seed(0)
        x = np.random.uniform(0.0, 3.0, (50, 1))
        y = _1d_test_function(x)
        xitp = np.random.uniform(0.0, 3.0, (50, 1))
        yitp1 = self.build(x, y)(xitp)
        yitp2 = RBFInterpolator(x, y)(xitp)
        assert_allclose(yitp1, yitp2, atol=1e-8)
