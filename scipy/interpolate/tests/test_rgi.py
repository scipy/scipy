import itertools

import pytest
import numpy as np

from numpy.testing import (assert_allclose, assert_equal, assert_warns,
                           assert_array_almost_equal, assert_array_equal)
from pytest import raises as assert_raises

from scipy.interpolate import (RegularGridInterpolator, interpn,
                               RectBivariateSpline,
                               NearestNDInterpolator, LinearNDInterpolator,
                               CloughTocher2DInterpolator)

from scipy.sparse._sputils import matrix


class TestRegularGridInterpolator:
    def _get_sample_4d(self):
        # create a 4-D grid of 3 points in each dimension
        points = [(0., .5, 1.)] * 4
        values = np.asarray([0., .5, 1.])
        values0 = values[:, np.newaxis, np.newaxis, np.newaxis]
        values1 = values[np.newaxis, :, np.newaxis, np.newaxis]
        values2 = values[np.newaxis, np.newaxis, :, np.newaxis]
        values3 = values[np.newaxis, np.newaxis, np.newaxis, :]
        values = (values0 + values1 * 10 + values2 * 100 + values3 * 1000)
        return points, values

    def _get_sample_4d_2(self):
        # create another 4-D grid of 3 points in each dimension
        points = [(0., .5, 1.)] * 2 + [(0., 5., 10.)] * 2
        values = np.asarray([0., .5, 1.])
        values0 = values[:, np.newaxis, np.newaxis, np.newaxis]
        values1 = values[np.newaxis, :, np.newaxis, np.newaxis]
        values2 = values[np.newaxis, np.newaxis, :, np.newaxis]
        values3 = values[np.newaxis, np.newaxis, np.newaxis, :]
        values = (values0 + values1 * 10 + values2 * 100 + values3 * 1000)
        return points, values

    def _get_sample_4d_3(self):
        # create another 4-D grid of 7 points in each dimension
        points = [(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)] * 4
        values = np.asarray([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
        values0 = values[:, np.newaxis, np.newaxis, np.newaxis]
        values1 = values[np.newaxis, :, np.newaxis, np.newaxis]
        values2 = values[np.newaxis, np.newaxis, :, np.newaxis]
        values3 = values[np.newaxis, np.newaxis, np.newaxis, :]
        values = (values0 + values1 * 10 + values2 * 100 + values3 * 1000)
        return points, values

    def _get_sample_4d_4(self):
        # create another 4-D grid of 2 points in each dimension
        points = [(0.0, 1.0)] * 4
        values = np.asarray([0.0, 1.0])
        values0 = values[:, np.newaxis, np.newaxis, np.newaxis]
        values1 = values[np.newaxis, :, np.newaxis, np.newaxis]
        values2 = values[np.newaxis, np.newaxis, :, np.newaxis]
        values3 = values[np.newaxis, np.newaxis, np.newaxis, :]
        values = (values0 + values1 * 10 + values2 * 100 + values3 * 1000)
        return points, values

    def test_list_input(self):
        points, values = self._get_sample_4d_3()

        sample = np.asarray([[0.1, 0.1, 1., .9], [0.2, 0.1, .45, .8],
                             [0.5, 0.5, .5, .5]])

        for method in ['linear', 'nearest', 'slinear', 'cubic', 'quintic']:
            interp = RegularGridInterpolator(points,
                                             values.tolist(),
                                             method=method)
            v1 = interp(sample.tolist())
            interp = RegularGridInterpolator(points,
                                             values,
                                             method=method)
            v2 = interp(sample)
            assert_allclose(v1, v2)

    def test_spline_dim_error(self):
        points, values = self._get_sample_4d_4()
        match = "points in dimension"

        # Check error raise when creating interpolator
        for method in ['cubic', 'quintic']:
            with pytest.raises(ValueError, match=match):
                RegularGridInterpolator(points, values, method=method)

        # Check error raise when creating interpolator
        interp = RegularGridInterpolator(points, values)
        sample = np.asarray([[0.1, 0.1, 1., .9], [0.2, 0.1, .45, .8],
                             [0.5, 0.5, .5, .5]])
        for method in ['cubic', 'quintic']:
            with pytest.raises(ValueError, match=match):
                interp(sample, method=method)

    def test_linear_and_slinear_close_1(self):
        points, values = self._get_sample_4d()
        sample = np.asarray([[0.1, 0.1, 1., .9], [0.2, 0.1, .45, .8],
                             [0.5, 0.5, .5, .5]])
        self._assert_linear_and_slinear_close(points, sample, values)

    def test_linear_and_slinear_close_2(self):
        points, values = self._get_sample_4d_2()
        sample = np.asarray([0.1, 0.1, 10., 9.])
        self._assert_linear_and_slinear_close(points, sample, values)

    def _assert_linear_and_slinear_close(self, points, sample, values):
        interp = RegularGridInterpolator(points, values, method="linear")
        v1 = interp(sample)
        interp = RegularGridInterpolator(points, values, method="slinear")
        v2 = interp(sample)
        assert_allclose(v1, v2)

    def test_complex(self):
        points, values = self._get_sample_4d_3()
        values = values - 2j*values
        sample = np.asarray([[0.1, 0.1, 1., .9], [0.2, 0.1, .45, .8],
                             [0.5, 0.5, .5, .5]])

        for method in ['linear', 'nearest', "slinear", "cubic", "quintic"]:
            interp = RegularGridInterpolator(points, values,
                                             method=method)
            rinterp = RegularGridInterpolator(points, values.real,
                                              method=method)
            iinterp = RegularGridInterpolator(points, values.imag,
                                              method=method)

            v1 = interp(sample)
            v2 = rinterp(sample) + 1j*iinterp(sample)
            assert_allclose(v1, v2)

    def test_linear_xi1d(self):
        points, values = self._get_sample_4d_2()
        interp = RegularGridInterpolator(points, values)
        sample = np.asarray([0.1, 0.1, 10., 9.])
        wanted = 1001.1
        assert_array_almost_equal(interp(sample), wanted)

    def test_linear_xi3d(self):
        points, values = self._get_sample_4d()
        interp = RegularGridInterpolator(points, values)
        sample = np.asarray([[0.1, 0.1, 1., .9], [0.2, 0.1, .45, .8],
                             [0.5, 0.5, .5, .5]])
        wanted = np.asarray([1001.1, 846.2, 555.5])
        assert_array_almost_equal(interp(sample), wanted)

    def test_nearest(self):
        points, values = self._get_sample_4d()
        interp = RegularGridInterpolator(points, values, method="nearest")
        sample = np.asarray([0.1, 0.1, .9, .9])
        wanted = 1100.
        assert_array_almost_equal(interp(sample), wanted)
        sample = np.asarray([0.1, 0.1, 0.1, 0.1])
        wanted = 0.
        assert_array_almost_equal(interp(sample), wanted)
        sample = np.asarray([0., 0., 0., 0.])
        wanted = 0.
        assert_array_almost_equal(interp(sample), wanted)
        sample = np.asarray([1., 1., 1., 1.])
        wanted = 1111.
        assert_array_almost_equal(interp(sample), wanted)
        sample = np.asarray([0.1, 0.4, 0.6, 0.9])
        wanted = 1055.
        assert_array_almost_equal(interp(sample), wanted)

    def test_linear_edges(self):
        points, values = self._get_sample_4d()
        interp = RegularGridInterpolator(points, values)
        sample = np.asarray([[0., 0., 0., 0.], [1., 1., 1., 1.]])
        wanted = np.asarray([0., 1111.])
        assert_array_almost_equal(interp(sample), wanted)

    def test_valid_create(self):
        # create a 2-D grid of 3 points in each dimension
        points = [(0., .5, 1.), (0., 1., .5)]
        values = np.asarray([0., .5, 1.])
        values0 = values[:, np.newaxis]
        values1 = values[np.newaxis, :]
        values = (values0 + values1 * 10)
        assert_raises(ValueError, RegularGridInterpolator, points, values)
        points = [((0., .5, 1.), ), (0., .5, 1.)]
        assert_raises(ValueError, RegularGridInterpolator, points, values)
        points = [(0., .5, .75, 1.), (0., .5, 1.)]
        assert_raises(ValueError, RegularGridInterpolator, points, values)
        points = [(0., .5, 1.), (0., .5, 1.), (0., .5, 1.)]
        assert_raises(ValueError, RegularGridInterpolator, points, values)
        points = [(0., .5, 1.), (0., .5, 1.)]
        assert_raises(ValueError, RegularGridInterpolator, points, values,
                      method="undefmethod")

    def test_valid_call(self):
        points, values = self._get_sample_4d()
        interp = RegularGridInterpolator(points, values)
        sample = np.asarray([[0., 0., 0., 0.], [1., 1., 1., 1.]])
        assert_raises(ValueError, interp, sample, "undefmethod")
        sample = np.asarray([[0., 0., 0.], [1., 1., 1.]])
        assert_raises(ValueError, interp, sample)
        sample = np.asarray([[0., 0., 0., 0.], [1., 1., 1., 1.1]])
        assert_raises(ValueError, interp, sample)

    def test_out_of_bounds_extrap(self):
        points, values = self._get_sample_4d()
        interp = RegularGridInterpolator(points, values, bounds_error=False,
                                         fill_value=None)
        sample = np.asarray([[-.1, -.1, -.1, -.1], [1.1, 1.1, 1.1, 1.1],
                             [21, 2.1, -1.1, -11], [2.1, 2.1, -1.1, -1.1]])
        wanted = np.asarray([0., 1111., 11., 11.])
        assert_array_almost_equal(interp(sample, method="nearest"), wanted)
        wanted = np.asarray([-111.1, 1222.1, -11068., -1186.9])
        assert_array_almost_equal(interp(sample, method="linear"), wanted)

    def test_out_of_bounds_extrap2(self):
        points, values = self._get_sample_4d_2()
        interp = RegularGridInterpolator(points, values, bounds_error=False,
                                         fill_value=None)
        sample = np.asarray([[-.1, -.1, -.1, -.1], [1.1, 1.1, 1.1, 1.1],
                             [21, 2.1, -1.1, -11], [2.1, 2.1, -1.1, -1.1]])
        wanted = np.asarray([0., 11., 11., 11.])
        assert_array_almost_equal(interp(sample, method="nearest"), wanted)
        wanted = np.asarray([-12.1, 133.1, -1069., -97.9])
        assert_array_almost_equal(interp(sample, method="linear"), wanted)

    def test_out_of_bounds_fill(self):
        points, values = self._get_sample_4d()
        interp = RegularGridInterpolator(points, values, bounds_error=False,
                                         fill_value=np.nan)
        sample = np.asarray([[-.1, -.1, -.1, -.1], [1.1, 1.1, 1.1, 1.1],
                             [2.1, 2.1, -1.1, -1.1]])
        wanted = np.asarray([np.nan, np.nan, np.nan])
        assert_array_almost_equal(interp(sample, method="nearest"), wanted)
        assert_array_almost_equal(interp(sample, method="linear"), wanted)
        sample = np.asarray([[0.1, 0.1, 1., .9], [0.2, 0.1, .45, .8],
                             [0.5, 0.5, .5, .5]])
        wanted = np.asarray([1001.1, 846.2, 555.5])
        assert_array_almost_equal(interp(sample), wanted)

    def test_nearest_compare_qhull(self):
        points, values = self._get_sample_4d()
        interp = RegularGridInterpolator(points, values, method="nearest")
        points_qhull = itertools.product(*points)
        points_qhull = [p for p in points_qhull]
        points_qhull = np.asarray(points_qhull)
        values_qhull = values.reshape(-1)
        interp_qhull = NearestNDInterpolator(points_qhull, values_qhull)
        sample = np.asarray([[0.1, 0.1, 1., .9], [0.2, 0.1, .45, .8],
                             [0.5, 0.5, .5, .5]])
        assert_array_almost_equal(interp(sample), interp_qhull(sample))

    def test_linear_compare_qhull(self):
        points, values = self._get_sample_4d()
        interp = RegularGridInterpolator(points, values)
        points_qhull = itertools.product(*points)
        points_qhull = [p for p in points_qhull]
        points_qhull = np.asarray(points_qhull)
        values_qhull = values.reshape(-1)
        interp_qhull = LinearNDInterpolator(points_qhull, values_qhull)
        sample = np.asarray([[0.1, 0.1, 1., .9], [0.2, 0.1, .45, .8],
                             [0.5, 0.5, .5, .5]])
        assert_array_almost_equal(interp(sample), interp_qhull(sample))

    def test_duck_typed_values(self):
        x = np.linspace(0, 2, 5)
        y = np.linspace(0, 1, 7)

        values = MyValue((5, 7))

        for method in ('nearest', 'linear'):
            interp = RegularGridInterpolator((x, y), values,
                                             method=method)
            v1 = interp([0.4, 0.7])

            interp = RegularGridInterpolator((x, y), values._v,
                                             method=method)
            v2 = interp([0.4, 0.7])
            assert_allclose(v1, v2)

    def test_invalid_fill_value(self):
        np.random.seed(1234)
        x = np.linspace(0, 2, 5)
        y = np.linspace(0, 1, 7)
        values = np.random.rand(5, 7)

        # integers can be cast to floats
        RegularGridInterpolator((x, y), values, fill_value=1)

        # complex values cannot
        assert_raises(ValueError, RegularGridInterpolator,
                      (x, y), values, fill_value=1+2j)

    def test_fillvalue_type(self):
        # from #3703; test that interpolator object construction succeeds
        values = np.ones((10, 20, 30), dtype='>f4')
        points = [np.arange(n) for n in values.shape]
        # xi = [(1, 1, 1)]
        RegularGridInterpolator(points, values)
        RegularGridInterpolator(points, values, fill_value=0.)

    def test_length_one_axis(self):
        # gh-5890, gh-9524 : length-1 axis is legal for method='linear'.
        # Along the axis it's linear interpolation; away from the length-1
        # axis, it's an extrapolation, so fill_value should be used.
        def f(x, y):
            return x + y
        x = np.linspace(1, 1, 1)
        y = np.linspace(1, 10, 10)
        data = f(*np.meshgrid(x, y, indexing="ij", sparse=True))

        interp = RegularGridInterpolator((x, y), data, method="linear",
                                         bounds_error=False, fill_value=101)

        # check values at the grid
        assert_allclose(interp(np.array([[1, 1], [1, 5], [1, 10]])),
                        [2, 6, 11],
                        atol=1e-14)

        # check off-grid interpolation is indeed linear
        assert_allclose(interp(np.array([[1, 1.4], [1, 5.3], [1, 10]])),
                        [2.4, 6.3, 11],
                        atol=1e-14)

        # check exrapolation w/ fill_value
        assert_allclose(interp(np.array([1.1, 2.4])),
                        interp.fill_value,
                        atol=1e-14)

        # check extrapolation: linear along the `y` axis, const along `x`
        interp.fill_value = None
        assert_allclose(interp([[1, 0.3], [1, 11.5]]),
                        [1.3, 12.5], atol=1e-15)

        assert_allclose(interp([[1.5, 0.3], [1.9, 11.5]]),
                        [1.3, 12.5], atol=1e-15)

        # extrapolation with method='nearest'
        interp = RegularGridInterpolator((x, y), data, method="nearest",
                                         bounds_error=False, fill_value=None)
        assert_allclose(interp([[1.5, 1.8], [-4, 5.1]]),
                        [3, 6],
                        atol=1e-15)

    @pytest.mark.parametrize("fill_value", [None, np.nan, np.pi])
    @pytest.mark.parametrize("method", ['linear', 'nearest'])
    def test_length_one_axis2(self, fill_value, method):
        options = {"fill_value": fill_value, "bounds_error": False,
                   "method": method}

        x = np.linspace(0, 2*np.pi, 20)
        z = np.sin(x)

        fa = RegularGridInterpolator((x,), z[:], **options)
        fb = RegularGridInterpolator((x, [0]), z[:, None], **options)

        x1a = np.linspace(-1, 2*np.pi+1, 100)
        za = fa(x1a)

        # evaluated at provided y-value, fb should behave exactly as fa
        y1b = np.zeros(100)
        zb = fb(np.vstack([x1a, y1b]).T)
        assert_allclose(zb, za)

        # evaluated at a different y-value, fb should return fill value
        y1b = np.ones(100)
        zb = fb(np.vstack([x1a, y1b]).T)
        if fill_value is None:
            assert_allclose(zb, za)
        else:
            assert_allclose(zb, fill_value)

    @pytest.mark.parametrize("method", ['nearest', 'linear'])
    def test_nan_x_1d(self, method):
        # gh-6624 : if x is nan, result should be nan
        f = RegularGridInterpolator(([1, 2, 3],), [10, 20, 30], fill_value=1,
                                    bounds_error=False, method=method)
        assert np.isnan(f([np.nan]))

        # test arbitrary nan pattern
        rng = np.random.default_rng(8143215468)
        x = rng.random(size=100)*4
        i = rng.random(size=100) > 0.5
        x[i] = np.nan
        with np.errstate(invalid='ignore'):
            # out-of-bounds comparisons, `out_of_bounds += x < grid[0]`,
            # generate numpy warnings if `x` contains nans.
            # These warnings should propagate to user (since `x` is user
            # input) and we simply filter them out.
            res = f(x)

        assert_equal(res[i], np.nan)
        assert_equal(res[~i], f(x[~i]))

    @pytest.mark.parametrize("method", ['nearest', 'linear'])
    def test_nan_x_2d(self, method):
        x, y = np.array([0, 1, 2]), np.array([1, 3, 7])

        def f(x, y):
            return x**2 + y**2

        xg, yg = np.meshgrid(x, y, indexing='ij', sparse=True)
        data = f(xg, yg)
        interp = RegularGridInterpolator((x, y), data,
                                         method=method, bounds_error=False)

        with np.errstate(invalid='ignore'):
            res = interp([[1.5, np.nan], [1, 1]])
        assert_allclose(res[1], 2, atol=1e-14)
        assert np.isnan(res[0])

        # test arbitrary nan pattern
        rng = np.random.default_rng(8143215468)
        x = rng.random(size=100)*4-1
        y = rng.random(size=100)*8
        i1 = rng.random(size=100) > 0.5
        i2 = rng.random(size=100) > 0.5
        i = i1 | i2
        x[i1] = np.nan
        y[i2] = np.nan
        z = np.array([x, y]).T
        with np.errstate(invalid='ignore'):
            # out-of-bounds comparisons, `out_of_bounds += x < grid[0]`,
            # generate numpy warnings if `x` contains nans.
            # These warnings should propagate to user (since `x` is user
            # input) and we simply filter them out.
            res = interp(z)

        assert_equal(res[i], np.nan)
        assert_equal(res[~i], interp(z[~i]))

    def test_broadcastable_input(self):
        # input data
        np.random.seed(0)
        x = np.random.random(10)
        y = np.random.random(10)
        z = np.hypot(x, y)

        # x-y grid for interpolation
        X = np.linspace(min(x), max(x))
        Y = np.linspace(min(y), max(y))
        X, Y = np.meshgrid(X, Y)
        XY = np.vstack((X.ravel(), Y.ravel())).T

        for interpolator in (NearestNDInterpolator, LinearNDInterpolator,
                             CloughTocher2DInterpolator):
            interp = interpolator(list(zip(x, y)), z)
            # single array input
            interp_points0 = interp(XY)
            # tuple input
            interp_points1 = interp((X, Y))
            interp_points2 = interp((X, 0.0))
            # broadcastable input
            interp_points3 = interp(X, Y)
            interp_points4 = interp(X, 0.0)

            assert_equal(interp_points0.size ==
                         interp_points1.size ==
                         interp_points2.size ==
                         interp_points3.size ==
                         interp_points4.size, True)

    def test_read_only(self):
        # input data
        np.random.seed(0)
        xy = np.random.random((10, 2))
        x, y = xy[:, 0], xy[:, 1]
        z = np.hypot(x, y)

        # interpolation points
        XY = np.random.random((50, 2))

        xy.setflags(write=False)
        z.setflags(write=False)
        XY.setflags(write=False)

        for interpolator in (NearestNDInterpolator, LinearNDInterpolator,
                             CloughTocher2DInterpolator):
            interp = interpolator(xy, z)
            interp(XY)

    def test_descending_points(self):
        def val_func_3d(x, y, z):
            return 2 * x ** 3 + 3 * y ** 2 - z

        x = np.linspace(1, 4, 11)
        y = np.linspace(4, 7, 22)
        z = np.linspace(7, 9, 33)
        points = (x, y, z)
        values = val_func_3d(
            *np.meshgrid(*points, indexing='ij', sparse=True))
        my_interpolating_function = RegularGridInterpolator(points,
                                                            values)
        pts = np.array([[2.1, 6.2, 8.3], [3.3, 5.2, 7.1]])
        correct_result = my_interpolating_function(pts)

        # descending data
        x_descending = x[::-1]
        y_descending = y[::-1]
        z_descending = z[::-1]
        points_shuffled = (x_descending, y_descending, z_descending)
        values_shuffled = val_func_3d(
            *np.meshgrid(*points_shuffled, indexing='ij', sparse=True))
        my_interpolating_function = RegularGridInterpolator(
            points_shuffled, values_shuffled)
        test_result = my_interpolating_function(pts)

        assert_array_equal(correct_result, test_result)

    def test_invalid_points_order(self):
        def val_func_2d(x, y):
            return 2 * x ** 3 + 3 * y ** 2

        x = np.array([.5, 2., 0., 4., 5.5])  # not ascending or descending
        y = np.array([.5, 2., 3., 4., 5.5])
        points = (x, y)
        values = val_func_2d(*np.meshgrid(*points, indexing='ij',
                                          sparse=True))
        match = "must be strictly ascending or descending"
        with pytest.raises(ValueError, match=match):
            RegularGridInterpolator(points, values)


class MyValue:
    """
    Minimal indexable object
    """

    def __init__(self, shape):
        self.ndim = 2
        self.shape = shape
        self._v = np.arange(np.prod(shape)).reshape(shape)

    def __getitem__(self, idx):
        return self._v[idx]

    def __array_interface__(self):
        return None

    def __array__(self):
        raise RuntimeError("No array representation")


class TestInterpN:
    def _sample_2d_data(self):
        x = np.arange(1, 6)
        x = np.array([.5, 2., 3., 4., 5.5])
        y = np.arange(1, 6)
        y = np.array([.5, 2., 3., 4., 5.5])
        z = np.array([[1, 2, 1, 2, 1], [1, 2, 1, 2, 1], [1, 2, 3, 2, 1],
                      [1, 2, 2, 2, 1], [1, 2, 1, 2, 1]])
        return x, y, z

    def test_spline_2d(self):
        x, y, z = self._sample_2d_data()
        lut = RectBivariateSpline(x, y, z)

        xi = np.array([[1, 2.3, 5.3, 0.5, 3.3, 1.2, 3],
                       [1, 3.3, 1.2, 4.0, 5.0, 1.0, 3]]).T
        assert_array_almost_equal(interpn((x, y), z, xi, method="splinef2d"),
                                  lut.ev(xi[:, 0], xi[:, 1]))

    def test_list_input(self):
        x, y, z = self._sample_2d_data()
        xi = np.array([[1, 2.3, 5.3, 0.5, 3.3, 1.2, 3],
                       [1, 3.3, 1.2, 4.0, 5.0, 1.0, 3]]).T

        for method in ['nearest', 'linear', 'splinef2d']:
            v1 = interpn((x, y), z, xi, method=method)
            v2 = interpn((x.tolist(), y.tolist()), z.tolist(),
                         xi.tolist(), method=method)
            assert_allclose(v1, v2, err_msg=method)

    def test_spline_2d_outofbounds(self):
        x = np.array([.5, 2., 3., 4., 5.5])
        y = np.array([.5, 2., 3., 4., 5.5])
        z = np.array([[1, 2, 1, 2, 1], [1, 2, 1, 2, 1], [1, 2, 3, 2, 1],
                      [1, 2, 2, 2, 1], [1, 2, 1, 2, 1]])
        lut = RectBivariateSpline(x, y, z)

        xi = np.array([[1, 2.3, 6.3, 0.5, 3.3, 1.2, 3],
                       [1, 3.3, 1.2, -4.0, 5.0, 1.0, 3]]).T
        actual = interpn((x, y), z, xi, method="splinef2d",
                         bounds_error=False, fill_value=999.99)
        expected = lut.ev(xi[:, 0], xi[:, 1])
        expected[2:4] = 999.99
        assert_array_almost_equal(actual, expected)

        # no extrapolation for splinef2d
        assert_raises(ValueError, interpn, (x, y), z, xi, method="splinef2d",
                      bounds_error=False, fill_value=None)

    def _sample_4d_data(self):
        points = [(0., .5, 1.)] * 2 + [(0., 5., 10.)] * 2
        values = np.asarray([0., .5, 1.])
        values0 = values[:, np.newaxis, np.newaxis, np.newaxis]
        values1 = values[np.newaxis, :, np.newaxis, np.newaxis]
        values2 = values[np.newaxis, np.newaxis, :, np.newaxis]
        values3 = values[np.newaxis, np.newaxis, np.newaxis, :]
        values = (values0 + values1 * 10 + values2 * 100 + values3 * 1000)
        return points, values

    def test_linear_4d(self):
        # create a 4-D grid of 3 points in each dimension
        points, values = self._sample_4d_data()
        interp_rg = RegularGridInterpolator(points, values)
        sample = np.asarray([[0.1, 0.1, 10., 9.]])
        wanted = interpn(points, values, sample, method="linear")
        assert_array_almost_equal(interp_rg(sample), wanted)

    def test_4d_linear_outofbounds(self):
        # create a 4-D grid of 3 points in each dimension
        points, values = self._sample_4d_data()
        sample = np.asarray([[0.1, -0.1, 10.1, 9.]])
        wanted = 999.99
        actual = interpn(points, values, sample, method="linear",
                         bounds_error=False, fill_value=999.99)
        assert_array_almost_equal(actual, wanted)

    def test_nearest_4d(self):
        # create a 4-D grid of 3 points in each dimension
        points, values = self._sample_4d_data()
        interp_rg = RegularGridInterpolator(points, values, method="nearest")
        sample = np.asarray([[0.1, 0.1, 10., 9.]])
        wanted = interpn(points, values, sample, method="nearest")
        assert_array_almost_equal(interp_rg(sample), wanted)

    def test_4d_nearest_outofbounds(self):
        # create a 4-D grid of 3 points in each dimension
        points, values = self._sample_4d_data()
        sample = np.asarray([[0.1, -0.1, 10.1, 9.]])
        wanted = 999.99
        actual = interpn(points, values, sample, method="nearest",
                         bounds_error=False, fill_value=999.99)
        assert_array_almost_equal(actual, wanted)

    def test_xi_1d(self):
        # verify that 1-D xi works as expected
        points, values = self._sample_4d_data()
        sample = np.asarray([0.1, 0.1, 10., 9.])
        v1 = interpn(points, values, sample, bounds_error=False)
        v2 = interpn(points, values, sample[None,:], bounds_error=False)
        assert_allclose(v1, v2)

    def test_xi_nd(self):
        # verify that higher-d xi works as expected
        points, values = self._sample_4d_data()

        np.random.seed(1234)
        sample = np.random.rand(2, 3, 4)

        v1 = interpn(points, values, sample, method='nearest',
                     bounds_error=False)
        assert_equal(v1.shape, (2, 3))

        v2 = interpn(points, values, sample.reshape(-1, 4),
                     method='nearest', bounds_error=False)
        assert_allclose(v1, v2.reshape(v1.shape))

    def test_xi_broadcast(self):
        # verify that the interpolators broadcast xi
        x, y, values = self._sample_2d_data()
        points = (x, y)

        xi = np.linspace(0, 1, 2)
        yi = np.linspace(0, 3, 3)

        for method in ['nearest', 'linear', 'splinef2d']:
            sample = (xi[:,None], yi[None,:])
            v1 = interpn(points, values, sample, method=method,
                         bounds_error=False)
            assert_equal(v1.shape, (2, 3))

            xx, yy = np.meshgrid(xi, yi)
            sample = np.c_[xx.T.ravel(), yy.T.ravel()]

            v2 = interpn(points, values, sample,
                         method=method, bounds_error=False)
            assert_allclose(v1, v2.reshape(v1.shape))

    def test_nonscalar_values(self):
        # Verify that non-scalar valued values also works
        points, values = self._sample_4d_data()

        np.random.seed(1234)
        values = np.random.rand(3, 3, 3, 3, 6)
        sample = np.random.rand(7, 11, 4)

        for method in ['nearest', 'linear']:
            v = interpn(points, values, sample, method=method,
                        bounds_error=False)
            assert_equal(v.shape, (7, 11, 6), err_msg=method)

            vs = [interpn(points, values[...,j], sample, method=method,
                          bounds_error=False)
                  for j in range(6)]
            v2 = np.array(vs).transpose(1, 2, 0)

            assert_allclose(v, v2, err_msg=method)

        # Vector-valued splines supported with fitpack
        assert_raises(ValueError, interpn, points, values, sample,
                      method='splinef2d')

    def test_complex(self):
        x, y, values = self._sample_2d_data()
        points = (x, y)
        values = values - 2j*values

        sample = np.array([[1, 2.3, 5.3, 0.5, 3.3, 1.2, 3],
                           [1, 3.3, 1.2, 4.0, 5.0, 1.0, 3]]).T

        for method in ['linear', 'nearest']:
            v1 = interpn(points, values, sample, method=method)
            v2r = interpn(points, values.real, sample, method=method)
            v2i = interpn(points, values.imag, sample, method=method)
            v2 = v2r + 1j*v2i
            assert_allclose(v1, v2)

        # Complex-valued data not supported by spline2fd
        assert_warns(np.ComplexWarning, interpn, points, values,
                     sample, method='splinef2d')

    def test_duck_typed_values(self):
        x = np.linspace(0, 2, 5)
        y = np.linspace(0, 1, 7)

        values = MyValue((5, 7))

        for method in ('nearest', 'linear'):
            v1 = interpn((x, y), values, [0.4, 0.7], method=method)
            v2 = interpn((x, y), values._v, [0.4, 0.7], method=method)
            assert_allclose(v1, v2)

    def test_matrix_input(self):
        x = np.linspace(0, 2, 5)
        y = np.linspace(0, 1, 7)

        values = matrix(np.random.rand(5, 7))

        sample = np.random.rand(3, 7, 2)

        for method in ('nearest', 'linear', 'splinef2d'):
            v1 = interpn((x, y), values, sample, method=method)
            v2 = interpn((x, y), np.asarray(values), sample, method=method)
            assert_allclose(v1, v2)

    def test_length_one_axis(self):
        # gh-5890, gh-9524 : length-1 axis is legal for method='linear'.
        # Along the axis it's linear interpolation; away from the length-1
        # axis, it's an extrapolation, so fill_value should be used.

        values = np.array([[0.1, 1, 10]])
        xi = np.array([[1, 2.2], [1, 3.2], [1, 3.8]])

        res = interpn(([1], [2, 3, 4]), values, xi)
        wanted = [0.9*0.2 + 0.1,   # on [2, 3) it's 0.9*(x-2) + 0.1
                  9*0.2 + 1,       # on [3, 4] it's 9*(x-3) + 1
                  9*0.8 + 1]

        assert_allclose(res, wanted, atol=1e-15)

        # check extrapolation
        xi = np.array([[1.1, 2.2], [1.5, 3.2], [-2.3, 3.8]])
        res = interpn(([1], [2, 3, 4]), values, xi,
                      bounds_error=False, fill_value=None)

        assert_allclose(res, wanted, atol=1e-15)

    def test_descending_points(self):
        def value_func_4d(x, y, z, a):
            return 2 * x ** 3 + 3 * y ** 2 - z - a

        x1 = np.array([0, 1, 2, 3])
        x2 = np.array([0, 10, 20, 30])
        x3 = np.array([0, 10, 20, 30])
        x4 = np.array([0, .1, .2, .30])
        points = (x1, x2, x3, x4)
        values = value_func_4d(
            *np.meshgrid(*points, indexing='ij', sparse=True))
        pts = (0.1, 0.3, np.transpose(np.linspace(0, 30, 4)),
               np.linspace(0, 0.3, 4))
        correct_result = interpn(points, values, pts)

        x1_descend = x1[::-1]
        x2_descend = x2[::-1]
        x3_descend = x3[::-1]
        x4_descend = x4[::-1]
        points_shuffled = (x1_descend, x2_descend, x3_descend, x4_descend)
        values_shuffled = value_func_4d(
            *np.meshgrid(*points_shuffled, indexing='ij', sparse=True))
        test_result = interpn(points_shuffled, values_shuffled, pts)

        assert_array_equal(correct_result, test_result)

    def test_invalid_points_order(self):
        x = np.array([.5, 2., 0., 4., 5.5])  # not ascending or descending
        y = np.array([.5, 2., 3., 4., 5.5])
        z = np.array([[1, 2, 1, 2, 1], [1, 2, 1, 2, 1], [1, 2, 3, 2, 1],
                      [1, 2, 2, 2, 1], [1, 2, 1, 2, 1]])
        xi = np.array([[1, 2.3, 6.3, 0.5, 3.3, 1.2, 3],
                       [1, 3.3, 1.2, -4.0, 5.0, 1.0, 3]]).T

        match = "must be strictly ascending or descending"
        with pytest.raises(ValueError, match=match):
            interpn((x, y), z, xi)
