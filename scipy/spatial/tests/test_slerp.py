from __future__ import division, absolute_import, print_function

import math
import numpy as np
from numpy.testing import assert_equal, assert_allclose

import pytest
from scipy.spatial import geometric_slerp

def _generate_spherical_points(ndim=3, n_pts=2):
    # generate uniform points on sphere
    # see: https://stackoverflow.com/a/23785326
    # tentatively extended to arbitrary dims
    # for 0-sphere it will always produce antipodes
    np.random.seed(123)
    points = np.random.normal(size=(n_pts, ndim))
    points /= np.linalg.norm(points, axis=1)[:, np.newaxis]
    return points[0], points[1]

class TestGeometricSlerp(object):
    # Test various properties of the geometric slerp code

    @pytest.mark.parametrize("n_dims", [1, 2, 3, 5, 7, 9])
    @pytest.mark.parametrize("n_pts", [3, 17])
    def test_shape_property(self, n_dims, n_pts):
        # geometric_slerp output shape should match
        # input dimensionality & requested number
        # of interpolation points
        if n_dims > 1:
            start, end = _generate_spherical_points(n_dims, 2)
        else:
            start, end = (np.array([1]), np.array([0]))

        actual = geometric_slerp(start_coord=start,
                                 end_coord=end,
                                 t_values=np.linspace(0, 1, n_pts))

        assert actual.shape == (n_pts, n_dims)

    @pytest.mark.parametrize("n_dims", [1, 2, 3, 5, 7, 9])
    @pytest.mark.parametrize("n_pts", [3, 17])
    def test_include_ends(self, n_dims, n_pts):
        # geometric_slerp should return a data structure
        # that includes the start and end coordinates
        # when t includes 0 and 1 ends
        # this is convenient for plotting surfaces represented
        # by interpolations for example

        # the generator doesn't work so well for the unit
        # sphere (it always produces antipodes), so use
        # custom values there
        if n_dims > 1:
            start, end = _generate_spherical_points(n_dims, 2)
        else:
            start, end = (np.array([1]), np.array([0]))

        actual = geometric_slerp(start_coord=start,
                                 end_coord=end,
                                 t_values=np.linspace(0, 1, n_pts))

        assert_equal(actual[0], start)
        assert_equal(actual[-1], end)

    @pytest.mark.parametrize("start, end", [
        # both arrays are not flat
        (np.zeros((1,3)), np.ones((1,3))),
        # only start array is not flat
        (np.zeros((1,3)), np.ones(3)),
        # only end array is not flat
        (np.zeros(1), np.ones((3, 1))),
        ])
    def test_input_shape_flat(self, start, end):
        # geometric_slerp should handle input arrays that are
        # not flat appropriately
        with pytest.raises(ValueError, match='flat'):
            geometric_slerp(start_coord=start,
                            end_coord=end,
                            t_values=np.linspace(0, 1, 10))

    @pytest.mark.parametrize("start, end", [
        # 7-D and 3-D ends
        (np.zeros(7), np.ones(3)),
        # 2-D and 1-D ends
        (np.zeros(2), np.ones(1)),
        # empty, "3D" will also get caught this way
        (np.array([]), np.ones(3)),
        ])
    def test_input_dim_mismatch(self, start, end):
        # geometric_slerp must appropriately handle cases where
        # an interpolation is attempted across two different
        # dimensionalities
        with pytest.raises(ValueError, match='dimensions'):
            geometric_slerp(start_coord=start,
                            end_coord=end,
                            t_values=np.linspace(0, 1, 10))

    @pytest.mark.parametrize("start, end", [
        # both empty
        (np.array([]), np.array([])),
        ])
    def test_input_at_least1d(self, start, end):
        # empty inputs to geometric_slerp must
        # be handled appropriately when not detected
        # by mismatch
        with pytest.raises(ValueError, match='empty'):
            geometric_slerp(start_coord=start,
                            end_coord=end,
                            t_values=np.linspace(0, 1, 10))

    @pytest.mark.parametrize("start, end, expected", [
        # North and South Poles are definitely antipodes
        (np.array([0, 0, 1.0]), np.array([0, 0, -1.0]), "error"),
        # this case will error; North Pole was rotated very slightly
        # using r = R.from_euler('x', 0.035, degrees=True)
        # to achieve Euclidean distance offset from diameter by
        # 9.328908379124812e-08, within the default tol
        (np.array([0.00000000e+00,
                  -6.10865200e-04,
                  9.99999813e-01]), np.array([0, 0, -1.0]), "error"),
        # this case should succeed because a sufficiently large
        # rotation was applied to North Pole point to shift it
        # to a Euclidean distance of 2.3036691931821451e-07
        # from South Pole, which is larger than tol
        (np.array([0.00000000e+00,
                  -9.59930941e-04,
                  9.99999539e-01]), np.array([0, 0, -1.0]), "success"),
        ])
    def test_handle_antipodes(self, start, end, expected):
        # antipodal points must be handled appropriately;
        # there are an infinite number of possible geodesic
        # interpolations between them in higher dims
        if expected == "error":
            with pytest.raises(ValueError, match='antipodes'):
                geometric_slerp(start_coord=start,
                                end_coord=end,
                                t_values=np.linspace(0, 1, 10))
        else:
            geometric_slerp(start_coord=start,
                            end_coord=end,
                            t_values=np.linspace(0, 1, 10))

    @pytest.mark.parametrize("start, end, expected", [
        # 1-D is the 0-sphere, which is not path-connected
        # we effectively follow a 1-D coordinate along
        # a "circumference"
        (np.array([1]),
         np.array([0]),
         np.array([[1],
                   [math.sqrt(3) / 2],
                   [0.5],
                   [0]])),
        # 2-D with n_pts=4 (two new interpolation points)
        # this is an actual circle
        (np.array([1, 0]),
         np.array([0, 1]),
         np.array([[1, 0],
                   [math.sqrt(3) / 2, 0.5], # 30 deg on unit circle
                   [0.5, math.sqrt(3) / 2], # 60 deg on unit circle
                   [0, 1]])),
        # likewise for 3-D (add z = 0 plane)
        # this is an ordinary sphere
        (np.array([1, 0, 0]),
         np.array([0, 1, 0]),
         np.array([[1, 0, 0],
                   [math.sqrt(3) / 2, 0.5, 0],
                   [0.5, math.sqrt(3) / 2, 0],
                   [0, 1, 0]])),
        # for 5-D, pad more columns with constants
        # zeros are easiest--non-zero values on unit
        # circle are more difficult to reason about
        # at higher dims
        (np.array([1, 0, 0, 0, 0]),
         np.array([0, 1, 0, 0, 0]),
         np.array([[1, 0, 0, 0, 0],
                   [math.sqrt(3) / 2, 0.5, 0, 0, 0],
                   [0.5, math.sqrt(3) / 2, 0, 0, 0],
                   [0, 1, 0, 0, 0]])),

    ])
    def test_straightforward_examples(self, start, end, expected):
        # some straightforward interpolation tests, sufficiently
        # simple to use the unit circle to deduce expected values;
        # for larger dimensions, pad with constants so that the
        # data is N-D but simpler to reason about
        actual = geometric_slerp(start_coord=start,
                                 end_coord=end,
                                 t_values=np.linspace(0, 1, 4))
        assert_allclose(actual, expected)

    @pytest.mark.parametrize("t_values", [
        # both interval ends clearly violate limits
        np.linspace(-20, 20, 300),
        # only one interval end violating limit slightly
        np.linspace(-0.0001, 0.0001, 17),
        ])
    def test_t_values_limits(self, t_values):
        # geometric_slerp() should appropriately handle
        # interpolation parameters < 0 and > 1
        with pytest.raises(ValueError, match='interpolation parameter'):
            actual = geometric_slerp(start_coord=np.array([1]),
                                     end_coord=np.array([0]),
                                     t_values=t_values)

    @pytest.mark.parametrize("tol", [
        # an integer currently raises
        5,
        # string raises
        "7", 
        # list and arrays also raise
        [5, 6, 7], np.array(9.0),
        ])
    def test_tol_type(self, tol):
        # geometric_slerp() should raise if tol is not
        # a suitable float type
        with pytest.raises(ValueError, match='must be a float'):
            actual = geometric_slerp(start_coord=np.array([1]),
                                     end_coord=np.array([0]),
                                     t_values=np.linspace(0, 1, 5),
                                     tol=tol)

    @pytest.mark.parametrize("tol", [
        -5e-6,
        -7e-10,
        ])
    def test_tol_sign(self, tol):
        # geometric_slerp() currently handles negative
        # tol values, as long as they are floats
        actual = geometric_slerp(start_coord=np.array([1]),
                                 end_coord=np.array([0]),
                                 t_values=np.linspace(0, 1, 5),
                                 tol=tol)

    @pytest.mark.parametrize("start, end", [
        # 0-sphere with two points too far away
        (np.array([2]), np.array([0])),
        # 1-sphere (circle) with one point at origin
        # and the other on the circle
        (np.array([1, 0]), np.array([0, 0])),
        # 2-sphere (normal sphere) with both points
        # just slightly off sphere by the same amount
        # in different directions
        (np.array([1 + 1e-6, 0, 0]),
         np.array([0, 1 - 1e-6, 0])),
        # same thing in 4-D
        (np.array([1 + 1e-6, 0, 0, 0]),
         np.array([0, 1 - 1e-6, 0, 0])),
        ])
    def test_unit_sphere_enforcement(self, start, end):
        # geometric_slerp() should raise on input that clearly
        # cannot be on an n-sphere of radius 1
        with pytest.raises(ValueError, match='unit n-sphere'):
            geometric_slerp(start_coord=start,
                            end_coord=end,
                            t_values=np.linspace(0, 1, 5))

    @pytest.mark.parametrize("start, end", [
        # 0-sphere 90 degree case
        (np.array([1]), np.array([0])),
        # 1-sphere 45 degree case
        (np.array([1, 0]),
         np.array([math.sqrt(2)/ 2.,
                   math.sqrt(2)/ 2.])),
        # 2-sphere 135 degree case
        (np.array([1, 0]),
         np.array([-math.sqrt(2)/ 2.,
                   math.sqrt(2)/ 2.])),
        ])
    @pytest.mark.parametrize("t_func", [
        np.linspace, np.logspace])
    def test_order_handling(self, start, end, t_func):
        # geometric_slerp() should handle scenarios with
        # ascending and descending t value arrays gracefully;
        # results should simply be reversed

        # for scrambled / unsorted parameters, the same values
        # should be returned, just in scrambled order

        num_t_vals = 20
        np.random.seed(789)
        forward_t_vals = t_func(0, 10, num_t_vals)
        # normalize to max of 1
        forward_t_vals /= forward_t_vals.max()
        reverse_t_vals = np.flipud(forward_t_vals)
        shuffled_indices = np.arange(num_t_vals)
        np.random.shuffle(shuffled_indices)
        scramble_t_vals = forward_t_vals.copy()[shuffled_indices]

        forward_results = geometric_slerp(start_coord=start,
                                          end_coord=end,
                                          t_values=forward_t_vals)
        reverse_results = geometric_slerp(start_coord=start,
                                          end_coord=end,
                                          t_values=reverse_t_vals)
        scrambled_results = geometric_slerp(start_coord=start,
                                            end_coord=end,
                                            t_values=scramble_t_vals)

        # check fidelity to input order
        assert_allclose(forward_results, np.flipud(reverse_results))
        assert_allclose(forward_results[shuffled_indices],
                        scrambled_results)

    @pytest.mark.parametrize("t_values", [
        # string:
        "15, 5, 7",
        # complex numbers currently produce a warning
        # but not sure we need to worry about it too much:
        # [3 + 1j, 5 + 2j],
        # empty list:
        [],
        ])
    def test_t_values_conversion(self, t_values):
        with pytest.raises(ValueError):
            scrambled_results = geometric_slerp(start_coord=np.array([1]),
                                                end_coord=np.array([0]),
                                                t_values=t_values)
