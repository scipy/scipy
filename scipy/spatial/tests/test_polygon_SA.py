from __future__ import print_function
import numpy as np
from scipy.spatial import _polygon_surface_area as psa
from scipy.spatial.distance import pdist
from numpy.testing import (assert_equal, assert_raises,
                           assert_allclose,
                           assert_raises_regex)
from hypothesis import given
from hypothesis.strategies import floats

class TestDegenerateInput(object):
    # tests for problematic input
    def setUp(self):
        self.cython = None

    @given(floats(min_value=1e-20,
                  max_value=1e20),
           floats(min_value=1e-22,
                 max_value=1e-2))
    def test_duplicate_filter(self, radius,
                              threshold):
        # check that a ValueError is raised
        # for duplicate polygon vertices
        # within a user-specified threshold
        # ensure that the error message string
        # contains the correct values as well
        vertices = np.array([[-1,0,0],
                             [1,0,0],
                             [0,0,1]]) * radius
        vertices = np.concatenate((vertices,
                                   vertices))
        min_dist = pdist(vertices).min()
        expected_str = '''Duplicate vertices detected based on minimum
                     distance {min_dist} and threshold value
                     {threshold}.'''.format(min_dist=min_dist,
                                            threshold=threshold)
        with assert_raises_regex(ValueError, expected_str):
            psa.poly_area(vertices,
                          radius,
                          threshold,
                          cython=self.cython)

class TestDegenerateInputCython(TestDegenerateInput):
    def setUp(self):
        self.cython = 1

class TestSimpleAreas(object):
    # test polygon surface area calculations
    # for known / straightforward cases
    def setUp(self):
        self.cython = None

    # property test random subset of all floats 
    # for the radius value
    @given(floats(min_value=1e-20,
                  max_value=1e20)) 
    def test_half_hemisphere_area(self, radius):
        # the area of half a hemisphere should
        # be 1/4 the area of the entire sphere
        vertices = np.array([[-1,0,0],
                             [1,0,0],
                             [0,0,1]]) * radius
        expected_area = np.pi * (radius ** 2)
        actual_area = psa.poly_area(vertices=vertices,
                                    radius=radius,
                                    cython=self.cython)
        assert_allclose(actual_area, expected_area)

    @given(floats(min_value=1e-20,
                  max_value=1e20)) 
    def test_half_hemisphere_area_reverse_order(self, radius):
        # the area of half a hemisphere should
        # be 1/4 the area of the entire sphere
        # reverse order of vertex sorting
        vertices = np.array([[0,0,1],
                             [1,0,0],
                             [-1,0,0]]) * radius
        expected_area = np.pi * (radius ** 2)
        actual_area = psa.poly_area(vertices=vertices,
                                    radius=radius,
                                    cython=self.cython)
        assert_allclose(actual_area, expected_area)

    @given(floats(min_value=1e-20,
                  max_value=1e20))
    def test_quarter_hemisphere_area(self, radius):
        # the area of 1/4 of a hemisphere should
        # be 1/8 the area of the entire sphere
        vertices = np.array([[-1,0,0],
                             [0,1,0],
                             [0,0,1]]) * radius
        expected_area = (np.pi * (radius ** 2)) / 2.
        actual_area = psa.poly_area(vertices=vertices,
                                    radius=radius,
                                    cython=self.cython)
        assert_allclose(actual_area, expected_area)

    @given(floats(min_value=-1e20,
                  max_value=0.0)) 
    def test_zero_radius_area(self, radius):
        # an appropriate exception should be raised
        # for r <= 0.0
        vertices = np.array([[-1,0,0],
                             [1,0,0],
                             [0,0,1]]) * radius
        with assert_raises(ValueError):
            psa.poly_area(vertices=vertices,
                          radius=radius,
                          cython=self.cython)

    @given(floats(min_value=1e-20, max_value=1e20),
           floats(min_value=1e-20, max_value=1e20))
    def test_planar_triangle_area(self, base, height):
        # simple triangle area test
        # confirm that base * height / 2 result
        # is respected for a variety of base and
        # height values
        triangle_vertices = np.array([[0,0,0],
                                      [base,0,0],
                                      [base / 2.,height,0]])
        expected = 0.5 * base * height
        actual = psa.poly_area(vertices=triangle_vertices,
                               cython=self.cython)
        assert_allclose(actual, expected)

    @given(floats(min_value=1e-20, max_value=1e20),
           floats(min_value=1e-20, max_value=1e20))
    def test_planar_triangle_area_reverse(self, base, height):
        # simple triangle area test
        # confirm that base * height / 2 result
        # is respected for a variety of base and
        # height values
        # reverse vertex sort order
        triangle_vertices = np.array([[base / 2.,height,0],
                                      [base,0,0],
                                      [0,0,0]])
        expected = 0.5 * base * height
        actual = psa.poly_area(vertices=triangle_vertices,
                               cython=self.cython)
        assert_allclose(actual, expected)

class TestSimpleAreasCython(TestSimpleAreas):
    def setUp(self):
        self.cython=1

class TestRadianAreas(object):
    # compare spherical polygon surface areas
    # with values calculated using another well-known
    # approach
    # see: Weisstein, Eric W. "Spherical Polygon." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/SphericalPolygon.html
    # which cites: Beyer, W. H. CRC Standard Mathematical Tables, 28th ed. Boca Raton, FL: CRC Press, p. 131, 1987.
    # see also: Bevis and Cambereri (1987) Mathematical Geology 19: 335-346.
    # the above authors cite the angle-based equation used for reference
    # calcs here as well-known and stated frequently in handbooks of mathematics
    def setUp(self):
        self.cython=None

    def _angle_area(self, sum_radian_angles, n_vertices, radius):
        # suface area of spherical polygon using
        # alternate approach
        area = (sum_radian_angles - (n_vertices - 2) * np.pi) * (radius ** 2)
        return area

    @given(floats(min_value=1e-20,
                  max_value=1e20)) 
    def test_double_octant_area_both_orders(self, radius):
        # the octant of a sphere (1/8 of total area;
        # 1/4 of a hemisphere) is known to have 3
        # right angles
        # if we extend the octant to include both poles
        # of the sphere as vertices we end up with a 
        # diamond-shaped area with angles that can
        # be logically deduced
        # equatorial vertices should have angles that double
        # to pi
        # polar vertices should retain right angles at pi/2

        expected_area = self._angle_area((2 * np.pi) + (np.pi),
                                          4,
                                          radius)

        sample_vertices = np.array([[-1,0,0],
                                    [0,0,1],
                                    [0,1,0],
                                    [0,0,-1]]) * radius

        # check cw and ccw vertex sorting
        actual_area = psa.poly_area(vertices=sample_vertices,
                                    radius=radius,
                                    cython=self.cython)
        actual_area_reverse = psa.poly_area(vertices=sample_vertices[::-1],
                                    radius=radius,
                                    cython=self.cython)
        assert_allclose(actual_area, expected_area)
        assert_allclose(actual_area_reverse, expected_area)

class TestRadianAreasCython(TestRadianAreas):
    def setUp(self):
        self.cython = 1

class TestConvolutedAreas(object):
    # test more convoluted / tricky shapes
    # as input for surface area calculations

    def setUp(self):
        self.cython = None

    @given(floats(min_value=1e-20,
                  max_value=1e20))
    def test_spherical_three_sixteen(self, radius):
        # a 5 vertex spherical polygon consisting of
        # an octant (1/8 total sphere area) and
        # half an octant (1/16 total sphere area)
        expected_area = 4. * np.pi * (radius ** 2) * (3./16.)
        sample_vertices = np.array([[0,0,1],
                                    [-1,0,0],
                                    [0,0,-1],
                                    [0,1,0],
                                    [-0.5,0.5,0]]) * radius
        # check cw and ccw vertex sorting
        actual_area = psa.poly_area(vertices=sample_vertices,
                                    radius=radius,
                                    cython=self.cython)
        actual_area_reverse = psa.poly_area(vertices=sample_vertices[::-1],
                                    radius=radius,
                                    cython=self.cython)
        assert_allclose(actual_area, expected_area)
        assert_allclose(actual_area_reverse, expected_area)

class TestConvolutedAreasCython(TestConvolutedAreas):
    def setUp(self):
        self.cython = 1
