from __future__ import print_function
import numpy as np
from scipy.spatial import _polygon_surface_area as psa
from numpy.testing import assert_equal, assert_raises
from hypothesis import given
from hypothesis.strategies import floats

class TestSimpleAreas(object):
    # test polygon surface area calculations
    # for known / straightforward cases

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
                                    radius=radius)
        assert_equal(actual_area, expected_area)

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
                                    radius=radius)
        assert_equal(actual_area, expected_area)

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
                          radius=radius)

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
        actual = psa.poly_area(vertices=triangle_vertices)
        assert_equal(actual, expected)
