import numpy as np
from scipy.spatial import _polygon_surface_area as psa
from scipy.spatial.distance import pdist
from numpy.testing import (assert_equal, assert_raises,
                           assert_allclose,
                           assert_raises_regex)
import split_spherical_triangle
import pytest

@pytest.mark.parametrize("cython, radius, threshold", [
    (None, 1e-20, 1e-22),
    (1, 1e-20, 1e-22),
    (None, 1e-10, 1e-12),
    (1, 1e-10, 1e-12),
    (None, 1e10, 1e-2),
    (1, 1e10, 1e-2),
    ])
def test_duplicate_filter(cython,
                          radius,
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
    with pytest.raises(ValueError) as excinfo:
        psa.poly_area(vertices,
                      radius,
                      threshold,
                      cython=cython)
        assert str(excinfo.value) == expected_str

# property test random subset of floats 
# for the radius value
@pytest.mark.parametrize("cython, radius", [
(None, 1e-20),
(1, 1e-20),
(None, 1e20),
(1, 1e20),
(None, 1e10),
(1, 1e10),
(None, 0.05),
(1, 0.05),
])
class TestSimpleAreas(object):
    # test polygon surface area calculations
    # for known / straightforward cases
    def test_half_hemisphere_area(self, cython, radius):
        # the area of half a hemisphere should
        # be 1/4 the area of the entire sphere
        vertices = np.array([[-1,0,0],
                             [1,0,0],
                             [0,0,1]]) * radius
        expected_area = np.pi * (radius ** 2)
        actual_area = psa.poly_area(vertices=vertices,
                                    radius=radius,
                                    cython=cython,
                                    discretizations=7000)
        assert_allclose(actual_area, expected_area)

    def test_half_hemisphere_area_reverse_order(self, cython, radius):
        # the area of half a hemisphere should
        # be 1/4 the area of the entire sphere
        # reverse order of vertex sorting
        vertices = np.array([[0,0,1],
                             [1,0,0],
                             [-1,0,0]]) * radius
        expected_area = np.pi * (radius ** 2)
        actual_area = psa.poly_area(vertices=vertices,
                                    radius=radius,
                                    cython=cython,
                                    discretizations=9000)
        assert_allclose(actual_area, expected_area)

    def test_quarter_hemisphere_area(self, cython, radius):
        # the area of 1/4 of a hemisphere should
        # be 1/8 the area of the entire sphere
        vertices = np.array([[-1,0,0],
                             [0,1,0],
                             [0,0,1]]) * radius
        expected_area = (np.pi * (radius ** 2)) / 2.
        actual_area = psa.poly_area(vertices=vertices,
                                    radius=radius,
                                    cython=cython,
                                    discretizations=9000)
        assert_allclose(actual_area, expected_area)

@pytest.mark.parametrize("cython, radius", [
    (None, 0),
    (1, 0),
    (None, -1e-10),
    (1, -1e-10),
    (None, -5),
    (1, -5),
    ])
def test_zero_radius_area(cython, radius):
    # an appropriate exception should be raised
    # for r <= 0.0
    vertices = np.array([[-1,0,0],
                         [1,0,0],
                         [0,0,1]]) * radius
    with pytest.raises(ValueError):
        psa.poly_area(vertices=vertices,
                      radius=radius,
                      cython=cython)

@pytest.mark.parametrize("cython, base, height", [
# stochastic failures observed with really high
# magnitudes here (e20), so scaled back to
# more reasonable magnitudes (e8)
(None, 1e-8, 1e-8),
(1, 1e-8, 1e-8),
(None, 1e8, 1e8),
(1, 1e8, 1e8),
(None, 0.5, 0.5),
(1, 0.5, 0.5),
])
class TestSimplePlanarTri(object):

    def test_planar_triangle_area(self, cython, base, height):
        # simple triangle area test
        # confirm that base * height / 2 result
        # is respected for a variety of base and
        # height values
        triangle_vertices = np.array([[0,0,0],
                                      [base,0,0],
                                      [base / 2.,height,0]])
        expected = 0.5 * base * height
        actual = psa.poly_area(vertices=triangle_vertices,
                               cython=cython)
        assert_allclose(actual, expected)

    def test_planar_triangle_area_reverse(self, cython, base, height):
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
                               cython=cython)
        assert_allclose(actual, expected)

@pytest.mark.parametrize("cython, radius", [
(None, 1e-20),
(1, 1e-20),
(None, 1e20),
(1, 1e20),
(None, 0.5),
(1, 0.5),
])
class TestRadianAreas(object):
    # compare spherical polygon surface areas
    # with values calculated using another well-known
    # approach
    # see: Weisstein, Eric W. "Spherical Polygon." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/SphericalPolygon.html
    # which cites: Beyer, W. H. CRC Standard Mathematical Tables, 28th ed. Boca Raton, FL: CRC Press, p. 131, 1987.
    # see also: Bevis and Cambereri (1987) Mathematical Geology 19: 335-346.
    # the above authors cite the angle-based equation used for reference
    # calcs here as well-known and stated frequently in handbooks of mathematics
    def _angle_area(self, sum_radian_angles, n_vertices, radius):
        # suface area of spherical polygon using
        # alternate approach
        area = (sum_radian_angles - (n_vertices - 2) * np.pi) * (radius ** 2)
        return area

    def test_double_octant_area_both_orders(self, cython, radius):
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
                                    cython=cython,
                                    discretizations=7000)
        actual_area_reverse = psa.poly_area(vertices=sample_vertices[::-1],
                                    radius=radius,
                                    cython=cython,
                                    discretizations=7000)
        assert_allclose(actual_area, expected_area)
        assert_allclose(actual_area_reverse, expected_area)

@pytest.mark.parametrize("cython, radius", [
(None, 1e-20),
(1, 1e-20),
(None, 1e20),
(1, 1e20),
(None, 0.5),
(1, 0.5),
])
class TestConvolutedAreas(object):
    # test more convoluted / tricky shapes
    # as input for surface area calculations

    def test_spherical_three_sixteen(self, cython, radius):
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
                                    cython=cython,
                                    discretizations=9000)
        actual_area_reverse = psa.poly_area(vertices=sample_vertices[::-1],
                                    radius=radius,
                                    cython=cython,
                                    discretizations=9000)
        assert_allclose(actual_area, expected_area)
        assert_allclose(actual_area_reverse, expected_area)

@pytest.mark.parametrize("cython", [None, 1])
class TestSplittingTriangles(object):
    # tests that leverage the splitting of
    # spherical triangles into 3 equal
    # area subtriangles

    def test_subtriangles_simple(self, cython):
        # start off with a spherical
        # triangle that covers 1/4
        # of a hemisphere and split into
        # 3 equal area subtriangles

        # verify that each subtriangle
        # has exactly 1/3 original triangle
        # expected area
        radius = 1.0
        vertices = np.array([[-1,0,0],
                             [0,1,0],
                             [0,0,1]])
        # divide by six for expected subtriangle areas
        # compare to expected area in test_quarter_hemisphere_area
        expected_sub_area = (np.pi * (radius ** 2)) / 6.

        # the original spherical triangle area
        # (before splitting into 3 subtriangles)
        original_tri_area = psa.poly_area(vertices=vertices,
                                    radius=radius,
                                    cython=cython,
                                    discretizations=9000)

        # find the central point D that splits the
        # spherical triangle into 3 equal area
        # subtriangles
        D = split_spherical_triangle.find_ternary_split_point(vertices,
                                                              1.0,
                                                              original_tri_area)

        vertices_subtriangle_1 = np.array([vertices[0],
                                           vertices[1],
                                           D])

        vertices_subtriangle_2 = np.array([vertices[1],
                                           vertices[2],
                                           D])

        vertices_subtriangle_3 = np.array([vertices[2],
                                           vertices[0],
                                           D])

        for subtriangle_vertices in [vertices_subtriangle_1,
                                     vertices_subtriangle_2,
                                     vertices_subtriangle_3]:
            actual_subtriangle_area = psa.poly_area(subtriangle_vertices,
                                                      1.0,
                                                      cython=cython,
                                                      discretizations=7000)
            assert_allclose(actual_subtriangle_area, expected_sub_area)

    @pytest.mark.parametrize("radius, iterations", [
    # NOTE: ~10 subdivisions seems to be
    # current tolerance limit for area calc
    (1.0, 1),
    (200.2, 1),
    (1.0, 10),
    (200.2, 10),
    (3.0, 7),
    (390.9872, 5),
    ])
    def test_subtriangles_iterative(self,
                                    radius,
                                    iterations,
                                    cython):
        # iteratively subdivide a simple
        # initial spherical triangle into 3 equal
        # area parts & verify expected
        # area of the final subtriangle
        vertices = np.array([[-radius,0,0],
                             [0,radius,0],
                             [0,0,radius]])

        original_tri_area = psa.poly_area(vertices=vertices,
                                    radius=radius,
                                    cython=cython,
                                    discretizations=12000)
        new_tri_area = original_tri_area

        for _ in range(iterations):
            D = split_spherical_triangle.find_ternary_split_point(vertices,
                                                                  radius,
                                                                  new_tri_area)
            new_vertices = np.array([vertices[0],
                                     vertices[1],
                                     D])
            vertices = new_vertices
            new_tri_area = psa.poly_area(vertices=vertices,
                                        radius=radius,
                                        cython=cython,
                                        discretizations=9000)

        expected_sub_area = original_tri_area / (3 ** iterations)
        actual_subtriangle_area = psa.poly_area(vertices=new_vertices,
                                                radius=radius,
                                                cython=cython,
                                                discretizations=9000)
        assert_allclose(actual_subtriangle_area,
                        expected_sub_area)

    @pytest.mark.parametrize("radius, iterations", [
    (1.0, 3),
    (200.2, 2),
    (1.0, 10),
    (200.2, 8),
    (3.0, 7),
    (390.9872, 5),
    ])
    def test_convoluted_polygons(self,
                                 radius,
                                 iterations,
                                 cython):
        # use spherical triangle 3-way splitting
        # to generate convoluted multi-vertex
        # spherical polygons
        vertices = np.array([[-radius,0,0],
                             [0,radius,0],
                             [0,0,radius]])

        original_area = psa.poly_area(vertices=vertices,
                                      radius=radius,
                                      cython=cython,
                                      discretizations=12000)

        # build up convoluted spherical polygons by
        # finding the 3-way triangle split point and
        # add this as a vertex to a newly-formed polygon

        expected_area = original_area
        new_tri_area = original_area
        new_subtri = vertices
        for iteration in range(iterations):
            # always split the 'most recent' subtriangle
            D = split_spherical_triangle.find_ternary_split_point(new_subtri,
                                                                  radius,
                                                                  new_tri_area)
            # the new polygon has an extra vertex between the
            # final and initial vertices on the first iteration,
            # and then always after the third vertex
            if iteration == 0:
                vertices = np.concatenate((vertices, D.reshape((1,3))))
            else:
                vertices = np.concatenate((vertices[:3,...],
                                           D.reshape((1,3)),
                                           vertices[3:,...]))

            new_subtri = np.array([vertices[1],
                                   vertices[2],
                                   D])

            # because we add the 3-way split point of the latest
            # spherical subtriangle, we adjust the expected area
            # of the growing convoluted spherical polygon to drop
            # by 1/3 of the surface area of the most recently-divided
            # subtriangle
            expected_area -= (new_tri_area * (1. / 3))

            new_tri_area = psa.poly_area(vertices=new_subtri,
                                        radius=radius,
                                        cython=cython,
                                        discretizations=92000)

        actual_polygon_area = psa.poly_area(vertices=vertices,
                                            radius=radius,
                                            cython=cython,
                                            discretizations=92000)

        assert_allclose(actual_polygon_area,
                        expected_area)
