from __future__ import print_function
import numpy as np
from numpy.testing import (TestCase,
                           assert_array_equal,
                           assert_array_almost_equal)
from scipy.spatial import spherical_voronoi


class TestCircumcenters(TestCase):

    def test_circumcenters(self):
        tetrahedrons = np.array([
            [[1, 2, 3],
             [-1.1, -2.1, -3.1],
             [-1.2, 2.2, 3.2],
             [-1.3, -2.3, 3.3]],
            [[10, 20, 30],
             [-10.1, -20.1, -30.1],
             [-10.2, 20.2, 30.2],
             [-10.3, -20.3, 30.3]]
        ])

        result = spherical_voronoi.calc_circumcenters(tetrahedrons)

        expected = [
            [-0.5680861153262529, -0.133279590288315, 0.1843323216995444],
            [-0.5965330784014926, -0.1480377040397778, 0.1981967854886021]
        ]

        assert_array_almost_equal(result, expected)


class TestProjectToSphere(TestCase):

    def test_unit_sphere(self):
        points = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        center = np.array([0, 0, 0])
        radius = 1
        projected = spherical_voronoi.project_to_sphere(points, center, radius)
        assert_array_almost_equal(points, projected)

    def test_scaled_points(self):
        points = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        center = np.array([0, 0, 0])
        radius = 1
        scaled = points * 2
        projected = spherical_voronoi.project_to_sphere(scaled, center, radius)
        assert_array_almost_equal(points, projected)

    def test_translated_sphere(self):
        points = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        center = np.array([1, 2, 3])
        translated = points + center
        radius = 1
        projected = spherical_voronoi.project_to_sphere(translated, center,
                                                        radius)
        assert_array_almost_equal(translated, projected)


class TestSphericalVoronoi(TestCase):

    points_unsymmetric = np.array([
        [-0.78928481, -0.16341094, 0.59188373],
        [-0.66839141, 0.73309634, 0.12578818],
        [0.32535778, -0.92476944, -0.19734181],
        [-0.90177102, -0.03785291, -0.43055335],
        [0.71781344, 0.68428936, 0.12842096],
        [-0.96064876, 0.23492353, -0.14820556],
        [0.73181537, -0.22025898, -0.6449281],
        [0.79979205, 0.54555747, 0.25039913]]
    )

    def test_constructor(self):
        generators = TestSphericalVoronoi.points_unsymmetric
        center = np.array([1, 2, 3])
        radius = 2
        s1 = spherical_voronoi.SphericalVoronoi(generators)
        s2 = spherical_voronoi.SphericalVoronoi(generators, radius)
        s3 = spherical_voronoi.SphericalVoronoi(generators, None, center)
        s4 = spherical_voronoi.SphericalVoronoi(generators, radius, center)
        assert_array_equal(s1.center, np.array([0, 0, 0]))
        self.assertEqual(s1.radius, 1)
        assert_array_equal(s2.center, np.array([0, 0, 0]))
        self.assertEqual(s2.radius, 2)
        assert_array_equal(s3.center, center)
        self.assertEqual(s3.radius, 1)
        assert_array_equal(s4.center, center)
        self.assertEqual(s4.radius, radius)

    def test_vertices_regions_translation_invariance(self):
        points = TestSphericalVoronoi.points_unsymmetric
        sv_origin = spherical_voronoi.SphericalVoronoi(points)
        center = np.array([1, 1, 1])
        points += center
        sv_translated = spherical_voronoi.SphericalVoronoi(points, None, center)
        assert_array_equal(sv_origin.regions, sv_translated.regions)
        assert_array_almost_equal(sv_origin.vertices + center,
                                  sv_translated.vertices)

    def test_vertices_regions_scaling_invariance(self):
        points = TestSphericalVoronoi.points_unsymmetric
        sv_unit = spherical_voronoi.SphericalVoronoi(points)
        points *= 2
        sv_scaled = spherical_voronoi.SphericalVoronoi(points, 2)
        assert_array_equal(sv_unit.regions, sv_scaled.regions)
        assert_array_almost_equal(sv_unit.vertices * 2,
                                  sv_scaled.vertices)

    def test_sort_vertices_of_regions(self):
        generators = TestSphericalVoronoi.points_unsymmetric
        sv = spherical_voronoi.SphericalVoronoi(generators)
        unsorted_regions = sv.regions
        sv.sort_vertices_of_regions()
        assert_array_equal(sorted(sv.regions), sorted(unsorted_regions))
