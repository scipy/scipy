from __future__ import print_function
import numpy as np
import itertools
from numpy.testing import (assert_equal,
                           assert_almost_equal,
                           assert_array_equal,
                           assert_array_almost_equal)
import warnings
from pytest import raises as assert_raises
from pytest import warns as assert_warns
from scipy.spatial import SphericalVoronoi, distance
from scipy.spatial import _spherical_voronoi as spherical_voronoi
from scipy._lib._numpy_compat import suppress_warnings


class TestSphericalVoronoi(object):

    def setup_method(self):
        self.points = np.array([
            [-0.78928481, -0.16341094, 0.59188373],
            [-0.66839141, 0.73309634, 0.12578818],
            [0.32535778, -0.92476944, -0.19734181],
            [-0.90177102, -0.03785291, -0.43055335],
            [0.71781344, 0.68428936, 0.12842096],
            [-0.96064876, 0.23492353, -0.14820556],
            [0.73181537, -0.22025898, -0.6449281],
            [0.79979205, 0.54555747, 0.25039913]]
        )

        # Issue #9386
        self.hemisphere_points = np.array([
            [0.88610999, -0.42383021, 0.18755541],
            [0.51980039, -0.72622668, 0.4498915],
            [0.56540011, -0.81629197, -0.11827989],
            [0.69659682, -0.69972598, 0.15854467]])

        # Issue #8859
        phi = np.linspace(0, 2 * np.pi, 10, endpoint=False)    # azimuth angle
        theta = np.linspace(0.001, np.pi * 0.4, 5)    # polar angle
        theta = theta[np.newaxis, :].T

        phiv, thetav = np.meshgrid(phi, theta)
        phiv = np.reshape(phiv, (50, 1))
        thetav = np.reshape(thetav, (50, 1))

        x = np.cos(phiv) * np.sin(thetav)
        y = np.sin(phiv) * np.sin(thetav)
        z = np.cos(thetav)
        self.hemisphere_points2 = np.concatenate([x, y, z], axis=1)

    def test_constructor(self):
        center = np.array([1, 2, 3])
        radius = 2
        s1 = SphericalVoronoi(self.points)
        # user input checks in SphericalVoronoi now require
        # the radius / center to match the generators so adjust
        # accordingly here
        s2 = SphericalVoronoi(self.points * radius, radius)
        s3 = SphericalVoronoi(self.points + center, center=center)
        s4 = SphericalVoronoi(self.points * radius + center, radius, center)
        assert_array_equal(s1.center, np.array([0, 0, 0]))
        assert_equal(s1.radius, 1)
        assert_array_equal(s2.center, np.array([0, 0, 0]))
        assert_equal(s2.radius, 2)
        assert_array_equal(s3.center, center)
        assert_equal(s3.radius, 1)
        assert_array_equal(s4.center, center)
        assert_equal(s4.radius, radius)

    def test_vertices_regions_translation_invariance(self):
        sv_origin = SphericalVoronoi(self.points)
        center = np.array([1, 1, 1])
        sv_translated = SphericalVoronoi(self.points + center, center=center)
        assert_array_equal(sv_origin.regions, sv_translated.regions)
        assert_array_almost_equal(sv_origin.vertices + center,
                                  sv_translated.vertices)

    def test_vertices_regions_scaling_invariance(self):
        sv_unit = SphericalVoronoi(self.points)
        sv_scaled = SphericalVoronoi(self.points * 2, 2)
        assert_array_equal(sv_unit.regions, sv_scaled.regions)
        assert_array_almost_equal(sv_unit.vertices * 2,
                                  sv_scaled.vertices)

    def test_old_radius_api(self):
        sv_unit = SphericalVoronoi(self.points, radius=1)
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning, "`radius` is `None`")
            sv = SphericalVoronoi(self.points, None)
            assert_array_almost_equal(sv_unit.vertices, sv.vertices)

    def test_old_radius_api_warning(self):
        with assert_warns(DeprecationWarning):
            sv = SphericalVoronoi(self.points, None)

    def test_old_center_api(self):
        sv_unit = SphericalVoronoi(self.points, radius=1, center=(0, 0, 0))
        with suppress_warnings() as sup:
            sup.filter(DeprecationWarning, "`radius` is `None`")
            sup.filter(DeprecationWarning, "`center` is `None`")
            sv = SphericalVoronoi(self.points, None, None)
            assert_array_almost_equal(sv_unit.vertices, sv.vertices)

    def test_old_center_api_warning(self):
        with assert_warns(DeprecationWarning):
            sv = SphericalVoronoi(self.points, None, None)

    def test_sort_vertices_of_regions(self):
        sv = SphericalVoronoi(self.points)
        unsorted_regions = sv.regions
        sv.sort_vertices_of_regions()
        assert_array_equal(sorted(sv.regions), sorted(unsorted_regions))

    def test_sort_vertices_of_regions_flattened(self):
        expected = sorted([[0, 6, 5, 2, 3], [2, 3, 10, 11, 8, 7], [0, 6, 4, 1], [4, 8,
            7, 5, 6], [9, 11, 10], [2, 7, 5], [1, 4, 8, 11, 9], [0, 3, 10, 9,
                1]])
        expected = list(itertools.chain(*sorted(expected)))
        sv = SphericalVoronoi(self.points)
        sv.sort_vertices_of_regions()
        actual = list(itertools.chain(*sorted(sv.regions)))
        assert_array_equal(actual, expected)

    def test_num_vertices(self):
        # for any n >= 3, a spherical Voronoi diagram has 2n - 4
        # vertices; this is a direct consequence of Euler's formula
        # as explained by Dinis and Mamede (2010) Proceedings of the
        # 2010 International Symposium on Voronoi Diagrams in Science
        # and Engineering
        sv = SphericalVoronoi(self.points)
        expected = self.points.shape[0] * 2 - 4
        actual = sv.vertices.shape[0]
        assert_equal(actual, expected)

    def test_voronoi_circles(self):
        sv = spherical_voronoi.SphericalVoronoi(self.points)
        for vertex in sv.vertices:
            distances = distance.cdist(sv.points,np.array([vertex]))
            closest = np.array(sorted(distances)[0:3])
            assert_almost_equal(closest[0], closest[1], 7, str(vertex))
            assert_almost_equal(closest[0], closest[2], 7, str(vertex))

    def test_duplicate_point_handling(self):
        # an exception should be raised for degenerate generators
        # related to Issue# 7046
        self.degenerate = np.concatenate((self.points, self.points))
        with assert_raises(ValueError):
            sv = spherical_voronoi.SphericalVoronoi(self.degenerate)

    def test_incorrect_radius_handling(self):
        # an exception should be raised if the radius provided
        # cannot possibly match the input generators
        with assert_raises(ValueError):
            sv = spherical_voronoi.SphericalVoronoi(self.points,
                                                    radius=0.98)

    def test_incorrect_center_handling(self):
        # an exception should be raised if the center provided
        # cannot possibly match the input generators
        with assert_raises(ValueError):
            sv = spherical_voronoi.SphericalVoronoi(self.points,
                                                    center=[0.1,0,0])

    def test_single_hemisphere_handling(self):
        # Test solution of Issues #9386, #8859

        for points in [self.hemisphere_points, self.hemisphere_points2]:
            sv = SphericalVoronoi(points)
            triangles = sv._tri.points[sv._tri.simplices]
            dots = np.einsum('ij,ij->i', sv.vertices, triangles[:, 0])
            circumradii = np.arccos(np.clip(dots, -1, 1))
            assert np.max(circumradii) > np.pi / 2
