from __future__ import division, absolute_import, print_function

import numpy as np

try:
    from scipy.spatial import (cKDTree, KDTree, SphericalVoronoi, distance,
    ConvexHull, Voronoi)
except ImportError:
    pass

from .common import Benchmark


class Build(Benchmark):
    params = [
        [(3,10000,1000), (8,10000,1000), (16,10000,1000)],
        ['KDTree', 'cKDTree'],
    ]
    param_names = ['(m, n, r)', 'class']

    def setup(self, mnr, cls_name):
        self.cls = KDTree if cls_name == 'KDTree' else cKDTree
        m, n, r = mnr

        np.random.seed(1234)
        self.data = np.concatenate((np.random.randn(n//2,m),
                                    np.random.randn(n-n//2,m)+np.ones(m)))

        self.queries = np.concatenate((np.random.randn(r//2,m),
                                       np.random.randn(r-r//2,m)+np.ones(m)))

    def time_build(self, mnr, cls_name):
        """
        Constructing kd-tree
        =======================
        dim | # points |  time
        """
        m, n, r = mnr
        if cls_name == 'cKDTree_flat':
            self.T = self.cls(self.data, leafsize=n)
        else:
            self.cls(self.data)

LEAF_SIZES = [8, 128]
BOX_SIZES = [None, 0.0, 1.0]

class Query(Benchmark):
    params = [
        [(3,10000,1000), (8,10000,1000), (16,10000,1000)],
        [1, 2, np.inf],
        BOX_SIZES, LEAF_SIZES,
    ]
    param_names = ['(m, n, r)', 'p', 'boxsize', 'leafsize']

    @staticmethod
    def do_setup(self, mnr, p, boxsize, leafsize):
        m, n, r = mnr

        np.random.seed(1234)

        self.data = np.random.uniform(size=(n, m))
        self.queries = np.random.uniform(size=(r, m))

        self.T = cKDTree(self.data, leafsize=leafsize, boxsize=boxsize)

    def setup(self, mnr, p, boxsize, leafsize):
        Query.do_setup(self, mnr, p, boxsize, leafsize)

    def time_query(self, mnr, p, boxsize, leafsize):
        """
        Querying kd-tree
        dim | # points | # queries |  KDTree  | cKDTree | flat cKDTree
        """
        self.T.query(self.queries, p=p)


class Radius(Benchmark):
    params = [
        [(3,10000,1000)],
        [1, 2, np.inf],
        [0.2, 0.5],
        BOX_SIZES, LEAF_SIZES,
    ]
    param_names = ['(m, n, r)', 'p', 'probe radius', 'boxsize', 'leafsize']

    def __init__(self):
        self.time_query_pairs.__func__.params = list(self.params)
        self.time_query_pairs.__func__.params[0] = [(3,1000,30),
                                                    (8,1000,30),
                                                    (16,1000,30)]

    def setup(self, mnr, p, probe_radius, boxsize, leafsize):
        Query.do_setup(self, mnr, p, boxsize, leafsize)

    def time_query_ball_point(self, mnr, p, probe_radius, boxsize, leafsize):
        self.T.query_ball_point(self.queries, probe_radius, p=p)

    def time_query_pairs(self, mnr, p, probe_radius, boxsize, leafsize):
        self.T.query_pairs(probe_radius, p=p)


class Neighbors(Benchmark):
    params = [
        [(3,1000,1000),
         (8,1000,1000),
         (16,1000,1000)],
        [1, 2, np.inf],
        [0.2, 0.5],
        BOX_SIZES, LEAF_SIZES,
        ['cKDTree', 'cKDTree_weighted'],
    ]
    param_names = ['(m, n1, n2)', 'p', 'probe radius', 'boxsize', 'leafsize', 'cls']

    def setup(self, mn1n2, p, probe_radius, boxsize, leafsize, cls):
        m, n1, n2 = mn1n2

        self.data1 = np.random.uniform(size=(n1, m))
        self.data2 = np.random.uniform(size=(n2, m))

        self.w1 = np.ones(n1)
        self.w2 = np.ones(n2)

        self.T1 = cKDTree(self.data1, boxsize=boxsize, leafsize=leafsize)
        self.T2 = cKDTree(self.data2, boxsize=boxsize, leafsize=leafsize)

    def time_sparse_distance_matrix(self, mn1n2, p, probe_radius, boxsize, leafsize, cls):
        self.T1.sparse_distance_matrix(self.T2, probe_radius, p=p)

    def time_count_neighbors(self, mn1n2, p, probe_radius, boxsize, leafsize, cls):
        """
        Count neighbors kd-tree
        dim | # points T1 | # points T2 | p | probe radius |  BoxSize | LeafSize | cls
        """

        if cls != 'cKDTree_weighted':
            self.T1.count_neighbors(self.T2, probe_radius, p=p)
        else:
            self.T1.count_neighbors(self.T2, probe_radius, weights=(self.w1, self.w2), p=p)

class CNeighbors(Benchmark):
    params = [
        [
          (2,1000,1000),
          (8,1000,1000),
          (16,1000,1000)
        ],
        [2, 10, 100, 400, 1000],
    ]
    param_names = ['(m, n1, n2)', 'Nr']

    def setup(self, mn1n2, Nr):
        m, n1, n2 = mn1n2

        data1 = np.random.uniform(size=(n1, m))
        data2 = np.random.uniform(size=(n2, m))
        self.w1 = np.ones(len(data1))
        self.w2 = np.ones(len(data2))
 
        self.T1d = cKDTree(data1, leafsize=1)
        self.T2d = cKDTree(data2, leafsize=1)
        self.T1s = cKDTree(data1, leafsize=8)
        self.T2s = cKDTree(data2, leafsize=8)
        self.r = np.linspace(0, 0.5, Nr)

    def time_count_neighbors_deep(self, mn1n2, Nr):
        """
        Count neighbors for a very deep kd-tree
        dim | # points T1 | # points T2 | Nr
        """
        self.T1d.count_neighbors(self.T2d, self.r)

    def time_count_neighbors_shallow(self, mn1n2, Nr):
        """
        Count neighbors for a shallow kd-tree
        dim | # points T1 | # points T2 | Nr
        """
        self.T1s.count_neighbors(self.T2s, self.r)

def generate_spherical_points(num_points):
        # generate uniform points on sphere (see:
        # http://stackoverflow.com/a/23785326/2942522)
        np.random.seed(123)
        points = np.random.normal(size=(num_points, 3))
        points /= np.linalg.norm(points, axis=1)[:, np.newaxis]
        return points

class SphericalVor(Benchmark):
    params = [10, 100, 1000, 5000, 10000]
    param_names = ['num_points']

    def setup(self, num_points):
        self.points = generate_spherical_points(num_points)

    def time_spherical_voronoi_calculation(self, num_points):
        """Perform spherical Voronoi calculation, but not the sorting of
        vertices in the Voronoi polygons.
        """
        SphericalVoronoi(self.points, radius=1, center=np.zeros(3))

class SphericalVorSort(Benchmark):
    params = [10, 100, 1000, 5000, 10000]
    param_names = ['num_points']

    def setup(self, num_points):
        self.points = generate_spherical_points(num_points)
        self.sv = SphericalVoronoi(self.points, radius=1,
                                   center=np.zeros(3))

    def time_spherical_polygon_vertex_sorting(self, num_points):
        """Time the vertex sorting operation in the Spherical Voronoi
        code.
        """
        self.sv.sort_vertices_of_regions()

class Cdist(Benchmark):
    params = ([10, 100, 1000], ['euclidean', 'minkowski', 'cityblock',
    'seuclidean', 'sqeuclidean', 'cosine', 'correlation', 'hamming', 'jaccard',
    'chebyshev', 'canberra', 'braycurtis', 'mahalanobis', 'yule', 'dice',
    'kulsinski', 'rogerstanimoto', 'russellrao', 'sokalmichener',
    'sokalsneath', 'wminkowski'])
    param_names = ['num_points', 'metric']

    def setup(self, num_points, metric):
        np.random.seed(123)
        self.points = np.random.random_sample((num_points, 3))

    def time_cdist(self, num_points, metric):
        """Time scipy.spatial.distance.cdist over a range of input data
        sizes and metrics.
        """
        distance.cdist(self.points, self.points, metric)


class ConvexHullBench(Benchmark):
    params = ([10, 100, 1000, 5000], [True, False])
    param_names = ['num_points', 'incremental']

    def setup(self, num_points, incremental):
        np.random.seed(123)
        self.points = np.random.random_sample((num_points, 3))

    def time_convex_hull(self, num_points, incremental):
        """Time scipy.spatial.ConvexHull over a range of input data sizes
        and settings.
        """
        ConvexHull(self.points, incremental)


class VoronoiBench(Benchmark):
    params = ([10, 100, 1000, 5000, 10000], [False, True])
    param_names = ['num_points', 'furthest_site']

    def setup(self, num_points, furthest_site):
        np.random.seed(123)
        self.points = np.random.random_sample((num_points, 3))

    def time_voronoi_calculation(self, num_points, furthest_site):
        """Time conventional Voronoi diagram calculation."""
        Voronoi(self.points, furthest_site=furthest_site)

class Hausdorff(Benchmark):
    params = [10, 100, 1000]
    param_names = ['num_points']

    def setup(self, num_points):
        np.random.seed(123)
        self.points1 = np.random.random_sample((num_points, 3))
        np.random.seed(71890)
        self.points2 = np.random.random_sample((num_points, 3))

    def time_directed_hausdorff(self, num_points):
        # time directed_hausdorff code in 3 D
        distance.directed_hausdorff(self.points1, self.points2)

