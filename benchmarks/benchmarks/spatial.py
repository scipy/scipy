from __future__ import division, absolute_import, print_function

import numpy as np

try:
    from scipy.spatial import cKDTree, KDTree
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
BOX_SIZES = [None, 1.0]

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
    ]
    param_names = ['(m, n1, n2)', 'p', 'probe radius', 'boxsize', 'leafsize']

    def setup(self, mn1n2, p, probe_radius, boxsize, leafsize):
        m, n1, n2 = mn1n2

        self.data1 = np.random.uniform(size=(n1, m))
        self.data2 = np.random.uniform(size=(n2, m))

        self.T1 = cKDTree(self.data1, boxsize=boxsize, leafsize=leafsize)
        self.T2 = cKDTree(self.data2, boxsize=boxsize, leafsize=leafsize)

    def time_sparse_distance_matrix(self, mn1n2, p, probe_radius, boxsize, leafsize):
        self.T1.sparse_distance_matrix(self.T2, probe_radius, p=p)

    def time_count_neighbors(self, mn1n2, p, probe_radius, boxsize, leafsize):
        """
        Count neighbors kd-tree
        dim | # points T1 | # points T2 | p | probe radius |  BoxSize | LeafSize
        """
        self.T1.count_neighbors(self.T2, probe_radius, p=p)
