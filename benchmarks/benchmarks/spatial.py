from __future__ import division, absolute_import, print_function

import sys
from scipy.spatial import cKDTree, KDTree
import numpy as np


class Build(object):
    params = [
        [(3,10000,1000), (8,10000,1000), (16,10000,1000)],
        ['KDTree', 'cKDTree'],
    ]
    param_names = ['(m, n, r)', 'class']
    goal_time = 0.5

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


class Query(object):
    params = [
        [(3,10000,1000), (8,10000,1000), (16,10000,1000)],
        ['KDTree', 'cKDTree', 'cKDTree_flat'],
    ]
    param_names = ['(m, n, r)', 'class']
    goal_time = 0.5

    @staticmethod
    def do_setup(self, mnr, cls_name):
        self.cls = KDTree if cls_name == 'KDTree' else cKDTree
        m, n, r = mnr

        np.random.seed(1234)
        self.data = np.concatenate((np.random.randn(n//2,m),
                                    np.random.randn(n-n//2,m)+np.ones(m)))

        self.queries = np.concatenate((np.random.randn(r//2,m),
                                       np.random.randn(r-r//2,m)+np.ones(m)))

        if cls_name == 'cKDTree_flat':
            self.T = self.cls(self.data, leafsize=n)
        else:
            self.T = self.cls(self.data)

    def setup(self, mnr, cls_name):
        Query.do_setup(self, mnr, cls_name)

    def time_query(self, mnr, cls_name):
        """
        Querying kd-tree
        dim | # points | # queries |  KDTree  | cKDTree | flat cKDTree
        """
        self.T.query(self.queries)


class Radius(object):
    params = [
        [(3,10000,1000)],
        [0.2, 0.5],
        ['KDTree', 'cKDTree', 'cKDTree_flat'],
    ]
    param_names = ['(m, n, r)', 'probe radius', 'class']
    goal_time = 0.5

    def __init__(self):
        self.time_query_pairs.__func__.params = list(self.params)
        self.time_query_pairs.__func__.params[0] = [(3,1000,30),
                                                    (8,1000,30),
                                                    (16,1000,30)]

    def setup(self, mnr, probe_radius, cls_name):
        Query.do_setup(self, mnr, cls_name)

    def time_query_ball_point(self, mnr, probe_radius, cls_name):
        self.T.query_ball_point(self.queries, probe_radius)

    def time_query_pairs(self, mnr, probe_radius, cls_name):
        self.T.query_pairs(probe_radius)


class Neighbors(object):
    params = [
        [(3,1000,1000),
         (8,1000,1000),
         (16,1000,1000)],
        [0.2, 0.5],
        ['KDTree', 'cKDTree'],
    ]
    param_names = ['(m, n1, n2)', 'probe radius', 'class']
    goal_time = 0.5
    timeout = 120

    def setup(self, mn1n2, probe_radius, cls_str):
        m, n1, n2 = mn1n2

        cls = KDTree if cls_str == 'KDTree' else cKDTree

        data1 = np.concatenate((np.random.randn(n1//2,m),
                                np.random.randn(n1-n1//2,m)+np.ones(m)))
        data2 = np.concatenate((np.random.randn(n2//2,m),
                                np.random.randn(n2-n2//2,m)+np.ones(m)))

        self.T1 = cls(data1)
        self.T2 = cls(data2)

    def time_sparse_distance_matrix(self, mn1n2, probe_radius, cls_str):
        self.T1.sparse_distance_matrix(self.T2, probe_radius)

    def time_count_neighbors(self, mn1n2, probe_radius, cls_str):
        """
        Count neighbors kd-tree
        dim | # points T1 | # points T2 | probe radius |  KDTree  | cKDTree
        """
        self.T1.count_neighbors(self.T2, probe_radius)
