#! /usr/bin/env python

# David Cournapeau
# Last Change: Wed Nov 05 07:00 PM 2008 J
from __future__ import division, print_function, absolute_import

import os.path
import warnings

import numpy as np
from numpy.testing import (assert_array_equal, assert_array_almost_equal,
    TestCase, run_module_suite, assert_raises, assert_allclose, assert_equal,
    assert_)

from scipy.cluster.vq import (kmeans, kmeans2, py_vq, py_vq2, vq, whiten,
    ClusterError)
from scipy.cluster import _vq

# Optional:
# import modules that are located in the same directory as this file.
DATAFILE1 = os.path.join(os.path.dirname(__file__), "data.txt")

# Global data
X = np.array([[3.0, 3], [4, 3], [4, 2],
               [9, 2], [5, 1], [6, 2], [9, 4],
               [5, 2], [5, 4], [7, 4], [6, 5]])

CODET1 = np.array([[3.0000, 3.0000],
                   [6.2000, 4.0000],
                   [5.8000, 1.8000]])

CODET2 = np.array([[11.0/3, 8.0/3],
                   [6.7500, 4.2500],
                   [6.2500, 1.7500]])

LABEL1 = np.array([0, 1, 2, 2, 2, 2, 1, 2, 1, 1, 1])


class TestWhiten(TestCase):
    def test_whiten(self):
        obs = np.array([[0.98744510, 0.82766775],
                        [0.62093317, 0.19406729],
                        [0.87545741, 0.00735733],
                        [0.85124403, 0.26499712],
                        [0.45067590, 0.45464607]])
        result = np.array([[5.08738849, 2.97091878],
                           [3.19909255, 0.69660580],
                           [4.51041982, 0.02640918],
                           [4.38567074, 0.95120889],
                           [2.32191480, 1.63195503]])
        assert_allclose(whiten(obs), result, rtol=1e-5)

    def test_whiten_zero_std(self):
        obs = np.array([[0., 1., 0.74109533],
                        [0., 1., 0.34243798],
                        [0., 1., 0.96785929]])
        result = np.array([[0., 1.0, 2.86666544],
                           [0., 1.0, 1.32460034],
                           [0., 1.0, 3.74382172]])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            assert_allclose(whiten(obs), result, rtol=1e-5)
            assert_equal(len(w), 1)
            assert_(issubclass(w[-1].category, RuntimeWarning))


class TestVq(TestCase):
    def test_py_vq(self):
        initc = np.concatenate(([[X[0]], [X[1]], [X[2]]]))
        label1 = py_vq(X, initc)[0]
        assert_array_equal(label1, LABEL1)

    def test_py_vq2(self):
        initc = np.concatenate(([[X[0]], [X[1]], [X[2]]]))
        label1 = py_vq2(X, initc)[0]
        assert_array_equal(label1, LABEL1)

    def test_vq(self):
        initc = np.concatenate(([[X[0]], [X[1]], [X[2]]]))
        label1, dist = _vq.vq(X, initc)
        assert_array_equal(label1, LABEL1)
        tlabel1, tdist = vq(X, initc)

    # def test_py_vq_1d(self):
    #     """Test special rank 1 vq algo, python implementation."""
    #     data = X[:, 0]
    #     initc = data[:3]
    #     a, b = _py_vq_1d(data, initc)
    #     ta, tb = py_vq(data[:, np.newaxis], initc[:, np.newaxis])
    #     assert_array_equal(a, ta)
    #     assert_array_equal(b, tb)

    def test_vq_1d(self):
        # Test special rank 1 vq algo, python implementation.
        data = X[:, 0]
        initc = data[:3]
        a, b = _vq.vq(data, initc)
        ta, tb = py_vq(data[:, np.newaxis], initc[:, np.newaxis])
        assert_array_equal(a, ta)
        assert_array_equal(b, tb)

    def test__vq_sametype(self):
        a = np.array([1.0, 2.0], dtype=np.float64)
        b = a.astype(np.float32)
        assert_raises(TypeError, _vq.vq, a, b)

    def test__vq_invalid_type(self):
        a = np.array([1, 2], dtype=np.int)
        assert_raises(TypeError, _vq.vq, a, a)

    def test_vq_large_nfeat(self):
        X = np.random.rand(20, 20)
        code_book = np.random.rand(3, 20)

        codes0, dis0 = _vq.vq(X, code_book)
        codes1, dis1 = py_vq(X, code_book)
        assert_allclose(dis0, dis1, 1e-5)
        assert_array_equal(codes0, codes1)

        X = X.astype(np.float32)
        code_book = code_book.astype(np.float32)

        codes0, dis0 = _vq.vq(X, code_book)
        codes1, dis1 = py_vq(X, code_book)
        assert_allclose(dis0, dis1, 1e-5)
        assert_array_equal(codes0, codes1)

    def test_vq_large_features(self):
        X = np.random.rand(10, 5) * 1000000
        code_book = np.random.rand(2, 5) * 1000000

        codes0, dis0 = _vq.vq(X, code_book)
        codes1, dis1 = py_vq(X, code_book)
        assert_allclose(dis0, dis1, 1e-5)
        assert_array_equal(codes0, codes1)


class TestKMean(TestCase):
    def test_large_features(self):
        # Generate a data set with large values, and run kmeans on it to
        # (regression for 1077).
        d = 300
        n = 100

        m1 = np.random.randn(d)
        m2 = np.random.randn(d)
        x = 10000 * np.random.randn(n, d) - 20000 * m1
        y = 10000 * np.random.randn(n, d) + 20000 * m2

        data = np.empty((x.shape[0] + y.shape[0], d), np.double)
        data[:x.shape[0]] = x
        data[x.shape[0]:] = y

        kmeans(data, 2)

    def test_kmeans_simple(self):
        initc = np.concatenate(([[X[0]], [X[1]], [X[2]]]))
        code = initc.copy()
        code1 = kmeans(X, code, iter=1)[0]

        assert_array_almost_equal(code1, CODET2)

    def test_kmeans_lost_cluster(self):
        # This will cause kmean to have a cluster with no points.
        data = np.fromfile(DATAFILE1, sep=", ")
        data = data.reshape((200, 2))
        initk = np.array([[-1.8127404, -0.67128041],
                         [2.04621601, 0.07401111],
                         [-2.31149087,-0.05160469]])

        kmeans(data, initk)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            kmeans2(data, initk, missing='warn')

        assert_raises(ClusterError, kmeans2, data, initk, missing='raise')

    def test_kmeans2_simple(self):
        initc = np.concatenate(([[X[0]], [X[1]], [X[2]]]))
        code = initc.copy()
        code1 = kmeans2(X, code, iter=1)[0]
        code2 = kmeans2(X, code, iter=2)[0]

        assert_array_almost_equal(code1, CODET1)
        assert_array_almost_equal(code2, CODET2)

    def test_kmeans2_rank1(self):
        data = np.fromfile(DATAFILE1, sep=", ")
        data = data.reshape((200, 2))
        data1 = data[:, 0]

        initc = data1[:3]
        code = initc.copy()
        kmeans2(data1, code, iter=1)[0]
        kmeans2(data1, code, iter=2)[0]

    def test_kmeans2_rank1_2(self):
        data = np.fromfile(DATAFILE1, sep=", ")
        data = data.reshape((200, 2))
        data1 = data[:, 0]
        kmeans2(data1, 2, iter=1)

    def test_kmeans2_init(self):
        data = np.fromfile(DATAFILE1, sep=", ")
        data = data.reshape((200, 2))

        kmeans2(data, 3, minit='points')
        kmeans2(data[:, :1], 3, minit='points')  # special case (1-D)

        # minit='random' can give warnings, filter those
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore',
                        message="One of the clusters is empty. Re-run")
            kmeans2(data, 3, minit='random')
            kmeans2(data[:, :1], 3, minit='random')  # special case (1-D)

    def test_kmeans2_empty(self):
        # Regression test for gh-1032.
        assert_raises(ValueError, kmeans2, [], 2)

    def test_kmeans_0k(self):
        # Regression test for gh-1073: fail when k arg is 0.
        assert_raises(ValueError, kmeans, X, 0)
        assert_raises(ValueError, kmeans2, X, 0)
        assert_raises(ValueError, kmeans2, X, np.array([]))


if __name__ == "__main__":
    run_module_suite()
