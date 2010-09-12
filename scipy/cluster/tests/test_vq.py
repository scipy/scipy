#! /usr/bin/env python

# David Cournapeau
# Last Change: Wed Nov 05 07:00 PM 2008 J

import os.path
import warnings

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, \
        TestCase, run_module_suite

from scipy.cluster.vq import kmeans, kmeans2, py_vq, py_vq2, vq, ClusterError
try:
    from scipy.cluster import _vq
    TESTC=True
except ImportError:
    print "== Error while importing _vq, not testing C imp of vq =="
    TESTC=False

#Optional:
# import modules that are located in the same directory as this file.
DATAFILE1 = os.path.join(os.path.dirname(__file__), "data.txt")

# Global data
X   = np.array([[3.0, 3], [4, 3], [4, 2],
               [9, 2], [5, 1], [6, 2], [9, 4],
               [5, 2], [5, 4], [7, 4], [6, 5]])

CODET1  = np.array([[3.0000, 3.0000],
                   [6.2000, 4.0000],
                   [5.8000, 1.8000]])

CODET2  = np.array([[11.0/3, 8.0/3],
                   [6.7500, 4.2500],
                   [6.2500, 1.7500]])

LABEL1  = np.array([0, 1, 2, 2, 2, 2, 1, 2, 1, 1, 1])

class TestVq(TestCase):
    def test_py_vq(self):
        initc = np.concatenate(([[X[0]], [X[1]], [X[2]]]))
        code = initc.copy()
        label1 = py_vq(X, initc)[0]
        assert_array_equal(label1, LABEL1)

    def test_py_vq2(self):
        initc = np.concatenate(([[X[0]], [X[1]], [X[2]]]))
        code = initc.copy()
        label1 = py_vq2(X, initc)[0]
        assert_array_equal(label1, LABEL1)

    def test_vq(self):
        initc = np.concatenate(([[X[0]], [X[1]], [X[2]]]))
        code = initc.copy()
        if TESTC:
            label1, dist = _vq.vq(X, initc)
            assert_array_equal(label1, LABEL1)
            tlabel1, tdist = vq(X, initc)
        else:
            print "== not testing C imp of vq =="

    #def test_py_vq_1d(self):
    #    """Test special rank 1 vq algo, python implementation."""
    #    data = X[:, 0]
    #    initc = data[:3]
    #    code = initc.copy()
    #    a, b = _py_vq_1d(data, initc)
    #    ta, tb = py_vq(data[:, np.newaxis], initc[:, np.newaxis])
    #    assert_array_equal(a, ta)
    #    assert_array_equal(b, tb)

    def test_vq_1d(self):
        """Test special rank 1 vq algo, python implementation."""
        data = X[:, 0]
        initc = data[:3]
        code = initc.copy()
        if TESTC:
            a, b = _vq.vq(data, initc)
            ta, tb = py_vq(data[:, np.newaxis], initc[:, np.newaxis])
            assert_array_equal(a, ta)
            assert_array_equal(b, tb)
        else:
            print "== not testing C imp of vq (rank 1) =="

class TestKMean(TestCase):
    def test_large_features(self):
        # Generate a data set with large values, and run kmeans on it to
        # (regression for 1077).
        d = 300
        n = 1e2

        m1 = np.random.randn(d)
        m2 = np.random.randn(d)
        x = 10000 * np.random.randn(n, d) - 20000 * m1
        y = 10000 * np.random.randn(n, d) + 20000 * m2

        data = np.empty((x.shape[0] + y.shape[0], d), np.double)
        data[:x.shape[0]] = x
        data[x.shape[0]:] = y

        res = kmeans(data, 2)
    def test_kmeans_simple(self):
        initc = np.concatenate(([[X[0]], [X[1]], [X[2]]]))
        code = initc.copy()
        code1 = kmeans(X, code, iter = 1)[0]

        assert_array_almost_equal(code1, CODET2)

    def test_kmeans_lost_cluster(self):
        """This will cause kmean to have a cluster with no points."""
        data = np.fromfile(open(DATAFILE1), sep = ", ")
        data = data.reshape((200, 2))
        initk = np.array([[-1.8127404, -0.67128041],
                         [ 2.04621601, 0.07401111],
                         [-2.31149087,-0.05160469]])

        res = kmeans(data, initk)
        warnings.simplefilter('ignore', UserWarning)
        try:
            res = kmeans2(data, initk, missing = 'warn')
        finally:
            warnings.simplefilter('default', UserWarning)

        try :
            res = kmeans2(data, initk, missing = 'raise')
            raise AssertionError("Exception not raised ! Should not happen")
        except ClusterError, e:
            pass

    def test_kmeans2_simple(self):
        """Testing simple call to kmeans2 and its results."""
        initc = np.concatenate(([[X[0]], [X[1]], [X[2]]]))
        code = initc.copy()
        code1 = kmeans2(X, code, iter = 1)[0]
        code2 = kmeans2(X, code, iter = 2)[0]

        assert_array_almost_equal(code1, CODET1)
        assert_array_almost_equal(code2, CODET2)

    def test_kmeans2_rank1(self):
        """Testing simple call to kmeans2 with rank 1 data."""
        data = np.fromfile(open(DATAFILE1), sep = ", ")
        data = data.reshape((200, 2))
        data1 = data[:, 0]
        data2 = data[:, 1]

        initc = data1[:3]
        code = initc.copy()
        code1 = kmeans2(data1, code, iter = 1)[0]
        code2 = kmeans2(data1, code, iter = 2)[0]

    def test_kmeans2_rank1_2(self):
        """Testing simple call to kmeans2 with rank 1 data."""
        data = np.fromfile(open(DATAFILE1), sep = ", ")
        data = data.reshape((200, 2))
        data1 = data[:, 0]

        code1 = kmeans2(data1, 2, iter = 1)

    def test_kmeans2_init(self):
        """Testing that kmeans2 init methods work."""
        data = np.fromfile(open(DATAFILE1), sep = ", ")
        data = data.reshape((200, 2))

        kmeans2(data, 3, minit = 'random')
        kmeans2(data, 3, minit = 'points')

        # Check special case 1d
        data = data[:, :1]
        kmeans2(data, 3, minit = 'random')
        kmeans2(data, 3, minit = 'points')

    def test_kmeans2_empty(self):
        """Ticket #505."""
        try:
            kmeans2([], 2)
            raise AssertionError("This should not succeed.")
        except ValueError, e:
            # OK, that's what we expect
            pass

    def test_kmeans_0k(self):
        """Regression test for #546: fail when k arg is 0."""
        try:
            kmeans(X, 0)
            raise AssertionError("kmeans with 0 clusters should fail.")
        except ValueError:
            pass

        try:
            kmeans2(X, 0)
            raise AssertionError("kmeans2 with 0 clusters should fail.")
        except ValueError:
            pass

        try:
            kmeans2(X, np.array([]))
            raise AssertionError("kmeans2 with 0 clusters should fail.")
        except ValueError:
            pass

if __name__ == "__main__":
    run_module_suite()
