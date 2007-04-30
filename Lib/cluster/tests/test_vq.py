#! /usr/bin/env python

# David Cournapeau
# Last Change: Thu Apr 26 09:00 PM 2007 J

# For now, just copy the tests from sandbox.pyem, so we can check that
# kmeans works OK for trivial examples.

import sys
from numpy.testing import *

import numpy as N

set_package_path()
from cluster.vq import kmeans, kmeans2, py_vq, py_vq2, _py_vq_1d
restore_path()

#Optional:
set_local_path()
# import modules that are located in the same directory as this file.
import os.path
DATAFILE1 = os.path.join(sys.path[0], "data.txt")
restore_path()

# Global data
X   = N.array([[3.0, 3], [4, 3], [4, 2],
               [9, 2], [5, 1], [6, 2], [9, 4],
               [5, 2], [5, 4], [7, 4], [6, 5]])

CODET1  = N.array([[3.0000, 3.0000],
                   [6.2000, 4.0000],
                   [5.8000, 1.8000]])

CODET2  = N.array([[11.0/3, 8.0/3],
                   [6.7500, 4.2500],
                   [6.2500, 1.7500]])

LABEL1  = N.array([0, 1, 2, 2, 2, 2, 1, 2, 1, 1, 1])

class test_vq(NumpyTestCase):
    def check_py_vq(self, level=1):
        initc = N.concatenate(([[X[0]], [X[1]], [X[2]]]))
        code = initc.copy()
        label1 = py_vq(X, initc)[0]
        assert_array_equal(label1, LABEL1)

    def check_py_vq2(self, level=1):
        initc = N.concatenate(([[X[0]], [X[1]], [X[2]]]))
        code = initc.copy()
        label1 = py_vq2(X, initc)[0]
        assert_array_equal(label1, LABEL1)

    def check_vq(self, level=1):
        initc = N.concatenate(([[X[0]], [X[1]], [X[2]]]))
        code = initc.copy()
        try:
            import _vq
            label1 = _vq.double_vq(X, initc)[0]
            assert_array_equal(label1, LABEL1)
        except ImportError:
            print "== Error while importing _vq, not testing C imp of vq =="

    #def check_vq_1d(self, level=1):
    #    data = X[:, 0]
    #    initc = data[:3]
    #    code = initc.copy()
    #    print _py_vq_1d(data, initc)

class test_kmean(NumpyTestCase):
    def check_kmeans_simple(self, level=1):
        initc = N.concatenate(([[X[0]], [X[1]], [X[2]]]))
        code = initc.copy()
        code1 = kmeans(X, code, iter = 1)[0]

        assert_array_almost_equal(code1, CODET2)

    def check_kmeans_lost_cluster(self, level=1):
        """This will cause kmean to have a cluster with no points."""
        data = N.fromfile(open(DATAFILE1), sep = ", ")
        data = data.reshape((200, 2))
        initk = N.array([[-1.8127404, -0.67128041],
                         [ 2.04621601, 0.07401111],
                         [-2.31149087,-0.05160469]])

        res = kmeans(data, initk)

    def check_kmeans2_simple(self, level=1):
        """Testing simple call to kmeans2 and its results."""
        initc = N.concatenate(([[X[0]], [X[1]], [X[2]]]))
        code = initc.copy()
        code1 = kmeans2(X, code, niter = 1)[0]
        code2 = kmeans2(X, code, niter = 2)[0]

        assert_array_almost_equal(code1, CODET1)
        assert_array_almost_equal(code2, CODET2)

    #def check_kmeans2_rank1(self, level=1):
    #    """Testing simple call to kmeans2 with rank 1 data."""
    #    data = N.fromfile(open(DATAFILE1), sep = ", ")
    #    data = data.reshape((200, 2))
    #    data1 = data[:, 0]
    #    data2 = data[:, 1]

    #    initc = data1[:3]
    #    code = initc.copy()
    #    print _py_vq_1d(data1, code)
    #    code1 = kmeans2(data1, code, niter = 1)[0]
    #    code2 = kmeans2(data1, code, niter = 2)[0]

    def check_kmeans2_init(self, level = 1):
        """Testing that kmeans2 init methods work."""
        data = N.fromfile(open(DATAFILE1), sep = ", ")
        data = data.reshape((200, 2))

        kmeans2(data, 3, minit = 'random')
        kmeans2(data, 3, minit = 'points')

if __name__ == "__main__":
    NumpyTest().run()
