#! /usr/bin/env python

# David Cournapeau
# Last Change: Thu Apr 26 05:00 PM 2007 J

# For now, just copy the tests from sandbox.pyem, so we can check that
# kmeans works OK for trivial examples.

import sys
from numpy.testing import *

import numpy as N

set_package_path()
from cluster.vq import kmeans, kmeans_, py_vq, py_vq2
restore_path()

#Optional:
set_local_path()
# import modules that are located in the same directory as this file.
import os.path
DATAFILE1   = os.path.join(sys.path[0], "data.txt")
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
        initc   = N.concatenate(([[X[0]], [X[1]], [X[2]]])) 
        code    = initc.copy()
        label1  = py_vq(X, initc)[0]
        assert_array_equal(label1, LABEL1)

    def check_py_vq2(self, level=1):
        initc   = N.concatenate(([[X[0]], [X[1]], [X[2]]])) 
        code    = initc.copy()
        label1  = py_vq2(X, initc)[0]
        assert_array_equal(label1, LABEL1)

    def check_vq(self, level=1):
        initc   = N.concatenate(([[X[0]], [X[1]], [X[2]]])) 
        code    = initc.copy()
        try:
            import _vq
            label1  = _vq.double_vq(X, initc)[0]
            assert_array_equal(label1, LABEL1)
        except ImportError:
            print "== Error while importing _vq, not testing C imp of vq =="

class test_kmean(NumpyTestCase):
    def check_kmeans(self, level=1):
        initc   = N.concatenate(([[X[0]], [X[1]], [X[2]]])) 
        code    = initc.copy()
        #code1   = kmeans(X, code, iter = 1)[0]

        #assert_array_almost_equal(code1, CODET2)

    def check_kmeans_lost_cluster(self, level=1):
        """This will cause kmean to have a cluster with no points."""
        data    = N.fromfile(open(DATAFILE1), sep = ", ")
        data    = data.reshape((200, 2))
        initk   = N.array([[-1.8127404, -0.67128041], [ 2.04621601, 0.07401111], 
                    [-2.31149087,-0.05160469]])

        res     = kmeans(data, initk)

if __name__ == "__main__":
    NumpyTest().run()
