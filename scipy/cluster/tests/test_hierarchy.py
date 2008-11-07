#! /usr/bin/env python
#
# Author: Damian Eads
# Date: April 17, 2008
#
# Copyright (C) 2008 Damian Eads
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#
# 3. The name of the author may not be used to endorse or promote
#    products derived from this software without specific prior
#    written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os.path

import numpy as np
from numpy.testing import *

from scipy.cluster.hierarchy import linkage, from_mlab_linkage, to_mlab_linkage, numobs_linkage, inconsistent, cophenet, from_mlab_linkage, fclusterdata, fcluster, is_isomorphic, single, complete, average, weighted, centroid, median, ward
from scipy.spatial.distance import squareform, pdist, numobs_dm, numobs_y

_tdist = np.array([[0,    662,  877,  255,  412,  996],
                   [662,  0,    295,  468,  268,  400],
                   [877,  295,  0,    754,  564,  138],
                   [255,  468,  754,  0,    219,  869],
                   [412,  268,  564,  219,  0,    669],
                   [996,  400,  138,  869,  669,  0  ]], dtype='double')

_ytdist = squareform(_tdist)


eo = {}

_filenames = ["iris.txt",
              "Q-X.txt",
              "fclusterdata-maxclusts-2.txt",
              "fclusterdata-maxclusts-3.txt",
              "fclusterdata-maxclusts-4.txt",
              "linkage-single-tdist.txt",
              "linkage-complete-tdist.txt",
              "linkage-average-tdist.txt",
              "linkage-weighted-tdist.txt",
              "inconsistent-Q-single-1.txt",
              "inconsistent-Q-single-2.txt",
              "inconsistent-Q-single-3.txt",
              "inconsistent-Q-single-4.txt",
              "inconsistent-Q-single-5.txt",
              "inconsistent-Q-single-6.txt",
              "inconsistent-complete-tdist-depth-1.txt",
              "inconsistent-complete-tdist-depth-2.txt",
              "inconsistent-complete-tdist-depth-3.txt",
              "inconsistent-complete-tdist-depth-4.txt",
              "inconsistent-single-tdist-depth-0.txt",
              "inconsistent-single-tdist-depth-1.txt",
              "inconsistent-single-tdist-depth-2.txt",
              "inconsistent-single-tdist-depth-3.txt",
              "inconsistent-single-tdist-depth-4.txt",
              "inconsistent-single-tdist-depth-5.txt",
              "inconsistent-single-tdist.txt",
              "inconsistent-weighted-tdist-depth-1.txt",
              "inconsistent-weighted-tdist-depth-2.txt",
              "inconsistent-weighted-tdist-depth-3.txt",
              "inconsistent-weighted-tdist-depth-4.txt",
              "linkage-Q-average.txt",
              "linkage-Q-complete.txt",
              "linkage-Q-single.txt",
              "linkage-Q-weighted.txt",
              "linkage-Q-centroid.txt",
              "linkage-Q-median.txt",
              "linkage-Q-ward.txt"
              ]

def load_testing_files():
    for fn in _filenames:
        name = fn.replace(".txt", "").replace("-ml", "")
        fqfn = os.path.join(os.path.dirname(__file__), fn)
        eo[name] = np.loadtxt(open(fqfn))
        #print "%s: %s   %s" % (name, str(eo[name].shape), str(eo[name].dtype))
    #eo['pdist-boolean-inp'] = np.bool_(eo['pdist-boolean-inp'])

load_testing_files()

class TestNumObs(TestCase):

    ############## numobs_dm
    def test_numobs_dm_multi_matrix(self):
        "Tests numobs_dm with observation matrices of multiple sizes."
        for n in xrange(1, 10):
            X = np.random.rand(n, 4)
            Y = pdist(X)
            A = squareform(Y)
            if verbose >= 3:
                print A.shape, Y.shape
            self.failUnless(numobs_dm(A) == n)

    def test_numobs_y_multi_matrix(self):
        "Tests numobs_y with observation matrices of multiple sizes."
        for n in xrange(2, 10):
            X = np.random.rand(n, 4)
            Y = pdist(X)
            #print A.shape, Y.shape, Yr.shape
            self.failUnless(numobs_y(Y) == n)

    def test_numobs_linkage_multi_matrix(self):
        "Tests numobs_linkage with observation matrices of multiple sizes."
        for n in xrange(2, 10):
            X = np.random.rand(n, 4)
            Y = pdist(X)
            Z = linkage(Y)
            #print Z
            #print A.shape, Y.shape, Yr.shape
            self.failUnless(numobs_linkage(Z) == n)

class TestLinkage(TestCase):

    ################### linkage
    def test_linkage_single_tdist(self):
        "Tests linkage(Y, 'single') on the tdist data set."
        Z = linkage(_ytdist, 'single')
        Zmlab = eo['linkage-single-tdist']
        eps = 1e-10
        expectedZ = from_mlab_linkage(Zmlab)
        self.failUnless(within_tol(Z, expectedZ, eps))

    def test_linkage_complete_tdist(self):
        "Tests linkage(Y, 'complete') on the tdist data set."
        Z = linkage(_ytdist, 'complete')
        Zmlab = eo['linkage-complete-tdist']
        eps = 1e-10
        expectedZ = from_mlab_linkage(Zmlab)
        self.failUnless(within_tol(Z, expectedZ, eps))

    def test_linkage_average_tdist(self):
        "Tests linkage(Y, 'average') on the tdist data set."
        Z = linkage(_ytdist, 'average')
        Zmlab = eo['linkage-average-tdist']
        eps = 1e-05
        expectedZ = from_mlab_linkage(Zmlab)
        #print Z, expectedZ, np.abs(Z - expectedZ).max()
        self.failUnless(within_tol(Z, expectedZ, eps))

    def test_linkage_weighted_tdist(self):
        "Tests linkage(Y, 'weighted') on the tdist data set."
        Z = linkage(_ytdist, 'weighted')
        Zmlab = eo['linkage-weighted-tdist']
        eps = 1e-10
        expectedZ = from_mlab_linkage(Zmlab)
        #print Z, expectedZ, np.abs(Z - expectedZ).max()
        self.failUnless(within_tol(Z, expectedZ, eps))

    ################### linkage on Q
    def test_linkage_single_q(self):
        "Tests linkage(Y, 'single') on the Q data set."
        X = eo['Q-X']
        Z = single(X)
        Zmlab = eo['linkage-Q-single']
        eps = 1e-06
        expectedZ = from_mlab_linkage(Zmlab)
        print abs(Z-expectedZ).max()
        self.failUnless(within_tol(Z, expectedZ, eps))

    def test_linkage_complete_q(self):
        "Tests linkage(Y, 'complete') on the Q data set."
        X = eo['Q-X']
        Z = complete(X)
        Zmlab = eo['linkage-Q-complete']
        eps = 1e-07
        expectedZ = from_mlab_linkage(Zmlab)
        print abs(Z-expectedZ).max()
        self.failUnless(within_tol(Z, expectedZ, eps))

    def test_linkage_centroid_q(self):
        "Tests linkage(Y, 'centroid') on the Q data set."
        X = eo['Q-X']
        Z = centroid(X)
        Zmlab = eo['linkage-Q-centroid']
        eps = 1e-07
        expectedZ = from_mlab_linkage(Zmlab)
        print abs(Z-expectedZ).max()
        self.failUnless(within_tol(Z, expectedZ, eps))

    def test_linkage_weighted_q(self):
        "Tests linkage(Y, 'weighted') on the Q data set."
        X = eo['Q-X']
        Z = weighted(X)
        Zmlab = eo['linkage-Q-weighted']
        eps = 1e-07
        expectedZ = from_mlab_linkage(Zmlab)
        print abs(Z-expectedZ).max()
        self.failUnless(within_tol(Z, expectedZ, eps))

class TestInconsistent(TestCase):

    def test_single_inconsistent_tdist_1(self):
        "Testing inconsistency matrix calculation (depth=1) on a single linkage."
        Y = squareform(_tdist)
        Z = linkage(Y, 'single')
        R = inconsistent(Z, 1)
        Rright = eo['inconsistent-single-tdist-depth-1']
        eps = 1e-15
        print np.abs(R - Rright).max()
        self.failUnless(within_tol(R, Rright, eps))

    def test_single_inconsistent_tdist_2(self):
        "Testing inconsistency matrix calculation (depth=2) on a single linkage."
        Y = squareform(_tdist)
        Z = linkage(Y, 'single')
        R = inconsistent(Z, 2)
        Rright = eo['inconsistent-single-tdist-depth-2']
        eps = 1e-05
        print np.abs(R - Rright).max()
        self.failUnless(within_tol(R, Rright, eps))

    def test_single_inconsistent_tdist_3(self):
        "Testing inconsistency matrix calculation (depth=3) on a single linkage."
        Y = squareform(_tdist)
        Z = linkage(Y, 'single')
        R = inconsistent(Z, 3)
        Rright = eo['inconsistent-single-tdist-depth-3']
        eps = 1e-05
        print np.abs(R - Rright).max()
        self.failUnless(within_tol(R, Rright, eps))

    def test_single_inconsistent_tdist_4(self):
        "Testing inconsistency matrix calculation (depth=4) on a single linkage."
        Y = squareform(_tdist)
        Z = linkage(Y, 'single')
        R = inconsistent(Z, 4)
        Rright = eo['inconsistent-single-tdist-depth-4']
        eps = 1e-05
        print np.abs(R - Rright).max()
        self.failUnless(within_tol(R, Rright, eps))

    # with complete linkage...

    def test_complete_inconsistent_tdist_1(self):
        "Testing inconsistency matrix calculation (depth=1) on a complete linkage."
        Y = squareform(_tdist)
        Z = linkage(Y, 'complete')
        R = inconsistent(Z, 1)
        Rright = eo['inconsistent-complete-tdist-depth-1']
        eps = 1e-15
        print np.abs(R - Rright).max()
        self.failUnless(within_tol(R, Rright, eps))

    def test_complete_inconsistent_tdist_2(self):
        "Testing inconsistency matrix calculation (depth=2) on a complete linkage."
        Y = squareform(_tdist)
        Z = linkage(Y, 'complete')
        R = inconsistent(Z, 2)
        Rright = eo['inconsistent-complete-tdist-depth-2']
        eps = 1e-05
        print np.abs(R - Rright).max()
        self.failUnless(within_tol(R, Rright, eps))

    def test_complete_inconsistent_tdist_3(self):
        "Testing inconsistency matrix calculation (depth=3) on a complete linkage."
        Y = squareform(_tdist)
        Z = linkage(Y, 'complete')
        R = inconsistent(Z, 3)
        Rright = eo['inconsistent-complete-tdist-depth-3']
        eps = 1e-05
        print np.abs(R - Rright).max()
        self.failUnless(within_tol(R, Rright, eps))

    def test_complete_inconsistent_tdist_4(self):
        "Testing inconsistency matrix calculation (depth=4) on a complete linkage."
        Y = squareform(_tdist)
        Z = linkage(Y, 'complete')
        R = inconsistent(Z, 4)
        Rright = eo['inconsistent-complete-tdist-depth-4']
        eps = 1e-05
        print np.abs(R - Rright).max()
        self.failUnless(within_tol(R, Rright, eps))

    # with single linkage and Q data set

    def test_single_inconsistent_tdist_1(self):
        "Testing inconsistency matrix calculation (depth=1) on a weighted linkage."
        X = eo['Q-X']
        Z = linkage(X, 'single', 'euclidean')
        R = inconsistent(Z, 1)
        Rright = eo['inconsistent-Q-single-1']
        eps = 1e-06
        print np.abs(R - Rright).max()
        self.failUnless(within_tol(R, Rright, eps))

    def test_single_inconsistent_tdist_2(self):
        "Testing inconsistency matrix calculation (depth=2) on a weighted linkage."
        X = eo['Q-X']
        Z = linkage(X, 'single', 'euclidean')
        R = inconsistent(Z, 2)
        Rright = eo['inconsistent-Q-single-2']
        eps = 1e-06
        print np.abs(R - Rright).max()
        self.failUnless(within_tol(R, Rright, eps))

    def test_single_inconsistent_tdist_3(self):
        "Testing inconsistency matrix calculation (depth=3) on a weighted linkage."
        X = eo['Q-X']
        Z = linkage(X, 'single', 'euclidean')
        R = inconsistent(Z, 3)
        Rright = eo['inconsistent-Q-single-3']
        eps = 1e-05
        print np.abs(R - Rright).max()
        self.failUnless(within_tol(R, Rright, eps))

    def test_single_inconsistent_tdist_4(self):
        "Testing inconsistency matrix calculation (depth=4) on a weighted linkage."
        X = eo['Q-X']
        Z = linkage(X, 'single', 'euclidean')
        R = inconsistent(Z, 4)
        Rright = eo['inconsistent-Q-single-4']
        eps = 1e-05
        print np.abs(R - Rright).max()
        self.failUnless(within_tol(R, Rright, eps))

class TestCopheneticDistance(TestCase):

    def test_linkage_cophenet_tdist_Z(self):
        "Testing cophenet(Z) on tdist data set."
        expectedM = np.array([268, 295, 255, 255, 295, 295, 268, 268, 295, 295, 295, 138, 219, 295, 295]);
        Z = linkage(_ytdist, 'single')
        M = cophenet(Z)
        eps = 1e-10
        self.failUnless(within_tol(M, expectedM, eps))

    def test_linkage_cophenet_tdist_Z_Y(self):
        "Testing cophenet(Z, Y) on tdist data set."
        Z = linkage(_ytdist, 'single')
        c = cophenet(Z, _ytdist)
        expectedc = 0.639931296433393415057366837573
        eps = 1e-10
        self.failUnless(np.abs(c - expectedc) <= eps)

    def test_linkage_cophenet_tdist_Z_Y_EL(self):
        "Testing cophenet(Z, Y, []) on tdist data set."
        Z = linkage(_ytdist, 'single')
        (c, M) = cophenet(Z, _ytdist, [])
        eps = 1e-10
        expectedM = np.array([268, 295, 255, 255, 295, 295, 268, 268, 295, 295, 295, 138, 219, 295, 295]);
        expectedc = 0.639931296433393415057366837573
        self.failUnless(np.abs(c - expectedc) <= eps)
        self.failUnless(within_tol(M, expectedM, eps))

class TestFromMLabLinkage(TestCase):

    def test_from_mlab_linkage_empty(self):
        "Testing from_mlab_linkage on empty linkage array."
        X = np.asarray([])
        R = from_mlab_linkage([])
        self.failUnless((R == X).all())

    def test_from_mlab_linkage_single_row(self):
        "Testing from_mlab_linkage on linkage array with single row."
        expectedZP = np.asarray([[ 0.,  1.,  3.,  2.]])
        Z = [[1,2,3]]
        ZP = from_mlab_linkage(Z)
        return self.failUnless((ZP == expectedZP).all())

    def test_from_mlab_linkage_multiple_rows(self):
        "Testing from_mlab_linkage on linkage array with multiple rows."
        Z = np.asarray([[3, 6, 138], [4, 5, 219],
                        [1, 8, 255], [2, 9, 268], [7, 10, 295]])
        expectedZS = np.array([[   2.,    5.,  138.,    2.],
                               [   3.,    4.,  219.,    2.],
                               [   0.,    7.,  255.,    3.],
                               [   1.,    8.,  268.,    4.],
                               [   6.,    9.,  295.,    6.]],
                              dtype=np.double)
        ZS = from_mlab_linkage(Z)
        print expectedZS, ZS
        self.failUnless((expectedZS == ZS).all())


class TestToMLabLinkage(TestCase):

    def test_to_mlab_linkage_empty(self):
        "Testing to_mlab_linkage on empty linkage array."
        X = np.asarray([])
        R = to_mlab_linkage([])
        self.failUnless((R == X).all())

    def test_to_mlab_linkage_single_row(self):
        "Testing to_mlab_linkage on linkage array with single row."
        Z = np.asarray([[ 0.,  1.,  3.,  2.]])
        expectedZP = np.asarray([[1,2,3]])
        ZP = to_mlab_linkage(Z)
        return self.failUnless((ZP == expectedZP).all())

    def test_from_mlab_linkage_multiple_rows(self):
        "Testing to_mlab_linkage on linkage array with multiple rows."
        expectedZM = np.asarray([[3, 6, 138], [4, 5, 219],
                        [1, 8, 255], [2, 9, 268], [7, 10, 295]])
        Z = np.array([[   2.,    5.,  138.,    2.],
                      [   3.,    4.,  219.,    2.],
                      [   0.,    7.,  255.,    3.],
                      [   1.,    8.,  268.,    4.],
                      [   6.,    9.,  295.,    6.]],
                     dtype=np.double)
        ZM = to_mlab_linkage(Z)
        print expectedZM, ZM
        self.failUnless((expectedZM == ZM).all())

class TestFcluster(TestCase):

    def test_fclusterdata_maxclusts_2(self):
        "Tests fclusterdata(X, criterion='maxclust', t=2) on a random 3-cluster data set."
        expectedT = eo['fclusterdata-maxclusts-2']
        X = eo['Q-X']
        T = fclusterdata(X, criterion='maxclust', t=2)
        self.failUnless(is_isomorphic(T, expectedT))

    def test_fclusterdata_maxclusts_3(self):
        "Tests fclusterdata(X, criterion='maxclust', t=3) on a random 3-cluster data set."
        expectedT = eo['fclusterdata-maxclusts-3']
        X = eo['Q-X']
        T = fclusterdata(X, criterion='maxclust', t=3)
        self.failUnless(is_isomorphic(T, expectedT))

    def test_fclusterdata_maxclusts_4(self):
        "Tests fclusterdata(X, criterion='maxclust', t=4) on a random 3-cluster data set."
        expectedT = eo['fclusterdata-maxclusts-4']
        X = eo['Q-X']
        T = fclusterdata(X, criterion='maxclust', t=4)
        self.failUnless(is_isomorphic(T, expectedT))

    def test_fcluster_maxclusts_2(self):
        "Tests fcluster(Z, criterion='maxclust', t=2) on a random 3-cluster data set."
        expectedT = eo['fclusterdata-maxclusts-2']
        X = eo['Q-X']
        Y = pdist(X)
        Z = linkage(Y)
        T = fcluster(Z, criterion='maxclust', t=2)
        self.failUnless(is_isomorphic(T, expectedT))

    def test_fcluster_maxclusts_3(self):
        "Tests fcluster(Z, criterion='maxclust', t=3) on a random 3-cluster data set."
        expectedT = eo['fclusterdata-maxclusts-3']
        X = eo['Q-X']
        Y = pdist(X)
        Z = linkage(Y)
        T = fcluster(Z, criterion='maxclust', t=3)
        self.failUnless(is_isomorphic(T, expectedT))

    def test_fcluster_maxclusts_4(self):
        "Tests fcluster(Z, criterion='maxclust', t=4) on a random 3-cluster data set."
        expectedT = eo['fclusterdata-maxclusts-4']
        X = eo['Q-X']
        Y = pdist(X)
        Z = linkage(Y)
        T = fcluster(Z, criterion='maxclust', t=4)
        self.failUnless(is_isomorphic(T, expectedT))

def help_single_inconsistent_depth(self, i):
    Y = squareform(_tdist)
    Z = linkage(Y, 'single')
    R = inconsistent(Z, i)
    Rright = eo['inconsistent-single-tdist-depth-' + str(i)]
    eps = 1e-05
    print np.abs(R - Rright).max()
    self.failUnless(within_tol(R, Rright, eps))

def within_tol(a, b, tol):
    return np.abs(a - b).max() < tol

if __name__ == "__main__":
    run_module_suite()
