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

import sys
import os.path
import numpy

from scipy.testing import *
from scipy.cluster.hierarchy import squareform, linkage, from_mlab_linkage, numobs_dm, numobs_y, numobs_linkage
from scipy.cluster.distance import pdist, matching, jaccard, dice, sokalsneath, rogerstanimoto, russellrao, yule

_tdist = numpy.array([[0,    662,  877,  255,  412,  996],
                      [662,  0,    295,  468,  268,  400],
                      [877,  295,  0,    754,  564,  138],
                      [255,  468,  754,  0,    219,  869],
                      [412,  268,  564,  219,  0,    669],
                      [996,  400,  138,  869,  669,  0  ]], dtype='double')

_ytdist = squareform(_tdist)


eo = {}

_filenames = ["iris.txt",
              "linkage-single-tdist.txt",
              "linkage-complete-tdist.txt",
              "linkage-average-tdist.txt",
              "linkage-weighted-tdist.txt",
              "random-bool-data.txt"]


def load_testing_files():
    for fn in _filenames:
        name = fn.replace(".txt", "").replace("-ml", "")
        fqfn = os.path.join(os.path.dirname(__file__), fn)
        eo[name] = numpy.loadtxt(open(fqfn))
        #print "%s: %s   %s" % (name, str(eo[name].shape), str(eo[name].dtype))
    #eo['pdist-boolean-inp'] = numpy.bool_(eo['pdist-boolean-inp'])

load_testing_files()

class TestSquareForm(TestCase):

    ################### squareform
    def test_squareform_empty_matrix(self):
        "Tests squareform on an empty matrix."
        A = numpy.zeros((0,0))
        rA = squareform(numpy.array(A, dtype='double'))
        self.failUnless(rA.shape == (0,))

    def test_squareform_empty_vector(self):
        v = numpy.zeros((0,))
        rv = squareform(numpy.array(v, dtype='double'))
        self.failUnless(rv.shape == (1,1))
        self.failUnless(rv[0, 0] == 0)

    def test_squareform_1by1_matrix(self):
        "Tests squareform on a 1x1 matrix."
        A = numpy.zeros((1,1))
        rA = squareform(numpy.array(A, dtype='double'))
        self.failUnless(rA.shape == (0,))

    def test_squareform_one_vector(self):
        "Tests squareform on a 1-D array, length=1."
        v = numpy.ones((1,)) * 8.3
        rv = squareform(numpy.array(v, dtype='double'))
        self.failUnless(rv.shape == (2,2))
        self.failUnless(rv[0,1] == 8.3)
        self.failUnless(rv[1,0] == 8.3)

    def test_squareform_2by2_matrix(self):
        "Tests squareform on a 2x2 matrix."
        A = numpy.zeros((2,2))
        A[0,1]=0.8
        A[1,0]=0.8
        rA = squareform(numpy.array(A, dtype='double'))
        self.failUnless(rA.shape == (1,))
        self.failUnless(rA[0] == 0.8)

    def test_squareform_multi_matrix(self):
        "Tests squareform on a square matrices of multiple sizes."
        for n in xrange(2, 5):
            yield self.check_squareform_multi_matrix(n)

    def check_squareform_multi_matrix(self, n):
        X = numpy.random.rand(n, 4)
        Y = pdist(X)
        self.failUnless(len(Y.shape) == 1)
        A = squareform(Y)
        Yr = squareform(A)
        s = A.shape
        k = 0
        print A.shape, Y.shape, Yr.shape
        self.failUnless(len(s) == 2)
        self.failUnless(len(Yr.shape) == 1)
        self.failUnless(s[0] == s[1])
        for i in xrange(0, s[0]):
            for j in xrange(i+1, s[1]):
                if i != j:
                    #print i, j, k, A[i, j], Y[k]
                    self.failUnless(A[i, j] == Y[k])
                    k += 1
                else:
                    self.failUnless(A[i, j] == 0)

class TestNumObs(TestCase):

    ############## numobs_dm
    def test_numobs_dm_multi_matrix(self):
        "Tests numobs_dm with observation matrices of multiple sizes."
        for n in xrange(1, 10):
            X = numpy.random.rand(n, 4)
            Y = pdist(X)
            A = squareform(Y)
            print A.shape, Y.shape
            self.failUnless(numobs_dm(A) == n)

    def test_numobs_y_multi_matrix(self):
        "Tests numobs_y with observation matrices of multiple sizes."
        for n in xrange(2, 10):
            X = numpy.random.rand(n, 4)
            Y = pdist(X)
            #print A.shape, Y.shape, Yr.shape
            self.failUnless(numobs_y(Y) == n)

    def test_numobs_linkage_multi_matrix(self):
        "Tests numobs_linkage with observation matrices of multiple sizes."
        for n in xrange(2, 10):
            X = numpy.random.rand(n, 4)
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
        #print Z, expectedZ, numpy.abs(Z - expectedZ).max()
        self.failUnless(within_tol(Z, expectedZ, eps))

    def test_linkage_weighted_tdist(self):
        "Tests linkage(Y, 'weighted') on the tdist data set."
        Z = linkage(_ytdist, 'weighted')
        Zmlab = eo['linkage-weighted-tdist']
        eps = 1e-10
        expectedZ = from_mlab_linkage(Zmlab)
        #print Z, expectedZ, numpy.abs(Z - expectedZ).max()
        self.failUnless(within_tol(Z, expectedZ, eps))

def within_tol(a, b, tol):
    return numpy.abs(a - b).max() < tol

if __name__ == "__main__":
    nose.run(argv=['', __file__])

