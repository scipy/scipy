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
from scipy.testing import *
from scipy.cluster.hierarchy import pdist, squareform, linkage, from_mlab_linkage, numobs_dm, numobs_y, numobs_linkage, matching, jaccard, dice, sokalsneath, rogerstanimoto, russellrao, yule

import numpy
#import math

#from scipy.cluster.hierarchy import pdist, euclidean

_filenames = ["iris.txt",
              "pdist-hamming-ml.txt",
              "pdist-boolean-inp.txt",
              "pdist-jaccard-ml.txt",
              "pdist-cityblock-ml-iris.txt",
              "pdist-minkowski-3.2-ml-iris.txt",
              "pdist-cityblock-ml.txt",
              "pdist-correlation-ml-iris.txt",
              "pdist-minkowski-5.8-ml-iris.txt",
              "pdist-correlation-ml.txt",
              "pdist-minkowski-3.2-ml.txt",
              "pdist-cosine-ml-iris.txt",
              "pdist-seuclidean-ml-iris.txt",
              "pdist-cosine-ml.txt",
              "pdist-seuclidean-ml.txt",
              "pdist-double-inp.txt",
              "pdist-spearman-ml.txt",
              "pdist-euclidean-ml.txt",
              "pdist-euclidean-ml-iris.txt",
              "pdist-chebychev-ml.txt",
              "pdist-chebychev-ml-iris.txt",
              "linkage-single-tdist.txt",
              "linkage-complete-tdist.txt",
              "linkage-average-tdist.txt",
              "linkage-weighted-tdist.txt"]

_tdist = numpy.array([[0,    662,  877,  255,  412,  996],
                      [662,  0,    295,  468,  268,  400],
                      [877,  295,  0,    754,  564,  138],
                      [255,  468,  754,  0,    219,  869],
                      [412,  268,  564,  219,  0,    669],
                      [996,  400,  138,  869,  669,  0  ]], dtype='double')

_ytdist = squareform(_tdist)

# A hashmap of expected output arrays for the tests. These arrays
# come from a list of text files, which are read prior to testing.

eo = {}

def load_testing_files():
    for fn in _filenames:
        name = fn.replace(".txt", "").replace("-ml", "")
        fqfn = os.path.join(os.path.dirname(__file__), fn)
        eo[name] = numpy.loadtxt(open(fqfn))
        #print "%s: %s   %s" % (name, str(eo[name].shape), str(eo[name].dtype))
    eo['pdist-boolean-inp'] = numpy.bool_(eo['pdist-boolean-inp'])

load_testing_files()

#print eo.keys()


#print numpy.abs(Y_test2 - Y_right).max()
#print numpy.abs(Y_test1 - Y_right).max()

class TestPdist(TestCase):

    def test_pdist_raises_type_error_float32(self):
        "Testing whether passing a float32 observation array generates an exception."
        X = numpy.zeros((10, 10), dtype=numpy.float32)
        try:
            pdist(X, 'euclidean')
        except TypeError:
            pass
        except:
            self.fail("float32 observation matrices should generate an error in pdist.")

    def test_pdist_raises_type_error_longdouble(self):
        "Testing whether passing a longdouble observation array generates an exception."
        X = numpy.zeros((10, 10), dtype=numpy.longdouble)
        try:
            pdist(X, 'euclidean')
        except TypeError:
            pass
        except:
            self.fail("longdouble observation matrices should generate an error in pdist.")

    def test_pdist_var_raises_type_error_float32(self):
        "Testing whether passing a float32 variance matrix generates an exception."
        X = numpy.zeros((10, 10))
        V = numpy.zeros((10, 10), dtype=numpy.float32)
        try:
            pdist(X, 'seuclidean', V=V)
        except TypeError:
            pass
        except:
            self.fail("float32 V matrices should generate an error in pdist('seuclidean').")

    def test_pdist_var_raises_type_error_longdouble(self):
        "Testing whether passing a longdouble variance matrix generates an exception."
        X = numpy.zeros((10, 10))
        V = numpy.zeros((10, 10), dtype=numpy.longdouble)

        try:
            pdist(X, 'seuclidean', V=V)
        except TypeError:
            pass
        except:
            self.fail("longdouble matrices should generate an error in pdist('seuclidean').")

    def test_pdist_ivar_raises_type_error_float32(self):
        "Testing whether passing a float32 variance matrix generates an exception."
        X = numpy.zeros((10, 10))
        VI = numpy.zeros((10, 10), dtype=numpy.float32)
        try:
            pdist(X, 'mahalanobis', VI=VI)
        except TypeError:
            pass
        except:
            self.fail("float32 matrices should generate an error in pdist('mahalanobis').")

    def test_pdist_ivar_raises_type_error_longdouble(self):
        "Testing whether passing a longdouble variance matrix generates an exception."
        X = numpy.zeros((10, 10))
        VI = numpy.zeros((10, 10), dtype=numpy.longdouble)

        try:
            pdist(X, 'mahalanobis', VI=VI)
        except TypeError:
            pass
        except:
            self.fail("longdouble matrices should generate an error in pdist('mahalanobis').")

    ################### pdist: euclidean
    def test_pdist_euclidean_random(self):
        "Tests pdist(X, 'euclidean') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-euclidean']

        Y_test1 = pdist(X, 'euclidean')
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_euclidean_random_nonC(self):
        "Tests pdist(X, 'test_euclidean') [the non-C implementation] on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-euclidean']
        Y_test2 = pdist(X, 'test_euclidean')
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    def test_pdist_euclidean_iris(self):
        "Tests pdist(X, 'euclidean') on the Iris data set."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-euclidean-iris']

        Y_test1 = pdist(X, 'euclidean')
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_euclidean_iris_nonC(self):
        "Tests pdist(X, 'test_euclidean') [the non-C implementation] on the Iris data set."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-euclidean-iris']
        Y_test2 = pdist(X, 'test_euclidean')
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: seuclidean
    def test_pdist_seuclidean_random(self):
        "Tests pdist(X, 'seuclidean') on random data."
        eps = 1e-05
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-seuclidean']

        Y_test1 = pdist(X, 'seuclidean')
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_seuclidean_random_nonC(self):
        "Tests pdist(X, 'test_sqeuclidean') [the non-C implementation] on random data."
        eps = 1e-05
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-seuclidean']
        Y_test2 = pdist(X, 'test_sqeuclidean')
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    def test_pdist_seuclidean_iris(self):
        "Tests pdist(X, 'seuclidean') on the Iris data set."
        eps = 1e-05
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-seuclidean-iris']

        Y_test1 = pdist(X, 'seuclidean')
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_seuclidean_iris_nonC(self):
        "Tests pdist(X, 'test_seuclidean') [the non-C implementation] on the Iris data set."
        eps = 1e-05
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-seuclidean-iris']
        Y_test2 = pdist(X, 'test_sqeuclidean')
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: cosine
    def test_pdist_cosine_random(self):
        "Tests pdist(X, 'cosine') on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-cosine']

        Y_test1 = pdist(X, 'cosine')
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_cosine_random_nonC(self):
        "Tests pdist(X, 'test_cosine') [the non-C implementation] on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-cosine']
        Y_test2 = pdist(X, 'test_cosine')
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    def test_pdist_cosine_iris(self):
        "Tests pdist(X, 'cosine') on the Iris data set."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-cosine-iris']

        Y_test1 = pdist(X, 'cosine')
        self.failUnless(within_tol(Y_test1, Y_right, eps))
        #print "cosine-iris", numpy.abs(Y_test1 - Y_right).max()

    def test_pdist_cosine_iris_nonC(self):
        "Tests pdist(X, 'test_cosine') [the non-C implementation] on the Iris data set."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-cosine-iris']
        Y_test2 = pdist(X, 'test_cosine')
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: cityblock
    def test_pdist_cityblock_random(self):
        "Tests pdist(X, 'cityblock') on random data."
        eps = 1e-06
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-cityblock']

        Y_test1 = pdist(X, 'cityblock')
        #print "cityblock", numpy.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_cityblock_random_nonC(self):
        "Tests pdist(X, 'test_cityblock') [the non-C implementation] on random data."
        eps = 1e-06
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-cityblock']
        Y_test2 = pdist(X, 'test_cityblock')
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    def test_pdist_cityblock_iris(self):
        "Tests pdist(X, 'cityblock') on the Iris data set."
        eps = 1e-14
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-cityblock-iris']

        Y_test1 = pdist(X, 'cityblock')
        self.failUnless(within_tol(Y_test1, Y_right, eps))
        #print "cityblock-iris", numpy.abs(Y_test1 - Y_right).max()

    def test_pdist_cityblock_iris_nonC(self):
        "Tests pdist(X, 'test_cityblock') [the non-C implementation] on the Iris data set."
        eps = 1e-14
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-cityblock-iris']
        Y_test2 = pdist(X, 'test_cityblock')
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: correlation
    def test_pdist_correlation_random(self):
        "Tests pdist(X, 'correlation') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-correlation']

        Y_test1 = pdist(X, 'correlation')
        #print "correlation", numpy.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_correlation_random_nonC(self):
        "Tests pdist(X, 'test_correlation') [the non-C implementation] on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-correlation']
        Y_test2 = pdist(X, 'test_correlation')
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    def test_pdist_correlation_iris(self):
        "Tests pdist(X, 'correlation') on the Iris data set."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-correlation-iris']

        Y_test1 = pdist(X, 'correlation')
        #print "correlation-iris", numpy.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_correlation_iris_nonC(self):
        "Tests pdist(X, 'test_correlation') [the non-C implementation] on the Iris data set."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-correlation-iris']
        Y_test2 = pdist(X, 'test_correlation')
        #print "test-correlation-iris", numpy.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################# minkowski

    def test_pdist_minkowski_random(self):
        "Tests pdist(X, 'minkowski') on random data."
        eps = 1e-05
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-minkowski-3.2']

        Y_test1 = pdist(X, 'minkowski', 3.2)
        #print "minkowski", numpy.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_minkowski_random_nonC(self):
        "Tests pdist(X, 'test_minkowski') [the non-C implementation] on random data."
        eps = 1e-05
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-minkowski-3.2']
        Y_test2 = pdist(X, 'test_minkowski', 3.2)
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    def test_pdist_minkowski_iris(self):
        "Tests pdist(X, 'minkowski') on iris data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-minkowski-3.2-iris']

        Y_test1 = pdist(X, 'minkowski', 3.2)
        #print "minkowski-iris-3.2", numpy.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_minkowski_iris_nonC(self):
        "Tests pdist(X, 'test_minkowski') [the non-C implementation] on iris data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-minkowski-3.2-iris']
        Y_test2 = pdist(X, 'test_minkowski', 3.2)
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    def test_pdist_minkowski_iris(self):
        "Tests pdist(X, 'minkowski') on iris data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-minkowski-5.8-iris']

        Y_test1 = pdist(X, 'minkowski', 5.8)
        #print "minkowski-iris-5.8", numpy.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_minkowski_iris_nonC(self):
        "Tests pdist(X, 'test_minkowski') [the non-C implementation] on iris data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-minkowski-5.8-iris']
        Y_test2 = pdist(X, 'test_minkowski', 5.8)
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: hamming
    def test_pdist_hamming_random(self):
        "Tests pdist(X, 'hamming') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-hamming']

        Y_test1 = pdist(X, 'hamming')
        #print "hamming", numpy.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_hamming_random_nonC(self):
        "Tests pdist(X, 'test_hamming') [the non-C implementation] on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-hamming']
        Y_test2 = pdist(X, 'test_hamming')
        #print "test-hamming", numpy.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: hamming (double)
    def test_pdist_dhamming_random(self):
        "Tests pdist(X, 'hamming') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = numpy.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']

        Y_test1 = pdist(X, 'hamming')
        #print "hamming", numpy.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_dhamming_random_nonC(self):
        "Tests pdist(X, 'test_hamming') [the non-C implementation] on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = numpy.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']
        Y_test2 = pdist(X, 'test_hamming')
        #print "test-hamming", numpy.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: jaccard
    def test_pdist_jaccard_random(self):
        "Tests pdist(X, 'jaccard') on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-jaccard']

        Y_test1 = pdist(X, 'jaccard')
        #print "jaccard", numpy.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_jaccard_random_nonC(self):
        "Tests pdist(X, 'test_jaccard') [the non-C implementation] on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-jaccard']
        Y_test2 = pdist(X, 'test_jaccard')
        #print "test-jaccard", numpy.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: jaccard (double)
    def test_pdist_djaccard_random(self):
        "Tests pdist(X, 'jaccard') on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = numpy.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']

        Y_test1 = pdist(X, 'jaccard')
        #print "jaccard", numpy.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_djaccard_random_nonC(self):
        "Tests pdist(X, 'test_jaccard') [the non-C implementation] on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = numpy.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']
        Y_test2 = pdist(X, 'test_jaccard')
        #print "test-jaccard", numpy.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: chebychev
    def test_pdist_chebychev_random(self):
        "Tests pdist(X, 'chebychev') on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-chebychev']

        Y_test1 = pdist(X, 'chebychev')
        #print "chebychev", numpy.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_chebychev_random_nonC(self):
        "Tests pdist(X, 'test_chebychev') [the non-C implementation] on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-chebychev']
        Y_test2 = pdist(X, 'test_chebychev')
        #print "test-chebychev", numpy.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    def test_pdist_chebychev_iris(self):
        "Tests pdist(X, 'chebychev') on the Iris data set."
        eps = 1e-15
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-chebychev-iris']

        Y_test1 = pdist(X, 'chebychev')
        #print "chebychev-iris", numpy.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_chebychev_iris_nonC(self):
        "Tests pdist(X, 'test_chebychev') [the non-C implementation] on the Iris data set."
        eps = 1e-15
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-chebychev-iris']
        Y_test2 = pdist(X, 'test_chebychev')
        #print "test-chebychev-iris", numpy.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    def test_pdist_chebychev_iris_nonC(self):
        "Tests pdist(X, 'test_chebychev') [the non-C implementation] on the Iris data set."
        eps = 1e-15
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-chebychev-iris']
        Y_test2 = pdist(X, 'test_chebychev')
        #print "test-chebychev-iris", numpy.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    def test_pdist_matching_mtica1(self):
        "Tests matching(*,*) with mtica example #1."
        m = matching(numpy.array([1, 0, 1, 1, 0]),
                     numpy.array([1, 1, 0, 1, 1]))
        m2 = matching(numpy.array([1, 0, 1, 1, 0], dtype=numpy.bool),
                      numpy.array([1, 1, 0, 1, 1], dtype=numpy.bool))
        self.failUnless(numpy.abs(m - 0.6) <= 1e-10)
        self.failUnless(numpy.abs(m2 - 0.6) <= 1e-10)

    def test_pdist_matching_mtica2(self):
        "Tests matching(*,*) with mtica example #2."
        m = matching(numpy.array([1, 0, 1]),
                     numpy.array([1, 1, 0]))
        m2 = matching(numpy.array([1, 0, 1], dtype=numpy.bool),
                      numpy.array([1, 1, 0], dtype=numpy.bool))
        self.failUnless(numpy.abs(m - (2.0/3.0)) <= 1e-10)
        self.failUnless(numpy.abs(m2 - (2.0/3.0)) <= 1e-10)

    def test_pdist_jaccard_mtica1(self):
        "Tests jaccard(*,*) with mtica example #1."
        m = jaccard(numpy.array([1, 0, 1, 1, 0]),
                    numpy.array([1, 1, 0, 1, 1]))
        m2 = jaccard(numpy.array([1, 0, 1, 1, 0], dtype=numpy.bool),
                     numpy.array([1, 1, 0, 1, 1], dtype=numpy.bool))
        self.failUnless(numpy.abs(m - 0.6) <= 1e-10)
        self.failUnless(numpy.abs(m2 - 0.6) <= 1e-10)

    def test_pdist_jaccard_mtica2(self):
        "Tests jaccard(*,*) with mtica example #2."
        m = jaccard(numpy.array([1, 0, 1]),
                    numpy.array([1, 1, 0]))
        m2 = jaccard(numpy.array([1, 0, 1], dtype=numpy.bool),
                     numpy.array([1, 1, 0], dtype=numpy.bool))
        self.failUnless(numpy.abs(m - (2.0/3.0)) <= 1e-10)
        self.failUnless(numpy.abs(m2 - (2.0/3.0)) <= 1e-10)

    def test_pdist_yule_mtica1(self):
        "Tests yule(*,*) with mtica example #1."
        m = yule(numpy.array([1, 0, 1, 1, 0]),
                 numpy.array([1, 1, 0, 1, 1]))
        m2 = yule(numpy.array([1, 0, 1, 1, 0], dtype=numpy.bool),
                  numpy.array([1, 1, 0, 1, 1], dtype=numpy.bool))
        print m
        self.failUnless(numpy.abs(m - 2.0) <= 1e-10)
        self.failUnless(numpy.abs(m2 - 2.0) <= 1e-10)

    def test_pdist_yule_mtica2(self):
        "Tests yule(*,*) with mtica example #2."
        m = yule(numpy.array([1, 0, 1]),
                 numpy.array([1, 1, 0]))
        m2 = yule(numpy.array([1, 0, 1], dtype=numpy.bool),
                  numpy.array([1, 1, 0], dtype=numpy.bool))
        print m
        self.failUnless(numpy.abs(m - 2.0) <= 1e-10)
        self.failUnless(numpy.abs(m2 - 2.0) <= 1e-10)

    def test_pdist_dice_mtica1(self):
        "Tests dice(*,*) with mtica example #1."
        m = dice(numpy.array([1, 0, 1, 1, 0]),
                 numpy.array([1, 1, 0, 1, 1]))
        m2 = dice(numpy.array([1, 0, 1, 1, 0], dtype=numpy.bool),
                  numpy.array([1, 1, 0, 1, 1], dtype=numpy.bool))
        print m
        self.failUnless(numpy.abs(m - (3.0/7.0)) <= 1e-10)
        self.failUnless(numpy.abs(m2 - (3.0/7.0)) <= 1e-10)

    def test_pdist_dice_mtica2(self):
        "Tests dice(*,*) with mtica example #2."
        m = dice(numpy.array([1, 0, 1]),
                 numpy.array([1, 1, 0]))
        m2 = dice(numpy.array([1, 0, 1], dtype=numpy.bool),
                  numpy.array([1, 1, 0], dtype=numpy.bool))
        print m
        self.failUnless(numpy.abs(m - 0.5) <= 1e-10)
        self.failUnless(numpy.abs(m2 - 0.5) <= 1e-10)

    def test_pdist_sokalsneath_mtica1(self):
        "Tests sokalsneath(*,*) with mtica example #1."
        m = sokalsneath(numpy.array([1, 0, 1, 1, 0]),
                        numpy.array([1, 1, 0, 1, 1]))
        m2 = sokalsneath(numpy.array([1, 0, 1, 1, 0], dtype=numpy.bool),
                         numpy.array([1, 1, 0, 1, 1], dtype=numpy.bool))
        print m
        self.failUnless(numpy.abs(m - (3.0/4.0)) <= 1e-10)
        self.failUnless(numpy.abs(m2 - (3.0/4.0)) <= 1e-10)

    def test_pdist_sokalsneath_mtica2(self):
        "Tests sokalsneath(*,*) with mtica example #2."
        m = sokalsneath(numpy.array([1, 0, 1]),
                        numpy.array([1, 1, 0]))
        m2 = sokalsneath(numpy.array([1, 0, 1], dtype=numpy.bool),
                         numpy.array([1, 1, 0], dtype=numpy.bool))
        print m
        self.failUnless(numpy.abs(m - (4.0/5.0)) <= 1e-10)
        self.failUnless(numpy.abs(m2 - (4.0/5.0)) <= 1e-10)

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

