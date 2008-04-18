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
from scipy.cluster.hierarchy import pdist

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
              "pdist-chebychev-ml-iris.txt"]

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
        X = numpy.zeros((10, 10), dtype='float32')
        try:
            pdist(X, 'euclidean')
        except TypeError:
            pass
        except:
            self.fail("float32 observation matrices should generate an error in pdist.")

    def test_pdist_raises_type_error_float96(self):
        "Testing whether passing a float96 observation array generates an exception."
        X = numpy.zeros((10, 10), dtype='float96')
        try:
            pdist(X, 'euclidean')
        except TypeError:
            pass
        except:
            self.fail("float96 observation matrices should generate an error in pdist.")

    def test_pdist_var_raises_type_error_float32(self):
        "Testing whether passing a float32 variance matrix generates an exception."
        X = numpy.zeros((10, 10))
        V = numpy.zeros((10, 10), dtype='float32')
        try:
            pdist(X, 'seuclidean', V=V)
        except TypeError:
            pass
        except:
            self.fail("float32 V matrices should generate an error in pdist('seuclidean').")

    def test_pdist_var_raises_type_error_float96(self):
        "Testing whether passing a float96 variance matrix generates an exception."
        X = numpy.zeros((10, 10))
        V = numpy.zeros((10, 10), dtype='float96')

        try:
            pdist(X, 'seuclidean', V=V)
        except TypeError:
            pass
        except:
            self.fail("float96 matrices should generate an error in pdist('seuclidean').")

    def test_pdist_ivar_raises_type_error_float32(self):
        "Testing whether passing a float32 variance matrix generates an exception."
        X = numpy.zeros((10, 10))
        VI = numpy.zeros((10, 10), dtype='float32')
        try:
            pdist(X, 'mahalanobis', VI=VI)
        except TypeError:
            pass
        except:
            self.fail("float32 matrices should generate an error in pdist('mahalanobis').")

    def test_pdist_ivar_raises_type_error_float96(self):
        "Testing whether passing a float96 variance matrix generates an exception."
        X = numpy.zeros((10, 10))
        VI = numpy.zeros((10, 10), dtype='float96')

        try:
            pdist(X, 'mahalanobis', VI=VI)
        except TypeError:
            pass
        except:
            self.fail("float96 matrices should generate an error in pdist('mahalanobis').")

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

def within_tol(a, b, tol):
    return numpy.abs(a - b).max() < tol

if __name__ == "__main__":
    nose.run(argv=['', __file__])
