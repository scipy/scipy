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
from scipy.spatial.distance import squareform, pdist, cdist, matching, \
                                   jaccard, dice, sokalsneath, rogerstanimoto, \
                                   russellrao, yule, num_obs_y, num_obs_dm, \
                                   is_valid_dm, is_valid_y

_filenames = ["iris.txt",
              "cdist-X1.txt",
              "cdist-X2.txt",
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
              "random-bool-data.txt"]

_tdist = np.array([[0,    662,  877,  255,  412,  996],
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
    "Loading test data files for the scipy.spatial.distance tests."
    for fn in _filenames:
        name = fn.replace(".txt", "").replace("-ml", "")
        fqfn = os.path.join(os.path.dirname(__file__), fn)
        eo[name] = np.loadtxt(open(fqfn))
        #print "%s: %s   %s" % (name, str(eo[name].shape), str(eo[name].dtype))
    eo['pdist-boolean-inp'] = np.bool_(eo['pdist-boolean-inp'])

load_testing_files()

#print eo.keys()


#print np.abs(Y_test2 - Y_right).max()
#print np.abs(Y_test1 - Y_right).max()

class TestCdist(TestCase):
    """
    Test suite for the cdist function.
    """

    def test_cdist_euclidean_random(self):
        "Tests cdist(X, 'euclidean') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'euclidean')
        Y2 = cdist(X1, X2, 'test_euclidean')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))
        
    def test_cdist_euclidean_random_unicode(self):
        "Tests cdist(X, u'euclidean') using unicode metric string"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, u'euclidean')
        Y2 = cdist(X1, X2, u'test_euclidean')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))
        
    def test_cdist_sqeuclidean_random(self):
        "Tests cdist(X, 'sqeuclidean') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'sqeuclidean')
        Y2 = cdist(X1, X2, 'test_sqeuclidean')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_cityblock_random(self):
        "Tests cdist(X, 'cityblock') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'cityblock')
        Y2 = cdist(X1, X2, 'test_cityblock')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_hamming_double_random(self):
        "Tests cdist(X, 'hamming') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'hamming')
        Y2 = cdist(X1, X2, 'test_hamming')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_hamming_bool_random(self):
        "Tests cdist(X, 'hamming') on random boolean data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'hamming')
        Y2 = cdist(X1, X2, 'test_hamming')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_jaccard_double_random(self):
        "Tests cdist(X, 'jaccard') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'jaccard')
        Y2 = cdist(X1, X2, 'test_jaccard')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_jaccard_bool_random(self):
        "Tests cdist(X, 'jaccard') on random boolean data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'jaccard')
        Y2 = cdist(X1, X2, 'test_jaccard')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_chebychev_random(self):
        "Tests cdist(X, 'chebychev') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'chebychev')
        Y2 = cdist(X1, X2, 'test_chebychev')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_minkowski_random_p3d8(self):
        "Tests cdist(X, 'minkowski') on random data. (p=3.8)"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'minkowski', p=3.8)
        Y2 = cdist(X1, X2, 'test_minkowski', p=3.8)
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_minkowski_random_p4d6(self):
        "Tests cdist(X, 'minkowski') on random data. (p=4.6)"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'minkowski', p=4.6)
        Y2 = cdist(X1, X2, 'test_minkowski', p=4.6)
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_minkowski_random_p1d23(self):
        "Tests cdist(X, 'minkowski') on random data. (p=1.23)"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'minkowski', p=1.23)
        Y2 = cdist(X1, X2, 'test_minkowski', p=1.23)
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))


    def test_cdist_wminkowski_random_p3d8(self):
        "Tests cdist(X, 'wminkowski') on random data. (p=3.8)"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        w = 1.0 / X1.std(axis=0)
        Y1 = cdist(X1, X2, 'wminkowski', p=3.8, w=w)
        Y2 = cdist(X1, X2, 'test_wminkowski', p=3.8, w=w)
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_wminkowski_random_p4d6(self):
        "Tests cdist(X, 'wminkowski') on random data. (p=4.6)"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        w = 1.0 / X1.std(axis=0)
        Y1 = cdist(X1, X2, 'wminkowski', p=4.6, w=w)
        Y2 = cdist(X1, X2, 'test_wminkowski', p=4.6, w=w)
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_wminkowski_random_p1d23(self):
        "Tests cdist(X, 'wminkowski') on random data. (p=1.23)"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        w = 1.0 / X1.std(axis=0)
        Y1 = cdist(X1, X2, 'wminkowski', p=1.23, w=w)
        Y2 = cdist(X1, X2, 'test_wminkowski', p=1.23, w=w)
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))


    def test_cdist_seuclidean_random(self):
        "Tests cdist(X, 'seuclidean') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'seuclidean')
        Y2 = cdist(X1, X2, 'test_seuclidean')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_sqeuclidean_random(self):
        "Tests cdist(X, 'sqeuclidean') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'sqeuclidean')
        Y2 = cdist(X1, X2, 'test_sqeuclidean')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_cosine_random(self):
        "Tests cdist(X, 'cosine') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'cosine')
        Y2 = cdist(X1, X2, 'test_cosine')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_correlation_random(self):
        "Tests cdist(X, 'correlation') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'correlation')
        Y2 = cdist(X1, X2, 'test_correlation')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_mahalanobis_random(self):
        "Tests cdist(X, 'mahalanobis') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'mahalanobis')
        Y2 = cdist(X1, X2, 'test_mahalanobis')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_canberra_random(self):
        "Tests cdist(X, 'canberra') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'canberra')
        Y2 = cdist(X1, X2, 'test_canberra')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_braycurtis_random(self):
        "Tests cdist(X, 'braycurtis') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'braycurtis')
        Y2 = cdist(X1, X2, 'test_braycurtis')
        if verbose > 2:
            print Y1, Y2
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_yule_random(self):
        "Tests cdist(X, 'yule') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'yule')
        Y2 = cdist(X1, X2, 'test_yule')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_matching_random(self):
        "Tests cdist(X, 'matching') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'matching')
        Y2 = cdist(X1, X2, 'test_matching')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_kulsinski_random(self):
        "Tests cdist(X, 'kulsinski') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'kulsinski')
        Y2 = cdist(X1, X2, 'test_kulsinski')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_dice_random(self):
        "Tests cdist(X, 'dice') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'dice')
        Y2 = cdist(X1, X2, 'test_dice')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_rogerstanimoto_random(self):
        "Tests cdist(X, 'rogerstanimoto') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'rogerstanimoto')
        Y2 = cdist(X1, X2, 'test_rogerstanimoto')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_russellrao_random(self):
        "Tests cdist(X, 'russellrao') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'russellrao')
        Y2 = cdist(X1, X2, 'test_russellrao')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_sokalmichener_random(self):
        "Tests cdist(X, 'sokalmichener') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'sokalmichener')
        Y2 = cdist(X1, X2, 'test_sokalmichener')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

    def test_cdist_sokalsneath_random(self):
        "Tests cdist(X, 'sokalsneath') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'sokalsneath')
        Y2 = cdist(X1, X2, 'test_sokalsneath')
        if verbose > 2:
            print (Y1-Y2).max()
        self.failUnless(within_tol(Y1, Y2, eps))

class TestPdist(TestCase):
    """
    Test suite for the pdist function.
    """

    ################### pdist: euclidean
    def test_pdist_euclidean_random(self):
        "Tests pdist(X, 'euclidean') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-euclidean']

        Y_test1 = pdist(X, 'euclidean')
        self.failUnless(within_tol(Y_test1, Y_right, eps))
        
    def test_pdist_euclidean_random(self):
        "Tests pdist(X, 'euclidean') with unicode metric string"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-euclidean']

        Y_test1 = pdist(X, u'euclidean')
        self.failUnless(within_tol(Y_test1, Y_right, eps))
        
    def test_pdist_euclidean_random_float32(self):
        "Tests pdist(X, 'euclidean') on random data (float32)."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['pdist-double-inp'])
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

    def test_pdist_euclidean_iris_double(self):
        "Tests pdist(X, 'euclidean') on the Iris data set."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-euclidean-iris']

        Y_test1 = pdist(X, 'euclidean')
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_euclidean_iris_float32(self):
        "Tests pdist(X, 'euclidean') on the Iris data set. (float32)"
        eps = 1e-06
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-euclidean-iris']

        Y_test1 = pdist(X, 'euclidean')
        if verbose > 2:
            print np.abs(Y_right - Y_test1).max()
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

    def test_pdist_seuclidean_random_float32(self):
        "Tests pdist(X, 'seuclidean') on random data (float32)."
        eps = 1e-05
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['pdist-double-inp'])
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

    def test_pdist_seuclidean_iris_float32(self):
        "Tests pdist(X, 'seuclidean') on the Iris data set (float32)."
        eps = 1e-05
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['iris'])
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

    def test_pdist_cosine_random_float32(self):
        "Tests pdist(X, 'cosine') on random data. (float32)"
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['pdist-double-inp'])
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
        #print "cosine-iris", np.abs(Y_test1 - Y_right).max()

    def test_pdist_cosine_iris_float32(self):
        "Tests pdist(X, 'cosine') on the Iris data set."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-cosine-iris']

        Y_test1 = pdist(X, 'cosine')
        if verbose > 2:
            print np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))
        #print "cosine-iris", np.abs(Y_test1 - Y_right).max()

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
        #print "cityblock", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_cityblock_random_float32(self):
        "Tests pdist(X, 'cityblock') on random data. (float32)"
        eps = 1e-06
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-cityblock']
        Y_test1 = pdist(X, 'cityblock')
        #print "cityblock", np.abs(Y_test1 - Y_right).max()
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
        #print "cityblock-iris", np.abs(Y_test1 - Y_right).max()

    def test_pdist_cityblock_iris_float32(self):
        "Tests pdist(X, 'cityblock') on the Iris data set. (float32)"
        eps = 1e-06
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-cityblock-iris']

        Y_test1 = pdist(X, 'cityblock')
        if verbose > 2:
            print "cityblock-iris-float32", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

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
        #print "correlation", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_correlation_random_float32(self):
        "Tests pdist(X, 'correlation') on random data. (float32)"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-correlation']

        Y_test1 = pdist(X, 'correlation')
        #print "correlation", np.abs(Y_test1 - Y_right).max()
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
        #print "correlation-iris", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_correlation_iris_float32(self):
        "Tests pdist(X, 'correlation') on the Iris data set. (float32)"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = np.float32(eo['pdist-correlation-iris'])

        Y_test1 = pdist(X, 'correlation')
        if verbose > 2:
            print "correlation-iris", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_correlation_iris_nonC(self):
        "Tests pdist(X, 'test_correlation') [the non-C implementation] on the Iris data set."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-correlation-iris']
        Y_test2 = pdist(X, 'test_correlation')
        #print "test-correlation-iris", np.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################# minkowski

    def test_pdist_minkowski_random(self):
        "Tests pdist(X, 'minkowski') on random data."
        eps = 1e-05
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-minkowski-3.2']

        Y_test1 = pdist(X, 'minkowski', 3.2)
        #print "minkowski", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_minkowski_random_float32(self):
        "Tests pdist(X, 'minkowski') on random data. (float32)"
        eps = 1e-05
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-minkowski-3.2']

        Y_test1 = pdist(X, 'minkowski', 3.2)
        #print "minkowski", np.abs(Y_test1 - Y_right).max()
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
        #print "minkowski-iris-3.2", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_minkowski_iris_float32(self):
        "Tests pdist(X, 'minkowski') on iris data. (float32)"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-minkowski-3.2-iris']
        Y_test1 = pdist(X, 'minkowski', 3.2)
        #print "minkowski-iris-3.2", np.abs(Y_test1 - Y_right).max()
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
        #print "minkowski-iris-5.8", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_minkowski_iris_float32(self):
        "Tests pdist(X, 'minkowski') on iris data. (float32)"
        eps = 1e-06
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-minkowski-5.8-iris']

        Y_test1 = pdist(X, 'minkowski', 5.8)
        if verbose > 2:
            print "minkowski-iris-5.8", np.abs(Y_test1 - Y_right).max()
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
        #print "hamming", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_hamming_random_float32(self):
        "Tests pdist(X, 'hamming') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']

        Y_test1 = pdist(X, 'hamming')
        #print "hamming", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_hamming_random_nonC(self):
        "Tests pdist(X, 'test_hamming') [the non-C implementation] on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-hamming']
        Y_test2 = pdist(X, 'test_hamming')
        #print "test-hamming", np.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: hamming (double)
    def test_pdist_dhamming_random(self):
        "Tests pdist(X, 'hamming') on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = np.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']
        Y_test1 = pdist(X, 'hamming')
        #print "hamming", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_dhamming_random_float32(self):
        "Tests pdist(X, 'hamming') on random data. (float32)"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']
        Y_test1 = pdist(X, 'hamming')
        #print "hamming", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_dhamming_random_nonC(self):
        "Tests pdist(X, 'test_hamming') [the non-C implementation] on random data."
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = np.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']
        Y_test2 = pdist(X, 'test_hamming')
        #print "test-hamming", np.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: jaccard
    def test_pdist_jaccard_random(self):
        "Tests pdist(X, 'jaccard') on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-jaccard']

        Y_test1 = pdist(X, 'jaccard')
        #print "jaccard", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_jaccard_random_float32(self):
        "Tests pdist(X, 'jaccard') on random data. (float32)"
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']

        Y_test1 = pdist(X, 'jaccard')
        #print "jaccard", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_jaccard_random_nonC(self):
        "Tests pdist(X, 'test_jaccard') [the non-C implementation] on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-jaccard']
        Y_test2 = pdist(X, 'test_jaccard')
        #print "test-jaccard", np.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: jaccard (double)
    def test_pdist_djaccard_random(self):
        "Tests pdist(X, 'jaccard') on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = np.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']

        Y_test1 = pdist(X, 'jaccard')
        #print "jaccard", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_djaccard_random_float32(self):
        "Tests pdist(X, 'jaccard') on random data. (float32)"
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']

        Y_test1 = pdist(X, 'jaccard')
        #print "jaccard", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_djaccard_random_nonC(self):
        "Tests pdist(X, 'test_jaccard') [the non-C implementation] on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = np.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']
        Y_test2 = pdist(X, 'test_jaccard')
        #print "test-jaccard", np.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    ################### pdist: chebychev
    def test_pdist_chebychev_random(self):
        "Tests pdist(X, 'chebychev') on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-chebychev']

        Y_test1 = pdist(X, 'chebychev')
        #print "chebychev", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_chebychev_random_float32(self):
        "Tests pdist(X, 'chebychev') on random data. (float32)"
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-chebychev']

        Y_test1 = pdist(X, 'chebychev')
        if verbose > 2:
            print "chebychev", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_chebychev_random_nonC(self):
        "Tests pdist(X, 'test_chebychev') [the non-C implementation] on random data."
        eps = 1e-08
        # Get the data: the input matrix and the right output.
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-chebychev']
        Y_test2 = pdist(X, 'test_chebychev')
        #print "test-chebychev", np.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    def test_pdist_chebychev_iris(self):
        "Tests pdist(X, 'chebychev') on the Iris data set."
        eps = 1e-15
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-chebychev-iris']
        Y_test1 = pdist(X, 'chebychev')
        #print "chebychev-iris", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_chebychev_iris_float32(self):
        "Tests pdist(X, 'chebychev') on the Iris data set. (float32)"
        eps = 1e-06
        # Get the data: the input matrix and the right output.
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-chebychev-iris']
        Y_test1 = pdist(X, 'chebychev')
        if verbose > 2:
            print "chebychev-iris", np.abs(Y_test1 - Y_right).max()
        self.failUnless(within_tol(Y_test1, Y_right, eps))

    def test_pdist_chebychev_iris_nonC(self):
        "Tests pdist(X, 'test_chebychev') [the non-C implementation] on the Iris data set."
        eps = 1e-15
        # Get the data: the input matrix and the right output.
        X = eo['iris']
        Y_right = eo['pdist-chebychev-iris']
        Y_test2 = pdist(X, 'test_chebychev')
        #print "test-chebychev-iris", np.abs(Y_test2 - Y_right).max()
        self.failUnless(within_tol(Y_test2, Y_right, eps))

    def test_pdist_matching_mtica1(self):
        "Tests matching(*,*) with mtica example #1 (nums)."
        m = matching(np.array([1, 0, 1, 1, 0]),
                     np.array([1, 1, 0, 1, 1]))
        m2 = matching(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                      np.array([1, 1, 0, 1, 1], dtype=np.bool))
        self.failUnless(np.abs(m - 0.6) <= 1e-10)
        self.failUnless(np.abs(m2 - 0.6) <= 1e-10)

    def test_pdist_matching_mtica2(self):
        "Tests matching(*,*) with mtica example #2."
        m = matching(np.array([1, 0, 1]),
                     np.array([1, 1, 0]))
        m2 = matching(np.array([1, 0, 1], dtype=np.bool),
                      np.array([1, 1, 0], dtype=np.bool))
        self.failUnless(np.abs(m - (2.0/3.0)) <= 1e-10)
        self.failUnless(np.abs(m2 - (2.0/3.0)) <= 1e-10)

    def test_pdist_matching_match(self):
        "Tests pdist(X, 'matching') to see if the two implementations match on random boolean input data."
        D = eo['random-bool-data']
        B = np.bool_(D)
        if verbose > 2:
            print B.shape, B.dtype
        eps = 1e-10
        y1 = pdist(B, "matching")
        y2 = pdist(B, "test_matching")
        y3 = pdist(D, "test_matching")
        if verbose > 2:
            print np.abs(y1-y2).max()
            print np.abs(y1-y3).max()
        self.failUnless(within_tol(y1, y2, eps))
        self.failUnless(within_tol(y2, y3, eps))

    def test_pdist_jaccard_mtica1(self):
        "Tests jaccard(*,*) with mtica example #1."
        m = jaccard(np.array([1, 0, 1, 1, 0]),
                    np.array([1, 1, 0, 1, 1]))
        m2 = jaccard(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                     np.array([1, 1, 0, 1, 1], dtype=np.bool))
        self.failUnless(np.abs(m - 0.6) <= 1e-10)
        self.failUnless(np.abs(m2 - 0.6) <= 1e-10)

    def test_pdist_jaccard_mtica2(self):
        "Tests jaccard(*,*) with mtica example #2."
        m = jaccard(np.array([1, 0, 1]),
                    np.array([1, 1, 0]))
        m2 = jaccard(np.array([1, 0, 1], dtype=np.bool),
                     np.array([1, 1, 0], dtype=np.bool))
        self.failUnless(np.abs(m - (2.0/3.0)) <= 1e-10)
        self.failUnless(np.abs(m2 - (2.0/3.0)) <= 1e-10)

    def test_pdist_jaccard_match(self):
        "Tests pdist(X, 'jaccard') to see if the two implementations match on random double input data."
        D = eo['random-bool-data']
        if verbose > 2:
            print D.shape, D.dtype
        eps = 1e-10
        y1 = pdist(D, "jaccard")
        y2 = pdist(D, "test_jaccard")
        y3 = pdist(np.bool_(D), "test_jaccard")
        if verbose > 2:
            print np.abs(y1-y2).max()
            print np.abs(y2-y3).max()
        self.failUnless(within_tol(y1, y2, eps))
        self.failUnless(within_tol(y2, y3, eps))

    def test_pdist_yule_mtica1(self):
        "Tests yule(*,*) with mtica example #1."
        m = yule(np.array([1, 0, 1, 1, 0]),
                 np.array([1, 1, 0, 1, 1]))
        m2 = yule(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                  np.array([1, 1, 0, 1, 1], dtype=np.bool))
        if verbose > 2:
            print m
        self.failUnless(np.abs(m - 2.0) <= 1e-10)
        self.failUnless(np.abs(m2 - 2.0) <= 1e-10)

    def test_pdist_yule_mtica2(self):
        "Tests yule(*,*) with mtica example #2."
        m = yule(np.array([1, 0, 1]),
                 np.array([1, 1, 0]))
        m2 = yule(np.array([1, 0, 1], dtype=np.bool),
                  np.array([1, 1, 0], dtype=np.bool))
        if verbose > 2:
            print m
        self.failUnless(np.abs(m - 2.0) <= 1e-10)
        self.failUnless(np.abs(m2 - 2.0) <= 1e-10)

    def test_pdist_yule_match(self):
        "Tests pdist(X, 'yule') to see if the two implementations match on random double input data."
        D = eo['random-bool-data']
        if verbose > 2:
            print D.shape, D.dtype
        eps = 1e-10
        y1 = pdist(D, "yule")
        y2 = pdist(D, "test_yule")
        y3 = pdist(np.bool_(D), "test_yule")
        if verbose > 2:
            print np.abs(y1-y2).max()
            print np.abs(y2-y3).max()
        self.failUnless(within_tol(y1, y2, eps))
        self.failUnless(within_tol(y2, y3, eps))

    def test_pdist_dice_mtica1(self):
        "Tests dice(*,*) with mtica example #1."
        m = dice(np.array([1, 0, 1, 1, 0]),
                 np.array([1, 1, 0, 1, 1]))
        m2 = dice(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                  np.array([1, 1, 0, 1, 1], dtype=np.bool))
        if verbose > 2:
            print m
        self.failUnless(np.abs(m - (3.0/7.0)) <= 1e-10)
        self.failUnless(np.abs(m2 - (3.0/7.0)) <= 1e-10)

    def test_pdist_dice_mtica2(self):
        "Tests dice(*,*) with mtica example #2."
        m = dice(np.array([1, 0, 1]),
                 np.array([1, 1, 0]))
        m2 = dice(np.array([1, 0, 1], dtype=np.bool),
                  np.array([1, 1, 0], dtype=np.bool))
        if verbose > 2:
            print m
        self.failUnless(np.abs(m - 0.5) <= 1e-10)
        self.failUnless(np.abs(m2 - 0.5) <= 1e-10)

    def test_pdist_dice_match(self):
        "Tests pdist(X, 'dice') to see if the two implementations match on random double input data."
        D = eo['random-bool-data']
        if verbose > 2:
            print D.shape, D.dtype
        eps = 1e-10
        y1 = pdist(D, "dice")
        y2 = pdist(D, "test_dice")
        y3 = pdist(D, "test_dice")
        if verbose > 2:
            print np.abs(y1-y2).max()
            print np.abs(y2-y3).max()
        self.failUnless(within_tol(y1, y2, eps))
        self.failUnless(within_tol(y2, y3, eps))

    def test_pdist_sokalsneath_mtica1(self):
        "Tests sokalsneath(*,*) with mtica example #1."
        m = sokalsneath(np.array([1, 0, 1, 1, 0]),
                        np.array([1, 1, 0, 1, 1]))
        m2 = sokalsneath(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                         np.array([1, 1, 0, 1, 1], dtype=np.bool))
        if verbose > 2:
            print m
        self.failUnless(np.abs(m - (3.0/4.0)) <= 1e-10)
        self.failUnless(np.abs(m2 - (3.0/4.0)) <= 1e-10)

    def test_pdist_sokalsneath_mtica2(self):
        "Tests sokalsneath(*,*) with mtica example #2."
        m = sokalsneath(np.array([1, 0, 1]),
                        np.array([1, 1, 0]))
        m2 = sokalsneath(np.array([1, 0, 1], dtype=np.bool),
                         np.array([1, 1, 0], dtype=np.bool))
        if verbose > 2:
            print m
        self.failUnless(np.abs(m - (4.0/5.0)) <= 1e-10)
        self.failUnless(np.abs(m2 - (4.0/5.0)) <= 1e-10)

    def test_pdist_sokalsneath_match(self):
        "Tests pdist(X, 'sokalsneath') to see if the two implementations match on random double input data."
        D = eo['random-bool-data']
        if verbose > 2:
            print D.shape, D.dtype
        eps = 1e-10
        y1 = pdist(D, "sokalsneath")
        y2 = pdist(D, "test_sokalsneath")
        y3 = pdist(np.bool_(D), "test_sokalsneath")
        if verbose > 2:
            print np.abs(y1-y2).max()
            print np.abs(y2-y3).max()
        self.failUnless(within_tol(y1, y2, eps))
        self.failUnless(within_tol(y2, y3, eps))

    def test_pdist_rogerstanimoto_mtica1(self):
        "Tests rogerstanimoto(*,*) with mtica example #1."
        m = rogerstanimoto(np.array([1, 0, 1, 1, 0]),
                           np.array([1, 1, 0, 1, 1]))
        m2 = rogerstanimoto(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                            np.array([1, 1, 0, 1, 1], dtype=np.bool))
        if verbose > 2:
            print m
        self.failUnless(np.abs(m - (3.0/4.0)) <= 1e-10)
        self.failUnless(np.abs(m2 - (3.0/4.0)) <= 1e-10)

    def test_pdist_rogerstanimoto_mtica2(self):
        "Tests rogerstanimoto(*,*) with mtica example #2."
        m = rogerstanimoto(np.array([1, 0, 1]),
                           np.array([1, 1, 0]))
        m2 = rogerstanimoto(np.array([1, 0, 1], dtype=np.bool),
                            np.array([1, 1, 0], dtype=np.bool))
        if verbose > 2:
            print m
        self.failUnless(np.abs(m - (4.0/5.0)) <= 1e-10)
        self.failUnless(np.abs(m2 - (4.0/5.0)) <= 1e-10)

    def test_pdist_rogerstanimoto_match(self):
        "Tests pdist(X, 'rogerstanimoto') to see if the two implementations match on random double input data."
        D = eo['random-bool-data']
        if verbose > 2:
            print D.shape, D.dtype
        eps = 1e-10
        y1 = pdist(D, "rogerstanimoto")
        y2 = pdist(D, "test_rogerstanimoto")
        y3 = pdist(np.bool_(D), "test_rogerstanimoto")
        if verbose > 2:
            print np.abs(y1-y2).max()
            print np.abs(y2-y3).max()
        self.failUnless(within_tol(y1, y2, eps))
        self.failUnless(within_tol(y2, y3, eps))

    def test_pdist_russellrao_mtica1(self):
        "Tests russellrao(*,*) with mtica example #1."
        m = russellrao(np.array([1, 0, 1, 1, 0]),
                       np.array([1, 1, 0, 1, 1]))
        m2 = russellrao(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                        np.array([1, 1, 0, 1, 1], dtype=np.bool))
        if verbose > 2:
            print m
        self.failUnless(np.abs(m - (3.0/5.0)) <= 1e-10)
        self.failUnless(np.abs(m2 - (3.0/5.0)) <= 1e-10)

    def test_pdist_russellrao_mtica2(self):
        "Tests russellrao(*,*) with mtica example #2."
        m = russellrao(np.array([1, 0, 1]),
                       np.array([1, 1, 0]))
        m2 = russellrao(np.array([1, 0, 1], dtype=np.bool),
                        np.array([1, 1, 0], dtype=np.bool))
        if verbose > 2:
            print m
        self.failUnless(np.abs(m - (2.0/3.0)) <= 1e-10)
        self.failUnless(np.abs(m2 - (2.0/3.0)) <= 1e-10)

    def test_pdist_russellrao_match(self):
        "Tests pdist(X, 'russellrao') to see if the two implementations match on random double input data."
        D = eo['random-bool-data']
        if verbose > 2:
            print D.shape, D.dtype
        eps = 1e-10
        y1 = pdist(D, "russellrao")
        y2 = pdist(D, "test_russellrao")
        y3 = pdist(np.bool_(D), "test_russellrao")
        if verbose > 2:
            print np.abs(y1-y2).max()
            print np.abs(y2-y3).max()
        self.failUnless(within_tol(y1, y2, eps))
        self.failUnless(within_tol(y2, y3, eps))

    def test_pdist_sokalmichener_match(self):
        "Tests pdist(X, 'sokalmichener') to see if the two implementations match on random double input data."
        D = eo['random-bool-data']
        if verbose > 2:
            print D.shape, D.dtype
        eps = 1e-10
        y1 = pdist(D, "sokalmichener")
        y2 = pdist(D, "test_sokalmichener")
        y3 = pdist(np.bool_(D), "test_sokalmichener")
        if verbose > 2:
            print np.abs(y1-y2).max()
            print np.abs(y2-y3).max()
        self.failUnless(within_tol(y1, y2, eps))
        self.failUnless(within_tol(y2, y3, eps))

    def test_pdist_kulsinski_match(self):
        "Tests pdist(X, 'kulsinski') to see if the two implementations match on random double input data."
        D = eo['random-bool-data']
        if verbose > 2:
            print D.shape, D.dtype
        eps = 1e-10
        y1 = pdist(D, "kulsinski")
        y2 = pdist(D, "test_kulsinski")
        y3 = pdist(np.bool_(D), "test_kulsinski")
        if verbose > 2:
            print np.abs(y1-y2).max()
        self.failUnless(within_tol(y1, y2, eps))

    def test_pdist_canberra_match(self):
        "Tests pdist(X, 'canberra') to see if the two implementations match on the Iris data set."
        D = eo['iris']
        if verbose > 2:
            print D.shape, D.dtype
        eps = 1e-10
        y1 = pdist(D, "canberra")
        y2 = pdist(D, "test_canberra")
        if verbose > 2:
            print np.abs(y1-y2).max()
        self.failUnless(within_tol(y1, y2, eps))

    def test_pdist_canberra_ticket_711(self):
        "Tests pdist(X, 'canberra') to see if Canberra gives the right result as reported in Scipy bug report 711."
        eps = 1e-8
        pdist_y = pdist(([3.3], [3.4]), "canberra")
        right_y = 0.01492537
        if verbose > 2:
            print np.abs(pdist_y-right_y).max()
        self.failUnless(within_tol(pdist_y, right_y, eps))

def within_tol(a, b, tol):
    return np.abs(a - b).max() < tol


class TestSquareForm(TestCase):

    ################### squareform
    def test_squareform_empty_matrix(self):
        "Tests squareform on an empty matrix."
        A = np.zeros((0,0))
        rA = squareform(np.array(A, dtype='double'))
        self.failUnless(rA.shape == (0,))

    def test_squareform_empty_vector(self):
        "Tests squareform on an empty vector."
        v = np.zeros((0,))
        rv = squareform(np.array(v, dtype='double'))
        self.failUnless(rv.shape == (1,1))
        self.failUnless(rv[0, 0] == 0)

    def test_squareform_1by1_matrix(self):
        "Tests squareform on a 1x1 matrix."
        A = np.zeros((1,1))
        rA = squareform(np.array(A, dtype='double'))
        self.failUnless(rA.shape == (0,))

    def test_squareform_one_vector(self):
        "Tests squareform on a 1-D array, length=1."
        v = np.ones((1,)) * 8.3
        rv = squareform(np.array(v, dtype='double'))
        self.failUnless(rv.shape == (2,2))
        self.failUnless(rv[0,1] == 8.3)
        self.failUnless(rv[1,0] == 8.3)

    def test_squareform_2by2_matrix(self):
        "Tests squareform on a 2x2 matrix."
        A = np.zeros((2,2))
        A[0,1]=0.8
        A[1,0]=0.8
        rA = squareform(np.array(A, dtype='double'))
        self.failUnless(rA.shape == (1,))
        self.failUnless(rA[0] == 0.8)

    def test_squareform_multi_matrix(self):
        "Tests squareform on a square matrices of multiple sizes."
        for n in xrange(2, 5):
            yield self.check_squareform_multi_matrix(n)

    def check_squareform_multi_matrix(self, n):
        X = np.random.rand(n, 4)
        Y = pdist(X)
        self.failUnless(len(Y.shape) == 1)
        A = squareform(Y)
        Yr = squareform(A)
        s = A.shape
        k = 0
        if verbose >= 3:
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

class TestNumObsY(TestCase):

    def test_num_obs_y_multi_matrix(self):
        "Tests num_obs_y with observation matrices of multiple sizes."
        for n in xrange(2, 10):
            X = np.random.rand(n, 4)
            Y = pdist(X)
            #print A.shape, Y.shape, Yr.shape
            self.failUnless(num_obs_y(Y) == n)

    def test_num_obs_y_1(self):
        "Tests num_obs_y(y) on a condensed distance matrix over 1 observations. Expecting exception."
        self.failUnlessRaises(ValueError, self.check_y, 1)

    def test_num_obs_y_2(self):
        "Tests num_obs_y(y) on a condensed distance matrix over 2 observations."
        self.failUnless(self.check_y(2))

    def test_num_obs_y_3(self):
        "Tests num_obs_y(y) on a condensed distance matrix over 3 observations."
        self.failUnless(self.check_y(3))

    def test_num_obs_y_4(self):
        "Tests num_obs_y(y) on a condensed distance matrix over 4 observations."
        self.failUnless(self.check_y(4))

    def test_num_obs_y_5_10(self):
        "Tests num_obs_y(y) on a condensed distance matrix between 5 and 15 observations."
        for i in xrange(5, 16):
            self.minit(i)

    def test_num_obs_y_2_100(self):
        "Tests num_obs_y(y) on 100 improper condensed distance matrices. Expecting exception."
        a = set([])
        for n in xrange(2, 16):
            a.add(n*(n-1)/2)
        for i in xrange(5, 105):
            if i not in a:
                self.failUnlessRaises(ValueError, self.bad_y, i)

    def minit(self, n):
        self.failUnless(self.check_y(n))

    def bad_y(self, n):
        y = np.random.rand(n)
        return num_obs_y(y)

    def check_y(self, n):
        return num_obs_y(self.make_y(n)) == n

    def make_y(self, n):
        return np.random.rand((n*(n-1)/2))

class TestNumObsDM(TestCase):

    ############## num_obs_dm
    def test_num_obs_dm_multi_matrix(self):
        "Tests num_obs_dm with observation matrices of multiple sizes."
        for n in xrange(1, 10):
            X = np.random.rand(n, 4)
            Y = pdist(X)
            A = squareform(Y)
            if verbose >= 3:
                print A.shape, Y.shape
            self.failUnless(num_obs_dm(A) == n)

    def test_num_obs_dm_0(self):
        "Tests num_obs_dm(D) on a 0x0 distance matrix. Expecting exception."
        self.failUnless(self.check_D(0))

    def test_num_obs_dm_1(self):
        "Tests num_obs_dm(D) on a 1x1 distance matrix."
        self.failUnless(self.check_D(1))

    def test_num_obs_dm_2(self):
        "Tests num_obs_dm(D) on a 2x2 distance matrix."
        self.failUnless(self.check_D(2))

    def test_num_obs_dm_3(self):
        "Tests num_obs_dm(D) on a 3x3 distance matrix."
        self.failUnless(self.check_D(2))

    def test_num_obs_dm_4(self):
        "Tests num_obs_dm(D) on a 4x4 distance matrix."
        self.failUnless(self.check_D(4))

    def check_D(self, n):
        return num_obs_dm(self.make_D(n)) == n

    def make_D(self, n):
        return np.random.rand(n, n)

def is_valid_dm_throw(D):
    return is_valid_dm(D, throw=True)

class TestIsValidDM(TestCase):

    def test_is_valid_dm_int16_array_E(self):
        "Tests is_valid_dm(*) on an int16 array. Exception expected."
        D = np.zeros((5, 5), dtype='i')
        self.failUnlessRaises(TypeError, is_valid_dm_throw, (D))

    def test_is_valid_dm_int16_array_F(self):
        "Tests is_valid_dm(*) on an int16 array. False expected."
        D = np.zeros((5, 5), dtype='i')
        self.failUnless(is_valid_dm(D) == False)

    def test_is_valid_dm_improper_shape_1D_E(self):
        "Tests is_valid_dm(*) on a 1D array. Exception expected."
        D = np.zeros((5,), dtype=np.double)
        self.failUnlessRaises(ValueError, is_valid_dm_throw, (D))

    def test_is_valid_dm_improper_shape_1D_F(self):
        "Tests is_valid_dm(*) on a 1D array. False expected."
        D = np.zeros((5,), dtype=np.double)
        self.failUnless(is_valid_dm(D) == False)

    def test_is_valid_dm_improper_shape_3D_E(self):
        "Tests is_valid_dm(*) on a 3D array. Exception expected."
        D = np.zeros((3,3,3), dtype=np.double)
        self.failUnlessRaises(ValueError, is_valid_dm_throw, (D))

    def test_is_valid_dm_improper_shape_3D_F(self):
        "Tests is_valid_dm(*) on a 3D array. False expected."
        D = np.zeros((3,3,3), dtype=np.double)
        self.failUnless(is_valid_dm(D) == False)

    def test_is_valid_dm_nonzero_diagonal_E(self):
        "Tests is_valid_dm(*) on a distance matrix with a nonzero diagonal. Exception expected."
        y = np.random.rand(10)
        D = squareform(y)
        for i in xrange(0, 5):
            D[i, i] = 2.0
        self.failUnlessRaises(ValueError, is_valid_dm_throw, (D))

    def test_is_valid_dm_nonzero_diagonal_F(self):
        "Tests is_valid_dm(*) on a distance matrix with a nonzero diagonal. False expected."
        y = np.random.rand(10)
        D = squareform(y)
        for i in xrange(0, 5):
            D[i, i] = 2.0
        self.failUnless(is_valid_dm(D) == False)

    def test_is_valid_dm_assymetric_E(self):
        "Tests is_valid_dm(*) on an assymetric distance matrix. Exception expected."
        y = np.random.rand(10)
        D = squareform(y)
        D[1,3] = D[3,1] + 1
        self.failUnlessRaises(ValueError, is_valid_dm_throw, (D))

    def test_is_valid_dm_assymetric_F(self):
        "Tests is_valid_dm(*) on an assymetric distance matrix. False expected."
        y = np.random.rand(10)
        D = squareform(y)
        D[1,3] = D[3,1] + 1
        self.failUnless(is_valid_dm(D) == False)

    def test_is_valid_dm_correct_1_by_1(self):
        "Tests is_valid_dm(*) on a correct 1x1. True expected."
        D = np.zeros((1,1), dtype=np.double)
        self.failUnless(is_valid_dm(D) == True)

    def test_is_valid_dm_correct_2_by_2(self):
        "Tests is_valid_dm(*) on a correct 2x2. True expected."
        y = np.random.rand(1)
        D = squareform(y)
        self.failUnless(is_valid_dm(D) == True)

    def test_is_valid_dm_correct_3_by_3(self):
        "Tests is_valid_dm(*) on a correct 3x3. True expected."
        y = np.random.rand(3)
        D = squareform(y)
        self.failUnless(is_valid_dm(D) == True)

    def test_is_valid_dm_correct_4_by_4(self):
        "Tests is_valid_dm(*) on a correct 4x4. True expected."
        y = np.random.rand(6)
        D = squareform(y)
        self.failUnless(is_valid_dm(D) == True)

    def test_is_valid_dm_correct_5_by_5(self):
        "Tests is_valid_dm(*) on a correct 5x5. True expected."
        y = np.random.rand(10)
        D = squareform(y)
        self.failUnless(is_valid_dm(D) == True)

def is_valid_y_throw(y):
    return is_valid_y(y, throw=True)

class TestIsValidY(TestCase):

    def test_is_valid_y_int16_array_E(self):
        "Tests is_valid_y(*) on an int16 array. Exception expected."
        y = np.zeros((10,), dtype='i')
        self.failUnlessRaises(TypeError, is_valid_y_throw, (y))

    def test_is_valid_y_int16_array_F(self):
        "Tests is_valid_y(*) on an int16 array. False expected."
        y = np.zeros((10,), dtype='i')
        self.failUnless(is_valid_y(y) == False)

    def test_is_valid_y_improper_shape_2D_E(self):
        "Tests is_valid_y(*) on a 2D array. Exception expected."
        y = np.zeros((3,3,), dtype=np.double)
        self.failUnlessRaises(ValueError, is_valid_y_throw, (y))

    def test_is_valid_y_improper_shape_2D_F(self):
        "Tests is_valid_y(*) on a 2D array. False expected."
        y = np.zeros((3,3,), dtype=np.double)
        self.failUnless(is_valid_y(y) == False)

    def test_is_valid_y_improper_shape_3D_E(self):
        "Tests is_valid_y(*) on a 3D array. Exception expected."
        y = np.zeros((3,3,3), dtype=np.double)
        self.failUnlessRaises(ValueError, is_valid_y_throw, (y))

    def test_is_valid_y_improper_shape_3D_F(self):
        "Tests is_valid_y(*) on a 3D array. False expected."
        y = np.zeros((3,3,3), dtype=np.double)
        self.failUnless(is_valid_y(y) == False)

    def test_is_valid_y_correct_2_by_2(self):
        "Tests is_valid_y(*) on a correct 2x2 condensed. True expected."
        y = self.correct_n_by_n(2)
        self.failUnless(is_valid_y(y) == True)

    def test_is_valid_y_correct_3_by_3(self):
        "Tests is_valid_y(*) on a correct 3x3 condensed. True expected."
        y = self.correct_n_by_n(3)
        self.failUnless(is_valid_y(y) == True)

    def test_is_valid_y_correct_4_by_4(self):
        "Tests is_valid_y(*) on a correct 4x4 condensed. True expected."
        y = self.correct_n_by_n(4)
        self.failUnless(is_valid_y(y) == True)

    def test_is_valid_y_correct_5_by_5(self):
        "Tests is_valid_y(*) on a correct 5x5 condensed. True expected."
        y = self.correct_n_by_n(5)
        self.failUnless(is_valid_y(y) == True)

    def test_is_valid_y_2_100(self):
        "Tests is_valid_y(*) on 100 improper condensed distance matrices. Expecting exception."
        a = set([])
        for n in xrange(2, 16):
            a.add(n*(n-1)/2)
        for i in xrange(5, 105):
            if i not in a:
                self.failUnlessRaises(ValueError, self.bad_y, i)

    def bad_y(self, n):
        y = np.random.rand(n)
        return is_valid_y(y, throw=True)

    def correct_n_by_n(self, n):
        y = np.random.rand(n*(n-1)/2)
        return y

if __name__=="__main__":
    run_module_suite()
