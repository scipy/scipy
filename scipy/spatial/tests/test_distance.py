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

from __future__ import division, print_function, absolute_import

import os.path
from scipy.lib.six import xrange

import numpy as np
from numpy.linalg import norm
from numpy.testing import (verbose, TestCase, run_module_suite, assert_,
        assert_raises, assert_array_equal, assert_equal, assert_almost_equal,
        assert_allclose)

from scipy.lib.six import u

from scipy.spatial.distance import (squareform, pdist, cdist, matching,
        jaccard, dice, sokalsneath, rogerstanimoto, russellrao, yule,
        num_obs_y, num_obs_dm, is_valid_dm, is_valid_y, minkowski, wminkowski,
        euclidean, sqeuclidean, cosine, correlation, hamming, mahalanobis,
        canberra, braycurtis, sokalmichener, _validate_vector)


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

_tdist = np.array([[0, 662, 877, 255, 412, 996],
                      [662, 0, 295, 468, 268, 400],
                      [877, 295, 0, 754, 564, 138],
                      [255, 468, 754, 0, 219, 869],
                      [412, 268, 564, 219, 0, 669],
                      [996, 400, 138, 869, 669, 0]], dtype='double')

_ytdist = squareform(_tdist)

# A hashmap of expected output arrays for the tests. These arrays
# come from a list of text files, which are read prior to testing.
# Each test loads inputs and outputs from this dictionary.
eo = {}


def load_testing_files():
    for fn in _filenames:
        name = fn.replace(".txt", "").replace("-ml", "")
        fqfn = os.path.join(os.path.dirname(__file__), 'data', fn)
        fp = open(fqfn)
        eo[name] = np.loadtxt(fp)
        fp.close()
    eo['pdist-boolean-inp'] = np.bool_(eo['pdist-boolean-inp'])

load_testing_files()


class TestCdist(TestCase):

    def test_cdist_euclidean_random(self):
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'euclidean')
        Y2 = cdist(X1, X2, 'test_euclidean')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_euclidean_random_unicode(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, u('euclidean'))
        Y2 = cdist(X1, X2, u('test_euclidean'))
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_sqeuclidean_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'sqeuclidean')
        Y2 = cdist(X1, X2, 'test_sqeuclidean')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_cityblock_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'cityblock')
        Y2 = cdist(X1, X2, 'test_cityblock')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_hamming_double_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'hamming')
        Y2 = cdist(X1, X2, 'test_hamming')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_hamming_bool_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'hamming')
        Y2 = cdist(X1, X2, 'test_hamming')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_jaccard_double_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'jaccard')
        Y2 = cdist(X1, X2, 'test_jaccard')
        if verbose > 2:
            print((Y1-Y2).max())
        self.assertTrue(within_tol(Y1, Y2, eps))

    def test_cdist_jaccard_bool_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'jaccard')
        Y2 = cdist(X1, X2, 'test_jaccard')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_chebychev_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'chebychev')
        Y2 = cdist(X1, X2, 'test_chebychev')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_minkowski_random_p3d8(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'minkowski', p=3.8)
        Y2 = cdist(X1, X2, 'test_minkowski', p=3.8)
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_minkowski_random_p4d6(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'minkowski', p=4.6)
        Y2 = cdist(X1, X2, 'test_minkowski', p=4.6)
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_minkowski_random_p1d23(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'minkowski', p=1.23)
        Y2 = cdist(X1, X2, 'test_minkowski', p=1.23)
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_wminkowski_random_p3d8(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        w = 1.0 / X1.std(axis=0)
        Y1 = cdist(X1, X2, 'wminkowski', p=3.8, w=w)
        Y2 = cdist(X1, X2, 'test_wminkowski', p=3.8, w=w)
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_wminkowski_random_p4d6(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        w = 1.0 / X1.std(axis=0)
        Y1 = cdist(X1, X2, 'wminkowski', p=4.6, w=w)
        Y2 = cdist(X1, X2, 'test_wminkowski', p=4.6, w=w)
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_wminkowski_random_p1d23(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        w = 1.0 / X1.std(axis=0)
        Y1 = cdist(X1, X2, 'wminkowski', p=1.23, w=w)
        Y2 = cdist(X1, X2, 'test_wminkowski', p=1.23, w=w)
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_seuclidean_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'seuclidean')
        Y2 = cdist(X1, X2, 'test_seuclidean')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_cosine_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'cosine')
        Y2 = cdist(X1, X2, 'test_cosine')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_correlation_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'correlation')
        Y2 = cdist(X1, X2, 'test_correlation')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_mahalanobis_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'mahalanobis')
        Y2 = cdist(X1, X2, 'test_mahalanobis')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_canberra_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'canberra')
        Y2 = cdist(X1, X2, 'test_canberra')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_braycurtis_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'braycurtis')
        Y2 = cdist(X1, X2, 'test_braycurtis')
        if verbose > 2:
            print(Y1, Y2)
            print((Y1-Y2).max())
        _assert_within_tol(Y1, Y2, eps)

    def test_cdist_yule_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'yule')
        Y2 = cdist(X1, X2, 'test_yule')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_matching_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'matching')
        Y2 = cdist(X1, X2, 'test_matching')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_kulsinski_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'kulsinski')
        Y2 = cdist(X1, X2, 'test_kulsinski')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_dice_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'dice')
        Y2 = cdist(X1, X2, 'test_dice')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_rogerstanimoto_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'rogerstanimoto')
        Y2 = cdist(X1, X2, 'test_rogerstanimoto')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_russellrao_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'russellrao')
        Y2 = cdist(X1, X2, 'test_russellrao')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_sokalmichener_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'sokalmichener')
        Y2 = cdist(X1, X2, 'test_sokalmichener')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_sokalsneath_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = cdist(X1, X2, 'sokalsneath')
        Y2 = cdist(X1, X2, 'test_sokalsneath')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)


class TestPdist(TestCase):

    def test_pdist_euclidean_random(self):
        eps = 1e-07
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-euclidean']

        Y_test1 = pdist(X, 'euclidean')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_euclidean_random_u(self):
        eps = 1e-07
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-euclidean']

        Y_test1 = pdist(X, u('euclidean'))
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_euclidean_random_float32(self):
        eps = 1e-07
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-euclidean']

        Y_test1 = pdist(X, 'euclidean')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_euclidean_random_nonC(self):
        eps = 1e-07
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-euclidean']
        Y_test2 = pdist(X, 'test_euclidean')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_euclidean_iris_double(self):
        eps = 1e-07
        X = eo['iris']
        Y_right = eo['pdist-euclidean-iris']

        Y_test1 = pdist(X, 'euclidean')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_euclidean_iris_float32(self):
        eps = 1e-06
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-euclidean-iris']

        Y_test1 = pdist(X, 'euclidean')
        if verbose > 2:
            print(np.abs(Y_right - Y_test1).max())
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_euclidean_iris_nonC(self):
        # Test pdist(X, 'test_euclidean') [the non-C implementation] on the
        # Iris data set.
        eps = 1e-07
        X = eo['iris']
        Y_right = eo['pdist-euclidean-iris']
        Y_test2 = pdist(X, 'test_euclidean')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_seuclidean_random(self):
        eps = 1e-05
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-seuclidean']

        Y_test1 = pdist(X, 'seuclidean')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_seuclidean_random_float32(self):
        eps = 1e-05
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-seuclidean']

        Y_test1 = pdist(X, 'seuclidean')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_seuclidean_random_nonC(self):
        # Test pdist(X, 'test_sqeuclidean') [the non-C implementation]
        eps = 1e-05
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-seuclidean']
        Y_test2 = pdist(X, 'test_sqeuclidean')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_seuclidean_iris(self):
        eps = 1e-05
        X = eo['iris']
        Y_right = eo['pdist-seuclidean-iris']

        Y_test1 = pdist(X, 'seuclidean')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_seuclidean_iris_float32(self):
        # Tests pdist(X, 'seuclidean') on the Iris data set (float32).
        eps = 1e-05
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-seuclidean-iris']

        Y_test1 = pdist(X, 'seuclidean')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_seuclidean_iris_nonC(self):
        # Test pdist(X, 'test_seuclidean') [the non-C implementation] on the
        # Iris data set.
        eps = 1e-05
        X = eo['iris']
        Y_right = eo['pdist-seuclidean-iris']
        Y_test2 = pdist(X, 'test_sqeuclidean')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_cosine_random(self):
        eps = 1e-08
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-cosine']
        Y_test1 = pdist(X, 'cosine')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_cosine_random_float32(self):
        eps = 1e-08
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-cosine']

        Y_test1 = pdist(X, 'cosine')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_cosine_random_nonC(self):
        # Test pdist(X, 'test_cosine') [the non-C implementation]
        eps = 1e-08
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-cosine']
        Y_test2 = pdist(X, 'test_cosine')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_cosine_iris(self):
        eps = 1e-08
        X = eo['iris']
        Y_right = eo['pdist-cosine-iris']

        Y_test1 = pdist(X, 'cosine')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_cosine_iris_float32(self):
        eps = 1e-07
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-cosine-iris']

        Y_test1 = pdist(X, 'cosine')
        if verbose > 2:
            print(np.abs(Y_test1 - Y_right).max())
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_cosine_iris_nonC(self):
        eps = 1e-08
        X = eo['iris']
        Y_right = eo['pdist-cosine-iris']
        Y_test2 = pdist(X, 'test_cosine')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_cityblock_random(self):
        eps = 1e-06
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-cityblock']
        Y_test1 = pdist(X, 'cityblock')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_cityblock_random_float32(self):
        eps = 1e-06
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-cityblock']
        Y_test1 = pdist(X, 'cityblock')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_cityblock_random_nonC(self):
        eps = 1e-06
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-cityblock']
        Y_test2 = pdist(X, 'test_cityblock')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_cityblock_iris(self):
        eps = 1e-14
        X = eo['iris']
        Y_right = eo['pdist-cityblock-iris']

        Y_test1 = pdist(X, 'cityblock')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_cityblock_iris_float32(self):
        eps = 1e-06
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-cityblock-iris']

        Y_test1 = pdist(X, 'cityblock')
        if verbose > 2:
            print("cityblock-iris-float32", np.abs(Y_test1 - Y_right).max())
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_cityblock_iris_nonC(self):
        # Test pdist(X, 'test_cityblock') [the non-C implementation] on the
        # Iris data set.
        eps = 1e-14
        X = eo['iris']
        Y_right = eo['pdist-cityblock-iris']
        Y_test2 = pdist(X, 'test_cityblock')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_correlation_random(self):
        eps = 1e-07
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-correlation']

        Y_test1 = pdist(X, 'correlation')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_correlation_random_float32(self):
        eps = 1e-07
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-correlation']

        Y_test1 = pdist(X, 'correlation')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_correlation_random_nonC(self):
        eps = 1e-07
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-correlation']
        Y_test2 = pdist(X, 'test_correlation')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_correlation_iris(self):
        eps = 1e-08
        X = eo['iris']
        Y_right = eo['pdist-correlation-iris']

        Y_test1 = pdist(X, 'correlation')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_correlation_iris_float32(self):
        eps = 1e-07
        X = eo['iris']
        Y_right = np.float32(eo['pdist-correlation-iris'])

        Y_test1 = pdist(X, 'correlation')
        if verbose > 2:
            print("correlation-iris", np.abs(Y_test1 - Y_right).max())
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_correlation_iris_nonC(self):
        eps = 1e-08
        X = eo['iris']
        Y_right = eo['pdist-correlation-iris']
        Y_test2 = pdist(X, 'test_correlation')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_minkowski_random(self):
        eps = 1e-05
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-minkowski-3.2']

        Y_test1 = pdist(X, 'minkowski', 3.2)
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_minkowski_random_float32(self):
        eps = 1e-05
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-minkowski-3.2']

        Y_test1 = pdist(X, 'minkowski', 3.2)
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_minkowski_random_nonC(self):
        eps = 1e-05
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-minkowski-3.2']
        Y_test2 = pdist(X, 'test_minkowski', 3.2)
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_minkowski_3_2_iris(self):
        eps = 1e-07
        X = eo['iris']
        Y_right = eo['pdist-minkowski-3.2-iris']
        Y_test1 = pdist(X, 'minkowski', 3.2)
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_minkowski_3_2_iris_float32(self):
        eps = 1e-06
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-minkowski-3.2-iris']
        Y_test1 = pdist(X, 'minkowski', 3.2)
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_minkowski_3_2_iris_nonC(self):
        eps = 1e-07
        X = eo['iris']
        Y_right = eo['pdist-minkowski-3.2-iris']
        Y_test2 = pdist(X, 'test_minkowski', 3.2)
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_minkowski_5_8_iris(self):
        eps = 1e-07
        X = eo['iris']
        Y_right = eo['pdist-minkowski-5.8-iris']
        Y_test1 = pdist(X, 'minkowski', 5.8)
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_minkowski_5_8_iris_float32(self):
        eps = 1e-06
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-minkowski-5.8-iris']

        Y_test1 = pdist(X, 'minkowski', 5.8)
        if verbose > 2:
            print("minkowski-iris-5.8", np.abs(Y_test1 - Y_right).max())
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_minkowski_5_8_iris_nonC(self):
        eps = 1e-07
        X = eo['iris']
        Y_right = eo['pdist-minkowski-5.8-iris']
        Y_test2 = pdist(X, 'test_minkowski', 5.8)
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_wminkowski(self):
        x = np.array([[0.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [1.0, 1.0, 1.0]])

        p2_expected = [1.0, 1.0, np.sqrt(3),
                       np.sqrt(2), np.sqrt(2),
                       np.sqrt(2)]
        p1_expected = [0.5, 1.0, 3.5,
                       1.5, 3.0,
                       2.5]
        dist = pdist(x, metric=wminkowski, w=[1.0, 1.0, 1.0])
        assert_allclose(dist, p2_expected, rtol=1e-14)

        dist = pdist(x, metric=wminkowski, w=[0.5, 1.0, 2.0], p=1)
        assert_allclose(dist, p1_expected, rtol=1e-14)

        dist = pdist(x, metric='wminkowski', w=[1.0, 1.0, 1.0])
        assert_allclose(dist, p2_expected, rtol=1e-14)

        dist = pdist(x, metric='wminkowski', w=[0.5, 1.0, 2.0], p=1)
        assert_allclose(dist, p1_expected, rtol=1e-14)

    def test_pdist_hamming_random(self):
        eps = 1e-07
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-hamming']

        Y_test1 = pdist(X, 'hamming')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_hamming_random_float32(self):
        eps = 1e-07
        X = np.float32(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']

        Y_test1 = pdist(X, 'hamming')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_hamming_random_nonC(self):
        eps = 1e-07
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-hamming']
        Y_test2 = pdist(X, 'test_hamming')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_dhamming_random(self):
        eps = 1e-07
        X = np.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']
        Y_test1 = pdist(X, 'hamming')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_dhamming_random_float32(self):
        eps = 1e-07
        X = np.float32(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']
        Y_test1 = pdist(X, 'hamming')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_dhamming_random_nonC(self):
        eps = 1e-07
        X = np.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']
        Y_test2 = pdist(X, 'test_hamming')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_jaccard_random(self):
        eps = 1e-08
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-jaccard']

        Y_test1 = pdist(X, 'jaccard')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_jaccard_random_float32(self):
        eps = 1e-08
        X = np.float32(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']

        Y_test1 = pdist(X, 'jaccard')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_jaccard_random_nonC(self):
        eps = 1e-08
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-jaccard']
        Y_test2 = pdist(X, 'test_jaccard')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_djaccard_random(self):
        eps = 1e-08
        X = np.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']

        Y_test1 = pdist(X, 'jaccard')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_djaccard_random_float32(self):
        eps = 1e-08
        X = np.float32(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']

        Y_test1 = pdist(X, 'jaccard')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_djaccard_random_nonC(self):
        eps = 1e-08
        X = np.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']
        Y_test2 = pdist(X, 'test_jaccard')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_chebychev_random(self):
        eps = 1e-08
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-chebychev']

        Y_test1 = pdist(X, 'chebychev')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_chebychev_random_float32(self):
        eps = 1e-07
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-chebychev']

        Y_test1 = pdist(X, 'chebychev')
        if verbose > 2:
            print("chebychev", np.abs(Y_test1 - Y_right).max())
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_chebychev_random_nonC(self):
        eps = 1e-08
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-chebychev']
        Y_test2 = pdist(X, 'test_chebychev')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_chebychev_iris(self):
        eps = 1e-15
        X = eo['iris']
        Y_right = eo['pdist-chebychev-iris']
        Y_test1 = pdist(X, 'chebychev')
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_chebychev_iris_float32(self):
        eps = 1e-06
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-chebychev-iris']
        Y_test1 = pdist(X, 'chebychev')
        if verbose > 2:
            print("chebychev-iris", np.abs(Y_test1 - Y_right).max())
        self.assertTrue(within_tol(Y_test1, Y_right, eps))

    def test_pdist_chebychev_iris_nonC(self):
        eps = 1e-15
        X = eo['iris']
        Y_right = eo['pdist-chebychev-iris']
        Y_test2 = pdist(X, 'test_chebychev')
        self.assertTrue(within_tol(Y_test2, Y_right, eps))

    def test_pdist_matching_mtica1(self):
        # Test matching(*,*) with mtica example #1 (nums).
        m = matching(np.array([1, 0, 1, 1, 0]),
                     np.array([1, 1, 0, 1, 1]))
        m2 = matching(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                      np.array([1, 1, 0, 1, 1], dtype=np.bool))
        self.assertTrue(np.abs(m - 0.6) <= 1e-10)
        self.assertTrue(np.abs(m2 - 0.6) <= 1e-10)

    def test_pdist_matching_mtica2(self):
        # Test matching(*,*) with mtica example #2.
        m = matching(np.array([1, 0, 1]),
                     np.array([1, 1, 0]))
        m2 = matching(np.array([1, 0, 1], dtype=np.bool),
                      np.array([1, 1, 0], dtype=np.bool))
        self.assertTrue(np.abs(m - (2.0/3.0)) <= 1e-10)
        self.assertTrue(np.abs(m2 - (2.0/3.0)) <= 1e-10)

    def test_pdist_matching_match(self):
        # Test pdist(X, 'matching') to see if the two implementations match on
        # random boolean input data.
        D = eo['random-bool-data']
        B = np.bool_(D)
        if verbose > 2:
            print(B.shape, B.dtype)
        eps = 1e-10
        y1 = pdist(B, "matching")
        y2 = pdist(B, "test_matching")
        y3 = pdist(D, "test_matching")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y1-y3).max())
        self.assertTrue(within_tol(y1, y2, eps))
        self.assertTrue(within_tol(y2, y3, eps))

    def test_pdist_jaccard_mtica1(self):
        m = jaccard(np.array([1, 0, 1, 1, 0]),
                    np.array([1, 1, 0, 1, 1]))
        m2 = jaccard(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                     np.array([1, 1, 0, 1, 1], dtype=np.bool))
        self.assertTrue(np.abs(m - 0.6) <= 1e-10)
        self.assertTrue(np.abs(m2 - 0.6) <= 1e-10)

    def test_pdist_jaccard_mtica2(self):
        m = jaccard(np.array([1, 0, 1]),
                    np.array([1, 1, 0]))
        m2 = jaccard(np.array([1, 0, 1], dtype=np.bool),
                     np.array([1, 1, 0], dtype=np.bool))
        self.assertTrue(np.abs(m - (2.0/3.0)) <= 1e-10)
        self.assertTrue(np.abs(m2 - (2.0/3.0)) <= 1e-10)

    def test_pdist_jaccard_match(self):
        # Test pdist(X, 'jaccard') to see if the two implementations match on
        # random double input data.
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = pdist(D, "jaccard")
        y2 = pdist(D, "test_jaccard")
        y3 = pdist(np.bool_(D), "test_jaccard")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        self.assertTrue(within_tol(y1, y2, eps))
        self.assertTrue(within_tol(y2, y3, eps))

    def test_pdist_yule_mtica1(self):
        m = yule(np.array([1, 0, 1, 1, 0]),
                 np.array([1, 1, 0, 1, 1]))
        m2 = yule(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                  np.array([1, 1, 0, 1, 1], dtype=np.bool))
        if verbose > 2:
            print(m)
        self.assertTrue(np.abs(m - 2.0) <= 1e-10)
        self.assertTrue(np.abs(m2 - 2.0) <= 1e-10)

    def test_pdist_yule_mtica2(self):
        m = yule(np.array([1, 0, 1]),
                 np.array([1, 1, 0]))
        m2 = yule(np.array([1, 0, 1], dtype=np.bool),
                  np.array([1, 1, 0], dtype=np.bool))
        if verbose > 2:
            print(m)
        self.assertTrue(np.abs(m - 2.0) <= 1e-10)
        self.assertTrue(np.abs(m2 - 2.0) <= 1e-10)

    def test_pdist_yule_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = pdist(D, "yule")
        y2 = pdist(D, "test_yule")
        y3 = pdist(np.bool_(D), "test_yule")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        self.assertTrue(within_tol(y1, y2, eps))
        self.assertTrue(within_tol(y2, y3, eps))

    def test_pdist_dice_mtica1(self):
        m = dice(np.array([1, 0, 1, 1, 0]),
                 np.array([1, 1, 0, 1, 1]))
        m2 = dice(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                  np.array([1, 1, 0, 1, 1], dtype=np.bool))
        if verbose > 2:
            print(m)
        self.assertTrue(np.abs(m - (3.0/7.0)) <= 1e-10)
        self.assertTrue(np.abs(m2 - (3.0/7.0)) <= 1e-10)

    def test_pdist_dice_mtica2(self):
        m = dice(np.array([1, 0, 1]),
                 np.array([1, 1, 0]))
        m2 = dice(np.array([1, 0, 1], dtype=np.bool),
                  np.array([1, 1, 0], dtype=np.bool))
        if verbose > 2:
            print(m)
        self.assertTrue(np.abs(m - 0.5) <= 1e-10)
        self.assertTrue(np.abs(m2 - 0.5) <= 1e-10)

    def test_pdist_dice_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = pdist(D, "dice")
        y2 = pdist(D, "test_dice")
        y3 = pdist(D, "test_dice")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        self.assertTrue(within_tol(y1, y2, eps))
        self.assertTrue(within_tol(y2, y3, eps))

    def test_pdist_sokalsneath_mtica1(self):
        m = sokalsneath(np.array([1, 0, 1, 1, 0]),
                        np.array([1, 1, 0, 1, 1]))
        m2 = sokalsneath(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                         np.array([1, 1, 0, 1, 1], dtype=np.bool))
        if verbose > 2:
            print(m)
        self.assertTrue(np.abs(m - (3.0/4.0)) <= 1e-10)
        self.assertTrue(np.abs(m2 - (3.0/4.0)) <= 1e-10)

    def test_pdist_sokalsneath_mtica2(self):
        m = sokalsneath(np.array([1, 0, 1]),
                        np.array([1, 1, 0]))
        m2 = sokalsneath(np.array([1, 0, 1], dtype=np.bool),
                         np.array([1, 1, 0], dtype=np.bool))
        if verbose > 2:
            print(m)
        self.assertTrue(np.abs(m - (4.0/5.0)) <= 1e-10)
        self.assertTrue(np.abs(m2 - (4.0/5.0)) <= 1e-10)

    def test_pdist_sokalsneath_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = pdist(D, "sokalsneath")
        y2 = pdist(D, "test_sokalsneath")
        y3 = pdist(np.bool_(D), "test_sokalsneath")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        self.assertTrue(within_tol(y1, y2, eps))
        self.assertTrue(within_tol(y2, y3, eps))

    def test_pdist_rogerstanimoto_mtica1(self):
        m = rogerstanimoto(np.array([1, 0, 1, 1, 0]),
                           np.array([1, 1, 0, 1, 1]))
        m2 = rogerstanimoto(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                            np.array([1, 1, 0, 1, 1], dtype=np.bool))
        if verbose > 2:
            print(m)
        self.assertTrue(np.abs(m - (3.0/4.0)) <= 1e-10)
        self.assertTrue(np.abs(m2 - (3.0/4.0)) <= 1e-10)

    def test_pdist_rogerstanimoto_mtica2(self):
        m = rogerstanimoto(np.array([1, 0, 1]),
                           np.array([1, 1, 0]))
        m2 = rogerstanimoto(np.array([1, 0, 1], dtype=np.bool),
                            np.array([1, 1, 0], dtype=np.bool))
        if verbose > 2:
            print(m)
        self.assertTrue(np.abs(m - (4.0/5.0)) <= 1e-10)
        self.assertTrue(np.abs(m2 - (4.0/5.0)) <= 1e-10)

    def test_pdist_rogerstanimoto_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = pdist(D, "rogerstanimoto")
        y2 = pdist(D, "test_rogerstanimoto")
        y3 = pdist(np.bool_(D), "test_rogerstanimoto")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        self.assertTrue(within_tol(y1, y2, eps))
        self.assertTrue(within_tol(y2, y3, eps))

    def test_pdist_russellrao_mtica1(self):
        m = russellrao(np.array([1, 0, 1, 1, 0]),
                       np.array([1, 1, 0, 1, 1]))
        m2 = russellrao(np.array([1, 0, 1, 1, 0], dtype=np.bool),
                        np.array([1, 1, 0, 1, 1], dtype=np.bool))
        if verbose > 2:
            print(m)
        self.assertTrue(np.abs(m - (3.0/5.0)) <= 1e-10)
        self.assertTrue(np.abs(m2 - (3.0/5.0)) <= 1e-10)

    def test_pdist_russellrao_mtica2(self):
        m = russellrao(np.array([1, 0, 1]),
                       np.array([1, 1, 0]))
        m2 = russellrao(np.array([1, 0, 1], dtype=np.bool),
                        np.array([1, 1, 0], dtype=np.bool))
        if verbose > 2:
            print(m)
        self.assertTrue(np.abs(m - (2.0/3.0)) <= 1e-10)
        self.assertTrue(np.abs(m2 - (2.0/3.0)) <= 1e-10)

    def test_pdist_russellrao_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = pdist(D, "russellrao")
        y2 = pdist(D, "test_russellrao")
        y3 = pdist(np.bool_(D), "test_russellrao")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        self.assertTrue(within_tol(y1, y2, eps))
        self.assertTrue(within_tol(y2, y3, eps))

    def test_pdist_sokalmichener_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = pdist(D, "sokalmichener")
        y2 = pdist(D, "test_sokalmichener")
        y3 = pdist(np.bool_(D), "test_sokalmichener")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        self.assertTrue(within_tol(y1, y2, eps))
        self.assertTrue(within_tol(y2, y3, eps))

    def test_pdist_kulsinski_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = pdist(D, "kulsinski")
        y2 = pdist(D, "test_kulsinski")
        y3 = pdist(np.bool_(D), "test_kulsinski")
        if verbose > 2:
            print(np.abs(y1-y2).max())
        self.assertTrue(within_tol(y1, y2, eps))
        self.assertTrue(within_tol(y2, y3, eps))

    def test_pdist_canberra_match(self):
        D = eo['iris']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = pdist(D, "canberra")
        y2 = pdist(D, "test_canberra")
        if verbose > 2:
            print(np.abs(y1-y2).max())
        self.assertTrue(within_tol(y1, y2, eps))

    def test_pdist_canberra_ticket_711(self):
        # Test pdist(X, 'canberra') to see if Canberra gives the right result
        # as reported on gh-1238.
        eps = 1e-8
        pdist_y = pdist(([3.3], [3.4]), "canberra")
        right_y = 0.01492537
        if verbose > 2:
            print(np.abs(pdist_y-right_y).max())
        self.assertTrue(within_tol(pdist_y, right_y, eps))


def within_tol(a, b, tol):
    return np.abs(a - b).max() < tol


def _assert_within_tol(a, b, atol, verbose_=False):
    if verbose_:
        print(np.abs(a-b).max())
    assert_allclose(a, b, rtol=0, atol=atol)


class TestSomeDistanceFunctions(TestCase):

    def setUp(self):
        # 1D arrays
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([1.0, 1.0, 5.0])
        # 3x1 arrays
        x31 = x[:,np.newaxis]
        y31 = y[:,np.newaxis]
        # 1x3 arrays
        x13 = x31.T
        y13 = y31.T

        self.cases = [(x,y), (x31, y31), (x13, y13)]

    def test_minkowski(self):
        for x, y in self.cases:
            dist1 = minkowski(x, y, p=1)
            assert_almost_equal(dist1, 3.0)
            dist1p5 = minkowski(x, y, p=1.5)
            assert_almost_equal(dist1p5, (1.0+2.0**1.5)**(2./3))
            dist2 = minkowski(x, y, p=2)
            assert_almost_equal(dist2, np.sqrt(5))

    def test_wminkowski(self):
        w = np.array([1.0, 2.0, 0.5])
        for x, y in self.cases:
            dist1 = wminkowski(x, y, p=1, w=w)
            assert_almost_equal(dist1, 3.0)
            dist1p5 = wminkowski(x, y, p=1.5, w=w)
            assert_almost_equal(dist1p5, (2.0**1.5+1.0)**(2./3))
            dist2 = wminkowski(x, y, p=2, w=w)
            assert_almost_equal(dist2, np.sqrt(5))

    def test_euclidean(self):
        for x, y in self.cases:
            dist = euclidean(x, y)
            assert_almost_equal(dist, np.sqrt(5))

    def test_sqeuclidean(self):
        for x, y in self.cases:
            dist = sqeuclidean(x, y)
            assert_almost_equal(dist, 5.0)

    def test_cosine(self):
        for x, y in self.cases:
            dist = cosine(x, y)
            assert_almost_equal(dist, 1.0 - 18.0/(np.sqrt(14)*np.sqrt(27)))

    def test_correlation(self):
        xm = np.array([-1.0, 0, 1.0])
        ym = np.array([-4.0/3, -4.0/3, 5.0-7.0/3])
        for x, y in self.cases:
            dist = correlation(x, y)
            assert_almost_equal(dist, 1.0 - np.dot(xm, ym)/(norm(xm)*norm(ym)))

    def test_mahalanobis(self):
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([1.0, 1.0, 5.0])
        vi = np.array([[2.0, 1.0, 0.0],[1.0, 2.0, 1.0], [0.0, 1.0, 2.0]])
        for x, y in self.cases:
            dist = mahalanobis(x, y, vi)
            assert_almost_equal(dist, np.sqrt(6.0))


class TestSquareForm(TestCase):

    def test_squareform_empty_matrix(self):
        A = np.zeros((0,0))
        rA = squareform(np.array(A, dtype='double'))
        self.assertTrue(rA.shape == (0,))

    def test_squareform_empty_vector(self):
        v = np.zeros((0,))
        rv = squareform(np.array(v, dtype='double'))
        self.assertTrue(rv.shape == (1,1))
        self.assertTrue(rv[0, 0] == 0)

    def test_squareform_1by1_matrix(self):
        A = np.zeros((1,1))
        rA = squareform(np.array(A, dtype='double'))
        self.assertTrue(rA.shape == (0,))

    def test_squareform_one_vector(self):
        v = np.ones((1,)) * 8.3
        rv = squareform(np.array(v, dtype='double'))
        self.assertTrue(rv.shape == (2,2))
        self.assertTrue(rv[0,1] == 8.3)
        self.assertTrue(rv[1,0] == 8.3)

    def test_squareform_one_binary_vector(self):
        # Tests squareform on a 1x1 binary matrix; conversion to double was
        # causing problems (see pull request 73).
        v = np.ones((1,), dtype=np.bool)
        rv = squareform(v)
        self.assertTrue(rv.shape == (2,2))
        self.assertTrue(rv[0,1])

    def test_squareform_2by2_matrix(self):
        A = np.zeros((2,2))
        A[0,1] = 0.8
        A[1,0] = 0.8
        rA = squareform(np.array(A, dtype='double'))
        self.assertTrue(rA.shape == (1,))
        self.assertTrue(rA[0] == 0.8)

    def test_squareform_multi_matrix(self):
        for n in xrange(2, 5):
            yield self.check_squareform_multi_matrix(n)

    def check_squareform_multi_matrix(self, n):
        X = np.random.rand(n, 4)
        Y = pdist(X)
        self.assertTrue(len(Y.shape) == 1)
        A = squareform(Y)
        Yr = squareform(A)
        s = A.shape
        k = 0
        if verbose >= 3:
            print(A.shape, Y.shape, Yr.shape)
        self.assertTrue(len(s) == 2)
        self.assertTrue(len(Yr.shape) == 1)
        self.assertTrue(s[0] == s[1])
        for i in xrange(0, s[0]):
            for j in xrange(i+1, s[1]):
                if i != j:
                    self.assertTrue(A[i, j] == Y[k])
                    k += 1
                else:
                    self.assertTrue(A[i, j] == 0)


class TestNumObsY(TestCase):

    def test_num_obs_y_multi_matrix(self):
        for n in xrange(2, 10):
            X = np.random.rand(n, 4)
            Y = pdist(X)
            self.assertTrue(num_obs_y(Y) == n)

    def test_num_obs_y_1(self):
        # Tests num_obs_y(y) on a condensed distance matrix over 1
        # observations. Expecting exception.
        assert_raises(ValueError, self.check_y, 1)

    def test_num_obs_y_2(self):
        # Tests num_obs_y(y) on a condensed distance matrix over 2
        # observations.
        self.assertTrue(self.check_y(2))

    def test_num_obs_y_3(self):
        self.assertTrue(self.check_y(3))

    def test_num_obs_y_4(self):
        self.assertTrue(self.check_y(4))

    def test_num_obs_y_5_10(self):
        for i in xrange(5, 16):
            self.minit(i)

    def test_num_obs_y_2_100(self):
        # Tests num_obs_y(y) on 100 improper condensed distance matrices.
        # Expecting exception.
        a = set([])
        for n in xrange(2, 16):
            a.add(n*(n-1)/2)
        for i in xrange(5, 105):
            if i not in a:
                assert_raises(ValueError, self.bad_y, i)

    def minit(self, n):
        self.assertTrue(self.check_y(n))

    def bad_y(self, n):
        y = np.random.rand(n)
        return num_obs_y(y)

    def check_y(self, n):
        return num_obs_y(self.make_y(n)) == n

    def make_y(self, n):
        return np.random.rand((n * (n - 1)) // 2)


class TestNumObsDM(TestCase):

    def test_num_obs_dm_multi_matrix(self):
        for n in xrange(1, 10):
            X = np.random.rand(n, 4)
            Y = pdist(X)
            A = squareform(Y)
            if verbose >= 3:
                print(A.shape, Y.shape)
            self.assertTrue(num_obs_dm(A) == n)

    def test_num_obs_dm_0(self):
        # Tests num_obs_dm(D) on a 0x0 distance matrix. Expecting exception.
        self.assertTrue(self.check_D(0))

    def test_num_obs_dm_1(self):
        # Tests num_obs_dm(D) on a 1x1 distance matrix.
        self.assertTrue(self.check_D(1))

    def test_num_obs_dm_2(self):
        self.assertTrue(self.check_D(2))

    def test_num_obs_dm_3(self):
        self.assertTrue(self.check_D(2))

    def test_num_obs_dm_4(self):
        self.assertTrue(self.check_D(4))

    def check_D(self, n):
        return num_obs_dm(self.make_D(n)) == n

    def make_D(self, n):
        return np.random.rand(n, n)


def is_valid_dm_throw(D):
    return is_valid_dm(D, throw=True)


class TestIsValidDM(TestCase):

    def test_is_valid_dm_int16_array_E(self):
        # Tests is_valid_dm(*) on an int16 array. Exception expected.
        D = np.zeros((5, 5), dtype='i')
        assert_raises(TypeError, is_valid_dm_throw, (D))

    def test_is_valid_dm_int16_array_F(self):
        D = np.zeros((5, 5), dtype='i')
        self.assertTrue(is_valid_dm(D) == False)

    def test_is_valid_dm_improper_shape_1D_E(self):
        D = np.zeros((5,), dtype=np.double)
        assert_raises(ValueError, is_valid_dm_throw, (D))

    def test_is_valid_dm_improper_shape_1D_F(self):
        D = np.zeros((5,), dtype=np.double)
        self.assertTrue(is_valid_dm(D) == False)

    def test_is_valid_dm_improper_shape_3D_E(self):
        D = np.zeros((3,3,3), dtype=np.double)
        assert_raises(ValueError, is_valid_dm_throw, (D))

    def test_is_valid_dm_improper_shape_3D_F(self):
        D = np.zeros((3,3,3), dtype=np.double)
        self.assertTrue(is_valid_dm(D) == False)

    def test_is_valid_dm_nonzero_diagonal_E(self):
        y = np.random.rand(10)
        D = squareform(y)
        for i in xrange(0, 5):
            D[i, i] = 2.0
        assert_raises(ValueError, is_valid_dm_throw, (D))

    def test_is_valid_dm_nonzero_diagonal_F(self):
        y = np.random.rand(10)
        D = squareform(y)
        for i in xrange(0, 5):
            D[i, i] = 2.0
        self.assertTrue(is_valid_dm(D) == False)

    def test_is_valid_dm_asymmetric_E(self):
        y = np.random.rand(10)
        D = squareform(y)
        D[1,3] = D[3,1] + 1
        assert_raises(ValueError, is_valid_dm_throw, (D))

    def test_is_valid_dm_asymmetric_F(self):
        y = np.random.rand(10)
        D = squareform(y)
        D[1,3] = D[3,1] + 1
        self.assertTrue(is_valid_dm(D) == False)

    def test_is_valid_dm_correct_1_by_1(self):
        D = np.zeros((1,1), dtype=np.double)
        self.assertTrue(is_valid_dm(D) == True)

    def test_is_valid_dm_correct_2_by_2(self):
        y = np.random.rand(1)
        D = squareform(y)
        self.assertTrue(is_valid_dm(D) == True)

    def test_is_valid_dm_correct_3_by_3(self):
        y = np.random.rand(3)
        D = squareform(y)
        self.assertTrue(is_valid_dm(D) == True)

    def test_is_valid_dm_correct_4_by_4(self):
        y = np.random.rand(6)
        D = squareform(y)
        self.assertTrue(is_valid_dm(D) == True)

    def test_is_valid_dm_correct_5_by_5(self):
        y = np.random.rand(10)
        D = squareform(y)
        self.assertTrue(is_valid_dm(D) == True)


def is_valid_y_throw(y):
    return is_valid_y(y, throw=True)


class TestIsValidY(TestCase):
    # If test case name ends on "_E" then an exception is expected for the
    # given input, if it ends in "_F" then False is expected for the is_valid_y
    # check.  Otherwise the input is expected to be valid.

    def test_is_valid_y_int16_array_E(self):
        y = np.zeros((10,), dtype='i')
        assert_raises(TypeError, is_valid_y_throw, (y))

    def test_is_valid_y_int16_array_F(self):
        y = np.zeros((10,), dtype='i')
        self.assertTrue(is_valid_y(y) == False)

    def test_is_valid_y_improper_shape_2D_E(self):
        y = np.zeros((3,3,), dtype=np.double)
        assert_raises(ValueError, is_valid_y_throw, (y))

    def test_is_valid_y_improper_shape_2D_F(self):
        y = np.zeros((3,3,), dtype=np.double)
        self.assertTrue(is_valid_y(y) == False)

    def test_is_valid_y_improper_shape_3D_E(self):
        y = np.zeros((3,3,3), dtype=np.double)
        assert_raises(ValueError, is_valid_y_throw, (y))

    def test_is_valid_y_improper_shape_3D_F(self):
        y = np.zeros((3,3,3), dtype=np.double)
        self.assertTrue(is_valid_y(y) == False)

    def test_is_valid_y_correct_2_by_2(self):
        y = self.correct_n_by_n(2)
        self.assertTrue(is_valid_y(y) == True)

    def test_is_valid_y_correct_3_by_3(self):
        y = self.correct_n_by_n(3)
        self.assertTrue(is_valid_y(y) == True)

    def test_is_valid_y_correct_4_by_4(self):
        y = self.correct_n_by_n(4)
        self.assertTrue(is_valid_y(y) == True)

    def test_is_valid_y_correct_5_by_5(self):
        y = self.correct_n_by_n(5)
        self.assertTrue(is_valid_y(y) == True)

    def test_is_valid_y_2_100(self):
        a = set([])
        for n in xrange(2, 16):
            a.add(n*(n-1)/2)
        for i in xrange(5, 105):
            if i not in a:
                assert_raises(ValueError, self.bad_y, i)

    def bad_y(self, n):
        y = np.random.rand(n)
        return is_valid_y(y, throw=True)

    def correct_n_by_n(self, n):
        y = np.random.rand((n * (n - 1)) // 2)
        return y


def test_bad_p():
    # Raise ValueError if p < 1.
    p = 0.5
    assert_raises(ValueError, minkowski, [1, 2], [3, 4], p)
    assert_raises(ValueError, wminkowski, [1, 2], [3, 4], p, [1, 1])


def test_sokalsneath_all_false():
    # Regression test for ticket #876
    assert_raises(ValueError, sokalsneath, [False, False, False], [False, False, False])


def test_canberra():
    # Regression test for ticket #1430.
    assert_equal(canberra([1,2,3], [2,4,6]), 1)
    assert_equal(canberra([1,1,0,0], [1,0,1,0]), 2)


def test_braycurtis():
    # Regression test for ticket #1430.
    assert_almost_equal(braycurtis([1,2,3], [2,4,6]), 1./3, decimal=15)
    assert_almost_equal(braycurtis([1,1,0,0], [1,0,1,0]), 0.5, decimal=15)


def test_euclideans():
    # Regression test for ticket #1328.
    x1 = np.array([1, 1, 1])
    x2 = np.array([0, 0, 0])

    # Basic test of the calculation.
    assert_almost_equal(sqeuclidean(x1, x2), 3.0, decimal=14)
    assert_almost_equal(euclidean(x1, x2), np.sqrt(3), decimal=14)

    # Check flattening for (1, N) or (N, 1) inputs
    assert_almost_equal(euclidean(x1[np.newaxis, :], x2[np.newaxis, :]),
                        np.sqrt(3), decimal=14)
    assert_almost_equal(sqeuclidean(x1[np.newaxis, :], x2[np.newaxis, :]),
                        3.0, decimal=14)
    assert_almost_equal(sqeuclidean(x1[:, np.newaxis], x2[:, np.newaxis]),
                        3.0, decimal=14)

    # Distance metrics only defined for vectors (= 1-D)
    x = np.arange(4).reshape(2, 2)
    assert_raises(ValueError, euclidean, x, x)
    assert_raises(ValueError, sqeuclidean, x, x)

    # Another check, with random data.
    rs = np.random.RandomState(1234567890)
    x = rs.rand(10)
    y = rs.rand(10)
    d1 = euclidean(x, y)
    d2 = sqeuclidean(x, y)
    assert_almost_equal(d1**2, d2, decimal=14)


def test_hamming_unequal_length():
    # Regression test for gh-4290.
    x = [0, 0, 1]
    y = [1, 0, 1, 0]
    # Used to give an AttributeError from ndarray.mean called on bool
    assert_raises(ValueError, hamming, x, y)


def test_hamming_string_array():
    # https://github.com/scikit-learn/scikit-learn/issues/4014
    a = np.array(['eggs', 'spam', 'spam', 'eggs', 'spam', 'spam', 'spam',
                  'spam', 'spam', 'spam', 'spam', 'eggs', 'eggs', 'spam',
                  'eggs', 'eggs', 'eggs', 'eggs', 'eggs', 'spam'],
                  dtype='|S4')
    b = np.array(['eggs', 'spam', 'spam', 'eggs', 'eggs', 'spam', 'spam',
                  'spam', 'spam', 'eggs', 'spam', 'eggs', 'spam', 'eggs',
                  'spam', 'spam', 'eggs', 'spam', 'spam', 'eggs'],
                  dtype='|S4')
    desired = 0.45
    assert_allclose(hamming(a, b), desired)


def test_sqeuclidean_dtypes():
    # Assert that sqeuclidean returns the right types of values.
    # Integer types should be converted to floating for stability.
    # Floating point types should be the same as the input.
    x = [1, 2, 3]
    y = [4, 5, 6]

    for dtype in [np.int8, np.int16, np.int32, np.int64]:
        d = sqeuclidean(np.asarray(x, dtype=dtype), np.asarray(y, dtype=dtype))
        assert_(np.issubdtype(d.dtype, np.floating))

    for dtype in [np.uint8, np.uint16, np.uint32, np.uint64]:
        d1 = sqeuclidean([0], np.asarray([-1], dtype=dtype))
        d2 = sqeuclidean(np.asarray([-1], dtype=dtype), [0])

        assert_equal(d1, d2)
        assert_equal(d1, np.float64(np.iinfo(dtype).max) ** 2)

    dtypes = [np.float32, np.float64, np.complex64, np.complex128]
    for dtype in ['float16', 'float128']:
        # These aren't present in older numpy versions; float128 may also not
        # be present on all platforms.
        if hasattr(np, dtype):
            dtypes.append(getattr(np, dtype))

    for dtype in dtypes:
        d = sqeuclidean(np.asarray(x, dtype=dtype), np.asarray(y, dtype=dtype))
        assert_equal(d.dtype, dtype)


def test_sokalmichener():
    # Test that sokalmichener has the same result for bool and int inputs.
    p = [True, True, False]
    q = [True, False, True]
    x = [int(b) for b in p]
    y = [int(b) for b in q]
    dist1 = sokalmichener(p, q)
    dist2 = sokalmichener(x, y)
    # These should be exactly the same.
    assert_equal(dist1, dist2)


def test__validate_vector():
    x = [1, 2, 3]
    y = _validate_vector(x)
    assert_array_equal(y, x)

    y = _validate_vector(x, dtype=np.float64)
    assert_array_equal(y, x)
    assert_equal(y.dtype, np.float64)

    x = [1]
    y = _validate_vector(x)
    assert_equal(y.ndim, 1)
    assert_equal(y, x)

    x = 1
    y = _validate_vector(x)
    assert_equal(y.ndim, 1)
    assert_equal(y, [x])

    x = np.arange(5).reshape(1, -1, 1)
    y = _validate_vector(x)
    assert_equal(y.ndim, 1)
    assert_array_equal(y, x[0, :, 0])

    x = [[1, 2], [3, 4]]
    assert_raises(ValueError, _validate_vector, x)


if __name__ == "__main__":
    run_module_suite()
