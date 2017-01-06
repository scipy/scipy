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

import warnings
import os.path
from functools import wraps, partial
from scipy._lib.six import xrange, u
from scipy.stats.stats import _chk_weights

import numpy as np
from numpy.linalg import norm
from numpy.testing import (verbose, TestCase, run_module_suite, assert_,
        assert_raises, assert_array_equal, assert_equal, assert_almost_equal,
        assert_allclose)

from scipy.spatial.distance import (squareform, pdist, cdist,
        jaccard, dice, sokalsneath, rogerstanimoto, russellrao, yule,
        num_obs_y, num_obs_dm, is_valid_dm, is_valid_y, minkowski,
        euclidean, sqeuclidean, cosine, correlation, hamming, mahalanobis,
        canberra, braycurtis, sokalmichener, _validate_vector,
        chebyshev, cityblock, kulsinski, seuclidean)
matching = hamming  # suppresses warning messages in testing

_filenames = [
              "cdist-X1.txt",
              "cdist-X2.txt",
              "iris.txt",
              "pdist-boolean-inp.txt",
              "pdist-chebychev-ml-iris.txt",
              "pdist-chebychev-ml.txt",
              "pdist-cityblock-ml-iris.txt",
              "pdist-cityblock-ml.txt",
              "pdist-correlation-ml-iris.txt",
              "pdist-correlation-ml.txt",
              "pdist-cosine-ml-iris.txt",
              "pdist-cosine-ml.txt",
              "pdist-double-inp.txt",
              "pdist-euclidean-ml-iris.txt",
              "pdist-euclidean-ml.txt",
              "pdist-hamming-ml.txt",
              "pdist-jaccard-ml.txt",
              "pdist-minkowski-3.2-ml-iris.txt",
              "pdist-minkowski-3.2-ml.txt",
              "pdist-minkowski-5.8-ml-iris.txt",
              "pdist-seuclidean-ml-iris.txt",
              "pdist-seuclidean-ml.txt",
              "pdist-spearman-ml.txt",
              "random-bool-data.txt",
              "random-double-data.txt",
              "random-int-data.txt",
              "random-uint-data.txt",
]

_tdist = np.array([[0, 662, 877, 255, 412, 996],
                      [662, 0, 295, 468, 268, 400],
                      [877, 295, 0, 754, 564, 138],
                      [255, 468, 754, 0, 219, 869],
                      [412, 268, 564, 219, 0, 669],
                      [996, 400, 138, 869, 669, 0]], dtype='double')

_ytdist = squareform(_tdist)

_metrics = ['braycurtis', 'canberra', 'chebyshev', 'cityblock',
            'correlation', 'cosine', 'dice', 'euclidean', 'hamming',
            'jaccard', 'kulsinski', 'mahalanobis', 'matching',
            'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean',
            'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']

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
    eo['random-bool-data'] = np.bool_(eo['random-bool-data'])
    eo['random-float32-data'] = np.float32(eo['random-double-data'])
    eo['random-int-data'] = np.int_(eo['random-int-data'])
    eo['random-uint-data'] = np.uint(eo['random-uint-data'])

load_testing_files()


def within_tol(a, b, tol):
    return np.abs(a - b).max() < tol


def _assert_within_tol(a, b, atol=0, rtol=0, verbose_=False):
    if verbose_:
        print(np.abs(a - b).max())
    assert_allclose(a, b, rtol=rtol, atol=atol)


def _rand_split(arrays, weights, axis, split_per, seed=None):
    # inverse operation for stats.collapse_weights
    weights = np.array(weights)  # modified inplace; need a copy
    seeded_rand = np.random.RandomState(seed)

    def mytake(a, ix, axis):
        record = np.asanyarray(np.take(a, ix, axis=axis))
        return record.reshape([a.shape[i] if i != axis else 1
                               for i in range(a.ndim)])

    n_obs = arrays[0].shape[axis]
    assert all(a.shape[axis] == n_obs for a in arrays), "data must be aligned on sample axis"
    for i in range(int(split_per) * n_obs):
        split_ix = seeded_rand.randint(n_obs + i)
        prev_w = weights[split_ix]
        q = seeded_rand.rand()
        weights[split_ix] = q * prev_w
        weights = np.append(weights, (1-q) * prev_w)
        arrays = [np.append(a, mytake(a, split_ix, axis=axis),
                            axis=axis) for a in arrays]
    return arrays, weights

def _rough_check(a, b, compare_assert=partial(assert_allclose, atol=1e-5),
                  key=lambda x: x, w=None):
    check_a = key(a)
    check_b = key(b)
    try:
        if np.array(check_a != check_b).any():  # try strict equality for string types
            compare_assert(check_a, check_b)
    except AttributeError:  # masked array
        compare_assert(check_a, check_b)
    except (TypeError, ValueError):  # nested data structure
        for a_i, b_i in zip(check_a, check_b):
            _rough_check(a_i, b_i, compare_assert=compare_assert)

# diff from test_stats:
#  n_args=2, weight_arg='w', default_axis=None
#  ma_safe = False, nan_safe = False
def _weight_checked(fn, n_args=2, default_axis=None, key=lambda x: x, weight_arg='w',
                    squeeze=True, silent=False,
                    ones_test=True, const_test=True, dup_test=True,
                    split_test=True, dud_test=True, ma_safe=False, ma_very_safe=False, nan_safe=False,
                    split_per=1.0, seed=0, compare_assert=partial(assert_allclose, atol=1e-5)):
    """runs fn on its arguments 2 or 3 ways, checks that the results are the same,
       then returns the same thing it would have returned before"""
    @wraps(fn)
    def wrapped(*args, **kwargs):
        result = fn(*args, **kwargs)

        arrays = args[:n_args]
        rest = args[n_args:]
        weights = kwargs.get(weight_arg, None)
        axis = kwargs.get('axis', default_axis)

        chked = _chk_weights(arrays, weights=weights, axis=axis, force_weights=True, mask_screen=True)
        arrays, weights, axis = chked[:-2], chked[-2], chked[-1]
        if squeeze:
            arrays = [np.atleast_1d(a.squeeze()) for a in arrays]

        try:
            # WEIGHTS CHECK 1: EQUAL WEIGHTED OBESERVATIONS
            args = tuple(arrays) + rest
            if ones_test:
                kwargs[weight_arg] = weights
                _rough_check(result, fn(*args, **kwargs), key=key)
            if const_test:
                kwargs[weight_arg] = weights * 101.0
                _rough_check(result, fn(*args, **kwargs), key=key)
                kwargs[weight_arg] = weights * 0.101
                try:
                    _rough_check(result, fn(*args, **kwargs), key=key)
                except Exception as e:
                    raise type(e)((e, arrays, weights))

            # WEIGHTS CHECK 2: ADDL 0-WEIGHTED OBS
            if dud_test:
                # add randomly resampled rows, weighted at 0
                dud_arrays, dud_weights = _rand_split(arrays, weights, axis, split_per=split_per, seed=seed)
                dud_weights[:weights.size] = weights  # not exactly 1 because of masked arrays
                dud_weights[weights.size:] = 0
                dud_args = tuple(dud_arrays) + rest
                kwargs[weight_arg] = dud_weights
                _rough_check(result, fn(*dud_args, **kwargs), key=key)
                # increase the value of those 0-weighted rows
                for a in dud_arrays:
                    indexer = [slice(None)] * a.ndim
                    indexer[axis] = slice(weights.size, None)
                    a[indexer] = a[indexer] * 101 
                dud_args = tuple(dud_arrays) + rest
                _rough_check(result, fn(*dud_args, **kwargs), key=key)
                # set those 0-weighted rows to NaNs
                for a in dud_arrays:
                    indexer = [slice(None)] * a.ndim
                    indexer[axis] = slice(weights.size, None)
                    a[indexer] = a[indexer] * np.nan
                if kwargs.get("nan_policy", None) == "omit" and nan_safe:
                    dud_args = tuple(dud_arrays) + rest
                    _rough_check(result, fn(*dud_args, **kwargs), key=key)
                # mask out those nan values
                if ma_safe:
                    dud_arrays = [np.ma.masked_invalid(a) for a in dud_arrays]
                    dud_args = tuple(dud_arrays) + rest
                    _rough_check(result, fn(*dud_args, **kwargs), key=key)
                    if ma_very_safe:
                        kwargs[weight_arg] = None
                        _rough_check(result, fn(*dud_args, **kwargs), key=key)
                del dud_arrays, dud_args, dud_weights

            # WEIGHTS CHECK 3: DUPLICATE DATA (DUMB SPLITTING)
            if dup_test:
                dup_arrays = [np.append(a, a, axis=axis) for a in arrays]
                dup_weights = np.append(weights, weights) / 2.0
                dup_args = tuple(dup_arrays) + rest
                kwargs[weight_arg] = dup_weights
                _rough_check(result, fn(*dup_args, **kwargs), key=key)
                del dup_args, dup_arrays, dup_weights

            # WEIGHT CHECK 3: RANDOM SPLITTING
            if split_test and split_per > 0:
                split_arrays, split_weights = _rand_split(arrays, weights, axis, split_per=split_per, seed=seed) 
                split_args = tuple(split_arrays) + rest
                kwargs[weight_arg] = split_weights
                _rough_check(result, fn(*split_args, **kwargs), key=key)
        except NotImplementedError as e:
            # when some combination of arguments makes weighting impossible,
            #  this is the desired response
            if not silent:
                warnings.warn("%s NotImplemented weights: %s" % (fn.__name__, e))
        return result
    return wrapped


wcdist = _weight_checked(cdist, default_axis=1, squeeze=False)
wcdist_no_const = _weight_checked(cdist, default_axis=1, squeeze=False, const_test=False)
wpdist = _weight_checked(pdist, default_axis=1, squeeze=False, n_args=1)
wpdist_no_const = _weight_checked(pdist, default_axis=1, squeeze=False, const_test=False, n_args=1)
wrogerstanimoto = _weight_checked(rogerstanimoto)
wmatching = whamming = _weight_checked(hamming, dud_test=False)
wyule = _weight_checked(yule)
wdice = _weight_checked(dice)
wcosine = _weight_checked(cosine)
wcorrelation = _weight_checked(correlation)
wminkowski = _weight_checked(minkowski, const_test=False)
wjaccard = _weight_checked(jaccard)
weuclidean = _weight_checked(euclidean, const_test=False)
wsqeuclidean = _weight_checked(sqeuclidean, const_test=False)
wbraycurtis = _weight_checked(braycurtis)
wcanberra = _weight_checked(canberra, const_test=False)
wsokalsneath = _weight_checked(sokalsneath)
wsokalmichener = _weight_checked(sokalmichener)
wrussellrao = _weight_checked(russellrao)


class TestCdist(TestCase):
    def setUp(self):
        self.rnd_eo_names = ['random-float32-data', 'random-int-data',
                             'random-uint-data', 'random-double-data',
                             'random-bool-data']
        self.valid_upcasts = {'bool': [np.uint, np.int_, np.float32, np.double],
                              'uint': [np.int_, np.float32, np.double],
                              'int': [np.float32, np.double],
                              'float32': [np.double]}

    def test_cdist_euclidean_random(self):
        eps = 1e-07
        # Get the data: the input matrix and the right output.
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = wcdist_no_const(X1, X2, 'euclidean')
        Y2 = wcdist_no_const(X1, X2, 'test_euclidean')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_euclidean_random_unicode(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = wcdist_no_const(X1, X2, u('euclidean'))
        Y2 = wcdist_no_const(X1, X2, u('test_euclidean'))
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_sqeuclidean_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = wcdist_no_const(X1, X2, 'sqeuclidean')
        Y2 = wcdist_no_const(X1, X2, 'test_sqeuclidean')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_cityblock_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = wcdist_no_const(X1, X2, 'cityblock')
        Y2 = wcdist_no_const(X1, X2, 'test_cityblock')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_hamming_double_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = wcdist(X1, X2, 'hamming')
        Y2 = cdist(X1, X2, 'test_hamming')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_hamming_bool_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = wcdist(X1, X2, 'hamming')
        Y2 = cdist(X1, X2, 'test_hamming')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_jaccard_double_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = wcdist(X1, X2, 'jaccard')
        Y2 = wcdist(X1, X2, 'test_jaccard')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_jaccard_bool_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = wcdist(X1, X2, 'jaccard')
        Y2 = wcdist(X1, X2, 'test_jaccard')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_chebychev_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = wcdist(X1, X2, 'chebychev')
        Y2 = wcdist(X1, X2, 'test_chebychev')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_minkowski_random_p3d8(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = wcdist_no_const(X1, X2, 'minkowski', p=3.8)
        Y2 = wcdist_no_const(X1, X2, 'test_minkowski', p=3.8)
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_minkowski_random_p4d6(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = wcdist_no_const(X1, X2, 'minkowski', p=4.6)
        Y2 = wcdist_no_const(X1, X2, 'test_minkowski', p=4.6)
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_minkowski_random_p1d23(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = wcdist_no_const(X1, X2, 'minkowski', p=1.23)
        Y2 = wcdist_no_const(X1, X2, 'test_minkowski', p=1.23)
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
        Y1 = wcdist(X1, X2, 'cosine')

        # Naive implementation
        def norms(X):
            return np.linalg.norm(X, axis=1).reshape(-1, 1)

        Y2 = 1 - np.dot((X1 / norms(X1)), (X2 / norms(X2)).T)

        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_correlation_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = wcdist(X1, X2, 'correlation')
        Y2 = cdist(X1, X2, 'test_correlation')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_mahalanobis_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1']
        X2 = eo['cdist-X2']
        Y1 = cdist(X1, X2, 'mahalanobis')
        Y2 = cdist(X1, X2, 'test_mahalanobis')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_mahalanobis(self):
        # 1-dimensional observations
        x1 = np.array([[2], [3]])
        x2 = np.array([[2], [5]])
        dist = cdist(x1, x2, metric='mahalanobis')
        assert_allclose(dist, [[0.0, np.sqrt(4.5)], [np.sqrt(0.5), np.sqrt(2)]])

        # 2-dimensional observations
        x1 = np.array([[0, 0], [-1, 0]])
        x2 = np.array([[0, 2], [1, 0], [0, -2]])
        dist = cdist(x1, x2, metric='mahalanobis')
        rt2 = np.sqrt(2)
        assert_allclose(dist, [[rt2, rt2, rt2], [2, 2*rt2, 2]])

        # Too few observations
        assert_raises(ValueError,
                      cdist, [[0, 1]], [[2, 3]], metric='mahalanobis')

    def test_cdist_canberra_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = wcdist_no_const(X1, X2, 'canberra')
        Y2 = wcdist_no_const(X1, X2, 'test_canberra')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_braycurtis_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = wcdist(X1, X2, 'braycurtis')
        Y2 = cdist(X1, X2, 'test_braycurtis')
        if verbose > 2:
            print(Y1, Y2)
            print((Y1-Y2).max())
        _assert_within_tol(Y1, Y2, eps)

    def test_cdist_yule_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = wcdist(X1, X2, 'yule')
        Y2 = cdist(X1, X2, 'test_yule')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_matching_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            Y1 = wcdist(X1, X2, 'matching')
            Y2 = cdist(X1, X2, 'test_matching')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_kulsinski_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = wcdist(X1, X2, 'kulsinski')
        Y2 = cdist(X1, X2, 'test_kulsinski')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_dice_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = wcdist(X1, X2, 'dice')
        Y2 = cdist(X1, X2, 'test_dice')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_rogerstanimoto_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = wcdist(X1, X2, 'rogerstanimoto')
        Y2 = cdist(X1, X2, 'test_rogerstanimoto')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_russellrao_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = wcdist(X1, X2, 'russellrao')
        Y2 = cdist(X1, X2, 'test_russellrao')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_sokalmichener_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = wcdist(X1, X2, 'sokalmichener')
        Y2 = cdist(X1, X2, 'test_sokalmichener')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_sokalsneath_random(self):
        eps = 1e-07
        X1 = eo['cdist-X1'] < 0.5
        X2 = eo['cdist-X2'] < 0.5
        Y1 = wcdist(X1, X2, 'sokalsneath')
        Y2 = cdist(X1, X2, 'test_sokalsneath')
        _assert_within_tol(Y1, Y2, eps, verbose > 2)

    def test_cdist_custom_notdouble(self):
        class myclass(object):
            pass

        def dummy_metric(x, y):
            if not isinstance(x[0], myclass) or not isinstance(y[0], myclass):
                raise ValueError("Type has been changed")
            return 1.123
        data = np.array([[myclass()]], dtype=object)
        cdist_y = cdist(data, data, metric=dummy_metric)
        right_y = 1.123
        assert_equal(cdist_y, right_y, verbose=verbose > 2)

    def _check_calling_conventions(self, X1, X2, metric, eo_name, eps=1e-04, **kwargs):
        # helper function for test_cdist_calling_conventions
        try:
            y1 = cdist(X1, X2, metric=metric, **kwargs)
            y2 = cdist(X1, X2, metric=eval(metric), **kwargs)
            y3 = cdist(X1, X2, metric="test_" + metric, **kwargs)
        except Exception as e:
            e_cls = e.__class__
            if verbose > 2:
                print(e_cls.__name__)
                print(e)
            try:
                assert_raises(Exception, cdist, X1, X2, metric=metric, **kwargs)
                assert_raises(Exception, cdist, X1, X2, metric=eval(metric), **kwargs)
                assert_raises(Exception, cdist, X1, X2, metric="test_" + metric, **kwargs)
            except Exception as e2:
                raise type(e)((e2, e, metric, eo_name))
        else:
            try:
                _assert_within_tol(y1, y2, rtol=eps, verbose_=verbose > 2)
                _assert_within_tol(y2, y3, rtol=eps, verbose_=verbose > 2)
            except Exception as e:
                raise type(e)((e, metric, eo_name))

    def test_cdist_calling_conventions(self):
        # Ensures that specifying the metric with a str or scipy function
        # gives the same behaviour (i.e. same result or same exception).
        # NOTE: The correctness should be checked within each metric tests.
        for eo_name in self.rnd_eo_names:
            X1 = eo[eo_name][:, ::-1]
            X2 = eo[eo_name][:-3:2]
            for metric in _metrics:
                if verbose > 2:
                    print("testing: ", metric, " with: ", eo_name)
                if metric in {'dice', 'yule', 'kulsinski', 'matching', 'rogerstanimoto',
                              'russellrao', 'sokalmichener', 'sokalsneath'} and 'bool' not in eo_name:
                    # python version permits non-bools e.g. for fuzzy logic
                    continue
                self._check_calling_conventions(X1, X2, metric, eo_name)

                # Testing built-in metrics with extra args
                if metric == "seuclidean":
                    X12 = np.vstack([X1, X2]).astype(np.double)
                    V = np.var(X12, axis=0, ddof=1)
                    self._check_calling_conventions(X1, X2, metric, eo_name, V=V)
                elif metric == "mahalanobis":
                    X12 = np.vstack([X1, X2]).astype(np.double)
                    V = np.atleast_2d(np.cov(X12.T))
                    VI = np.array(np.linalg.inv(V).T)
                    self._check_calling_conventions(X1, X2, metric, eo_name, VI=VI)

    def test_cdist_dtype_equivalence(self):
        # Tests that the result is not affected by type up-casting
        eps = 1e-07
        tests = [(eo['random-bool-data'], self.valid_upcasts['bool']),
                 (eo['random-uint-data'], self.valid_upcasts['uint']),
                 (eo['random-int-data'], self.valid_upcasts['int']),
                 (eo['random-float32-data'], self.valid_upcasts['float32'])]
        for metric in _metrics:
            for test in tests:
                X1 = test[0]
                try:
                    y1 = cdist(X1, metric=metric)
                except Exception as e:
                    e_cls = e.__class__
                    if verbose > 2:
                        print(e_cls.__name__)
                        print(e)
                    for new_type in test[1]:
                        X2 = new_type(test[0])
                        assert_raises(e_cls, cdist, X2, metric=metric)
                else:
                    for new_type in test[1]:
                        y2 = cdist(new_type(X1), metric=metric)
                        _assert_within_tol(y1, y2, eps, verbose > 2)


class TestPdist(TestCase):
    def setUp(self):
        self.rnd_eo_names = ['random-float32-data', 'random-int-data',
                             'random-uint-data', 'random-double-data',
                             'random-bool-data']
        self.valid_upcasts = {'bool': [np.uint, np.int_, np.float32, np.double],
                              'uint': [np.int_, np.float32, np.double],
                              'int': [np.float32, np.double],
                              'float32': [np.double]}

    def test_pdist_euclidean_random(self):
        eps = 1e-07
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-euclidean']
        Y_test1 = wpdist_no_const(X, 'euclidean')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_euclidean_random_u(self):
        eps = 1e-07
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-euclidean']
        Y_test1 = wpdist_no_const(X, u('euclidean'))
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_euclidean_random_float32(self):
        eps = 1e-07
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-euclidean']
        Y_test1 = wpdist_no_const(X, 'euclidean')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_euclidean_random_nonC(self):
        eps = 1e-07
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-euclidean']
        Y_test2 = wpdist_no_const(X, 'test_euclidean')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_euclidean_iris_double(self):
        eps = 1e-07
        X = eo['iris']
        Y_right = eo['pdist-euclidean-iris']
        Y_test1 = wpdist_no_const(X, 'euclidean')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_euclidean_iris_float32(self):
        eps = 1e-06
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-euclidean-iris']
        Y_test1 = wpdist_no_const(X, 'euclidean')
        _assert_within_tol(Y_test1, Y_right, eps, verbose > 2)

    def test_pdist_euclidean_iris_nonC(self):
        # Test pdist(X, 'test_euclidean') [the non-C implementation] on the
        # Iris data set.
        eps = 1e-07
        X = eo['iris']
        Y_right = eo['pdist-euclidean-iris']
        Y_test2 = wpdist_no_const(X, 'test_euclidean')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_seuclidean_random(self):
        eps = 1e-05
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-seuclidean']
        Y_test1 = pdist(X, 'seuclidean')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_seuclidean_random_float32(self):
        eps = 1e-05
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-seuclidean']
        Y_test1 = pdist(X, 'seuclidean')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_seuclidean_random_nonC(self):
        # Test pdist(X, 'test_sqeuclidean') [the non-C implementation]
        eps = 1e-05
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-seuclidean']
        Y_test2 = pdist(X, 'test_seuclidean')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_seuclidean_iris(self):
        eps = 1e-05
        X = eo['iris']
        Y_right = eo['pdist-seuclidean-iris']
        Y_test1 = pdist(X, 'seuclidean')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_seuclidean_iris_float32(self):
        # Tests pdist(X, 'seuclidean') on the Iris data set (float32).
        eps = 1e-05
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-seuclidean-iris']
        Y_test1 = pdist(X, 'seuclidean')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_seuclidean_iris_nonC(self):
        # Test pdist(X, 'test_sqeuclidean') [the non-C implementation] on the
        # Iris data set.
        eps = 1e-05
        X = eo['iris']
        Y_right = eo['pdist-seuclidean-iris']
        Y_test2 = pdist(X, 'test_seuclidean')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_cosine_random(self):
        eps = 1e-08
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-cosine']
        Y_test1 = wpdist(X, 'cosine')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_cosine_random_float32(self):
        eps = 1e-08
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-cosine']
        Y_test1 = wpdist(X, 'cosine')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_cosine_random_nonC(self):
        # Test pdist(X, 'test_cosine') [the non-C implementation]
        eps = 1e-08
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-cosine']
        Y_test2 = wpdist(X, 'test_cosine')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_cosine_iris(self):
        eps = 1e-08
        X = eo['iris']
        Y_right = eo['pdist-cosine-iris']
        Y_test1 = wpdist(X, 'cosine')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_cosine_iris_float32(self):
        eps = 1e-07
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-cosine-iris']
        Y_test1 = wpdist(X, 'cosine')
        _assert_within_tol(Y_test1, Y_right, eps, verbose > 2)

    def test_pdist_cosine_iris_nonC(self):
        eps = 1e-08
        X = eo['iris']
        Y_right = eo['pdist-cosine-iris']
        Y_test2 = wpdist(X, 'test_cosine')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_cosine_bounds(self):
        # Test adapted from @joernhees's example at gh-5208: case where
        # cosine distance used to be negative. XXX: very sensitive to the
        # specific norm computation.
        x = np.abs(np.random.RandomState(1337).rand(91))
        X = np.vstack([x, x])
        assert_(wpdist(X, 'cosine')[0] >= 0,
                msg='cosine distance should be non-negative')

    def test_pdist_cityblock_random(self):
        eps = 1e-06
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-cityblock']
        Y_test1 = wpdist_no_const(X, 'cityblock')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_cityblock_random_float32(self):
        eps = 1e-06
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-cityblock']
        Y_test1 = wpdist_no_const(X, 'cityblock')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_cityblock_random_nonC(self):
        eps = 1e-06
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-cityblock']
        Y_test2 = wpdist_no_const(X, 'test_cityblock')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_cityblock_iris(self):
        eps = 1e-14
        X = eo['iris']
        Y_right = eo['pdist-cityblock-iris']
        Y_test1 = wpdist_no_const(X, 'cityblock')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_cityblock_iris_float32(self):
        eps = 1e-06
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-cityblock-iris']
        Y_test1 = wpdist_no_const(X, 'cityblock')
        _assert_within_tol(Y_test1, Y_right, eps, verbose > 2)

    def test_pdist_cityblock_iris_nonC(self):
        # Test pdist(X, 'test_cityblock') [the non-C implementation] on the
        # Iris data set.
        eps = 1e-14
        X = eo['iris']
        Y_right = eo['pdist-cityblock-iris']
        Y_test2 = wpdist_no_const(X, 'test_cityblock')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_correlation_random(self):
        eps = 1e-07
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-correlation']
        Y_test1 = wpdist(X, 'correlation')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_correlation_random_float32(self):
        eps = 1e-07
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-correlation']
        Y_test1 = wpdist(X, 'correlation')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_correlation_random_nonC(self):
        eps = 1e-07
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-correlation']
        Y_test2 = wpdist(X, 'test_correlation')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_correlation_iris(self):
        eps = 1e-08
        X = eo['iris']
        Y_right = eo['pdist-correlation-iris']
        Y_test1 = wpdist(X, 'correlation')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_correlation_iris_float32(self):
        eps = 1e-07
        X = eo['iris']
        Y_right = np.float32(eo['pdist-correlation-iris'])
        Y_test1 = wpdist(X, 'correlation')
        _assert_within_tol(Y_test1, Y_right, eps, verbose > 2)

    def test_pdist_correlation_iris_nonC(self):
        eps = 1e-08
        X = eo['iris']
        Y_right = eo['pdist-correlation-iris']
        Y_test2 = wpdist(X, 'test_correlation')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_minkowski_random(self):
        eps = 1e-05
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-minkowski-3.2']
        Y_test1 = wpdist_no_const(X, 'minkowski', 3.2)
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_minkowski_random_float32(self):
        eps = 1e-05
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-minkowski-3.2']
        Y_test1 = wpdist_no_const(X, 'minkowski', 3.2)
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_minkowski_random_nonC(self):
        eps = 1e-05
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-minkowski-3.2']
        Y_test2 = wpdist_no_const(X, 'test_minkowski', 3.2)
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_minkowski_3_2_iris(self):
        eps = 1e-07
        X = eo['iris']
        Y_right = eo['pdist-minkowski-3.2-iris']
        Y_test1 = wpdist_no_const(X, 'minkowski', 3.2)
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_minkowski_3_2_iris_float32(self):
        eps = 1e-06
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-minkowski-3.2-iris']
        Y_test1 = wpdist_no_const(X, 'minkowski', 3.2)
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_minkowski_3_2_iris_nonC(self):
        eps = 1e-07
        X = eo['iris']
        Y_right = eo['pdist-minkowski-3.2-iris']
        Y_test2 = wpdist_no_const(X, 'test_minkowski', 3.2)
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_minkowski_5_8_iris(self):
        eps = 1e-07
        X = eo['iris']
        Y_right = eo['pdist-minkowski-5.8-iris']
        Y_test1 = wpdist_no_const(X, 'minkowski', 5.8)
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_minkowski_5_8_iris_float32(self):
        eps = 1e-06
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-minkowski-5.8-iris']
        Y_test1 = wpdist_no_const(X, 'minkowski', 5.8)
        _assert_within_tol(Y_test1, Y_right, eps, verbose > 2)

    def test_pdist_minkowski_5_8_iris_nonC(self):
        eps = 1e-07
        X = eo['iris']
        Y_right = eo['pdist-minkowski-5.8-iris']
        Y_test2 = wpdist_no_const(X, 'test_minkowski', 5.8)
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_mahalanobis(self):
        # 1-dimensional observations
        x = np.array([2.0, 2.0, 3.0, 5.0]).reshape(-1, 1)
        dist = pdist(x, metric='mahalanobis')
        assert_allclose(dist, [0.0, np.sqrt(0.5), np.sqrt(4.5),
                               np.sqrt(0.5), np.sqrt(4.5), np.sqrt(2.0)])

        # 2-dimensional observations
        x = np.array([[0, 0], [-1, 0], [0, 2], [1, 0], [0, -2]])
        dist = pdist(x, metric='mahalanobis')
        rt2 = np.sqrt(2)
        assert_allclose(dist, [rt2, rt2, rt2, rt2, 2, 2*rt2, 2, 2, 2*rt2, 2])

        # Too few observations
        assert_raises(ValueError,
                      wpdist, [[0, 1], [2, 3]], metric='mahalanobis')

    def test_pdist_hamming_random(self):
        eps = 1e-07
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-hamming']
        Y_test1 = wpdist(X, 'hamming')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_hamming_random_float32(self):
        eps = 1e-07
        X = np.float32(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']
        Y_test1 = wpdist(X, 'hamming')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_hamming_random_nonC(self):
        eps = 1e-07
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-hamming']
        Y_test2 = wpdist(X, 'test_hamming')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_dhamming_random(self):
        eps = 1e-07
        X = np.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']
        Y_test1 = wpdist(X, 'hamming')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_dhamming_random_float32(self):
        eps = 1e-07
        X = np.float32(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']
        Y_test1 = wpdist(X, 'hamming')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_dhamming_random_nonC(self):
        eps = 1e-07
        X = np.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-hamming']
        Y_test2 = wpdist(X, 'test_hamming')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_jaccard_random(self):
        eps = 1e-08
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-jaccard']
        Y_test1 = wpdist(X, 'jaccard')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_jaccard_random_float32(self):
        eps = 1e-08
        X = np.float32(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']
        Y_test1 = wpdist(X, 'jaccard')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_jaccard_random_nonC(self):
        eps = 1e-08
        X = eo['pdist-boolean-inp']
        Y_right = eo['pdist-jaccard']
        Y_test2 = wpdist(X, 'test_jaccard')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_djaccard_random(self):
        eps = 1e-08
        X = np.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']
        Y_test1 = wpdist(X, 'jaccard')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_djaccard_random_float32(self):
        eps = 1e-08
        X = np.float32(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']
        Y_test1 = wpdist(X, 'jaccard')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_djaccard_random_nonC(self):
        eps = 1e-08
        X = np.float64(eo['pdist-boolean-inp'])
        Y_right = eo['pdist-jaccard']
        Y_test2 = wpdist(X, 'test_jaccard')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_chebychev_random(self):
        eps = 1e-08
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-chebychev']
        Y_test1 = wpdist(X, 'chebychev')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_chebychev_random_float32(self):
        eps = 1e-07
        X = np.float32(eo['pdist-double-inp'])
        Y_right = eo['pdist-chebychev']
        Y_test1 = wpdist(X, 'chebychev')
        _assert_within_tol(Y_test1, Y_right, eps, verbose > 2)

    def test_pdist_chebychev_random_nonC(self):
        eps = 1e-08
        X = eo['pdist-double-inp']
        Y_right = eo['pdist-chebychev']
        Y_test2 = wpdist(X, 'test_chebychev')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_chebychev_iris(self):
        eps = 1e-15
        X = eo['iris']
        Y_right = eo['pdist-chebychev-iris']
        Y_test1 = wpdist(X, 'chebychev')
        _assert_within_tol(Y_test1, Y_right, eps)

    def test_pdist_chebychev_iris_float32(self):
        eps = 1e-06
        X = np.float32(eo['iris'])
        Y_right = eo['pdist-chebychev-iris']
        Y_test1 = wpdist(X, 'chebychev')
        _assert_within_tol(Y_test1, Y_right, eps, verbose > 2)

    def test_pdist_chebychev_iris_nonC(self):
        eps = 1e-15
        X = eo['iris']
        Y_right = eo['pdist-chebychev-iris']
        Y_test2 = wpdist(X, 'test_chebychev')
        _assert_within_tol(Y_test2, Y_right, eps)

    def test_pdist_matching_mtica1(self):
        # Test matching(*,*) with mtica example #1 (nums).
        m = wmatching(np.array([1, 0, 1, 1, 0]),
                      np.array([1, 1, 0, 1, 1]))
        m2 = wmatching(np.array([1, 0, 1, 1, 0], dtype=bool),
                       np.array([1, 1, 0, 1, 1], dtype=bool))
        assert_allclose(m, 0.6, rtol=0, atol=1e-10)
        assert_allclose(m2, 0.6, rtol=0, atol=1e-10)

    def test_pdist_matching_mtica2(self):
        # Test matching(*,*) with mtica example #2.
        m = wmatching(np.array([1, 0, 1]),
                      np.array([1, 1, 0]))
        m2 = wmatching(np.array([1, 0, 1], dtype=bool),
                       np.array([1, 1, 0], dtype=bool))
        assert_allclose(m, 2/3, rtol=0, atol=1e-10)
        assert_allclose(m2, 2/3, rtol=0, atol=1e-10)

    def test_pdist_matching_match(self):
        # Test pdist(X, 'matching') to see if the two implementations match on
        # random boolean input data.
        D = eo['random-bool-data']
        B = np.bool_(D)
        if verbose > 2:
            print(B.shape, B.dtype)
        eps = 1e-10
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            y1 = wpdist(B, "matching")
            y2 = wpdist(B, "test_matching")
            y3 = wpdist(D, "test_matching")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y1-y3).max())
        _assert_within_tol(y1, y2, eps)
        _assert_within_tol(y2, y3, eps)

    def test_pdist_jaccard_mtica1(self):
        m = wjaccard(np.array([1, 0, 1, 1, 0]),
                     np.array([1, 1, 0, 1, 1]))
        m2 = wjaccard(np.array([1, 0, 1, 1, 0], dtype=bool),
                      np.array([1, 1, 0, 1, 1], dtype=bool))
        assert_allclose(m, 0.6, rtol=0, atol=1e-10)
        assert_allclose(m2, 0.6, rtol=0, atol=1e-10)

    def test_pdist_jaccard_mtica2(self):
        m = wjaccard(np.array([1, 0, 1]),
                     np.array([1, 1, 0]))
        m2 = wjaccard(np.array([1, 0, 1], dtype=bool),
                      np.array([1, 1, 0], dtype=bool))
        assert_allclose(m, 2/3, rtol=0, atol=1e-10)
        assert_allclose(m2, 2/3, rtol=0, atol=1e-10)

    def test_pdist_jaccard_match(self):
        # Test pdist(X, 'jaccard') to see if the two implementations match on
        # random double input data.
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = wpdist(D, "jaccard")
        y2 = wpdist(D, "test_jaccard")
        y3 = wpdist(np.bool_(D), "test_jaccard")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        _assert_within_tol(y1, y2, eps)
        _assert_within_tol(y2, y3, eps)

    def test_pdist_yule_mtica1(self):
        m = wyule(np.array([1, 0, 1, 1, 0]),
                  np.array([1, 1, 0, 1, 1]))
        m2 = wyule(np.array([1, 0, 1, 1, 0], dtype=bool),
                   np.array([1, 1, 0, 1, 1], dtype=bool))
        if verbose > 2:
            print(m)
        assert_allclose(m, 2, rtol=0, atol=1e-10)
        assert_allclose(m2, 2, rtol=0, atol=1e-10)

    def test_pdist_yule_mtica2(self):
        m = wyule(np.array([1, 0, 1]),
                  np.array([1, 1, 0]))
        m2 = wyule(np.array([1, 0, 1], dtype=bool),
                   np.array([1, 1, 0], dtype=bool))
        if verbose > 2:
            print(m)
        assert_allclose(m, 2, rtol=0, atol=1e-10)
        assert_allclose(m2, 2, rtol=0, atol=1e-10)

    def test_pdist_yule_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = wpdist(D, "yule")
        y2 = wpdist(D, "test_yule")
        y3 = wpdist(np.bool_(D), "test_yule")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        _assert_within_tol(y1, y2, eps)
        _assert_within_tol(y2, y3, eps)

    def test_pdist_dice_mtica1(self):
        m = wdice(np.array([1, 0, 1, 1, 0]),
                  np.array([1, 1, 0, 1, 1]))
        m2 = wdice(np.array([1, 0, 1, 1, 0], dtype=bool),
                   np.array([1, 1, 0, 1, 1], dtype=bool))
        if verbose > 2:
            print(m)
        assert_allclose(m, 3/7, rtol=0, atol=1e-10)
        assert_allclose(m2, 3/7, rtol=0, atol=1e-10)

    def test_pdist_dice_mtica2(self):
        m = wdice(np.array([1, 0, 1]),
                  np.array([1, 1, 0]))
        m2 = wdice(np.array([1, 0, 1], dtype=bool),
                   np.array([1, 1, 0], dtype=bool))
        if verbose > 2:
            print(m)
        assert_allclose(m, 0.5, rtol=0, atol=1e-10)
        assert_allclose(m2, 0.5, rtol=0, atol=1e-10)

    def test_pdist_dice_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = wpdist(D, "dice")
        y2 = wpdist(D, "test_dice")
        y3 = wpdist(D, "test_dice")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        _assert_within_tol(y1, y2, eps)
        _assert_within_tol(y2, y3, eps)

    def test_pdist_sokalsneath_mtica1(self):
        m = wsokalsneath(np.array([1, 0, 1, 1, 0]),
                         np.array([1, 1, 0, 1, 1]))
        m2 = wsokalsneath(np.array([1, 0, 1, 1, 0], dtype=bool),
                          np.array([1, 1, 0, 1, 1], dtype=bool))
        if verbose > 2:
            print(m)
        assert_allclose(m, 3/4, rtol=0, atol=1e-10)
        assert_allclose(m2, 3/4, rtol=0, atol=1e-10)

    def test_pdist_sokalsneath_mtica2(self):
        m = wsokalsneath(np.array([1, 0, 1]),
                         np.array([1, 1, 0]))
        m2 = wsokalsneath(np.array([1, 0, 1], dtype=bool),
                          np.array([1, 1, 0], dtype=bool))
        if verbose > 2:
            print(m)
        assert_allclose(m, 4/5, rtol=0, atol=1e-10)
        assert_allclose(m2, 4/5, rtol=0, atol=1e-10)

    def test_pdist_sokalsneath_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = wpdist(D, "sokalsneath")
        y2 = wpdist(D, "test_sokalsneath")
        y3 = wpdist(np.bool_(D), "test_sokalsneath")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        _assert_within_tol(y1, y2, eps)
        _assert_within_tol(y2, y3, eps)

    def test_pdist_rogerstanimoto_mtica1(self):
        m = wrogerstanimoto(np.array([1, 0, 1, 1, 0]),
                            np.array([1, 1, 0, 1, 1]))
        m2 = wrogerstanimoto(np.array([1, 0, 1, 1, 0], dtype=bool),
                             np.array([1, 1, 0, 1, 1], dtype=bool))
        if verbose > 2:
            print(m)
        assert_allclose(m, 3/4, rtol=0, atol=1e-10)
        assert_allclose(m2, 3/4, rtol=0, atol=1e-10)

    def test_pdist_rogerstanimoto_mtica2(self):
        m = wrogerstanimoto(np.array([1, 0, 1]),
                            np.array([1, 1, 0]))
        m2 = wrogerstanimoto(np.array([1, 0, 1], dtype=bool),
                             np.array([1, 1, 0], dtype=bool))
        if verbose > 2:
            print(m)
        assert_allclose(m, 4/5, rtol=0, atol=1e-10)
        assert_allclose(m2, 4/5, rtol=0, atol=1e-10)

    def test_pdist_rogerstanimoto_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = wpdist(D, "rogerstanimoto")
        y2 = wpdist(D, "test_rogerstanimoto")
        y3 = wpdist(np.bool_(D), "test_rogerstanimoto")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        _assert_within_tol(y1, y2, eps)
        _assert_within_tol(y2, y3, eps)

    def test_pdist_russellrao_mtica1(self):
        m = wrussellrao(np.array([1, 0, 1, 1, 0]),
                        np.array([1, 1, 0, 1, 1]))
        m2 = wrussellrao(np.array([1, 0, 1, 1, 0], dtype=bool),
                         np.array([1, 1, 0, 1, 1], dtype=bool))
        if verbose > 2:
            print(m)
        assert_allclose(m, 3/5, rtol=0, atol=1e-10)
        assert_allclose(m2, 3/5, rtol=0, atol=1e-10)

    def test_pdist_russellrao_mtica2(self):
        m = wrussellrao(np.array([1, 0, 1]),
                        np.array([1, 1, 0]))
        m2 = wrussellrao(np.array([1, 0, 1], dtype=bool),
                         np.array([1, 1, 0], dtype=bool))
        if verbose > 2:
            print(m)
        assert_allclose(m, 2/3, rtol=0, atol=1e-10)
        assert_allclose(m2, 2/3, rtol=0, atol=1e-10)

    def test_pdist_russellrao_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = wpdist(D, "russellrao")
        y2 = wpdist(D, "test_russellrao")
        y3 = wpdist(np.bool_(D), "test_russellrao")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        _assert_within_tol(y1, y2, eps)
        _assert_within_tol(y2, y3, eps)

    def test_pdist_sokalmichener_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = wpdist(D, "sokalmichener")
        y2 = wpdist(D, "test_sokalmichener")
        y3 = wpdist(np.bool_(D), "test_sokalmichener")
        if verbose > 2:
            print(np.abs(y1-y2).max())
            print(np.abs(y2-y3).max())
        _assert_within_tol(y1, y2, eps)
        _assert_within_tol(y2, y3, eps)

    def test_pdist_kulsinski_match(self):
        D = eo['random-bool-data']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = wpdist(D, "kulsinski")
        y2 = wpdist(D, "test_kulsinski")
        y3 = wpdist(np.bool_(D), "test_kulsinski")
        _assert_within_tol(y1, y2, eps, verbose > 2)
        _assert_within_tol(y2, y3, eps)

    def test_pdist_canberra_match(self):
        D = eo['iris']
        if verbose > 2:
            print(D.shape, D.dtype)
        eps = 1e-10
        y1 = wpdist_no_const(D, "canberra")
        y2 = wpdist_no_const(D, "test_canberra")
        _assert_within_tol(y1, y2, eps, verbose > 2)

    def test_pdist_canberra_ticket_711(self):
        # Test pdist(X, 'canberra') to see if Canberra gives the right result
        # as reported on gh-1238.
        eps = 1e-8
        pdist_y = wpdist_no_const(([3.3], [3.4]), "canberra")
        right_y = 0.01492537
        _assert_within_tol(pdist_y, right_y, eps, verbose > 2)

    def _check_calling_conventions(self, X, metric, eps=1e-07, **kwargs):
        # helper function for test_pdist_calling_conventions
        try:
            y1 = pdist(X, metric=metric, **kwargs)
            y2 = pdist(X, metric=eval(metric), **kwargs)
            y3 = pdist(X, metric="test_" + metric, **kwargs)
        except Exception as e:
            e_cls = e.__class__
            if verbose > 2:
                print(e_cls.__name__)
                print(e)
            try:
                assert_raises(Exception, pdist, X, metric=metric, **kwargs)
                assert_raises(Exception, pdist, X, metric=eval(metric), **kwargs)
                assert_raises(Exception, pdist, X, metric="test_" + metric, **kwargs)
            except Exception as e2:
                raise type(e2)((e2, e, metric))
        else:
            try:
                _assert_within_tol(y1, y2, rtol=eps, verbose_=verbose > 2)
                _assert_within_tol(y2, y3, rtol=eps, verbose_=verbose > 2)
            except Exception as e:
                raise type(e)((e, metric))

    def test_pdist_calling_conventions(self):
        # Ensures that specifying the metric with a str or scipy function
        # gives the same behaviour (i.e. same result or same exception).
        # NOTE: The correctness should be checked within each metric tests.
        # NOTE: Extra args should be checked with a dedicated test
        eps = 1e-07
        for eo_name in self.rnd_eo_names:
            X = eo[eo_name][::2]
            for metric in _metrics:
                if verbose > 2:
                    print("testing: ", metric, " with: ", eo_name)
                if metric in {'dice', 'yule', 'kulsinski', 'matching', 'rogerstanimoto',
                              'russellrao', 'sokalmichener', 'sokalsneath'} and 'bool' not in eo_name:
                    # python version permits non-bools e.g. for fuzzy logic
                    continue
                self._check_calling_conventions(X, metric)

                # Testing built-in metrics with extra args
                if metric == "seuclidean":
                    V = np.var(X.astype(np.double), axis=0, ddof=1).astype(np.double)
                    self._check_calling_conventions(X, metric, V=V)
                elif metric == "mahalanobis":
                    V = np.atleast_2d(np.cov(X.astype(np.double).T))
                    VI = np.array(np.linalg.inv(V).T)
                    self._check_calling_conventions(X, metric, VI=VI)

    def test_pdist_dtype_equivalence(self):
        # Tests that the result is not affected by type up-casting
        eps = 1e-07
        tests = [(eo['random-bool-data'], self.valid_upcasts['bool']),
                 (eo['random-uint-data'], self.valid_upcasts['uint']),
                 (eo['random-int-data'], self.valid_upcasts['int']),
                 (eo['random-float32-data'], self.valid_upcasts['float32'])]
        for metric in _metrics:
            for test in tests:
                X1 = test[0]
                try:
                    y1 = pdist(X1, metric=metric)
                except Exception as e:
                    e_cls = e.__class__
                    if verbose > 2:
                        print(e_cls.__name__)
                        print(e)
                    for new_type in test[1]:
                        X2 = new_type(test[0])
                        assert_raises(e_cls, pdist, X2, metric=metric)
                else:
                    for new_type in test[1]:
                        y2 = pdist(new_type(X1), metric=metric)
                        _assert_within_tol(y1, y2, eps, verbose > 2)


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
            dist1 = wminkowski(x, y, p=1)
            assert_almost_equal(dist1, 3.0)
            dist1p5 = wminkowski(x, y, p=1.5)
            assert_almost_equal(dist1p5, (1.0+2.0**1.5)**(2./3))
            dist2 = wminkowski(x, y, p=2)
            assert_almost_equal(dist2, np.sqrt(5))

    def test_euclidean(self):
        for x, y in self.cases:
            dist = weuclidean(x, y)
            assert_almost_equal(dist, np.sqrt(5))

    def test_sqeuclidean(self):
        for x, y in self.cases:
            dist = wsqeuclidean(x, y)
            assert_almost_equal(dist, 5.0)

    def test_cosine(self):
        for x, y in self.cases:
            dist = wcosine(x, y)
            assert_almost_equal(dist, 1.0 - 18.0/(np.sqrt(14)*np.sqrt(27)))

    def test_correlation(self):
        xm = np.array([-1.0, 0, 1.0])
        ym = np.array([-4.0/3, -4.0/3, 5.0-7.0/3])
        for x, y in self.cases:
            dist = wcorrelation(x, y)
            assert_almost_equal(dist, 1.0 - np.dot(xm, ym)/(norm(xm)*norm(ym)))

    def test_mahalanobis(self):
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([1.0, 1.0, 5.0])
        vi = np.array([[2.0, 1.0, 0.0],[1.0, 2.0, 1.0], [0.0, 1.0, 2.0]])
        for x, y in self.cases:
            dist = mahalanobis(x, y, vi)
            assert_almost_equal(dist, np.sqrt(6.0))


class TestSquareForm(object):
    checked_dtypes = [np.float64, np.float32, np.int32, np.int8, bool]

    def test_squareform_matrix(self):
        for dtype in self.checked_dtypes:
            yield self.check_squareform_matrix, dtype

    def test_squareform_vector(self):
        for dtype in self.checked_dtypes:
            yield self.check_squareform_vector, dtype

    def check_squareform_matrix(self, dtype):
        A = np.zeros((0,0), dtype=dtype)
        rA = squareform(A)
        assert_equal(rA.shape, (0,))
        assert_equal(rA.dtype, dtype)

        A = np.zeros((1,1), dtype=dtype)
        rA = squareform(A)
        assert_equal(rA.shape, (0,))
        assert_equal(rA.dtype, dtype)

        A = np.array([[0,4.2],[4.2,0]], dtype=dtype)
        rA = squareform(A)
        assert_equal(rA.shape, (1,))
        assert_equal(rA.dtype, dtype)
        assert_array_equal(rA, np.array([4.2], dtype=dtype))

    def check_squareform_vector(self, dtype):
        v = np.zeros((0,), dtype=dtype)
        rv = squareform(v)
        assert_equal(rv.shape, (1,1))
        assert_equal(rv.dtype, dtype)
        assert_array_equal(rv, [[0]])

        v = np.array([8.3], dtype=dtype)
        rv = squareform(v)
        assert_equal(rv.shape, (2,2))
        assert_equal(rv.dtype, dtype)
        assert_array_equal(rv, np.array([[0,8.3],[8.3,0]], dtype=dtype))

    def test_squareform_multi_matrix(self):
        for n in xrange(2, 5):
            yield self.check_squareform_multi_matrix, n

    def check_squareform_multi_matrix(self, n):
        X = np.random.rand(n, 4)
        Y = wpdist_no_const(X)
        assert_equal(len(Y.shape), 1)
        A = squareform(Y)
        Yr = squareform(A)
        s = A.shape
        k = 0
        if verbose >= 3:
            print(A.shape, Y.shape, Yr.shape)
        assert_equal(len(s), 2)
        assert_equal(len(Yr.shape), 1)
        assert_equal(s[0], s[1])
        for i in xrange(0, s[0]):
            for j in xrange(i+1, s[1]):
                if i != j:
                    assert_equal(A[i, j], Y[k])
                    k += 1
                else:
                    assert_equal(A[i, j], 0)


class TestNumObsY(TestCase):

    def test_num_obs_y_multi_matrix(self):
        for n in xrange(2, 10):
            X = np.random.rand(n, 4)
            Y = wpdist_no_const(X)
            assert_equal(num_obs_y(Y), n)

    def test_num_obs_y_1(self):
        # Tests num_obs_y(y) on a condensed distance matrix over 1
        # observations. Expecting exception.
        assert_raises(ValueError, self.check_y, 1)

    def test_num_obs_y_2(self):
        # Tests num_obs_y(y) on a condensed distance matrix over 2
        # observations.
        assert_(self.check_y(2))

    def test_num_obs_y_3(self):
        assert_(self.check_y(3))

    def test_num_obs_y_4(self):
        assert_(self.check_y(4))

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
        assert_(self.check_y(n))

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
            Y = wpdist_no_const(X)
            A = squareform(Y)
            if verbose >= 3:
                print(A.shape, Y.shape)
            assert_equal(num_obs_dm(A), n)

    def test_num_obs_dm_0(self):
        # Tests num_obs_dm(D) on a 0x0 distance matrix. Expecting exception.
        assert_(self.check_D(0))

    def test_num_obs_dm_1(self):
        # Tests num_obs_dm(D) on a 1x1 distance matrix.
        assert_(self.check_D(1))

    def test_num_obs_dm_2(self):
        assert_(self.check_D(2))

    def test_num_obs_dm_3(self):
        assert_(self.check_D(2))

    def test_num_obs_dm_4(self):
        assert_(self.check_D(4))

    def check_D(self, n):
        return num_obs_dm(self.make_D(n)) == n

    def make_D(self, n):
        return np.random.rand(n, n)


def is_valid_dm_throw(D):
    return is_valid_dm(D, throw=True)


class TestIsValidDM(TestCase):

    def test_is_valid_dm_improper_shape_1D_E(self):
        D = np.zeros((5,), dtype=np.double)
        assert_raises(ValueError, is_valid_dm_throw, (D))

    def test_is_valid_dm_improper_shape_1D_F(self):
        D = np.zeros((5,), dtype=np.double)
        assert_equal(is_valid_dm(D), False)

    def test_is_valid_dm_improper_shape_3D_E(self):
        D = np.zeros((3,3,3), dtype=np.double)
        assert_raises(ValueError, is_valid_dm_throw, (D))

    def test_is_valid_dm_improper_shape_3D_F(self):
        D = np.zeros((3,3,3), dtype=np.double)
        assert_equal(is_valid_dm(D), False)

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
        assert_equal(is_valid_dm(D), False)

    def test_is_valid_dm_asymmetric_E(self):
        y = np.random.rand(10)
        D = squareform(y)
        D[1,3] = D[3,1] + 1
        assert_raises(ValueError, is_valid_dm_throw, (D))

    def test_is_valid_dm_asymmetric_F(self):
        y = np.random.rand(10)
        D = squareform(y)
        D[1,3] = D[3,1] + 1
        assert_equal(is_valid_dm(D), False)

    def test_is_valid_dm_correct_1_by_1(self):
        D = np.zeros((1,1), dtype=np.double)
        assert_equal(is_valid_dm(D), True)

    def test_is_valid_dm_correct_2_by_2(self):
        y = np.random.rand(1)
        D = squareform(y)
        assert_equal(is_valid_dm(D), True)

    def test_is_valid_dm_correct_3_by_3(self):
        y = np.random.rand(3)
        D = squareform(y)
        assert_equal(is_valid_dm(D), True)

    def test_is_valid_dm_correct_4_by_4(self):
        y = np.random.rand(6)
        D = squareform(y)
        assert_equal(is_valid_dm(D), True)

    def test_is_valid_dm_correct_5_by_5(self):
        y = np.random.rand(10)
        D = squareform(y)
        assert_equal(is_valid_dm(D), True)


def is_valid_y_throw(y):
    return is_valid_y(y, throw=True)


class TestIsValidY(TestCase):
    # If test case name ends on "_E" then an exception is expected for the
    # given input, if it ends in "_F" then False is expected for the is_valid_y
    # check.  Otherwise the input is expected to be valid.

    def test_is_valid_y_improper_shape_2D_E(self):
        y = np.zeros((3,3,), dtype=np.double)
        assert_raises(ValueError, is_valid_y_throw, (y))

    def test_is_valid_y_improper_shape_2D_F(self):
        y = np.zeros((3,3,), dtype=np.double)
        assert_equal(is_valid_y(y), False)

    def test_is_valid_y_improper_shape_3D_E(self):
        y = np.zeros((3,3,3), dtype=np.double)
        assert_raises(ValueError, is_valid_y_throw, (y))

    def test_is_valid_y_improper_shape_3D_F(self):
        y = np.zeros((3,3,3), dtype=np.double)
        assert_equal(is_valid_y(y), False)

    def test_is_valid_y_correct_2_by_2(self):
        y = self.correct_n_by_n(2)
        assert_equal(is_valid_y(y), True)

    def test_is_valid_y_correct_3_by_3(self):
        y = self.correct_n_by_n(3)
        assert_equal(is_valid_y(y), True)

    def test_is_valid_y_correct_4_by_4(self):
        y = self.correct_n_by_n(4)
        assert_equal(is_valid_y(y), True)

    def test_is_valid_y_correct_5_by_5(self):
        y = self.correct_n_by_n(5)
        assert_equal(is_valid_y(y), True)

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
    assert_raises(ValueError, wminkowski, [1, 2], [3, 4], p)
    assert_raises(ValueError, wminkowski, [1, 2], [3, 4], p, [1, 1])


def test_canberra():
    # Regression test for ticket #1430.
    assert_equal(wcanberra([1,2,3], [2,4,6]), 1)
    assert_equal(wcanberra([1,1,0,0], [1,0,1,0]), 2)


def test_braycurtis():
    # Regression test for ticket #1430.
    assert_almost_equal(wbraycurtis([1,2,3], [2,4,6]), 1./3, decimal=15)
    assert_almost_equal(wbraycurtis([1,1,0,0], [1,0,1,0]), 0.5, decimal=15)


def test_euclideans():
    # Regression test for ticket #1328.
    x1 = np.array([1, 1, 1])
    x2 = np.array([0, 0, 0])

    # Basic test of the calculation.
    assert_almost_equal(wsqeuclidean(x1, x2), 3.0, decimal=14)
    assert_almost_equal(weuclidean(x1, x2), np.sqrt(3), decimal=14)

    # Check flattening for (1, N) or (N, 1) inputs
    assert_almost_equal(weuclidean(x1[np.newaxis, :], x2[np.newaxis, :]),
                        np.sqrt(3), decimal=14)
    assert_almost_equal(wsqeuclidean(x1[np.newaxis, :], x2[np.newaxis, :]),
                        3.0, decimal=14)
    assert_almost_equal(wsqeuclidean(x1[:, np.newaxis], x2[:, np.newaxis]),
                        3.0, decimal=14)

    # Distance metrics only defined for vectors (= 1-D)
    x = np.arange(4).reshape(2, 2)
    assert_raises(ValueError, weuclidean, x, x)
    assert_raises(ValueError, wsqeuclidean, x, x)

    # Another check, with random data.
    rs = np.random.RandomState(1234567890)
    x = rs.rand(10)
    y = rs.rand(10)
    d1 = weuclidean(x, y)
    d2 = wsqeuclidean(x, y)
    assert_almost_equal(d1**2, d2, decimal=14)


def test_hamming_unequal_length():
    # Regression test for gh-4290.
    x = [0, 0, 1]
    y = [1, 0, 1, 0]
    # Used to give an AttributeError from ndarray.mean called on bool
    assert_raises(ValueError, whamming, x, y)


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
    assert_allclose(whamming(a, b), desired)


def test_sqeuclidean_dtypes():
    # Assert that sqeuclidean returns the right types of values.
    # Integer types should be converted to floating for stability.
    # Floating point types should be the same as the input.
    x = [1, 2, 3]
    y = [4, 5, 6]

    for dtype in [np.int8, np.int16, np.int32, np.int64]:
        d = wsqeuclidean(np.asarray(x, dtype=dtype), np.asarray(y, dtype=dtype))
        assert_(np.issubdtype(d.dtype, np.floating))

    for dtype in [np.uint8, np.uint16, np.uint32, np.uint64]:
        d1 = wsqeuclidean([0], np.asarray([-1], dtype=dtype))
        d2 = wsqeuclidean(np.asarray([-1], dtype=dtype), [0])

        assert_equal(d1, d2)
        assert_equal(d1, np.float64(np.iinfo(dtype).max) ** 2)

    dtypes = [np.float32, np.float64, np.complex64, np.complex128]
    for dtype in ['float16', 'float128']:
        # These aren't present in older numpy versions; float128 may also not
        # be present on all platforms.
        if hasattr(np, dtype):
            dtypes.append(getattr(np, dtype))

    for dtype in dtypes:
        d = wsqeuclidean(np.asarray(x, dtype=dtype), np.asarray(y, dtype=dtype))
        assert_equal(d.dtype, dtype)


def test_sokalmichener():
    # Test that sokalmichener has the same result for bool and int inputs.
    p = [True, True, False]
    q = [True, False, True]
    x = [int(b) for b in p]
    y = [int(b) for b in q]
    dist1 = wsokalmichener(p, q)
    dist2 = wsokalmichener(x, y)
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
