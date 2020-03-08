from __future__ import division, print_function, absolute_import

import itertools

import numpy as np
from numpy.testing import assert_equal
from pytest import raises as assert_raises

from scipy.sparse.linalg import utils


def test_make_system_bad_shape():
    assert_raises(ValueError, utils.make_system, np.zeros((5,3)), None, np.zeros(4), np.zeros(4))


def test_column_vector_preservation():
    # Check that make_system preserves the kind of column vector
    # passed as input.

    arr_types = (np.matrix, np.array)

    vecs = [
        np.zeros((5,)),
        np.zeros((5, 1)),
        np.matrix(np.zeros((5, 1)))
    ]

    bad_vecs = [
        np.zeros((1, 5)),
        np.matrix(np.zeros((1, 5)))
    ]

    for t1, t2 in itertools.product(arr_types, arr_types):
        A = t1(np.zeros((5, 5)))
        M = t2(np.zeros((5, 5)))

        for b, x in itertools.product(vecs, [None] + vecs):
            Av, Mv, xv, bv, postprocess = utils.make_system(A, M, x, b)
            r1 = postprocess(Av.dot(bv) - xv)
            r2 = postprocess(Mv.dot(bv) - xv)
            assert_equal(r1.shape, b.shape)
            assert_equal(r2.shape, b.shape)
            assert_equal(type(r1), type(b))
            assert_equal(type(r2), type(b))

        for b, x in itertools.product(bad_vecs, [None] + vecs + bad_vecs):
            assert_raises(ValueError, utils.make_system, A, M, x, b)

        for b, x in itertools.product(vecs + bad_vecs, bad_vecs):
            assert_raises(ValueError, utils.make_system, A, M, x, b)
