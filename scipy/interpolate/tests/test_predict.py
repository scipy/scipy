from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (assert_almost_equal, assert_,
        run_module_suite, assert_raises)

import scipy.interpolate.predict as predict


def test_cache_bigger_zero():
    c = predict.Cache(0)
    assert_raises(RuntimeError, c.add(1))


def test_linear_always_stable_predictor():
    for length in range(20, 50, 5):
        a = [float(item) for item in range(length)]
        a = np.array(a)
        res1 = predict.always_stable_predictor(a)
        assert_almost_equal(res1, length)

        p = predict.Cache(length)
        for item in a:
            p.add(item)
        res2 = p.predict()
        assert_almost_equal(res2, length)


def test_square_always_stable_predictor():
    for length in range(20, 50, 5):
        a = [float(item*item) for item in range(length)]
        a = np.array(a)
        res1 = predict.always_stable_predictor(a)
        assert_(res1 < length*length)
        assert_(0.99*length*length < res1)

        p = predict.Cache(length)
        for item in a:
            p.add(item)
        res2 = p.predict()
        assert_(res2 < length*length)
        assert_(0.99*length*length < res2)


if __name__ == "__main__":
    run_module_suite()
