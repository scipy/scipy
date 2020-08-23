import numpy as np
from numpy.testing import assert_allclose, assert_equal
import pytest

import scipy.special as sc

class TestInverseErrorFunction:
    def test_compliment(self):
        # Test erfcinv(1 - x) == erfinv(x)
        x = np.linspace(-1, 1, 101)
        assert_allclose(sc.erfcinv(1 - x), sc.erfinv(x), rtol=0, atol=1e-15)

    def test_literal_values(self):
        # calculated via https://keisan.casio.com/exec/system/1180573448
        # for y = 0, 0.1, ... , 0.9
        actual = sc.erfinv(np.linspace(0, 0.9, 10))
        expected = [
            0,
            0.08885599049425768701574,
            0.1791434546212916764928,
            0.27246271472675435562,
            0.3708071585935579290583,
            0.4769362762044698733814,
            0.5951160814499948500193,
            0.7328690779592168522188,
            0.9061938024368232200712,
            1.163087153676674086726,
        ]
        assert_allclose(actual, expected, rtol=0, atol=1e-15)

    @pytest.mark.parametrize(
        'f, x, y',
        [
            (sc.erfinv, -1, -np.inf),
            (sc.erfinv, 0, 0),
            (sc.erfinv, 1, np.inf),
            (sc.erfinv, -100, np.nan),
            (sc.erfinv, 100, np.nan),
            (sc.erfcinv, 0, np.inf),
            (sc.erfcinv, 1, -0.0),
            (sc.erfcinv, 2, -np.inf),
            (sc.erfcinv, -100, np.nan),
            (sc.erfcinv, 100, np.nan),
        ],
        ids=[
            'erfinv at lower bound',
            'erfinv at midpoint',
            'erfinv at upper bound',
            'erfinv below lower bound',
            'erfinv above upper bound',
            'erfcinv at lower bound',
            'erfcinv at midpoint',
            'erfcinv at upper bound',
            'erfcinv below lower bound',
            'erfcinv above upper bound',
        ]
    )
    def test_domain_bounds(self, f, x, y):
        assert_equal(f(x), y)

    def test_erfinv_asympt(self):
        # regression test for gh-12758: erfinv(x) loses precision at small x
        # This test checks erfinv against the Taylor expansion:
        # erf(x) = 2/\sqrt{\pi} (x - x^3 / 3 + O(x^5)),    x\to 0 
        # where we only retain the linear term.
        x = np.array([1e-20, 1e-15, 1e-14, 1e-10, 1e-8, 0.9e-7, 1.1e-7, 1e-6])
        assert_allclose(sc.erfinv(x),
                        np.sqrt(np.pi)/2 * x,
                        rtol=1e-10)

        # also test the roundtrip consistency
        assert_allclose(sc.erf(sc.erfinv(x)),
                        x,
                        rtol=1e-10)

