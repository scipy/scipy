import pytest

import numpy as np
from numpy.testing import assert_allclose
import scipy.special as sc


class TestExp1:

    def test_branch_cut(self):
        assert np.isnan(sc.exp1(-1))
        assert sc.exp1(complex(-1, 0)).imag == (
            -sc.exp1(complex(-1, -0.0)).imag
        )

        assert_allclose(
            sc.exp1(complex(-1, 0)),
            sc.exp1(-1 + 1e-20j),
            atol=0,
            rtol=1e-15
        )
        assert_allclose(
            sc.exp1(complex(-1, -0.0)),
            sc.exp1(-1 - 1e-20j),
            atol=0,
            rtol=1e-15
        )

    def test_834(self):
        # Regression test for #834
        a = sc.exp1(-complex(19.9999990))
        b = sc.exp1(-complex(19.9999991))
        assert_allclose(a.imag, b.imag, atol=0, rtol=1e-15)


class TestLogExp1:

    # The expected values were computed with mpmath, e.g.
    #    from mpmath import mp
    #    mp.dps = 75
    #    x = 1000
    #    y = mp.log(mp.expint(1, x))
    #    print(float(y))
    # prints
    #    -1006.9087537832978
    @pytest.mark.parametrize('x, expected',
                             [(2**-50, 3.528714908615874),
                              (0.10, 0.6004417824948862),
                              (0.25, 0.0433301754689424),
                              (0.996, -1.5102201140777374),
                              (1.000, -1.5169319590020456),
                              (1.004, -1.5236351344283192),
                              (2.5, -3.6922885436511623),
                              (25, -28.256715322371637),
                              (200, -205.3032803974004),
                              (450, -456.11146244470496),
                              (600, -606.3985921751421),
                              (1000, -1006.9087537832978),
                              (1e80, -1e80)])
    def test_logexpint1_basic(self, x, expected):
        y = sc.logexp1(x)
        assert_allclose(y, expected, rtol=5e-15)


class TestExpi:

    @pytest.mark.parametrize('result', [
        sc.expi(complex(-1, 0)),
        sc.expi(complex(-1, -0.0)),
        sc.expi(-1)
    ])
    def test_branch_cut(self, result):
        desired = -0.21938393439552027368  # Computed using Mpmath
        assert_allclose(result, desired, atol=0, rtol=1e-14)

    def test_near_branch_cut(self):
        lim_from_above = sc.expi(-1 + 1e-20j)
        lim_from_below = sc.expi(-1 - 1e-20j)
        assert_allclose(
            lim_from_above.real,
            lim_from_below.real,
            atol=0,
            rtol=1e-15
        )
        assert_allclose(
            lim_from_above.imag,
            -lim_from_below.imag,
            atol=0,
            rtol=1e-15
        )

    def test_continuity_on_positive_real_axis(self):
        assert_allclose(
            sc.expi(complex(1, 0)),
            sc.expi(complex(1, -0.0)),
            atol=0,
            rtol=1e-15
        )


class TestExpn:

    def test_out_of_domain(self):
        assert all(np.isnan([sc.expn(-1, 1.0), sc.expn(1, -1.0)]))
