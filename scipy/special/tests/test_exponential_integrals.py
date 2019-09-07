from numpy.testing import assert_allclose
import scipy.special as sc


class TestExp1(object):

    def test_branch_cut(self):
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

    def test_exp1_834(self):
        # Regression test for #834
        a = sc.exp1(-complex(19.9999990))
        b = sc.exp1(-complex(19.9999991))
        assert_allclose(a.imag, b.imag, atol=0, rtol=1e-15)
