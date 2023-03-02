import pytest
import numpy as np
from numpy.testing import assert_equal
from mpmath import mp
from .mpmath_distribution import Normal, SkewNormal, Beta

class TestInit:
    def test_mpf_params(self):
        # ensure that parameters are of type `mpf`
        dist = SkewNormal(a=[1, 2, 3])
        assert type(dist.a[0]) == type(mp.mpf(0))

    def test_invalid_param(self):
        with pytest.raises(ValueError, match="Invalid parameters {'a'}."):
            Normal(a=[1, 2, 3])

        with pytest.raises(ValueError, match="Invalid parameters..."):
            SkewNormal(b=[1, 2, 3])

        with pytest.raises(ValueError, match="Invalid parameters..."):
            Beta(a=0.5, b=0.5, c=0.5)

    def test_missing_param(self):
        with pytest.raises(ValueError, match="Parameters {'a'} not..."):
            SkewNormal()

        with pytest.raises(ValueError, match="Parameters {'a'} not..."):
            Beta(b=0.5)

    def test_incompatible_shapes(self):
        with pytest.raises(ValueError, match="shape mismatch"):
            Beta(a=[1, 2], b=[1, 2, 3])

    def test_domain(self):
        dist = SkewNormal(a=-np.inf)
        assert_equal(dist._in_domain, False)

        dist = SkewNormal(a=[-np.inf, -10, 0, 10, np.inf])
        assert_equal(dist._in_domain, [False, True, True, True, False])

        dist = Beta(a=[0, 0.25, 0.5, 0.75, 1], b=[[0], [0.5], [1]])
        assert_equal(dist._in_domain, [[False, False, False, False, False],
                                       [False, True, True, True, False],
                                       [False, False, False, False, False]])

    class TestSupport:
        def setup_method(self):
            self.rng = np.random.default_rng(5203023286415660796)

        def test_supported(self):
            dist = SkewNormal(a=1)
            assert_equal(dist.supported(10), True)

            a = self.rng.random(size=2)
            dist = SkewNormal(a=a)
            assert_equal(dist.supported(10), [True, True])

            x = self.rng.random(size=(3, 1))
            dist = SkewNormal(a=[1, 2])
            res = dist.supported(x)
            assert_equal(res, True)
            assert res.shape == (3, 2)
