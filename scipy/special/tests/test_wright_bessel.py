import pytest
import numpy as np
from numpy.testing import assert_equal, assert_allclose

import scipy.special as sc
from scipy.special import wright_bessel


@pytest.mark.parametrize('a', [-0.9, 0, 0.5, 1, 10])
@pytest.mark.parametrize('b_real', [-10., -0.5, 0, 0.5, 10.5])
@pytest.mark.parametrize('b_img', [-10., -0.5, 0, 0.5, 10.5])
def test_wright_bessel_zero(a, b_real, b_img):
    if b_img == 0:
        b = b_real
    else:
        b = complex(b_real, b_img)
    assert_equal(wright_bessel(a, b, 0.), 0.)


@pytest.mark.parametrize('v',      [-10, 1, -0.9, -0.1, 0., 0.1, 0.9, 1, 10])
@pytest.mark.parametrize('z_real', [-10, 1, -0.9, -0.1, 0., 0.1, 0.9, 1, 10])
@pytest.mark.parametrize('z_img',  [-10, 1, -0.9, -0.1, 0., 0.1, 0.9, 1, 10])
def test_wright_bessel_iv(v, z_real, z_img):
    """Test relation of wright_bessel and modified bessel function iv.

    iv(z) = (1/2*⁢z)**v *⁢ Phi⁡(1, v+1; 1/4⁢*z**2).
    See https://dlmf.nist.gov/10.46.E2
    """
    z = complex(z_real, z_img)

    if z != 0:
        wb = wright_bessel(1, v+1, z**2/4., tol=1e-12)
        assert_allclose(np.power(z/2., v) * wb,
                        sc.iv(v, z),
                        rtol=1e-12, atol=1e-12)
