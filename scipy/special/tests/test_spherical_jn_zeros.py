import numpy as np
from scipy.special import spherical_jn_zeros

def test_spherical_jn_zeros():
    # Test with  l=0 and n_zeros=3
    zeros = spherical_jn_zeros(0, 3)
    assert len(zeros) == 3
    assert np.allclose(zeros, [np.pi, 2*np.pi, 3*np.pi], atol=1e-6)

    # Test with l=1 and n_zeros=2
    zeros = spherical_jn_zeros(1, 2)
    assert len(zeros) == 2
    assert np.allclose(zeros, [4.49340945790906, 7.72525183693770], atol=1e-6)

    # Test with invalid values
    try:
        spherical_jn_zeros(-1, 3)
        assert False, "ValueError: `l` and `n_zeros` have to be natural integers"
    except ValueError:
        pass

    try:
        spherical_jn_zeros(0, 0)
        assert False, "ValueError: `l` and `n_zeros` have to be natural integers"
    except ValueError:
        pass
