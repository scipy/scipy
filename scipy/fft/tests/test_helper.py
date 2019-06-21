from scipy.fft import next_fast_len
from numpy.testing import assert_equal

_5_smooth_numbers = [
    2, 3, 4, 5, 6, 8, 9, 10,
    2 * 3 * 5,
    2**3 * 3**5,
    2**3 * 3**3 * 5**2,
]

def test_next_fast_len():
    for n in _5_smooth_numbers:
        assert_equal(next_fast_len(n), n)
