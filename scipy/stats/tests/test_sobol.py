
from numpy.testing import assert_equal

from scipy.stats._sobol import low_0_bit


class TestLow0Bit:
    def test_examples(self):
        test_vector = [
            # from low_0_bit's docstring
            (0b0000, 1),
            (0b0001, 2),
            (0b0010, 1),
            (0b0101, 2),
            (0b0111, 4),
            # gh-23409
            (2 ** 32 - 1, 33),
            (2 ** 32, 1),
            (2 ** 33 - 1, 34),
            (2 ** 64 - 1, 65),
        ]
        for in_, out in test_vector:
            # Try both specializations of this function if the argument fits
            if in_ <= 2 ** 32 - 1:
                assert_equal(low_0_bit['uint32_t'](in_), out)
            assert_equal(low_0_bit['uint64_t'](in_), out)

