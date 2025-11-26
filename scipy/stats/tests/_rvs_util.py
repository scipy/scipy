import numpy as np


class TailBinner:
    """
    Aggregate values in the tail of a decreasing array of expected frequencies.

    This is a utility class used in the tests of zipfian.rvs().
    It "bins" or aggregates the tail of the array into a single element
    so that no values in the binned array are less than minfreq.

    It provides the method `tail_bin(x)` to aggregate another array using
    that same tail size as was found for `expected` in the constructor.
    This allows the binned expected array and the binned input array `x`
    to be used in a chi-squared or similar test.
    """
    def __init__(self, expected, minfreq):
        """
        The values in `expected` must be positive and decreasing, with
        ``sum(expected) >= minfreq``.
        """
        expected = np.array(expected)
        if expected[0] < minfreq:
            expected2 = np.sum(expected, keepdims=True)
            if expected2[0] < minfreq:
                raise ValueError('must have sum(expected) >= minfreq')
            bin_start = 0
        else:
            n = len(expected)
            if expected[-1] >= minfreq:
                expected2 = expected
                bin_start = n
            else:
                # k is the index of the first element in expected that
                # is less than minfreq
                bin_start = np.argmax(expected < minfreq)
                tail_sum = np.sum(expected[bin_start:])
                # If tail_sum >= minfreq, we bin them together.
                # If tail_sum < minfreq, we combine them with the previous index.
                if tail_sum < minfreq:
                    bin_start -= 1
                    tail_sum += expected[bin_start]
                expected2 = np.concatenate((expected[:bin_start], [tail_sum]))

        self.binned_expected = expected2
        self.bin_start = bin_start
        self.size = len(expected)
        self.total = np.sum(expected)

    def tail_bin(self, x):
        if len(x) != self.size:
            raise ValueError('input x is not the same size as the original '
                             'array that was tail binned.')
        if not np.isclose(np.sum(x), self.total, rtol=1e-8):
            raise ValueError('sum(x) does not match sum(expected)')
        if self.bin_start == self.size:
            return x
        return np.concatenate((x[:self.bin_start],
                               np.sum(x[self.bin_start:], keepdims=True)))
