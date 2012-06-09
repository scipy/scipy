
import numpy as np
from numpy.testing import TestCase, run_module_suite, assert_equal, \
        assert_array_equal

from scipy.stats import rankdata, tiecorrect


class TestTieCorrect(TestCase):

    def test_empty(self):
        """An empty array requires no correction, should return 1.0."""
        ranks = np.array([], dtype=np.float64)
        c = tiecorrect(ranks)
        assert_equal(c, 1.0)

    def test_one(self):
        """A single element requires no correction, should return 1.0."""
        ranks = np.array([1.0], dtype=np.float64)
        c = tiecorrect(ranks)
        assert_equal(c, 1.0)

    def test_no_correction(self):
        """Arrays with no ties require no correction."""
        ranks = np.arange(2.0)
        c = tiecorrect(ranks)
        assert_equal(c, 1.0)
        ranks = np.arange(3.0)
        c = tiecorrect(ranks)
        assert_equal(c, 1.0)

    def test_basic(self):
        """Check a few basic examples of the tie correction factor."""
        # One tie of two elements
        ranks = np.array([1.0, 2.5, 2.5])
        c = tiecorrect(ranks)
        T = 2.0
        N = ranks.size
        expected = 1.0 - (T**3 - T) / (N**3 - N)
        assert_equal(c, expected)

        # One tie of two elements (same as above, but tie is not at the end)
        ranks = np.array([1.5, 1.5, 3.0])
        c = tiecorrect(ranks)
        T = 2.0
        N = ranks.size
        expected = 1.0 - (T**3 - T) / (N**3 - N)
        assert_equal(c, expected)

        # One tie of three elements
        ranks = np.array([1.0, 3.0, 3.0, 3.0])
        c = tiecorrect(ranks)
        T = 3.0
        N = ranks.size
        expected = 1.0 - (T**3 - T) / (N**3 - N)
        assert_equal(c, expected)

        # Two ties, lengths 2 and 3.
        ranks = np.array([1.5, 1.5, 4.0, 4.0, 4.0])
        c = tiecorrect(ranks)
        T1 = 2.0
        T2 = 3.0
        N = ranks.size
        expected = 1.0 - ((T1**3 - T1) + (T2**3 - T2)) / (N**3 - N)
        assert_equal(c, expected)


class TestRankData(TestCase):

    def test_empty(self):
        """stats.rankdata([]) should return an empty array."""
        a = np.array([], dtype=np.int)
        r = rankdata(a)
        assert_array_equal(r, np.array([], dtype=np.float64))
        r = rankdata([])
        assert_array_equal(r, np.array([], dtype=np.float64))

    def test_one(self):
        """Check stats.rankdata with an array of length 1."""
        data = [100]
        a = np.array(data, dtype=np.int)
        r = rankdata(a)
        assert_array_equal(r, np.array([1.0], dtype=np.float64))
        r = rankdata(data)
        assert_array_equal(r, np.array([1.0], dtype=np.float64))        

    def test_basic(self):
        """Basic tests of stats.rankdata."""
        data = [100, 10, 50]
        expected = np.array([3.0, 1.0, 2.0], dtype=np.float64)
        a = np.array(data, dtype=np.int)
        r = rankdata(a)
        assert_array_equal(r, expected)
        r = rankdata(data)
        assert_array_equal(r, expected)

        data = [40, 10, 30, 10, 50]
        expected = np.array([4.0, 1.5, 3.0, 1.5, 5.0], dtype=np.float64)
        a = np.array(data, dtype=np.int)
        r = rankdata(a)
        assert_array_equal(r, expected)
        r = rankdata(data)
        assert_array_equal(r, expected)

        data = [20, 20, 20, 10, 10, 10]
        expected = np.array([5.0, 5.0, 5.0, 2.0, 2.0, 2.0], dtype=np.float64)
        a = np.array(data, dtype=np.int)
        r = rankdata(a)
        assert_array_equal(r, expected)
        r = rankdata(data)
        assert_array_equal(r, expected)
        # The docstring states explicitly that the argument is flattened.
        a2d = a.reshape(2, 3)
        r = rankdata(a2d)
        assert_array_equal(r, expected)        

    def test_large_int(self):
        data = np.array([2**60, 2**60+1], dtype=np.uint64)
        r = rankdata(data)
        assert_array_equal(r, [1.0, 2.0])

        data = np.array([2**60, 2**60+1], dtype=np.int64)
        r = rankdata(data)
        assert_array_equal(r, [1.0, 2.0])

        data = np.array([2**60, -2**60+1], dtype=np.int64)
        r = rankdata(data)
        assert_array_equal(r, [2.0, 1.0])


if __name__ == "__main__":
    run_module_suite()
