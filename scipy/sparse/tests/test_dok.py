from numpy.testing import assert_array_equal
from scipy.sparse import dok_array


def test_fromkeys_default():
    # test with default value
    edges = [(0, 2), (1, 0), (2, 1)]
    Xdok = dok_array.fromkeys(edges)
    assert_array_equal(
        Xdok.toarray(),
        [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
    )

def test_fromkeys_positional():
    # test with positional value
    edges = [(0, 2), (1, 0), (2, 1)]
    Xdok = dok_array.fromkeys(edges, -1)
    assert_array_equal(
        Xdok.toarray(),
        [[0, 0, -1], [-1, 0, 0], [0, -1, 0]],
    )

def test_fromkeys_iterator():
    it = ((a, a % 2) for a in range(4))
    Xdok = dok_array.fromkeys(it)
    assert_array_equal(
        Xdok.toarray(),
        [[1, 0], [0, 1], [1, 0], [0, 1]],
    )
