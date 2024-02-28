from numpy.testing import assert_array_almost_equal
from scipy.sparse import dok_matrix


def test_fromkeys_default():
    # test default value value
    edges = [(0,2), (1,0), (2,1)]
    Xdok = dok_matrix.fromkeys(edges)
    assert_array_almost_equal(
        [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
        Xdok.todense()
    )

def test_fromkeys_positional():
    # test with value (positional)
    edges = [(0,2), (1,0), (2,1)]
    Xdok = dok_matrix.fromkeys(edges, -1)
    assert_array_almost_equal(
        [[0, 0, -1], [-1, 0, 0], [0, -1, 0]],
        Xdok.todense()
    )

def test_fromkeys_kwargs():
    # test with value (kwarg)
    edges = [(0,2), (1,0), (2,1)]
    Xdok = dok_matrix.fromkeys(edges, value=-1)
    assert_array_almost_equal(
        [[0, 0, -1], [-1, 0, 0], [0, -1, 0]],
        Xdok.todense()
    )
