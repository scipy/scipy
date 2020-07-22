import numpy as np
from numpy.testing import assert_array_equal
from scipy.sparse.csgraph import DisjointSet


def test_linear_union_sequence():
    n = 10
    dis = DisjointSet(n)
    assert dis.nc == n

    for i in range(n - 1):
        assert dis.union(i, i + 1)
        assert dis.nc == n - 1 - i

    roots = [dis.find(i) for i in range(n)]
    assert_array_equal(0, roots)
    assert not dis.union(0, n - 1)


def test_self_unions():
    n = 10
    dis = DisjointSet(n)

    for i in range(n):
        assert not dis.union(i, i)
    assert dis.nc == n

    roots = [dis.find(i) for i in range(n)]
    assert_array_equal(np.arange(n), roots)


def test_equal_size_ordering():
    n = 10
    dis = DisjointSet(n)

    for i in range(0, n, 2):
        assert dis.union(i, i + 1)
        assert dis.find(i) == i
        assert dis.find(i + 1) == i


def test_binary_tree():
    kmax = 4
    n = 2**kmax
    dis = DisjointSet(n)

    for k in 2**np.arange(kmax):
        for i in range(0, n, 2 * k):
            assert dis.union(i, i + k)

        roots = [dis.find(i) for i in range(n)]
        expected = np.arange(n) - np.arange(n) % (2 * k)
        assert_array_equal(roots, expected)
