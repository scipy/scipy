import numpy as np
from numpy.testing import assert_array_equal
from scipy.sparse.csgraph import DisjointSet


def test_linear_union_sequence():
    n = 10
    dis = DisjointSet(n)
    assert dis.n_components == n

    for i in range(n - 1):
        assert dis.union(i, i + 1)
        assert dis.n_components == n - 1 - i

    roots = [dis.find(i) for i in range(n)]
    assert_array_equal(0, roots)
    assert not dis.union(0, n - 1)


def test_self_unions():
    n = 10
    dis = DisjointSet(n)

    for i in range(n):
        assert not dis.union(i, i)
    assert dis.n_components == n

    roots = [dis.find(i) for i in range(n)]
    assert_array_equal(np.arange(n), roots)


def test_equal_size_ordering():
    n = 20
    dis = DisjointSet(n)
    rng = np.random.RandomState(seed=0)

    for i in range(0, n, 2):
        indices = [i, i + 1]
        rng.shuffle(indices)
        assert dis.union(indices[0], indices[1])
        assert dis.find(i) == i
        assert dis.find(i + 1) == i


def test_binary_tree():
    kmax = 10
    n = 2**kmax
    dis = DisjointSet(n)
    rng = np.random.RandomState(seed=0)

    for k in 2**np.arange(kmax):
        for i in range(0, n, 2 * k):
            r1, r2 = rng.randint(0, k, size=2)
            assert dis.union(i + r1, i + k + r2)

        roots = [dis.find(i) for i in range(n)]
        expected = np.arange(n) - np.arange(n) % (2 * k)
        assert_array_equal(roots, expected)
