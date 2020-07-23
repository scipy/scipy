import pytest
import numpy as np
from numpy.testing import assert_array_equal
from scipy.sparse.csgraph import DisjointSet


@pytest.mark.parametrize("nodes", [list(range(10)),
                                   [1, 2, 3, 6, 7, 8, 10, 11]])
def test_linear_union_sequence(nodes):
    dis = DisjointSet()

    for i in range(len(nodes) - 1):
        assert dis.union(nodes[i], nodes[i + 1])
        assert dis.n_components == 1

    roots = [dis.find(i) for i in nodes]
    assert_array_equal(nodes[0], roots)
    assert not dis.union(nodes[0], nodes[-1])


@pytest.mark.parametrize("nodes", [list(range(10)),
                                   [1, 2, 3, 6, 7, 8, 10, 11]])
def test_self_unions(nodes):
    dis = DisjointSet()

    for i in nodes:
        assert not dis.union(i, i)
    assert dis.n_components == len(nodes)

    roots = [dis.find(i) for i in nodes]
    assert_array_equal(nodes, roots)


def test_equal_size_ordering():
    n = 20
    dis = DisjointSet()
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
    dis = DisjointSet()
    rng = np.random.RandomState(seed=0)

    for k in 2**np.arange(kmax):
        for i in range(0, n, 2 * k):
            r1, r2 = rng.randint(0, k, size=2)
            assert dis.union(i + r1, i + k + r2)

        roots = [dis.find(i) for i in range(n)]
        expected = np.arange(n) - np.arange(n) % (2 * k)
        assert_array_equal(roots, expected)


def test_input_validation():
    dis = DisjointSet()

    with pytest.raises(TypeError, match="integer"):
        dis.find(1.0)

    with pytest.raises(TypeError, match="integer"):
        dis.find("123")

    with pytest.raises(TypeError, match="integer"):
        dis.union(1, "123")

    with pytest.raises(TypeError, match="integer"):
        dis.union("123", 1)

    with pytest.raises(TypeError, match="integer"):
        dis.union(2.0, "123")

    with pytest.raises(TypeError, match="integer"):
        dis.union("123", 2.0)
