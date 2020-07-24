import pytest
import numpy as np
from numpy.testing import assert_array_equal
from scipy.sparse.csgraph import DisjointSet


def get_nodes(n, shuffle):
    nodes = np.arange(n)
    if shuffle:
        rng = np.random.RandomState(seed=0)
        rng.shuffle(nodes)
    return nodes


@pytest.mark.parametrize("n", [10, 100])
@pytest.mark.parametrize("shuffle", [False, True])
def test_linear_union_sequence(n, shuffle):
    nodes = get_nodes(n, shuffle)
    dis = DisjointSet()

    for i in range(len(nodes) - 1):
        assert dis.union(nodes[i], nodes[i + 1])
        assert dis.n_components == 1

    roots = [dis.find(i) for i in nodes]
    assert_array_equal(nodes[0], roots)
    assert not dis.union(nodes[0], nodes[-1])


@pytest.mark.parametrize("n", [10, 100])
@pytest.mark.parametrize("shuffle", [False, True])
def test_self_unions(n, shuffle):
    nodes = get_nodes(n, shuffle)
    dis = DisjointSet()

    for i in nodes:
        assert not dis.union(i, i)
    assert dis.n_components == len(nodes)

    roots = [dis.find(i) for i in nodes]
    assert_array_equal(nodes, roots)

@pytest.mark.parametrize("n", [10, 100])
@pytest.mark.parametrize("shuffle", [False, True])
def test_equal_size_ordering(n, shuffle):
    nodes = get_nodes(n, shuffle)
    dis = DisjointSet()
    for i in range(0, len(nodes), 2):
        a, b = nodes[i], nodes[i + 1]
        assert dis.union(a, b)
        assert dis.find(a) == a
        assert dis.find(b) == a


@pytest.mark.parametrize("kmax", [5, 10])
@pytest.mark.parametrize("shuffle", [False, True])
def test_binary_tree(kmax, shuffle):
    n = 2**kmax
    nodes = get_nodes(n, shuffle)
    dis = DisjointSet()
    rng = np.random.RandomState(seed=0)

    for k in 2**np.arange(kmax):
        for i in range(0, n, 2 * k):
            r1, r2 = rng.randint(0, k, size=2)
            assert dis.union(nodes[i + r1], nodes[i + k + r2])

        roots = [dis.find(i) for i in nodes]
        expected_indices = np.arange(n) - np.arange(n) % (2 * k)
        expected = [nodes[i] for i in expected_indices]
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
