import pytest
import numpy as np
from scipy.sparse.csgraph import DisjointSet
import string


def generate_random_token():
    k = len(string.ascii_letters)
    tokens = list(range(k))
    tokens += list([float(i) for i in range(k)])
    tokens += list(string.ascii_letters)
    tokens += [None for i in range(k)]
    rng = np.random.RandomState(seed=0)

    while 1:
        size = rng.randint(1, 3)
        node = rng.choice(tokens, size)
        if size == 1:
            yield node[0]
        else:
            yield tuple(node)


def get_nodes(n, shuffle):
    nodes = set()
    for node in generate_random_token():
        if node not in nodes:
            nodes.add(node)
            if len(nodes) >= n:
                break
    return list(nodes)


@pytest.mark.parametrize("n", [10, 100])
@pytest.mark.parametrize("shuffle", [False, True])
def test_linear_union_sequence(n, shuffle):
    nodes = get_nodes(n, shuffle)
    dis = DisjointSet()

    for i in range(len(nodes) - 1):
        assert dis.union(nodes[i], nodes[i + 1])
        assert dis.n_components == 1

    roots = [dis.find(i) for i in nodes]
    assert all([nodes[0] == r for r in roots])
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
    assert nodes == roots


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
        assert roots == expected
