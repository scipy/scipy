import pytest
import numpy as np
from scipy.sparse.csgraph import DisjointSet
import string


def generate_random_token():
    k = len(string.ascii_letters)
    tokens = list(np.arange(k, dtype=int))
    tokens += list(np.arange(k, dtype=float))
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


def get_nodes(n):
    nodes = set()
    for node in generate_random_token():
        if node not in nodes:
            nodes.add(node)
            if len(nodes) >= n:
                break
    return list(nodes)


@pytest.mark.parametrize("n", [10, 100])
def test_linear_union_sequence(n):
    nodes = get_nodes(n)
    dis = DisjointSet()

    for i in range(len(nodes) - 1):
        assert not dis.connected(nodes[i], nodes[i + 1])
        assert dis.merge(nodes[i], nodes[i + 1])
        assert dis.connected(nodes[i], nodes[i + 1])
        assert dis.n_components == 1

    assert nodes == list(dis)
    roots = [dis[i] for i in nodes]
    assert all([nodes[0] == r for r in roots])
    assert not dis.merge(nodes[0], nodes[-1])


@pytest.mark.parametrize("n", [10, 100])
def test_self_unions(n):
    nodes = get_nodes(n)
    dis = DisjointSet()

    for x in nodes:
        assert dis.connected(x, x)
        assert not dis.merge(x, x)
        assert dis.connected(x, x)
    assert dis.n_components == len(nodes)

    assert nodes == list(dis)
    roots = [dis[x] for x in nodes]
    assert nodes == roots


@pytest.mark.parametrize("kmax", [5, 10])
def test_binary_tree(kmax):
    n = 2**kmax
    nodes = get_nodes(n)
    dis = DisjointSet()
    rng = np.random.RandomState(seed=0)

    for k in 2**np.arange(kmax):
        for i in range(0, n, 2 * k):
            r1, r2 = rng.randint(0, k, size=2)
            a, b = nodes[i + r1], nodes[i + k + r2]
            assert not dis.connected(a, b)
            assert dis.merge(a, b)
            assert dis.connected(a, b)

        assert nodes == list(dis)
        roots = [dis[i] for i in nodes]
        expected_indices = np.arange(n) - np.arange(n) % (2 * k)
        expected = [nodes[i] for i in expected_indices]
        assert roots == expected
