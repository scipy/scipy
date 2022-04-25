import pytest
import numpy as np

from .. import _construct, array, coo_matrix, _arrays


@pytest.mark.parametrize(
    ("fn", "args"),
    [
        ("diags", ([0, 1, 2],)),
        ("eye", (3,)),
        ("spdiags", ([1, 2, 3], 0, 3, 3)),
        ("identity", (3,)),
        # kron with dense B
        ("kron", (np.array([[1, 2], [3, 4]]), np.array([[4, 3], [2, 1]]))),
        # kron with sparse B
        ("kron", (np.array([[1, 2], [3, 4]]), np.array([[1, 0], [0, 0]]))),
        ("kronsum", (np.array([[1, 0], [0, 1]]), np.array([[0, 1], [1, 0]]))),
        ("random", (3, 3)),
    ],
)
def test_default_container_format(fn, args):
    """Default container format for array creation functions should be the same
    as their matrix counterparts."""
    ary_fn, mat_fn = (getattr(mod, fn) for mod in (array, _construct))
    a = ary_fn(*args)
    assert isinstance(a, _arrays._sparray)
    assert a.format == mat_fn(*args).format


@pytest.mark.parametrize("fn", ("hstack", "vstack"))
def test_stacks_default_container_format(fn):
    """Same idea as `test_default_container_format`, but for the stacking
    creation functions."""
    A = coo_matrix(np.eye(2))
    B = coo_matrix([[0, 1], [1, 0]])

    ary_fn, mat_fn = (getattr(mod, fn) for mod in (array, _construct))
    a = ary_fn([A, B])
    assert isinstance(a, _arrays._sparray)
    assert a.format == mat_fn([A, B]).format


def test_blocks_default_container_format():
    """Same idea as `test_default_container_format`, but for the block
    block creation funtions."""
    A = coo_matrix(np.eye(2))
    B = coo_matrix([[2], [0]])
    C = coo_matrix([[3]])

    # block diag
    a = array.block_diag((A, B, C))
    m = _construct.block_diag((A, B, C))
    assert isinstance(a, _arrays._sparray)
    assert a.format == m.format

    # bmat
    a = array.bmat([[A, None], [None, C]])
    m = _construct.bmat([[A, None], [None, C]])
    assert isinstance(a, _arrays._sparray)
    assert a.format == m.format
