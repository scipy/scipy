import tracemalloc

import numpy as np
from numpy.testing import assert_equal
import pytest
from scipy.sparse import coo_array, lil_array

pytestmark = pytest.mark.thread_unsafe


def _assert_rhs_not_densified(A, key, rhs_sp, dense_nbytes):
    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        A[key] = rhs_sp
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < dense_nbytes // 4, (
        f"Sparse assignment to A[{key[0]}, {key[1]}]\n"
        f"  allocated {peak} bytes (>= {dense_nbytes // 4} bytes)\n"
        f"  dense size={dense_nbytes} bytes"
    )


_N_1D = 180_000
_DENSE_1D_NBYTES = _N_1D * 8
_IDX_1D = np.arange(100, 100 + _N_1D, dtype=np.intp)
_RHS_1D_ROW = coo_array(([7.0, 8.0], ([0, 0], [10, 50_000])), shape=(1, _N_1D))
_RHS_1D_COL = coo_array(([7.0, 8.0], ([10, 50_000], [0, 0])), shape=(_N_1D, 1))
_RHS_ZERO_SCALAR = coo_array((1, 1))


@pytest.mark.parametrize(
    ["is_wide", "key", "rhs_sp", "checks"],
    [
        # _set_intXslice_sparse
        (True, (2, slice(100, 100 + _N_1D)), _RHS_1D_ROW,
         [((2, 110), 7.0), ((2, 500), 0.0)]),
        (True, (-1, slice(100 + _N_1D, 100, -1)), _RHS_1D_ROW,
         [((9, 100 + _N_1D - 10), 7.0)]),
        (True, (2, slice(100, 100 + _N_1D)), _RHS_ZERO_SCALAR,
         [((2, 500), 0.0), ((2, 100_000), 0.0)]),
        # _set_sliceXint_sparse
        (False, (slice(100, 100 + _N_1D), 3), _RHS_1D_COL,
         [((110, 3), 7.0), ((500, 3), 0.0)]),
        (False, (slice(100 + _N_1D, 100, -1), -2), _RHS_1D_COL,
         [((100 + _N_1D - 10, 8), 7.0)]),
        (False, (slice(100, 100 + _N_1D), 3), _RHS_ZERO_SCALAR,
         [((500, 3), 0.0), ((100_000, 3), 0.0)]),
        # _set_intXarray_sparse
        (True, (2, _IDX_1D), _RHS_1D_ROW,
         [((2, 110), 7.0), ((2, 500), 0.0)]),
        (True, (2, _IDX_1D), _RHS_ZERO_SCALAR,
         [((2, 500), 0.0), ((2, 100_000), 0.0)]),
        # _set_arrayXint_sparse
        (False, (_IDX_1D, 3), _RHS_1D_COL,
         [((110, 3), 7.0), ((500, 3), 0.0)]),
        (False, (_IDX_1D, 3), _RHS_1D_ROW,
         [((110, 3), 7.0), ((500, 3), 0.0)]),
        (False, (_IDX_1D, 3), _RHS_ZERO_SCALAR,
         [((500, 3), 0.0), ((100_000, 3), 0.0)]),
    ],
)
def test_1d_sparse_assignment_not_densified(is_wide, key, rhs_sp, checks):
    if is_wide:
        A = coo_array(([1.0, 2.0], ([2, 2], [500, 100_000])), shape=(10, 200_000))
    else:
        A = coo_array(([1.0, 2.0], ([500, 100_000], [3, 3])), shape=(200_000, 10))
    A = lil_array(A)

    _assert_rhs_not_densified(A, key, rhs_sp, _DENSE_1D_NBYTES)
    for (r, c), expected in checks:
        assert A[r, c] == expected


_DENSE_2D_NBYTES = 1000 * 1000 * 8
_IDX_2D = np.arange(100, 1100, dtype=np.intp)
_RHS_2D_MAP = {
    "full": coo_array(([3.0, 4.0], ([5, 500], [10, 600])), shape=(1000, 1000)),
    "broadcast_row": coo_array(([3.0, 4.0], ([0, 0], [10, 600])), shape=(1, 1000)),
    "broadcast_col": coo_array(([3.0, 4.0], ([5, 500], [0, 0])), shape=(1000, 1)),
}


@pytest.mark.parametrize("rhs_kind", ["full", "broadcast_row", "broadcast_col"])
@pytest.mark.parametrize(
    "key",
    [
        # _set_sliceXslice_sparse (contiguous, positive step, negative step)
        (slice(100, 1100), slice(100, 1100)),
        (slice(100, 2100, 2), slice(100, 2100, 2)),
        (slice(2100, 100, -2), slice(2100, 100, -2)),
        # _set_arrayXslice_sparse (contiguous, negative step)
        (_IDX_2D, slice(100, 1100)),
        (_IDX_2D, slice(1100, 100, -1)),
        # _set_sliceXarray_sparse (contiguous, negative step)
        (slice(100, 1100), _IDX_2D),
        (slice(1100, 100, -1), _IDX_2D),
        # _set_columnXarray_sparse (outer and transposed outer indexing)
        (_IDX_2D[:, None], _IDX_2D),
        (_IDX_2D[:, None], _IDX_2D[None, :]),
        (_IDX_2D[None, :], _IDX_2D[:, None]),
    ],
)
def test_2d_outer_sparse_assignment_not_densified(key, rhs_kind):
    A = lil_array(coo_array(([1., 2.], ([150, 800], [200, 900])), shape=(2200, 2200)))
    _assert_rhs_not_densified(A, key, _RHS_2D_MAP[rhs_kind], _DENSE_2D_NBYTES)


_DENSE_INNER_NBYTES = 400 * 500 * 8
_R_INNER = np.arange(200_000, dtype=np.intp).reshape(400, 500) % 1400
_C_INNER = (np.arange(200_000, dtype=np.intp).reshape(400, 500) * 7) % 1400
_RHS_INNER_MAP = {
    "full": coo_array(([3.0, 4.0], ([5, 200], [10, 300])), shape=(400, 500)),
    "broadcast_row": coo_array(([3.0, 4.0], ([0, 0], [10, 300])), shape=(1, 500)),
    "broadcast_col": coo_array(([3.0, 4.0], ([5, 200], [0, 0])), shape=(400, 1)),
}


@pytest.mark.parametrize("rhs_kind", ["full", "broadcast_row", "broadcast_col"])
def test_2d_inner_sparse_assignment_not_densified(rhs_kind):
    A = lil_array(coo_array(([1., 2.], ([150, 800], [200, 900])), shape=(1500, 1500)))
    _assert_rhs_not_densified(
        A, (_R_INNER, _C_INNER), _RHS_INNER_MAP[rhs_kind], _DENSE_INNER_NBYTES
    )


@pytest.mark.parametrize(
    "key",
    [
        (slice(1, 3), slice(1, 4)),
        (np.array([1, 3]), slice(1, 4)),
        (slice(1, 3), np.array([1, 2, 4])),
        (np.array([[1], [3]]), np.array([1, 2, 4])),
    ],
)
def test_sparse_assignment_explicit_zeros(key):
    B = np.ones((6, 6), dtype=np.float32)
    A = lil_array(B)
    rhs_explicit_zero = coo_array(
        (np.array([0.0, 9.0], dtype=np.float32), ([0, 0], [0, 1])),
        shape=(2, 3),
        dtype=np.float32,
    )

    A[key] = rhs_explicit_zero
    B[key] = rhs_explicit_zero.toarray()
    assert_equal(A.toarray(), B)

