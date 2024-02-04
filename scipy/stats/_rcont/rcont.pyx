cimport numpy as np
import numpy as np

np.import_array()

from numpy.random cimport bitgen_t
from cpython.pycapsule cimport PyCapsule_GetPointer, PyCapsule_IsValid

ctypedef np.int64_t tab_t

cdef extern from "./_rcont.h":

    void rcont1_init(tab_t*, int, const tab_t*)
    void rcont1(tab_t*, int, const tab_t*, int, const tab_t*,
                tab_t, tab_t*, bitgen_t*)
    void rcont2(tab_t*, int, const tab_t*, int, const tab_t*,
                tab_t, bitgen_t*)


cdef bitgen_t* get_bitgen(random_state):
    if isinstance(random_state, np.random.RandomState):
        bg = random_state._bit_generator
    elif isinstance(random_state, np.random.Generator):
        bg = random_state.bit_generator
    else:
        raise ValueError('random_state is not RandomState or Generator')
    capsule = bg.capsule

    cdef:
        const char *capsule_name = "BitGenerator"

    if not PyCapsule_IsValid(capsule, capsule_name):
        raise ValueError("invalid pointer to anon_func_state")

    return <bitgen_t *> PyCapsule_GetPointer(capsule, capsule_name)


def rvs_rcont1(tab_t[::1] row, tab_t[::1] col, tab_t ntot,
               int size, random_state):

    cdef:
        bitgen_t *rstate = get_bitgen(random_state)
        int nr = row.shape[0]
        int nc = col.shape[0]

    cdef np.ndarray[tab_t, ndim=3, mode="c"] result = np.zeros(
        (size, nr, nc), dtype=np.int64
    )

    cdef np.ndarray[tab_t, ndim=1, mode="c"] work = np.empty(
        ntot, dtype=np.int64
    )

    if nc == 0 or nr == 0 or ntot == 0:
        return result

    rcont1_init(&work[0], nc, &col[0])

    for i in range(size):
        rcont1(&result[i, 0, 0], nr, &row[0], nc, &col[0], ntot,
               &work[0], rstate)

    return result


def rvs_rcont2(tab_t[::1] row, tab_t[::1] col, tab_t ntot,
               int size, random_state):
    cdef:
        bitgen_t *rstate = get_bitgen(random_state)
        int nr = row.shape[0]
        int nc = col.shape[0]

    cdef np.ndarray[tab_t, ndim=3, mode="c"] result = np.zeros(
        (size, nr, nc), dtype=np.int64
    )

    if nc == 0 or nr == 0 or ntot == 0:
        return result

    for i in range(size):
        rcont2(&result[i, 0, 0], nr, &row[0], nc, &col[0], ntot,
               rstate)

    return result
