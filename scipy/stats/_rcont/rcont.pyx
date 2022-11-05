cimport numpy as np
import numpy as np

np.import_array()

from numpy.random cimport bitgen_t
from cpython.pycapsule cimport PyCapsule_GetPointer, PyCapsule_IsValid


cdef extern from "./_rcont.h":
    void rcont1_init(int*, int, const double*)
    void rcont1(double*, int, const double*, int, const double*,
                double, int*, bitgen_t*)
    void rcont2(double*, int, const double*, int, const double*,
                double, bitgen_t*)


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


def rvs_rcont1(double[::1] row, double[::1] col, double ntot,
               int size, random_state):

    cdef:
        bitgen_t *rstate = get_bitgen(random_state)
        int nr = row.shape[0]
        int nc = col.shape[0]

    cdef np.ndarray[double, ndim=3, mode="c"] result = np.zeros(
        (size, nr, nc), dtype=np.double
    )

    cdef np.ndarray[int, ndim=1, mode="c"] work = np.empty(
        <int>ntot, dtype=np.intc
    )

    rcont1_init(&work[0], nc, &col[0])

    for i in range(size):
        rcont1(&result[i, 0, 0], nr, &row[0], nc, &col[0], ntot,
               &work[0], rstate)

    return result


def rvs_rcont2(double[::1] row, double[::1] col, double ntot,
               int size, random_state):
    cdef:
        bitgen_t *rstate = get_bitgen(random_state)
        int nr = row.shape[0]
        int nc = col.shape[0]

    cdef np.ndarray[double, ndim=3, mode="c"] result = np.zeros(
        (size, nr, nc), dtype=np.double
    )

    for i in range(size):
        rcont2(&result[i, 0, 0], nr, &row[0], nc, &col[0], ntot,
               rstate)

    return result
