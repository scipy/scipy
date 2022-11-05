cimport numpy as np
import numpy as np

np.import_array()

from numpy.random cimport bitgen_t
from cpython.pycapsule cimport PyCapsule_GetPointer, PyCapsule_IsValid


cdef extern from "./rcont.c":
    int rcont1(double*, int, const double*, int, const double*,
               int**, bitgen_t*)
    int rcont2(double*, int, const double*, int, const double*,
               double* ntot, bitgen_t*)


cdef bitgen_t* get_bitgen(random_state):
    if isinstance(random_state, np.random.RandomState):
        bg = random_state._bit_generator
    elif isinstance(random_state, np.random.Generator):
        bg = random_state.bit_generator
    else:
        raise ValueError('random_state is not one of None, int, RandomState, Generator')
    capsule = bg.capsule

    cdef:
        const char *capsule_name = "BitGenerator"

    if not PyCapsule_IsValid(capsule, capsule_name):
        raise ValueError("Invalid pointer to anon_func_state.")

    return <bitgen_t *> PyCapsule_GetPointer(capsule, capsule_name)


cdef rvs_rcont1(double[:] row, double[:] col, int size, double ntot, random_state):

    cdef:
        bitgen_t *rstate
        int nr = row.shape[0]
        int nc = col.shape[0]
        int** work = NULL

    rstate = get_bitgen(random_state)

    result = np.zeros((size, nr, nc), dtype=np.double)

    cdef double [:,:,:] result_view = result

    cdef int error = 0
    for i in range(size):
        error = rcont1(<double*>result_view[i, :, :],
                       nr, <const double*>row,
                       nc, <const double*>col,
                       work, rstate)
        if error != 0:
            break
    
    return error, result