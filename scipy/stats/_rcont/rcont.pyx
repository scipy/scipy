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


cdef get_bitgen(random_state):
    if isinstance(random_state, np.random.RandomState):
        bg = random_state._bit_generator
    elif isinstance(random_state, np.random.Generator):
        bg = random_state.bit_generator
    else:
        raise ValueError('random_state is not one of None, int, RandomState, Generator')
    capsule = bg.capsule

    if not PyCapsule_IsValid(capsule, capsule_name):
        raise ValueError("Invalid pointer to anon_func_state.")

    cdef:
        bitgen_t *bitgen
        const char *capsule_name = "BitGenerator"

    bitgen = <bitgen_t *> PyCapsule_GetPointer(capsule, capsule_name)

    return bitgen


cdef rvs_rcont1(row, col, int size, random_state):

    cdef:
        bitgen_t *bitgen
        Py_ssize_t nr = row.shape[0]
        Py_ssize_t nc = col.shape[0]

    bitgen = get_bitgen(random_state)

    result = np.zeros((size, nr, nc), dtype=np.int64)

    cdef np.intc 

    for i in range(size):
        