# Author: Eric Larson
# 2014

import numpy as np
cimport numpy as np
cimport cython

np.import_array()

# Fast inner loop of max_len_seq.
@cython.cdivision(True)  # faster modulo
@cython.boundscheck(False)  # designed to stay within bounds
@cython.wraparound(False)  # we don't use negative indexing
def _max_len_seq_inner(Py_ssize_t[::1] taps,
                       np.int8_t[::1] state,
                       Py_ssize_t nbits, Py_ssize_t length,
                       np.int8_t[::1] seq):
    # Here we compute MLS using a shift register, indexed using a ring buffer
    # technique (faster than using something like np.roll to shift)
    cdef Py_ssize_t n_taps = taps.shape[0]
    cdef Py_ssize_t idx = 0
    cdef np.int8_t feedback
    cdef Py_ssize_t i
    for i in range(length):
        feedback = state[idx]
        seq[i] = feedback
        for ti in range(n_taps):
            feedback ^= state[(taps[ti] + idx) % nbits]
        state[idx] = feedback
        idx = (idx + 1) % nbits
    # state must be rolled s.t. next run, when idx==0, it's in the right place
    return np.roll(state, -idx, axis=0)
