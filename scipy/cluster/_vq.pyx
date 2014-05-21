"""
Cython rewrite of the vector quantization module, originally written
in C at src/vq.c and the wrapper at src/vq_module.c. This should be
easier to maintain than old SWIG output.

Original C version by Damian Eads.
Translated to Cython by David Warde-Farley, October 2009.
"""

import numpy as np
cimport numpy as np
from cblas cimport *

cdef extern from "math.h":
    float sqrtf(float num)
    double sqrt(double num)

cdef extern from "numpy/arrayobject.h":
    object PyArray_EMPTY(int, np.npy_intp*, int, int)
    cdef enum:
        PyArray_INTP

cdef extern from "numpy/npy_math.h":
    cdef enum:
        NPY_INFINITY

# C types
ctypedef np.float64_t float64_t
ctypedef np.float32_t float32_t
ctypedef np.int32_t int32_t

# Initialize the NumPy C API
np.import_array()


cdef void float32_tvq(float32_t *obs, float32_t *code_book,
                      int ncodes, int nfeat, int nobs,
                      np.npy_intp *codes, float32_t *low_dist):
    # Naive algorithm is prefered when nfeat is small
    if nfeat < 5:
        float32_tvq_small_nf(obs, code_book, ncodes, nfeat, nobs,
                             codes, low_dist)
        return

    cdef int obs_index, code_index
    cdef float32_t *p_obs, *p_codes, dist_sqr
    cdef np.ndarray[float32_t, ndim=1] obs_sqr, codes_sqr
    cdef np.ndarray[float32_t, ndim=2] M

    obs_sqr = np.ndarray(nobs, np.float32)
    codes_sqr = np.ndarray(ncodes, np.float32)
    M = np.ndarray((nobs, ncodes), np.float32)

    p_obs = obs
    for obs_index in range(nobs):
        # obs_sqr[i] is the inner product of the i-th observation with itself
        obs_sqr[obs_index] = cblas_sdot(nfeat, p_obs, 1, p_obs, 1)
        p_obs += nfeat

    p_codes = code_book
    for code_index in range(ncodes):
        # codes_sqr[i] is the inner product of the i-th code with itself
        codes_sqr[code_index] = cblas_sdot(nfeat, p_codes, 1, p_codes, 1)
        p_codes += nfeat

    # M[i][j] is the inner product of the i-th obs and j-th code
    # M = Obs * Codes.T
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                nobs, ncodes, nfeat, -2.0, obs, nfeat, code_book, nfeat,
                0.0, <float32_t *>M.data, ncodes)

    for obs_index in range(nobs):
        low_dist[obs_index] = NPY_INFINITY
        for code_index in range(ncodes):
            dist_sqr = (M[obs_index, code_index] +
                    obs_sqr[obs_index] + codes_sqr[code_index])
            if dist_sqr < low_dist[obs_index]:
                codes[obs_index] = code_index
                low_dist[obs_index] = dist_sqr

        # dist_sqr may be negative due to float point errors
        if low_dist[obs_index] > 0:
            low_dist[obs_index] = sqrtf(low_dist[obs_index])
        else:
            low_dist[obs_index] = 0


cdef void float32_tvq_small_nf(float32_t *obs, float32_t *code_book,
                               int ncodes, int nfeat, int nobs,
                               np.npy_intp *codes, float32_t *low_dist):
    """
    Vector quantization for float32 using naive algorithm.
    This is perfered when nfeat is small.
    """
    # Temporary variables
    cdef float32_t dist_sqr, diff
    cdef int32_t obs_index, code_index, feature
    cdef int32_t offset = 0

    # Index and pointer to keep track of the current position in
    # both arrays so that we don't have to always do index * nfeat.
    cdef int codebook_pos
    cdef float32_t *current_obs

    for obs_index in range(nobs):
        codebook_pos = 0
        low_dist[obs_index] = NPY_INFINITY
        for code_index in range(ncodes):
            dist_sqr = 0

            # Distance between code_book[code_index] and obs[obs_index]
            for feature in range(nfeat):

                # Use current_obs pointer and codebook_pos to minimize
                # pointer arithmetic necessary (i.e. no multiplications)
                current_obs = &(obs[offset])
                diff = code_book[codebook_pos] - current_obs[feature]
                dist_sqr += diff * diff
                codebook_pos += 1

            # Replace the code assignment and record distance if necessary
            if dist_sqr < low_dist[obs_index]:
                codes[obs_index] = code_index
                low_dist[obs_index] = dist_sqr

        low_dist[obs_index] = sqrtf(low_dist[obs_index])

        # Update the offset of the current observation
        offset += nfeat


cdef void float64_tvq(float64_t *obs, float64_t *code_book,
                      int ncodes, int nfeat, int nobs,
                      np.npy_intp *codes, float64_t *low_dist):
    # Naive algorithm is prefered when nfeat is small
    if nfeat < 5:
        float64_tvq_small_nf(obs, code_book, ncodes, nfeat, nobs,
                             codes, low_dist)
        return

    cdef int obs_index, code_index
    cdef float64_t *p_obs, *p_codes, dist_sqr
    cdef np.ndarray[float64_t, ndim=1] obs_sqr, codes_sqr
    cdef np.ndarray[float64_t, ndim=2] M

    obs_sqr = np.ndarray(nobs, np.float64)
    codes_sqr = np.ndarray(ncodes, np.float64)
    M = np.ndarray((nobs, ncodes), np.float64)

    p_obs = obs
    for obs_index in range(nobs):
        # obs_sqr[i] is the inner product of the i-th observation with itself
        obs_sqr[obs_index] = cblas_ddot(nfeat, p_obs, 1, p_obs, 1)
        p_obs += nfeat

    p_codes = code_book
    for code_index in range(ncodes):
        # codes_sqr[i] is the inner product of the i-th code with itself
        codes_sqr[code_index] = cblas_ddot(nfeat, p_codes, 1, p_codes, 1)
        p_codes += nfeat

    # M[i][j] is the inner product of the i-th obs and j-th code
    # M = Obs * Codes.T
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                nobs, ncodes, nfeat, -2.0, obs, nfeat, code_book, nfeat,
                0.0, <float64_t *>M.data, ncodes)

    for obs_index in range(nobs):
        low_dist[obs_index] = NPY_INFINITY
        for code_index in range(ncodes):
            dist_sqr = (M[obs_index, code_index] +
                    obs_sqr[obs_index] + codes_sqr[code_index])
            if dist_sqr < low_dist[obs_index]:
                codes[obs_index] = code_index
                low_dist[obs_index] = dist_sqr

        # dist_sqr may be negative due to float point errors
        if low_dist[obs_index] > 0:
            low_dist[obs_index] = sqrt(low_dist[obs_index])
        else:
            low_dist[obs_index] = 0


cdef void float64_tvq_small_nf(float64_t *obs, float64_t *code_book,
                               int ncodes, int nfeat, int nobs,
                               np.npy_intp *codes, float64_t *low_dist):
    """
    Vector quantization for float64 using naive algorithm.
    This is perfered when nfeat is small.
    """
    # Temporary variables
    cdef float64_t dist_sqr, diff
    cdef int32_t obs_index, code_index, feature
    cdef int32_t offset = 0

    # Index and pointer to keep track of the current position in
    # both arrays so that we don't have to always do index * nfeat.
    cdef int codebook_pos
    cdef float64_t *current_obs

    for obs_index in range(nobs):
        codebook_pos = 0
        low_dist[obs_index] = NPY_INFINITY
        for code_index in range(ncodes):
            dist_sqr = 0

            # Distance between code_book[code_index] and obs[obs_index]
            for feature in range(nfeat):

                # Use current_obs pointer and codebook_pos to minimize
                # pointer arithmetic necessary (i.e. no multiplications)
                current_obs = &(obs[offset])
                diff = code_book[codebook_pos] - current_obs[feature]
                dist_sqr += diff * diff
                codebook_pos += 1

            # Replace the code assignment and record distance if necessary
            if dist_sqr < low_dist[obs_index]:
                codes[obs_index] = code_index
                low_dist[obs_index] = dist_sqr

        low_dist[obs_index] = sqrt(low_dist[obs_index])

        # Update the offset of the current observation
        offset += nfeat


def vq(np.ndarray obs, np.ndarray codes):
    """
    vq(obs, codes)

    Vector quantization ndarray wrapper.
    """
    cdef np.npy_intp nobs, ncodes, nfeat
    cdef np.ndarray obs_a, codes_a
    cdef np.ndarray outcodes, outdists
    cdef int flags = np.NPY_CONTIGUOUS | np.NPY_NOTSWAPPED | np.NPY_ALIGNED

    # Ensure the arrays are contiguous - done in C to minimize overhead.
    obs_a = np.PyArray_FROM_OF(obs, flags)
    codes_a = np.PyArray_FROM_OF(codes, flags)

    if obs_a.ndim != codes_a.ndim:
        raise ValueError('observation and code should have same rank')

    if obs_a.ndim == 1:
        nfeat = 1
        nobs = obs_a.shape[0]
        ncodes = codes_a.shape[0]
    elif obs_a.ndim == 2:
        nfeat = obs_a.shape[1]
        nobs = obs_a.shape[0]
        ncodes = codes_a.shape[0]
        if nfeat != codes_a.shape[1]:
            raise ValueError('obs and code should have same number of '
                             'features (columns)')
    else:
        raise ValueError('rank different than 1 or 2 are not supported')

    # We create this with the C API so that we can be sure that
    # the resulting array has elements big enough to store indices
    # on that platform. Hence, PyArray_INTP.
    outcodes = PyArray_EMPTY(1, &nobs, PyArray_INTP, 0)

    # This we just want to match the dtype of the input, so np.empty is fine.
    outdists = np.empty((nobs,), dtype=obs.dtype)

    if obs.dtype == np.float32:
        float32_tvq(<float32_t *>obs_a.data, <float32_t *>codes_a.data,
                    ncodes, nfeat, nobs, <np.npy_intp *>outcodes.data,
                    <float32_t *>outdists.data)
    elif obs.dtype == np.float64:
        float64_tvq(<float64_t *>obs_a.data, <float64_t *>codes_a.data,
                    ncodes, nfeat, nobs, <np.npy_intp *>outcodes.data,
                    <float64_t *>outdists.data)
    else:
        raise ValueError('type other than float or double not supported')

    return outcodes, outdists

