"""
Cython rewrite of the vector quantization module, originally written
in C at src/vq.c and the wrapper at src/vq_module.c. This should be
easier to maintain than old SWIG output.

Original C version by Damian Eads.
Translated to Cython by David Warde-Farley, October 2009.
"""

import numpy as np
cimport numpy as np
from cluster_blas cimport *

cdef extern from "math.h":
    float sqrtf(float num)
    double sqrt(double num)

ctypedef np.float64_t float64_t
ctypedef np.float32_t float32_t
ctypedef np.int32_t int32_t

# Use Cython fused types for templating
# Define supported data types as vq_type
ctypedef fused vq_type:
    float32_t
    float64_t

# When the number of features is less than this number,
# switch back to the naive algorithm to avoid high overhead.
DEF NFEATURES_CUTOFF=5

# Initialize the NumPy C API
np.import_array()


cdef inline vq_type _sqrt(vq_type x):
    if vq_type is float32_t:
        return sqrtf(x)
    else:
        return sqrt(x)


cdef inline vq_type vec_sqr(int n, vq_type *p):
    cdef vq_type result = 0.0
    cdef int i
    for i in range(n):
        result += p[i] * p[i]
    return result


cdef inline void cal_M(int nobs, int ncodes, int nfeat, vq_type *obs,
                       vq_type *code_book, vq_type *M):
    """
    Calculate M = obs * code_book.T
    """
    cdef vq_type alpha = -2.0, beta = 0.0

    # Call BLAS functions with Fortran ABI
    # Note that BLAS Fortran ABI uses column-major order
    if vq_type is float32_t:
        f_sgemm("T", "N", &ncodes, &nobs, &nfeat,
                &alpha, code_book, &nfeat, obs, &nfeat, &beta, M, &ncodes)
    else:
        f_dgemm("T", "N", &ncodes, &nobs, &nfeat,
                &alpha, code_book, &nfeat, obs, &nfeat, &beta, M, &ncodes)


cdef void _vq(vq_type *obs, vq_type *code_book,
              int ncodes, int nfeat, int nobs,
              int32_t *codes, vq_type *low_dist):
    # Naive algorithm is prefered when nfeat is small
    if nfeat < NFEATURES_CUTOFF:
        _vq_small_nf(obs, code_book, ncodes, nfeat, nobs, codes, low_dist)
        return

    cdef int obs_index, code_index
    cdef vq_type *p_obs
    cdef vq_type *p_codes
    cdef vq_type dist_sqr
    cdef np.ndarray[vq_type, ndim=1] obs_sqr, codes_sqr
    cdef np.ndarray[vq_type, ndim=2] M

    if vq_type is float32_t:
        obs_sqr = np.ndarray(nobs, np.float32)
        codes_sqr = np.ndarray(ncodes, np.float32)
        M = np.ndarray((nobs, ncodes), np.float32)
    else:
        obs_sqr = np.ndarray(nobs, np.float64)
        codes_sqr = np.ndarray(ncodes, np.float64)
        M = np.ndarray((nobs, ncodes), np.float64)

    p_obs = obs
    for obs_index in range(nobs):
        # obs_sqr[i] is the inner product of the i-th observation with itself
        obs_sqr[obs_index] = vec_sqr(nfeat, p_obs)
        p_obs += nfeat

    p_codes = code_book
    for code_index in range(ncodes):
        # codes_sqr[i] is the inner product of the i-th code with itself
        codes_sqr[code_index] = vec_sqr(nfeat, p_codes)
        p_codes += nfeat

    # M[i][j] is the inner product of the i-th obs and j-th code
    # M = obs * codes.T
    cal_M(nobs, ncodes, nfeat, obs, code_book, <vq_type *>M.data)

    for obs_index in range(nobs):
        for code_index in range(ncodes):
            dist_sqr = (M[obs_index, code_index] +
                    obs_sqr[obs_index] + codes_sqr[code_index])
            if dist_sqr < low_dist[obs_index]:
                codes[obs_index] = code_index
                low_dist[obs_index] = dist_sqr

        # dist_sqr may be negative due to float point errors
        if low_dist[obs_index] > 0:
            low_dist[obs_index] = _sqrt(low_dist[obs_index])
        else:
            low_dist[obs_index] = 0


cdef void _vq_small_nf(vq_type *obs, vq_type *code_book,
                       int ncodes, int nfeat, int nobs,
                       int32_t *codes, vq_type *low_dist):
    """
    Vector quantization for float32 using naive algorithm.
    This is prefered when nfeat is small.
    """
    # Temporary variables
    cdef vq_type dist_sqr, diff
    cdef int32_t obs_index, code_index, feature
    cdef int32_t offset = 0

    # Index and pointer to keep track of the current position in
    # both arrays so that we don't have to always do index * nfeat.
    cdef int codebook_pos
    cdef vq_type *current_obs

    for obs_index in range(nobs):
        codebook_pos = 0
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

        low_dist[obs_index] = _sqrt(low_dist[obs_index])

        # Update the offset of the current observation
        offset += nfeat


def vq(np.ndarray obs, np.ndarray codes):
    """
    Vector quantization ndarray wrapper. Only support float32 and float64.

    Parameters
    ----------
    obs : ndarray
        The observation matrix. Each row is an observation.
    codes : ndarray
        The code book matrix.

    Notes
    -----
    The observation matrix and code book matrix should have same ndim and
    same number of columns (features). Only 1-dimensional and 2-dimensional
    arrays are supported.
    """
    cdef int nobs, ncodes, nfeat
    cdef np.ndarray obs_a, codes_a
    cdef np.ndarray outcodes, outdists
    cdef int flags = np.NPY_CONTIGUOUS | np.NPY_NOTSWAPPED | np.NPY_ALIGNED

    # Ensure the arrays are contiguous - done in C to minimize overhead.
    obs_a = np.PyArray_FROM_OF(obs, flags)
    codes_a = np.PyArray_FROM_OF(codes, flags)

    if obs.dtype != codes.dtype:
        raise TypeError('observation and code should have same dtype')
    if obs.dtype not in (np.float32, np.float64):
        raise TypeError('type other than float or double not supported')
    if obs_a.ndim != codes_a.ndim:
        raise ValueError(
            'observation and code should have same number of dimensions')

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
        raise ValueError('ndim different than 1 or 2 are not supported')

    # Initialize outdists and outcodes array.
    # Outdists should be initialized as INF.
    outdists = np.empty((nobs,), dtype=obs.dtype)
    outcodes = np.empty((nobs,), dtype=np.int32)
    outdists.fill(np.inf)

    if obs.dtype.type is np.float32:
        _vq(<float32_t *>obs_a.data, <float32_t *>codes_a.data,
            ncodes, nfeat, nobs, <int32_t *>outcodes.data,
            <float32_t *>outdists.data)
    elif obs.dtype.type is np.float64:
        _vq(<float64_t *>obs_a.data, <float64_t *>codes_a.data,
            ncodes, nfeat, nobs, <int32_t *>outcodes.data,
            <float64_t *>outdists.data)

    return outcodes, outdists

