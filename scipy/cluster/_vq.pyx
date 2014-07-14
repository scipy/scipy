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
    """
    The underlying function (template) of _vq.vq.

    Parameters
    ----------
    obs : vq_type*
        The pointer to the observation matrix.
    code_book : vq_type*
        The pointer to the code book matrix.
    ncodes : int
        The number of centroids (codes).
    nfeat : int
        The number of features of each observation.
    nobs : int
        The number of observations.
    codes : vq_type*
        The pointer to the new codes array.
    low_dist : vq_type*
        low_dist[i] is the Euclidean distance from obs[i] to the corresponding
        centroid.
    """
    # Naive algorithm is prefered when nfeat is small
    if nfeat < NFEATURES_CUTOFF:
        _vq_small_nf(obs, code_book, ncodes, nfeat, nobs, codes, low_dist)
        return

    cdef np.npy_intp i, j
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
    for i in range(nobs):
        # obs_sqr[i] is the inner product of the i-th observation with itself
        obs_sqr[i] = vec_sqr(nfeat, p_obs)
        p_obs += nfeat

    p_codes = code_book
    for i in range(ncodes):
        # codes_sqr[i] is the inner product of the i-th code with itself
        codes_sqr[i] = vec_sqr(nfeat, p_codes)
        p_codes += nfeat

    # M[i][j] is the inner product of the i-th obs and j-th code
    # M = obs * codes.T
    cal_M(nobs, ncodes, nfeat, obs, code_book, <vq_type *>M.data)

    for i in range(nobs):
        for j in range(ncodes):
            dist_sqr = (M[i, j] +
                    obs_sqr[i] + codes_sqr[j])
            if dist_sqr < low_dist[i]:
                codes[i] = j
                low_dist[i] = dist_sqr

        # dist_sqr may be negative due to float point errors
        if low_dist[i] > 0:
            low_dist[i] = _sqrt(low_dist[i])
        else:
            low_dist[i] = 0


cdef void _vq_small_nf(vq_type *obs, vq_type *code_book,
                       int ncodes, int nfeat, int nobs,
                       int32_t *codes, vq_type *low_dist):
    """
    Vector quantization using naive algorithm.
    This is prefered when nfeat is small.
    The parameters are the same as those of _vq.
    """
    # Temporary variables
    cdef vq_type dist_sqr, diff
    cdef np.npy_intp i, j, k, obs_offset = 0, code_offset

    # Index and pointer to keep track of the current position in
    # both arrays so that we don't have to always do index * nfeat.
    cdef vq_type *current_obs
    cdef vq_type *current_code

    for i in range(nobs):
        code_offset = 0
        current_obs = &(obs[obs_offset])

        for j in range(ncodes):
            dist_sqr = 0
            current_code = &(code_book[code_offset])

            # Distance between code_book[j] and obs[i]
            for k in range(nfeat):
                diff = current_code[k] - current_obs[k]
                dist_sqr += diff * diff
            code_offset += nfeat

            # Replace the code assignment and record distance if necessary
            if dist_sqr < low_dist[i]:
                codes[i] = j
                low_dist[i] = dist_sqr

        low_dist[i] = _sqrt(low_dist[i])

        # Update the offset of the current observation
        obs_offset += nfeat


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
    cdef np.ndarray outcodes, outdists

    # Ensure the arrays are contiguous
    obs = np.ascontiguousarray(obs)
    codes = np.ascontiguousarray(codes)

    if obs.dtype != codes.dtype:
        raise TypeError('observation and code should have same dtype')
    if obs.dtype not in (np.float32, np.float64):
        raise TypeError('type other than float or double not supported')
    if obs.ndim != codes.ndim:
        raise ValueError(
            'observation and code should have same number of dimensions')

    if obs.ndim == 1:
        nfeat = 1
        nobs = obs.shape[0]
        ncodes = codes.shape[0]
    elif obs.ndim == 2:
        nfeat = obs.shape[1]
        nobs = obs.shape[0]
        ncodes = codes.shape[0]
        if nfeat != codes.shape[1]:
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
        _vq(<float32_t *>obs.data, <float32_t *>codes.data,
            ncodes, nfeat, nobs, <int32_t *>outcodes.data,
            <float32_t *>outdists.data)
    elif obs.dtype.type is np.float64:
        _vq(<float64_t *>obs.data, <float64_t *>codes.data,
            ncodes, nfeat, nobs, <int32_t *>outcodes.data,
            <float64_t *>outdists.data)

    return outcodes, outdists


cdef np.ndarray _update_cluster_means(vq_type *obs, int32_t *labels,
                                      vq_type *cb, int nobs, int nc, int nfeat):
    """
    The underlying function (template) of _vq.update_cluster_means.

    Parameters
    ----------
    obs : vq_type*
        The pointer to the observation matrix.
    labels : int32_t*
        The pointer to the array of the labels (codes) of the observations.
    cb : vq_type*
        The pointer to the new code book matrix.
    nobs : int
        The number of observations.
    nc : int
        The number of centroids (codes).
    nfeat : int
        The number of features of each observation.

    Returns
    -------
    has_members : ndarray
        A boolean array indicating which clusters have members.
    """
    cdef np.npy_intp i, j, cluster_size, label
    cdef vq_type *obs_p
    cdef vq_type *cb_p
    cdef np.ndarray[int, ndim=1] obs_count

    # Calculate the sums the numbers of obs in each cluster
    obs_count = np.zeros(nc, np.int32)
    obs_p = obs
    for i in range(nobs):
        label = labels[i]
        cb_p = cb + nfeat * label

        for j in range(nfeat):
            cb_p[j] += obs_p[j]

        # Count the obs in each cluster
        obs_count[label] += 1
        obs_p += nfeat

    cb_p = cb
    for i in range(nc):
        cluster_size = obs_count[i]

        if cluster_size > 0:
            # Calculate the centroid of each cluster
            for j in range(nfeat):
                cb_p[j] /= cluster_size

        cb_p += nfeat

    # Return a boolean array indicating which clusters have members
    return obs_count > 0


def update_cluster_means(np.ndarray obs, np.ndarray labels, int nc):
    """
    The update-step of K-means. Calculate the mean of observations in each
    cluster.

    Parameters
    ----------
    obs : ndarray
        The observation matrix. Each row is an observation. Its dtype must be
        float32 or float64.
    labels : ndarray
        The label of each observation. Must be an 1d array.
    nc : int
        The number of centroids.

    Returns
    -------
    cb : ndarray
        The new code book.
    has_members : ndarray
        A boolean array indicating which clusters have members.

    Notes
    -----
    The empty clusters will be set to all zeros and the curresponding elements
    in `has_members` will be `False`. The upper level function should decide
    how to deal with them.
    """
    cdef np.ndarray has_members, cb
    cdef int nfeat

    # Ensure the arrays are contiguous
    obs = np.ascontiguousarray(obs)
    labels = np.ascontiguousarray(labels)

    if obs.dtype not in (np.float32, np.float64):
        raise TypeError('type other than float or double not supported')
    if labels.dtype.type is not np.int32:
        labels = labels.astype(np.int32)
    if labels.ndim != 1:
        raise ValueError('labels must be an 1d array')

    if obs.ndim == 1:
        nfeat = 1
        cb = np.zeros(nc, dtype=obs.dtype)
    elif obs.ndim == 2:
        nfeat = obs.shape[1]
        cb = np.zeros((nc, nfeat), dtype=obs.dtype)
    else:
        raise ValueError('ndim different than 1 or 2 are not supported')

    if obs.dtype.type is np.float32:
        has_members = _update_cluster_means(<float32_t *>obs.data,
                                            <int32_t *>labels.data,
                                            <float32_t *>cb.data,
                                            obs.shape[0], nc, nfeat)
    elif obs.dtype.type is np.float64:
        has_members = _update_cluster_means(<float64_t *>obs.data,
                                            <int32_t *>labels.data,
                                            <float64_t *>cb.data,
                                            obs.shape[0], nc, nfeat)

    return cb, has_members

