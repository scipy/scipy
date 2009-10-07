"""
Cython rewrite of the vector quantization module, originally written
in C at src/vq.c and the wrapper at src/vq_module.c. This should be
easier to maintain than old SWIG output.

Original C version by Damian Eads. 
Translated to Cython by David Warde-Farley, October 2009.
"""

import numpy as np
cimport numpy as np

cimport cython

cdef extern from "math.h":
    float sqrtf(float num)
    double sqrt(double num)

#cdef double rbig = 1e100

# Python/NumPy types
FLOAT64 = np.float64
FLOAT32 = np.float32
INT32 = np.int32

# C types
ctypedef np.float64_t FLOAT64_t
ctypedef np.float32_t FLOAT32_t
ctypedef np.int32_t INT32_t

@cython.boundscheck(False)
@cython.wraparound(False)
def float_tvq(np.ndarray[FLOAT32_t, ndim=2] obs, 
                   np.ndarray[FLOAT32_t, ndim=2] code_book,
                   np.ndarray[INT32_t, ndim=1] codes,
                   np.ndarray[FLOAT32_t, ndim=1] low_dist):
    """
    Quantize Nobs observations to their nearest codebook 
    entry (single-precision version).
    """
    # Temporary variables
    cdef float dist, diff
    cdef int obs_index, code_index, feature
    
    # Loop limits
    cdef int ncodes = code_book.shape[0]
    cdef int nfeat = code_book.shape[1]
    cdef int nobs = obs.shape[0]
    
    for obs_index in range(nobs):
        
        low_dist[obs_index] = np.inf

        for code_index in range(ncodes):
            dist = 0
            
            # Distance between code_book[code_index] and obs[obs_index]
            for feature in range(nfeat):
                diff = code_book[code_index, feature] - obs[obs_index, feature]
                dist += diff * diff
            dist = sqrtf(dist)
            
            # Replace the code assignment and record distance if necessary
            if dist < low_dist[obs_index]:
                codes[obs_index] = code_index
                low_dist[obs_index] = dist

def vq(np.ndarray obs, np.ndarray codes):
    """
    Testing    
    """
    cdef np.npy_intp n_obs, n_codes, n_features
    cdef np.ndarray outcodes, outdists
    obs = np.atleast_2d(np.ascontiguousarray(obs))
    codes = np.atleast_2d(np.ascontiguousarray(codes))

    if obs.ndim > 2 or codes.ndim > 2:
        raise ValueError("rank > 2 for arguments not supported")
    
    if obs.dtype != codes.dtype:
        raise ValueError("obs and codes must be of same dtype")
    
    if obs.shape[1] != codes.shape[1]:
        raise ValueError("obs and codes should have same number of " + \
                         "features (columns)")

    nobs = obs.shape[0]
    outcodes = np.zeros((nobs,), dtype=np.int32)
    outdists = np.zeros((nobs,), dtype=obs.dtype)
    
    if obs.dtype == np.float32:
        float_tvq(obs, codes, outcodes, outdists)

    return outcodes, outdists
