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


cdef extern from "Python.h":
    ctypedef int Py_intptr_t

cdef extern from "math.h":
    float sqrtf(float num)
    double sqrt(double num)

cdef extern from "numpy/arrayobject.h":
    ctypedef Py_intptr_t npy_intp
    cdef npy_intp PyArray_TYPE(object arr)
    cdef npy_intp PyArray_DIM(object arr, int n)
    cdef npy_intp PyArray_NDIM(object arr)

cdef extern from "numpy/npy_math.h":
    cdef enum:
        NPY_INFINITY
    
# Python/NumPy types
float64 = np.float64
float32 = np.float32
int32 = np.int32

# C types
ctypedef np.float64_t float64_t
ctypedef np.float32_t float32_t
ctypedef np.int32_t int32_t

cdef int float_tvq(float32_t *obs, float32_t *code_book,
                   int ncodes, int nfeat, int nobs,
                   int32_t *codes, float32_t *low_dist):
    """
    Quantize Nobs observations to their nearest codebook 
    entry (single-precision version).
    """
    # Temporary variables
    cdef float32_t dist, diff
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
            dist = 0
            # Distance between code_book[code_index] and obs[obs_index]
            for feature in range(nfeat):
                # Use current_obs pointer and codebook_pos to minimize
                # pointer arithmetic necessary (i.e. no multiplications)
                current_obs = &(obs[offset])
                diff = code_book[codebook_pos] - current_obs[feature]
                dist += diff * diff
                codebook_pos += 1
            dist = sqrtf(dist)
            
            # Replace the code assignment and record distance if necessary
            if dist < low_dist[obs_index]:
                codes[obs_index] = code_index
                low_dist[obs_index] = dist
        
        # Update the offset of the current observation
        offset += nfeat
    return 0

def vq(np.ndarray obs, np.ndarray codes):
    """
    Vector quantization ndarray wrapper.
    """
    cdef np.npy_intp nobs, ncodes, nfeat, nfeat_codes
    cdef np.ndarray outcodes, outdists
    obs = np.ascontiguousarray(obs)
    codes = np.ascontiguousarray(codes)

    if obs.ndim == 2:
        nobs = obs.shape[0]
        nfeat = obs.shape[1]
    elif obs.ndim == 1:
        nfeat = obs.shape[0]
        nobs = 1
    else:
        raise ValueError('obs must have 0 < obs.ndim <= 2')
    
    if codes.ndim == 2:
        ncodes = codes.shape[0]
        nfeat_codes = codes.shape[1]
    elif obs.ndim == 1:
        nfeat_codes = codes.shape[0]
        ncodes = 1
    else:
        raise ValueError('codes must have 0 < codes.ndim <= 2')

    if nfeat_codes != nfeat:
        raise ValueError('obs and codes must have same # of ' + \
                         'features (columns)')
    
    outcodes = np.empty((nobs,), dtype=np.int32)
    outdists = np.empty((nobs,), dtype=obs.dtype)
    
    if obs.dtype == np.float32:
        float_tvq(<float32_t *>obs.data, <float32_t *>codes.data, 
                  ncodes, nfeat, nobs, <int32_t *>outcodes.data, 
                  <float32_t *>outdists.data)
     
    return outcodes, outdists

