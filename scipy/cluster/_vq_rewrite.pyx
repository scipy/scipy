"""
Cython rewrite of the vector quantization module, originally written
in C at src/vq.c and the wrapper at src/vq_module.c. This should be
easier to maintain than old SWIG output.

Original C version by Damian Eads. 
Translated to Cython by David Warde-Farley, October 2009.
"""

import numpy as np
cimport numpy as np

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

cdef void float_tvq(float32_t *obs, float32_t *code_book,
                   int ncodes, int nfeat, int nobs,
                   np.npy_intp *codes, float32_t *low_dist):
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

def vq(np.ndarray obs, np.ndarray codes):
    """
    vq(obs, codes)

    Vector quantization ndarray wrapper.
    """
    cdef np.npy_intp nobs, ncodes, nfeat, nfeat_codes
    cdef np.ndarray obs_a, codes_a
    cdef np.ndarray outcodes, outdists
    cdef int flags = np.NPY_CONTIGUOUS | np.NPY_NOTSWAPPED | np.NPY_ALIGNED

    # Ensure the arrays are contiguous - done in C to minimize overhead.
    obs_a = np.PyArray_FROM_OF(obs, flags)
    codes_a = np.PyArray_FROM_OF(codes, flags)
     
    if obs_a.ndim == 2:
        nobs = obs_a.shape[0]
        nfeat = obs_a.shape[1]
    elif obs_a.ndim == 1:
        nfeat = obs_a.shape[0]
        nobs = 1
    else:
        raise ValueError('obs must have 0 < obs.ndim <= 2')
    
    if codes_a.ndim == 2:
        ncodes = codes_a.shape[0]
        nfeat_codes = codes_a.shape[1]
    elif codes.ndim == 1:
        # Treat one dimensional arrays as row vectors.
        nfeat_codes = codes_a.shape[0]
        ncodes = 1
    else:
        raise ValueError('codes must have 0 < codes.ndim <= 2')
    
    # Ensure the two arrays have the same number of features (columns).
    if nfeat_codes != nfeat:
        raise ValueError('obs and codes must have same # of ' + \
                         'features (columns)')
    
    # We create this with the C API so that we can be sure that
    # the resulting array has elements big enough to store indices
    # on that platform. Hence, PyArray_INTP.
    outcodes = PyArray_EMPTY(1, &nobs, PyArray_INTP, 0)
    
    # This we just want to match the dtype of the input, so np.empty is fine.
    outdists = np.empty((nobs,), dtype=obs.dtype)
    
    if obs.dtype == np.float32:
        float_tvq(<float32_t *>obs_a.data, <float32_t *>codes_a.data, 
                  ncodes, nfeat, nobs, <np.npy_intp *>outcodes.data, 
                  <float32_t *>outdists.data)
    
    return outcodes, outdists

