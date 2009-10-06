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

cdef double rbig = 1e100

# Python/NumPy types
DTYPE_FLOAT64 = np.float64
DTYPE_FLOAT32 = np.float32

# C types
ctypedef np.double_t float64_t
ctypedef np.float_t float_t

cdef int float_vq_obs(float *obs,       # FIXME: should be ndarray
                      float *code_book, # FIXME: should be ndarray
                      int Ncodes,
                      int Nfeatures,
                      int *code,
                      float *lowest_dist):
    """
    Quantize a single observation 'obs' to its nearest codebook 
    entry (single-precision version).
    """
    cdef float dist, diff
    cdef int i, j, k
    for i in range(Ncodes):
        dist = 0
        for j in range(Nfeatures):
            diff = code_book[k] - obs[j]
            dist += diff * diff
            k += 1
        dist = sqrtf(dist)
        if dist < lowest_dist[0]:
            code[0] = i
            lowest_dist[0] = dist
    return 0

cdef int double_vq_obs(double *obs,       # FIXME: should be ndarray
                       double *code_book, # FIXME: should be ndarray
                       int Ncodes, 
                       int Nfeatures,
                       int *code,
                       double *lowest_dist):
    """
    Quantize a single observation 'obs' to its nearest codebook 
    entry (double-precision version).
    """
    cdef double dist, diff
    cdef int i, j, k
    k = 0
    lowest_dist[0] = rbig 
    for i in range(Ncodes):
        dist = 0
        for j in range(Nfeatures):
            diff = code_book[k] - obs[j]
            dist += diff * diff
            k += 1
        dist = sqrt(dist)
        if dist < lowest_dist[0]:
            code[0]= i
            lowest_dist[0] = dist
    return 0

cdef int float_tvq(float *obs,        # FIXME: should be ndarray
                   float *code_book,  # FIXME: should be ndarray
                   int Nobs, int Ncodes, int Nfeatures,
                   int *codes,        # FIXME: should be ndarray
                   float *lowest_dist # FIXME: should be ndarray
                  ):
    """
    Quantize Nobs observations to their nearest codebook 
    entry (single-precision version).
    """
    cdef int i
    for i in range(Nobs):
        float_vq_obs(&(obs[i * Nfeatures]), code_book,
                     Ncodes, Nfeatures, &(codes[i]),
                     &(lowest_dist[i]))

cdef int double_tvq(double *obs,        # FIXME: should be ndarray
                    double *code_book,  # FIXME: should be ndarray
                    int Nobs, int Ncodes, int Nfeatures,
                    int *codes,         # FIXME: should be ndarray
                    double *lowest_dist # FIXME: should be ndarray
                   ):
    """
    Quantize Nobs observations to their nearest codebook 
    entry (double-precision version).
    """
    cdef int i
    for i in range(Nobs):
        double_vq_obs(&(obs[i * Nfeatures]), 
                      code_book, Ncodes, Nfeatures,  
                      &(codes[i]), &(lowest_dist[i]))
    return 0

