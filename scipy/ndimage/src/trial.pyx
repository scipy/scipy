#Input arrays nrrd to be matched from the original function
#Byteswapping needs to be used.

######################################################################
# Cython version of ndimage.morphology file
######################################################################
cimport cython
from cython cimport sizeof
import numpy as np
cimport numpy as np
from libc.math cimport sqrt

np.import_array()

DEF DISTANCE_EUCLIDIAN = 1
DEF DISTANCE_CITY_BLOCK = 2
DEF NI_DISTANCE_CHESSBOARD = 3

cdef extern from 'float.h':
    double DBL_MAX

cdef extern from "limits.h":
    long long LONG_MAX
    unsigned int UINT_MAX

cdef extern from *:
   ctypedef int Py_intptr_t

cdef extern from "numpy/arrayobject.h" nogil:
    ctypedef struct PyArrayIterObject:
        np.intp_t *coordinates
        char *dataptr
        np.npy_bool contiguous

    void *PyDataMem_NEW(size_t)
    void PyDataMem_FREE(void *)

def fun(np.ndarray a = None):
    if a is not None:
        print 'yes'
