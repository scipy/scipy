''' Cython mio5 utility routines (-*- python -*- like)

'''

import sys

from copy import copy as pycopy

from libc.stdlib cimport calloc, free
from libc.string cimport strcmp, strlen

from cpython cimport Py_INCREF, Py_DECREF
from cpython cimport PyObject

cdef extern from "Python.h":
    ctypedef struct PyTypeObject:
        pass

from cpython cimport PyBytes_Size, PyBytes_FromString, \
    PyBytes_FromStringAndSize

import numpy as np
from numpy.compat import asbytes, asstr
cimport numpy as cnp

cdef extern from "numpy/arrayobject.h":
    PyTypeObject PyArray_Type
    cnp.ndarray PyArray_NewFromDescr(PyTypeObject *subtype,
                                     cnp.dtype newdtype,
                                     int nd,
                                     cnp.npy_intp* dims,
                                     cnp.npy_intp* strides,
                                     void* data,
                                     int flags,
                                     object parent)

cdef extern from "numpy_rephrasing.h":
    void PyArray_Set_BASE(cnp.ndarray arr, object obj)

# Numpy must be initialized before any code using the numpy C-API
# directly
cnp.import_array()

# Constant from numpy - max number of array dimensions
DEF _MAT_MAXDIMS = 32
# max number of integer indices of matlab data types (miINT8 etc)
DEF _N_MIS = 20
# max number of integer indices of matlab class types (mxINT8_CLASS etc)
DEF _N_MXS = 20

cimport streams
import scipy.io.matlab.miobase as miob
from scipy.io.matlab.mio_utils import squeeze_element, chars_to_strings
import scipy.io.matlab.mio5_params as mio5p
import scipy.sparse


cdef enum:
    miINT8 = 1
    miUINT8 = 2
    miINT16 = 3
    miUINT16 = 4
    miINT32 = 5
    miUINT32 = 6
    miSINGLE = 7
    miDOUBLE = 9
    miINT64 = 12
    miUINT64 = 13
    miMATRIX = 14
    miCOMPRESSED = 15
    miUTF8 = 16
    miUTF16 = 17
    miUTF32 = 18

cdef enum: # see comments in mio5_params
    mxCELL_CLASS = 1
    mxSTRUCT_CLASS = 2
    mxOBJECT_CLASS = 3
    mxCHAR_CLASS = 4
    mxSPARSE_CLASS = 5
    mxDOUBLE_CLASS = 6
    mxSINGLE_CLASS = 7
    mxINT8_CLASS = 8
    mxUINT8_CLASS = 9
    mxINT16_CLASS = 10
    mxUINT16_CLASS = 11
    mxINT32_CLASS = 12
    mxUINT32_CLASS = 13
    mxINT64_CLASS = 14
    mxUINT64_CLASS = 15
    mxFUNCTION_CLASS = 16
    mxOPAQUE_CLASS = 17 # This appears to be a function workspace
    mxOBJECT_CLASS_FROM_MATRIX_H = 18

sys_is_le = sys.byteorder == 'little'
native_code = sys_is_le and '<' or '>'
swapped_code = sys_is_le and '>' or '<'

cdef cnp.dtype OPAQUE_DTYPE = mio5p.OPAQUE_DTYPE


cpdef cnp.uint32_t byteswap_u4(cnp.uint32_t u4):
    return ((u4 << 24) |
           ((u4 << 8) & 0xff0000U) |
           ((u4 >> 8 & 0xff00u)) |
           (u4 >> 24))


cdef class VarHeader5:
    cdef readonly object name
    cdef readonly int mclass
    cdef readonly object dims
    cdef cnp.int32_t dims_ptr[_MAT_MAXDIMS]
    cdef int n_dims
    cdef int check_stream_limit
    cdef int is_complex
    cdef readonly int is_logical
    cdef public int is_global
    cdef readonly size_t nzmax

    def set_dims(self, dims):
        """ Allow setting of dimensions from python

        This is for constructing headers for tests
        """
        self.dims = dims
        self.n_dims = len(dims)
        for i, dim in enumerate(dims):
            self.dims_ptr[i] = <cnp.int32_t>int(dim)
