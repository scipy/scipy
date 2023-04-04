# distutils: language=c++
# cython: language_level=3

from libcpp.vector cimport vector

from .HConst cimport MatrixFormat


cdef extern from "HighsSparseMatrix.h" nogil:
    cdef cppclass HighsSparseMatrix:
        MatrixFormat format_
        int num_col_
        int num_row_
        vector[int] start_
        vector[int] index_
        vector[double] value_
