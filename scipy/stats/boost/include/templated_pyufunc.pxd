# distutils: language = c++
# cython: language_level=3

from numpy cimport npy_intp

cdef extern from "Templated_PyUFunc.hpp" namespace "" nogil:
    void PyUFunc_T[T, N](char** args, npy_intp *dimensions, npy_intp *steps, void *data) except +
