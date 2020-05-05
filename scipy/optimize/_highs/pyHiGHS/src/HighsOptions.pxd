# distutils: language=c++
# cython: language_level=3

from libc.stdio cimport FILE

from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "HighsOptions.h" nogil:
    cdef cppclass HighsOptions:
        FILE * output
        int message_level
        bool mip
        string solution_file
        bool write_solution_to_file
        bool write_solution_pretty
