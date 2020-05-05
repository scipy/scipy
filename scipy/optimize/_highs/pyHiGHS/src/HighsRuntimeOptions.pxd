# distutils: language=c++
# cython: language_level=3

from libcpp cimport bool

from HighsOptions cimport HighsOptions

cdef extern from "HighsRuntimeOptions.h" nogil:
    bool loadOptions(int argc, char** argv, HighsOptions& options)
