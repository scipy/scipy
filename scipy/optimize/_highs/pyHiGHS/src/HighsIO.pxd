# distutils: language=c++
# cython: language_level=3

from libc.stdio cimport FILE

cdef extern from "HighsIO.h" nogil:
    void HighsPrintMessage(FILE* pass_output, const int message_level, const int level, const char* format, ...)
