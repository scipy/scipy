# distutils: language=c++
# cython: language_level=3

from libc.stdio cimport FILE

cdef extern from "HighsIO.h" nogil:
    cdef enum HighsLogType:
        kInfo "HighsLogType::kInfo" = 1
        kDetailed "HighsLogType::kDetailed"
        kVerbose "HighsLogType::kVerbose"
        kWarning "HighsLogType::kWarning"
        kError "HighsLogType::kError"
