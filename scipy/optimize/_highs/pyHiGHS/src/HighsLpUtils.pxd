# distutils: language=c++
# cython: language_level=3

from libcpp cimport bool

from .HighsStatus cimport HighsStatus
from .HighsLp cimport HighsLp
from .HighsOptions cimport HighsOptions

cdef extern from "HighsLpUtils.h" nogil:
    HighsStatus assessLp(HighsLp& lp, const HighsOptions& options, const bool normalise)
