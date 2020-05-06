# distutils: language=c++
# cython: language_level=3

from HighsOptions cimport HighsOptions
from HighsLp cimport HighsLp
from HighsStatus cimport HighsStatus

cdef extern from "LoadProblem.h" nogil:
    HighsStatus loadLpFromFile(const HighsOptions& options, HighsLp& lp)
