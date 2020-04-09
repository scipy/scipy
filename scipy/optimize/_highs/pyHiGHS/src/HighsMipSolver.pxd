# distutils: language=c++
# cython: language_level=3

from HighsOptions cimport HighsOptions
from HighsLp cimport HighsLp

cdef extern from "HighsMipSolver.h" nogil:
    cdef enum HighsMipStatus:
        HighsMipStatuskOptimal "HighsMipStatus::kOptimal"
        HighsMipStatuskTimeout "HighsMipStatus::kTimeout"
        HighsMipStatuskError "HighsMipStatus::kError"
        HighsMipStatuskNodeOptimal "HighsMipStatus::kNodeOptimal"
        HighsMipStatuskNodeInfeasible "HighsMipStatus::kNodeInfeasible"
        HighsMipStatuskNodeUnbounded "HighsMipStatus::kNodeUnbounded"
        HighsMipStatuskNodeNotOptimal "HighsMipStatus::kNodeNotOptimal"
        HighsMipStatuskNodeError "HighsMipStatus::kNodeError"
        HighsMipStatuskRootNodeOptimal "HighsMipStatus::kRootNodeOptimal"
        HighsMipStatuskRootNodeNotOptimal "HighsMipStatus::kRootNodeNotOptimal"
        HighsMipStatuskRootNodeError "HighsMipStatus::kRootNodeError"
        HighsMipStatuskMaxNodeReached "HighsMipStatus::kMaxNodeReached"
        HighsMipStatuskUnderDevelopment "HighsMipStatus::kUnderDevelopment"
        HighsMipStatuskTreeExhausted "HighsMipStatus::kTreeExhausted"

    cdef cppclass HighsMipSolver:
        HighsMipSolver(const HighsOptions& options, const HighsLp& lp) except +
        HighsMipStatus runMipSolver()
