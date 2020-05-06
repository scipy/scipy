# distutils: language=c++
# cython: language_level=3

cdef extern from "HConst.h" nogil:

    const int HIGHS_CONST_I_INF = 2147483647
    const double HIGHS_CONST_INF = 1e200
    const double HIGHS_CONST_TINY = 1e-14
    const double HIGHS_CONST_ZERO = 1e-50
    const int HIGHS_THREAD_LIMIT = 8

    cdef enum HighsPrintMessageLevel:
        ML_MIN = 0
        ML_NONE = ML_MIN
        ML_VERBOSE = 1
        ML_DETAILED = 2
        ML_MINIMAL = 4
        ML_ALWAYS = ML_VERBOSE | ML_DETAILED | ML_MINIMAL
        ML_MAX = ML_ALWAYS

    cdef enum HighsBasisStatus:
        LOWER "HighsBasisStatus::LOWER" = 0, # (slack) variable is at its lower bound [including fixed variables]
        BASIC "HighsBasisStatus::BASIC" # (slack) variable is basic
        UPPER "HighsBasisStatus::UPPER" # (slack) variable is at its upper bound
        ZERO "HighsBasisStatus::ZERO" # free variable is non-basic and set to zero
        NONBASIC "HighsBasisStatus::NONBASIC" # nonbasic with no specific bound information - useful for users and postsolve
        SUPER "HighsBasisStatus::SUPER"

    cdef enum SolverOption:
        SOLVER_OPTION_SIMPLEX "SolverOption::SOLVER_OPTION_SIMPLEX" = -1
        SOLVER_OPTION_CHOOSE "SolverOption::SOLVER_OPTION_CHOOSE"
        SOLVER_OPTION_IPM "SolverOption::SOLVER_OPTION_IPM"

    cdef enum PrimalDualStatus:
        PrimalDualStatusSTATUS_NOT_SET "PrimalDualStatus::STATUS_NOT_SET" = -1
        PrimalDualStatusSTATUS_NO_SOLUTION "PrimalDualStatus::STATUS_NO_SOLUTION"
        PrimalDualStatusSTATUS_UNKNOWN "PrimalDualStatus::STATUS_UNKOWN"
        PrimalDualStatusSTATUS_INFEASIBLE_POINT "PrimalDualStatus::STATUS_INFEASIBLE_POINT"
        PrimalDualStatusSTATUS_FEASIBLE_POINT "PrimalDualStatus::STATUS_FEASIBLE_POINT"
