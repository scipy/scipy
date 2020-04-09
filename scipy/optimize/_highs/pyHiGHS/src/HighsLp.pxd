# distutils: language=c++
# cython: language_level=3

from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from HConst cimport HighsBasisStatus

cdef extern from "HighsLp.h" nogil:
    cdef cppclass HighsLp:
        int numCol_
        int numRow_
        int numInt_
        int nnz_

        string model_name_

        vector[int] Astart_
        vector[int] Aindex_
        vector[double] Avalue_
        vector[double] colCost_
        vector[double] colLower_
        vector[double] colUpper_
        vector[double] rowLower_
        vector[double] rowUpper_

    cdef enum HighsModelStatus:
        HighsModelStatusNOTSET "HighsModelStatus::NOTSET"
        HighsModelStatusLOAD_ERROR "HighsModelStatus::LOAD_ERROR"
        HighsModelStatusMODEL_ERROR "HighsModelStatus::MODEL_ERROR"
        HighsModelStatusMODEL_EMPTY "HighsModelStatus::MODEL_EMPTY"
        HighsModelStatusPRESOLVE_ERROR "HighsModelStatus::PRESOLVE_ERROR"
        HighsModelStatusSOLVE_ERROR "HighsModelStatus::SOLVE_ERROR"
        HighsModelStatusPOSTSOLVE_ERROR "HighsModelStatus::POSTSOLVE_ERROR"
        HighsModelStatusPRIMAL_INFEASIBLE "HighsModelStatus::PRIMAL_INFEASIBLE"
        HighsModelStatusPRIMAL_UNBOUNDED "HighsModelStatus::PRIMAL_UNBOUNDED"
        HighsModelStatusOPTIMAL "HighsModelStatus::OPTIMAL"
        HighsModelStatusREACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND "HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND"
        HighsModelStatusREACHED_TIME_LIMIT "HighsModelStatus::REACHED_TIME_LIMIT"
        HighsModelStatusREACHED_ITERATION_LIMIT "HighsModelStatus::REACHED_ITERATION_LIMIT"

    cdef cppclass HighsSolution:
        vector[double] col_value
        vector[double] col_dual
        vector[double] row_value
        vector[double] row_dual

    cdef cppclass HighsBasis:
        bool valid_
        vector[HighsBasisStatus] col_status
        vector[HighsBasisStatus] row_status
