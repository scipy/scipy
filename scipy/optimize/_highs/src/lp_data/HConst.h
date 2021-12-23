/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HConst.h
 * @brief Constants for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HCONST_H_
#define LP_DATA_HCONST_H_

#include <limits>
#include <string>

const int HIGHS_CONST_I_INF = std::numeric_limits<int>::max();
const double HIGHS_CONST_INF = std::numeric_limits<double>::infinity();
const double HIGHS_CONST_TINY = 1e-14;
const double HIGHS_CONST_ZERO = 1e-50;
const std::string off_string = "off";
const std::string choose_string = "choose";
const std::string on_string = "on";
const int HIGHS_THREAD_LIMIT = 8;  // 32;

enum HighsDebugLevel {
  HIGHS_DEBUG_LEVEL_MIN = 0,
  HIGHS_DEBUG_LEVEL_NONE = HIGHS_DEBUG_LEVEL_MIN,  // 0
  HIGHS_DEBUG_LEVEL_CHEAP,                         // 1
  HIGHS_DEBUG_LEVEL_COSTLY,                        // 2
  HIGHS_DEBUG_LEVEL_EXPENSIVE,                     // 3
  HIGHS_DEBUG_LEVEL_MAX = HIGHS_DEBUG_LEVEL_EXPENSIVE
};

enum class HighsDebugStatus {
  NOT_CHECKED = -1,
  OK,
  SMALL_ERROR,
  WARNING,
  LARGE_ERROR,
  ERROR,
  EXCESSIVE_ERROR,
  LOGICAL_ERROR,
};

enum class HighsOptionType { BOOL = 0, INT, DOUBLE, STRING };

enum class HighsInfoType { INT = 1, DOUBLE };

enum OptionOffChooseOn { OPTION_OFF = -1, OPTION_CHOOSE, OPTION_ON };

/** SCIP/HiGHS Objective sense */
enum class ObjSense { MINIMIZE = 1, MAXIMIZE = -1 };

enum SolverOption {
  SOLVER_OPTION_SIMPLEX = -1,
  SOLVER_OPTION_CHOOSE,
  SOLVER_OPTION_IPM
};

enum PrimalDualStatus {
  STATUS_NOTSET = -1,
  STATUS_MIN = STATUS_NOTSET,
  STATUS_NO_SOLUTION,
  STATUS_UNKNOWN,
  STATUS_INFEASIBLE_POINT,
  STATUS_FEASIBLE_POINT,
  STATUS_MAX = STATUS_FEASIBLE_POINT
};

const std::string FILENAME_DEFAULT = "";

// Need to allow infinite costs to pass SCIP LPI unit tests
const bool allow_infinite_costs = true;

// Primal/dual statuses and corresponding HighsModelStatus
// values. Note that if dual infeasibility is identified, then the
// prototype primal code is used to distinguish PRIMAL_DUAL_INFEASIBLE
// from PRIMAL_UNBOUNDED. If this fails, then HiGHS may just return
// DUAL_INFEASIBLE
//
//           | Du Infeas    | Du Feas   | Du UnBd
// Pr Infeas | PR_DU_INFEAS | PR_INFEAS | PR_INFEAS
// Pr Feas   | PR_UNBD      | OPTIMAL   |   N/A
// Pr Unbd   | PR_UNBD      |     N/A   |   N/A
//
// Dual infeasibility is recognised by infeasibility at dual phase 1 optimality
// (and implied by primal unboundedness)
//
// Dual feasibility is recognised by feasibility at dual phase 1 optimality or
// primal phase 2 optimality
//
// Dual unboundedness is recognised by unboundedness in dual phase 2
//
// Primal infeasibility is recognised by infeasibility at primal phase 1
// optimality (and implied by dual unboundedness)
//
// Primal feasibility is recognised by feasibility at primal phase 1 optimality
// or dual phase 2 optimality
//
// Primal unboundedness is recognised by unboundedness in primal phase 2
//

enum class HighsModelStatus {
  // NB Add new status values to the end so that int cast of status
  // values is unchanged, since enums are not preserved in some
  // interfaces
  NOTSET = 0,
  HIGHS_MODEL_STATUS_MIN = NOTSET,
  LOAD_ERROR,
  MODEL_ERROR,
  PRESOLVE_ERROR,
  SOLVE_ERROR,
  POSTSOLVE_ERROR,
  MODEL_EMPTY,
  PRIMAL_INFEASIBLE,
  PRIMAL_UNBOUNDED,
  OPTIMAL,
  REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND,
  REACHED_TIME_LIMIT,
  REACHED_ITERATION_LIMIT,
  PRIMAL_DUAL_INFEASIBLE,
  DUAL_INFEASIBLE,
  HIGHS_MODEL_STATUS_MAX = DUAL_INFEASIBLE
};

/** SCIP/CPLEX-like HiGHS basis status for columns and rows. */
enum class HighsBasisStatus {
  LOWER =
      0,  // (slack) variable is at its lower bound [including fixed variables]
  BASIC,  // (slack) variable is basic
  UPPER,  // (slack) variable is at its upper bound
  ZERO,   // free variable is non-basic and set to zero
  NONBASIC,  // nonbasic with no specific bound information - useful for users
             // and postsolve
  SUPER      // Super-basic variable: non-basic and either free and
             // nonzero or not at a bound. No SCIP equivalent
};

#endif /* LP_DATA_HCONST_H_ */
