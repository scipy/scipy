/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/SimplexConst.h
 * @brief Constants for HiGHS simplex solvers
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_SIMPLEXCONST_H_
#define SIMPLEX_SIMPLEXCONST_H_

enum class SimplexSolutionStatus {
  UNSET = -1,
  OPTIMAL,
  PRIMAL_FEASIBLE,
  DUAL_FEASIBLE,
  INFEASIBLE,
  UNBOUNDED,
  SINGULAR,
  FAILED,
  REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND,
  OUT_OF_TIME
};

enum class SimplexAlgorithm { PRIMAL = 0, DUAL };

enum SimplexStrategy {
  SIMPLEX_STRATEGY_MIN = 0,
  SIMPLEX_STRATEGY_CHOOSE = SIMPLEX_STRATEGY_MIN,
  SIMPLEX_STRATEGY_DUAL,
  SIMPLEX_STRATEGY_DUAL_PLAIN = SIMPLEX_STRATEGY_DUAL,
  SIMPLEX_STRATEGY_DUAL_TASKS,
  SIMPLEX_STRATEGY_DUAL_MULTI,
  SIMPLEX_STRATEGY_PRIMAL,
  SIMPLEX_STRATEGY_MAX = SIMPLEX_STRATEGY_PRIMAL,
  SIMPLEX_STRATEGY_NUM
};

enum SimplexSolvePhase {
  SOLVE_PHASE_MIN = -3,
  SOLVE_PHASE_ERROR = SOLVE_PHASE_MIN,  // -3
  SOLVE_PHASE_EXIT,                     // -2,
  SOLVE_PHASE_UNKNOWN,                  // -1
  SOLVE_PHASE_OPTIMAL,                  // 0
  SOLVE_PHASE_1,                        // 1
  SOLVE_PHASE_2,                        // 2
  SOLVE_PHASE_CLEANUP = 4,
  SOLVE_PHASE_MAX = SOLVE_PHASE_CLEANUP
};

enum DualSimplexCleanupStrategy {
  DUAL_SIMPLEX_CLEANUP_STRATEGY_MIN = 0,
  DUAL_SIMPLEX_CLEANUP_STRATEGY_NONE = DUAL_SIMPLEX_CLEANUP_STRATEGY_MIN,
  DUAL_SIMPLEX_CLEANUP_STRATEGY_HPRIMAL,
  DUAL_SIMPLEX_CLEANUP_STRATEGY_HQPRIMAL,
  DUAL_SIMPLEX_CLEANUP_STRATEGY_MAX = DUAL_SIMPLEX_CLEANUP_STRATEGY_HQPRIMAL
};

enum SimplexScaleStrategy {
  SIMPLEX_SCALE_STRATEGY_MIN = 0,
  SIMPLEX_SCALE_STRATEGY_OFF = SIMPLEX_SCALE_STRATEGY_MIN,
  SIMPLEX_SCALE_STRATEGY_HIGHS,
  SIMPLEX_SCALE_STRATEGY_HIGHS_FORCED,
  SIMPLEX_SCALE_STRATEGY_015,
  SIMPLEX_SCALE_STRATEGY_0157,
  SIMPLEX_SCALE_STRATEGY_MAX = SIMPLEX_SCALE_STRATEGY_0157
};

enum SimplexCrashStrategy {
  SIMPLEX_CRASH_STRATEGY_MIN = 0,
  SIMPLEX_CRASH_STRATEGY_OFF = SIMPLEX_CRASH_STRATEGY_MIN,
  SIMPLEX_CRASH_STRATEGY_LTSSF_K,
  SIMPLEX_CRASH_STRATEGY_LTSSF = SIMPLEX_CRASH_STRATEGY_LTSSF_K,
  SIMPLEX_CRASH_STRATEGY_BIXBY,
  SIMPLEX_CRASH_STRATEGY_LTSSF_PRI,
  SIMPLEX_CRASH_STRATEGY_LTSF_K,
  SIMPLEX_CRASH_STRATEGY_LTSF_PRI,
  SIMPLEX_CRASH_STRATEGY_LTSF,
  SIMPLEX_CRASH_STRATEGY_BIXBY_NO_NONZERO_COL_COSTS,
  SIMPLEX_CRASH_STRATEGY_BASIC,
  SIMPLEX_CRASH_STRATEGY_TEST_SING,
  SIMPLEX_CRASH_STRATEGY_MAX = SIMPLEX_CRASH_STRATEGY_TEST_SING
};

enum SimplexDualEdgeWeightStrategy {
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_MIN = -1,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_CHOOSE =
      SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_MIN,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL,
  SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_MAX =
      SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL
};

enum SimplexPrimalEdgeWeightStrategy {
  SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_MIN = -1,
  SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_CHOOSE =
      SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_MIN,
  SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DANTZIG,
  SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX,
  SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_MAX =
      SIMPLEX_PRIMAL_EDGE_WEIGHT_STRATEGY_DEVEX
};

enum SimplexPriceStrategy {
  SIMPLEX_PRICE_STRATEGY_MIN = 0,
  SIMPLEX_PRICE_STRATEGY_COL = SIMPLEX_PRICE_STRATEGY_MIN,
  SIMPLEX_PRICE_STRATEGY_ROW,
  SIMPLEX_PRICE_STRATEGY_ROW_SWITCH,
  SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH,
  SIMPLEX_PRICE_STRATEGY_MAX = SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH
};

enum SimplexDualChuzcStrategy {
  SIMPLEX_DUAL_CHUZC_STRATEGY_MIN = 0,
  SIMPLEX_DUAL_CHUZC_STRATEGY_CHOOSE = SIMPLEX_DUAL_CHUZC_STRATEGY_MIN,
  SIMPLEX_DUAL_CHUZC_STRATEGY_QUAD,
  SIMPLEX_DUAL_CHUZC_STRATEGY_HEAP,
  SIMPLEX_DUAL_CHUZC_STRATEGY_BOTH,
  SIMPLEX_DUAL_CHUZC_STRATEGY_MAX = SIMPLEX_DUAL_CHUZC_STRATEGY_BOTH
};

// Not an enum class since invert_hint is used in so many places
enum InvertHint {
  INVERT_HINT_NO = 0,
  INVERT_HINT_UPDATE_LIMIT_REACHED,
  INVERT_HINT_SYNTHETIC_CLOCK_SAYS_INVERT,
  INVERT_HINT_POSSIBLY_OPTIMAL,
  INVERT_HINT_POSSIBLY_PRIMAL_UNBOUNDED,
  INVERT_HINT_POSSIBLY_DUAL_UNBOUNDED,
  INVERT_HINT_POSSIBLY_SINGULAR_BASIS,
  INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX,
  INVERT_HINT_CHOOSE_COLUMN_FAIL,
  INVERT_HINT_Count
};

enum class DualEdgeWeightMode { DANTZIG = 0, DEVEX, STEEPEST_EDGE, Count };

enum class PriceMode { ROW = 0, COL };

const int PARALLEL_THREADS_DEFAULT = 8;
const int DUAL_TASKS_MIN_THREADS = 3;
const int DUAL_MULTI_MIN_THREADS = 1;  // 2;

// TODO: Set this false tactically to make mip interface more
// efficient by preventing reinversion on optimality in phase 1 or
// phase 2
const bool invert_if_row_out_negative = true;

/** Simplex nonbasicFlag status for columns and rows. Don't use enum
    class since they are used as int to replace conditional statements
    by multiplication */
const int NONBASIC_FLAG_TRUE = 1;   // Nonbasic
const int NONBASIC_FLAG_FALSE = 0;  // Basic

/** Simplex nonbasicMove status for columns and rows. Don't use enum
    class since they are used in conditional statements */
const int NONBASIC_MOVE_UP = 1;   // Free to move (only) up
const int NONBASIC_MOVE_DN = -1;  // Free to move (only) down
const int NONBASIC_MOVE_ZE = 0;   // Fixed or free to move up and down
//
// Relation between HiGHS basis and Simplex basis
//
// Data structures
// ===============
//
// HiGHS basis consists of vectors
//
// * col_status[numCol]
// * row_status[numRow]
//
// Simplex basis consists of vectors
//
// * nonbasicMove[numTot]
// * basicIndex[numRow]
// * nonbasicFlag[numTot]
//
// where nonbasicFlag is duplicate information but is used to identify
// whether a particular variable is basic or nonbasic.
//
// Basic variables
// ===============
//
// Highs: *_status value of BASIC
//
// <=>
//
// Simplex: nonbasicFlag value of NONBASIC_FLAG_FALSE
//
// Nonbasic variables
// ==================
//
// Relations complicated by the fact that
//
// * HiGHS   rows have bounds [ l,  u]
// * Simplex rows have bounds [-u, -l]
//
// Nonbasic columns
// ================
//
// Highs: col_status value of LOWER - at lower bound
// <=>
// Simplex: nonbasicMove value of NONBASIC_MOVE_UP - [l, Inf] column free to
// move up and negative dual
//
// Highs: col_status value of ZERO - at zero
// =>
// Simplex: nonbasicMove value of NONBASIC_MOVE_ZE - free variable treated
// specially in simplex
//
// Highs: col_status value of UPPER - at upper bound
// =>
// Simplex: Either
// * nonbasicMove value of NONBASIC_MOVE_DN - [-Inf, u] column free to move down
// and positive dual
// * nonbasicMove value of NONBASIC_MOVE_ZE - [   l, u] column ?? and free dual
//
// Simplex: nonbasicMove value of NONBASIC_MOVE_DN - [-Inf, u] column free to
// move down and positive dual
// =>
// Highs: col_status value of UPPER - at upper bound
//
// Simplex: nonbasicMove value of NONBASIC_MOVE_ZE - [l, u] column ?? and free
// dual
// =>
// Highs: Either
// * col_status value of UPPER - at upper bound
// * col_status value of ZERO - at zero
//
// Nonbasic rows
// =============
//
#endif /* SIMPLEX_SIMPLEXCONST_H_ */
