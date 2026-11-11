/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/SimplexConst.h
 * @brief Constants for HiGHS simplex solvers
 */
#ifndef SIMPLEX_SIMPLEXCONST_H_
#define SIMPLEX_SIMPLEXCONST_H_

#include "util/HighsType.h"

enum class SimplexAlgorithm { kNone = 0, kPrimal, kDual };

enum SimplexStrategy {
  kSimplexStrategyMin = 0,
  kSimplexStrategyChoose = kSimplexStrategyMin,      // 0
  kSimplexStrategyDual,                              // 1
  kSimplexStrategyDualPlain = kSimplexStrategyDual,  // 1
  kSimplexStrategyDualTasks,                         // 2
  kSimplexStrategyDualMulti,                         // 3
  kSimplexStrategyPrimal,                            // 4
  kSimplexStrategyMax = kSimplexStrategyPrimal,
  kSimplexStrategyNum
};

enum SimplexUnscaledSolutionStrategy {
  kSimplexUnscaledSolutionStrategyMin = 0,
  kSimplexUnscaledSolutionStrategyNone =
      kSimplexUnscaledSolutionStrategyMin,  // 0
  kSimplexUnscaledSolutionStrategyRefine,   // 1
  kSimplexUnscaledSolutionStrategyDirect,   // 2
  kSimplexUnscaledSolutionStrategyMax = kSimplexUnscaledSolutionStrategyDirect,
  kSimplexUnscaledSolutionStrategyNum
};

enum SimplexSolvePhase {
  kSolvePhaseMin = -3,
  kSolvePhaseError = kSolvePhaseMin,   // -3
  kSolvePhaseExit,                     // -2,
  kSolvePhaseUnknown,                  // -1
  kSolvePhaseOptimal,                  // 0
  kSolvePhase1,                        // 1
  kSolvePhase2,                        // 2
  kSolvePhasePrimalInfeasibleCleanup,  // 3
  kSolvePhaseOptimalCleanup,           // 4
  kSolvePhaseTabooBasis,               // 5
  kSolvePhaseMax = kSolvePhaseTabooBasis
};

enum SimplexCrashStrategy {
  kSimplexCrashStrategyMin = 0,
  kSimplexCrashStrategyOff = kSimplexCrashStrategyMin,
  kSimplexCrashStrategyLtssfK,
  kSimplexCrashStrategyLtssf = kSimplexCrashStrategyLtssfK,
  kSimplexCrashStrategyBixby,
  kSimplexCrashStrategyLtssfPri,
  kSimplexCrashStrategyLtsfK,
  kSimplexCrashStrategyLtsfPri,
  kSimplexCrashStrategyLtsf,
  kSimplexCrashStrategyBixbyNoNonzeroColCosts,
  kSimplexCrashStrategyBasic,
  kSimplexCrashStrategyTestSing,
  kSimplexCrashStrategyMax = kSimplexCrashStrategyTestSing
};

enum SimplexEdgeWeightStrategy {
  kSimplexEdgeWeightStrategyMin = -1,
  kSimplexEdgeWeightStrategyChoose = kSimplexEdgeWeightStrategyMin,
  kSimplexEdgeWeightStrategyDantzig,
  kSimplexEdgeWeightStrategyDevex,
  kSimplexEdgeWeightStrategySteepestEdge,
  kSimplexEdgeWeightStrategyMax = kSimplexEdgeWeightStrategySteepestEdge
};

enum SimplexPriceStrategy {
  kSimplexPriceStrategyMin = 0,
  kSimplexPriceStrategyCol = kSimplexPriceStrategyMin,
  kSimplexPriceStrategyRow,
  kSimplexPriceStrategyRowSwitch,
  kSimplexPriceStrategyRowSwitchColSwitch,
  kSimplexPriceStrategyMax = kSimplexPriceStrategyRowSwitchColSwitch
};

enum SimplexPivotalRowRefinementStrategy {
  kSimplexInfeasibilityProofRefinementMin = 0,
  kSimplexInfeasibilityProofRefinementNo =
      kSimplexInfeasibilityProofRefinementMin,
  kSimplexInfeasibilityProofRefinementUnscaledLp,    // 1
  kSimplexInfeasibilityProofRefinementAlsoScaledLp,  // 2
  kSimplexInfeasibilityProofRefinementMax =
      kSimplexInfeasibilityProofRefinementAlsoScaledLp
};

enum SimplexPrimalCorrectionStrategy {
  kSimplexPrimalCorrectionStrategyNone = 0,
  kSimplexPrimalCorrectionStrategyInRebuild,
  kSimplexPrimalCorrectionStrategyAlways,
  //  kSimplexPrimalCorrectionStrategyRefined
};

// Not an enum class since rebuild_reason is used in so many places
enum RebuildReason {
  kRebuildReasonCleanup = -1,
  kRebuildReasonNo = 0,
  kRebuildReasonUpdateLimitReached,               // 1
  kRebuildReasonSyntheticClockSaysInvert,         // 2
  kRebuildReasonPossiblyOptimal,                  // 3
  kRebuildReasonPossiblyPhase1Feasible,           // 4
  kRebuildReasonPossiblyPrimalUnbounded,          // 5
  kRebuildReasonPossiblyDualUnbounded,            // 6
  kRebuildReasonPossiblySingularBasis,            // 7
  kRebuildReasonPrimalInfeasibleInPrimalSimplex,  // 8
  kRebuildReasonChooseColumnFail,                 // 9
  kRebuildReasonForceRefactor,                    // 10
  kRebuildReasonExcessivePrimalValue,             // 11
  kRebuildReasonCount
};

enum SimplexNlaOperation {
  kSimplexNlaNull = -1,
  kSimplexNlaBtranFull = 0,
  kSimplexNlaPriceFull,
  kSimplexNlaBtranBasicFeasibilityChange,
  kSimplexNlaPriceBasicFeasibilityChange,
  kSimplexNlaBtranEp,
  kSimplexNlaPriceAp,
  kSimplexNlaFtran,
  kSimplexNlaFtranBfrt,
  kSimplexNlaFtranDse,
  kSimplexNlaBtranPse,
  kNumSimplexNlaOperation
};

enum class EdgeWeightMode { kDantzig = 0, kDevex, kSteepestEdge, kCount };

const HighsInt kDualTasksMinConcurrency = 3;
const HighsInt kDualMultiMinConcurrency = 1;  // 2;

// Simplex nonbasicFlag status for columns and rows. Don't use enum
// class since they are used as HighsInt to replace conditional
// statements by multiplication
const int8_t kNonbasicFlagTrue = 1;   // Nonbasic
const int8_t kNonbasicFlagFalse = 0;  // Basic
const int8_t kIllegalFlagValue =
    -99;  // Used to see whether valid flag value has been set

// Simplex nonbasicMove status for columns and rows. Don't use enum
// class since they are used in conditional statements
const int8_t kNonbasicMoveUp = 1;   // Free to move (only) up
const int8_t kNonbasicMoveDn = -1;  // Free to move (only) down
const int8_t kNonbasicMoveZe = 0;   // Fixed or free to move up and down
const int8_t kIllegalMoveValue =
    -99;  // Used to see whether valid move value has been set

// Threshold for accepting updated DSE weight
const double kAcceptDseWeightThreshold = 0.25;

const double kMinDualSteepestEdgeWeight = 1e-4;

const HighsInt kNoRowSought = -2;
const HighsInt kNoRowChosen = -1;

const HighsInt kNoRayIndex = -1;
const HighsInt kNoRaySign = 0;

// Switch to use code to check that, unless the basis supplied by the
// MIP solver was alien, the simplex solver starts from dual
// feasibility.
const bool kDebugMipNodeDualFeasible = false;

enum class LpAction {
  kScale = 0,
  kNewCosts,
  kNewBounds,
  kNewBasis,
  kNewCols,
  kNewRows,
  kDelCols,
  kDelNonbasicCols,
  kDelRows,
  kDelRowsBasisOk,
  kScaledCol,
  kScaledRow,
  kBacktracking
};

enum class BadBasisChangeReason {
  kAll = 0,
  kSingular,
  kCycling,
  kFailedInfeasibilityProof
};

const HighsInt kAllowedNumBadDevexWeight = 3;
const double kBadDevexWeightFactor = 3;

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
// Simplex: nonbasicFlag value of kNonbasicFlagFalse
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
// Simplex: nonbasicMove value of kNonbasicMoveUp - [l, Inf] column free to
// move up and negative dual
//
// Highs: col_status value of ZERO - at zero
// =>
// Simplex: nonbasicMove value of kNonbasicMoveZe - free variable treated
// specially in simplex
//
// Highs: col_status value of UPPER - at upper bound
// =>
// Simplex: Either
// * nonbasicMove value of kNonbasicMoveDn - [-Inf, u] column free to move down
// and positive dual
// * nonbasicMove value of kNonbasicMoveZe - [   l, u] column ?? and free dual
//
// Simplex: nonbasicMove value of kNonbasicMoveDn - [-Inf, u] column free to
// move down and positive dual
// =>
// Highs: col_status value of UPPER - at upper bound
//
// Simplex: nonbasicMove value of kNonbasicMoveZe - [l, u] column ?? and free
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
