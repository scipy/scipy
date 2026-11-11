/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HConst.h
 * @brief Constants for HiGHS
 */
#ifndef LP_DATA_HCONST_H_
#define LP_DATA_HCONST_H_

#include <cmath>
#include <limits>
#include <string>

#include "util/HighsType.h"

const std::string kHighsCopyrightStatement =
    "Copyright (c) 2026 under MIT licence terms";

const size_t kHighsSize_tInf = (std::numeric_limits<size_t>::max)();
const HighsInt kHighsIInf = (std::numeric_limits<HighsInt>::max)();
const HighsInt kHighsIInf32 = (std::numeric_limits<int>::max)();
const double kHighsInf = std::numeric_limits<double>::infinity();
const double kHighsUndefined = kHighsInf;
const double kHighsTiny = 1e-14;
const double kHighsMacheps = std::ldexp(1, -52);
// 2^52: doubles of larger magnitude cannot represent fractional parts
const double kHighsInvMacheps = 1 / kHighsMacheps;
const double kHighsZero = 1e-50;
const std::string kHighsOffString = "off";
const std::string kHighsChooseString = "choose";
const std::string kHighsOnString = "on";
const std::string kHighsSparseString = "sparse";
const std::string kHighsPartialString = "partial";
const HighsInt kHighsMaxStringLength = 512;
const HighsInt kSimplexConcurrencyLimit = 8;
const double kRunningAverageMultiplier = 0.05;

const double kExcessivelySmallObjectiveCoefficient = 1e-4;
const double kExcessivelyLargeObjectiveCoefficient = 1e6;
const double kExcessivelySmallBoundValue = 1e-4;
const double kExcessivelyLargeBoundValue = 1e6;
const double kExcessivelyLargeIntegerBoundValue = 1e5;

const HighsInt kNoThreadInstance = -1;
const bool kAllowDeveloperAssert = false;
const bool kExtendInvertWhenAddingRows = false;

enum class HighsLogType { kInfo = 1, kDetailed, kVerbose, kWarning, kError };

enum SimplexScaleStrategy {
  kSimplexScaleStrategyMin = 0,
  kSimplexScaleStrategyOff = kSimplexScaleStrategyMin,  // 0
  kSimplexScaleStrategyChoose,                          // 1
  kSimplexScaleStrategyEquilibration,                   // 2
  kSimplexScaleStrategyForcedEquilibration,             // 3
  kSimplexScaleStrategyMaxValue,                        // 4
  kSimplexScaleStrategyMaxValue015 = kSimplexScaleStrategyMaxValue,
  kSimplexScaleStrategyMaxValue0157 = kSimplexScaleStrategyMaxValue,
  kSimplexScaleStrategyMax = kSimplexScaleStrategyMaxValue
};

enum HighsDebugLevel {
  kHighsDebugLevelNone = 0,
  kHighsDebugLevelCheap,
  kHighsDebugLevelCostly,
  kHighsDebugLevelExpensive,
  kHighsDebugLevelMin = kHighsDebugLevelNone,
  kHighsDebugLevelMax = kHighsDebugLevelExpensive
};

enum class HighsDebugStatus {
  kNotChecked = -1,
  kOk,
  kSmallError,
  kWarning,
  kLargeError,
  kError,
  kExcessiveError,
  kLogicalError,
};

enum HighsAnalysisLevel {
  kHighsAnalysisLevelNone = 0,
  kHighsAnalysisLevelModelData = 1,
  kHighsAnalysisLevelSolverSummaryData = 2,
  kHighsAnalysisLevelSolverRuntimeData = 4,
  kHighsAnalysisLevelSolverTime = 8,
  kHighsAnalysisLevelNlaData = 16,
  kHighsAnalysisLevelNlaTime = 32,
  kHighsAnalysisLevelMipData = 64,
  kHighsAnalysisLevelMipTime = 128,
  kHighsAnalysisLevelPresolveTime = 256,
  kHighsAnalysisLevelMin = kHighsAnalysisLevelNone,
  kHighsAnalysisLevelMax =
      kHighsAnalysisLevelModelData + kHighsAnalysisLevelSolverSummaryData +
      kHighsAnalysisLevelSolverRuntimeData + kHighsAnalysisLevelSolverTime +
      kHighsAnalysisLevelNlaData + kHighsAnalysisLevelNlaTime +
      kHighsAnalysisLevelMipData + kHighsAnalysisLevelMipTime +
      kHighsAnalysisLevelPresolveTime
};

enum class HighsVarType : uint8_t {
  kContinuous = 0,
  kInteger = 1,
  kSemiContinuous = 2,
  kSemiInteger = 3,
  kImplicitInteger = 4,
};

enum class HighsOptionType { kBool = 0, kInt, kDouble, kString };

enum class HighsInfoType { kInt64 = -1, kInt = 1, kDouble };

enum class HighsRunDataType { kInt64 = -1, kInt = 1, kDouble };

enum OptionOffChooseOn {
  kHighsOptionOff = -1,
  kHighsOptionChoose,
  kHighsOptionOn
};

enum IpxDualizeStrategy {
  kIpxDualizeStrategyOff = kHighsOptionOff,
  kIpxDualizeStrategyChoose = kHighsOptionChoose,
  kIpxDualizeStrategyOn = kHighsOptionOn,
  kIpxDualizeStrategyLukas,
  kIpxDualizeStrategyFilippo,
  kIpxDualizeStrategyMin = kIpxDualizeStrategyOff,
  kIpxDualizeStrategyMax = kIpxDualizeStrategyFilippo,
};

/** SCIP/HiGHS Objective sense */
enum class ObjSense { kMinimize = 1, kMaximize = -1 };

enum class MatrixFormat { kColwise = 1, kRowwise, kRowwisePartitioned };

enum class HessianFormat { kTriangular = 1, kSquare };

enum SolutionStatus {
  kSolutionStatusNone = 0,
  kSolutionStatusInfeasible,
  kSolutionStatusFeasible,
  kSolutionStatusMin = kSolutionStatusNone,
  kSolutionStatusMax = kSolutionStatusFeasible
};

enum BasisValidity {
  kBasisValidityInvalid = 0,
  kBasisValidityValid,
  kBasisValidityMin = kBasisValidityInvalid,
  kBasisValidityMax = kBasisValidityValid
};

enum CentringStatus {
  kCentringStatusNone = 0,
  kCentringStatusFailed,
  kCentringStatusSuccess,
  kCentringStatusMin = kCentringStatusNone,
  kCentringStatusMax = kCentringStatusSuccess
};

const std::string kHighsBasisFileV1 = "v1";  // Deprecated
const std::string kHighsBasisFileV2 = "v2";

enum SolutionStyle {
  kSolutionStyleOldRaw = -1,
  kSolutionStyleRaw = 0,
  kSolutionStylePretty,        // 1;
  kSolutionStyleGlpsolRaw,     // 2;
  kSolutionStyleGlpsolPretty,  // 3;
  kSolutionStyleSparse,        // 4;
  kSolutionStylePartial,       // 5;
  kSolutionStyleMin = kSolutionStyleOldRaw,
  kSolutionStyleMax = kSolutionStylePartial
};

enum GlpsolCostRowLocation {
  kGlpsolCostRowLocationLast = -2,
  kGlpsolCostRowLocationNone,         // -1
  kGlpsolCostRowLocationNoneIfEmpty,  // 0
  kGlpsolCostRowLocationMin = kGlpsolCostRowLocationLast
};

const std::string kHighsFilenameDefault = "";
const std::string kHighsMinimalColNamePrefix = "c";
const std::string kHighsMinimalrowNamePrefix = "r";
const std::string kHighsUniqueColNamePrefix = "c_ekk";
const std::string kHighsUniquerowNamePrefix = "r_ekk";

enum class HighsPresolveStatus {
  kNotPresolved = -1,
  kNotReduced,
  kInfeasible,
  kUnboundedOrInfeasible,
  kReduced,
  kReducedToEmpty,
  kTimeout,
  kNullError,     // V2.0: Delete since it's not used!
  kOptionsError,  // V2.0: Delete since it's not used!
  kNotSet,
  kOutOfMemory,  // V2.0: Move above kNotSet
  kException,    // V2.0: Move above kNotSet
  kMatrixValueError
};

enum class HighsPostsolveStatus {  // V2.0: Delete if not used!
  kNotPresolved = -1,
  kNoPrimalSolutionError,
  kSolutionRecovered,
  kBasisError
};

enum class HighsModelStatus {
  // NB Add new status values to the end so that int cast of status
  // values is unchanged, since enums are not preserved in some
  // interfaces
  kNotset = 0,
  kLoadError,  // V2.0: Delete since it's not used!
  kModelError,
  kPresolveError,
  kSolveError,
  kPostsolveError,  // V2.0: Delete if not used! Add to documentation if used
  kModelEmpty,
  kOptimal,
  kInfeasible,
  kUnboundedOrInfeasible,
  kUnbounded,
  kObjectiveBound,
  kObjectiveTarget,
  kTimeLimit,
  kIterationLimit,
  // V2.0: put kUnknown after kSolutionLimit and kInterrupt - and then
  // modify kMax and highs_c_api.h, highs_csharp_api.cs,
  // highspy/highs_bindings.cpp
  kUnknown,
  kSolutionLimit,
  kInterrupt,
  kMemoryLimit,
  kHighsInterrupt,
  kMin = kNotset,
  kMax = kHighsInterrupt
};

enum HighsCallbackType : int {
  kCallbackMin = 0,
  kCallbackLogging = kCallbackMin,    // 0
  kCallbackSimplexInterrupt,          // 1
  kCallbackIpmInterrupt,              // 2
  kCallbackMipSolution,               // 3
  kCallbackMipImprovingSolution,      // 4
  kCallbackMipLogging,                // 5
  kCallbackMipInterrupt,              // 6
  kCallbackMipGetCutPool,             // 7
  kCallbackMipDefineLazyConstraints,  // 8
  kCallbackMipUserSolution,           // 9
  kCallbackQpFirstFeasiblePoint,      // 10
  kCallbackQpInterrupt,               // 11
  kCallbackMax = kCallbackQpInterrupt,
  kNumCallbackType
};

/** SCIP/CPLEX-like HiGHS basis status for columns and rows. */
enum class HighsBasisStatus : uint8_t {
  kLower =
      0,   // (slack) variable is at its lower bound [including fixed variables]
  kBasic,  // (slack) variable is basic
  kUpper,  // (slack) variable is at its upper bound
  kZero,   // free variable is nonbasic and set to zero
  kNonbasic  // nonbasic with no specific bound information - useful for users
             // and postsolve
};

// Types of presolve modes
enum PresolveMode : uint8_t {
  kPresolveModeMin = 0,
  kPresolveModeDefault = kPresolveModeMin,
  kPresolveModeInitialSweep,
  kPresolveModeLight,
  kPresolveModeFull,
  kPresolveModeMax = kPresolveModeFull
};

// Types of presolve rules
enum PresolveRuleType : int {
  kPresolveRuleIllegal = -1,
  kPresolveRuleMin = 0,
  kPresolveRuleEmptyRow = kPresolveRuleMin,
  kPresolveRuleSingletonRow,
  kPresolveRuleRedundantRow,
  kPresolveRuleEmptyCol,
  kPresolveRuleFixedCol,
  kPresolveRuleDominatedCol,
  kPresolveRuleForcingRow,
  kPresolveRuleForcingCol,
  kPresolveRuleFreeColSubstitution,
  kPresolveRuleDoubletonEquation,
  kPresolveRuleDependentEquations,
  kPresolveRuleDependentFreeCols,
  kPresolveRuleAggregator,
  kPresolveRuleParallelRowsAndCols,
  kPresolveRuleSparsify,
  kPresolveRuleProbing,
  kPresolveRuleEnumeration,
  kPresolveRuleDualFixing,
  kPresolveRuleZeroCostSingleton,
  kPresolveRuleColStuffing,
  kPresolveRuleDualFixProbing,
  kPresolveRuleInitialSweep,
  kPresolveRuleFourierMotzkin,
  kPresolveRuleWeaklyDominatedCol,
  kPresolveRuleMax = kPresolveRuleWeaklyDominatedCol,
  kPresolveRuleLastAllowOff = kPresolveRuleMax,
  kPresolveRuleCount
};

enum IisStrategy : int {
  kIisStrategyMin = 0,
  kIisStrategyLight = kIisStrategyMin,  // 0
  kIisStrategyFromRay = 1,
  kIisStrategyFromLp = 2,
  kIisStrategyIrreducible = 4,
  kIisStrategyColPriority = 8,
  kIisStrategyRelaxation = 16,
  kIisStrategyDefault = kIisStrategyLight,
  kIisStrategyMax = kIisStrategyFromRay + kIisStrategyFromLp +
                    kIisStrategyIrreducible + kIisStrategyColPriority +
                    kIisStrategyRelaxation
};

enum IisStatus : int {
  kIisStatusMin = -1,
  kIisStatusNotInConflict = kIisStatusMin,  // -1
  kIisStatusMaybeInConflict,                // 0
  kIisStatusInConflict,                     // 1
  kIisStatusMax = kIisStatusInConflict
};

enum HessianOracleCallType : int {
  kHessianOracleCallTypeMin = 0,
  kHessianOracleCallTypeEntry = kHessianOracleCallTypeMin,
  kHessianOracleCallTypeColumn,
  kHessianOracleCallTypeProduct,
  kHessianOracleCallTypeMax = kHessianOracleCallTypeProduct
};

enum MipChooseSubMipRecord : int {
  kMipRecord = -1,
  kChooseRecord,
  kSubMipRecord
};

enum PresolveSolvePostsolveIndex : int {
  kPresolveTime = 0,
  kSolveTime,
  kPostsolveTime,
  kToPresolveSolvePostsolve
};

enum SubSolverIndex : int {
  kFromSubSolver = kToPresolveSolvePostsolve,
  kSubSolverMip = kFromSubSolver,
  kSubSolverDuSimplexBasis,
  kSubSolverDuSimplexNoBasis,
  kSubSolverPrSimplexBasis,
  kSubSolverPrSimplexNoBasis,
  kSubSolverHipo,
  kSubSolverIpx,
  kSubSolverHipoAc,
  kSubSolverIpxAc,
  kSubSolverPdlp,
  kSubSolverQpAsm,
  kSubSolverSubMip,
  kLastSubSolver = kSubSolverSubMip,
  kToSubSolver = kLastSubSolver + 1
};

// Minimum and default KKT tolerance
const double kMinimumKktTolerance = 1e-10;
const double kDefaultKktTolerance = 1e-7;

// Minimum and default MIP tolerance
const double kMinimumMipTolerance = 1e-10;
const double kDefaultMipTolerance = 1e-6;

// Minimum and default IPM optimality tolerance
const double kMinimumIpmTolerance = 1e-12;
const double kDefaultIpmTolerance = 1e-1 * kDefaultKktTolerance;

// Default QP Hessian regularization value
const double kHessianRegularizationValue = 1e-7;

// Default and max allowed power-of-two matrix scale factor
const HighsInt kDefaultAllowedMatrixPow2Scale = 20;
const HighsInt kMaxAllowedMatrixPow2Scale = 30;

// Illegal values of num/max/sum infeasibility - used to indicate that true
// values aren't known
const double kHighsIllegalInfeasibilityMeasure = kHighsInf;
const HighsInt kHighsIllegalInfeasibilityCount = -1;

// Illegal values of num/max/sum residual - used to indicate that true
// values aren't known
const double kHighsIllegalResidualMeasure = kHighsInf;
const HighsInt kHighsIllegalResidualCount = -1;

// Illegal values for HighsError - used to indicate that true
// values aren't known
const double kHighsIllegalErrorValue = kHighsInf;
const HighsInt kHighsIllegalErrorIndex = -1;

// Illegal values for complementarity violations used to indicate that true
// values aren't known
const double kHighsIllegalComplementarityViolation = kHighsInf;
const HighsInt kHighsIllegalComplementarityCount = -1;

const double kHighsIllegalDoubleMeasure = -kHighsInf;
const HighsInt kHighsIllegalIntMeasure = -1;

// Tolerance on asymmetry in square Hessians
const double kSquareHessianAsymmetryTolerance = 1e-10;

// Maximum upper bound on semi-variables
const double kMaxSemiVariableUpper = 1e5;

// Limit on primal values being realistic
const double kExcessivePrimalValue = 1e25;

// Tolerance values for highsDoubleToString
const double kModelValueToStringTolerance = 1e-15;
const double kRangingValueToStringTolerance = 1e-13;
const double kHighsSolutionValueToStringTolerance = 1e-13;
const double kGlpsolSolutionValueToStringTolerance = 1e-12;

// Termination link in linked lists
const HighsInt kNoLink = -1;

const int8_t kPivotIllegal = -1;
const int8_t kPivotLogical = 0;
const int8_t kPivotUnit = 1;
const int8_t kPivotRowSingleton = 2;
const int8_t kPivotColSingleton = 3;
const int8_t kPivotMarkowitz = 4;

// For converting general LPs to form for PDLP
//
// Requires non-conforming names for cuPDLP-C
enum ConstraintType { EQ = 0, LEQ, GEQ, BOUND, FREE };

// Mask for switching off PDLP features
enum PdlpFeaturesOff {
  kPdlpAllFeaturesOn = 0,
  kPdlpScalingOff = 1,
  kPdlpRestartOff = 2,
  kPdlpAdaptiveStepSizeOff = 4,
  kPdlpAllFeaturesOff =
      kPdlpScalingOff + kPdlpRestartOff + kPdlpAdaptiveStepSizeOff
};

enum PdlpScalingBit {
  kPdlpScalingMin = 0,
  kPdlpScalingRuiz = 1,
  kPdlpScalingL2 = 2,
  kPdlpScalingPC = 4,
  kPdlpScalingMax = kPdlpScalingRuiz + kPdlpScalingL2 + kPdlpScalingPC
};

enum PdlpStepSizeStrategy {
  kPdlpStepSizeStrategyFixed = 0,
  kPdlpStepSizeStrategyMin = kPdlpStepSizeStrategyFixed,
  kPdlpStepSizeStrategyAdaptive,
  kPdlpStepSizeStrategyMalitskyPock,
  kPdlpStepSizeStrategyPid,
  kPdlpStepSizeStrategyMax = kPdlpStepSizeStrategyPid
};

enum PdlpRestartStrategy {
  kPdlpRestartStrategyOff = 0,
  kPdlpRestartStrategyMin = kPdlpRestartStrategyOff,
  kPdlpRestartStrategyFixed,
  kPdlpRestartStrategyAdaptive,
  kPdlpRestartStrategyHalpern,
  kPdlpRestartStrategyMax = kPdlpRestartStrategyHalpern
};

namespace hipo {
enum class ParallelTechnique {
  kMin = 0,
  kAnalyse = kMin,
  kOrderNE,
  kOrderAS,
  kNEStruct,
  kNEValues,
  kTree,
  kNode,
  kForwardSolve,
  kDiagonalSolve,
  kBackwardSolve,
  kCount,
  kMaxSum = (1 << kCount) - 1
};
}

#endif /* LP_DATA_HCONST_H_ */
