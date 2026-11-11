/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HFactorConst.h
 * @brief Constants for basis matrix factorization, update and solves for HiGHS
 */
#ifndef HFACTORCONST_H_
#define HFACTORCONST_H_

#include "util/HighsType.h"

enum UPDATE_METHOD {
  kUpdateMethodFt = 1,
  kUpdateMethodPf = 2,
  kUpdateMethodMpf = 3,
  kUpdateMethodApf = 4
};
/**
 * Limits and default value of pivoting threshold
 */
const double kMinPivotThreshold = 8e-4;
const double kDefaultPivotThreshold = 0.1;
const double kPivotThresholdChangeFactor = 5.0;
const double kMaxPivotThreshold = 0.5;
/**
 * Limits and default value of minimum absolute pivot
 */
const double kMinPivotTolerance = 0;
const double kDefaultPivotTolerance = 1e-10;
const double kMaxPivotTolerance = 1.0;
/**
 * Necessary thresholds for expected density to trigger
 * hyper-sparse TRANs,
 */
const double kHyperFtranL = 0.15;
const double kHyperFtranU = 0.10;
const double kHyperBtranL = 0.10;
const double kHyperBtranU = 0.15;
/**
 * Necessary threshold for RHS density to trigger hyper-sparse TRANs,
 */
const double kHyperCancel = 0.05;
/**
 * Threshold for result density for it to be considered as
 * hyper-sparse - only for reporting
 */
const double kHyperResult = 0.10;

/**
 * Parameters for reinversion on synthetic clock
 */
const double kMultiBuildSyntheticTickMu = 1.0;
const double kNumericalTroubleTolerance = 1e-7;
const double kMultiNumericalTroubleTolerance = 1e-7;

const HighsInt kSyntheticTickReinversionMinUpdateCount = 50;
const HighsInt kMultiSyntheticTickReinversionMinUpdateCount =
    kSyntheticTickReinversionMinUpdateCount;

// Constants defining the space available for dimension-related
// identifiers like starts, and multipliers (of
// basis_matrix_limit_size, the basis matrix limit size) for
// fill-related identifiers like indices/values in Markowitz, and
// update.
const HighsInt kMCExtraEntriesMultiplier = 2;
const HighsInt kMRExtraEntriesMultiplier = 2;
const HighsInt kLFactorExtraEntriesMultiplier = 3;
const HighsInt kUFactorExtraVectors = 1000;
const HighsInt kUFactorExtraEntriesMultiplier = 3;
const HighsInt kPFFPivotEntries = 1000;
const HighsInt kPFVectors = 2000;
const HighsInt kPFEntriesMultiplier = 4;
const HighsInt kNewLRRowsExtraNz = 100;

enum ReportLuOption { kReportLuJustL = 1, kReportLuJustU, kReportLuBoth };

#endif /* HFACTORCONST_H_ */
