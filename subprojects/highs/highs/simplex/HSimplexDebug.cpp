/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplexDebug.cpp
 * @brief
 */

#include "simplex/HSimplexDebug.h"

#include <string>

#include "lp_data/HighsDebug.h"

void debugDualChuzcFailNorms(
    const HighsInt workCount,
    const std::vector<std::pair<HighsInt, double>>& workData,
    double& workDataNorm, const HighsInt numVar, const double* workDual,
    double& workDualNorm) {
  workDataNorm = 0;
  for (HighsInt i = 0; i < workCount; i++) {
    double value = workData[i].second;
    workDataNorm += value * value;
  }
  workDataNorm = sqrt(workDataNorm);
  workDualNorm = 0;
  for (HighsInt iVar = 0; iVar < numVar; iVar++) {
    double value = workDual[iVar];
    workDualNorm += value * value;
  }
  workDualNorm = sqrt(workDualNorm);
}

HighsDebugStatus debugDualChuzcFailQuad0(
    const HighsOptions& options, const HighsInt workCount,
    const std::vector<std::pair<HighsInt, double>>& workData,
    const HighsInt numVar, const double* workDual, const double selectTheta,
    const double remainTheta, const bool force) {
  // Non-trivially expensive assessment of CHUZC failure
  if (options.highs_debug_level < kHighsDebugLevelCostly && !force)
    return HighsDebugStatus::kNotChecked;

  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     No change in loop 2 so return error\n");
  double workDataNorm;
  double workDualNorm;
  debugDualChuzcFailNorms(workCount, workData, workDataNorm, numVar, workDual,
                          workDualNorm);
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     workCount = %" HIGHSINT_FORMAT
              "; selectTheta=%g; remainTheta=%g\n",
              workCount, selectTheta, remainTheta);
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     workDataNorm = %g; workDualNorm = %g\n",
              workDataNorm, workDualNorm);
  return HighsDebugStatus::kOk;
}

HighsDebugStatus debugDualChuzcFailQuad1(
    const HighsOptions& options, const HighsInt workCount,
    const std::vector<std::pair<HighsInt, double>>& workData,
    const HighsInt numVar, const double* workDual, const double selectTheta,
    const bool force) {
  // Non-trivially expensive assessment of CHUZC failure
  if (options.highs_debug_level < kHighsDebugLevelCostly && !force)
    return HighsDebugStatus::kNotChecked;

  highsLogDev(
      options.log_options, HighsLogType::kInfo,
      "DualChuzC:     No group identified in quad search so return error\n");
  double workDataNorm;
  double workDualNorm;
  debugDualChuzcFailNorms(workCount, workData, workDataNorm, numVar, workDual,
                          workDualNorm);
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     workCount = %" HIGHSINT_FORMAT
              "; selectTheta=%g\n",
              workCount, selectTheta);
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     workDataNorm = %g; workDualNorm = %g\n",
              workDataNorm, workDualNorm);
  return HighsDebugStatus::kOk;
}

HighsDebugStatus debugDualChuzcFailHeap(
    const HighsOptions& options, const HighsInt workCount,
    const std::vector<std::pair<HighsInt, double>>& workData,
    const HighsInt numVar, const double* workDual, const double selectTheta,
    const bool force) {
  // Non-trivially expensive assessment of CHUZC failure
  if (options.highs_debug_level < kHighsDebugLevelCostly && !force)
    return HighsDebugStatus::kNotChecked;

  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     No entries in heap so return error\n");
  double workDataNorm;
  double workDualNorm;
  debugDualChuzcFailNorms(workCount, workData, workDataNorm, numVar, workDual,
                          workDualNorm);
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     workCount = %" HIGHSINT_FORMAT
              "; selectTheta=%g\n",
              workCount, selectTheta);
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "DualChuzC:     workDataNorm = %g; workDualNorm = %g\n",
              workDataNorm, workDualNorm);
  return HighsDebugStatus::kOk;
}

HighsDebugStatus debugNonbasicFlagConsistent(const HighsOptions& options,
                                             const HighsLp& lp,
                                             const SimplexBasis& basis) {
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  HighsInt numTot = lp.num_col_ + lp.num_row_;
  const bool right_size = (HighsInt)basis.nonbasicFlag_.size() == numTot;
  if (!right_size) {
    highsLogDev(options.log_options, HighsLogType::kError,
                "nonbasicFlag size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  HighsInt num_basic_variables = 0;
  for (HighsInt var = 0; var < numTot; var++) {
    if (basis.nonbasicFlag_[var] == kNonbasicFlagFalse) {
      num_basic_variables++;
    } else {
      assert(basis.nonbasicFlag_[var] == kNonbasicFlagTrue);
    }
  }
  bool right_num_basic_variables = num_basic_variables == lp.num_row_;
  if (!right_num_basic_variables) {
    highsLogDev(options.log_options, HighsLogType::kError,
                "nonbasicFlag has %" HIGHSINT_FORMAT ", not %" HIGHSINT_FORMAT
                " basic variables\n",
                num_basic_variables, lp.num_row_);
    assert(right_num_basic_variables);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}
