/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplexDebug.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "simplex/HSimplexDebug.h"

#include <string>

#include "lp_data/HighsDebug.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsSolutionDebug.h"
#include "simplex/HDualRow.h"
#include "simplex/HFactorDebug.h"
#include "simplex/HSimplex.h"
#include "simplex/SimplexTimer.h"

const double excessive_absolute_primal_norm = 1e12;
const double excessive_relative_primal_norm = 1e6;
const double large_absolute_primal_norm = sqrt(excessive_absolute_primal_norm);
const double large_relative_primal_norm = sqrt(excessive_relative_primal_norm);

const double excessive_absolute_nonbasic_dual_norm = 1e12;
const double excessive_relative_nonbasic_dual_norm = 1e6;
const double large_absolute_nonbasic_dual_norm =
    sqrt(excessive_absolute_nonbasic_dual_norm);
const double large_relative_nonbasic_dual_norm =
    sqrt(excessive_relative_nonbasic_dual_norm);

const double large_absolute_basic_dual_norm = 1e-12;
const double large_relative_basic_dual_norm = 1e-14;
const double excessive_absolute_basic_dual_norm =
    sqrt(large_absolute_basic_dual_norm);
const double excessive_relative_basic_dual_norm =
    sqrt(large_relative_basic_dual_norm);

const double computed_primal_excessive_absolute_norm =
    excessive_absolute_primal_norm;
const double computed_primal_excessive_relative_norm =
    excessive_relative_primal_norm;
const double computed_primal_large_absolute_norm = large_absolute_primal_norm;
const double computed_primal_large_relative_norm = large_relative_primal_norm;

const double computed_dual_excessive_absolute_nonbasic_dual_norm =
    excessive_absolute_nonbasic_dual_norm;
const double computed_dual_excessive_relative_nonbasic_dual_norm =
    excessive_relative_nonbasic_dual_norm;
const double computed_dual_large_absolute_nonbasic_dual_norm =
    large_absolute_nonbasic_dual_norm;
const double computed_dual_large_relative_nonbasic_dual_norm =
    large_relative_nonbasic_dual_norm;

const double computed_dual_excessive_absolute_basic_dual_norm =
    excessive_absolute_basic_dual_norm;
const double computed_dual_excessive_relative_basic_dual_norm =
    excessive_relative_basic_dual_norm;
const double computed_dual_large_absolute_basic_dual_norm =
    large_absolute_basic_dual_norm;
const double computed_dual_large_relative_basic_dual_norm =
    large_relative_basic_dual_norm;

const double computed_dual_small_relative_nonbasic_dual_change_norm = 1e-12;
const double computed_dual_large_relative_nonbasic_dual_change_norm =
    sqrt(computed_dual_small_relative_nonbasic_dual_change_norm);
const double computed_dual_small_absolute_nonbasic_dual_change_norm = 1e-6;
const double computed_dual_large_absolute_nonbasic_dual_change_norm =
    sqrt(computed_dual_small_absolute_nonbasic_dual_change_norm);

const double updated_objective_small_relative_error = 1e-12;
const double updated_objective_large_relative_error =
    sqrt(updated_objective_small_relative_error);
const double updated_objective_small_absolute_error = 1e-6;
const double updated_objective_large_absolute_error =
    sqrt(updated_objective_small_absolute_error);

const double excessive_basis_condition = 1e16;
const double large_basis_condition = sqrt(excessive_basis_condition);
const double fair_basis_condition = sqrt(large_basis_condition);

const double cleanup_large_absolute_nonbasic_dual_change_norm = 1e-12;
const double cleanup_large_relative_nonbasic_dual_change_norm = 1e-6;
const double cleanup_excessive_absolute_nonbasic_dual_change_norm =
    sqrt(cleanup_large_absolute_nonbasic_dual_change_norm);
const double cleanup_excessive_relative_nonbasic_dual_change_norm =
    sqrt(cleanup_large_relative_nonbasic_dual_change_norm);

const double freelist_excessive_pct_num_entries = 25.0;
const double freelist_large_pct_num_entries = 10.0;
const double freelist_fair_pct_num_entries = 1.0;

HighsDebugStatus debugSimplexLp(const HighsModelObject& highs_model_object) {
  // Non-trivially expensive check that the .simplex_lp, if valid is .lp scaled
  // according to .scale
  const HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  if (!simplex_lp_status.valid ||
      highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  const HighsOptions& options = highs_model_object.options_;
  const HighsLp& lp = highs_model_object.lp_;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsScale& scale = highs_model_object.scale_;
  const HFactor& factor = highs_model_object.factor_;

  bool right_size = true;
  right_size = (int)scale.col_.size() == lp.numCol_ && right_size;
  right_size = (int)scale.row_.size() == lp.numRow_ && right_size;
  if (!right_size) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "scale size error");
    assert(right_size);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  // Take a copy of the original LP
  HighsLp check_lp = lp;
  if (applyScalingToLp(options, check_lp, scale) != HighsStatus::OK) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "debugSimplexLp: Error scaling check LP");
    return HighsDebugStatus::LOGICAL_ERROR;
  }
  const bool simplex_lp_data_ok = check_lp == simplex_lp;
  if (!simplex_lp_data_ok) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "debugSimplexLp: Check LP and simplex LP not equal");
    assert(simplex_lp_data_ok);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }

  if (simplex_lp_status.has_basis) {
    const bool simplex_basis_correct =
        debugDebugToHighsStatus(debugSimplexBasisCorrect(highs_model_object)) !=
        HighsStatus::Error;
    if (!simplex_basis_correct) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Supposed to be a Simplex basis, but incorrect");
      assert(simplex_basis_correct);
      return_status = HighsDebugStatus::LOGICAL_ERROR;
    }
  }

  if (simplex_lp_status.has_invert) {
    const bool invert_ok = debugDebugToHighsStatus(debugCheckInvert(
                               options, factor)) != HighsStatus::Error;
    if (!invert_ok) {
      HighsLogMessage(
          options.logfile, HighsMessageType::ERROR,
          "Supposed to be a Simplex basis inverse, but too inaccurate");
      assert(invert_ok);
      return_status = HighsDebugStatus::LOGICAL_ERROR;
    }
  }
  return return_status;
}

HighsDebugStatus debugBasisRightSize(const HighsOptions& options,
                                     const HighsLp& simplex_lp,
                                     const SimplexBasis& simplex_basis) {
  // Cheap analysis of a Simplex basis, checking vector sizes
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  bool right_size = isBasisRightSize(simplex_lp, simplex_basis);
  if (!right_size) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Simplex basis size error");
    assert(right_size);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  return return_status;
}

HighsDebugStatus debugSimplexInfoBasisRightSize(
    const HighsModelObject& highs_model_object) {
  // Trivially cheap check of dimensions and sizes
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;

  const HighsOptions& options = highs_model_object.options_;
  const HighsLp& lp = highs_model_object.lp_;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;

  int numCol = lp.numCol_;
  int numRow = lp.numRow_;
  int numTot = numCol + numRow;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  bool dimension_ok =
      numCol == simplex_lp.numCol_ && numRow == simplex_lp.numRow_;
  assert(dimension_ok);
  if (!dimension_ok) {
    HighsLogMessage(
        options.logfile, HighsMessageType::ERROR,
        "LP-SimplexLP dimension incompatibility (%d, %d) != (%d, %d)", numCol,
        simplex_lp.numCol_, numRow, simplex_lp.numRow_);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  bool right_size = true;
  right_size = (int)simplex_info.workCost_.size() == numTot && right_size;
  right_size = (int)simplex_info.workDual_.size() == numTot && right_size;
  right_size = (int)simplex_info.workShift_.size() == numTot && right_size;
  right_size = (int)simplex_info.workLower_.size() == numTot && right_size;
  right_size = (int)simplex_info.workUpper_.size() == numTot && right_size;
  right_size = (int)simplex_info.workRange_.size() == numTot && right_size;
  right_size = (int)simplex_info.workValue_.size() == numTot && right_size;
  if (!right_size) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "simplex_info work vector size error");
    assert(right_size);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  if (debugBasisRightSize(options, simplex_lp, simplex_basis) !=
      HighsDebugStatus::OK)
    return_status = HighsDebugStatus::LOGICAL_ERROR;

  return return_status;
}

HighsDebugStatus debugComputePrimal(const HighsModelObject& highs_model_object,
                                    const std::vector<double>& primal_rhs) {
  // Non-trivially expensive analysis of computed primal values.
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::NOT_CHECKED;
  const std::vector<double>& primal_value =
      highs_model_object.simplex_info_.baseValue_;

  int num_row = highs_model_object.simplex_lp_.numRow_;

  // Use the size of the RHS to determine whether to use it
  const bool have_primal_rhs = (int)primal_rhs.size() == num_row;

  double primal_rhs_norm = 0;
  if (have_primal_rhs) {
    for (int iRow = 0; iRow < num_row; iRow++)
      primal_rhs_norm += fabs(primal_rhs[iRow]);
  }
  double computed_absolute_primal_norm = 0;
  for (int iRow = 0; iRow < num_row; iRow++)
    computed_absolute_primal_norm += fabs(primal_value[iRow]);

  std::string value_adjective;
  int report_level;
  return_status = HighsDebugStatus::OK;
  double computed_relative_primal_norm;
  if (primal_rhs_norm) {
    computed_relative_primal_norm =
        computed_absolute_primal_norm / primal_rhs_norm;
  } else {
    computed_relative_primal_norm = -1;
  }
  if (computed_relative_primal_norm > computed_primal_excessive_relative_norm ||
      computed_absolute_primal_norm > computed_primal_excessive_absolute_norm) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::ERROR;
  } else if (computed_relative_primal_norm >
                 computed_primal_large_relative_norm ||
             computed_absolute_primal_norm >
                 computed_primal_large_absolute_norm) {
    value_adjective = "Large";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::WARNING;
  } else {
    value_adjective = "SMALL";
    report_level = ML_VERBOSE;
  }
  HighsPrintMessage(
      highs_model_object.options_.output,
      highs_model_object.options_.message_level, report_level,
      "ComputePrimal: %-9s absolute (%9.4g) or relative (%9.4g) norm of "
      "primal values\n",
      value_adjective.c_str(), computed_absolute_primal_norm,
      computed_relative_primal_norm);
  return return_status;
}
HighsDebugStatus debugComputeDual(const HighsModelObject& highs_model_object,
                                  const std::vector<double>& previous_dual,
                                  const std::vector<double>& basic_costs,
                                  const std::vector<double>& row_dual) {
  // Non-trivially expensive analysis of computed dual values.
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::NOT_CHECKED;
  const std::vector<double>& new_dual =
      highs_model_object.simplex_info_.workDual_;

  int num_row = highs_model_object.simplex_lp_.numRow_;
  int num_col = highs_model_object.simplex_lp_.numCol_;

  const bool have_basic_costs = (int)basic_costs.size() == num_row;
  const bool have_row_dual = (int)row_dual.size() == num_row;
  const bool have_previous_dual =
      (int)previous_dual.size() == num_col + num_row;

  double basic_costs_norm = 0;
  if (have_basic_costs) {
    for (int iRow = 0; iRow < num_row; iRow++)
      basic_costs_norm += fabs(basic_costs[iRow]);
  }
  double row_dual_norm = 0;
  if (have_row_dual) {
    for (int iRow = 0; iRow < num_row; iRow++)
      row_dual_norm += fabs(row_dual[iRow]);
  }
  double computed_dual_absolute_basic_dual_norm = 0;
  double computed_dual_absolute_nonbasic_dual_norm = 0;
  for (int iVar = 0; iVar < num_row + num_col; iVar++) {
    if (!highs_model_object.simplex_basis_.nonbasicFlag_[iVar]) {
      computed_dual_absolute_basic_dual_norm += fabs(new_dual[iVar]);
      continue;
    }
    computed_dual_absolute_nonbasic_dual_norm += fabs(new_dual[iVar]);
  }
  std::string value_adjective;
  int report_level;
  return_status = HighsDebugStatus::OK;
  // Comment on the norm of the basic costs being zero
  if (have_basic_costs && !basic_costs_norm) {
    HighsLogMessage(
        highs_model_object.options_.logfile, HighsMessageType::WARNING,
        "ComputeDual:   basic cost norm is = %9.4g", basic_costs_norm);
    return_status = HighsDebugStatus::WARNING;
  }
  // Comment on the norm of the nonbasic duals being zero
  if (!computed_dual_absolute_nonbasic_dual_norm) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::WARNING,
                    "ComputeDual:   nonbasic dual norm is = %9.4g",
                    computed_dual_absolute_nonbasic_dual_norm);
    return_status = HighsDebugStatus::WARNING;
  }

  // Comment on the norm of basic duals (relative to the norm of the
  // basic costs) which, as c_B-BB^{-1}c_B, should be zero
  double computed_dual_relative_basic_dual_norm;
  if (basic_costs_norm) {
    computed_dual_relative_basic_dual_norm =
        computed_dual_absolute_basic_dual_norm / basic_costs_norm;
  } else {
    computed_dual_relative_basic_dual_norm = -1;
  }
  if (computed_dual_relative_basic_dual_norm >
          computed_dual_excessive_relative_basic_dual_norm ||
      computed_dual_absolute_basic_dual_norm >
          computed_dual_excessive_absolute_basic_dual_norm) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::ERROR;
  } else if (computed_dual_relative_basic_dual_norm >
                 computed_dual_large_relative_basic_dual_norm ||
             computed_dual_absolute_basic_dual_norm >
                 computed_dual_large_absolute_basic_dual_norm) {
    value_adjective = "Large";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::WARNING;
  } else {
    value_adjective = "OK";
    report_level = ML_VERBOSE;
  }
  HighsPrintMessage(
      highs_model_object.options_.output,
      highs_model_object.options_.message_level, report_level,
      "ComputeDual:   %-9s absolute (%9.4g) or relative (%9.4g) norm of "
      "   basic dual values\n",
      value_adjective.c_str(), computed_dual_absolute_basic_dual_norm,
      computed_dual_relative_basic_dual_norm);
  // Comment on the norm of nonbasic duals relative to the norm of the
  // basic costs
  double computed_dual_relative_nonbasic_dual_norm;
  if (basic_costs_norm) {
    computed_dual_relative_nonbasic_dual_norm =
        computed_dual_absolute_nonbasic_dual_norm / basic_costs_norm;
  } else {
    computed_dual_relative_nonbasic_dual_norm = -1;
  }
  if (computed_dual_relative_nonbasic_dual_norm >
          computed_dual_excessive_relative_nonbasic_dual_norm ||
      computed_dual_absolute_nonbasic_dual_norm >
          computed_dual_excessive_absolute_nonbasic_dual_norm) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::ERROR;
  } else if (computed_dual_relative_nonbasic_dual_norm >
                 computed_dual_large_relative_nonbasic_dual_norm ||
             computed_dual_absolute_nonbasic_dual_norm >
                 computed_dual_large_absolute_nonbasic_dual_norm) {
    value_adjective = "Large";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::WARNING;
  } else {
    value_adjective = "OK";
    report_level = ML_VERBOSE;
  }
  HighsPrintMessage(
      highs_model_object.options_.output,
      highs_model_object.options_.message_level, report_level,
      "ComputeDual:   %-9s absolute (%9.4g) or relative (%9.4g) norm of "
      "nonbasic dual values\n",
      value_adjective.c_str(), computed_dual_absolute_nonbasic_dual_norm,
      computed_dual_relative_nonbasic_dual_norm);
  double report_basic_costs_norm = -1;
  if (basic_costs_norm) report_basic_costs_norm = basic_costs_norm;
  double report_row_dual_norm = -1;
  if (row_dual_norm) report_row_dual_norm = row_dual_norm;
  HighsPrintMessage(highs_model_object.options_.output,
                    highs_model_object.options_.message_level, report_level,
                    "ComputeDual:   B.pi=c_B has |c_B| = %9.4g; |pi| = %9.4g; "
                    "|pi^TA-c| = [basic %9.4g; nonbasic %9.4g]\n",
                    report_basic_costs_norm, report_row_dual_norm,
                    computed_dual_absolute_basic_dual_norm,
                    computed_dual_absolute_nonbasic_dual_norm);
  if (have_previous_dual) {
    // Comment on the change in the dual values
    std::string change_adjective;
    double computed_dual_absolute_nonbasic_dual_change_norm = 0;
    for (int iVar = 0; iVar < num_row + num_col; iVar++) {
      if (!highs_model_object.simplex_basis_.nonbasicFlag_[iVar]) continue;
      computed_dual_absolute_nonbasic_dual_change_norm +=
          fabs(new_dual[iVar] - previous_dual[iVar]);
    }
    double computed_dual_relative_nonbasic_dual_change_norm;
    if (computed_dual_absolute_nonbasic_dual_norm) {
      computed_dual_relative_nonbasic_dual_change_norm =
          computed_dual_absolute_nonbasic_dual_change_norm /
          computed_dual_absolute_nonbasic_dual_norm;
    } else {
      computed_dual_relative_nonbasic_dual_change_norm = -1;
    }
    if (computed_dual_relative_nonbasic_dual_change_norm >
            computed_dual_large_relative_nonbasic_dual_change_norm ||
        computed_dual_absolute_nonbasic_dual_change_norm >
            computed_dual_large_absolute_nonbasic_dual_change_norm) {
      change_adjective = "Large";
      report_level = ML_ALWAYS;
      return_status = HighsDebugStatus::WARNING;
    } else if (computed_dual_relative_nonbasic_dual_change_norm >
                   computed_dual_small_relative_nonbasic_dual_change_norm ||
               computed_dual_absolute_nonbasic_dual_change_norm >
                   computed_dual_small_absolute_nonbasic_dual_change_norm) {
      change_adjective = "Small";
      report_level = ML_DETAILED;
      return_status = HighsDebugStatus::WARNING;
    } else {
      change_adjective = "OK";
      report_level = ML_VERBOSE;
    }
    HighsPrintMessage(highs_model_object.options_.output,
                      highs_model_object.options_.message_level, report_level,
                      "ComputeDual:   %-9s absolute (%9.4g) or relative "
                      "(%9.4g) nonbasic dual change\n",
                      change_adjective.c_str(),
                      computed_dual_absolute_nonbasic_dual_change_norm,
                      computed_dual_relative_nonbasic_dual_change_norm);
  }
  return return_status;
}

HighsDebugStatus debugSimplexDualFeasibility(
    const HighsModelObject& highs_model_object, const std::string message,
    const bool force) {
  // Non-trivially expensive check of dual feasibility.
  if (highs_model_object.options_.highs_debug_level <
          HIGHS_DEBUG_LEVEL_COSTLY &&
      !force)
    return HighsDebugStatus::NOT_CHECKED;
  if (force)
    HighsPrintMessage(highs_model_object.options_.output, 1, 1,
                      "SmplxDuFeas:   Forcing debug\n");

  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  double scaled_dual_feasibility_tolerance =
      highs_model_object.scaled_solution_params_.dual_feasibility_tolerance;

  int num_dual_infeasibilities = 0;
  double sum_dual_infeasibilities = 0;
  double max_dual_infeasibility = 0;
  for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double dual = simplex_info.workDual_[iVar];
    const double lower = simplex_info.workLower_[iVar];
    const double upper = simplex_info.workUpper_[iVar];
    double dual_infeasibility = 0;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(dual);
    } else {
      // Not free: any dual infeasibility is given by the dual value
      // signed by nonbasicMove
      dual_infeasibility = -simplex_basis.nonbasicMove_[iVar] * dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibilities++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
  }
  if (num_dual_infeasibilities) {
    HighsPrintMessage(highs_model_object.options_.output,
                      highs_model_object.options_.message_level, ML_ALWAYS,
                      "SmplxDuFeas:   num/max/sum simplex dual infeasibilities "
                      "= %d / %g / %g - %s\n",
                      num_dual_infeasibilities, max_dual_infeasibility,
                      sum_dual_infeasibilities, message.c_str());
    return HighsDebugStatus::LOGICAL_ERROR;
  }
  return HighsDebugStatus::OK;
}

HighsDebugStatus debugUpdatedObjectiveValue(
    HighsModelObject& highs_model_object, const SimplexAlgorithm algorithm,
    const int phase, const std::string message, const bool force) {
  // Non-trivially expensive check of updated objective value. Computes the
  // exact objective value
  if (highs_model_object.options_.highs_debug_level <
          HIGHS_DEBUG_LEVEL_COSTLY &&
      !force)
    return HighsDebugStatus::NOT_CHECKED;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;

  static bool have_previous_exact_primal_objective_value;
  static double previous_exact_primal_objective_value;
  static double previous_updated_primal_objective_value;
  static double updated_primal_objective_correction;

  static bool have_previous_exact_dual_objective_value;
  static double previous_exact_dual_objective_value;
  static double previous_updated_dual_objective_value;
  static double updated_dual_objective_correction;
  if (phase < 0) {
    if (algorithm == SimplexAlgorithm::PRIMAL) {
      have_previous_exact_primal_objective_value = false;
    } else {
      have_previous_exact_dual_objective_value = false;
    }
    return HighsDebugStatus::OK;
  }
  double exact_objective_value;
  double updated_objective_value;
  bool have_previous_exact_objective_value;
  // Assign values to prevent compiler warning
  double previous_exact_objective_value = 0;
  double previous_updated_objective_value = 0;
  double updated_objective_correction = 0;
  std::string algorithm_name;
  if (algorithm == SimplexAlgorithm::PRIMAL) {
    algorithm_name = "primal";
    have_previous_exact_objective_value =
        have_previous_exact_primal_objective_value;
    if (have_previous_exact_objective_value) {
      previous_exact_objective_value = previous_exact_primal_objective_value;
      previous_updated_objective_value =
          previous_updated_primal_objective_value;
      updated_objective_correction = updated_primal_objective_correction;
    }
    updated_objective_value = simplex_info.updated_primal_objective_value;
    // Save the current objective value so that it can be recovered
    // after calling computePrimalObjectiveValue
    double save_objective_value = simplex_info.primal_objective_value;
    computePrimalObjectiveValue(highs_model_object);
    exact_objective_value = simplex_info.primal_objective_value;
    simplex_info.primal_objective_value = save_objective_value;
  } else {
    algorithm_name = "dual";
    have_previous_exact_objective_value =
        have_previous_exact_dual_objective_value;
    if (have_previous_exact_objective_value) {
      previous_exact_objective_value = previous_exact_dual_objective_value;
      previous_updated_objective_value = previous_updated_dual_objective_value;
      updated_objective_correction = updated_dual_objective_correction;
    }
    updated_objective_value = simplex_info.updated_dual_objective_value;
    // Save the current objective value so that it can be recovered
    // after calling computeDualObjectiveValue
    double save_objective_value = simplex_info.dual_objective_value;
    computeDualObjectiveValue(highs_model_object, phase);
    exact_objective_value = simplex_info.dual_objective_value;
    simplex_info.dual_objective_value = save_objective_value;
  }
  double change_in_objective_value = 0;
  double change_in_updated_objective_value = 0;
  if (have_previous_exact_objective_value) {
    change_in_objective_value =
        exact_objective_value - previous_exact_objective_value;
    change_in_updated_objective_value =
        updated_objective_value - previous_updated_objective_value;
    updated_objective_value += updated_objective_correction;
  } else {
    updated_objective_correction = 0;
  }
  const double updated_objective_error =
      exact_objective_value - updated_objective_value;
  const double updated_objective_absolute_error = fabs(updated_objective_error);
  const double updated_objective_relative_error =
      updated_objective_absolute_error / max(1.0, fabs(exact_objective_value));
  updated_objective_correction += updated_objective_error;

  // Now update the records of previous objective value
  if (algorithm == SimplexAlgorithm::PRIMAL) {
    have_previous_exact_primal_objective_value = true;
    previous_exact_primal_objective_value = exact_objective_value;
    previous_updated_primal_objective_value = updated_objective_value;
    updated_primal_objective_correction = updated_objective_correction;
  } else {
    have_previous_exact_dual_objective_value = true;
    previous_exact_dual_objective_value = exact_objective_value;
    previous_updated_dual_objective_value = updated_objective_value;
    updated_dual_objective_correction = updated_objective_correction;
  }

  // Now analyse the error
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  std::string error_adjective;
  int report_level;
  bool at_least_small_error =
      updated_objective_relative_error >
          updated_objective_small_relative_error ||
      updated_objective_absolute_error > updated_objective_small_absolute_error;
  if (!at_least_small_error) return return_status;
  if (updated_objective_relative_error >
          updated_objective_large_relative_error ||
      updated_objective_absolute_error >
          updated_objective_large_absolute_error) {
    error_adjective = "Large";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::LARGE_ERROR;
  } else if (updated_objective_relative_error >
                 updated_objective_small_relative_error ||
             updated_objective_absolute_error >
                 updated_objective_small_absolute_error) {
    error_adjective = "Small";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::SMALL_ERROR;
  } else {
    error_adjective = "OK";
    report_level = ML_VERBOSE;
    return_status = HighsDebugStatus::OK;
  }
  HighsPrintMessage(
      highs_model_object.options_.output,
      highs_model_object.options_.message_level, report_level,
      "UpdateObjVal:  %-9s absolute (%9.4g) or relative (%9.4g) error in "
      "updated %s objective value"
      " - objective change - exact (%9.4g) updated (%9.4g) | %s\n",
      error_adjective.c_str(), updated_objective_error,
      updated_objective_relative_error, algorithm_name.c_str(),
      change_in_objective_value, change_in_updated_objective_value,
      message.c_str());
  return return_status;
}

HighsDebugStatus debugUpdatedObjectiveValue(
    const HighsModelObject& highs_model_object,
    const SimplexAlgorithm algorithm) {
  // Cheap check of updated objective value - assumes that the
  // objective value computed directly is correct, so only call after
  // this has been done
  if (highs_model_object.options_.highs_debug_level == HIGHS_DEBUG_LEVEL_NONE)
    return HighsDebugStatus::NOT_CHECKED;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  std::string algorithm_name = "dual";
  if (algorithm == SimplexAlgorithm::PRIMAL) algorithm_name = "primal";
  double exact_objective_value;
  double updated_objective_value;
  if (algorithm == SimplexAlgorithm::PRIMAL) {
    assert(highs_model_object.simplex_lp_status_.has_primal_objective_value);
    exact_objective_value = simplex_info.primal_objective_value;
    updated_objective_value = simplex_info.updated_primal_objective_value;
  } else {
    assert(highs_model_object.simplex_lp_status_.has_dual_objective_value);
    exact_objective_value = simplex_info.dual_objective_value;
    updated_objective_value = simplex_info.updated_dual_objective_value;
  }
  const double updated_objective_error =
      exact_objective_value - updated_objective_value;
  const double updated_objective_absolute_error = fabs(updated_objective_error);
  const double updated_objective_relative_error =
      updated_objective_absolute_error / max(1.0, fabs(exact_objective_value));

  // Now analyse the error
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  std::string error_adjective;
  int report_level;
  if (updated_objective_relative_error >
          updated_objective_large_relative_error ||
      updated_objective_absolute_error >
          updated_objective_large_absolute_error) {
    error_adjective = "Large";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::LARGE_ERROR;
  } else if (updated_objective_relative_error >
                 updated_objective_small_relative_error ||
             updated_objective_absolute_error >
                 updated_objective_small_absolute_error) {
    error_adjective = "Small";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::SMALL_ERROR;
  } else {
    error_adjective = "OK";
    report_level = ML_VERBOSE;
    return_status = HighsDebugStatus::OK;
  }
  HighsPrintMessage(highs_model_object.options_.output,
                    highs_model_object.options_.message_level, report_level,
                    "UpdateObjVal:  %-9s large absolute (%9.4g) or relative "
                    "(%9.4g) error in updated %s objective value\n",
                    error_adjective.c_str(), updated_objective_error,
                    updated_objective_relative_error, algorithm_name.c_str());
  return return_status;
}

HighsDebugStatus debugFixedNonbasicMove(
    const HighsModelObject& highs_model_object) {
  // Non-trivially expensive check of nonbasicMove for fixed variables
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  int num_fixed_variable_move_errors = 0;
  for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    if (simplex_info.workLower_[iVar] == simplex_info.workUpper_[iVar] &&
        simplex_basis.nonbasicMove_[iVar])
      num_fixed_variable_move_errors++;
  }
  assert(num_fixed_variable_move_errors == 0);
  if (num_fixed_variable_move_errors) {
    HighsPrintMessage(highs_model_object.options_.output,
                      highs_model_object.options_.message_level, ML_ALWAYS,
                      "There are %d fixed nonbasicMove errors",
                      num_fixed_variable_move_errors);
    return HighsDebugStatus::LOGICAL_ERROR;
  }
  return HighsDebugStatus::OK;
}

HighsDebugStatus debugNonbasicMove(const HighsModelObject& highs_model_object) {
  // Non-trivially expensive check of NonbasicMove
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  const HighsOptions& options = highs_model_object.options_;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  int num_free_variable_move_errors = 0;
  int num_lower_bounded_variable_move_errors = 0;
  int num_upper_bounded_variable_move_errors = 0;
  int num_boxed_variable_move_errors = 0;
  int num_fixed_variable_move_errors = 0;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  bool right_size = (int)simplex_basis.nonbasicMove_.size() == numTot;
  // Check consistency of nonbasicMove
  if (!right_size) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "nonbasicMove size error");
    assert(right_size);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  double lower;
  double upper;

  for (int iVar = 0; iVar < numTot; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic variable
    if (iVar < simplex_lp.numCol_) {
      lower = simplex_lp.colLower_[iVar];
      upper = simplex_lp.colUpper_[iVar];
    } else {
      int iRow = iVar - simplex_lp.numCol_;
      lower = -simplex_lp.rowUpper_[iRow];
      upper = -simplex_lp.rowLower_[iRow];
    }

    if (highs_isInfinity(upper)) {
      if (highs_isInfinity(-lower)) {
        // Free
        if (simplex_basis.nonbasicMove_[iVar]) {
          num_free_variable_move_errors++;
        }
      } else {
        // Only lower bounded
        if (simplex_basis.nonbasicMove_[iVar] != NONBASIC_MOVE_UP) {
          num_lower_bounded_variable_move_errors++;
        }
      }
    } else {
      if (highs_isInfinity(-lower)) {
        // Only upper bounded
        if (simplex_basis.nonbasicMove_[iVar] != NONBASIC_MOVE_DN) {
          num_upper_bounded_variable_move_errors++;
        }
      } else {
        // Boxed or fixed
        if (lower != upper) {
          // Boxed
          if (!simplex_basis.nonbasicMove_[iVar]) {
            num_boxed_variable_move_errors++;
          }
        } else {
          // Fixed
          if (simplex_basis.nonbasicMove_[iVar]) {
            num_fixed_variable_move_errors++;
          }
        }
      }
    }
  }
  int num_errors =
      num_free_variable_move_errors + num_lower_bounded_variable_move_errors +
      num_upper_bounded_variable_move_errors + num_boxed_variable_move_errors +
      num_fixed_variable_move_errors;

  if (num_errors) {
    HighsLogMessage(
        options.logfile, HighsMessageType::ERROR,
        "There are %d nonbasicMove errors: %d free; %d lower; %d upper; %d "
        "boxed; %d fixed",
        num_errors, num_free_variable_move_errors,
        num_lower_bounded_variable_move_errors,
        num_upper_bounded_variable_move_errors, num_boxed_variable_move_errors,
        num_fixed_variable_move_errors);
    assert(num_errors == 0);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  return return_status;
}

HighsDebugStatus debugBasisCondition(const HighsModelObject& highs_model_object,
                                     const std::string message) {
  // Non-trivially expensive assessment of basis condition
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  double basis_condition = computeBasisCondition(highs_model_object);
  std::string value_adjective;
  int report_level;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  if (basis_condition > excessive_basis_condition) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::ERROR;
  } else if (basis_condition > large_basis_condition) {
    value_adjective = "Large";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::WARNING;
  } else if (basis_condition > fair_basis_condition) {
    value_adjective = "Fair";
    report_level = ML_VERBOSE;
    return_status = HighsDebugStatus::OK;
  } else {
    value_adjective = "OK";
    report_level = ML_VERBOSE;
    return_status = HighsDebugStatus::OK;
  }
  HighsPrintMessage(
      highs_model_object.options_.output,
      highs_model_object.options_.message_level, report_level,
      "BasisCond:     %-9s basis condition estimate (%9.4g) - %s\n",
      value_adjective.c_str(), basis_condition, message.c_str());
  return return_status;
}

HighsDebugStatus debugCleanup(HighsModelObject& highs_model_object,
                              const std::vector<double>& original_dual) {
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
#ifdef HiGHSDEV
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;
#endif
  // Make sure that the original_dual has been set up
  assert((int)original_dual.size() == simplex_lp.numCol_ + simplex_lp.numRow_);
  const std::vector<double>& new_dual = simplex_info.workDual_;

  const double dual_feasibility_tolerance =
      highs_model_object.scaled_solution_params_.dual_feasibility_tolerance;
  int num_dual_sign_change = 0;
  double cleanup_absolute_nonbasic_dual_change_norm = 0;
  double cleanup_absolute_nonbasic_dual_norm = 0;
  for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    cleanup_absolute_nonbasic_dual_norm += std::fabs(new_dual[iVar]);
#ifdef HiGHSDEV
    const double nonbasic_dual_change =
        std::fabs(new_dual[iVar] - original_dual[iVar]);
    updateValueDistribution(nonbasic_dual_change,
                            analysis.cleanup_dual_change_distribution);
    cleanup_absolute_nonbasic_dual_change_norm += nonbasic_dual_change;
#endif
    const double max_dual =
        std::max(std::fabs(new_dual[iVar]), std::fabs(original_dual[iVar]));
    if (max_dual > dual_feasibility_tolerance &&
        new_dual[iVar] * original_dual[iVar] < 0)
      num_dual_sign_change++;
  }
  // Comment on the norm of the nonbasic duals being zero
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  if (!cleanup_absolute_nonbasic_dual_norm) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::WARNING,
                    "DualCleanup:   dual norm is = %9.4g",
                    cleanup_absolute_nonbasic_dual_norm);
    return_status = HighsDebugStatus::WARNING;
  }
  // Comment on the norm of the change being zero
  if (!cleanup_absolute_nonbasic_dual_change_norm) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::WARNING,
                    "DualCleanup:   dual norm is = %9.4g",
                    cleanup_absolute_nonbasic_dual_change_norm);
    return_status = HighsDebugStatus::WARNING;
  }
  double cleanup_relative_nonbasic_dual_change_norm;
  if (cleanup_absolute_nonbasic_dual_norm) {
    cleanup_relative_nonbasic_dual_change_norm =
        cleanup_absolute_nonbasic_dual_change_norm /
        cleanup_absolute_nonbasic_dual_norm;
  } else {
    cleanup_relative_nonbasic_dual_change_norm = -1;
  }
  std::string value_adjective;
  int report_level;
  if (cleanup_absolute_nonbasic_dual_change_norm >
          cleanup_excessive_absolute_nonbasic_dual_change_norm ||
      cleanup_relative_nonbasic_dual_change_norm >
          cleanup_excessive_relative_nonbasic_dual_change_norm) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::ERROR;
  } else if (cleanup_absolute_nonbasic_dual_change_norm >
                 cleanup_large_absolute_nonbasic_dual_change_norm ||
             cleanup_relative_nonbasic_dual_change_norm >
                 cleanup_large_relative_nonbasic_dual_change_norm) {
    value_adjective = "Large";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::WARNING;
  } else {
    value_adjective = "OK";
    report_level = ML_VERBOSE;
    return_status = HighsDebugStatus::OK;
  }
  HighsPrintMessage(
      highs_model_object.options_.output,
      highs_model_object.options_.message_level, report_level,
      "DualCleanup:   %-9s absolute (%9.4g) or relative (%9.4g) dual change, "
      "with %d meaningful sign change(s)\n",
      value_adjective.c_str(), cleanup_absolute_nonbasic_dual_change_norm,
      cleanup_relative_nonbasic_dual_change_norm, num_dual_sign_change);
  return return_status;
}

HighsDebugStatus debugFreeListNumEntries(
    const HighsModelObject& highs_model_object, const std::set<int>& freeList) {
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;

  int freelist_num_entries = 0;
  if (freeList.size() > 0) {
    std::set<int>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++)
      freelist_num_entries++;
  }

  const int numTot = highs_model_object.simplex_lp_.numCol_ +
                     highs_model_object.simplex_lp_.numRow_;
  double pct_freelist_num_entries = (100.0 * freelist_num_entries) / numTot;

  std::string value_adjective;
  int report_level;
  HighsDebugStatus return_status = HighsDebugStatus::NOT_CHECKED;

  if (pct_freelist_num_entries > freelist_excessive_pct_num_entries) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
  } else if (pct_freelist_num_entries > freelist_large_pct_num_entries) {
    value_adjective = "Large";
    report_level = ML_DETAILED;
  } else if (pct_freelist_num_entries > freelist_fair_pct_num_entries) {
    value_adjective = "Fair";
    report_level = ML_VERBOSE;
  } else {
    value_adjective = "OK";
    if (freelist_num_entries) {
      report_level = ML_ALWAYS;
    } else {
      report_level = ML_VERBOSE;
    }
    return_status = HighsDebugStatus::OK;
  }

  HighsPrintMessage(
      highs_model_object.options_.output,
      highs_model_object.options_.message_level, report_level,
      "FreeList   :   %-9s percentage (%6.2g) of %d variables on free list\n",
      value_adjective.c_str(), pct_freelist_num_entries, numTot);

  return return_status;
}

HighsDebugStatus debugDualChuzcFail(
    const HighsOptions& options, const int workCount,
    const std::vector<std::pair<int, double>>& workData, const double* workDual,
    const double selectTheta, const double remainTheta) {
  // Non-trivially expensive assessment of basis condition
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;

  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "DualChuzC:     No change in loop 2 so return error\n");
  double workDataNorm = 0;
  double dualNorm = 0;
  for (int i = 0; i < workCount; i++) {
    int iCol = workData[i].first;
    double value = workData[i].second;
    workDataNorm += value * value;
    value = workDual[iCol];
    dualNorm += value * value;
  }
  workDataNorm += sqrt(workDataNorm);
  dualNorm += sqrt(dualNorm);
  HighsPrintMessage(
      options.output, options.message_level, ML_ALWAYS,
      "DualChuzC:     workCount = %d; selectTheta=%g; remainTheta=%g\n",
      workCount, selectTheta, remainTheta);
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "DualChuzC:     workDataNorm = %g; dualNorm = %g\n",
                    workDataNorm, dualNorm);
  return HighsDebugStatus::OK;
}

void debugDualChuzcWorkDataAndGroupReport(
    const HighsModelObject& highs_model_object, const double workDelta,
    const double workTheta, const std::string message,
    const int report_workCount,
    const std::vector<std::pair<int, double>>& report_workData,
    const std::vector<int>& report_workGroup) {
  const HighsOptions& options = highs_model_object.options_;
  const std::vector<int>& workMove =
      highs_model_object.simplex_basis_.nonbasicMove_;
  const std::vector<double>& workDual =
      highs_model_object.simplex_info_.workDual_;
  const std::vector<double>& workRange =
      highs_model_object.simplex_info_.workRange_;
  const double Td =
      highs_model_object.scaled_solution_params_.dual_feasibility_tolerance;
  double totalChange = initial_total_change;
  const double totalDelta = fabs(workDelta);
  HighsPrintMessage(
      options.output, options.message_level, ML_ALWAYS,
      "\n%s: totalDelta = %10.4g\nworkData\n  En iCol       Dual      Value    "
      "  Ratio     Change\n",
      message.c_str(), totalDelta);
  for (int i = 0; i < report_workCount; i++) {
    int iCol = report_workData[i].first;
    double value = report_workData[i].second;
    double dual = workMove[iCol] * workDual[iCol];
    totalChange += value * (workRange[iCol]);
    HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                      "%4d %4d %10.4g %10.4g %10.4g %10.4g\n", i, iCol, dual,
                      value, dual / value, totalChange);
  }
  double selectTheta = workTheta;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "workGroup\n  Ix:   selectTheta Entries\n");
  for (int group = 0; group < (int)report_workGroup.size() - 1; group++) {
    HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                      "%4d: selectTheta = %10.4g ", group, selectTheta);
    for (int en = report_workGroup[group]; en < report_workGroup[group + 1];
         en++) {
      HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                        "%4d ", en);
    }
    HighsPrintMessage(options.output, options.message_level, ML_ALWAYS, "\n");
    int en = report_workGroup[group + 1];
    int iCol = report_workData[en].first;
    double value = report_workData[en].second;
    double dual = workMove[iCol] * workDual[iCol];
    selectTheta = (dual + Td) / value;
  }
}

HighsDebugStatus debugDualChuzcWorkDataAndGroup(
    const HighsModelObject& highs_model_object, const double workDelta,
    const double workTheta, const int workCount, const int alt_workCount,
    const int breakIndex, const int alt_breakIndex,
    const std::vector<std::pair<int, double>>& workData,
    const std::vector<std::pair<int, double>>& sorted_workData,
    const std::vector<int>& workGroup, const std::vector<int>& alt_workGroup) {
  // Cheap comparison and possible non-trivially expensive reporting
  // of the two sorting methods for BFRT nodes in dual CHUZC
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  const HighsOptions& options = highs_model_object.options_;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  int workPivot = workData[breakIndex].first;
  int alt_workPivot = sorted_workData[alt_breakIndex].first;
  if (alt_workPivot != workPivot) {
    HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                      "Quad workPivot = %d; Heap workPivot = %d\n", workPivot,
                      alt_workPivot);
    return_status = HighsDebugStatus::WARNING;
    if (highs_model_object.options_.highs_debug_level <
        HIGHS_DEBUG_LEVEL_COSTLY)
      return return_status;
    debugDualChuzcWorkDataAndGroupReport(highs_model_object, workDelta,
                                         workTheta, "Original", workCount,
                                         workData, workGroup);
    debugDualChuzcWorkDataAndGroupReport(
        highs_model_object, workDelta, workTheta, "Heap-derived", alt_workCount,
        sorted_workData, alt_workGroup);
  }
  return return_status;
}

HighsDebugStatus debugSimplexBasicSolution(
    const string message, const HighsModelObject& highs_model_object) {
  // Non-trivially expensive analysis of a simplex basic solution, starting from
  // solution_params
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;

  if (highsStatusFromHighsModelStatus(
          highs_model_object.scaled_model_status_) != HighsStatus::OK)
    return HighsDebugStatus::OK;
  HighsDebugStatus return_status = HighsDebugStatus::NOT_CHECKED;

  const HighsLp& lp = highs_model_object.lp_;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsScale& scale = highs_model_object.scale_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;

  return_status = debugSimplexInfoBasisRightSize(highs_model_object);
  if (return_status == HighsDebugStatus::LOGICAL_ERROR) return return_status;

  // Determine a HiGHS basis from the simplex basis. Only basic/nonbasic is
  // needed
  HighsBasis basis;
  basis.col_status.resize(lp.numCol_);
  basis.row_status.resize(lp.numRow_);
  // Now scatter the indices of basic variables
  for (int iVar = 0; iVar < lp.numCol_ + lp.numRow_; iVar++) {
    if (iVar < lp.numCol_) {
      int iCol = iVar;
      if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE) {
        basis.col_status[iCol] = HighsBasisStatus::NONBASIC;
      } else {
        basis.col_status[iCol] = HighsBasisStatus::BASIC;
      }
    } else {
      int iRow = iVar - lp.numCol_;
      if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE) {
        basis.row_status[iRow] = HighsBasisStatus::NONBASIC;
      } else {
        basis.row_status[iRow] = HighsBasisStatus::BASIC;
      }
    }
  }
  basis.valid_ = true;
  // Possibly scaled model
  // Determine a HiGHS solution simplex solution
  HighsSolution solution;
  solution.col_value.resize(lp.numCol_);
  solution.col_dual.resize(lp.numCol_);
  solution.row_value.resize(lp.numRow_);
  solution.row_dual.resize(lp.numRow_);

  for (int iVar = 0; iVar < lp.numCol_ + lp.numRow_; iVar++) {
    if (iVar < lp.numCol_) {
      int iCol = iVar;
      solution.col_value[iCol] = simplex_info.workValue_[iVar];
      solution.col_dual[iCol] =
          (int)simplex_lp.sense_ * simplex_info.workDual_[iVar];
    } else {
      int iRow = iVar - lp.numCol_;
      solution.row_value[iRow] = -simplex_info.workValue_[iVar];
      solution.row_dual[iRow] =
          (int)simplex_lp.sense_ * simplex_info.workDual_[iVar];
    }
  }
  // Now insert the basic values
  for (int ix = 0; ix < lp.numRow_; ix++) {
    int iVar = simplex_basis.basicIndex_[ix];
    if (iVar < lp.numCol_) {
      solution.col_value[iVar] = simplex_info.baseValue_[ix];
      solution.col_dual[iVar] = 0;
    } else {
      int iRow = iVar - lp.numCol_;
      solution.row_value[iRow] = -simplex_info.baseValue_[ix];
      solution.row_dual[iRow] = 0;
    }
  }

  const std::string message_scaled = message + " - scaled";
  return_status = debugWorseStatus(
      debugHighsBasicSolution(message_scaled, highs_model_object.options_,
                              simplex_lp, basis, solution,
                              highs_model_object.scaled_solution_params_,
                              highs_model_object.scaled_model_status_),
      return_status);

  if (!highs_model_object.scale_.is_scaled_) return return_status;

  // Doesn't work if simplex LP has permuted columns
  assert(!highs_model_object.simplex_lp_status_.is_permuted);
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    solution.col_value[iCol] *= scale.col_[iCol];
    solution.col_dual[iCol] /= (scale.col_[iCol] / scale.cost_);
  }
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
    solution.row_value[iRow] /= scale.row_[iRow];
    solution.row_dual[iRow] *= (scale.row_[iRow] * scale.cost_);
  }
  // Cannot assume unscaled solution params or unscaled model status are known
  const std::string message_unscaled = message + " - unscaled";
  return_status = debugWorseStatus(
      debugHighsBasicSolution(message_unscaled, highs_model_object.options_, lp,
                              basis, solution),
      return_status);

  // Scaled model
  return return_status;
}

HighsDebugStatus debugSimplexHighsSolutionDifferences(
    const HighsModelObject& highs_model_object) {
  // Nontrivially expensive check of dimensions and sizes
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;

  const HighsOptions& options = highs_model_object.options_;
  const HighsSolution& solution = highs_model_object.solution_;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  const HighsScale& scale = highs_model_object.scale_;

  HighsDebugStatus return_status = HighsDebugStatus::NOT_CHECKED;

  // Go through the columns, finding the differences in nonbasic column values
  // and duals
  double max_nonbasic_col_value_difference = 0;
  double max_nonbasic_col_dual_difference = 0;
  for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
    int iVar = iCol;
    if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE) {
      // Consider this nonbasic column
      double local_col_value = simplex_info.workValue_[iVar] * scale.col_[iCol];
      double local_col_dual = (int)simplex_lp.sense_ *
                              simplex_info.workDual_[iVar] /
                              (scale.col_[iCol] / scale.cost_);
      double value_difference =
          fabs(local_col_value - solution.col_value[iCol]);
      double dual_difference = fabs(local_col_dual - solution.col_dual[iCol]);
      max_nonbasic_col_value_difference =
          std::max(value_difference, max_nonbasic_col_value_difference);
      max_nonbasic_col_dual_difference =
          std::max(dual_difference, max_nonbasic_col_dual_difference);
    }
  }
  // Go through the rows, finding the differences in nonbasic and
  // basic row values and duals, as well as differences in basic
  // column values and duals
  double max_nonbasic_row_value_difference = 0;
  double max_nonbasic_row_dual_difference = 0;
  double max_basic_col_value_difference = 0;
  double max_basic_col_dual_difference = 0;
  double max_basic_row_value_difference = 0;
  double max_basic_row_dual_difference = 0;

  for (int ix = 0; ix < simplex_lp.numRow_; ix++) {
    int iRow = ix;
    int iVar = simplex_lp.numCol_ + iRow;
    if (simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE) {
      // Consider this nonbasic row
      double local_row_value =
          -simplex_info.workValue_[iVar] / scale.row_[iRow];
      double local_row_dual = (int)simplex_lp.sense_ *
                              simplex_info.workDual_[iVar] *
                              (scale.row_[iRow] * scale.cost_);
      double value_difference =
          fabs(local_row_value - solution.row_value[iRow]);
      double dual_difference = fabs(local_row_dual - solution.row_dual[iRow]);
      max_nonbasic_row_value_difference =
          std::max(value_difference, max_nonbasic_row_value_difference);
      max_nonbasic_row_dual_difference =
          std::max(dual_difference, max_nonbasic_row_dual_difference);
    }
    // Consider the basic variable associated with this row index
    iVar = simplex_basis.basicIndex_[ix];
    if (iVar < simplex_lp.numCol_) {
      // Consider this basic column
      int iCol = iVar;
      double local_col_value = simplex_info.baseValue_[ix] * scale.col_[iCol];
      double local_col_dual = 0;
      double value_difference =
          fabs(local_col_value - solution.col_value[iCol]);
      double dual_difference = fabs(local_col_dual - solution.col_dual[iCol]);
      max_basic_col_value_difference =
          std::max(value_difference, max_basic_col_value_difference);
      max_basic_col_dual_difference =
          std::max(dual_difference, max_basic_col_dual_difference);
    } else {
      // Consider this basic row
      iRow = iVar - simplex_lp.numCol_;
      double local_row_value = -simplex_info.baseValue_[ix] / scale.row_[iRow];
      double local_row_dual = 0;
      double value_difference =
          fabs(local_row_value - solution.row_value[iRow]);
      double dual_difference = fabs(local_row_dual - solution.row_dual[iRow]);
      max_basic_row_value_difference =
          std::max(value_difference, max_basic_row_value_difference);
      max_basic_row_dual_difference =
          std::max(dual_difference, max_basic_row_dual_difference);
    }
  }

  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "\nHiGHS-simplex solution differences\n");
  std::string value_adjective;
  int report_level;
  return_status = HighsDebugStatus::OK;
  if (max_nonbasic_col_value_difference > 0) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = debugWorseStatus(HighsDebugStatus::ERROR, return_status);
    HighsPrintMessage(
        options.output, options.message_level, report_level,
        "HighsSimplexD: %-9s Nonbasic column value difference: %9.4g\n",
        value_adjective.c_str(), max_nonbasic_col_value_difference);
  }
  if (max_nonbasic_row_value_difference > 0) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = debugWorseStatus(HighsDebugStatus::ERROR, return_status);
    HighsPrintMessage(
        options.output, options.message_level, report_level,
        "HighsSimplexD: %-9s Nonbasic row    value difference: %9.4g\n",
        value_adjective.c_str(), max_nonbasic_row_value_difference);
  }

  return_status = debugWorseStatus(
      debugAssessSolutionNormDifference(options, "Basic   column value",
                                        max_basic_col_value_difference),
      return_status);
  return_status = debugWorseStatus(
      debugAssessSolutionNormDifference(options, "Basic      row value",
                                        max_basic_row_value_difference),
      return_status);
  return_status = debugWorseStatus(
      debugAssessSolutionNormDifference(options, "Nonbasic column dual",
                                        max_nonbasic_col_dual_difference),
      return_status);
  return_status = debugWorseStatus(
      debugAssessSolutionNormDifference(options, "Nonbasic    row dual",
                                        max_nonbasic_row_dual_difference),
      return_status);

  if (max_basic_col_dual_difference > 0) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = debugWorseStatus(HighsDebugStatus::ERROR, return_status);
    HighsPrintMessage(
        options.output, options.message_level, report_level,
        "HighsSimplexD: %-9s Basic    column dual difference: %9.4g\n",
        value_adjective.c_str(), max_basic_col_dual_difference);
  }
  if (max_basic_row_dual_difference > 0) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = debugWorseStatus(HighsDebugStatus::ERROR, return_status);
    HighsPrintMessage(
        options.output, options.message_level, report_level,
        "HighsSimplexD: %-9s Basic    row     dual difference: %9.4g\n",
        value_adjective.c_str(), max_basic_row_dual_difference);
  }

  return return_status;
}

HighsDebugStatus debugAssessSolutionNormDifference(const HighsOptions& options,
                                                   const std::string type,
                                                   const double difference) {
  const double small_difference = 1e-12;
  const double large_difference = 1e-8;
  const double excessive_difference = 1e-4;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  if (difference <= small_difference) return return_status;
  std::string value_adjective;
  int report_level;

  if (difference > excessive_difference) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::ERROR;
  } else if (difference > large_difference) {
    value_adjective = "Large";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::WARNING;
  } else {
    value_adjective = "OK";
    report_level = ML_VERBOSE;
  }
  HighsPrintMessage(options.output, options.message_level, report_level,
                    "HighsSimplexD: %-9s %s difference: %9.4g\n",
                    value_adjective.c_str(), type.c_str(), difference);
  return return_status;
}

HighsDebugStatus debugSimplexBasisCorrect(
    const HighsModelObject& highs_model_object) {
  // Nontrivially expensive analysis of a Simplex basis, checking
  // consistency, and then correctness of nonbasicMove
  const HighsOptions& options = highs_model_object.options_;
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  const bool consistent =
      debugBasisConsistent(options, simplex_lp, simplex_basis) !=
      HighsDebugStatus::LOGICAL_ERROR;
  if (!consistent) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Supposed to be a Simplex basis, but not consistent");
    assert(consistent);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return return_status;
  const bool correct_nonbasicMove =
      debugNonbasicMove(highs_model_object) != HighsDebugStatus::LOGICAL_ERROR;
  if (!correct_nonbasicMove) {
    HighsLogMessage(
        options.logfile, HighsMessageType::ERROR,
        "Supposed to be a Simplex basis, but nonbasicMove is incorrect");
    assert(correct_nonbasicMove);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  return return_status;
}

HighsDebugStatus debugBasisConsistent(const HighsOptions& options,
                                      const HighsLp& simplex_lp,
                                      const SimplexBasis& simplex_basis) {
  // Cheap analysis of a Simplex basis, checking vector sizes, numbers
  // of basic/nonbasic variables and non-repetition of basic variables
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  // Check consistency of nonbasicFlag
  if (debugNonbasicFlagConsistent(options, simplex_lp, simplex_basis) ==
      HighsDebugStatus::LOGICAL_ERROR) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "nonbasicFlag inconsistent");
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  const bool right_size =
      (int)simplex_basis.basicIndex_.size() == simplex_lp.numRow_;
  // Check consistency of basicIndex
  if (!right_size) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "basicIndex size error");
    assert(right_size);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  // Use localNonbasicFlag so that duplicate entries in basicIndex can
  // be spotted
  vector<int> localNonbasicFlag = simplex_basis.nonbasicFlag_;
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
    int iCol = simplex_basis.basicIndex_[iRow];
    int flag = localNonbasicFlag[iCol];
    // Indicate that this column has been found in basicIndex
    localNonbasicFlag[iCol] = -1;
    if (flag) {
      // Nonzero value for localNonbasicFlag entry means that column is either
      if (flag == NONBASIC_FLAG_TRUE) {
        // Nonbasic...
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "Entry basicIndex_[%d] = %d is not basic", iRow, iCol);
      } else {
        // .. or is -1 since it has already been found in basicIndex
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "Entry basicIndex_[%d] = %d is already basic", iRow,
                        iCol);
        assert(flag == -1);
      }
      assert(!flag);
      return_status = HighsDebugStatus::LOGICAL_ERROR;
    }
  }
  return return_status;
}

HighsDebugStatus debugNonbasicFlagConsistent(
    const HighsOptions& options, const HighsLp& simplex_lp,
    const SimplexBasis& simplex_basis) {
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  const bool right_size = (int)simplex_basis.nonbasicFlag_.size() == numTot;
  if (!right_size) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "nonbasicFlag size error");
    assert(right_size);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  int num_basic_variables = 0;
  for (int var = 0; var < numTot; var++) {
    if (simplex_basis.nonbasicFlag_[var] == NONBASIC_FLAG_FALSE) {
      num_basic_variables++;
    } else {
      assert(simplex_basis.nonbasicFlag_[var] == NONBASIC_FLAG_TRUE);
    }
  }
  bool right_num_basic_variables = num_basic_variables == simplex_lp.numRow_;
  if (!right_num_basic_variables) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "nonbasicFlag has %d, not %d basic variables",
                    num_basic_variables, simplex_lp.numRow_);
    assert(right_num_basic_variables);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  return return_status;
}

HighsDebugStatus debugOkForSolve(const HighsModelObject& highs_model_object,
                                 const int phase) {
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  const HighsDebugStatus return_status = HighsDebugStatus::OK;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  const HighsOptions& options = highs_model_object.options_;
  bool ok;
  // Minimal check - just look at flags. This means we trust them!
  ok = simplex_lp_status.has_basis && simplex_lp_status.has_matrix_col_wise &&
       simplex_lp_status.has_matrix_row_wise &&
       simplex_lp_status.has_factor_arrays &&
       simplex_lp_status.has_dual_steepest_edge_weights &&
       simplex_lp_status.has_invert;
  if (!ok) {
    if (!simplex_lp_status.has_basis)
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Not OK to solve since simplex_lp_status.has_basis = %d",
                      simplex_lp_status.has_basis);
    if (!simplex_lp_status.has_matrix_col_wise)
      HighsLogMessage(
          options.logfile, HighsMessageType::ERROR,
          "Not OK to solve since simplex_lp_status.has_matrix_col_wise "
          "= %d",
          simplex_lp_status.has_matrix_col_wise);
    if (!simplex_lp_status.has_matrix_row_wise)
      HighsLogMessage(
          options.logfile, HighsMessageType::ERROR,
          "Not OK to solve since simplex_lp_status.has_matrix_row_wise "
          "= %d",
          simplex_lp_status.has_matrix_row_wise);
    //    if (!simplex_lp_status.has_factor_arrays)
    //      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
    //                  "Not OK to solve since
    //      simplex_lp_status.has_factor_arrays = %d",
    //             simplex_lp_status.has_factor_arrays);
    if (!simplex_lp_status.has_dual_steepest_edge_weights)
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Not OK to solve since "
                      "simplex_lp_status.has_dual_steepest_edge_weights = %d",
                      simplex_lp_status.has_dual_steepest_edge_weights);
    if (!simplex_lp_status.has_invert)
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Not OK to solve since simplex_lp_status.has_invert = %d",
                      simplex_lp_status.has_invert);
  }
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return return_status;
  // Basis and data check
  if (debugBasisConsistent(highs_model_object.options_, simplex_lp,
                           highs_model_object.simplex_basis_) ==
      HighsDebugStatus::LOGICAL_ERROR)
    return HighsDebugStatus::LOGICAL_ERROR;
  if (!debugWorkArraysOk(highs_model_object, phase))
    return HighsDebugStatus::LOGICAL_ERROR;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int var = 0; var < numTot; ++var) {
    if (simplex_basis.nonbasicFlag_[var]) {
      // Nonbasic variable
      if (!debugOneNonbasicMoveVsWorkArraysOk(highs_model_object, var))
        return HighsDebugStatus::LOGICAL_ERROR;
    }
  }
  return return_status;
}

// Methods below are not called externally

bool debugWorkArraysOk(const HighsModelObject& highs_model_object,
                       const int phase) {
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const HighsOptions& options = highs_model_object.options_;
  bool ok = true;
  // Only check phase 2 bounds: others will have been set by solve() so can be
  // trusted
  if (phase == 2) {
    for (int col = 0; col < simplex_lp.numCol_; ++col) {
      int var = col;
      if (!highs_isInfinity(-simplex_info.workLower_[var])) {
        ok = simplex_info.workLower_[var] == simplex_lp.colLower_[col];
        if (!ok) {
          HighsLogMessage(
              options.logfile, HighsMessageType::ERROR,
              "For col %d, simplex_info.workLower_ should be %g but is %g", col,
              simplex_lp.colLower_[col], simplex_info.workLower_[var]);
          return ok;
        }
      }
      if (!highs_isInfinity(simplex_info.workUpper_[var])) {
        ok = simplex_info.workUpper_[var] == simplex_lp.colUpper_[col];
        if (!ok) {
          HighsLogMessage(
              options.logfile, HighsMessageType::ERROR,
              "For col %d, simplex_info.workUpper_ should be %g but is %g", col,
              simplex_lp.colUpper_[col], simplex_info.workUpper_[var]);
          return ok;
        }
      }
    }
    for (int row = 0; row < simplex_lp.numRow_; ++row) {
      int var = simplex_lp.numCol_ + row;
      if (!highs_isInfinity(-simplex_info.workLower_[var])) {
        ok = simplex_info.workLower_[var] == -simplex_lp.rowUpper_[row];
        if (!ok) {
          HighsLogMessage(
              options.logfile, HighsMessageType::ERROR,
              "For row %d, simplex_info.workLower_ should be %g but is %g", row,
              -simplex_lp.rowUpper_[row], simplex_info.workLower_[var]);
          return ok;
        }
      }
      if (!highs_isInfinity(simplex_info.workUpper_[var])) {
        ok = simplex_info.workUpper_[var] == -simplex_lp.rowLower_[row];
        if (!ok) {
          HighsLogMessage(
              options.logfile, HighsMessageType::ERROR,
              "For row %d, simplex_info.workUpper_ should be %g but is %g", row,
              -simplex_lp.rowLower_[row], simplex_info.workUpper_[var]);
          return ok;
        }
      }
    }
  }
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int var = 0; var < numTot; ++var) {
    ok = simplex_info.workRange_[var] ==
         (simplex_info.workUpper_[var] - simplex_info.workLower_[var]);
    if (!ok) {
      HighsLogMessage(
          options.logfile, HighsMessageType::ERROR,
          "For variable %d, simplex_info.workRange_ should be %g = %g - %g "
          "but is %g",
          var, simplex_info.workUpper_[var] - simplex_info.workLower_[var],
          simplex_info.workUpper_[var], simplex_info.workLower_[var],
          simplex_info.workRange_[var]);
      return ok;
    }
  }
  // Don't check perturbed costs: these will have been set by solve() so can be
  // trusted
  if (!simplex_info.costs_perturbed) {
    for (int col = 0; col < simplex_lp.numCol_; ++col) {
      int var = col;
      ok = simplex_info.workCost_[var] ==
           (int)simplex_lp.sense_ * simplex_lp.colCost_[col];
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "For col %d, simplex_info.workLower_ should be %g but is %g", col,
            simplex_lp.colLower_[col], simplex_info.workCost_[var]);
        return ok;
      }
    }
    for (int row = 0; row < simplex_lp.numRow_; ++row) {
      int var = simplex_lp.numCol_ + row;
      ok = simplex_info.workCost_[var] == 0.;
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "For row %d, simplex_info.workCost_ should be zero but is %g", row,
            simplex_info.workCost_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool debugOneNonbasicMoveVsWorkArraysOk(
    const HighsModelObject& highs_model_object, const int var) {
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  const HighsOptions& options = highs_model_object.options_;
  assert(var >= 0);
  assert(var < simplex_lp.numCol_ + simplex_lp.numRow_);
  // Make sure we're not checking a basic variable
  if (!simplex_basis.nonbasicFlag_[var]) return true;
  bool ok;
  if (!highs_isInfinity(-simplex_info.workLower_[var])) {
    if (!highs_isInfinity(simplex_info.workUpper_[var])) {
      // Finite lower and upper bounds so nonbasic move depends on whether they
      // are equal
      if (simplex_info.workLower_[var] == simplex_info.workUpper_[var]) {
        // Fixed variable
        ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE;
        if (!ok) {
          HighsLogMessage(
              options.logfile, HighsMessageType::ERROR,
              "Fixed variable %d (simplex_lp.numCol_ = %d) [%11g, %11g, "
              "%11g] so nonbasic "
              "move should be zero but is %d",
              var, simplex_lp.numCol_, simplex_info.workLower_[var],
              simplex_info.workValue_[var], simplex_info.workUpper_[var],
              simplex_basis.nonbasicMove_[var]);
          return ok;
        }
        ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
        if (!ok) {
          HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                          "Fixed variable %d (simplex_lp.numCol_ = %d) so "
                          "simplex_info.work value should be %g but "
                          "is %g",
                          var, simplex_lp.numCol_, simplex_info.workLower_[var],
                          simplex_info.workValue_[var]);
          return ok;
        }
      } else {
        // Boxed variable
        ok = (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) ||
             (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_DN);
        if (!ok) {
          HighsLogMessage(
              options.logfile, HighsMessageType::ERROR,
              "Boxed variable %d (simplex_lp.numCol_ = %d) [%11g, %11g, "
              "%11g] range %g so "
              "nonbasic move should be up/down but is  %d",
              var, simplex_lp.numCol_, simplex_info.workLower_[var],
              simplex_info.workValue_[var], simplex_info.workUpper_[var],
              simplex_info.workUpper_[var] - simplex_info.workLower_[var],
              simplex_basis.nonbasicMove_[var]);
          return ok;
        }
        if (simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
          ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
          if (!ok) {
            HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                            "Boxed variable %d (simplex_lp.numCol_ = %d) with "
                            "NONBASIC_MOVE_UP so work "
                            "value should be %g but is %g",
                            var, simplex_lp.numCol_,
                            simplex_info.workLower_[var],
                            simplex_info.workValue_[var]);
            return ok;
          }
        } else {
          ok = simplex_info.workValue_[var] == simplex_info.workUpper_[var];
          if (!ok) {
            HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                            "Boxed variable %d (simplex_lp.numCol_ = %d) with "
                            "NONBASIC_MOVE_DN so work "
                            "value should be %g but is %g",
                            var, simplex_lp.numCol_,
                            simplex_info.workUpper_[var],
                            simplex_info.workValue_[var]);
            return ok;
          }
        }
      }
    } else {
      // Infinite upper bound
      ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_UP;
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Finite lower bound and infinite upper bound variable %d "
            "(simplex_lp.numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be up=%2d but is  "
            "%d",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var], simplex_info.workUpper_[var],
            NONBASIC_MOVE_UP, simplex_basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == simplex_info.workLower_[var];
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Finite lower bound and infinite upper bound variable %d "
            "(simplex_lp.numCol_ = "
            "%d) so work value should be %g but is %g",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var]);
        return ok;
      }
    }
  } else {
    // Infinite lower bound
    if (!highs_isInfinity(simplex_info.workUpper_[var])) {
      ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_DN;
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Finite upper bound and infinite lower bound variable %d "
            "(simplex_lp.numCol_ = "
            "%d) [%11g, %11g, %11g] so nonbasic move should be down but is  "
            "%d",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var], simplex_info.workUpper_[var],
            simplex_basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == simplex_info.workUpper_[var];
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Finite upper bound and infinite lower bound variable %d "
            "(simplex_lp.numCol_ = "
            "%d) so work value should be %g but is %g",
            var, simplex_lp.numCol_, simplex_info.workUpper_[var],
            simplex_info.workValue_[var]);
        return ok;
      }
    } else {
      // Infinite upper bound
      ok = simplex_basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE;
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Free variable %d (simplex_lp.numCol_ = %d) [%11g, %11g, %11g] "
            "so nonbasic "
            "move should be zero but is  %d",
            var, simplex_lp.numCol_, simplex_info.workLower_[var],
            simplex_info.workValue_[var], simplex_info.workUpper_[var],
            simplex_basis.nonbasicMove_[var]);
        return ok;
      }
      ok = simplex_info.workValue_[var] == 0.0;
      if (!ok) {
        HighsLogMessage(
            options.logfile, HighsMessageType::ERROR,
            "Free variable %d (simplex_lp.numCol_ = %d) so work value should "
            "be zero but "
            "is %g",
            var, simplex_lp.numCol_, simplex_info.workValue_[var]);
        return ok;
      }
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

bool debugAllNonbasicMoveVsWorkArraysOk(
    const HighsModelObject& highs_model_object) {
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  //    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  const HighsOptions& options = highs_model_object.options_;
  bool ok;
  const int numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  for (int var = 0; var < numTot; ++var) {
    HighsLogMessage(
        options.logfile, HighsMessageType::ERROR,
        "NonbasicMoveVsWorkArrays: var = %2d; simplex_basis.nonbasicFlag_[var] "
        "= %2d",
        var, simplex_basis.nonbasicFlag_[var]);
    if (!simplex_basis.nonbasicFlag_[var]) continue;
    ok = debugOneNonbasicMoveVsWorkArraysOk(highs_model_object, var);
    if (!ok) {
      HighsLogMessage(
          options.logfile, HighsMessageType::ERROR,
          "Error in NonbasicMoveVsWorkArrays for nonbasic variable %d", var);
      assert(ok);
      return ok;
    }
  }
  // ok must be true if we reach here
  assert(ok);
  return ok;
}

void debugReportReinvertOnNumericalTrouble(
    const std::string method_name, const HighsModelObject& highs_model_object,
    const double numerical_trouble_measure, const double alpha_from_col,
    const double alpha_from_row, const double numerical_trouble_tolerance,
    const bool reinvert) {
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return;
  const double abs_alpha_from_col = fabs(alpha_from_col);
  const double abs_alpha_from_row = fabs(alpha_from_row);
  const double abs_alpha_diff = fabs(abs_alpha_from_col - abs_alpha_from_row);
  const int iteration_count = highs_model_object.iteration_counts_.simplex;
  const int update_count = highs_model_object.simplex_info_.update_count;
  const std::string model_name = highs_model_object.simplex_lp_.model_name_;

  const bool numerical_trouble =
      numerical_trouble_measure > numerical_trouble_tolerance;
  const bool near_numerical_trouble =
      10 * numerical_trouble_measure > numerical_trouble_tolerance;

  const bool wrong_sign = alpha_from_col * alpha_from_row <= 0;
  if (!near_numerical_trouble && !wrong_sign) return;
  std::string adjective;
  if (numerical_trouble) {
    adjective = "       exceeds";
  } else if (near_numerical_trouble) {
    adjective = "almost exceeds";
  } else {
    adjective = "clearly satisfies";
  }
  HighsLogMessage(highs_model_object.options_.logfile,
                  HighsMessageType::WARNING,
                  "%s (%s) [Iter %d; Update %d] Col: %11.4g; Row: %11.4g; Diff "
                  "= %11.4g: Measure %11.4g %s %11.4g",
                  method_name.c_str(), model_name.c_str(), iteration_count,
                  update_count, abs_alpha_from_col, abs_alpha_from_row,
                  abs_alpha_diff, numerical_trouble_measure, adjective.c_str(),
                  numerical_trouble_tolerance);
  if (wrong_sign) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::WARNING,
                    "   Incompatible signs for Col: %11.4g and Row: %11.4g",
                    alpha_from_col, alpha_from_row);
  }
  if ((numerical_trouble || wrong_sign) && !reinvert) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::WARNING,
                    "   Numerical trouble or wrong sign and not reinverting");
  }
}
