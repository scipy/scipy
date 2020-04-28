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
    return_status = HighsDebugStatus::WARNING;
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
  if (have_primal_rhs && !primal_rhs_norm) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::WARNING,
                    "ComputePrimal: |PrimalRHS| = %9.4g", primal_rhs_norm);
    return_status = HighsDebugStatus::WARNING;
  }
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
    return_status = HighsDebugStatus::WARNING;
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
    return_status = HighsDebugStatus::WARNING;
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

HighsDebugStatus debugUpdatedObjectiveValue(
    HighsModelObject& highs_model_object, const SimplexAlgorithm algorithm,
    const int phase, const std::string message) {
  // Non-trivially expensive check of updated objective value. Computes the
  // exact objective value
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
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
      "UpdateObjVal:  %-9s large absolute (%9.4g) or relative (%9.4g) error in "
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
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  int num_free_variable_move_errors = 0;
  int num_lower_bounded_variable_move_errors = 0;
  int num_upper_bounded_variable_move_errors = 0;
  int num_boxed_variable_move_errors = 0;
  int num_fixed_variable_move_errors = 0;
  for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double lower = simplex_info.workLower_[iVar];
    const double upper = simplex_info.workUpper_[iVar];

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
    HighsPrintMessage(
        highs_model_object.options_.output,
        highs_model_object.options_.message_level, ML_ALWAYS,
        "There are %d nonbasicMove errors: %d free; %d lower; %d upper; %d "
        "boxed; %d fixed",
        num_errors, num_free_variable_move_errors,
        num_lower_bounded_variable_move_errors,
        num_upper_bounded_variable_move_errors, num_boxed_variable_move_errors,
        num_fixed_variable_move_errors);
  }
  assert(num_errors == 0);
  if (num_errors) return HighsDebugStatus::LOGICAL_ERROR;
  return HighsDebugStatus::OK;
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
    return_status = HighsDebugStatus::WARNING;
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
    return_status = HighsDebugStatus::WARNING;
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
