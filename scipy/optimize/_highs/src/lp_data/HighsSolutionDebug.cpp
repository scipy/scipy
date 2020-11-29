/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolutionDebug.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsSolutionDebug.h"

#include <math.h>

#include <vector>

#include "lp_data/HighsDebug.h"
#include "lp_data/HighsModelUtils.h"
#include "util/HighsUtils.h"

const double large_relative_solution_param_error = 1e-12;
const double excessive_relative_solution_param_error =
    sqrt(large_relative_solution_param_error);

const double large_residual_error = 1e-12;
const double excessive_residual_error = sqrt(large_residual_error);

HighsDebugStatus debugBasisConsistent(const HighsOptions& options,
                                      const HighsLp lp,
                                      const HighsBasis& basis) {
  // Cheap analysis of a HiGHS basis, checking vector sizes, numbers
  // of basic/nonbasic variables
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  if (!basis.valid_) return return_status;
  bool consistent = isBasisConsistent(lp, basis);
  if (!consistent) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "HiGHS basis inconsistency");
    assert(consistent);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  return return_status;
}

HighsDebugStatus debugBasisRightSize(const HighsOptions& options,
                                     const HighsLp lp,
                                     const HighsBasis& basis) {
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  bool right_size = isBasisRightSize(lp, basis);
  if (!right_size) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "HiGHS basis size error");
    assert(right_size);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  return return_status;
}

HighsDebugStatus debugSolutionRightSize(const HighsOptions& options,
                                        const HighsLp lp,
                                        const HighsSolution& solution) {
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  bool right_size = isSolutionRightSize(lp, solution);
  if (!right_size) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "HiGHS solution size error");
    assert(right_size);
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  }
  return return_status;
}

HighsDebugStatus debugHighsBasicSolution(
    const string message, const HighsModelObject& highs_model_object) {
  // Non-trivially expensive analysis of a HiGHS basic solution, starting from
  // highs_model_object
  return debugHighsBasicSolution(
      message, highs_model_object.options_, highs_model_object.lp_,
      highs_model_object.basis_, highs_model_object.solution_,
      highs_model_object.unscaled_solution_params_,
      highs_model_object.unscaled_model_status_);
}

HighsDebugStatus debugHighsBasicSolution(
    const string message, const HighsOptions& options, const HighsLp& lp,
    const HighsBasis& basis, const HighsSolution& solution,
    const HighsInfo& info, const HighsModelStatus model_status) {
  // Non-trivially expensive analysis of a HiGHS basic solution, starting from
  // options and info
  //
  // Extract the solution_params from info and options
  HighsSolutionParams solution_params;
  solution_params.primal_feasibility_tolerance =
      options.primal_feasibility_tolerance;
  solution_params.dual_feasibility_tolerance =
      options.dual_feasibility_tolerance;
  solution_params.primal_status = info.primal_status;
  solution_params.dual_status = info.dual_status;
  solution_params.objective_function_value = info.objective_function_value;
  solution_params.num_primal_infeasibilities = info.num_primal_infeasibilities;
  solution_params.max_primal_infeasibility = info.max_primal_infeasibility;
  solution_params.sum_primal_infeasibilities = info.sum_primal_infeasibilities;
  solution_params.num_dual_infeasibilities = info.num_dual_infeasibilities;
  solution_params.max_dual_infeasibility = info.max_dual_infeasibility;
  solution_params.sum_dual_infeasibilities = info.sum_dual_infeasibilities;

  return debugHighsBasicSolution(message, options, lp, basis, solution,
                                 solution_params, model_status);
}

HighsDebugStatus debugHighsBasicSolution(const string message,
                                         const HighsOptions& options,
                                         const HighsLp& lp,
                                         const HighsBasis& basis,
                                         const HighsSolution& solution) {
  // Non-trivially expensive analysis of a HiGHS basic solution, starting from
  // options, assuming no knowledge of solution parameters or model status
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;

  // Check that there is a solution and valid basis to use
  if (debugHaveBasisAndSolutionData(lp, basis, solution) !=
      HighsDebugStatus::OK)
    return HighsDebugStatus::LOGICAL_ERROR;

  // Extract the solution_params from options
  HighsSolutionParams solution_params;
  solution_params.primal_feasibility_tolerance =
      options.primal_feasibility_tolerance;
  solution_params.dual_feasibility_tolerance =
      options.dual_feasibility_tolerance;

  double check_primal_objective_value;
  double check_dual_objective_value;
  // Get values for solution params from scratch. Also get primal/dual errors
  HighsPrimalDualErrors primal_dual_errors;
  // Get the primal and dual infeasibilities and errors
  debugHighsBasicSolutionPrimalDualInfeasibilitiesAndErrors(
      options, lp, basis, solution, check_primal_objective_value,
      check_dual_objective_value, solution_params, primal_dual_errors);

  HighsModelStatus model_status = HighsModelStatus::NOTSET;
  if (solution_params.num_primal_infeasibilities == 0 &&
      solution_params.num_dual_infeasibilities == 0)
    model_status = HighsModelStatus::OPTIMAL;

  debugReportHighsBasicSolution(message, options, solution_params,
                                model_status);
  return debugAnalysePrimalDualErrors(options, primal_dual_errors);
}

HighsDebugStatus debugHighsBasicSolution(
    const string message, const HighsOptions& options, const HighsLp& lp,
    const HighsBasis& basis, const HighsSolution& solution,
    const HighsSolutionParams& solution_params,
    const HighsModelStatus model_status) {
  // Non-trivially expensive analysis of a HiGHS basic solution, starting from
  // solution_params
  if (options.highs_debug_level < HIGHS_DEBUG_LEVEL_CHEAP)
    return HighsDebugStatus::NOT_CHECKED;
  // No basis to test if model status corresponds to warning or error
  if (highsStatusFromHighsModelStatus(model_status) != HighsStatus::OK)
    return HighsDebugStatus::OK;

  // No basis to test if model status is primal infeasible or unbounded
  if (model_status == HighsModelStatus::PRIMAL_INFEASIBLE ||
      model_status == HighsModelStatus::PRIMAL_UNBOUNDED)
    return HighsDebugStatus::OK;

  // Check that there is a solution and valid basis to use
  if (debugHaveBasisAndSolutionData(lp, basis, solution) !=
      HighsDebugStatus::OK)
    return HighsDebugStatus::LOGICAL_ERROR;

  HighsSolutionParams check_solution_params;
  double check_primal_objective_value;
  double check_dual_objective_value;
  // Extract the primal and dual feasibility tolerances and solution status
  check_solution_params.primal_feasibility_tolerance =
      solution_params.primal_feasibility_tolerance;
  check_solution_params.dual_feasibility_tolerance =
      solution_params.dual_feasibility_tolerance;
  check_solution_params.primal_status = solution_params.primal_status;
  check_solution_params.dual_status = solution_params.dual_status;
  // Get values for solution params from scratch. Also get primal/dual errors
  HighsPrimalDualErrors primal_dual_errors;
  // Get the primal and dual infeasibilities and errors
  debugHighsBasicSolutionPrimalDualInfeasibilitiesAndErrors(
      options, lp, basis, solution, check_primal_objective_value,
      check_dual_objective_value, check_solution_params, primal_dual_errors);
  check_solution_params.objective_function_value = check_primal_objective_value;

  HighsDebugStatus return_status = debugCompareSolutionParams(
      options, solution_params, check_solution_params);
  debugReportHighsBasicSolution(message, options, solution_params,
                                model_status);
  return_status = debugWorseStatus(
      debugAnalysePrimalDualErrors(options, primal_dual_errors), return_status);

  return return_status;
}

// Methods below are not called externally

HighsDebugStatus debugHaveBasisAndSolutionData(const HighsLp& lp,
                                               const HighsBasis& basis,
                                               const HighsSolution& solution) {
  if (!isSolutionRightSize(lp, solution))
    return HighsDebugStatus::LOGICAL_ERROR;
  if (!isBasisRightSize(lp, basis) && basis.valid_)
    return HighsDebugStatus::LOGICAL_ERROR;
  return HighsDebugStatus::OK;
}

void debugHighsBasicSolutionPrimalDualInfeasibilitiesAndErrors(
    const HighsOptions& options, const HighsLp& lp, const HighsBasis& basis,
    const HighsSolution& solution, double& primal_objective_value,
    double& dual_objective_value, HighsSolutionParams& solution_params,
    HighsPrimalDualErrors& primal_dual_errors) {
  double primal_feasibility_tolerance =
      solution_params.primal_feasibility_tolerance;
  double dual_feasibility_tolerance =
      solution_params.dual_feasibility_tolerance;

  // solution_params are the values computed in this method.
  int& num_primal_infeasibilities = solution_params.num_primal_infeasibilities;
  double& max_primal_infeasibility = solution_params.max_primal_infeasibility;
  double& sum_primal_infeasibilities =
      solution_params.sum_primal_infeasibilities;
  int& num_dual_infeasibilities = solution_params.num_dual_infeasibilities;
  double& max_dual_infeasibility = solution_params.max_dual_infeasibility;
  double& sum_dual_infeasibilities = solution_params.sum_dual_infeasibilities;

  num_primal_infeasibilities = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibilities = 0;
  num_dual_infeasibilities = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibilities = 0;

  std::vector<double> primal_activities;
  std::vector<double> dual_activities;
  primal_activities.assign(lp.numRow_, 0);
  dual_activities.resize(lp.numCol_);
  int num_non_basic_var = 0;
  int num_basic_var = 0;

  int& num_nonzero_basic_duals = primal_dual_errors.num_nonzero_basic_duals;
  int& num_large_nonzero_basic_duals =
      primal_dual_errors.num_large_nonzero_basic_duals;
  double& max_nonzero_basic_dual = primal_dual_errors.max_nonzero_basic_dual;
  double& sum_nonzero_basic_duals = primal_dual_errors.sum_nonzero_basic_duals;

  int& num_off_bound_nonbasic = primal_dual_errors.num_off_bound_nonbasic;
  double& max_off_bound_nonbasic = primal_dual_errors.max_off_bound_nonbasic;
  double& sum_off_bound_nonbasic = primal_dual_errors.sum_off_bound_nonbasic;

  int& num_primal_residual = primal_dual_errors.num_primal_residual;
  double& max_primal_residual = primal_dual_errors.max_primal_residual;
  double& sum_primal_residual = primal_dual_errors.sum_primal_residual;

  int& num_dual_residual = primal_dual_errors.num_dual_residual;
  double& max_dual_residual = primal_dual_errors.max_dual_residual;
  double& sum_dual_residual = primal_dual_errors.sum_dual_residual;

  num_nonzero_basic_duals = 0;
  num_large_nonzero_basic_duals = 0;
  max_nonzero_basic_dual = 0;
  sum_nonzero_basic_duals = 0;

  num_off_bound_nonbasic = 0;
  max_off_bound_nonbasic = 0;
  sum_off_bound_nonbasic = 0;
  num_primal_residual = 0;
  max_primal_residual = 0;
  sum_primal_residual = 0;
  num_dual_residual = 0;
  max_dual_residual = 0;
  sum_dual_residual = 0;

  // Initialise the objective value calculations. Done using
  // HighsSolution so offset is vanilla
  primal_objective_value = lp.offset_;
  dual_objective_value = lp.offset_;

  bool header_written = false;
  double off_bound_nonbasic;
  double primal_infeasibility;
  double dual_infeasibility;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    double lower = lp.colLower_[iCol];
    double upper = lp.colUpper_[iCol];
    double value = solution.col_value[iCol];
    double dual = solution.col_dual[iCol];
    HighsBasisStatus status = basis.col_status[iCol];
    primal_objective_value += lp.colCost_[iCol] * value;
    if (status != HighsBasisStatus::BASIC) dual_objective_value += value * dual;
    // Flip dual according to lp.sense_
    dual *= (int)lp.sense_;
    bool report = false;
    bool query = debugBasicSolutionVariable(
        report, primal_feasibility_tolerance, dual_feasibility_tolerance,
        status, lower, upper, value, dual, num_non_basic_var, num_basic_var,
        off_bound_nonbasic, primal_infeasibility, dual_infeasibility);
    if (off_bound_nonbasic > 0) num_off_bound_nonbasic++;
    max_off_bound_nonbasic =
        std::max(off_bound_nonbasic, max_off_bound_nonbasic);
    sum_off_bound_nonbasic += off_bound_nonbasic;
    if (primal_infeasibility > primal_feasibility_tolerance)
      num_primal_infeasibilities++;
    max_primal_infeasibility =
        std::max(primal_infeasibility, max_primal_infeasibility);
    sum_primal_infeasibilities += primal_infeasibility;
    if (status == HighsBasisStatus::BASIC) {
      double abs_basic_dual = dual_infeasibility;
      if (abs_basic_dual > 0) {
        num_nonzero_basic_duals++;
        if (abs_basic_dual > dual_feasibility_tolerance)
          num_large_nonzero_basic_duals++;
        max_nonzero_basic_dual =
            std::max(abs_basic_dual, max_nonzero_basic_dual);
        sum_nonzero_basic_duals += abs_basic_dual;
      }
    } else {
      if (dual_infeasibility > dual_feasibility_tolerance)
        num_dual_infeasibilities++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
    report =
        options.highs_debug_level > HIGHS_DEBUG_LEVEL_EXPENSIVE ||
        (options.highs_debug_level == HIGHS_DEBUG_LEVEL_EXPENSIVE && query);
    if (report) {
      if (!header_written) {
        printf(
            "\nColumns\nIndex NonBs Mv [          LB,           UB]       "
            "Primal         Dual    PrimalIfs      DualIfs\n");
        header_written = true;
      }
      printf("%5d %5d [%12g, %12g] %12g %12g", iCol, (int)status, lower, upper,
             value, dual);
      printf(" %12g %12g", primal_infeasibility, dual_infeasibility);
      debugBasicSolutionVariable(
          report, primal_feasibility_tolerance, dual_feasibility_tolerance,
          status, lower, upper, value, dual, num_non_basic_var, num_basic_var,
          off_bound_nonbasic, primal_infeasibility, dual_infeasibility);
      printf("\n");
    }
    dual_activities[iCol] = lp.colCost_[iCol];
    for (int el = lp.Astart_[iCol]; el < lp.Astart_[iCol + 1]; el++) {
      int iRow = lp.Aindex_[el];
      double Avalue = lp.Avalue_[el];
      primal_activities[iRow] += value * Avalue;
      dual_activities[iCol] += solution.row_dual[iRow] * Avalue;
    }
  }
  bool report = options.highs_debug_level > HIGHS_DEBUG_LEVEL_EXPENSIVE;
  header_written = false;
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    double primal_residual_error =
        std::fabs(primal_activities[iRow] - solution.row_value[iRow]);
    if (primal_residual_error > large_residual_error) {
      if (report) {
        if (!header_written) {
          printf(
              "\nRow primal residuals\nIndex     Activity     Solution     "
              "Residual\n");
          header_written = true;
        }
        printf("%5d %12g %12g %12g\n", iRow, primal_activities[iRow],
               solution.row_value[iRow], primal_residual_error);
      }
      num_primal_residual++;
    }
    max_primal_residual = std::max(primal_residual_error, max_primal_residual);
    sum_primal_residual += primal_residual_error;
  }
  header_written = false;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    double dual_residual_error =
        std::fabs(dual_activities[iCol] - solution.col_dual[iCol]);
    if (dual_residual_error > large_residual_error) {
      if (report) {
        if (!header_written) {
          printf(
              "\nRow dual residuals\nIndex     Activity     Solution     "
              "Residual\n");
          header_written = true;
        }
        printf("%5d %12g %12g %12g\n", iCol, dual_activities[iCol],
               solution.col_dual[iCol], dual_residual_error);
      }
      num_dual_residual++;
    }
    max_dual_residual = std::max(dual_residual_error, max_dual_residual);
    sum_dual_residual += dual_residual_error;
  }
  header_written = false;
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    double lower = lp.rowLower_[iRow];
    double upper = lp.rowUpper_[iRow];
    double value = solution.row_value[iRow];
    double dual = -solution.row_dual[iRow];
    HighsBasisStatus status = basis.row_status[iRow];
    if (status != HighsBasisStatus::BASIC) dual_objective_value += value * dual;
    // Flip dual according to lp.sense_
    dual *= (int)lp.sense_;
    bool report = false;
    bool query = debugBasicSolutionVariable(
        report, primal_feasibility_tolerance, dual_feasibility_tolerance,
        status, lower, upper, value, dual, num_non_basic_var, num_basic_var,
        off_bound_nonbasic, primal_infeasibility, dual_infeasibility);
    if (off_bound_nonbasic > 0) num_off_bound_nonbasic++;
    max_off_bound_nonbasic =
        std::max(off_bound_nonbasic, max_off_bound_nonbasic);
    sum_off_bound_nonbasic += off_bound_nonbasic;
    if (primal_infeasibility > primal_feasibility_tolerance)
      num_primal_infeasibilities++;
    max_primal_infeasibility =
        std::max(primal_infeasibility, max_primal_infeasibility);
    sum_primal_infeasibilities += primal_infeasibility;
    if (status == HighsBasisStatus::BASIC) {
      double abs_basic_dual = dual_infeasibility;
      if (abs_basic_dual > 0) {
        num_nonzero_basic_duals++;
        if (abs_basic_dual > dual_feasibility_tolerance)
          num_large_nonzero_basic_duals++;
        max_nonzero_basic_dual =
            std::max(abs_basic_dual, max_nonzero_basic_dual);
        sum_nonzero_basic_duals += abs_basic_dual;
      }
    } else {
      if (dual_infeasibility > dual_feasibility_tolerance)
        num_dual_infeasibilities++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
    report =
        options.highs_debug_level > HIGHS_DEBUG_LEVEL_EXPENSIVE ||
        (options.highs_debug_level == HIGHS_DEBUG_LEVEL_EXPENSIVE && query);
    if (report) {
      if (!header_written) {
        printf(
            "Rows\nIndex NonBs Mv [          LB,           UB]       Primal    "
            "     Dual    PrimalIfs      DualIfs\n");
        header_written = true;
      }
      printf("%5d %5d [%12g, %12g] %12g %12g", iRow, (int)status, lower, upper,
             value, dual);
      printf(" %12g %12g", primal_infeasibility, dual_infeasibility);
      debugBasicSolutionVariable(
          report, primal_feasibility_tolerance, dual_feasibility_tolerance,
          status, lower, upper, value, dual, num_non_basic_var, num_basic_var,
          off_bound_nonbasic, primal_infeasibility, dual_infeasibility);
      printf("\n");
    }
  }
}

bool debugBasicSolutionVariable(
    bool report, const double primal_feasibility_tolerance,
    const double dual_feasibility_tolerance, const HighsBasisStatus status,
    const double lower, const double upper, const double value,
    const double dual, int& num_non_basic_var, int& num_basic_var,
    double& off_bound_nonbasic, double& primal_infeasibility,
    double& dual_infeasibility) {
  double middle = (lower + upper) * 0.5;

  bool query = false;
  bool count = !report;
  off_bound_nonbasic = 0;
  double primal_residual = std::max(lower - value, value - upper);
  primal_infeasibility = std::max(primal_residual, 0.);
  // ToDo Strange: nonbasic_flag seems to be inverted???
  if (status == HighsBasisStatus::BASIC) {
    // Basic variable: look for primal infeasibility
    if (count) num_basic_var++;
    if (primal_infeasibility > primal_feasibility_tolerance) {
      // Outside a bound
      if (value < lower) {
        query = true;
        if (report)
          printf(": Basic below lower bound by %12g", primal_residual);
      } else {
        query = true;
        if (report)
          printf(": Basic above upper bound by %12g", primal_residual);
      }
    }
    dual_infeasibility = std::fabs(dual);
    if (dual_infeasibility > dual_feasibility_tolerance) {
      query = true;
      if (report) printf(": Dual infeasibility of %12g", dual_infeasibility);
    }
  } else {
    // Nonbasic variable: look for primal and dual infeasibility
    if (count) num_non_basic_var++;

    if (primal_infeasibility > primal_feasibility_tolerance) {
      // Outside a bound
      off_bound_nonbasic = primal_infeasibility;
      dual_infeasibility = 0;
      if (value < lower) {
        query = true;
        if (report)
          printf(": Nonbasic below lower bound by %12g", primal_residual);
      } else {
        query = true;
        if (report)
          printf(": Nonbasic above upper bound by %12g", primal_residual);
      }
    } else if (primal_residual >= -primal_feasibility_tolerance) {
      // At a bound: check for dual feasibility
      off_bound_nonbasic = std::fabs(primal_residual);
      if (lower < upper) {
        // Non-fixed variable
        if (value < middle) {
          // At lower
          dual_infeasibility = std::max(-dual, 0.);
          if (dual_infeasibility > dual_feasibility_tolerance) {
            // Dual infeasiblility
            query = true;
            if (report)
              printf(": Dual infeasibility of %12g", dual_infeasibility);
          }
        } else {
          // At Upper
          dual_infeasibility = std::max(dual, 0.);
          if (dual_infeasibility > dual_feasibility_tolerance) {
            // Dual infeasiblility
            query = true;
            if (report)
              printf(": Dual infeasibility of %12g", dual_infeasibility);
          }
        }
      } else {
        // Fixed variable
        dual_infeasibility = 0;
      }
    } else {
      // Between bounds (or free)
      if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
        // Free
        if (report) printf(": Nonbasic free");
      } else {
        query = true;
        if (report) printf(": Nonbasic off bound by %12g", -primal_residual);
        off_bound_nonbasic = -primal_residual;
      }
      dual_infeasibility = std::fabs(dual);
      if (dual_infeasibility > dual_feasibility_tolerance) {
        query = true;
        if (report) printf(": Dual infeasibility of %12g", dual_infeasibility);
      }
    }
  }
  query = false;
  return query;
}

HighsDebugStatus debugAnalysePrimalDualErrors(
    const HighsOptions& options, HighsPrimalDualErrors& primal_dual_errors) {
  std::string value_adjective;
  int report_level;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  const bool force_report =
      options.highs_debug_level >= HIGHS_DEBUG_LEVEL_COSTLY;
  if (primal_dual_errors.num_nonzero_basic_duals) {
    value_adjective = "Error";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  } else {
    value_adjective = "";
    report_level = ML_NONE;
    return_status = HighsDebugStatus::OK;
  }
  if (force_report) report_level = ML_ALWAYS;
  HighsPrintMessage(options.output, options.message_level, report_level,
                    "PrDuErrors : %-9s Nonzero basic duals:       num = %2d; "
                    "max = %9.4g; sum = %9.4g\n",
                    value_adjective.c_str(),
                    primal_dual_errors.num_nonzero_basic_duals,
                    primal_dual_errors.max_nonzero_basic_dual,
                    primal_dual_errors.sum_nonzero_basic_duals);

  if (primal_dual_errors.num_off_bound_nonbasic) {
    value_adjective = "Error";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::LOGICAL_ERROR;
  } else {
    value_adjective = "";
    report_level = ML_NONE;
    return_status = HighsDebugStatus::OK;
  }
  if (force_report) report_level = ML_ALWAYS;
  HighsPrintMessage(options.output, options.message_level, report_level,
                    "PrDuErrors : %-9s Off-bound nonbasic values: num = %2d; "
                    "max = %9.4g; sum = %9.4g\n",
                    value_adjective.c_str(),
                    primal_dual_errors.num_off_bound_nonbasic,
                    primal_dual_errors.max_off_bound_nonbasic,
                    primal_dual_errors.sum_off_bound_nonbasic);

  if (primal_dual_errors.max_primal_residual > excessive_residual_error) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::ERROR;
  } else if (primal_dual_errors.max_primal_residual > large_residual_error) {
    value_adjective = "Large";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::WARNING;
  } else {
    value_adjective = "";
    report_level = ML_VERBOSE;
    return_status = HighsDebugStatus::OK;
  }
  if (force_report) report_level = ML_ALWAYS;
  HighsPrintMessage(options.output, options.message_level, report_level,
                    "PrDuErrors : %-9s Primal residual:           num = %2d; "
                    "max = %9.4g; sum = %9.4g\n",
                    value_adjective.c_str(),
                    primal_dual_errors.num_primal_residual,
                    primal_dual_errors.max_primal_residual,
                    primal_dual_errors.sum_primal_residual);

  if (primal_dual_errors.max_dual_residual > excessive_residual_error) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::ERROR;
  } else if (primal_dual_errors.max_dual_residual > large_residual_error) {
    value_adjective = "Large";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::WARNING;
  } else {
    value_adjective = "";
    report_level = ML_VERBOSE;
    return_status = HighsDebugStatus::OK;
  }
  if (force_report) report_level = ML_ALWAYS;
  HighsPrintMessage(options.output, options.message_level, report_level,
                    "PrDuErrors : %-9s Dual residual:             num = %2d; "
                    "max = %9.4g; sum = %9.4g\n",
                    value_adjective.c_str(),
                    primal_dual_errors.num_dual_residual,
                    primal_dual_errors.max_dual_residual,
                    primal_dual_errors.sum_dual_residual);

  return return_status;
}

HighsDebugStatus debugCompareSolutionParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1) {
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  return_status =
      debugWorseStatus(debugCompareSolutionObjectiveParams(
                           options, solution_params0, solution_params1),
                       return_status);
  return_status =
      debugWorseStatus(debugCompareSolutionStatusParams(
                           options, solution_params0, solution_params1),
                       return_status);
  return_status =
      debugWorseStatus(debugCompareSolutionInfeasibilityParams(
                           options, solution_params0, solution_params1),
                       return_status);
  return return_status;
}

HighsDebugStatus debugCompareSolutionObjectiveParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1) {
  return debugCompareSolutionParamValue(
      "objective_function_value", options,
      solution_params0.objective_function_value,
      solution_params1.objective_function_value);
}

HighsDebugStatus debugCompareSolutionStatusParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1) {
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  return_status = debugWorseStatus(
      debugCompareSolutionParamInteger("primal_status", options,
                                       solution_params0.primal_status,
                                       solution_params1.primal_status),
      return_status);
  return_status =
      debugWorseStatus(debugCompareSolutionParamInteger(
                           "dual_status", options, solution_params0.dual_status,
                           solution_params1.dual_status),
                       return_status);
  return return_status;
}

HighsDebugStatus debugCompareSolutionInfeasibilityParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1) {
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  return_status =
      debugWorseStatus(debugCompareSolutionParamInteger(
                           "num_primal_infeasibilities", options,
                           solution_params0.num_primal_infeasibilities,
                           solution_params1.num_primal_infeasibilities),
                       return_status);
  return_status =
      debugWorseStatus(debugCompareSolutionParamValue(
                           "sum_primal_infeasibilities", options,
                           solution_params0.sum_primal_infeasibilities,
                           solution_params1.sum_primal_infeasibilities),
                       return_status);
  return_status = debugWorseStatus(
      debugCompareSolutionParamValue("max_primal_infeasibility", options,
                                     solution_params0.max_primal_infeasibility,
                                     solution_params1.max_primal_infeasibility),
      return_status);

  return_status =
      debugWorseStatus(debugCompareSolutionParamInteger(
                           "num_dual_infeasibilities", options,
                           solution_params0.num_dual_infeasibilities,
                           solution_params1.num_dual_infeasibilities),
                       return_status);
  return_status = debugWorseStatus(
      debugCompareSolutionParamValue("sum_dual_infeasibilities", options,
                                     solution_params0.sum_dual_infeasibilities,
                                     solution_params1.sum_dual_infeasibilities),
      return_status);
  return_status = debugWorseStatus(
      debugCompareSolutionParamValue("max_dual_infeasibility", options,
                                     solution_params0.max_dual_infeasibility,
                                     solution_params1.max_dual_infeasibility),
      return_status);
  return return_status;
}

HighsDebugStatus debugCompareSolutionParamValue(const string name,
                                                const HighsOptions& options,
                                                const double v0,
                                                const double v1) {
  if (v0 == v1) return HighsDebugStatus::OK;
  double delta = highsRelativeDifference(v0, v1);
  std::string value_adjective;
  int report_level;
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  if (delta > excessive_relative_solution_param_error) {
    value_adjective = "Excessive";
    report_level = ML_ALWAYS;
    return_status = HighsDebugStatus::ERROR;
  } else if (delta > large_relative_solution_param_error) {
    value_adjective = "Large";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::WARNING;
  } else {
    value_adjective = "OK";
    report_level = ML_VERBOSE;
  }
  HighsPrintMessage(options.output, options.message_level, report_level,
                    "SolutionPar:  %-9s relative difference of %9.4g for %s\n",
                    value_adjective.c_str(), delta, name.c_str());
  return return_status;
}

HighsDebugStatus debugCompareSolutionParamInteger(const string name,
                                                  const HighsOptions& options,
                                                  const int v0, const int v1) {
  if (v0 == v1) return HighsDebugStatus::OK;
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "SolutionPar:  difference of %d for %s\n", v1 - v0,
                    name.c_str());
  return HighsDebugStatus::LOGICAL_ERROR;
}

void debugReportHighsBasicSolution(const string message,
                                   const HighsOptions& options,
                                   const HighsSolutionParams& solution_params,
                                   const HighsModelStatus model_status) {
  HighsPrintMessage(options.output, options.message_level, ML_ALWAYS,
                    "\nHiGHS basic solution: %s\n", message.c_str());
  HighsPrintMessage(
      options.output, options.message_level, ML_ALWAYS,
      "Infeas:                Pr %d(Max %.4g, Sum %.4g); Du %d(Max %.4g, "
      "Sum %.4g); Status: %s\n",
      solution_params.num_primal_infeasibilities,
      solution_params.max_primal_infeasibility,
      solution_params.sum_primal_infeasibilities,
      solution_params.num_dual_infeasibilities,
      solution_params.max_dual_infeasibility,
      solution_params.sum_dual_infeasibilities,
      utilHighsModelStatusToString(model_status).c_str());
}
