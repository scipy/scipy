/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolutionDebug.cpp
 * @brief
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

// Called from HighsSolve - solveLp
HighsDebugStatus debugHighsLpSolution(
    const std::string& message, const HighsLpSolverObject& solver_object) {
  // Non-trivially expensive analysis of a solution to a model
  //
  // Called to check the unscaled model status and solution params
  const bool check_model_status_and_highs_info = true;
  // Define an empty Hessian
  HighsHessian hessian;
  return debugHighsSolution(message, solver_object.options_, solver_object.lp_,
                            hessian, solver_object.solution_,
                            solver_object.basis_, solver_object.model_status_,
                            solver_object.highs_info_,
                            check_model_status_and_highs_info);
}

HighsDebugStatus debugHighsSolution(const string& message,
                                    const HighsOptions& options,
                                    const HighsModel& model,
                                    const HighsSolution& solution,
                                    const HighsBasis& basis) {
  // Non-trivially expensive analysis of a solution to a model
  //
  // Called to report on KKT errors after solving a model when only
  // the solution (possibly only primal) and (possibly) basis are
  // known
  //
  // Set up a HighsModelStatus and HighsInfo just to
  // complete the parameter list.By setting
  // check_model_status_and_highs_info to be false they aren't
  // used.
  HighsModelStatus dummy_model_status;
  HighsInfo dummy_highs_info;
  // Call resetModelStatusAndSolutionParams to side-step compiler
  // warning.
  resetModelStatusAndHighsInfo(dummy_model_status, dummy_highs_info);
  const bool check_model_status_and_highs_info = false;
  return debugHighsSolution(
      message, options, model.lp_, model.hessian_, solution, basis,
      dummy_model_status, dummy_highs_info, check_model_status_and_highs_info);
}

HighsDebugStatus debugHighsSolution(
    const string& message, const HighsOptions& options, const HighsModel& model,
    const HighsSolution& solution, const HighsBasis& basis,
    const HighsModelStatus model_status, const HighsInfo& info) {
  // Non-trivially expensive analysis of a solution to a model
  //
  // Called to check the HiGHS model_status and info
  //
  // Refer to info as highs_info so general method can be used
  //
  const HighsInfo& highs_info = info;
  const bool check_model_status_and_highs_info = true;
  return debugHighsSolution(message, options, model.lp_, model.hessian_,
                            solution, basis, model_status, highs_info,
                            check_model_status_and_highs_info);
}

HighsDebugStatus debugHighsSolution(
    const std::string& message, const HighsOptions& options, const HighsLp& lp,
    const HighsHessian& hessian, const HighsSolution& solution,
    const HighsBasis& basis, const HighsModelStatus model_status,
    const HighsInfo& highs_info, const bool check_model_status_and_highs_info) {
  // Non-trivially expensive analysis of a solution to a model
  //
  // Called to possibly check the model_status and highs_info,
  // and then report on the KKT and model status, plus any errors
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status;
  // Use local_model_status to for checking - or if it's not known
  HighsModelStatus local_model_status = HighsModelStatus::kNotset;
  // Use local_highs_info to determine highs_info for
  // checking - or if it's not known
  HighsInfo local_highs_info;
  if (check_model_status_and_highs_info) {
    double local_objective_function_value = 0;
    if (solution.value_valid)
      local_objective_function_value =
          lp.objectiveValue(solution.col_value) +
          hessian.objectiveValue(solution.col_value);
    local_highs_info.objective_function_value = local_objective_function_value;
  }
  HighsPrimalDualErrors primal_dual_errors;
  // Determine the extent to which KKT conditions are not satisfied,
  // accumulating data on primal/dual errors relating to any basis
  // implications and excessive residuals
  const bool get_residuals =
      true;  // options.highs_debug_level >= kHighsDebugLevelCostly;

  vector<double> gradient;
  const bool is_qp = hessian.dim_ > 0;
  if (is_qp) {
    hessian.product(solution.col_value, gradient);
  } else {
    gradient.assign(lp.num_col_, 0);
  }
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    gradient[iCol] += lp.col_cost_[iCol];
  getKktFailures(options, is_qp, lp, gradient, solution, local_highs_info,
                 get_residuals);
  getPrimalDualBasisErrors(options, lp, solution, basis, primal_dual_errors);
  getPrimalDualGlpsolErrors(options, lp, gradient, solution,
                            primal_dual_errors);
  HighsInt& num_primal_infeasibility =
      local_highs_info.num_primal_infeasibilities;
  HighsInt& num_dual_infeasibility = local_highs_info.num_dual_infeasibilities;
  if (check_model_status_and_highs_info) {
    // Can assume that model_status and highs_info are known, so should be
    // checked
    local_model_status = model_status;
    // Check that highs_info is the same as when computed from scratch
    return_status =
        debugCompareHighsInfo(options, highs_info, local_highs_info);
    if (return_status != HighsDebugStatus::kOk) return return_status;
  } else {
    // Determine whether optimality can be reported
    if (num_primal_infeasibility == 0 && num_dual_infeasibility == 0)
      local_model_status = HighsModelStatus::kOptimal;
  }
  if (check_model_status_and_highs_info &&
      model_status == HighsModelStatus::kOptimal) {
    bool error_found = false;
    if (num_primal_infeasibility > 0) {
      error_found = true;
      highsLogDev(options.log_options, HighsLogType::kError,
                  "debugHighsLpSolution: %" HIGHSINT_FORMAT
                  " primal infeasibilities but model status is %s\n",
                  num_primal_infeasibility,
                  utilModelStatusToString(model_status).c_str());
    }
    if (num_dual_infeasibility > 0) {
      error_found = true;
      highsLogDev(options.log_options, HighsLogType::kError,
                  "debugHighsLpSolution: %" HIGHSINT_FORMAT
                  " dual infeasibilities but model status is %s\n",
                  num_dual_infeasibility,
                  utilModelStatusToString(model_status).c_str());
    }
    if (error_found) return HighsDebugStatus::kLogicalError;
  }
  // Report on the solution
  debugReportHighsSolution(message, options.log_options, local_highs_info,
                           local_model_status);
  // Analyse the primal and dual errors
  return_status = debugAnalysePrimalDualErrors(options, primal_dual_errors);
  return return_status;
}

void debugReportHighsSolution(const string& message,
                              const HighsLogOptions& log_options,
                              const HighsInfo& highs_info,
                              const HighsModelStatus model_status) {
  highsLogDev(log_options, HighsLogType::kInfo, "\nHiGHS solution: %s\n",
              message.c_str());
  if (highs_info.num_primal_infeasibilities >= 0 ||
      highs_info.num_dual_infeasibilities >= 0) {
    highsLogDev(log_options, HighsLogType::kInfo, "Infeas:                ");
  }
  if (highs_info.num_primal_infeasibilities >= 0) {
    highsLogDev(log_options, HighsLogType::kInfo,
                "Pr %" HIGHSINT_FORMAT "(Max %.4g, Sum %.4g); ",
                highs_info.num_primal_infeasibilities,
                highs_info.max_primal_infeasibility,
                highs_info.sum_primal_infeasibilities);
  }
  if (highs_info.num_dual_infeasibilities >= 0) {
    highsLogDev(log_options, HighsLogType::kInfo,
                "Du %" HIGHSINT_FORMAT "(Max %.4g, Sum %.4g); ",
                highs_info.num_dual_infeasibilities,
                highs_info.max_dual_infeasibility,
                highs_info.sum_dual_infeasibilities);
  }
  highsLogDev(log_options, HighsLogType::kInfo, "Status: %s\n",
              utilModelStatusToString(model_status).c_str());
}

HighsDebugStatus debugHighsBasisConsistent(const HighsOptions& options,
                                           const HighsLp& lp,
                                           const HighsBasis& basis) {
  // Cheap analysis of a HiGHS basis, checking vector sizes, numbers
  // of basic/nonbasic variables
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  if (!basis.valid) return return_status;
  bool consistent = isBasisConsistent(lp, basis);
  if (!consistent) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "HiGHS basis inconsistency\n");
    assert(consistent);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus debugBasisRightSize(const HighsOptions& options,
                                     const HighsLp& lp,
                                     const HighsBasis& basis) {
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  bool right_size = isBasisRightSize(lp, basis);
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "HiGHS basis size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus debugPrimalSolutionRightSize(const HighsOptions& options,
                                              const HighsLp& lp,
                                              const HighsSolution& solution) {
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  bool right_size = isPrimalSolutionRightSize(lp, solution);
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "HiGHS primal solution size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus debugDualSolutionRightSize(const HighsOptions& options,
                                            const HighsLp& lp,
                                            const HighsSolution& solution) {
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  bool right_size = isDualSolutionRightSize(lp, solution);
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "HiGHS dual solution size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

// Methods below are not called externally

HighsDebugStatus debugAnalysePrimalDualErrors(
    const HighsOptions& options, HighsPrimalDualErrors& primal_dual_errors) {
  std::string value_adjective;
  HighsLogType report_level;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const bool force_report = options.highs_debug_level >= kHighsDebugLevelCostly;
  if (primal_dual_errors.num_nonzero_basic_duals >= 0) {
    if (primal_dual_errors.num_nonzero_basic_duals > 0) {
      value_adjective = "Error";
      report_level = HighsLogType::kError;
      return_status = HighsDebugStatus::kLogicalError;
    } else {
      value_adjective = "";
      report_level = HighsLogType::kVerbose;
      return_status = HighsDebugStatus::kOk;
    }
    if (force_report) report_level = HighsLogType::kInfo;
    highsLogDev(
        options.log_options, report_level,
        "PrDuErrors : %-9s Nonzero basic duals:       num = %7" HIGHSINT_FORMAT
        "; "
        "max = %9.4g; sum = %9.4g\n",
        value_adjective.c_str(), primal_dual_errors.num_nonzero_basic_duals,
        primal_dual_errors.max_nonzero_basic_dual,
        primal_dual_errors.sum_nonzero_basic_duals);
  }
  if (primal_dual_errors.num_off_bound_nonbasic >= 0) {
    if (primal_dual_errors.num_off_bound_nonbasic > 0) {
      value_adjective = "Error";
      report_level = HighsLogType::kError;
      return_status = HighsDebugStatus::kLogicalError;
    } else {
      value_adjective = "";
      report_level = HighsLogType::kVerbose;
      return_status = HighsDebugStatus::kOk;
    }
    if (force_report) report_level = HighsLogType::kInfo;
    highsLogDev(
        options.log_options, report_level,
        "PrDuErrors : %-9s Off-bound nonbasic values: num = %7" HIGHSINT_FORMAT
        "; "
        "max = %9.4g; sum = %9.4g\n",
        value_adjective.c_str(), primal_dual_errors.num_off_bound_nonbasic,
        primal_dual_errors.max_off_bound_nonbasic,
        primal_dual_errors.sum_off_bound_nonbasic);
  }
  if (primal_dual_errors.glpsol_num_primal_residual_errors >= 0) {
    if (primal_dual_errors.glpsol_max_primal_residual.absolute_value >
        excessive_residual_error) {
      value_adjective = "Excessive";
      report_level = HighsLogType::kError;
      return_status = HighsDebugStatus::kError;
    } else if (primal_dual_errors.glpsol_max_primal_residual.absolute_value >
               large_residual_error) {
      value_adjective = "Large";
      report_level = HighsLogType::kDetailed;
      return_status = HighsDebugStatus::kWarning;
    } else {
      value_adjective = "";
      report_level = HighsLogType::kVerbose;
      return_status = HighsDebugStatus::kOk;
    }
    if (force_report) report_level = HighsLogType::kInfo;
    highsLogDev(
        options.log_options, report_level,
        "PrDuErrors : %-9s Primal residual:           num = %7" HIGHSINT_FORMAT
        "; "
        "max = %9.4g\n",
        value_adjective.c_str(),
        primal_dual_errors.glpsol_num_primal_residual_errors,
        primal_dual_errors.glpsol_max_primal_residual.absolute_value);
  }
  if (primal_dual_errors.glpsol_num_dual_residual_errors >= 0) {
    if (primal_dual_errors.glpsol_max_dual_residual.absolute_value >
        excessive_residual_error) {
      value_adjective = "Excessive";
      report_level = HighsLogType::kError;
      return_status = HighsDebugStatus::kError;
    } else if (primal_dual_errors.glpsol_max_dual_residual.absolute_value >
               large_residual_error) {
      value_adjective = "Large";
      report_level = HighsLogType::kDetailed;
      return_status = HighsDebugStatus::kWarning;
    } else {
      value_adjective = "";
      report_level = HighsLogType::kVerbose;
      return_status = HighsDebugStatus::kOk;
    }
    if (force_report) report_level = HighsLogType::kInfo;
    highsLogDev(
        options.log_options, report_level,
        "PrDuErrors : %-9s Dual residual:             num = %7" HIGHSINT_FORMAT
        "; "
        "max = %9.4g\n",
        value_adjective.c_str(),
        primal_dual_errors.glpsol_num_dual_residual_errors,
        primal_dual_errors.glpsol_max_dual_residual.absolute_value);
  }
  return return_status;
}

HighsDebugStatus debugCompareHighsInfo(const HighsOptions& options,
                                       const HighsInfo& highs_info0,
                                       const HighsInfo& highs_info1) {
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  return_status = debugWorseStatus(
      debugCompareHighsInfoObjective(options, highs_info0, highs_info1),
      return_status);
  return_status = debugWorseStatus(
      debugCompareHighsInfoStatus(options, highs_info0, highs_info1),
      return_status);
  return_status = debugWorseStatus(
      debugCompareHighsInfoInfeasibility(options, highs_info0, highs_info1),
      return_status);
  return return_status;
}

HighsDebugStatus debugCompareHighsInfoObjective(const HighsOptions& options,
                                                const HighsInfo& highs_info0,
                                                const HighsInfo& highs_info1) {
  return debugCompareHighsInfoDouble("objective_function_value", options,
                                     highs_info0.objective_function_value,
                                     highs_info1.objective_function_value);
}

HighsDebugStatus debugCompareHighsInfoStatus(const HighsOptions& options,
                                             const HighsInfo& highs_info0,
                                             const HighsInfo& highs_info1) {
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  return_status = debugWorseStatus(
      debugCompareHighsInfoInteger("primal_status", options,
                                   highs_info0.primal_solution_status,
                                   highs_info1.primal_solution_status),
      return_status);
  return_status = debugWorseStatus(
      debugCompareHighsInfoInteger("dual_status", options,
                                   highs_info0.dual_solution_status,
                                   highs_info1.dual_solution_status),
      return_status);
  return return_status;
}

HighsDebugStatus debugCompareHighsInfoInfeasibility(
    const HighsOptions& options, const HighsInfo& highs_info0,
    const HighsInfo& highs_info1) {
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  return_status = debugWorseStatus(
      debugCompareHighsInfoInteger("num_primal_infeasibility", options,
                                   highs_info0.num_primal_infeasibilities,
                                   highs_info1.num_primal_infeasibilities),
      return_status);
  return_status = debugWorseStatus(
      debugCompareHighsInfoDouble("sum_primal_infeasibility", options,
                                  highs_info0.sum_primal_infeasibilities,
                                  highs_info1.sum_primal_infeasibilities),
      return_status);
  return_status = debugWorseStatus(
      debugCompareHighsInfoDouble("max_primal_infeasibility", options,
                                  highs_info0.max_primal_infeasibility,
                                  highs_info1.max_primal_infeasibility),
      return_status);

  return_status = debugWorseStatus(
      debugCompareHighsInfoInteger("num_dual_infeasibility", options,
                                   highs_info0.num_dual_infeasibilities,
                                   highs_info1.num_dual_infeasibilities),
      return_status);
  return_status = debugWorseStatus(
      debugCompareHighsInfoDouble("sum_dual_infeasibility", options,
                                  highs_info0.sum_dual_infeasibilities,
                                  highs_info1.sum_dual_infeasibilities),
      return_status);
  return_status = debugWorseStatus(
      debugCompareHighsInfoDouble("max_dual_infeasibility", options,
                                  highs_info0.max_dual_infeasibility,
                                  highs_info1.max_dual_infeasibility),
      return_status);
  return return_status;
}

HighsDebugStatus debugCompareHighsInfoDouble(const string& name,
                                             const HighsOptions& options,
                                             const double v0, const double v1) {
  if (v0 == v1) return HighsDebugStatus::kOk;
  double delta = highsRelativeDifference(v0, v1);
  std::string value_adjective;
  HighsLogType report_level;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  if (delta > excessive_relative_solution_param_error) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kError;
  } else if (delta > large_relative_solution_param_error) {
    value_adjective = "Large";
    report_level = HighsLogType::kDetailed;
    return_status = HighsDebugStatus::kWarning;
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
  }
  highsLogDev(options.log_options, report_level,
              "SolutionPar:  %-9s relative difference of %9.4g for %s\n",
              value_adjective.c_str(), delta, name.c_str());
  return return_status;
}

HighsDebugStatus debugCompareHighsInfoInteger(const string& name,
                                              const HighsOptions& options,
                                              const HighsInt v0,
                                              const HighsInt v1) {
  if (v0 == v1) return HighsDebugStatus::kOk;
  highsLogDev(options.log_options, HighsLogType::kError,
              "SolutionPar:  difference of %" HIGHSINT_FORMAT " for %s\n",
              v1 - v0, name.c_str());
  return HighsDebugStatus::kLogicalError;
}
