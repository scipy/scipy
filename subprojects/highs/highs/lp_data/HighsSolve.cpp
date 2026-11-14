/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolve.cpp
 * @brief Class-independent utilities for HiGHS
 */

#include "HighsExternalApi.h"
#include "ipm/IpxWrapper.h"
#include "lp_data/HighsSolutionDebug.h"
#include "mip/HighsMipSolver.h"
#include "model/HighsHessianUtils.h"
#include "pdlp/CupdlpWrapper.h"
#include "pdlp/HiPdlpWrapper.h"
#include "qpsolver/QpAsmWrapper.h"
#include "simplex/HApp.h"

#define SOLVE_CATCH_CALL(Solve, SolveString)                                  \
  do {                                                                        \
    try {                                                                     \
      call_status = Solve;                                                    \
    } catch (const std::exception& exception) {                               \
      solver_object.model_status_ =                                           \
          handleExceptionIsOom(options.log_options, SolveString, exception)   \
              ? HighsModelStatus::kMemoryLimit                                \
              : HighsModelStatus::kSolveError;                                \
      call_status = HighsStatus::kError;                                      \
    } catch (const HighsTask::Interrupt&) {                                   \
      highsLogDev(options.log_options, HighsLogType::kError,                  \
                  "HighsTask interrupt when solving with %s\n", SolveString); \
      solver_object.model_status_ = HighsModelStatus::kSolveError;            \
      call_status = HighsStatus::kError;                                      \
    }                                                                         \
  } while (0)

// The method below runs the simplex, IPX, HiPO or PDLP solver on the LP
HighsStatus solveLp(HighsLpSolverObject& solver_object,
                    const std::string& message) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsOptions& options = solver_object.options_;
  HighsProfiling* profiling = solver_object.profiling_;
  // Reset unscaled model status and solution params - except for
  // iteration counts
  resetModelStatusAndHighsInfo(solver_object);
  highsLogUser(options.log_options, HighsLogType::kInfo,
               (message + "\n").c_str());
  if (options.highs_debug_level > kHighsDebugLevelMin) {
    // Shouldn't have to check validity of the LP since this is done when it is
    // loaded or modified
    call_status = assessLp(solver_object.lp_, options);
    // If any errors have been found or normalisation carried out,
    // call_status will be ERROR or WARNING, so only valid return is OK.
    assert(call_status == HighsStatus::kOk);
    return_status = interpretCallStatus(options.log_options, call_status,
                                        return_status, "assessLp");
    if (return_status == HighsStatus::kError) return return_status;
  }
  const bool use_only_ipm = useIpm(options.solver) || options.run_centring;
  bool use_hipo = useHipo(options, kSolverString, solver_object.lp_);

  const bool use_ipx = use_only_ipm && !use_hipo;
  // Now actually solve LPs!
  //
  // lambda for solving LP by simplex
  auto simplexSolve = [&]() -> HighsStatus {
    return_status = HighsStatus::kOk;
    SOLVE_CATCH_CALL(solveLpSimplex(solver_object), "Simplex");
    return_status = interpretCallStatus(options.log_options, call_status,
                                        return_status, "solveLpSimplex");
    if (return_status == HighsStatus::kError) return return_status;
    if (!isSolutionRightSize(solver_object.lp_, solver_object.solution_)) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Inconsistent solution returned from solver\n");
      return_status = HighsStatus::kError;
    }
    return return_status;
  };
  if (!solver_object.lp_.num_row_ || solver_object.lp_.numNz() == 0) {
    // LP is unconstrained due to having no rows or a zero constraint
    // matrix, so solve directly
    call_status = solveUnconstrainedLp(solver_object);
    return_status = interpretCallStatus(options.log_options, call_status,
                                        return_status, "solveUnconstrainedLp");
    if (return_status == HighsStatus::kError) return return_status;
  } else if (use_only_ipm || usePdlp(options.solver)) {
    // Use IPM or PDLP
    if (use_only_ipm) {
      // Use IPM to solve the LP
      if (use_hipo) {
        // Use HiPO to solve the LP
        SOLVE_CATCH_CALL(solveLpHipo(solver_object), "HiPO(LP)");
        return_status = interpretCallStatus(options.log_options, call_status,
                                            return_status, "solveLpHipo");
      } else if (use_ipx) {
        // Use IPX to solve the LP
        SOLVE_CATCH_CALL(solveLpIpx(solver_object), "IPX");
        return_status = interpretCallStatus(options.log_options, call_status,
                                            return_status, "solveLpIpx");
      }
    } else {
      // Use cuPDLP-C or HiPDLP to solve the LP
      profiling->start(kSubSolverPdlp);
      if (options.solver == kPdlpString) {
        SOLVE_CATCH_CALL(solveLpCupdlp(solver_object), "cuPDLP-C");
      } else {
        SOLVE_CATCH_CALL(solveLpHiPdlp(solver_object), "HiPDLP");
      }
      profiling->stop(kSubSolverPdlp);
      return_status = interpretCallStatus(options.log_options, call_status,
                                          return_status, "solveLp-Pdlp");
    }
    // Check for error return
    if (return_status == HighsStatus::kError) return return_status;

    // Non-error return requires a primal solution
    assert(solver_object.solution_.value_valid);

    if (useIpm(options.solver) || options.run_centring) {
      // Setting the IPM-specific values of (highs_)info_ has been done in
      // solveLpHipo/Ipx
      const bool unwelcome_ipx_status =
          solver_object.model_status_ == HighsModelStatus::kUnknown ||
          (solver_object.model_status_ ==
               HighsModelStatus::kUnboundedOrInfeasible &&
           !options.allow_unbounded_or_infeasible);
      if (unwelcome_ipx_status) {
        // When performing an analytic centre calculation, the setting
        // of options.run_crossover is ignored, so simplex clean-up is
        // not possible - or desirable, anyway!
        highsLogUser(
            options.log_options, HighsLogType::kWarning,
            "Unwelcome IPX status of %s: basis is %svalid; solution is "
            "%svalid; run_crossover is \"%s\"\n",
            utilModelStatusToString(solver_object.model_status_).c_str(),
            solver_object.basis_.valid ? "" : "not ",
            solver_object.solution_.value_valid ? "" : "not ",
            options.run_centring ? kHighsOffString.c_str()
                                 : options.run_crossover.c_str());
        const bool allow_simplex_cleanup =
            options.run_crossover != kHighsOffString && !options.run_centring;
        if (allow_simplex_cleanup) {
          // IPX has returned a model status that HiGHS would rather
          // avoid, so perform simplex clean-up if crossover was allowed.
          //
          // This is an unusual situation, and the cost will usually be
          // acceptable. Worst case is if crossover wasn't run, in which
          // case there's no basis to start simplex
          //
          // ToDo: Check whether simplex can exploit the primal solution
          // returned by HiPO/IPX
          highsLogUser(options.log_options, HighsLogType::kWarning,
                       "IPM solution is imprecise, so clean up with simplex\n");
          return_status = simplexSolve();
          if (return_status == HighsStatus::kError) return return_status;
        }  // options.run_crossover == kHighsOnString
           // clang-format off
      }  // unwelcome_ipx_status
      // clang-format on
    }
  } else {
    // Use Simplex
    return_status = simplexSolve();
    if (return_status == HighsStatus::kError) return return_status;
  }
  // Analyse the HiGHS (basic) solution
  if (debugHighsLpSolution(message, solver_object) ==
      HighsDebugStatus::kLogicalError)
    return_status = HighsStatus::kError;
  return return_status;
}

// Solves an unconstrained LP without scaling, setting HighsBasis, HighsSolution
// and HighsInfo
HighsStatus solveUnconstrainedLp(HighsLpSolverObject& solver_object) {
  return (solveUnconstrainedLp(solver_object.options_, solver_object.lp_,
                               solver_object.model_status_,
                               solver_object.highs_info_,
                               solver_object.solution_, solver_object.basis_));
}

// Solves an unconstrained LP without scaling, setting HighsBasis, HighsSolution
// and HighsInfo
HighsStatus solveUnconstrainedLp(const HighsOptions& options, const HighsLp& lp,
                                 HighsModelStatus& model_status,
                                 HighsInfo& highs_info, HighsSolution& solution,
                                 HighsBasis& basis) {
  // Aliase to model status and solution parameters
  resetModelStatusAndHighsInfo(model_status, highs_info);

  // Check that the LP really is unconstrained!
  assert(lp.num_row_ == 0 || lp.numNz() == 0);
  if (lp.num_row_ > 0) {
    // LP has rows, but should only be here if the constraint matrix
    // is zero
    if (lp.numNz() > 0) return HighsStatus::kError;
  }

  highsLogUser(options.log_options, HighsLogType::kInfo,
               "Solving an unconstrained LP with %" HIGHSINT_FORMAT
               " columns\n",
               lp.num_col_);
  solution.col_value.assign(lp.num_col_, 0);
  solution.col_dual.assign(lp.num_col_, 0);
  basis.col_status.assign(lp.num_col_, HighsBasisStatus::kNonbasic);
  // No rows for primal solution, dual solution or basis
  solution.row_value.clear();
  solution.row_dual.clear();
  basis.row_status.clear();

  double primal_feasibility_tolerance = options.primal_feasibility_tolerance;
  double dual_feasibility_tolerance = options.dual_feasibility_tolerance;

  // Initialise the objective value calculation. Done using
  // HighsSolution so offset is vanilla
  double objective = lp.offset_;

  highs_info.num_primal_infeasibilities = 0;
  highs_info.max_primal_infeasibility = 0;
  highs_info.sum_primal_infeasibilities = 0;
  highs_info.num_dual_infeasibilities = 0;
  highs_info.max_dual_infeasibility = 0;
  highs_info.sum_dual_infeasibilities = 0;

  if (lp.num_row_ > 0) {
    // Assign primal, dual and basis status for rows, checking for
    // infeasibility
    for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
      double primal_infeasibility = 0;
      double lower = lp.row_lower_[iRow];
      double upper = lp.row_upper_[iRow];
      if (lower > primal_feasibility_tolerance) {
        // Lower bound too large for zero activity
        primal_infeasibility = lower;
      } else if (upper < -primal_feasibility_tolerance) {
        // Upper bound too small for zero activity
        primal_infeasibility = -upper;
      }
      solution.row_value.push_back(0);
      solution.row_dual.push_back(0);
      basis.row_status.push_back(HighsBasisStatus::kBasic);
      if (primal_infeasibility > primal_feasibility_tolerance)
        highs_info.num_primal_infeasibilities++;
      highs_info.sum_primal_infeasibilities += primal_infeasibility;
      highs_info.max_primal_infeasibility =
          std::max(primal_infeasibility, highs_info.max_primal_infeasibility);
    }
  }

  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    double cost = lp.col_cost_[iCol];
    double dual = (HighsInt)lp.sense_ * cost;
    double lower = lp.col_lower_[iCol];
    double upper = lp.col_upper_[iCol];
    double value;
    double primal_infeasibility = 0;
    double dual_infeasibility = -1;
    HighsBasisStatus status = HighsBasisStatus::kNonbasic;
    if (lower > upper) {
      // Inconsistent bounds, so set the variable to lower bound,
      // unless it's infinite. Otherwise set the variable to upper
      // bound, unless it's infinite. Otherwise set the variable to
      // zero.
      if (highs_isInfinity(lower)) {
        // Lower bound of +inf
        if (highs_isInfinity(-upper)) {
          // Upper bound of -inf
          value = 0;
          status = HighsBasisStatus::kZero;
          primal_infeasibility = kHighsInf;
          dual_infeasibility = std::fabs(dual);
        } else {
          // Finite upper bound - since lower exceeds it
          value = upper;
          status = HighsBasisStatus::kUpper;
          primal_infeasibility = lower - value;
          dual_infeasibility = std::max(dual, 0.);
        }
      } else {
        // Finite lower bound
        value = lower;
        status = HighsBasisStatus::kLower;
        primal_infeasibility = value - upper;
        dual_infeasibility = std::max(-dual, 0.);
      }
    } else if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free column: set to zero and record dual infeasibility
      value = 0;
      status = HighsBasisStatus::kZero;
      dual_infeasibility = std::fabs(dual);
    } else if (dual >= dual_feasibility_tolerance) {
      // Column with sufficiently positive dual
      if (!highs_isInfinity(-lower)) {
        // Set to this finite lower bound
        value = lower;
        status = HighsBasisStatus::kLower;
        dual_infeasibility = 0;
      } else {
        // Infinite lower bound so set to upper bound and record dual
        // infeasibility
        value = upper;
        status = HighsBasisStatus::kUpper;
        dual_infeasibility = dual;
      }
    } else if (dual <= -dual_feasibility_tolerance) {
      // Column with sufficiently negative dual
      if (!highs_isInfinity(upper)) {
        // Set to this finite upper bound
        value = upper;
        status = HighsBasisStatus::kUpper;
        dual_infeasibility = 0;
      } else {
        // Infinite upper bound so set to lower bound and record dual
        // infeasibility
        value = lower;
        status = HighsBasisStatus::kLower;
        dual_infeasibility = -dual;
      }
    } else {
      // Column with sufficiently small dual: set to lower bound (if
      // finite) otherwise upper bound
      if (highs_isInfinity(-lower)) {
        value = upper;
        status = HighsBasisStatus::kUpper;
      } else {
        value = lower;
        status = HighsBasisStatus::kLower;
      }
      dual_infeasibility = std::fabs(dual);
    }
    assert(status != HighsBasisStatus::kNonbasic);
    assert(dual_infeasibility >= 0);
    solution.col_value[iCol] = value;
    solution.col_dual[iCol] = (HighsInt)lp.sense_ * dual;
    basis.col_status[iCol] = status;
    objective += value * cost;
    if (primal_infeasibility > primal_feasibility_tolerance)
      highs_info.num_primal_infeasibilities++;
    highs_info.sum_primal_infeasibilities += primal_infeasibility;
    highs_info.max_primal_infeasibility =
        std::max(primal_infeasibility, highs_info.max_primal_infeasibility);
    if (dual_infeasibility > dual_feasibility_tolerance)
      highs_info.num_dual_infeasibilities++;
    highs_info.sum_dual_infeasibilities += dual_infeasibility;
    highs_info.max_dual_infeasibility =
        std::max(dual_infeasibility, highs_info.max_dual_infeasibility);
  }
  highs_info.objective_function_value = objective;
  solution.value_valid = true;
  solution.dual_valid = true;
  basis.valid = true;
  basis.useful = true;
  highs_info.basis_validity = kBasisValidityValid;
  setSolutionStatus(highs_info);
  if (highs_info.num_primal_infeasibilities) {
    // Primal infeasible
    model_status = HighsModelStatus::kInfeasible;
  } else if (highs_info.num_dual_infeasibilities) {
    // Dual infeasible => primal unbounded for unconstrained LP
    model_status = HighsModelStatus::kUnbounded;
  } else {
    model_status = HighsModelStatus::kOptimal;
  }

  return HighsStatus::kOk;
}

// Assuming that any user scaling in user_scale_data has been applied,
// determine the model coefficient ranges, assess it for values
// outside the [small, large] range, and give appropriate scaling
// recommendations
void assessExcessiveObjectiveBoundScaling(const HighsOptions& options,
                                          const HighsModel& model,
                                          HighsUserScaleData& user_scale_data) {
  const HighsLogOptions& log_options = options.log_options;
  const HighsLp& lp = model.lp_;
  if (lp.num_col_ == 0 || lp.num_row_ == 0) return;
  std::vector<double> effective_cost = getEffectiveCosts(lp, options);
  const bool user_cost_or_bound_scale =
      user_scale_data.user_objective_scale || user_scale_data.user_bound_scale;
  const double small_objective_coefficient =
      kExcessivelySmallObjectiveCoefficient;
  const double large_objective_coefficient =
      kExcessivelyLargeObjectiveCoefficient;
  const double small_bound = kExcessivelySmallBoundValue;
  const double large_bound = kExcessivelyLargeBoundValue;
  const double large_integer_bound = kExcessivelyLargeIntegerBoundValue;
  std::stringstream message;
  if (user_cost_or_bound_scale) {
    if (user_scale_data.user_objective_scale)
      message << highsFormatToString(" user_objective_scale option value of %d",
                                     user_scale_data.user_objective_scale);
    if (user_scale_data.user_bound_scale) {
      if (user_scale_data.user_objective_scale) message << " and";
      message << highsFormatToString(" user_bound_scale option value of %d",
                                     user_scale_data.user_bound_scale);
    }
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Assessing costs and bounds after applying%s\n",
                 message.str().c_str());
  }
  // Lambda for assessing a finite nonzero
  auto assessFiniteNonzero = [&](const double value, double& min_value,
                                 double& max_value) {
    double abs_value = std::abs(value);
    if (abs_value > 0 && abs_value < kHighsInf) {
      min_value = std::min(abs_value, min_value);
      max_value = std::max(abs_value, max_value);
    }
  };
  double min_continuous_col_cost = kHighsInf;
  double min_noncontinuous_col_cost = kHighsInf;
  double max_continuous_col_cost = -kHighsInf;
  double max_noncontinuous_col_cost = -kHighsInf;
  double min_effective_continuous_col_cost = kHighsInf;
  double min_effective_noncontinuous_col_cost = kHighsInf;
  double max_effective_continuous_col_cost = -kHighsInf;
  double max_effective_noncontinuous_col_cost = -kHighsInf;
  double min_continuous_col_bound = kHighsInf;
  double min_noncontinuous_col_bound = kHighsInf;
  double max_continuous_col_bound = -kHighsInf;
  double max_noncontinuous_col_bound = -kHighsInf;
  double min_continuous_matrix_value = kHighsInf;
  double min_noncontinuous_matrix_value = kHighsInf;
  double max_continuous_matrix_value = -kHighsInf;
  double max_noncontinuous_matrix_value = -kHighsInf;
  const bool is_mip = lp.integrality_.size();
  HighsInt num_continuous_variable = 0;
  HighsInt num_noncontinuous_variable = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (is_mip && lp.integrality_[iCol] != HighsVarType::kContinuous) {
      assessFiniteNonzero(lp.col_cost_[iCol], min_noncontinuous_col_cost,
                          max_noncontinuous_col_cost);
      assessFiniteNonzero(effective_cost[iCol],
                          min_effective_noncontinuous_col_cost,
                          max_effective_noncontinuous_col_cost);
      assessFiniteNonzero(lp.col_lower_[iCol], min_noncontinuous_col_bound,
                          max_noncontinuous_col_bound);
      assessFiniteNonzero(lp.col_upper_[iCol], min_noncontinuous_col_bound,
                          max_noncontinuous_col_bound);
      num_noncontinuous_variable++;
    } else {
      assessFiniteNonzero(lp.col_cost_[iCol], min_continuous_col_cost,
                          max_continuous_col_cost);
      assessFiniteNonzero(effective_cost[iCol],
                          min_effective_continuous_col_cost,
                          max_effective_continuous_col_cost);
      assessFiniteNonzero(lp.col_lower_[iCol], min_continuous_col_bound,
                          max_continuous_col_bound);
      assessFiniteNonzero(lp.col_upper_[iCol], min_continuous_col_bound,
                          max_continuous_col_bound);
      num_continuous_variable++;
    }
  }
  double min_col_cost =
      std::min(min_continuous_col_cost, min_noncontinuous_col_cost);
  double max_col_cost =
      std::max(max_continuous_col_cost, max_noncontinuous_col_cost);

  double min_effective_col_cost = std::min(
      min_effective_continuous_col_cost, min_effective_noncontinuous_col_cost);
  double max_effective_col_cost = std::max(
      max_effective_continuous_col_cost, max_effective_noncontinuous_col_cost);

  double min_matrix_value = kHighsInf;
  double max_matrix_value = -kHighsInf;
  const HighsInt num_matrix_nz = lp.numNz();
  for (HighsInt iEl = 0; iEl < num_matrix_nz; iEl++)
    assessFiniteNonzero(lp.a_matrix_.value_[iEl], min_matrix_value,
                        max_matrix_value);

  double min_row_bound = kHighsInf;
  double max_row_bound = -kHighsInf;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    assessFiniteNonzero(lp.row_lower_[iRow], min_row_bound, max_row_bound);
    assessFiniteNonzero(lp.row_upper_[iRow], min_row_bound, max_row_bound);
  }

  double min_continuous_hessian_value = kHighsInf;
  double max_continuous_hessian_value = -kHighsInf;
  const HighsInt num_hessian_nz = model.hessian_.numNz();
  for (HighsInt iEl = 0; iEl < num_hessian_nz; iEl++)
    assessFiniteNonzero(model.hessian_.value_[iEl],
                        min_continuous_hessian_value,
                        max_continuous_hessian_value);

  // Determine the minimum and maximum overall bounds that can be
  // scaled with user_bound_scale before zeroing extrema due to
  // absence of finite nonzero bounds

  double min_scalable_bound = std::min(min_continuous_col_bound, min_row_bound);
  double max_scalable_bound = std::max(max_continuous_col_bound, max_row_bound);
  if (min_scalable_bound == kHighsInf) min_scalable_bound = 0;
  if (max_scalable_bound == -kHighsInf) max_scalable_bound = 0;

  if (min_col_cost == kHighsInf) min_col_cost = 0;
  if (max_col_cost == -kHighsInf) max_col_cost = 0;
  if (min_effective_col_cost == kHighsInf) min_effective_col_cost = 0;
  if (max_effective_col_cost == -kHighsInf) max_effective_col_cost = 0;
  if (min_continuous_col_bound == kHighsInf) min_continuous_col_bound = 0;
  if (max_continuous_col_bound == -kHighsInf) max_continuous_col_bound = 0;
  if (min_noncontinuous_col_bound == kHighsInf) min_noncontinuous_col_bound = 0;
  if (max_noncontinuous_col_bound == -kHighsInf)
    max_noncontinuous_col_bound = 0;
  if (min_row_bound == kHighsInf) min_row_bound = 0;
  if (max_row_bound == -kHighsInf) max_row_bound = 0;

  double min_hessian_value = min_continuous_hessian_value;
  double max_hessian_value = max_continuous_hessian_value;
  if (min_hessian_value == kHighsInf) min_hessian_value = 0;
  if (max_hessian_value == -kHighsInf) max_hessian_value = 0;

  // Report on the coefficient ranges
  highsLogUser(log_options, HighsLogType::kInfo, "Coefficient ranges:\n");
  if (num_matrix_nz)
    highsLogUser(log_options, HighsLogType::kInfo, "  Matrix  [%5.0e, %5.0e]\n",
                 min_matrix_value, max_matrix_value);
  if (lp.num_col_) {
    highsLogUser(log_options, HighsLogType::kInfo, "  Cost    [%5.0e, %5.0e]\n",
                 min_col_cost, max_col_cost);
    if (min_effective_col_cost != min_col_cost ||
        max_effective_col_cost != max_col_cost)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "  Cost    [%5.0e, %5.0e] (effective)\n",
                   min_effective_col_cost, max_effective_col_cost);
    if (num_hessian_nz)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "  Hessian [%5.0e, %5.0e]\n", min_hessian_value,
                   max_hessian_value);
    if (num_continuous_variable)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "  Bound   [%5.0e, %5.0e]%s\n", min_continuous_col_bound,
                   max_continuous_col_bound,
                   num_noncontinuous_variable > 0 ? " (continuous)" : "");
    if (num_noncontinuous_variable)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "  Bound   [%5.0e, %5.0e]%s\n", min_noncontinuous_col_bound,
                   max_noncontinuous_col_bound,
                   num_continuous_variable > 0 ? " (non-continuous)" : "");
  }
  if (lp.num_row_)
    highsLogUser(log_options, HighsLogType::kInfo, "  RHS     [%5.0e, %5.0e]\n",
                 min_row_bound, max_row_bound);

  // LPs with no columns or no finite nonzero costs will have
  // max_col_cost = 0
  assert(max_col_cost >= 0);
  // LPs with no columns or no finite nonzero bounds will have
  // max_continuous_col_bound = 0 and max_noncontinuous_col_bound = 0
  assert(max_continuous_col_bound >= 0);
  assert(max_noncontinuous_col_bound >= 0);
  // LPs with no rows or no finite nonzero bounds will have
  // max_row_bound = 0
  assert(max_row_bound >= 0);

  const std::string problem =
      user_cost_or_bound_scale ? "User-scaled problem" : "Problem";

  // If the effective costs aren't small or large, don't report costs
  // that are small or large
  if (0 < min_col_cost && min_col_cost < small_objective_coefficient &&
      0 < min_effective_col_cost &&
      min_effective_col_cost < small_objective_coefficient)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively small costs\n", problem.c_str());
  if (max_col_cost > large_objective_coefficient &&
      max_effective_col_cost > large_objective_coefficient)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively large costs\n", problem.c_str());

  // If extreme effective costs aren't different from extreme costs,
  // don't report extreme costs
  if (0 < min_effective_col_cost && min_effective_col_cost != min_col_cost &&
      min_effective_col_cost < small_objective_coefficient)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively small effective costs\n",
                 problem.c_str());
  if (max_effective_col_cost != max_col_cost &&
      max_effective_col_cost > large_objective_coefficient)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively large effective costs\n",
                 problem.c_str());

  if (0 < min_hessian_value && min_hessian_value < small_objective_coefficient)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively small Hessian values\n",
                 problem.c_str());
  if (max_hessian_value > large_objective_coefficient)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively large Hessian values\n",
                 problem.c_str());
  if (0 < min_continuous_col_bound && min_continuous_col_bound < small_bound)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively small bounds on%s variables\n",
                 problem.c_str(), is_mip ? " continuous" : "");
  if (max_continuous_col_bound > large_bound)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively large bounds on%s variables\n",
                 problem.c_str(), is_mip ? " continuous" : "");
  if (0 < min_noncontinuous_col_bound &&
      min_noncontinuous_col_bound < small_bound)
    highsLogUser(
        log_options, HighsLogType::kWarning,
        "%s has some excessively small bounds on non-continuous variables\n",
        problem.c_str());
  if (max_noncontinuous_col_bound > large_integer_bound)
    highsLogUser(
        log_options, HighsLogType::kWarning,
        "%s has some excessively large bounds on non-continuous variables\n",
        problem.c_str());
  if (0 < min_row_bound && min_row_bound < small_bound)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively small bounds on constraints\n",
                 problem.c_str());
  if (max_row_bound > large_bound)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively large bounds on constraints\n",
                 problem.c_str());

  // Lambda to determine recommended user scaling values
  auto suggestScaling = [&](double min_value, double max_value,
                            double small_value, double large_value) {
    double ratio = 1;
    if (max_value > large_value) {
      // Max scalable value is large, so suggest scaling values down
      // so that the max value is large_value
      ratio = large_value / max_value;
    } else if (0 < max_value && max_value < small_value) {
      // All scalable values are small, so suggest scaling them up so
      // the max value is small_value
      ratio = small_value / max_value;
    }
    assert(ratio);
    return ratio;
  };

  // Lambda to for outward rounding of log
  auto outerRoundedLog = [&](const double value, const HighsInt base) {
    assert(value > 0);
    assert(base == 2 || base == 10);
    HighsInt rounded_log = 0;
    if (base == 2) {
      rounded_log = value < 1 ? std::floor(std::log2(value))
                              : std::ceil(std::log2(value));
    } else {
      rounded_log = value < 1 ? std::floor(std::log10(value))
                              : std::ceil(std::log10(value));
    }
    return rounded_log;
  };
  double suggested_bound_scaling = suggestScaling(
      min_scalable_bound, max_scalable_bound, small_bound, large_bound);
  // Determine the suggested (new) value for user_bound_scale,
  // allowing for the fact that the current value has been applied
  HighsInt dl_user_bound_scale = outerRoundedLog(suggested_bound_scaling, 2);
  user_scale_data.suggested_user_bound_scale =
      user_scale_data.user_bound_scale + dl_user_bound_scale;
  // Determine the order of magnitude of the suggested bound scaling -
  // just for logging
  HighsInt suggested_bound_scale_order_of_magnitude =
      outerRoundedLog(suggested_bound_scaling, 10);
  // Applying the suggested bound scaling requires the costs and
  // matrix columns of non-continuous variables to be scaled, and any
  // Hessian entries are also scaled
  //
  // Determine the corresponding extreme non-continuous costs and
  // update the extreme costs so that objective scaling can be
  // suggested
  double suggested_user_bound_scale_value =
      pow(2.0, user_scale_data.suggested_user_bound_scale);
  min_effective_noncontinuous_col_cost *= suggested_user_bound_scale_value;
  max_effective_noncontinuous_col_cost *= suggested_user_bound_scale_value;
  min_hessian_value /= suggested_user_bound_scale_value;
  max_hessian_value /= suggested_user_bound_scale_value;

  min_effective_col_cost = std::min(min_effective_continuous_col_cost,
                                    min_effective_noncontinuous_col_cost);
  max_effective_col_cost = std::max(max_effective_continuous_col_cost,
                                    max_effective_noncontinuous_col_cost);
  double min_objective_coefficient =
      std::min(min_effective_col_cost, min_continuous_hessian_value);
  double max_objective_coefficient =
      std::max(max_effective_col_cost, max_continuous_hessian_value);
  if (min_objective_coefficient == kHighsInf) min_objective_coefficient = 0;
  if (max_objective_coefficient == -kHighsInf) max_objective_coefficient = 0;

  double suggested_objective_scaling =
      suggestScaling(min_objective_coefficient, max_objective_coefficient,
                     small_objective_coefficient, large_objective_coefficient);
  // Determine the suggested (new) value for user_objective_scale,
  // allowing for the fact that the current value has been applied
  HighsInt dl_user_objective_scale =
      outerRoundedLog(suggested_objective_scaling, 2);
  user_scale_data.suggested_user_objective_scale =
      user_scale_data.user_objective_scale + dl_user_objective_scale;
  // Determine the order of magnitude of the suggested objective scaling -
  // just for logging
  HighsInt suggested_objective_scale_order_of_magnitude =
      outerRoundedLog(suggested_objective_scaling, 10);

  bool warning_issued = false;
  // Only report the order of magnitude scaling if there is no user
  // scaling
  bool order_of_magnitude_message =
      suggested_objective_scale_order_of_magnitude &&
      !user_scale_data.user_objective_scale;
  message.str(std::string());
  if (order_of_magnitude_message)
    message << highsFormatToString(
        "   Consider scaling the objective by 1e%+1d",
        int(suggested_objective_scale_order_of_magnitude));
  if (dl_user_objective_scale) {
    if (!order_of_magnitude_message) {
      message << "   Consider";
    } else {
      message << ", or";
    }
    message << highsFormatToString(
        " setting the user_objective_scale option to %d",
        int(user_scale_data.suggested_user_objective_scale));
  }
  if (order_of_magnitude_message || dl_user_objective_scale) {
    highsLogUser(log_options, HighsLogType::kWarning, "%s\n",
                 message.str().c_str());
    warning_issued = true;
  }

  message.str(std::string());
  order_of_magnitude_message = suggested_bound_scale_order_of_magnitude &&
                               !user_scale_data.user_bound_scale;
  message.str(std::string());
  if (order_of_magnitude_message)
    message << highsFormatToString(
        "   Consider scaling the    bounds by 1e%+1d",
        int(suggested_bound_scale_order_of_magnitude));
  if (dl_user_bound_scale) {
    if (!order_of_magnitude_message) {
      message << "   Consider";
    } else {
      message << ", or";
    }
    message << highsFormatToString(
        " setting the user_bound_scale option to %d",
        int(user_scale_data.suggested_user_bound_scale));
  }
  if (order_of_magnitude_message || dl_user_bound_scale) {
    highsLogUser(log_options, HighsLogType::kWarning, "%s\n",
                 message.str().c_str());
    warning_issued = true;
  }

  // Assess the reationship between extreme objective/bound values and
  // feasibility tolerances
  const bool use_mip_feasibility_tolerances =
      model.isMip() && !options.solve_relaxation;
  const double use_primal_feasibility_tolerance =
      use_mip_feasibility_tolerances
          ? options.getMipPrimalFeasibilityTolerance()
          : options.primal_feasibility_tolerance;
  const double use_dual_feasibility_tolerance =
      use_mip_feasibility_tolerances ? options.getMipDualFeasibilityTolerance()
                                     : options.dual_feasibility_tolerance;
  const double max_objective =
      std::max(max_effective_col_cost, max_hessian_value);
  const double max_bound =
      std::max(max_continuous_col_bound,
               std::max(max_noncontinuous_col_bound, max_row_bound));
  const double primal_ratio = max_bound / use_primal_feasibility_tolerance;
  const double dual_ratio = max_objective / use_dual_feasibility_tolerance;
  const bool primal_fragile = primal_ratio > 1e-1 / kHighsMacheps;
  const bool dual_fragile = dual_ratio > 1e-1 / kHighsMacheps;
  const bool numerically_fragile = primal_fragile || dual_fragile;
  if (numerically_fragile) {
    std::string blanks(problem.length(), ' ');
    if (primal_fragile)
      highsLogUser(log_options, HighsLogType::kWarning,
                   "%s has excessively large ratio (%.1g) between max variable "
                   "or constraint bound"
                   " (%.1g) and primal feasibility tolerance (%.1g)\n",
                   problem.c_str(), primal_ratio, max_bound,
                   use_primal_feasibility_tolerance);
    if (dual_fragile)
      highsLogUser(log_options, HighsLogType::kWarning,
                   "%s has excessively large ratio (%.1g) between max "
                   "objective coefficient"
                   " (%.1g) and dual feasibility tolerance (%.1g)\n",
                   (primal_fragile ? blanks : problem).c_str(), dual_ratio,
                   max_objective, use_dual_feasibility_tolerance);
    highsLogUser(log_options, HighsLogType::kWarning,
                 "   Improve scaling and/or increase tolerance\n");
    warning_issued = true;
  }
  // Give a scary warning to users asking too much of HiGHS!
  if (warning_issued)
    highsLogUser(
        log_options, HighsLogType::kWarning,
        "%s is badly scaled, %swhich may compromise the speed, accuracy "
        "and reliability of solvers in HiGHS\n",
        problem.c_str(),
        numerically_fragile ? "particularly relative to tolerances, " : "");
}

// Decide whether to use the HiPO IPM solver
bool useHipo(const HighsOptions& options,
             const std::string& specific_solver_option, const HighsLp& lp,
             const bool logging) {
  // specific_solver_option wll be "solver", "mip_lp_solver" or
  // "mip_ipm_solver" according to context
  assert(specific_solver_option == kSolverString ||
         specific_solver_option == kMipLpSolverString ||
         specific_solver_option == kMipIpmSolverString);
  const std::string specific_solver_option_value =
      specific_solver_option == kSolverString        ? options.solver
      : specific_solver_option == kMipLpSolverString ? options.mip_lp_solver
                                                     : options.mip_ipm_solver;
  // In the MIP solver there are situations where IPM must be used
  const bool force_ipm = specific_solver_option == kMipIpmSolverString;

  // Initialize to value for valgrind.
  bool use_hipo = false;

  if (specific_solver_option_value == kIpxString) {
    use_hipo = false;
  } else if (specific_solver_option_value == kIpmString ||
             specific_solver_option_value == kHipoString || force_ipm) {
    use_hipo = HighsExternalApi::isAvailable<HighsExtras::hipo>();
  }
  if (options.run_centring) use_hipo = false;
  // Later decide between HiPO and IPX based on LP properties
  if (specific_solver_option == kMipIpmSolverString) return use_hipo;
  // Later decide between simplex, HiPO and IPX based on LP properties
  return use_hipo;
}

HighsStatus solveQp(HighsQpSolverObject& solver_object,
                    const std::string& message) {
  HighsModel& model_ = solver_object.model_;
  HighsBasis& basis = solver_object.basis_;
  HighsSolution& solution = solver_object.solution_;
  HighsInfo& info = solver_object.highs_info_;
  HighsOptions& options = solver_object.options_;
  HighsProfiling* profiling = solver_object.profiling_;
  HighsModelStatus& model_status = solver_object.model_status_;

  // Check that the model is column-wise
  HighsLp& lp = model_.lp_;
  assert(model_.lp_.a_matrix_.isColwise());
  HighsHessian& hessian = model_.hessian_;
  HighsInt dim = hessian.dim_;
  if (dim > 0) {
    assert(hessian.format_ == HessianFormat::kTriangular);
    assert(!hessian.isOracle());
  } else {
    assert(hessian.isOracle());
    assert(hessian.oracle_.isValid());
    dim = hessian.oracle_.dim_;
  }
  if (dim != lp.num_col_) {
    highsLogDev(
        options.log_options, HighsLogType::kError,
        "Hessian dimension = %d is incompatible with matrix dimension = %d\n",
        int(dim), int(lp.num_col_));
    model_status = HighsModelStatus::kModelError;
    solution.value_valid = false;
    solution.dual_valid = false;
    return HighsStatus::kError;
  }

  HighsStatus call_status;

  // Choose solver
  bool use_hipo =
      (options.solver == kHipoString || options.solver == kIpmString) &&
      HighsExternalApi::isAvailable<HighsExtras::hipo>();

  if (use_hipo) {
    // Need to convert oracle to an explicit Hessian without an oracle
    const bool was_oracle = hessian.isOracle();
    auto oracle_call = hessian.oracle_.call_;
    if (was_oracle) {
      hessian.formFromOracle();
      hessian.oracle_.call_ = nullptr;
    }
    assert(!hessian.isOracle());
    // Run HiPO
    if (profiling) profiling->start(kSubSolverHipo);
    SOLVE_CATCH_CALL(solveQpHipo(solver_object), "HiPO(QP)");
    if (profiling) profiling->stop(kSubSolverHipo);
    // Restore any oracle call;
    hessian.oracle_.call_ = oracle_call;
    assert(hessian.isOracle() == was_oracle);
  } else {
    // Run the active set QP solver
    if (profiling) profiling->start(kSubSolverQpAsm);
    SOLVE_CATCH_CALL(solveQpAsm(solver_object), "active set QP solver");
    if (profiling) profiling->stop(kSubSolverQpAsm);
  }
  HighsStatus return_status = call_status;
  if (return_status == HighsStatus::kError) return return_status;

  // Get the objective and any KKT failures
  info.objective_function_value = model_.objectiveValue(solution.col_value);
  getQpKktFailures(options, model_, solution, info);
  info.valid = true;
  if (model_status == HighsModelStatus::kOptimal)
    return checkOptimality("QP", options, info, model_status);
  return return_status;
}

HighsStatus checkOptimality(const std::string& solver_type,
                            const HighsOptions& options, const HighsInfo& info,
                            HighsModelStatus& model_status) {
  // Check for infeasibility measures incompatible with optimality
  assert(model_status == HighsModelStatus::kOptimal);
  // Cannot expect to have no dual_infeasibilities since the QP solver
  // (and, of course, the MIP solver) give no dual information
  if (info.num_primal_infeasibilities == 0 &&
      info.num_dual_infeasibilities <= 0) {
    // Consider semi-continuous infeasibilities
    if (info.num_semi_infeasibilities > 0) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "%s solver claims optimality, but with num/max/sum %d/%g/%g "
                   "semi-variable infeasibilities: consider solving with "
                   "smaller mip_feasibility_tolerance\n",
                   solver_type.c_str(), int(info.num_semi_infeasibilities),
                   info.max_semi_infeasibility, info.sum_semi_infeasibilities);
      model_status = HighsModelStatus::kSolveError;
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Setting model status to %s\n",
                   utilModelStatusToString(model_status).c_str());
      return HighsStatus::kError;
    }
    return HighsStatus::kOk;
  }
  model_status = HighsModelStatus::kSolveError;
  std::stringstream ss;
  ss.str(std::string());
  ss << highsFormatToString(
      "%s solver claims optimality, but with num/max/sum "
      "primal(%d/%g/%g)",
      solver_type.c_str(), int(info.num_primal_infeasibilities),
      info.max_primal_infeasibility, info.sum_primal_infeasibilities);
  if (info.num_dual_infeasibilities > 0)
    ss << highsFormatToString(
        "and dual(%d/%g/%g)", int(info.num_dual_infeasibilities),
        info.max_dual_infeasibility, info.sum_dual_infeasibilities);
  ss << " infeasibilities\n";
  const std::string report_string = ss.str();
  highsLogUser(options.log_options, HighsLogType::kError, "%s",
               report_string.c_str());
  highsLogUser(options.log_options, HighsLogType::kError,
               "Setting model status to %s\n",
               utilModelStatusToString(model_status).c_str());
  return HighsStatus::kError;
}

HighsStatus solveMip(HighsMipSolverObject& solver_object,
                     const string& message) {
  HighsLp& lp = solver_object.lp_;
  HighsSolution& solution = solver_object.solution_;
  std::vector<HighsObjectiveSolution>& saved_objective_and_solution =
      solver_object.saved_objective_and_solution_;
  HighsInfo& info = solver_object.highs_info_;
  HighsCallback& callback = solver_object.callback_;
  HighsOptions& options = solver_object.options_;
  HighsProfiling* profiling = solver_object.profiling_;
  HighsModelStatus& model_status = solver_object.model_status_;

  // Run the MIP solver
  HighsInt log_dev_level = options.log_dev_level;
  //  options.log_dev_level = kHighsLogDevLevelInfo;
  // Check that the model isn't row-wise
  assert(lp.a_matrix_.format_ != MatrixFormat::kRowwise);
  const bool has_semi_variables = lp.hasSemiVariables();
  HighsLp use_lp;
  if (has_semi_variables) {
    // Replace any semi-variables by a continuous/integer variable and
    // a (temporary) binary. Any initial solution must accommodate this.
    use_lp = withoutSemiVariables(lp, solution,
                                  options.primal_feasibility_tolerance);
  }
  HighsLp& mip = has_semi_variables ? use_lp : lp;
  HighsMipSolver solver(callback, options, mip, solution);
  solver.setProfiling(profiling);
  profiling->start(kSubSolverMip);
  try {
    solver.run();
  } catch (const std::exception& exception) {
    solver.modelstatus_ =
        handleExceptionIsOom(options.log_options, "MIP solver", exception)
            ? HighsModelStatus::kMemoryLimit
            : HighsModelStatus::kSolveError;
  } catch (const HighsTask::Interrupt&) {
    highsLogDev(options.log_options, HighsLogType::kError,
                "HighsTask interrupt when solving with MIP solver\n");
    solver.modelstatus_ = HighsModelStatus::kSolveError;
  }
  profiling->stop(kSubSolverMip);
  options.log_dev_level = log_dev_level;
  // Set the return_status, model status and, for completeness, scaled
  // model status
  HighsStatus return_status =
      highsStatusFromHighsModelStatus(solver.modelstatus_);
  model_status = solver.modelstatus_;
  // Extract the solution
  if (solver.solution_objective_ != kHighsInf) {
    // There is a primal solution
    //
    // If the original model has semi-variables, its solution is
    // (still) given by the first lp.num_col_ entries of the
    // solution from the MIP solver
    //
    // #2547 This resize is unnecessary
    //
    // solution.col_value.resize(lp.num_col_);
    solution.col_value = solver.solution_;
    saved_objective_and_solution = solver.saved_objective_and_solution_;
    lp.a_matrix_.productQuad(solution.row_value, solution.col_value);
    solution.value_valid = true;
  } else {
    // There is no primal solution: should be so by default
    assert(!solution.value_valid);
  }
  // Check that no modified upper bounds for semi-variables are active
  if (solution.value_valid &&
      activeModifiedUpperBounds(options, lp, solution.col_value)) {
    solution.value_valid = false;
    model_status = HighsModelStatus::kSolveError;
    return_status = HighsStatus::kError;
  }
  // There is no dual solution: should be so by default
  assert(!solution.dual_valid);
  HighsBasis basis;
  // Get the objective and any KKT failures
  info.objective_function_value = solver.solution_objective_;
  // Remember to judge primal feasibility according to
  // mip_feasibility_tolerance, so take a copy of the original
  // value...
  const double primal_feasibility_tolerance =
      options.getMipPrimalFeasibilityTolerance();
  options.primal_feasibility_tolerance = primal_feasibility_tolerance;
  // NB getKktFailures sets the primal and dual solution status
  getLpKktFailures(options, lp, solution, basis, info);
  // Recover the primal feasibility tolerance
  options.primal_feasibility_tolerance = primal_feasibility_tolerance;
  // Set the MIP-specific values of info
  info.mip_node_count = solver.node_count_;
  info.mip_dual_bound = solver.dual_bound_;
  info.mip_gap = solver.gap_;
  info.primal_dual_integral = solver.primal_dual_integral_;
  // Get the number of LP iterations, avoiding overflow if the int64_t
  // value is too large
  int64_t mip_total_lp_iterations = solver.total_lp_iterations_;
  info.simplex_iteration_count = mip_total_lp_iterations > kHighsIInf
                                     ? -1
                                     : HighsInt(mip_total_lp_iterations);
  info.valid = true;
  if (model_status == HighsModelStatus::kOptimal)
    return_status = checkOptimality("MIP", options, info, model_status);
  // Overwrite max infeasibility to include integrality if there is a solution
  if (solver.solution_objective_ != kHighsInf) {
    const double mip_max_bound_violation =
        std::max(solver.row_violation_, solver.bound_violation_);
    const double delta_max_bound_violation =
        std::abs(mip_max_bound_violation - info.max_primal_infeasibility);
    // Possibly report a mis-match between the max bound violation
    // returned by the MIP solver, and the value obtained from the
    // solution
    if (delta_max_bound_violation > 1e-12)
      highsLogDev(options.log_options, HighsLogType::kWarning,
                  "Inconsistent max bound violation: MIP solver (%10.4g); LP "
                  "(%10.4g); Difference of %10.4g\n",
                  mip_max_bound_violation, info.max_primal_infeasibility,
                  delta_max_bound_violation);
    info.max_integrality_violation = solver.integrality_violation_;
    if (info.max_integrality_violation > options.mip_feasibility_tolerance)
      info.primal_solution_status = kSolutionStatusInfeasible;
  }
  return return_status;
}
