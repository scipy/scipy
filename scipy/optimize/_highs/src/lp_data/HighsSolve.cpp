/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolve.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "lp_data/HighsInfo.h"
#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsSolution.h"
#include "simplex/HApp.h"
#include "util/HighsUtils.h"
#ifdef IPX_ON
#include "ipm/IpxWrapper.h"
#else
#include "ipm/IpxWrapperEmpty.h"
#endif

// Solves an unconstrained LP without scaling, setting HighsBasis, HighsSolution
// and HighsSolutionParams
HighsStatus solveUnconstrainedLp(HighsModelObject& highs_model_object) {
  // Reset unscaled and scaled model status and solution params - except for
  // iteration counts
  resetModelStatusAndSolutionParams(highs_model_object);

  // Aliases to unscaled model status and solution parameters
  HighsSolutionParams& unscaled_solution_params =
      highs_model_object.unscaled_solution_params_;
  //  HighsModelStatus& unscaled_model_status =
  //  highs_model_object.unscaled_model_status_;

  // Check that the LP really is unconstrained!
  const HighsLp& lp = highs_model_object.lp_;
  assert(lp.numRow_ == 0);
  if (lp.numRow_ != 0) return HighsStatus::Error;

  HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                  "Solving an unconstrained LP with %d columns", lp.numCol_);

  HighsSolution& solution = highs_model_object.solution_;
  HighsBasis& basis = highs_model_object.basis_;

  solution.col_value.assign(lp.numCol_, 0);
  solution.col_dual.assign(lp.numCol_, 0);
  basis.col_status.assign(lp.numCol_, HighsBasisStatus::NONBASIC);

  double primal_feasibility_tolerance =
      unscaled_solution_params.primal_feasibility_tolerance;
  double dual_feasibility_tolerance =
      unscaled_solution_params.dual_feasibility_tolerance;

  // Initialise the objective value calculation. Done using
  // HighsSolution so offset is vanilla
  double objective = lp.offset_;
  bool infeasible = false;
  bool unbounded = false;

  unscaled_solution_params.num_primal_infeasibilities = 0;
  unscaled_solution_params.num_dual_infeasibilities = 0;

  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    double cost = lp.colCost_[iCol];
    double dual = (int)lp.sense_ * cost;
    double lower = lp.colLower_[iCol];
    double upper = lp.colUpper_[iCol];
    double value;
    double primal_infeasibility = 0;
    HighsBasisStatus status;
    if (lower > upper) {
      // Inconsistent bounds, so set the variable to lower bound,
      // unless it's infinite. Otherwise set the variable to upper
      // bound, unless it's infinite. Otherwise set the variable to
      // zero.
      if (highs_isInfinity(lower)) {
        // Lower bound of +inf
        if (highs_isInfinity(-upper)) {
          // Unite upper bound of -inf
          value = 0;
          status = HighsBasisStatus::ZERO;
          primal_infeasibility = HIGHS_CONST_INF;
        } else {
          value = upper;
          status = HighsBasisStatus::UPPER;
          primal_infeasibility = lower - value;
        }
      } else {
        value = lower;
        status = HighsBasisStatus::LOWER;
        primal_infeasibility = value - upper;
      }
    } else if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free column: must have zero cost
      value = 0;
      status = HighsBasisStatus::ZERO;
      if (fabs(dual) > dual_feasibility_tolerance) unbounded = true;
    } else if (dual >= dual_feasibility_tolerance) {
      // Column with sufficiently positive dual: set to lower bound
      // and check for unboundedness
      if (highs_isInfinity(-lower)) unbounded = true;
      value = lower;
      status = HighsBasisStatus::LOWER;
    } else if (dual <= -dual_feasibility_tolerance) {
      // Column with sufficiently negative dual: set to upper bound
      // and check for unboundedness
      if (highs_isInfinity(upper)) unbounded = true;
      value = upper;
      status = HighsBasisStatus::UPPER;
    } else {
      // Column with sufficiently small dual: set to lower bound (if
      // finite) otherwise upper bound
      if (highs_isInfinity(-lower)) {
        value = upper;
        status = HighsBasisStatus::UPPER;
      } else {
        value = lower;
        status = HighsBasisStatus::LOWER;
      }
    }
    solution.col_value[iCol] = value;
    solution.col_dual[iCol] = (int)lp.sense_ * dual;
    basis.col_status[iCol] = status;
    objective += value * cost;
    unscaled_solution_params.sum_primal_infeasibilities += primal_infeasibility;
    if (primal_infeasibility > primal_feasibility_tolerance) {
      infeasible = true;
      unscaled_solution_params.num_primal_infeasibilities++;
      unscaled_solution_params.max_primal_infeasibility =
          max(primal_infeasibility,
              unscaled_solution_params.max_primal_infeasibility);
    }
  }
  unscaled_solution_params.objective_function_value = objective;
  basis.valid_ = true;

  if (infeasible) {
    highs_model_object.unscaled_model_status_ =
        HighsModelStatus::PRIMAL_INFEASIBLE;
    unscaled_solution_params.primal_status = STATUS_INFEASIBLE_POINT;
  } else {
    unscaled_solution_params.primal_status = STATUS_FEASIBLE_POINT;
    if (unbounded) {
      highs_model_object.unscaled_model_status_ =
          HighsModelStatus::PRIMAL_UNBOUNDED;
      unscaled_solution_params.dual_status = STATUS_UNKNOWN;
    } else {
      highs_model_object.unscaled_model_status_ = HighsModelStatus::OPTIMAL;
      unscaled_solution_params.dual_status = STATUS_FEASIBLE_POINT;
    }
  }
  highs_model_object.scaled_model_status_ =
      highs_model_object.unscaled_model_status_;
  return HighsStatus::OK;
}

// The method below runs simplex or ipx solver on the lp.
HighsStatus solveLp(HighsModelObject& model, const string message) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsOptions& options = model.options_;
  // Reset unscaled and scaled model status and solution params - except for
  // iteration counts
  resetModelStatusAndSolutionParams(model);
  HighsLogMessage(options.logfile, HighsMessageType::INFO, message.c_str());
#ifdef HIGHSDEV
  // Shouldn't have to check validity of the LP since this is done when it is
  // loaded or modified
  call_status = assessLp(model.lp_, options_);
  // If any errors have been found or normalisation carried out,
  // call_status will be ERROR or WARNING, so only valid return is OK.
  assert(call_status == HighsStatus::OK);
  return_status = interpretCallStatus(call_status, return_status, "assessLp");
  if (return_status == HighsStatus::Error) return return_status;
#endif
  if (!model.lp_.numRow_) {
    // Unconstrained LP so solve directly
    call_status = solveUnconstrainedLp(model);
    return_status =
        interpretCallStatus(call_status, return_status, "solveUnconstrainedLp");
    if (return_status == HighsStatus::Error) return return_status;
  } else if (options.solver == ipm_string) {
    // Use IPM
#ifdef IPX_ON
    bool imprecise_solution;
    call_status = solveLpIpx(
        options, model.timer_, model.lp_, imprecise_solution, model.basis_,
        model.solution_, model.iteration_counts_, model.unscaled_model_status_,
        model.unscaled_solution_params_);
    return_status =
        interpretCallStatus(call_status, return_status, "solveLpIpx");
    if (return_status == HighsStatus::Error) return return_status;
    if (imprecise_solution) {
      // IPX+crossover has not obtained a solution satisfying the tolerances.
      // Use the simplex method to clean up
      call_status = solveLpSimplex(model);
      return_status =
          interpretCallStatus(call_status, return_status, "solveLpSimplex");
      if (return_status == HighsStatus::Error) return return_status;

      if (!isSolutionRightSize(model.lp_, model.solution_)) {
        HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                        "Inconsistent solution returned from solver");
        return HighsStatus::Error;
      }
    } else {
      // Set the scaled model status and solution params for completeness
      model.scaled_model_status_ = model.unscaled_model_status_;
      model.scaled_solution_params_ = model.unscaled_solution_params_;
    }
#else
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Model cannot be solved with IPM");
    return HighsStatus::Error;
#endif
  } else {
    // Use Simplex
    call_status = solveLpSimplex(model);
    return_status =
        interpretCallStatus(call_status, return_status, "solveLpSimplex");
    if (return_status == HighsStatus::Error) return return_status;

    if (!isSolutionRightSize(model.lp_, model.solution_)) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Inconsistent solution returned from solver");
      return HighsStatus::Error;
    }
  }
  // Possibly analyse the HiGHS basic solution
  debugHighsBasicSolution(message, model);

  return return_status;
}
