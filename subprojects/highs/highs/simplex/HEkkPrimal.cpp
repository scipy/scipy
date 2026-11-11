/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkPrimal.cpp
 * @brief
 */
#include "simplex/HEkkPrimal.h"

#include "../extern/pdqsort/pdqsort.h"
#include "simplex/HEkkDual.h"
#include "simplex/SimplexTimer.h"
#include "util/HighsSort.h"

using std::min;

HighsStatus HEkkPrimal::solve(const bool pass_force_phase2) {
  // Initialise control data for a particular solve
  initialiseSolve();
  // Assumes that the LP has a positive number of rows
  if (ekk_instance_.isUnconstrainedLp())
    return ekk_instance_.returnFromSolve(HighsStatus::kError);

  HighsOptions& options = *ekk_instance_.options_;
  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsSimplexStatus& status = ekk_instance_.status_;

  if (!status.has_invert) {
    highsLogDev(options.log_options, HighsLogType::kError,
                "HEkkPrimal::solve called without INVERT\n");
    assert(status.has_fresh_invert);
    return ekk_instance_.returnFromSolve(HighsStatus::kError);
  }

  if (debugPrimalSimplex("Initialise", true) == HighsDebugStatus::kLogicalError)
    return ekk_instance_.returnFromSolve(HighsStatus::kError);

  // Get the nonbasic free column set
  getNonbasicFreeColumnSet();

  const bool primal_feasible_with_unperturbed_bounds =
      info.num_primal_infeasibilities == 0;
  const bool force_phase2 =
      pass_force_phase2 ||
      info.max_primal_infeasibility * info.max_primal_infeasibility <
          options.primal_feasibility_tolerance;
  // Determine whether the solution is near-optimal. Values 1000 and
  // 1e-3 (ensuring sum<1) are unimportant, as the sum of dual
  // infeasibilities for near-optimal solutions is typically many
  // orders of magnitude smaller than 1, and the sum of dual
  // infeasibilities will be very much larger for non-trivial LPs
  // that are primal feasible for a logical or crash basis.
  //
  // Consider there to be no primal infeasibilities if there are none,
  // or if phase 2 is forced, in which case any primal infeasibilities
  // will be shifted
  const bool no_simplex_primal_infeasibilities =
      primal_feasible_with_unperturbed_bounds || force_phase2;
  const bool near_optimal = info.num_dual_infeasibilities < 1000 &&
                            info.max_dual_infeasibility < 1e-3 &&
                            no_simplex_primal_infeasibilities;
  // For reporting, save primal infeasibility data for the LP without
  // bound perturbations
  const HighsInt unperturbed_num_infeasibilities =
      info.num_primal_infeasibilities;
  const double unperturbed_max_infeasibility = info.max_primal_infeasibility;
  const double unperturbed_sum_infeasibilities =
      info.sum_primal_infeasibilities;
  if (near_optimal)
    highsLogDev(options.log_options, HighsLogType::kDetailed,
                "Primal feasible and num / max / sum "
                "dual infeasibilities of "
                "%" HIGHSINT_FORMAT
                " / %g "
                "/ %g, so near-optimal\n",
                info.num_dual_infeasibilities, info.max_dual_infeasibility,
                info.sum_dual_infeasibilities);

  // Perturb bounds according to whether the solution is near-optimal
  const bool perturb_bounds = !near_optimal;
  if (!perturb_bounds)
    highsLogDev(options.log_options, HighsLogType::kDetailed,
                "Near-optimal, so don't use bound perturbation\n");
  if (perturb_bounds && info.primal_simplex_bound_perturbation_multiplier) {
    ekk_instance_.initialiseBound(SimplexAlgorithm::kPrimal, kSolvePhaseUnknown,
                                  perturb_bounds);
    ekk_instance_.initialiseNonbasicValueAndMove();
    ekk_instance_.computePrimal();
    ekk_instance_.computeSimplexPrimalInfeasible();
  }

  // Check whether the time/iteration limit has been reached. First
  // point at which a non-error return can occur
  if (ekk_instance_.bailout())
    return ekk_instance_.returnFromSolve(HighsStatus::kWarning);

  // Now to do some iterations!
  HighsInt num_primal_infeasibility =
      ekk_instance_.info_.num_primal_infeasibilities;
  solve_phase = num_primal_infeasibility > 0 ? kSolvePhase1 : kSolvePhase2;
  if (force_phase2) {
    // Dual infeasibilities without cost perturbation involved
    // fixed variables or were (at most) small, so can easily be
    // removed by flips for and fixed variables shifts for the rest
    solve_phase = kSolvePhase2;
    if (!pass_force_phase2) {
      const bool local_report = false;  // true;
      if (!primal_feasible_with_unperturbed_bounds && local_report) {
        printf(
            "Solve %d: Forcing phase 2 since near primal feasible with "
            "unperturbed "
            "costs\n"
            "num / max / sum primal infeasibilities\n"
            "%d / %11.4g / %11.4g (  perturbed bounds)\n"
            "%d / %11.4g / %11.4g (unperturbed bounds)\n",
            (int)ekk_instance_.debug_solve_call_num_,
            (int)info.num_primal_infeasibilities, info.max_primal_infeasibility,
            info.sum_primal_infeasibilities,
            (int)unperturbed_num_infeasibilities, unperturbed_max_infeasibility,
            unperturbed_sum_infeasibilities);
      }
    }
  }
  if (ekk_instance_.debugOkForSolve(algorithm, solve_phase) ==
      HighsDebugStatus::kLogicalError)
    return ekk_instance_.returnFromSolve(HighsStatus::kError);
  // Resize the copy of scattered edge weights for backtracking
  info.backtracking_basis_edge_weight_.resize(num_tot);

  // The major solving loop
  //
  // Possibly write out the column header for iteration reports, and
  // initialise records for primal correction reporting
  localReportIter(true);
  correctPrimal(true);
  while (solve_phase) {
    HighsInt it0 = ekk_instance_.iteration_count_;
    // When starting a new phase the (updated) primal objective function
    // value isn't known. Indicate this so that when the value
    // computed from scratch in rebuild() isn't checked against the
    // updated value
    status.has_primal_objective_value = false;
    if (solve_phase == kSolvePhaseUnknown) {
      // Determine the number of primal infeasibilities, and hence the solve
      // phase
      ekk_instance_.computeSimplexPrimalInfeasible();
      num_primal_infeasibility = ekk_instance_.info_.num_primal_infeasibilities;
      solve_phase = num_primal_infeasibility > 0 ? kSolvePhase1 : kSolvePhase2;
      if (info.backtracking_) {
        // Backtracking
        ekk_instance_.initialiseCost(SimplexAlgorithm::kPrimal, solve_phase);
        ekk_instance_.initialiseNonbasicValueAndMove();
        // Can now forget that we might have been backtracking
        info.backtracking_ = false;
      }
    }
    assert(solve_phase == kSolvePhase1 || solve_phase == kSolvePhase2);
    if (solve_phase == kSolvePhase1) {
      //
      // Phase 1
      //
      // solve_phase = kSolvePhase1 if the iteration or time limit has
      // been reached
      //
      // solve_phase = kSolvePhase2 if there are no primal infeasibilities
      //
      // solve_phase = kSolvePhaseUnknown if backtracking
      //
      // solve_phase = kSolvePhaseExit if primal infeasibility is
      // detected, in which case model_status_ =
      // HighsModelStatus::kInfeasible is set
      //
      // solve_phase = kSolvePhaseTabooBasis is set if only basis change is
      // taboo
      //
      // solve_phase = kSolvePhaseError is set if an error occurs
      solvePhase1();
      assert(solve_phase == kSolvePhase1 || solve_phase == kSolvePhase2 ||
             solve_phase == kSolvePhaseUnknown ||
             solve_phase == kSolvePhaseExit ||
             solve_phase == kSolvePhaseTabooBasis ||
             solve_phase == kSolvePhaseError);
      info.primal_phase1_iteration_count +=
          (ekk_instance_.iteration_count_ - it0);
    } else if (solve_phase == kSolvePhase2) {
      //
      // Phase 2
      //
      // solve_phase = kSolvePhaseOptimal if there are no dual
      // infeasibilities
      //
      // solve_phase = kSolvePhase1 if there are primal
      // infeasibilities
      //
      // solve_phase = kSolvePhase2 if the iteration or time limit has
      // been reached
      //
      // solve_phase = kSolvePhaseOptimalCleanup if, after removing bound
      // shifts, there are primal infeasibilities to clean up
      //
      // solve_phase = kSolvePhaseUnknown if backtracking
      //
      // solve_phase = kSolvePhaseExit if primal unboundedness is
      // detected, in which case model_status_ =
      // HighsModelStatus::kUnbounded is set
      //
      // solve_phase = kSolvePhaseTabooBasis is set if only basis change is
      // taboo
      //
      // solve_phase = kSolvePhaseError is set if an error occurs
      solvePhase2();
      assert(solve_phase == kSolvePhaseOptimal || solve_phase == kSolvePhase1 ||
             solve_phase == kSolvePhase2 ||
             solve_phase == kSolvePhaseOptimalCleanup ||
             solve_phase == kSolvePhaseUnknown ||
             solve_phase == kSolvePhaseExit ||
             solve_phase == kSolvePhaseTabooBasis ||
             solve_phase == kSolvePhaseError);
      assert(solve_phase != kSolvePhaseExit ||
             ekk_instance_.model_status_ == HighsModelStatus::kUnbounded);
      info.primal_phase2_iteration_count +=
          (ekk_instance_.iteration_count_ - it0);
    } else {
      // Should only be kSolvePhase1 or kSolvePhase2
      ekk_instance_.model_status_ = HighsModelStatus::kSolveError;
      return ekk_instance_.returnFromSolve(HighsStatus::kError);
    }
    // Return if bailing out from solve
    if (ekk_instance_.solve_bailout_)
      return ekk_instance_.returnFromSolve(HighsStatus::kWarning);
    // Can have all possible cases of solve_phase
    assert(solve_phase >= kSolvePhaseMin && solve_phase <= kSolvePhaseMax);
    // Look for scenarios when the major solving loop ends
    if (solve_phase == kSolvePhaseTabooBasis) {
      // Only basis change is taboo so return HighsStatus::kWarning
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "HEkkPrimal::solve Only basis change is taboo\n");
      ekk_instance_.model_status_ = HighsModelStatus::kUnknown;
      return ekk_instance_.returnFromSolve(HighsStatus::kWarning);
    }
    if (solve_phase == kSolvePhaseError) {
      // Solver error so return HighsStatus::kError
      ekk_instance_.model_status_ = HighsModelStatus::kSolveError;
      return ekk_instance_.returnFromSolve(HighsStatus::kError);
    }
    if (solve_phase == kSolvePhaseExit) {
      // LP identified as not having an optimal solution
      assert(ekk_instance_.model_status_ == HighsModelStatus::kInfeasible ||
             ekk_instance_.model_status_ == HighsModelStatus::kUnbounded);
      // If infeasible, save the primal phase 1 dual values before
      // they are overwritten with the duals for the original
      // objective
      if (ekk_instance_.model_status_ == HighsModelStatus::kInfeasible)
        ekk_instance_.primal_phase1_dual_ = ekk_instance_.info_.workDual_;
      break;
    }
    if (solve_phase == kSolvePhaseOptimalCleanup) {
      // Primal infeasibilities after phase 2. Dual feasible with
      // primal infeasibilities so use dual simplex to clean up
      break;
    }
    // If solve_phase == kSolvePhaseOptimal == 0 then major solving
    // loop ends naturally since solve_phase is false
  }
  // If bailing out, should have returned already
  assert(!ekk_instance_.solve_bailout_);
  // Should only have these cases
  assert(solve_phase ==
             kSolvePhaseExit ||  // solve_phase == kSolvePhaseUnknown ||
         solve_phase == kSolvePhaseOptimal ||  // solve_phase == kSolvePhase1 ||
         solve_phase == kSolvePhaseOptimalCleanup);
  if (solve_phase == kSolvePhaseOptimal)
    ekk_instance_.model_status_ = HighsModelStatus::kOptimal;

  if (solve_phase == kSolvePhaseOptimalCleanup) {
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "HEkkPrimal:: Using dual simplex to try to clean up num / "
                "max / sum = %" HIGHSINT_FORMAT
                " / %g / %g primal infeasibilities\n",
                info.num_primal_infeasibilities, info.max_primal_infeasibility,
                info.sum_primal_infeasibilities);
    ekk_instance_.computePrimalObjectiveValue();
    // Use dual to clean up. This almost always yields optimality,
    // and shouldn't yield infeasibility - since the current point
    // is dual feasible - but can yield
    // unboundedness. Time/iteration limit return is, of course,
    // possible, as are solver error
    HighsStatus return_status = HighsStatus::kOk;
    // Switch off any bound perturbation
    double save_dual_simplex_cost_perturbation_multiplier =
        info.dual_simplex_cost_perturbation_multiplier;
    info.dual_simplex_cost_perturbation_multiplier = 0;
    HighsInt simplex_strategy = info.simplex_strategy;
    info.simplex_strategy = kSimplexStrategyDualPlain;
    HEkkDual dual_solver(ekk_instance_);
    HighsStatus call_status = dual_solver.solve(true);
    // Restore any bound perturbation
    info.dual_simplex_cost_perturbation_multiplier =
        save_dual_simplex_cost_perturbation_multiplier;
    info.simplex_strategy = simplex_strategy;
    assert(ekk_instance_.called_return_from_solve_);
    return_status = interpretCallStatus(options.log_options, call_status,
                                        return_status, "HEkkDual::solve");
    // Reset called_return_from_solve_ to be false, since it's
    // called for this solve
    ekk_instance_.called_return_from_solve_ = false;
    if (return_status != HighsStatus::kOk)
      return ekk_instance_.returnFromSolve(return_status);
    if (ekk_instance_.model_status_ == HighsModelStatus::kOptimal &&
        info.num_primal_infeasibilities + info.num_dual_infeasibilities)
      highsLogDev(options.log_options, HighsLogType::kWarning,
                  "HEkkPrimal:: Dual simplex clean up yields  optimality, but "
                  "with %" HIGHSINT_FORMAT
                  " (max %g) primal infeasibilities and %" HIGHSINT_FORMAT
                  " (max %g) dual infeasibilities\n",
                  info.num_primal_infeasibilities,
                  info.max_primal_infeasibility, info.num_dual_infeasibilities,
                  info.max_dual_infeasibility);
  }
  if (ekk_instance_.debugOkForSolve(algorithm, solve_phase) ==
      HighsDebugStatus::kLogicalError)
    return ekk_instance_.returnFromSolve(HighsStatus::kError);
  return ekk_instance_.returnFromSolve(HighsStatus::kOk);
}

void HEkkPrimal::initialiseInstance() {
  // Called in constructor for HEkkPrimal class
  analysis = &ekk_instance_.analysis_;

  num_col = ekk_instance_.lp_.num_col_;
  num_row = ekk_instance_.lp_.num_row_;
  num_tot = num_col + num_row;

  // Setup local vectors
  col_aq.setup(num_row);
  row_ep.setup(num_row);
  row_ap.setup(num_col);
  col_basic_feasibility_change.setup(num_row);
  row_basic_feasibility_change.setup(num_col);
  col_steepest_edge.setup(num_row);

  ph1SorterR.reserve(num_row);
  ph1SorterT.reserve(num_row);

  num_free_col = 0;
  for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
    if (ekk_instance_.info_.workLower_[iCol] == -kHighsInf &&
        ekk_instance_.info_.workUpper_[iCol] == kHighsInf) {
      // Free column
      num_free_col++;
    }
  }
  // Set up the HSet instances, possibly using the internal error reporting and
  // debug option
  const bool debug =
      ekk_instance_.options_->highs_debug_level > kHighsDebugLevelCheap;
  if (num_free_col) {
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "HEkkPrimal:: LP has %" HIGHSINT_FORMAT " free columns\n",
                num_free_col);
    nonbasic_free_col_set.setup(
        num_free_col, num_tot, ekk_instance_.options_->output_flag,
        ekk_instance_.options_->log_options.log_stream, debug);
  }
  // Set up the hyper-sparse CHUZC data
  hyper_chuzc_candidate.resize(1 + max_num_hyper_chuzc_candidates);
  hyper_chuzc_measure.resize(1 + max_num_hyper_chuzc_candidates);
  hyper_chuzc_candidate_set.setup(
      max_num_hyper_chuzc_candidates, num_tot,
      ekk_instance_.options_->output_flag,
      ekk_instance_.options_->log_options.log_stream, debug);
}

void HEkkPrimal::initialiseSolve() {
  // Copy values of simplex solver options to dual simplex options
  primal_feasibility_tolerance =
      ekk_instance_.options_->primal_feasibility_tolerance;
  dual_feasibility_tolerance =
      ekk_instance_.options_->dual_feasibility_tolerance;
  objective_target = ekk_instance_.options_->objective_target;

  ekk_instance_.status_.has_primal_objective_value = false;
  ekk_instance_.status_.has_dual_objective_value = false;

  ekk_instance_.model_status_ = HighsModelStatus::kNotset;
  ekk_instance_.solve_bailout_ = false;
  ekk_instance_.called_return_from_solve_ = false;
  ekk_instance_.exit_algorithm_ = SimplexAlgorithm::kPrimal;

  rebuild_reason = kRebuildReasonNo;
  if (!ekk_instance_.status_.has_dual_steepest_edge_weights) {
    // No dual weights to maintain, so ensure that the vectors are
    // assigned since they are used around factorization and when
    // setting up the backtracking information. ToDo Eliminate this
    // opacity
    ekk_instance_.dual_edge_weight_.assign(num_row, 1.0);
    ekk_instance_.scattered_dual_edge_weight_.resize(num_tot);
  }
  const HighsInt edge_weight_strategy =
      ekk_instance_.options_->simplex_primal_edge_weight_strategy;
  if (edge_weight_strategy == kSimplexEdgeWeightStrategyChoose ||
      edge_weight_strategy == kSimplexEdgeWeightStrategyDevex) {
    // By default, use Devex
    edge_weight_mode = EdgeWeightMode::kDevex;
  } else if (edge_weight_strategy == kSimplexEdgeWeightStrategyDantzig) {
    edge_weight_mode = EdgeWeightMode::kDantzig;
  } else {
    assert(edge_weight_strategy == kSimplexEdgeWeightStrategySteepestEdge);
    edge_weight_mode = EdgeWeightMode::kSteepestEdge;
  }
  if (edge_weight_mode == EdgeWeightMode::kDantzig) {
    edge_weight_.assign(num_tot, 1.0);
  } else if (edge_weight_mode == EdgeWeightMode::kDevex) {
    initialiseDevexFramework();
  } else if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
    computePrimalSteepestEdgeWeights();
  }
}

void HEkkPrimal::solvePhase1() {
  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsSimplexStatus& status = ekk_instance_.status_;
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the
  // updated value
  status.has_primal_objective_value = false;
  status.has_dual_objective_value = false;
  // Possibly bail out immediately if iteration limit is current value
  if (ekk_instance_.bailout()) return;
  highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kDetailed,
              "primal-phase1-start\n");
  // If there's no backtracking basis, save the initial basis in case of
  // backtracking
  if (!info.valid_backtracking_basis_) ekk_instance_.putBacktrackingBasis();

  // Main solving structure
  for (;;) {
    //
    // Rebuild
    //
    // solve_phase = kSolvePhaseError is set if the basis matrix is singular
    rebuild();
    if (solve_phase == kSolvePhaseError) return;
    if (solve_phase == kSolvePhaseUnknown) return;
    if (ekk_instance_.bailout()) return;
    assert(solve_phase == kSolvePhase1 || solve_phase == kSolvePhase2);
    //
    // solve_phase = kSolvePhase2 is set if no primal infeasibilities
    // are found in rebuild(), in which case return for phase 2
    if (solve_phase == kSolvePhase2) break;

    for (;;) {
      iterate();
      if (ekk_instance_.bailout()) return;
      if (solve_phase == kSolvePhaseError) return;
      assert(solve_phase == kSolvePhase1);
      if (rebuild_reason) break;
    }
    // If the data are fresh from rebuild() and no flips have
    // occurred, possibly break out of the outer loop to see what's
    // occurred
    bool finished = status.has_fresh_rebuild && num_flip_since_rebuild == 0 &&
                    !ekk_instance_.rebuildRefactor(rebuild_reason);
    if (finished && ekk_instance_.tabooBadBasisChange()) {
      // A bad basis change has had to be made taboo without any other
      // basis changes or flips having been performed from a fresh
      // rebuild. In other words, the only basis change that could be
      // made is not permitted, so no definitive statement about the
      // LP can be made.
      solve_phase = kSolvePhaseTabooBasis;
      return;
    }
    if (finished) break;
  }
  // If bailing out, should have returned already
  assert(!ekk_instance_.solve_bailout_);
  // Will only have accurate simplex info if moving to phase 2 - but
  // should check primal feasibility and residual information if LP
  // is primal infeasible
  if (debugPrimalSimplex("End of solvePhase1") ==
      HighsDebugStatus::kLogicalError) {
    solve_phase = kSolvePhaseError;
    return;
  }
  if (solve_phase == kSolvePhase1) {
    // Determine whether primal infeasibility has been identified
    if (variable_in < 0) {
      // Optimal in phase 1, so should have primal infeasibilities
      assert(info.num_primal_infeasibilities > 0);
      if (ekk_instance_.info_.bounds_shifted ||
          ekk_instance_.info_.bounds_perturbed) {
        // Remove any bound shifts or perturbations and return to
        // phase 1
        cleanup();
      } else {
        ekk_instance_.model_status_ = HighsModelStatus::kInfeasible;
        solve_phase = kSolvePhaseExit;
      }
    }
  }
  if (solve_phase == kSolvePhase2) {
    // Moving to phase 2 so comment if bound perturbation is not permitted
    //
    // It may have been prevented to avoid cleanup-perturbation loops
    if (!info.allow_bound_perturbation)
      highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kWarning,
                  "Moving to phase 2, but not allowing bound perturbation\n");
  }
}

void HEkkPrimal::solvePhase2() {
  HighsOptions& options = *ekk_instance_.options_;
  HighsSimplexStatus& status = ekk_instance_.status_;
  HighsModelStatus& model_status = ekk_instance_.model_status_;
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the
  // updated value
  status.has_primal_objective_value = false;
  status.has_dual_objective_value = false;
  // Possibly bail out immediately if iteration limit is current value
  if (ekk_instance_.bailout()) return;
  highsLogDev(options.log_options, HighsLogType::kDetailed,
              "primal-phase2-start\n");
  phase2UpdatePrimal(true);

  // If there's no backtracking basis Save the initial basis in case of
  // backtracking
  if (!ekk_instance_.info_.valid_backtracking_basis_)
    ekk_instance_.putBacktrackingBasis();

  // Main solving structure
  for (;;) {
    //
    // Rebuild
    //
    // solve_phase = kSolvePhaseError is set if the basis matrix is singular
    rebuild();
    if (solve_phase == kSolvePhaseError) return;
    if (solve_phase == kSolvePhaseUnknown) return;
    if (ekk_instance_.bailout()) return;
    assert(solve_phase == kSolvePhase1 || solve_phase == kSolvePhase2);
    //
    // solve_phase = kSolvePhase1 is set if primal infeasibilities
    // are found in rebuild(), in which case return for phase 1
    if (solve_phase == kSolvePhase1) break;

    for (;;) {
      iterate();
      if (ekk_instance_.bailout()) return;
      if (solve_phase == kSolvePhaseError) return;
      assert(solve_phase == kSolvePhase2);
      if (rebuild_reason) break;
    }
    // If the data are fresh from rebuild() and no flips have
    // occurred, possibly break out of the outer loop to see what's
    // occurred
    bool finished = status.has_fresh_rebuild && num_flip_since_rebuild == 0 &&
                    !ekk_instance_.rebuildRefactor(rebuild_reason);
    if (finished && ekk_instance_.tabooBadBasisChange()) {
      // A bad basis change has had to be made taboo without any other
      // basis changes or flips having been performed from a fresh
      // rebuild. In other words, the only basis change that could be
      // made is not permitted, so no definitive statement about the
      // LP can be made.
      solve_phase = kSolvePhaseTabooBasis;
      return;
    }
    if (finished) break;
  }
  // If bailing out, should have returned already
  assert(!ekk_instance_.solve_bailout_);
  if (debugPrimalSimplex("End of solvePhase2") ==
      HighsDebugStatus::kLogicalError) {
    solve_phase = kSolvePhaseError;
    return;
  }
  if (solve_phase == kSolvePhase1) {
    highsLogDev(options.log_options, HighsLogType::kDetailed,
                "primal-return-phase1\n");
  } else if (variable_in == -1) {
    // There is no candidate in CHUZC, even after rebuild so probably optimal
    highsLogDev(options.log_options, HighsLogType::kDetailed,
                "primal-phase-2-optimal\n");
    // Remove any bound perturbations and see if basis is still primal feasible
    cleanup();
    if (ekk_instance_.info_.num_primal_infeasibilities > 0) {
      // There are primal infeasibilities, so consider performing dual
      // simplex iterations to get primal feasibility
      solve_phase = kSolvePhaseOptimalCleanup;
    } else {
      // There are no primal infeasibilities so optimal!
      solve_phase = kSolvePhaseOptimal;
      highsLogDev(options.log_options, HighsLogType::kDetailed,
                  "problem-optimal\n");
      model_status = HighsModelStatus::kOptimal;
      ekk_instance_.computeDualObjectiveValue();  // Why?
    }
  } else if (row_out == kNoRowSought) {
    // CHUZR has not been performed - because the chosen reduced cost
    // was unattractive when computed from scratch and no rebuild was
    // required. This is very rare and should be handled otherwise
    //
    printf("HEkkPrimal::solvePhase2 row_out = %d solve %d\n", (int)row_out,
           (int)ekk_instance_.debug_solve_call_num_);
    fflush(stdout);
    assert(row_out != kNoRowSought);
  } else {
    // No candidate in CHUZR
    if (row_out >= 0) {
      printf("HEkkPrimal::solvePhase2 row_out = %d solve %d\n", (int)row_out,
             (int)ekk_instance_.debug_solve_call_num_);
      fflush(stdout);
    }
    // Ensure that CHUZR was performed and found no row
    assert(row_out == kNoRowChosen);

    // There is no candidate in CHUZR, so probably primal unbounded
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "primal-phase-2-unbounded\n");
    if (ekk_instance_.info_.bounds_shifted ||
        ekk_instance_.info_.bounds_perturbed) {
      // If the bounds have been shifted or perturbed, clean up and
      // return
      cleanup();
      // If there are primal infeasibilities, go back to phase 1
      if (ekk_instance_.info_.num_primal_infeasibilities > 0)
        solve_phase = kSolvePhase1;
    } else {
      // The bounds have not been perturbed, so primal unbounded
      solve_phase = kSolvePhaseExit;
      // Primal unbounded, so save primal ray
      savePrimalRay();
      // Model status should be unset
      assert(model_status == HighsModelStatus::kNotset);
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "problem-primal-unbounded\n");
      model_status = HighsModelStatus::kUnbounded;
    }
  }
}

void HEkkPrimal::cleanup() {
  HighsSimplexInfo& info = ekk_instance_.info_;
  if (!info.bounds_shifted && !info.bounds_perturbed) return;
  highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kDetailed,
              "primal-cleanup-shift\n");
  // Remove perturbation and don't permit further perturbation
  ekk_instance_.initialiseBound(SimplexAlgorithm::kPrimal, solve_phase, false);
  ekk_instance_.initialiseNonbasicValueAndMove();
  info.allow_bound_perturbation = false;
  // Possibly take a copy of the original duals before recomputing them
  /*
  vector<double> original_baseValue;
  if (ekk_instance_.options_->highs_debug_level > kHighsDebugLevelCheap)
    original_baseValue = info.baseValue_;
  */
  // Compute the primal values
  ekk_instance_.computePrimal();
  // Possibly analyse the change in duals
  /*  debugCleanup(ekk_instance_, original_baseValue); */
  // Compute the primal infeasibilities
  ekk_instance_.computeSimplexPrimalInfeasible();

  // Compute the primal objective value
  ekk_instance_.computePrimalObjectiveValue();
  // Now that there's a new primal_objective_value, reset the updated
  // value
  info.updated_primal_objective_value = info.primal_objective_value;

  //  if (!info.run_quiet) {
  // Report the dual infeasibilities
  ekk_instance_.computeSimplexDualInfeasible();
  // In phase 1, report the simplex LP dual infeasibilities
  // In phase 2, report the simplex dual infeasibilities (known)
  //    if (solve_phase == kSolvePhase1)
  //    computeSimplexLpDualInfeasible(ekk_instance_);
  reportRebuild(kRebuildReasonCleanup);
  //  }
}

void HEkkPrimal::rebuild() {
  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsSimplexStatus& status = ekk_instance_.status_;
  // Clear taboo flag from any bad basis changes
  ekk_instance_.clearBadBasisChangeTabooFlag();

  // Record whether the update objective value should be tested. If
  // the objective value is known, then the updated objective value
  // should be correct - once the correction due to recomputing the
  // dual values has been applied.
  //
  // Note that computePrimalObjectiveValue sets
  // has_primal_objective_value
  //
  // Have to do this before INVERT, as this permutes the indices of
  // basic variables, and baseValue only corresponds to the new
  // ordering once computePrimal has been called
  const bool check_updated_objective_value = status.has_primal_objective_value;
  double previous_primal_objective_value = -kHighsInf;
  if (check_updated_objective_value) {
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase,
    //    "Before INVERT");
    previous_primal_objective_value = info.updated_primal_objective_value;
  } else {
    // Reset the knowledge of previous objective values
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, -1, "");
  }

  // Decide whether refactorization should be performed
  const bool refactor_basis_matrix =
      ekk_instance_.rebuildRefactor(rebuild_reason);

  // Take a local copy of the rebuild reason and then reset the global value
  const HighsInt local_rebuild_reason = rebuild_reason;
  rebuild_reason = kRebuildReasonNo;
  if (refactor_basis_matrix) {
    // Get a nonsingular inverse if possible. One of three things
    // happens: Current basis is nonsingular; Current basis is
    // singular and last nonsingular basis is refactorized as
    // nonsingular - or found singular. Latter is code failure.
    if (!ekk_instance_.getNonsingularInverse(solve_phase)) {
      solve_phase = kSolvePhaseError;
      return;
    }
    // Record the synthetic clock for INVERT, and zero it for UPDATE
    ekk_instance_.resetSyntheticClock();
  }
  if (!ekk_instance_.status_.has_ar_matrix) {
    // Don't have the row-wise matrix, so reinitialise it
    //
    // Should only happen when backtracking
    assert(info.backtracking_);
    ekk_instance_.initialisePartitionedRowwiseMatrix();
  }
  assert(ekk_instance_.ar_matrix_.debugPartitionOk(
      ekk_instance_.basis_.nonbasicFlag_.data()));

  if (info.backtracking_) {
    // If backtracking, may change phase, so drop out
    solve_phase = kSolvePhaseUnknown;
    return;
  }

  ekk_instance_.computePrimal();
  if (solve_phase == kSolvePhase2) {
    bool correct_primal_ok = correctPrimal();
    if (kAllowDeveloperAssert) {
      assert(correct_primal_ok);
    }
  }
  getBasicPrimalInfeasibility();
  if (info.num_primal_infeasibilities > 0) {
    // Primal infeasibilities so should be in phase 1
    if (solve_phase == kSolvePhase2) {
      highsLogDev(
          ekk_instance_.options_->log_options, HighsLogType::kWarning,
          "HEkkPrimal::rebuild switching back to phase 1 from phase 2\n");
      solve_phase = kSolvePhase1;
    }
    phase1ComputeDual();
  } else {
    // No primal infeasibilities so in phase 2. Reset costs if was
    // previously in phase 1
    if (solve_phase == kSolvePhase1) {
      ekk_instance_.initialiseCost(SimplexAlgorithm::kPrimal, solve_phase);
      solve_phase = kSolvePhase2;
    }
    ekk_instance_.computeDual();
  }
  ekk_instance_.computeSimplexDualInfeasible();
  ekk_instance_.computePrimalObjectiveValue();
  if (check_updated_objective_value) {
    // Apply the objective value correction due to computing primal
    // values from scratch.
    const double primal_objective_value_correction =
        info.primal_objective_value - previous_primal_objective_value;
    info.updated_primal_objective_value += primal_objective_value_correction;
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm);
  }
  // Now that there's a new dual_objective_value, reset the updated
  // value
  info.updated_primal_objective_value = info.primal_objective_value;

  reportRebuild(local_rebuild_reason);

  // Record the synthetic clock for INVERT, and zero it for UPDATE
  ekk_instance_.resetSyntheticClock();

  // Determine whether to use hyper-sparse CHUZC
  if (solve_phase == kSolvePhase1) {
    use_hyper_chuzc = false;
  } else {
    use_hyper_chuzc = false;  // true;
  }
  hyperChooseColumnClear();

  num_flip_since_rebuild = 0;
  // Data are fresh from rebuild
  status.has_fresh_rebuild = true;
  assert(solve_phase == kSolvePhase1 || solve_phase == kSolvePhase2);
}

void HEkkPrimal::iterate() {
  const HighsInt from_check_iter = 15;
  const HighsInt to_check_iter = from_check_iter + 10;
  if (ekk_instance_.debug_solve_report_) {
    ekk_instance_.debug_iteration_report_ =
        ekk_instance_.iteration_count_ >= from_check_iter &&
        ekk_instance_.iteration_count_ <= to_check_iter;
    if (ekk_instance_.debug_iteration_report_) {
      printf("HEkkPrimal::iterate Debug iteration %d\n",
             (int)ekk_instance_.iteration_count_);
    }
  }

  if (debugPrimalSimplex("Before iteration") ==
      HighsDebugStatus::kLogicalError) {
    solve_phase = kSolvePhaseError;
    return;
  }
  // Initialise row_out so that aborting iteration before CHUZR due to
  // numerical test of chosen reduced cost can be spotted - and
  // eliminates the unassigned read that can occur when the first
  // iteration is aborted in primal clean-up
  row_out = kNoRowSought;
  // Perform CHUZC
  //
  chuzc();
  if (variable_in == -1) {
    rebuild_reason = kRebuildReasonPossiblyOptimal;
    return;
  }

  // Perform FTRAN - and dual value cross-check to decide whether to use the
  // variable
  //
  // rebuild_reason = kRebuildReasonPossiblySingularBasis is set if
  // numerical trouble is detected
  if (!useVariableIn()) {
    if (rebuild_reason)
      assert(rebuild_reason == kRebuildReasonPossiblySingularBasis);
    return;
  }
  assert(!rebuild_reason);

  // Perform CHUZR
  if (solve_phase == kSolvePhase1) {
    phase1ChooseRow();
    assert(row_out != kNoRowSought);
    if (row_out == kNoRowChosen) {
      highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kError,
                  "Primal phase 1 choose row failed\n");
      solve_phase = kSolvePhaseError;
      return;
    }
  } else {
    chooseRow();
  }
  assert(!rebuild_reason);

  // Consider whether to perform a bound swap - either because it's
  // shorter than the pivoting step or, in the case of Phase 1,
  // because it's cheaper than pivoting - which may be questionable
  //
  // rebuild_reason = kRebuildReasonPossiblyPrimalUnbounded is set
  // in phase 2 if there's no pivot or bound swap. In phase 1 there is
  // always a pivot at this stage since row_out < 0 is trapped (above)
  // as an error.
  assert(solve_phase == kSolvePhase2 || row_out >= 0);
  considerBoundSwap();
  if (rebuild_reason == kRebuildReasonPossiblyPrimalUnbounded) return;
  assert(!rebuild_reason);

  if (row_out >= 0) {
    //
    // Perform unit BTRAN and PRICE to get pivotal row - and do a
    // numerical check.
    //
    // rebuild_reason = kRebuildReasonPossiblySingularBasis is set
    // if numerical trouble is detected
    assessPivot();
    if (rebuild_reason) {
      assert(rebuild_reason == kRebuildReasonPossiblySingularBasis);
      return;
    }
  }

  if (isBadBasisChange()) return;

  // Any pivoting is numerically acceptable, so perform update.
  //
  // rebuild_reason =
  // kRebuildReasonPrimalInfeasibleInPrimalSimplex is set if a
  // primal infeasibility is found in phase 2
  //
  // rebuild_reason = kRebuildReasonPossiblyPhase1Feasible is set in
  // phase 1 if the number of primal infeasibilities is reduced to
  // zero
  //
  // rebuild_reason = kRebuildReasonUpdateLimitReached is set in
  // either phase if the update count reaches the limit!
  //
  // rebuild_reason = kRebuildReasonSyntheticClockSaysInvert is
  // set in updateFactor() if it is considered to be more efficient to
  // reinvert.
  update();
  // Force rebuild if there are no infeasibilities in phase 1
  if (!ekk_instance_.info_.num_primal_infeasibilities &&
      solve_phase == kSolvePhase1)
    rebuild_reason = kRebuildReasonPossiblyPhase1Feasible;

  const bool ok_rebuild_reason =
      rebuild_reason == kRebuildReasonNo ||
      rebuild_reason == kRebuildReasonPossiblyPhase1Feasible ||
      rebuild_reason == kRebuildReasonPrimalInfeasibleInPrimalSimplex ||
      rebuild_reason == kRebuildReasonSyntheticClockSaysInvert ||
      rebuild_reason == kRebuildReasonUpdateLimitReached;
  if (!ok_rebuild_reason) {
    printf("HEkkPrimal::rebuild Solve %d; Iter %d: rebuild_reason = %d\n",
           (int)ekk_instance_.debug_solve_call_num_,
           (int)ekk_instance_.iteration_count_, (int)rebuild_reason);
    fflush(stdout);
  }
  assert(ok_rebuild_reason);
  assert(solve_phase == kSolvePhase1 || solve_phase == kSolvePhase2);
}

void HEkkPrimal::chuzc() {
  if (done_next_chuzc) assert(use_hyper_chuzc);
  vector<double>& workDual = ekk_instance_.info_.workDual_;
  ekk_instance_.applyTabooVariableIn(workDual, 0);
  if (use_hyper_chuzc) {
    // Perform hyper-sparse CHUZC and then check result using full CHUZC
    if (!done_next_chuzc) chooseColumn(true);
    const bool check_hyper_chuzc = true;
    if (check_hyper_chuzc) {
      HighsInt hyper_sparse_variable_in = variable_in;
      chooseColumn(false);
      double hyper_sparse_measure = 0;
      if (hyper_sparse_variable_in >= 0) {
        double squared_dual_infeasibility = workDual[hyper_sparse_variable_in] *
                                            workDual[hyper_sparse_variable_in];
        hyper_sparse_measure =
            squared_dual_infeasibility / edge_weight_[hyper_sparse_variable_in];
      }
      double measure = 0;
      if (variable_in >= 0) {
        double squared_dual_infeasibility =
            workDual[variable_in] * workDual[variable_in];
        measure = squared_dual_infeasibility / edge_weight_[variable_in];
      }
      double abs_measure_error = fabs(hyper_sparse_measure - measure);
      bool measure_error = abs_measure_error > 1e-12;
      // if (measure_error)
      //   printf("Iteration %" HIGHSINT_FORMAT
      //          ": Hyper-sparse CHUZC measure %g != %g = Full "
      //          "CHUZC measure (%" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT
      //          "): error %g\n",
      //          ekk_instance_.iteration_count_, hyper_sparse_measure, measure,
      //          hyper_sparse_variable_in, variable_in, abs_measure_error);

      // todo this fails on some rarer occasions, e.g. on glass4
      assert(!measure_error);
      variable_in = hyper_sparse_variable_in;
    }
  } else {
    chooseColumn(false);
  }
  ekk_instance_.unapplyTabooVariableIn(workDual);
}

void HEkkPrimal::chooseColumn(const bool hyper_sparse) {
  assert(!hyper_sparse || !done_next_chuzc);
  const vector<int8_t>& nonbasicMove = ekk_instance_.basis_.nonbasicMove_;
  const vector<double>& workDual = ekk_instance_.info_.workDual_;
  double best_measure = 0;
  variable_in = -1;

  const bool local_use_hyper_chuzc = hyper_sparse;
  // Consider nonbasic free columns first
  const HighsInt& num_nonbasic_free_col = nonbasic_free_col_set.count();
  if (local_use_hyper_chuzc) {
    if (!initialise_hyper_chuzc) hyperChooseColumn();
    if (initialise_hyper_chuzc) {
      analysis->simplexTimerStart(ChuzcHyperInitialiselClock);
      num_hyper_chuzc_candidates = 0;
      if (num_nonbasic_free_col) {
        const vector<HighsInt>& nonbasic_free_col_set_entry =
            nonbasic_free_col_set.entry();
        for (HighsInt ix = 0; ix < num_nonbasic_free_col; ix++) {
          HighsInt iCol = nonbasic_free_col_set_entry[ix];
          double dual_infeasibility = fabs(workDual[iCol]);
          if (dual_infeasibility > dual_feasibility_tolerance) {
            double measure =
                dual_infeasibility * dual_infeasibility / edge_weight_[iCol];
            addToDecreasingHeap(
                num_hyper_chuzc_candidates, max_num_hyper_chuzc_candidates,
                hyper_chuzc_measure, hyper_chuzc_candidate, measure, iCol);
          }
        }
      }
      // Now look at other columns
      for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
        double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
        if (dual_infeasibility > dual_feasibility_tolerance) {
          double measure =
              dual_infeasibility * dual_infeasibility / edge_weight_[iCol];
          addToDecreasingHeap(
              num_hyper_chuzc_candidates, max_num_hyper_chuzc_candidates,
              hyper_chuzc_measure, hyper_chuzc_candidate, measure, iCol);
        }
      }
      // Sort the heap
      sortDecreasingHeap(num_hyper_chuzc_candidates, hyper_chuzc_measure,
                         hyper_chuzc_candidate);
      initialise_hyper_chuzc = false;
      analysis->simplexTimerStop(ChuzcHyperInitialiselClock);
      // Choose the first entry - if there is one
      if (num_hyper_chuzc_candidates) {
        variable_in = hyper_chuzc_candidate[1];
        best_measure = hyper_chuzc_measure[1];
        max_hyper_chuzc_non_candidate_measure =
            hyper_chuzc_measure[num_hyper_chuzc_candidates];
        if (report_hyper_chuzc)
          printf(
              "Full CHUZC: Max         measure is %9.4g for column "
              "%4" HIGHSINT_FORMAT
              ", and "
              "max non-candidate measure of  %9.4g\n",
              best_measure, variable_in, max_hyper_chuzc_non_candidate_measure);
      }
    }
  } else {
    analysis->simplexTimerStart(ChuzcPrimalClock);
    // Choose any attractive nonbasic free column
    if (num_nonbasic_free_col) {
      const vector<HighsInt>& nonbasic_free_col_set_entry =
          nonbasic_free_col_set.entry();
      for (HighsInt ix = 0; ix < num_nonbasic_free_col; ix++) {
        HighsInt iCol = nonbasic_free_col_set_entry[ix];
        double dual_infeasibility = fabs(workDual[iCol]);
        if (dual_infeasibility > dual_feasibility_tolerance &&
            dual_infeasibility * dual_infeasibility >
                best_measure * edge_weight_[iCol]) {
          variable_in = iCol;
          best_measure =
              dual_infeasibility * dual_infeasibility / edge_weight_[iCol];
        }
      }
    }
    // Now look at other columns
    for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
      double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
      if (dual_infeasibility > dual_feasibility_tolerance &&
          dual_infeasibility * dual_infeasibility >
              best_measure * edge_weight_[iCol]) {
        variable_in = iCol;
        best_measure =
            dual_infeasibility * dual_infeasibility / edge_weight_[iCol];
      }
    }
    analysis->simplexTimerStop(ChuzcPrimalClock);
  }
  //  printf("ChooseColumn: Iteration %" HIGHSINT_FORMAT ", choose column %"
  //  HIGHSINT_FORMAT " with measure %g\n",
  //	 ekk_instance_.iteration_count_, variable_in, best_measure);
}

bool HEkkPrimal::useVariableIn() {
  // rebuild_reason = kRebuildReasonPossiblySingularBasis is set if
  // numerical trouble is detected
  HighsSimplexInfo& info = ekk_instance_.info_;
  vector<double>& workDual = info.workDual_;
  const vector<int8_t>& nonbasicMove = ekk_instance_.basis_.nonbasicMove_;
  const double updated_theta_dual = workDual[variable_in];
  // Determine the move direction - can't use nonbasicMove_[variable_in]
  // due to free columns
  move_in = updated_theta_dual > 0 ? -1 : 1;
  // Unless the variable is free, nonbasicMove[variable_in] should be the same
  // as move_in
  if (nonbasicMove[variable_in]) assert(nonbasicMove[variable_in] == move_in);
  //
  // FTRAN
  //
  // Compute pivot column
  ekk_instance_.pivotColumnFtran(variable_in, col_aq);
  // Compute the dual for the pivot column and compare it with the
  // updated value
  double computed_theta_dual =
      ekk_instance_.computeDualForTableauColumn(variable_in, col_aq);
  ekk_instance_.debugUpdatedDual(updated_theta_dual, computed_theta_dual);
  // Feed in the computed dual value.
  //
  // The sum of dual infeasibilities (and maybe max dual
  // infeasibility) will be wrong, but there's a big tolerance on
  // this in debugSimplex. Have to be careful (below) if the computed
  // dual value is no longer a dual infeasibility
  info.workDual_[variable_in] = computed_theta_dual;
  // Reassign theta_dual to be the computed value
  theta_dual = info.workDual_[variable_in];
  // Determine whether theta_dual is too small or has changed sign
  const bool theta_dual_small = fabs(theta_dual) <= dual_feasibility_tolerance;
  const bool theta_dual_sign_error =
      updated_theta_dual * computed_theta_dual <= 0;

  // If theta_dual is small, then it's no longer a dual infeasibility,
  // so reduce the number of dual infeasibilities. Otherwise an error
  // is identified in debugSimplex
  if (theta_dual_small) ekk_instance_.info_.num_dual_infeasibilities--;
  if (theta_dual_small || theta_dual_sign_error) {
    // The computed dual is small or has a sign error, so don't use it
    std::string theta_dual_size = "";
    if (theta_dual_small) theta_dual_size = "; too small";
    std::string theta_dual_sign = "";
    if (theta_dual_sign_error) theta_dual_sign = "; sign error";
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "Chosen entering variable %" HIGHSINT_FORMAT
                " (Iter = %" HIGHSINT_FORMAT "; Update = %" HIGHSINT_FORMAT
                ") has computed "
                "(updated) dual of %10.4g (%10.4g) so don't use it%s%s\n",
                variable_in, ekk_instance_.iteration_count_, info.update_count,
                computed_theta_dual, updated_theta_dual,
                theta_dual_size.c_str(), theta_dual_sign.c_str());
    // If a significant computed dual has sign error, consider reinverting
    if (!theta_dual_small && info.update_count > 0)
      rebuild_reason = kRebuildReasonPossiblySingularBasis;
    hyperChooseColumnClear();
    return false;
  }
  return true;
}

void HEkkPrimal::phase1ChooseRow() {
  const HighsSimplexInfo& info = ekk_instance_.info_;
  const vector<double>& baseLower = info.baseLower_;
  const vector<double>& baseUpper = info.baseUpper_;
  const vector<double>& baseValue = info.baseValue_;
  analysis->simplexTimerStart(Chuzr1Clock);
  // Collect phase 1 theta lists
  //

  const double dPivotTol = info.update_count < 10   ? 1e-9
                           : info.update_count < 20 ? 1e-8
                                                    : 1e-7;
  ph1SorterR.clear();
  ph1SorterT.clear();
  for (HighsInt i = 0; i < col_aq.count; i++) {
    HighsInt iRow = col_aq.index[i];
    double dAlpha = col_aq.array[iRow] * move_in;

    // When the basic variable x[i] decrease
    if (dAlpha > +dPivotTol) {
      // Whether it can become feasible by going below its upper bound
      if (baseValue[iRow] > baseUpper[iRow] + primal_feasibility_tolerance) {
        double dFeasTheta =
            (baseValue[iRow] - baseUpper[iRow] - primal_feasibility_tolerance) /
            dAlpha;
        ph1SorterR.push_back(std::make_pair(dFeasTheta, iRow));
        ph1SorterT.push_back(std::make_pair(dFeasTheta, iRow));
      }
      // Whether it can become infeasible (again) by going below its
      // lower bound
      if (baseValue[iRow] > baseLower[iRow] - primal_feasibility_tolerance &&
          baseLower[iRow] > -kHighsInf) {
        double dRelaxTheta =
            (baseValue[iRow] - baseLower[iRow] + primal_feasibility_tolerance) /
            dAlpha;
        double dTightTheta = (baseValue[iRow] - baseLower[iRow]) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dRelaxTheta, iRow - num_row));
        ph1SorterT.push_back(std::make_pair(dTightTheta, iRow - num_row));
      }
    }

    // When the basic variable x[i] increase
    if (dAlpha < -dPivotTol) {
      // Whether it can become feasible by going above its lower bound
      if (baseValue[iRow] < baseLower[iRow] - primal_feasibility_tolerance) {
        double dFeasTheta =
            (baseValue[iRow] - baseLower[iRow] + primal_feasibility_tolerance) /
            dAlpha;
        ph1SorterR.push_back(std::make_pair(dFeasTheta, iRow - num_row));
        ph1SorterT.push_back(std::make_pair(dFeasTheta, iRow - num_row));
      }
      // Whether it can become infeasible (again) by going above its
      // upper bound
      if (baseValue[iRow] < baseUpper[iRow] + primal_feasibility_tolerance &&
          baseUpper[iRow] < +kHighsInf) {
        double dRelaxTheta =
            (baseValue[iRow] - baseUpper[iRow] - primal_feasibility_tolerance) /
            dAlpha;
        double dTightTheta = (baseValue[iRow] - baseUpper[iRow]) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dRelaxTheta, iRow));
        ph1SorterT.push_back(std::make_pair(dTightTheta, iRow));
      }
    }
  }

  analysis->simplexTimerStop(Chuzr1Clock);
  // When there are no candidates at all, we can leave it here
  if (ph1SorterR.empty()) {
    row_out = kNoRowChosen;
    variable_out = -1;
    return;
  }

  // Now sort the relaxed theta to find the final break point. TODO:
  // Consider partial sort. Or heapify [O(n)] and then pop k points
  // [kO(log(n))].

  analysis->simplexTimerStart(Chuzr2Clock);
  pdqsort(ph1SorterR.begin(), ph1SorterR.end());
  double dMaxTheta = ph1SorterR[0].first;
  double dGradient = fabs(theta_dual);
  for (const auto& candidate : ph1SorterR) {
    double dMyTheta = candidate.first;
    HighsInt index = candidate.second;
    HighsInt iRow = index >= 0 ? index : index + num_row;
    dGradient -= fabs(col_aq.array[iRow]);
    // Stop when the gradient start to decrease
    if (dGradient <= 0) {
      break;
    }
    dMaxTheta = dMyTheta;
  }

  // Find out the biggest possible alpha for pivot
  pdqsort(ph1SorterT.begin(), ph1SorterT.end());
  double dMaxAlpha = 0.0;
  size_t iLast = ph1SorterT.size();
  for (size_t i = 0; i < ph1SorterT.size(); i++) {
    double dMyTheta = ph1SorterT[i].first;
    HighsInt index = ph1SorterT[i].second;
    HighsInt iRow = index >= 0 ? index : index + num_row;
    double dAbsAlpha = fabs(col_aq.array[iRow]);
    // Stop when the theta is too large
    if (dMyTheta > dMaxTheta) {
      iLast = i;
      break;
    }
    // Update the maximal possible alpha
    if (dMaxAlpha < dAbsAlpha) {
      dMaxAlpha = dAbsAlpha;
    }
  }

  // Finally choose a pivot with good enough alpha, working backwards
  row_out = kNoRowChosen;
  variable_out = -1;
  move_out = 0;
  for (size_t i = iLast; i > 0; i--) {
    HighsInt index = ph1SorterT[i - 1].second;
    HighsInt iRow = index >= 0 ? index : index + num_row;
    double dAbsAlpha = fabs(col_aq.array[iRow]);
    if (dAbsAlpha > dMaxAlpha * 0.1) {
      row_out = iRow;
      move_out = index >= 0 ? 1 : -1;
      break;
    }
  }
  analysis->simplexTimerStop(Chuzr2Clock);
}

void HEkkPrimal::chooseRow() {
  HighsSimplexInfo& info = ekk_instance_.info_;
  const vector<double>& baseLower = info.baseLower_;
  const vector<double>& baseUpper = info.baseUpper_;
  const vector<double>& baseValue = info.baseValue_;
  analysis->simplexTimerStart(Chuzr1Clock);
  // Initialize
  row_out = kNoRowChosen;

  // Choose row pass 1
  double alphaTol = info.update_count < 10   ? 1e-9
                    : info.update_count < 20 ? 1e-8
                                             : 1e-7;

  double relaxTheta = 1e100;
  double relaxSpace;
  for (HighsInt i = 0; i < col_aq.count; i++) {
    HighsInt iRow = col_aq.index[i];
    double alpha = col_aq.array[iRow] * move_in;
    if (alpha > alphaTol) {
      relaxSpace =
          baseValue[iRow] - baseLower[iRow] + primal_feasibility_tolerance;
      if (relaxSpace < relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    } else if (alpha < -alphaTol) {
      relaxSpace =
          baseValue[iRow] - baseUpper[iRow] - primal_feasibility_tolerance;
      if (relaxSpace > relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    }
  }
  analysis->simplexTimerStop(Chuzr1Clock);

  analysis->simplexTimerStart(Chuzr2Clock);
  double bestAlpha = 0;
  for (HighsInt i = 0; i < col_aq.count; i++) {
    HighsInt iRow = col_aq.index[i];
    double alpha = col_aq.array[iRow] * move_in;
    if (alpha > alphaTol) {
      // Positive pivotal column entry
      double tightSpace = baseValue[iRow] - baseLower[iRow];
      if (tightSpace < relaxTheta * alpha) {
        if (bestAlpha < alpha) {
          bestAlpha = alpha;
          row_out = iRow;
        }
      }
    } else if (alpha < -alphaTol) {
      // Negative pivotal column entry
      double tightSpace = baseValue[iRow] - baseUpper[iRow];
      if (tightSpace > relaxTheta * alpha) {
        if (bestAlpha < -alpha) {
          bestAlpha = -alpha;
          row_out = iRow;
        }
      }
    }
  }
  analysis->simplexTimerStop(Chuzr2Clock);
}

void HEkkPrimal::considerBoundSwap() {
  const HighsSimplexInfo& info = ekk_instance_.info_;
  const vector<double>& workLower = info.workLower_;
  const vector<double>& workUpper = info.workUpper_;
  const vector<double>& baseLower = info.baseLower_;
  const vector<double>& baseUpper = info.baseUpper_;
  const vector<double>& workValue = info.workValue_;
  const vector<double>& baseValue = info.baseValue_;

  // Compute the primal theta and see if we should have done a bound
  // flip instead
  if (row_out == kNoRowChosen) {
    assert(solve_phase == kSolvePhase2);
    // No binding ratio in CHUZR, so flip or unbounded
    theta_primal = move_in * kHighsInf;
    move_out = 0;
  } else {
    assert(row_out >= 0);
    // Determine the step to the leaving bound
    //
    alpha_col = col_aq.array[row_out];
    // In Phase 1, move_out depends on whether the leaving variable is
    // becoming feasible - moves up to lower (down to upper) - or
    // remaining feasible - moves down to lower (up to upper) - so
    // can't be set so easily as in phase 2
    if (solve_phase == kSolvePhase2)
      move_out = alpha_col * move_in > 0 ? -1 : 1;
    theta_primal = 0;
    if (move_out == 1) {
      theta_primal = (baseValue[row_out] - baseUpper[row_out]) / alpha_col;
    } else {
      theta_primal = (baseValue[row_out] - baseLower[row_out]) / alpha_col;
    }
    assert(theta_primal > -kHighsInf && theta_primal < kHighsInf);
  }

  // Look to see if there is a bound flip
  bool flipped = false;
  double lower_in = workLower[variable_in];
  double upper_in = workUpper[variable_in];
  value_in = workValue[variable_in] + theta_primal;
  if (move_in > 0) {
    if (value_in > upper_in + primal_feasibility_tolerance) {
      flipped = true;
      row_out = kNoRowChosen;
      value_in = upper_in;
      theta_primal = upper_in - lower_in;
    }
  } else {
    if (value_in < lower_in - primal_feasibility_tolerance) {
      flipped = true;
      row_out = kNoRowChosen;
      value_in = lower_in;
      theta_primal = lower_in - upper_in;
    }
  }
  const bool pivot_or_flipped = row_out >= 0 || flipped;
  if (solve_phase == kSolvePhase2) {
    // Check for possible unboundedness
    if (!pivot_or_flipped) {
      rebuild_reason = kRebuildReasonPossiblyPrimalUnbounded;
      return;
    }
  }
  // Check for possible error
  assert(pivot_or_flipped);
  assert(flipped == (row_out == kNoRowChosen));
}

void HEkkPrimal::assessPivot() {
  assert(row_out >= 0);
  // Record the pivot entry
  alpha_col = col_aq.array[row_out];
  variable_out = ekk_instance_.basis_.basicIndex_[row_out];

  // Compute the tableau row
  //
  // BTRAN
  //
  // Compute unit BTran for tableau row and FT update
  ekk_instance_.unitBtran(row_out, row_ep);
  //
  // PRICE
  //
  const bool quad_precision = false;
  ekk_instance_.tableauRowPrice(quad_precision, row_ep, row_ap);

  // Checks row-wise pivot against column-wise pivot for
  // numerical trouble
  //
  // rebuild_reason = kRebuildReasonPossiblySingularBasis is set if
  // numerical trouble is detected
  updateVerify();
}

void HEkkPrimal::update() {
  // Perform update operations that are independent of phase
  HighsSimplexInfo& info = ekk_instance_.info_;
  assert(!rebuild_reason);
  bool flipped = row_out < 0;
  if (flipped) {
    variable_out = variable_in;
    alpha_col = 0;
    numericalTrouble = 0;
    info.workValue_[variable_in] = value_in;
    assert(ekk_instance_.basis_.nonbasicMove_[variable_in] == move_in);
    ekk_instance_.basis_.nonbasicMove_[variable_in] = -move_in;
  } else {
    // Adjust perturbation if leaving equation
    adjustPerturbedEquationOut();
  }

  // Start hyper-sparse CHUZC, that takes place through phase1Update()
  hyperChooseColumnStart();

  if (solve_phase == kSolvePhase1) {
    // Update primal values
    phase1UpdatePrimal();

    // Update the duals with respect to feasibility changes
    basicFeasibilityChangeUpdateDual();

    // For hyper-sparse CHUZC, analyse the duals that have just changed
    hyperChooseColumnBasicFeasibilityChange();

  } else {
    // Update primal values, and identify any infeasibilities
    //
    // rebuild_reason =
    // kRebuildReasonPrimalInfeasibleInPrimalSimplex is set if a
    // primal infeasibility is found
    phase2UpdatePrimal();
  }

  assert(rebuild_reason == kRebuildReasonNo ||
         rebuild_reason == kRebuildReasonPrimalInfeasibleInPrimalSimplex);

  if (flipped) {
    info.primal_bound_swap++;
    ekk_instance_.invalidateDualInfeasibilityRecord();
    iterationAnalysis();
    localReportIter();
    num_flip_since_rebuild++;
    // Update the synthetic clock for UPDATE
    ekk_instance_.total_synthetic_tick_ += col_aq.synthetic_tick;
    return;
  }

  assert(row_out >= 0);

  // Now set the value of the entering variable
  info.baseValue_[row_out] = value_in;
  // Consider whether the entering value is feasible and, if not, take
  // action
  //
  // rebuild_reason =
  // kRebuildReasonPrimalInfeasibleInPrimalSimplex is set in
  // phase 2 if a primal infeasibility is found
  considerInfeasibleValueIn();

  // Update the dual values
  theta_dual = info.workDual_[variable_in];
  updateDual();

  // Update any non-unit primal edge weights
  if (edge_weight_mode == EdgeWeightMode::kDevex) {
    updateDevex();
  } else if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
    debugPrimalSteepestEdgeWeights("before update");
    updatePrimalSteepestEdgeWeights();
  }

  // If entering column was nonbasic free, remove it from the set
  removeNonbasicFreeColumn();

  // For hyper-sparse CHUZC, analyse the duals and weights that have
  // just changed
  hyperChooseColumnDualChange();

  if (ekk_instance_.status_.has_dual_steepest_edge_weights) {
    ekk_instance_.devDebugDualSteepestEdgeWeights("before update");
    updateDualSteepestEdgeWeights();
  }
  // Perform pivoting
  //
  // Transform the vectors used in updateFactor if the simplex NLA involves
  // scaling
  ekk_instance_.transformForUpdate(&col_aq, &row_ep, variable_in, &row_out);
  //
  // Update the sets of indices of basic and nonbasic variables
  ekk_instance_.updatePivots(variable_in, row_out, move_out);
  //
  // Update the invertible representation of the basis matrix
  ekk_instance_.updateFactor(&col_aq, &row_ep, &row_out, &rebuild_reason);

  if (ekk_instance_.status_.has_dual_steepest_edge_weights)
    ekk_instance_.devDebugDualSteepestEdgeWeights("after  update");
  if (edge_weight_mode == EdgeWeightMode::kSteepestEdge)
    debugPrimalSteepestEdgeWeights("after update");
  //
  // Update the row-wise representation of the nonbasic columns
  ekk_instance_.updateMatrix(variable_in, variable_out);
  if (info.update_count >= info.update_limit)
    rebuild_reason = kRebuildReasonUpdateLimitReached;

  // Update the iteration count
  ekk_instance_.iteration_count_++;

  // Reset the devex when there are too many errors
  if (edge_weight_mode == EdgeWeightMode::kDevex &&
      num_bad_devex_weight_ > kAllowedNumBadDevexWeight)
    initialiseDevexFramework();

  // Report on the iteration
  iterationAnalysis();
  localReportIter();

  // Update the synthetic clock for UPDATE
  ekk_instance_.total_synthetic_tick_ += col_aq.synthetic_tick;
  ekk_instance_.total_synthetic_tick_ += row_ep.synthetic_tick;

  // Perform hyper-sparse CHUZC
  hyperChooseColumn();
}

void HEkkPrimal::hyperChooseColumn() {
  if (!use_hyper_chuzc) return;
  if (initialise_hyper_chuzc) return;
  analysis->simplexTimerStart(ChuzcHyperClock);
  const vector<int8_t>& nonbasicMove = ekk_instance_.basis_.nonbasicMove_;
  const vector<int8_t>& nonbasicFlag = ekk_instance_.basis_.nonbasicFlag_;
  const vector<double>& workDual = ekk_instance_.info_.workDual_;
  if (report_hyper_chuzc)
    printf(
        "H-S  CHUZC: Max changed measure is %9.4g for column %4" HIGHSINT_FORMAT
        "",
        max_changed_measure_value, max_changed_measure_column);
  double best_measure = max_changed_measure_value;
  variable_in = -1;
  if (max_changed_measure_column >= 0) {
    // Use max_changed_measure_column if it is well defined and has
    // nonzero dual. It may have been zeroed because it is taboo
    if (workDual[max_changed_measure_column])
      variable_in = max_changed_measure_column;
  }
  const bool consider_nonbasic_free_column =
      (nonbasic_free_col_set.count() != 0);
  if (num_hyper_chuzc_candidates) {
    for (HighsInt iEntry = 1; iEntry <= num_hyper_chuzc_candidates; iEntry++) {
      HighsInt iCol = hyper_chuzc_candidate[iEntry];
      if (nonbasicFlag[iCol] == kNonbasicFlagFalse) {
        assert(!nonbasicMove[iCol]);
        continue;
      }
      // Assess any dual infeasibility
      double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
      if (consider_nonbasic_free_column) {
        if (nonbasic_free_col_set.in(iCol))
          dual_infeasibility = fabs(workDual[iCol]);
      }
      if (dual_infeasibility > dual_feasibility_tolerance) {
        if (dual_infeasibility * dual_infeasibility >
            best_measure * edge_weight_[iCol]) {
          best_measure =
              dual_infeasibility * dual_infeasibility / edge_weight_[iCol];
          variable_in = iCol;
        }
      }
    }
  }
  if (variable_in != max_changed_measure_column) {
    if (report_hyper_chuzc)
      printf(
          ", and after HS CHUZC set it is now %9.4g for column "
          "%4" HIGHSINT_FORMAT "",
          best_measure, variable_in);
    max_hyper_chuzc_non_candidate_measure =
        max(max_changed_measure_value, max_hyper_chuzc_non_candidate_measure);
  }
  if (best_measure >= max_hyper_chuzc_non_candidate_measure) {
    // Candidate is at least as good as any unknown column, so accept it
    done_next_chuzc = true;
    if (report_hyper_chuzc)
      printf(", and no       has  measure >  %9.4g\n",
             max_hyper_chuzc_non_candidate_measure);
  } else {
    // Candidate isn't as good as best unknown column, so do a full CHUZC
    // Shouldn't claim to have done the next CHUZC
    assert(!done_next_chuzc);
    done_next_chuzc = false;
    initialise_hyper_chuzc = true;
    if (report_hyper_chuzc)
      printf(", but some may have measure >= %9.4g\n",
             max_hyper_chuzc_non_candidate_measure);
  }
  analysis->simplexTimerStop(ChuzcHyperClock);
}

void HEkkPrimal::hyperChooseColumnStart() {
  max_changed_measure_value = 0;
  max_changed_measure_column = -1;
  done_next_chuzc = false;
}

void HEkkPrimal::hyperChooseColumnClear() {
  initialise_hyper_chuzc = use_hyper_chuzc;
  max_hyper_chuzc_non_candidate_measure = -1;
  done_next_chuzc = false;
}

void HEkkPrimal::hyperChooseColumnChangedInfeasibility(
    const double infeasibility, const HighsInt iCol) {
  if (infeasibility * infeasibility >
      max_changed_measure_value * edge_weight_[iCol]) {
    max_hyper_chuzc_non_candidate_measure =
        max(max_changed_measure_value, max_hyper_chuzc_non_candidate_measure);
    max_changed_measure_value =
        infeasibility * infeasibility / edge_weight_[iCol];
    max_changed_measure_column = iCol;
  } else if (infeasibility * infeasibility >
             max_hyper_chuzc_non_candidate_measure * edge_weight_[iCol]) {
    max_hyper_chuzc_non_candidate_measure =
        infeasibility * infeasibility / edge_weight_[iCol];
  }
}

void HEkkPrimal::hyperChooseColumnBasicFeasibilityChange() {
  if (!use_hyper_chuzc) return;
  analysis->simplexTimerStart(ChuzcHyperBasicFeasibilityChangeClock);
  const vector<double>& workDual = ekk_instance_.info_.workDual_;
  const vector<int8_t>& nonbasicMove = ekk_instance_.basis_.nonbasicMove_;
  HighsInt to_entry;
  const bool use_row_indices = ekk_instance_.simplex_nla_.sparseLoopStyle(
      row_basic_feasibility_change.count, num_col, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iCol =
        use_row_indices ? row_basic_feasibility_change.index[iEntry] : iEntry;
    double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
    if (dual_infeasibility > dual_feasibility_tolerance)
      hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  const bool use_col_indices = ekk_instance_.simplex_nla_.sparseLoopStyle(
      col_basic_feasibility_change.count, num_row, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iRow =
        use_col_indices ? col_basic_feasibility_change.index[iEntry] : iEntry;
    HighsInt iCol = num_col + iRow;
    double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
    if (dual_infeasibility > dual_feasibility_tolerance)
      hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  // Any nonbasic free columns will be handled explicitly in
  // hyperChooseColumnDualChange, so only look at them here if not
  // flipping
  const HighsInt& num_nonbasic_free_col = nonbasic_free_col_set.count();
  if (row_out < 0 && num_nonbasic_free_col) {
    const vector<HighsInt>& nonbasic_free_col_set_entry =
        nonbasic_free_col_set.entry();
    for (HighsInt iEntry = 0; iEntry < num_nonbasic_free_col; iEntry++) {
      HighsInt iCol = nonbasic_free_col_set_entry[iEntry];
      double dual_infeasibility = fabs(workDual[iCol]);
      if (dual_infeasibility > dual_feasibility_tolerance)
        hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
    }
  }
  analysis->simplexTimerStop(ChuzcHyperBasicFeasibilityChangeClock);
}

void HEkkPrimal::hyperChooseColumnDualChange() {
  if (!use_hyper_chuzc) return;
  analysis->simplexTimerStart(ChuzcHyperDualClock);
  const vector<double>& workDual = ekk_instance_.info_.workDual_;
  const vector<int8_t>& nonbasicMove = ekk_instance_.basis_.nonbasicMove_;
  HighsInt to_entry;
  // Look at changes in the columns and assess any dual infeasibility
  const bool use_row_indices = ekk_instance_.simplex_nla_.sparseLoopStyle(
      row_ap.count, num_col, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iCol = use_row_indices ? row_ap.index[iEntry] : iEntry;
    double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
    if (iCol == check_column && ekk_instance_.iteration_count_ >= check_iter) {
      double measure =
          dual_infeasibility * dual_infeasibility / edge_weight_[iCol];
      if (report_hyper_chuzc) {
        printf("Changing column %" HIGHSINT_FORMAT ": measure = %g \n",
               check_column, measure);
      }
    }
    if (dual_infeasibility > dual_feasibility_tolerance)
      hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  // Look at changes in the rows and assess any dual infeasibility
  const bool use_col_indices = ekk_instance_.simplex_nla_.sparseLoopStyle(
      row_ep.count, num_row, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iRow = use_col_indices ? row_ep.index[iEntry] : iEntry;
    HighsInt iCol = iRow + num_col;
    double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
    if (iCol == check_column && ekk_instance_.iteration_count_ >= check_iter) {
      double measure =
          dual_infeasibility * dual_infeasibility / edge_weight_[iCol];
      if (report_hyper_chuzc) {
        printf("Changing column %" HIGHSINT_FORMAT ": measure = %g \n",
               check_column, measure);
      }
    }
    if (dual_infeasibility > dual_feasibility_tolerance)
      hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  // Look for measure changes in any nonbasic free columns and assess
  // any dual infeasibility
  const HighsInt& num_nonbasic_free_col = nonbasic_free_col_set.count();
  if (num_nonbasic_free_col) {
    const vector<HighsInt>& nonbasic_free_col_set_entry =
        nonbasic_free_col_set.entry();
    for (HighsInt iEntry = 0; iEntry < num_nonbasic_free_col; iEntry++) {
      HighsInt iCol = nonbasic_free_col_set_entry[iEntry];
      double dual_infeasibility = fabs(workDual[iCol]);
      if (dual_infeasibility > dual_feasibility_tolerance)
        hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
    }
  }
  // Assess any dual infeasibility for the leaving column - should be dual
  // feasible!
  HighsInt iCol = variable_out;
  double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
  if (dual_infeasibility > dual_feasibility_tolerance) {
    printf("Dual infeasibility %g for leaving column!\n", dual_infeasibility);
    assert(dual_infeasibility <= dual_feasibility_tolerance);
    hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  analysis->simplexTimerStop(ChuzcHyperDualClock);
}

void HEkkPrimal::updateDual() {
  analysis->simplexTimerStart(UpdateDualClock);
  assert(alpha_col);
  assert(row_out >= 0);
  vector<double>& workDual = ekk_instance_.info_.workDual_;
  //  const vector<HighsInt>& nonbasicMove =
  //  ekk_instance_.basis_.nonbasicMove_;
  // Update the duals
  theta_dual = workDual[variable_in] / alpha_col;
  for (HighsInt iEl = 0; iEl < row_ap.count; iEl++) {
    HighsInt iCol = row_ap.index[iEl];
    workDual[iCol] -= theta_dual * row_ap.array[iCol];
  }
  for (HighsInt iEl = 0; iEl < row_ep.count; iEl++) {
    HighsInt iRow = row_ep.index[iEl];
    HighsInt iCol = iRow + num_col;
    workDual[iCol] -= theta_dual * row_ep.array[iRow];
  }
  // Dual for the pivot
  workDual[variable_in] = 0;
  workDual[variable_out] = -theta_dual;

  ekk_instance_.invalidateDualInfeasibilityRecord();
  // After dual update in primal simplex the dual objective value is not known
  ekk_instance_.status_.has_dual_objective_value = false;
  analysis->simplexTimerStop(UpdateDualClock);
}

void HEkkPrimal::phase1ComputeDual() {
  HighsSimplexInfo& info = ekk_instance_.info_;
  const vector<int8_t>& nonbasicFlag = ekk_instance_.basis_.nonbasicFlag_;

  HVector buffer;
  buffer.setup(num_row);
  buffer.clear();
  buffer.count = 0;
  // Accumulate costs for checking
  info.workCost_.assign(num_tot, 0);
  // Zero the dual values
  info.workDual_.assign(num_tot, 0);
  // Determine the base value for cost perturbation
  const double base =
      info.primal_simplex_phase1_cost_perturbation_multiplier * 5e-7;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    const double value = info.baseValue_[iRow];
    const double lower = info.baseLower_[iRow];
    const double upper = info.baseUpper_[iRow];
    HighsInt bound_violated = 0;
    if (value < lower - primal_feasibility_tolerance) {
      bound_violated = -1;
    } else if (value > upper + primal_feasibility_tolerance) {
      bound_violated = 1;
    }
    if (!bound_violated) continue;
    double cost = bound_violated;
    if (base) cost *= 1 + base * info.numTotRandomValue_[iRow];
    buffer.array[iRow] = cost;
    buffer.index[buffer.count++] = iRow;
  }
  if (buffer.count <= 0) {
    // Strange, should be a non-trivial RHS
    assert(buffer.count > 0);
    return;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    info.workCost_[ekk_instance_.basis_.basicIndex_[iRow]] = buffer.array[iRow];
  //
  // Full BTRAN
  //
  ekk_instance_.fullBtran(buffer);
  //
  // Full PRICE
  //
  HVector bufferLong;
  bufferLong.setup(num_col);
  ekk_instance_.fullPrice(buffer, bufferLong);

  for (HighsInt iCol = 0; iCol < num_col; iCol++)
    info.workDual_[iCol] = -nonbasicFlag[iCol] * bufferLong.array[iCol];
  for (HighsInt iRow = 0, iCol = num_col; iRow < num_row; iRow++, iCol++)
    info.workDual_[iCol] = -nonbasicFlag[iCol] * buffer.array[iRow];
}

void HEkkPrimal::phase1UpdatePrimal() {
  analysis->simplexTimerStart(UpdatePrimalClock);
  HighsSimplexInfo& info = ekk_instance_.info_;
  col_basic_feasibility_change.clear();
  //
  // Update basic primal values, identifying all the feasibility
  // changes giving a value to col_basic_feasibility_change so that the duals
  // can be updated.
  //
  // Determine the base value for cost perturbation
  const double base =
      info.primal_simplex_phase1_cost_perturbation_multiplier * 5e-7;
  //  if (ekk_instance_.sparseLoopStyle(col_aq.count, num_row, to_entry)) {
  for (HighsInt iEl = 0; iEl < col_aq.count; iEl++) {
    HighsInt iRow = col_aq.index[iEl];
    info.baseValue_[iRow] -= theta_primal * col_aq.array[iRow];
    HighsInt iCol = ekk_instance_.basis_.basicIndex_[iRow];
    double was_cost = info.workCost_[iCol];
    const double value = info.baseValue_[iRow];
    const double lower = info.baseLower_[iRow];
    const double upper = info.baseUpper_[iRow];
    HighsInt bound_violated = 0;
    if (value < lower - primal_feasibility_tolerance) {
      bound_violated = -1.0;
    } else if (value > upper + primal_feasibility_tolerance) {
      bound_violated = 1.0;
    }
    double cost = bound_violated;
    if (base) cost *= 1 + base * info.numTotRandomValue_[iRow];
    info.workCost_[iCol] = cost;
    if (was_cost) {
      if (!cost) info.num_primal_infeasibilities--;
    } else {
      if (cost) info.num_primal_infeasibilities++;
    }
    double delta_cost = cost - was_cost;
    if (delta_cost) {
      col_basic_feasibility_change.array[iRow] = delta_cost;
      col_basic_feasibility_change.index[col_basic_feasibility_change.count++] =
          iRow;
      if (iCol >= num_col) info.workDual_[iCol] += delta_cost;
    }
  }
  // Don't set baseValue[row_out] yet so that dual update due to
  // feasibility changes is done correctly
  ekk_instance_.invalidatePrimalMaxSumInfeasibilityRecord();
  analysis->simplexTimerStop(UpdatePrimalClock);
}

void HEkkPrimal::considerInfeasibleValueIn() {
  assert(row_out >= 0);
  HighsSimplexInfo& info = ekk_instance_.info_;
  // Determine the base value for cost perturbation
  const double base =
      info.primal_simplex_phase1_cost_perturbation_multiplier * 5e-7;
  const double lower = info.workLower_[variable_in];
  const double upper = info.workUpper_[variable_in];
  HighsInt bound_violated = 0;
  if (value_in < lower - primal_feasibility_tolerance) {
    bound_violated = -1;
  } else if (value_in > upper + primal_feasibility_tolerance) {
    bound_violated = 1;
  }
  if (!bound_violated) return;
  // The primal value of the entering variable is not feasible
  if (solve_phase == kSolvePhase1) {
    info.num_primal_infeasibilities++;
    double cost = bound_violated;
    if (base) cost *= 1 + base * info.numTotRandomValue_[row_out];
    info.workCost_[variable_in] = cost;
    info.workDual_[variable_in] += cost;
  } else if (primal_correction_strategy ==
             kSimplexPrimalCorrectionStrategyNone) {
    // @primal_infeasibility calculation
    double primal_infeasibility;
    if (bound_violated < 0) {
      primal_infeasibility = lower - value_in;
    } else {
      primal_infeasibility = value_in - upper;
    }
    info.num_primal_infeasibilities++;
    highsLogDev(
        ekk_instance_.options_->log_options, HighsLogType::kWarning,
        "Entering variable has primal infeasibility of %g for [%g, %g, %g]\n",
        primal_infeasibility, lower, value_in, upper);
    rebuild_reason = kRebuildReasonPrimalInfeasibleInPrimalSimplex;
  } else {
    double bound_shift;
    if (bound_violated > 0) {
      // Perturb the upper bound to accommodate the infeasibility
      shiftBound(false, variable_in, value_in,
                 info.numTotRandomValue_[variable_in],
                 info.workUpper_[variable_in], bound_shift);
      info.workUpperShift_[variable_in] += bound_shift;
    } else {
      // Perturb the lower bound to accommodate the infeasibility
      shiftBound(true, variable_in, value_in,
                 info.numTotRandomValue_[variable_in],
                 info.workLower_[variable_in], bound_shift);
      info.workLowerShift_[variable_in] += bound_shift;
    }
    // Surely better to record this
    //
    //    info.bounds_shifted = true;
    info.bounds_perturbed = true;
  }
  ekk_instance_.invalidatePrimalMaxSumInfeasibilityRecord();
}

void HEkkPrimal::phase2UpdatePrimal(const bool initialise) {
  if (initialise) {
    max_max_local_primal_infeasibility_ = 0;
    max_max_ignored_violation_ = 0;
    return;
  }
  analysis->simplexTimerStart(UpdatePrimalClock);
  HighsSimplexInfo& info = ekk_instance_.info_;
  bool primal_infeasible = false;
  double max_local_primal_infeasibility = 0;
  double max_ignored_violation = 0;
  // If shifts are only identified in rebuild() the bounds can be
  // ignored. If they aren't ignored, then violations lead to either
  // identification of infeasibilities (and return to Phase 1) or
  // shifting of bounds to accommodate them.
  //
  // primal_correction_strategy is defined (const) as
  // kSimplexPrimalCorrectionStrategyAlways
  const bool ignore_bounds =
      primal_correction_strategy == kSimplexPrimalCorrectionStrategyInRebuild;
  assert(primal_correction_strategy == kSimplexPrimalCorrectionStrategyAlways);
  HighsInt to_entry;
  const bool use_col_indices = ekk_instance_.simplex_nla_.sparseLoopStyle(
      col_aq.count, num_row, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iRow = use_col_indices ? col_aq.index[iEntry] : iEntry;
    info.baseValue_[iRow] -= theta_primal * col_aq.array[iRow];
    //    if (ignore_bounds) continue;
    // Determine whether a bound is violated and take action
    double lower = info.baseLower_[iRow];
    double upper = info.baseUpper_[iRow];
    double value = info.baseValue_[iRow];
    HighsInt bound_violated = 0;
    if (value < lower - primal_feasibility_tolerance) {
      bound_violated = -1;
    } else if (value > upper + primal_feasibility_tolerance) {
      bound_violated = 1;
    }
    if (!bound_violated) continue;
    // A bound is violated
    if (primal_correction_strategy == kSimplexPrimalCorrectionStrategyNone) {
      // Not used
      assert(111 == 222);
      // @primal_infeasibility calculation
      double primal_infeasibility;
      if (bound_violated < 0) {
        primal_infeasibility = lower - value;
      } else {
        primal_infeasibility = value - upper;
      }
      max_local_primal_infeasibility =
          max(primal_infeasibility, max_local_primal_infeasibility);
      if (primal_infeasibility > primal_feasibility_tolerance) {
        info.num_primal_infeasibilities++;
        primal_infeasible = true;
      }
    } else if (ignore_bounds) {
      // Not used
      assert(111 == 333);
      double ignored_violation;
      if (bound_violated < 0) {
        ignored_violation = lower - value;
      } else {
        ignored_violation = value - upper;
      }
      max_ignored_violation = max(ignored_violation, max_ignored_violation);
    } else {
      HighsInt iCol = ekk_instance_.basis_.basicIndex_[iRow];
      double bound_shift;
      if (bound_violated > 0) {
        // Perturb the upper bound to accommodate the infeasibility
        shiftBound(false, iCol, info.baseValue_[iRow],
                   info.numTotRandomValue_[iCol], info.workUpper_[iCol],
                   bound_shift);
        info.baseUpper_[iRow] = info.workUpper_[iCol];
        info.workUpperShift_[iCol] += bound_shift;
      } else {
        // Perturb the lower bound to accommodate the infeasibility
        shiftBound(true, iCol, info.baseValue_[iRow],
                   info.numTotRandomValue_[iCol], info.workLower_[iCol],
                   bound_shift);
        info.baseLower_[iRow] = info.workLower_[iCol];
        info.workLowerShift_[iCol] += bound_shift;
      }
      info.bounds_shifted = true;
      assert(bound_shift > 0);
    }
  }
  if (primal_infeasible) {
    rebuild_reason = kRebuildReasonPrimalInfeasibleInPrimalSimplex;
    if (max_local_primal_infeasibility >
        max_max_local_primal_infeasibility_ * 2) {
      max_max_local_primal_infeasibility_ = max_local_primal_infeasibility;
      printf("phase2UpdatePrimal: max_local_primal_infeasibility = %g\n",
             max_local_primal_infeasibility);
    }
    ekk_instance_.invalidatePrimalMaxSumInfeasibilityRecord();
  }
  if (max_ignored_violation > max_max_ignored_violation_ * 2) {
    max_max_ignored_violation_ = max_ignored_violation;
    printf("phase2UpdatePrimal: max_ignored_violation = %g\n",
           max_ignored_violation);
  }
  info.updated_primal_objective_value +=
      info.workDual_[variable_in] * theta_primal;

  analysis->simplexTimerStop(UpdatePrimalClock);
}

bool HEkkPrimal::correctPrimal(const bool initialise) {
  if (primal_correction_strategy == kSimplexPrimalCorrectionStrategyNone)
    return true;
  if (initialise) {
    max_max_primal_correction_ = 0;
    return true;
  }
  assert(solve_phase == kSolvePhase2);
  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsInt num_primal_correction = 0;
  double max_primal_correction = 0;
  double sum_primal_correction = 0;
  HighsInt num_primal_correction_skipped = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    double lower = info.baseLower_[iRow];
    double upper = info.baseUpper_[iRow];
    double value = info.baseValue_[iRow];
    HighsInt bound_violated = 0;
    if (value < lower - primal_feasibility_tolerance) {
      bound_violated = -1;
    } else if (value > upper + primal_feasibility_tolerance) {
      bound_violated = 1;
    }
    if (bound_violated) {
      if (info.allow_bound_perturbation) {
        HighsInt iCol = ekk_instance_.basis_.basicIndex_[iRow];
        double bound_shift;
        if (bound_violated > 0) {
          // Perturb the upper bound to accommodate the infeasibility
          shiftBound(false, iCol, info.baseValue_[iRow],
                     info.numTotRandomValue_[iCol], info.workUpper_[iCol],
                     bound_shift);
          info.baseUpper_[iRow] = info.workUpper_[iCol];
          info.workUpperShift_[iCol] += bound_shift;
        } else {
          // Perturb the lower bound to accommodate the infeasibility
          shiftBound(true, iCol, info.baseValue_[iRow],
                     info.numTotRandomValue_[iCol], info.workLower_[iCol],
                     bound_shift);
          info.baseLower_[iRow] = info.workLower_[iCol];
          info.workLowerShift_[iCol] += bound_shift;
        }
        assert(bound_shift > 0);
        num_primal_correction++;
        max_primal_correction = max(bound_shift, max_primal_correction);
        sum_primal_correction += bound_shift;
        // Surely better to record this
        //
        // info.bounds_shifted = true;
        info.bounds_perturbed = true;
      } else {
        // Bound perturbation is not permitted
        num_primal_correction_skipped++;
      }
    }
  }
  if (num_primal_correction_skipped) {
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kError,
                "correctPrimal: Missed %d bound shifts\n",
                num_primal_correction_skipped);
    if (kAllowDeveloperAssert) {
      assert(!num_primal_correction_skipped);
    }
    return false;
  }
  if (max_primal_correction > 2 * max_max_primal_correction_) {
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "phase2CorrectPrimal: num / max / sum primal corrections = "
                "%" HIGHSINT_FORMAT
                " / %g / "
                "%g\n",
                num_primal_correction, max_primal_correction,
                sum_primal_correction);
    max_max_primal_correction_ = max_primal_correction;
  }
  return true;
}

void HEkkPrimal::basicFeasibilityChangeUpdateDual() {
  analysis->simplexTimerStart(UpdateDualBasicFeasibilityChangeClock);
  HighsSimplexInfo& info = ekk_instance_.info_;
  // For basic logicals, the change in the basic cost will be a
  // component in col_basic_feasibility_change. This will lead to it being
  // subtracted from workDual in the loop below over the
  // nonzeros in col_basic_feasibility_change, so add it in now. For basic
  // structurals, there will be no corresponding component in
  // row_basic_feasibility_change, since only the nonbasic components are
  // computed (avoided using row pricing, and basic components
  // zeroed after column pricing). Hence there will be no
  // subtraction in the loop below over the nonzeros in
  // row_basic_feasibility_change. Hence, only add in the basic cost change
  // for logicals.
  //
  // Assumes that row_basic_feasibility_change has been set up in
  // phase1UpdatePrimal()

  basicFeasibilityChangeBtran();
  basicFeasibilityChangePrice();
  HighsInt to_entry;
  const bool use_row_indices = ekk_instance_.simplex_nla_.sparseLoopStyle(
      row_basic_feasibility_change.count, num_col, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iCol =
        use_row_indices ? row_basic_feasibility_change.index[iEntry] : iEntry;
    info.workDual_[iCol] -= row_basic_feasibility_change.array[iCol];
  }
  const bool use_col_indices = ekk_instance_.simplex_nla_.sparseLoopStyle(
      col_basic_feasibility_change.count, num_row, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iRow =
        use_col_indices ? col_basic_feasibility_change.index[iEntry] : iEntry;
    HighsInt iCol = num_col + iRow;
    info.workDual_[iCol] -= col_basic_feasibility_change.array[iRow];
  }
  ekk_instance_.invalidateDualInfeasibilityRecord();
  analysis->simplexTimerStop(UpdateDualBasicFeasibilityChangeClock);
}

void HEkkPrimal::basicFeasibilityChangeBtran() {
  // Performs BTRAN on col_basic_feasibility_change. Make sure that
  // col_basic_feasibility_change.count is large (>lp_.num_row_ to be
  // sure) rather than 0 if the indices of the RHS (and true value of
  // col_basic_feasibility_change.count) isn't known.
  analysis->simplexTimerStart(BtranBasicFeasibilityChangeClock);
  const HighsInt solver_num_row = ekk_instance_.lp_.num_row_;
  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordBefore(
        kSimplexNlaBtranBasicFeasibilityChange, col_basic_feasibility_change,
        ekk_instance_.info_.col_basic_feasibility_change_density);
  ekk_instance_.simplex_nla_.btran(
      col_basic_feasibility_change,
      ekk_instance_.info_.col_basic_feasibility_change_density,
      analysis->pointer_serial_factor_clocks);

  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordAfter(kSimplexNlaBtranBasicFeasibilityChange,
                                   col_basic_feasibility_change);
  const double local_col_basic_feasibility_change_density =
      (double)col_basic_feasibility_change.count / solver_num_row;
  ekk_instance_.updateOperationResultDensity(
      local_col_basic_feasibility_change_density,
      ekk_instance_.info_.col_basic_feasibility_change_density);
  analysis->simplexTimerStop(BtranBasicFeasibilityChangeClock);
}

void HEkkPrimal::basicFeasibilityChangePrice() {
  analysis->simplexTimerStart(PriceBasicFeasibilityChangeClock);
  const HighsSimplexInfo& info = ekk_instance_.info_;
  const double local_density =
      1.0 * col_basic_feasibility_change.count / num_row;
  bool use_col_price;
  bool use_row_price_w_switch;
  ekk_instance_.choosePriceTechnique(info.price_strategy, local_density,
                                     use_col_price, use_row_price_w_switch);
  if (analysis->analyse_simplex_summary_data) {
    if (use_col_price) {
      const double expected_density = 1;
      analysis->operationRecordBefore(kSimplexNlaPriceBasicFeasibilityChange,
                                      col_basic_feasibility_change,
                                      expected_density);
      analysis->num_col_price++;
    } else if (use_row_price_w_switch) {
      analysis->operationRecordBefore(
          kSimplexNlaPriceBasicFeasibilityChange, col_basic_feasibility_change,
          ekk_instance_.info_.col_basic_feasibility_change_density);
      analysis->num_row_price_with_switch++;
    } else {
      analysis->operationRecordBefore(
          kSimplexNlaPriceBasicFeasibilityChange, col_basic_feasibility_change,
          ekk_instance_.info_.col_basic_feasibility_change_density);
      analysis->num_row_price++;
    }
  }
  row_basic_feasibility_change.clear();
  const bool quad_precision = false;
  if (use_col_price) {
    // Perform column-wise PRICE
    ekk_instance_.lp_.a_matrix_.priceByColumn(quad_precision,
                                              row_basic_feasibility_change,
                                              col_basic_feasibility_change);
  } else if (use_row_price_w_switch) {
    // Perform hyper-sparse row-wise PRICE, but switch if the density of
    // row_basic_feasibility_change becomes extreme
    //
    const double switch_density = kHyperPriceDensity;
    ekk_instance_.ar_matrix_.priceByRowWithSwitch(
        quad_precision, row_basic_feasibility_change,
        col_basic_feasibility_change, info.row_basic_feasibility_change_density,
        0, switch_density);
  } else {
    // Perform hyper-sparse row-wise PRICE
    ekk_instance_.ar_matrix_.priceByRow(quad_precision,
                                        row_basic_feasibility_change,
                                        col_basic_feasibility_change);
  }
  if (use_col_price) {
    // Column-wise PRICE computes components corresponding to basic
    // variables, so zero these by exploiting the fact that, for basic
    // variables, nonbasicFlag[*]=0
    const std::vector<int8_t>& nonbasicFlag =
        ekk_instance_.basis_.nonbasicFlag_;
    for (HighsInt iCol = 0; iCol < num_col; iCol++)
      row_basic_feasibility_change.array[iCol] *= nonbasicFlag[iCol];
  }
  // Update the record of average row_basic_feasibility_change density
  const double local_row_basic_feasibility_change_density =
      (double)row_basic_feasibility_change.count / num_col;
  ekk_instance_.updateOperationResultDensity(
      local_row_basic_feasibility_change_density,
      ekk_instance_.info_.row_basic_feasibility_change_density);
  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordAfter(kSimplexNlaPriceBasicFeasibilityChange,
                                   row_basic_feasibility_change);
  analysis->simplexTimerStop(PriceBasicFeasibilityChangeClock);
}

void HEkkPrimal::initialiseDevexFramework() {
  edge_weight_.assign(num_tot, 1.0);
  devex_index_.assign(num_tot, 0);
  for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
    const HighsInt nonbasicFlag = ekk_instance_.basis_.nonbasicFlag_[iCol];
    devex_index_[iCol] = nonbasicFlag * nonbasicFlag;
  }
  num_devex_iterations_ = 0;
  num_bad_devex_weight_ = 0;
  if (report_hyper_chuzc) printf("initialiseDevexFramework\n");
  hyperChooseColumnClear();
}

void HEkkPrimal::updateDevex() {
  analysis->simplexTimerStart(DevexUpdateWeightClock);
  // Compute the pivot weight from the reference set
  double dPivotWeight = 0.0;
  HighsInt to_entry;
  const bool use_col_indices = ekk_instance_.simplex_nla_.sparseLoopStyle(
      col_aq.count, num_row, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iRow = use_col_indices ? col_aq.index[iEntry] : iEntry;
    HighsInt iCol = ekk_instance_.basis_.basicIndex_[iRow];
    double dAlpha = devex_index_[iCol] * col_aq.array[iRow];
    dPivotWeight += dAlpha * dAlpha;
  }
  dPivotWeight += devex_index_[variable_in] * 1.0;

  // Check if the saved weight is too large
  if (edge_weight_[variable_in] > kBadDevexWeightFactor * dPivotWeight)
    num_bad_devex_weight_++;

  // Update the devex weight for all
  double dPivot = col_aq.array[row_out];
  dPivotWeight /= (dPivot * dPivot);

  for (HighsInt iEl = 0; iEl < row_ap.count; iEl++) {
    HighsInt iCol = row_ap.index[iEl];
    double alpha = row_ap.array[iCol];
    double devex = dPivotWeight * alpha * alpha;
    devex += devex_index_[iCol] * 1.0;
    if (edge_weight_[iCol] < devex) {
      edge_weight_[iCol] = devex;
    }
  }
  for (HighsInt iEl = 0; iEl < row_ep.count; iEl++) {
    HighsInt iRow = row_ep.index[iEl];
    HighsInt iCol = iRow + num_col;
    double alpha = row_ep.array[iRow];
    double devex = dPivotWeight * alpha * alpha;
    devex += devex_index_[iCol] * 1.0;
    if (edge_weight_[iCol] < devex) {
      edge_weight_[iCol] = devex;
    }
  }
  // Update devex weight for the pivots
  edge_weight_[variable_out] = max(1.0, dPivotWeight);
  edge_weight_[variable_in] = 1.0;
  num_devex_iterations_++;
  analysis->simplexTimerStop(DevexUpdateWeightClock);
}

void HEkkPrimal::computePrimalSteepestEdgeWeights() {
  const HighsInt report_var = -16;
  edge_weight_.resize(num_tot);
  if (ekk_instance_.logicalBasis()) {
    HighsSparseMatrix& a_matrix = ekk_instance_.lp_.a_matrix_;
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      edge_weight_[iCol] = 1;
      for (HighsInt iEl = a_matrix.start_[iCol];
           iEl < a_matrix.start_[iCol + 1]; iEl++)
        edge_weight_[iCol] += a_matrix.value_[iEl] * a_matrix.value_[iEl];
    }
  } else {
    HVector local_col_aq;
    local_col_aq.setup(num_row);
    for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
      if (ekk_instance_.basis_.nonbasicFlag_[iVar]) {
        edge_weight_[iVar] =
            computePrimalSteepestEdgeWeight(iVar, local_col_aq);
        if (iVar == report_var) {
          printf("Tableau column %d\nRow       Value\n", (int)report_var);
          for (HighsInt iRow = 0; iRow < num_row; iRow++) {
            if (local_col_aq.array[iRow])
              printf("%3d  %10.7g\n", (int)iRow, local_col_aq.array[iRow]);
          }
        }
      }
    }
  }
}

double HEkkPrimal::computePrimalSteepestEdgeWeight(const HighsInt iVar,
                                                   HVector& local_col_aq) {
  local_col_aq.clear();
  ekk_instance_.lp_.a_matrix_.collectAj(local_col_aq, iVar, 1);
  local_col_aq.packFlag = false;
  ekk_instance_.simplex_nla_.ftran(
      local_col_aq, ekk_instance_.info_.col_aq_density,
      ekk_instance_.analysis_.pointer_serial_factor_clocks);
  const double local_col_aq_density =
      (1.0 * local_col_aq.count) / ekk_instance_.lp_.num_row_;
  ekk_instance_.updateOperationResultDensity(
      local_col_aq_density, ekk_instance_.info_.col_aq_density);
  return 1 + local_col_aq.norm2();
}

void HEkkPrimal::updatePrimalSteepestEdgeWeights() {
  // Compute, for all j
  //
  // lambda_j = hat{a}_{pj} / hat{a}_{pq}
  //
  // Then contribution to updated weight j is
  //
  // 1 + (hat{a}_j - lambda_j*hat{a}_q).(hat{a}_j - lambda_j*hat{a}_q)
  //
  // = 1 + hat{a}_j.hat{a}_j - 2*lambda_j*hat{a}_qB^{-1}a_j +
  // lambda_j^2*hat{a}_q.hat{a}_q
  //
  // = w_j - 2*lambda_j*(hat{a}_qB^{-1})a_j + lambda_j^2*||hat{a}_q||^2
  //
  // So need to compute mu = B^{-T}hat{a}_q
  //
  // Note that hat{a}_pj - lambda_j*hat{a}_pq is zero, but the updated
  // tableau entry is lambda_j, so have to add lambda_j*lambda_j
  HighsSparseMatrix& a_matrix = ekk_instance_.lp_.a_matrix_;
  col_steepest_edge.copy(&col_aq);
  updateBtranPSE(col_steepest_edge);
  const double col_aq_squared_2norm = col_aq.norm2();
  const bool report_col_aq = false;
  if (report_col_aq) {
    printf(
        "updatePrimalSteepestEdgeWeights: in = %d; out = %d; ||col_aq||^2 = "
        "%g\n",
        (int)variable_in, (int)variable_out, col_aq_squared_2norm);
    printf("Pivotal column %d\nRow       Value\n", (int)variable_in);
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (col_aq.array[iRow])
        printf("%3d  %10.7g\n", (int)iRow, col_aq.array[iRow]);
    }
  }
  assert(ekk_instance_.basis_.nonbasicFlag_[variable_in]);
  HighsInt iVar;
  double pivotal_row_value;
  for (HighsInt iX = 0; iX < row_ap.count + row_ep.count; iX++) {
    if (iX < row_ap.count) {
      iVar = row_ap.index[iX];
      pivotal_row_value = row_ap.array[iVar];
    } else {
      HighsInt iRow = row_ep.index[iX - row_ap.count];
      iVar = num_col + iRow;
      pivotal_row_value = row_ep.array[iRow];
    }
    if (iVar == variable_in) continue;
    if (!ekk_instance_.basis_.nonbasicFlag_[iVar]) continue;
    const double lambda = pivotal_row_value / alpha_col;
    double mu_aj = 0;
    if (iVar < num_col) {
      for (HighsInt iEl = a_matrix.start_[iVar];
           iEl < a_matrix.start_[iVar + 1]; iEl++)
        mu_aj += col_steepest_edge.array[a_matrix.index_[iEl]] *
                 a_matrix.value_[iEl];
    } else {
      mu_aj = col_steepest_edge.array[iVar - num_col];
    }
    const double min_weight = 1 + lambda * lambda;
    edge_weight_[iVar] +=
        (lambda * lambda * col_aq_squared_2norm - 2 * lambda * mu_aj);
    edge_weight_[iVar] += lambda * lambda;
    if (edge_weight_[iVar] < min_weight) {
      //      printf("Augmenting weight(%2d)=%10.4g to %10.4g\n", (int)iVar,
      //      edge_weight_[iVar], min_weight);
      edge_weight_[iVar] = min_weight;
    }
  }
  // The tableau column for the variable leaving the basis is the
  // pivotal column, divided through by the pivot, except for the
  // value in the pivotal location, which is 1/pivot
  //
  // We have col_aq_squared_2norm = s^2 + pivot^2, where s^2 is the
  // sum of squares of the non-pivotal entries
  //
  // The new weight is s^2/pivot^2 + 1/pivot^2 + 1
  //
  // = (s^2 + pivot^2)/pivot^2 + 1/pivot^2
  //
  // = col_aq_squared_2norm/pivot^2 + 1/pivot^2
  //
  // = (col_aq_squared_2norm + 1) / pivot^2
  edge_weight_[variable_out] =
      (1 + col_aq_squared_2norm) / (alpha_col * alpha_col);
  edge_weight_[variable_in] = 0;
}

void HEkkPrimal::updateDualSteepestEdgeWeights() {
  col_steepest_edge.copy(&row_ep);
  updateFtranDSE(col_steepest_edge);
  std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
  // Compute the weight from row_ep and over-write the updated weight
  if (ekk_instance_.simplex_in_scaled_space_) {
    edge_weight[row_out] = row_ep.norm2();
  } else {
    edge_weight[row_out] =
        ekk_instance_.simplex_nla_.rowEp2NormInScaledSpace(row_out, row_ep);
  }
  const double pivot_in_scaled_space =
      ekk_instance_.simplex_nla_.pivotInScaledSpace(&col_aq, variable_in,
                                                    row_out);
  if (ekk_instance_.simplex_in_scaled_space_)
    assert(pivot_in_scaled_space == alpha_col);
  const double new_pivotal_edge_weight =
      edge_weight[row_out] / (pivot_in_scaled_space * pivot_in_scaled_space);
  const double Kai = -2 / pivot_in_scaled_space;
  ekk_instance_.updateDualSteepestEdgeWeights(row_out, variable_in, &col_aq,
                                              new_pivotal_edge_weight, Kai,
                                              col_steepest_edge.array.data());
  edge_weight[row_out] = new_pivotal_edge_weight;
}

void HEkkPrimal::updateFtranDSE(HVector& col_steepest_edge) {
  // For comments on scaling actions, see HEkkDual::updateFtranDSE
  analysis->simplexTimerStart(FtranDseClock);
  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordBefore(kSimplexNlaFtranDse, col_steepest_edge,
                                    ekk_instance_.info_.row_DSE_density);
  // Apply R{-1}
  ekk_instance_.simplex_nla_.unapplyBasisMatrixRowScale(col_steepest_edge);

  // Perform FTRAN DSE
  ekk_instance_.simplex_nla_.ftranInScaledSpace(
      col_steepest_edge, ekk_instance_.info_.row_DSE_density,
      analysis->pointer_serial_factor_clocks);
  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordAfter(kSimplexNlaFtranDse, col_steepest_edge);
  analysis->simplexTimerStop(FtranDseClock);
  const double local_row_DSE_density =
      (1.0 * col_steepest_edge.count) / num_row;
  ekk_instance_.updateOperationResultDensity(
      local_row_DSE_density, ekk_instance_.info_.row_DSE_density);
}

void HEkkPrimal::updateBtranPSE(HVector& col_steepest_edge) {
  analysis->simplexTimerStart(BtranPseClock);
  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordBefore(
        kSimplexNlaBtranPse, col_steepest_edge,
        ekk_instance_.info_.col_steepest_edge_density);
  // Perform BTRAN PSE
  ekk_instance_.simplex_nla_.btran(
      col_steepest_edge, ekk_instance_.info_.col_steepest_edge_density,
      analysis->pointer_serial_factor_clocks);
  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordAfter(kSimplexNlaBtranPse, col_steepest_edge);
  analysis->simplexTimerStop(BtranPseClock);
  const double local_col_steepest_edge_density =
      (1.0 * col_steepest_edge.count) / num_row;
  ekk_instance_.updateOperationResultDensity(
      local_col_steepest_edge_density,
      ekk_instance_.info_.col_steepest_edge_density);
}

void HEkkPrimal::updateVerify() {
  // updateVerify for primal
  const HighsSimplexInfo& info = ekk_instance_.info_;
  const double numerical_trouble_tolerance = 1e-7;
  numericalTrouble = 0;
  double abs_alpha_from_col = fabs(alpha_col);
  std::string alpha_row_source;
  if (variable_in < num_col) {
    alpha_row = row_ap.array[variable_in];
    alpha_row_source = "Col";
  } else {
    alpha_row = row_ep.array[variable_in - num_col];
    alpha_row_source = "Row";
  }
  double abs_alpha_from_row = fabs(alpha_row);
  double abs_alpha_diff = fabs(abs_alpha_from_col - abs_alpha_from_row);
  double min_abs_alpha = min(abs_alpha_from_col, abs_alpha_from_row);
  numericalTrouble = abs_alpha_diff / min_abs_alpha;
  if (numericalTrouble > numerical_trouble_tolerance)
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "Numerical check: Iter %4" HIGHSINT_FORMAT
                ": alpha_col = %12g, (From %3s alpha_row = "
                "%12g), aDiff = %12g: measure = %12g\n",
                ekk_instance_.iteration_count_, alpha_col,
                alpha_row_source.c_str(), alpha_row, abs_alpha_diff,
                numericalTrouble);
  if (kAllowDeveloperAssert) {
    assert(numericalTrouble < 1e-3);
  }
  // Reinvert if the relative difference is large enough, and updates have been
  // performed
  //
  if (numericalTrouble > 1e-7 && info.update_count > 0)
    rebuild_reason = kRebuildReasonPossiblySingularBasis;
}

void HEkkPrimal::iterationAnalysisData() {
  // Possibly compute the infeasibility data
  if (analysis->analyse_simplex_runtime_data)
    ekk_instance_.computeInfeasibilitiesForReporting(SimplexAlgorithm::kPrimal);
  HighsSimplexInfo& info = ekk_instance_.info_;
  analysis->simplex_strategy = kSimplexStrategyPrimal;
  analysis->edge_weight_mode = edge_weight_mode;
  analysis->solve_phase = solve_phase;
  analysis->simplex_iteration_count = ekk_instance_.iteration_count_;
  analysis->devex_iteration_count = num_devex_iterations_;
  analysis->pivotal_row_index = row_out;
  analysis->leaving_variable = variable_out;
  analysis->entering_variable = variable_in;
  analysis->rebuild_reason = rebuild_reason;
  analysis->reduced_rhs_value = 0;
  analysis->reduced_cost_value = 0;
  analysis->edge_weight = 0;
  analysis->primal_delta = 0;
  analysis->primal_step = theta_primal;
  analysis->dual_step = theta_dual;
  analysis->pivot_value_from_column = alpha_col;
  analysis->pivot_value_from_row = alpha_row;
  analysis->numerical_trouble = numericalTrouble;
  analysis->edge_weight_error = ekk_instance_.edge_weight_error_;
  analysis->objective_value = info.updated_primal_objective_value;
  analysis->num_primal_infeasibility = info.num_primal_infeasibilities;
  analysis->num_dual_infeasibility = info.num_dual_infeasibilities;
  analysis->sum_primal_infeasibility = info.sum_primal_infeasibilities;
  analysis->sum_dual_infeasibility = info.sum_dual_infeasibilities;
  if ((edge_weight_mode == EdgeWeightMode::kDevex) &&
      (num_devex_iterations_ == 0))
    analysis->num_devex_framework++;
  analysis->col_aq_density = info.col_aq_density;
  analysis->row_ep_density = info.row_ep_density;
  analysis->row_ap_density = info.row_ap_density;
  analysis->row_DSE_density = info.row_DSE_density;
  analysis->col_steepest_edge_density = info.col_steepest_edge_density;
  analysis->col_basic_feasibility_change_density =
      info.col_basic_feasibility_change_density;
  analysis->row_basic_feasibility_change_density =
      info.row_basic_feasibility_change_density;
  analysis->col_BFRT_density = info.col_BFRT_density;
  analysis->primal_col_density = info.primal_col_density;
  analysis->dual_col_density = info.dual_col_density;
  analysis->num_costly_DSE_iteration = info.num_costly_DSE_iteration;
  analysis->costly_DSE_measure = info.costly_DSE_measure;
}

void HEkkPrimal::iterationAnalysis() {
  iterationAnalysisData();
  analysis->iterationReport();
  if (analysis->analyse_simplex_summary_data) analysis->iterationRecord();
}

void HEkkPrimal::localReportIterHeader() {
  printf(" Iter ColIn Row_Out ColOut\n");
}

void HEkkPrimal::localReportIter(const bool header) {
  if (!report_hyper_chuzc) return;
  const HighsSimplexInfo& info = ekk_instance_.info_;
  HighsInt iteration_count = ekk_instance_.iteration_count_;
  if (header) {
    localReportIterHeader();
    last_header_iteration_count_ = iteration_count;
  } else {
    if (ekk_instance_.iteration_count_ > last_header_iteration_count_ + 10) {
      localReportIterHeader();
      last_header_iteration_count_ = iteration_count;
    }
    if (row_out >= 0) {
      printf("%5" HIGHSINT_FORMAT " %5" HIGHSINT_FORMAT "  %5" HIGHSINT_FORMAT
             "  %5" HIGHSINT_FORMAT "",
             iteration_count, variable_in, row_out, variable_out);
    } else {
      printf("%5" HIGHSINT_FORMAT " %5" HIGHSINT_FORMAT " Bound flip   ",
             iteration_count, variable_in);
    }
    if (check_column >= 0 && iteration_count >= check_iter) {
      HighsInt flag = ekk_instance_.basis_.nonbasicFlag_[check_column];
      HighsInt move = ekk_instance_.basis_.nonbasicMove_[check_column];
      double lower = info.workLower_[check_column];
      double upper = info.workUpper_[check_column];
      double value;
      if (flag == kNonbasicFlagTrue) {
        value = info.workValue_[check_column];
      } else {
        HighsInt iRow;
        for (iRow = 0; iRow < num_row; iRow++) {
          if (ekk_instance_.basis_.basicIndex_[iRow] == check_column) break;
        }
        assert(iRow < num_row);
        value = info.baseValue_[iRow];
      }
      printf(": Var %2" HIGHSINT_FORMAT " (%1" HIGHSINT_FORMAT
             ", %2" HIGHSINT_FORMAT ") [%9.4g, %9.4g, %9.4g]",
             check_column, flag, move, lower, value, upper);
      if (flag == kNonbasicFlagTrue) {
        double dual = info.workDual_[check_column];
        double weight = edge_weight_[check_column];
        double infeasibility = -move * dual;
        if (lower == -kHighsInf && upper == kHighsInf)
          infeasibility = fabs(dual);
        if (infeasibility < dual_feasibility_tolerance) infeasibility = 0;
        double measure = infeasibility * infeasibility / weight;
        printf(" Du = %9.4g; Wt = %9.4g; Ms = %9.4g", dual, weight, measure);
      }
    }
    printf("\n");
  }
}

void HEkkPrimal::reportRebuild(const HighsInt reason_for_rebuild) {
  analysis->simplexTimerStart(ReportRebuildClock);
  iterationAnalysisData();
  analysis->rebuild_reason = reason_for_rebuild;
  analysis->rebuild_reason_string =
      ekk_instance_.rebuildReason(reason_for_rebuild);
  if (ekk_instance_.options_->output_flag) analysis->invertReport();
  analysis->simplexTimerStop(ReportRebuildClock);
}

void HEkkPrimal::getNonbasicFreeColumnSet() {
  if (!num_free_col) return;
  assert(num_free_col > 0);
  const HighsSimplexInfo& info = ekk_instance_.info_;
  const SimplexBasis& basis = ekk_instance_.basis_;
  nonbasic_free_col_set.clear();
  for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
    bool nonbasic_free = basis.nonbasicFlag_[iCol] == kNonbasicFlagTrue &&
                         info.workLower_[iCol] <= -kHighsInf &&
                         info.workUpper_[iCol] >= kHighsInf;
    if (nonbasic_free) nonbasic_free_col_set.add(iCol);
  }
  //  nonbasic_free_col_set.print();
}

void HEkkPrimal::removeNonbasicFreeColumn() {
  bool remove_nonbasic_free_column =
      ekk_instance_.basis_.nonbasicMove_[variable_in] == 0;
  if (remove_nonbasic_free_column) {
    bool removed_nonbasic_free_column =
        nonbasic_free_col_set.remove(variable_in);
    if (!removed_nonbasic_free_column) {
      highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kError,
                  "HEkkPrimal::phase1update failed to remove nonbasic free "
                  "column %" HIGHSINT_FORMAT "\n",
                  variable_in);
      assert(removed_nonbasic_free_column);
    }
  }
}

void HEkkPrimal::adjustPerturbedEquationOut() {
  if (!ekk_instance_.info_.bounds_perturbed) return;
  const HighsLp& lp = ekk_instance_.lp_;
  HighsSimplexInfo& info = ekk_instance_.info_;
  double lp_lower;
  double lp_upper;
  if (variable_out < num_col) {
    lp_lower = lp.col_lower_[variable_out];
    lp_upper = lp.col_upper_[variable_out];
  } else {
    lp_lower = -lp.row_upper_[variable_out - num_col];
    lp_upper = -lp.row_lower_[variable_out - num_col];
  }
  if (lp_lower < lp_upper) return;
  // Leaving variable is fixed
  //  double save_theta_primal = theta_primal;
  double true_fixed_value = lp_lower;
  // Modify theta_primal so that variable leaves at true fixed value
  theta_primal = (info.baseValue_[row_out] - true_fixed_value) / alpha_col;
  /*
    printf("For equation %4" HIGHSINT_FORMAT " to be nonbasic at RHS %10.4g
    requires theta_primal to change by %10.4g from %10.4g to %10.4g\n",
    variable_out, true_fixed_value, theta_primal-save_theta_primal,
    save_theta_primal, theta_primal);
  */
  info.workLower_[variable_out] = true_fixed_value;
  info.workUpper_[variable_out] = true_fixed_value;
  info.workRange_[variable_out] = 0;
  value_in = info.workValue_[variable_in] + theta_primal;
}

void HEkkPrimal::getBasicPrimalInfeasibility() {
  // Gets the num/max/sum of basic primal infeasibilities,
  analysis->simplexTimerStart(ComputePrIfsClock);
  const double primal_feasibility_tolerance =
      ekk_instance_.options_->primal_feasibility_tolerance;
  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsInt& num_primal_infeasibility = info.num_primal_infeasibilities;
  double& max_primal_infeasibility = info.max_primal_infeasibility;
  double& sum_primal_infeasibility = info.sum_primal_infeasibilities;
  const HighsInt updated_num_primal_infeasibility = num_primal_infeasibility;
  num_primal_infeasibility = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibility = 0;

  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    double value = info.baseValue_[iRow];
    double lower = info.baseLower_[iRow];
    double upper = info.baseUpper_[iRow];
    // @primal_infeasibility calculation
    double primal_infeasibility = 0;
    if (value < lower - primal_feasibility_tolerance) {
      primal_infeasibility = lower - value;
    } else if (value > upper + primal_feasibility_tolerance) {
      primal_infeasibility = value - upper;
    }
    if (primal_infeasibility > 0) {
      if (primal_infeasibility > primal_feasibility_tolerance)
        num_primal_infeasibility++;
      max_primal_infeasibility =
          std::max(primal_infeasibility, max_primal_infeasibility);
      sum_primal_infeasibility += primal_infeasibility;
    }
  }
  if (updated_num_primal_infeasibility >= 0) {
    // The number of primal infeasibilities should be correct
    bool num_primal_infeasibility_ok =
        num_primal_infeasibility == updated_num_primal_infeasibility;
    // if (!num_primal_infeasibility_ok)
    //   printf("In iteration %" HIGHSINT_FORMAT
    //          ": num_primal_infeasibility = %" HIGHSINT_FORMAT
    //          " != %" HIGHSINT_FORMAT
    //          " = "
    //          "updated_num_primal_infeasibility\n",
    //          ekk_instance_.iteration_count_, num_primal_infeasibility,
    //          updated_num_primal_infeasibility);

    assert(num_primal_infeasibility_ok);
  }
  analysis->simplexTimerStop(ComputePrIfsClock);
}

void HEkkPrimal::shiftBound(const bool lower, const HighsInt iVar,
                            const double value, const double random_value,
                            double& bound, double& shift) {
  // If infeasibility is very large, then adding feasibility may not
  // yield a new value (see #1144) so new_infeasibility < 0 is false,
  // tripping the old assert
  //
  // Ambros proposed adding feasibility
  // *(infeasibility/scale_threshold) when
  // infeasibility/scale_threshold > 1, but this could lead to a large
  // value for new_infeasibility. Sounds better to accept degeneracy
  // in this edge case.
  double feasibility = (1 + random_value) * primal_feasibility_tolerance;
  double old_bound = bound;
  std::string type;
  double infeasibility;
  double new_infeasibility;
  if (lower) {
    // Bound to shift is lower
    type = "lower";
    assert(value < bound - primal_feasibility_tolerance);
    infeasibility = bound - value;
    assert(infeasibility > 0);
    // Determine the amount by which value will be feasible - so that
    // (ideally) it's not degenerate
    shift = infeasibility + feasibility;
    bound -= shift;
    new_infeasibility = bound - value;
  } else {
    // Bound to shift is upper
    type = "upper";
    assert(value > bound + primal_feasibility_tolerance);
    infeasibility = value - bound;
    assert(infeasibility > 0);
    // Determine the amount by which value will be feasible - so that
    // (ideally) it's not degenerate
    shift = infeasibility + feasibility;
    bound += shift;
    new_infeasibility = value - bound;
  }
  if (new_infeasibility > 0) {
    // new_infeasibility should be non-positive, and negative unless
    // bound is excessively large, whereas feasibility is positive
    double error = std::fabs(new_infeasibility + feasibility);
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "HEkkPrimal::shiftBound Value(%4d) = %10.4g exceeds %s: "
                "old_bound = %g; random_value = %g; value = %g; "
                "feasibility = %g; infeasibility = %g; shift = %g; bound = %g; "
                "new_infeasibility = %g with error %g\n",
                int(iVar), value, type.c_str(), old_bound, random_value, value,
                feasibility, infeasibility, shift, bound, new_infeasibility,
                error);
    fflush(stdout);
  }
  assert(new_infeasibility <= 0);
}

void HEkkPrimal::savePrimalRay() {
  assert(variable_in >= 0);
  assert(move_in != kNoRaySign);
  ekk_instance_.primal_ray_record_.clear();
  ekk_instance_.primal_ray_record_.index = variable_in;
  ekk_instance_.primal_ray_record_.sign = -move_in;
}

HighsDebugStatus HEkkPrimal::debugPrimalSimplex(const std::string& message,
                                                const bool initialise) {
  HighsDebugStatus return_status =
      ekk_instance_.debugSimplex(message, algorithm, solve_phase, initialise);
  if (return_status == HighsDebugStatus::kLogicalError) return return_status;
  if (initialise) return return_status;
  return_status = ekk_instance_.debugNonbasicFreeColumnSet(
      num_free_col, nonbasic_free_col_set);
  if (return_status == HighsDebugStatus::kLogicalError) return return_status;
  return HighsDebugStatus::kOk;
}

HighsDebugStatus HEkkPrimal::debugPrimalSteepestEdgeWeights(
    const std::string& message) {
  // Possibly force the expensive check for development work
  const bool check_primal_edge_weights = true;
  if (check_primal_edge_weights) {
    const bool check_all_primal_edge_weights = false;
    const HighsInt alt_debug_level = check_all_primal_edge_weights
                                         ? (HighsInt)kHighsDebugLevelExpensive
                                         : (HighsInt)kHighsDebugLevelCostly;
    //    printf("\nPerforming level %1d check %s for primal steepest edge
    //    weights\n", (int)alt_debug_level, message.c_str());
    return debugPrimalSteepestEdgeWeights(alt_debug_level);
  } else {
    return debugPrimalSteepestEdgeWeights();
  }
}

HighsDebugStatus HEkkPrimal::debugPrimalSteepestEdgeWeights(
    const HighsInt alt_debug_level) {
  const HighsInt use_debug_level =
      alt_debug_level >= 0 ? alt_debug_level
                           : ekk_instance_.options_->highs_debug_level;
  if (use_debug_level < kHighsDebugLevelCostly)
    return HighsDebugStatus::kNotChecked;
  const HighsLp& lp = ekk_instance_.lp_;
  const HighsInt num_row = lp.num_row_;
  const std::vector<int8_t> nonbasic_flag = ekk_instance_.basis_.nonbasicFlag_;
  double primal_steepest_edge_weight_norm = 0;
  double primal_steepest_edge_weight_error = 0;
  HighsInt num_check_weight;
  HVector local_col_aq;
  local_col_aq.setup(num_row);
  if (use_debug_level < kHighsDebugLevelExpensive) {
    for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
      primal_steepest_edge_weight_norm +=
          std::fabs(nonbasic_flag[iVar] * edge_weight_[iVar]);
    }
    // Just check a few weights
    num_check_weight =
        std::max((HighsInt)1, std::min((HighsInt)10, num_tot / 10));
    for (HighsInt iCheck = 0; iCheck < num_check_weight; iCheck++) {
      HighsInt iVar;
      for (;;) {
        iVar = random_.integer(num_tot);
        if (nonbasic_flag[iVar]) break;
      }
      const double true_weight =
          computePrimalSteepestEdgeWeight(iVar, local_col_aq);
      primal_steepest_edge_weight_error +=
          std::fabs(edge_weight_[iVar] - true_weight);
    }
  } else {
    // Check all weights
    num_check_weight = num_col;
    std::vector<double> updated_primal_edge_weight = edge_weight_;
    computePrimalSteepestEdgeWeights();
    for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
      if (!nonbasic_flag[iVar]) continue;
      primal_steepest_edge_weight_norm += std::fabs(edge_weight_[iVar]);
      const double error =
          std::fabs(updated_primal_edge_weight[iVar] - edge_weight_[iVar]);
      if (error > 1e-4)
        printf(
            "debugPrimalSteepestEdgeWeights: var = %2d; weight (true = %10.4g; "
            "updated = %10.4g) error = %10.4g\n",
            (int)iVar, edge_weight_[iVar], updated_primal_edge_weight[iVar],
            error);
      primal_steepest_edge_weight_error += error;
    }
    edge_weight_ = updated_primal_edge_weight;
  }
  // Now assess the relative error
  assert(primal_steepest_edge_weight_norm > 0);
  double relative_primal_steepest_edge_weight_error =
      primal_steepest_edge_weight_error / primal_steepest_edge_weight_norm;
  const double large_relative_primal_steepest_edge_weight_error = 1e-3;
  if (relative_primal_steepest_edge_weight_error >
      10 * debug_max_relative_primal_steepest_edge_weight_error) {
    printf(
        "HEkk::debugPrimalSteepestEdgeWeights Iteration %5d: Checked %2d "
        "weights: "
        "error = %10.4g; norm = %10.4g; relative error = %10.4g\n",
        (int)ekk_instance_.iteration_count_, (int)num_check_weight,
        primal_steepest_edge_weight_error, primal_steepest_edge_weight_norm,
        relative_primal_steepest_edge_weight_error);
    debug_max_relative_primal_steepest_edge_weight_error =
        relative_primal_steepest_edge_weight_error;
    if (relative_primal_steepest_edge_weight_error >
        large_relative_primal_steepest_edge_weight_error)
      return HighsDebugStatus::kLargeError;
  }
  return HighsDebugStatus::kOk;
}

bool HEkkPrimal::isBadBasisChange() {
  return ekk_instance_.isBadBasisChange(SimplexAlgorithm::kPrimal, variable_in,
                                        row_out, rebuild_reason);
}
