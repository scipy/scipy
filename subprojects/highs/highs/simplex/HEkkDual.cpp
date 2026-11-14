/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkDual.cpp
 * @brief
 */
#include "simplex/HEkkDual.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <set>

#include "lp_data/HighsLpUtils.h"
#include "parallel/HighsParallel.h"
#include "simplex/HEkkPrimal.h"
#include "simplex/SimplexTimer.h"

using std::fabs;

HighsStatus HEkkDual::solve(const bool pass_force_phase2) {
  // Initialise control data for a particular solve
  initialiseSolve();

  if (debugDualSimplex("Initialise", true) == HighsDebugStatus::kLogicalError)
    return ekk_instance_.returnFromSolve(HighsStatus::kError);
  // Assumes that the LP has a positive number of rows
  if (ekk_instance_.isUnconstrainedLp())
    return ekk_instance_.returnFromSolve(HighsStatus::kError);

  HighsOptions& options = *ekk_instance_.options_;
  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsSimplexStatus& status = ekk_instance_.status_;
  HighsModelStatus& model_status = ekk_instance_.model_status_;

  if (!dualInfoOk(ekk_instance_.lp_)) {
    highsLogDev(options.log_options, HighsLogType::kError,
                "HPrimalDual::solve has error in dual information\n");
    return ekk_instance_.returnFromSolve(HighsStatus::kError);
  }

  // Possibly use Li dual steepest edge weights by not storing squared
  // primal infeasibilities
  possiblyUseLiDualSteepestEdge();

  assert(status.has_invert);
  if (!status.has_invert) {
    highsLogDev(options.log_options, HighsLogType::kError,
                "HDual:: Should enter solve with INVERT\n");
    return ekk_instance_.returnFromSolve(HighsStatus::kError);
  }

  // Determine the duals without cost perturbation
  ekk_instance_.initialiseCost(SimplexAlgorithm::kDual, kSolvePhaseUnknown);
  ekk_instance_.computeDual();
  ekk_instance_.computeSimplexDualInfeasible();
  // Record whether the solution with unperturbed costs is dual feasible
  const bool dual_feasible_with_unperturbed_costs =
      info.num_dual_infeasibilities == 0;
  // Force phase 2 if dual infeasibilities without cost perturbation
  // involved fixed variables or were (at most) small
  force_phase2 = pass_force_phase2 ||
                 info.max_dual_infeasibility * info.max_dual_infeasibility <
                     ekk_instance_.options_->dual_feasibility_tolerance;
  // Within the MIP solver, unless the basis supplied was alien, the
  // simplex solver should be able to start from dual feasibility, so
  // possibly debug this property. Note that debug_dual_feasible is
  // set in HEkk::setBasis, and is false if kDebugMipNodeDualFeasible
  // is false
  if (ekk_instance_.debug_dual_feasible &&
      !dual_feasible_with_unperturbed_costs) {
    SimplexBasis& basis = ekk_instance_.basis_;
    highsLogDev(
        options.log_options, HighsLogType::kWarning,
        "Basis should be dual feasible, but duals without cost perturbation "
        "have num / max / sum = %4d / %g / %g infeasibilities",
        (int)info.num_dual_infeasibilities, info.max_dual_infeasibility,
        info.sum_dual_infeasibilities);
    if (!force_phase2) {
      highsLogDev(options.log_options, HighsLogType::kWarning,
                  " !!Not forcing phase 2!! basis Id = %d; update count = %d; "
                  "name = %s\n",
                  (int)basis.debug_id, (int)basis.debug_update_count,
                  basis.debug_origin_name.c_str());
    } else {
      highsLogDev(options.log_options, HighsLogType::kWarning, "\n");
    }
  }
  // Determine whether the solution is near-optimal. Values 1000 and
  // 1e-3 (ensuring sum<1) are unimportant, as the sum of primal
  // infeasibilities for near-optimal solutions is typically many
  // orders of magnitude smaller than 1, and the sum of primal
  // infeasibilities will be very much larger for non-trivial LPs
  // that are dual feasible for a logical or crash basis.
  //
  // Consider there to be no dual infeasibilities if there are none,
  // or if phase 2 is forced, in which case any dual infeasibilities
  // will be shifted
  const bool no_simplex_dual_infeasibilities =
      dual_feasible_with_unperturbed_costs || force_phase2;
  const bool near_optimal = no_simplex_dual_infeasibilities &&
                            info.num_primal_infeasibilities < 1000 &&
                            info.max_primal_infeasibility < 1e-3;
  if (near_optimal)
    highsLogDev(options.log_options, HighsLogType::kDetailed,
                "Dual feasible with unperturbed costs and num / max / sum "
                "primal infeasibilities of "
                "%" HIGHSINT_FORMAT
                " / %g "
                "/ %g, so near-optimal\n",
                info.num_primal_infeasibilities, info.max_primal_infeasibility,
                info.sum_primal_infeasibilities);

  // Perturb costs according to whether the solution is near-optimal
  const bool perturb_costs = !near_optimal;
  if (!perturb_costs)
    highsLogDev(options.log_options, HighsLogType::kDetailed,
                "Near-optimal, so don't use cost perturbation\n");
  ekk_instance_.initialiseCost(SimplexAlgorithm::kDual, kSolvePhaseUnknown,
                               perturb_costs);
  // Check whether the time/iteration limit has been reached. First
  // point at which a non-error return can occur
  if (ekk_instance_.bailout())
    return ekk_instance_.returnFromSolve(HighsStatus::kWarning);

  // Consider initialising edge weights
  if (status.has_dual_steepest_edge_weights) {
    // Dual steepest edge weights are known, so possibly check
    assert((HighsInt)ekk_instance_.dual_edge_weight_.size() >= solver_num_row);
    assert((HighsInt)ekk_instance_.scattered_dual_edge_weight_.size() >=
           solver_num_tot);
    ekk_instance_.devDebugDualSteepestEdgeWeights("before solve");
  } else {
    // Set up edge weights
    //
    // Assign unit weights - necessary for Dantzig and Devex, and
    // correct for steepest edge when B=I so, for clarity, do it for
    // all. Also ensure that the scattering vector is the right size
    ekk_instance_.dual_edge_weight_.assign(solver_num_row, 1.0);
    ekk_instance_.scattered_dual_edge_weight_.resize(solver_num_tot);
    if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
      // Intending to using dual steepest edge weights
      //
      // Exact DSE weights need to be computed if the basis contains
      // structurals
      if (ekk_instance_.logicalBasis()) {
        // Unit weights already set up for B=I
        status.has_dual_steepest_edge_weights = true;
      } else {
        // Non-logical basis
        if (near_optimal) {
          // Use Devex rather than compute steepest edge weights
          highsLogDev(
              options.log_options, HighsLogType::kDetailed,
              "Basis is not logical, but near-optimal, so use Devex rather "
              "than compute steepest edge weights\n");
          edge_weight_mode = EdgeWeightMode::kDevex;
          assert(!status.has_dual_steepest_edge_weights);
        } else {
          // Compute steepest edge weights
          highsLogDev(
              options.log_options, HighsLogType::kDetailed,
              "Basis is not logical, so compute steepest edge weights\n");
          ekk_instance_.computeDualSteepestEdgeWeights(true);
          status.has_dual_steepest_edge_weights = true;
        }
      }
    }
    if (edge_weight_mode == EdgeWeightMode::kDevex) initialiseDevexFramework();
    // Check on consistency between edge_weight_mode and
    // status.has_dual_steepest_edge_weights
    if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
      assert(status.has_dual_steepest_edge_weights);
    } else {
      assert(!status.has_dual_steepest_edge_weights);
    }
  }
  // Resize the copy of scattered edge weights for backtracking
  info.backtracking_basis_edge_weight_.resize(solver_num_tot);

  if (perturb_costs) {
    // Compute the dual values with perturbed costs
    ekk_instance_.computeDual();
    // Determine the number of dual infeasibilities after fixed
    // variable flips
    computeDualInfeasibilitiesWithFixedVariableFlips();
    dualInfeasCount = info.num_dual_infeasibilities;
  }

  // Determine the solve phase
  if (force_phase2) {
    // Dual infeasibilities without cost perturbation involved
    // fixed variables or were (at most) small, so can easily be
    // removed by flips for fixed variables and shifts for the rest
    solve_phase = kSolvePhase2;
  } else {
    // Phase depends on the number of dual infeasibilities after fixed
    // variable flips
    solve_phase = dualInfeasCount > 0 ? kSolvePhase1 : kSolvePhase2;
  }
  if (ekk_instance_.debugOkForSolve(SimplexAlgorithm::kDual, solve_phase) ==
      HighsDebugStatus::kLogicalError)
    return ekk_instance_.returnFromSolve(HighsStatus::kError);
  //
  // The major solving loop
  //
  while (solve_phase) {
    HighsInt it0 = ekk_instance_.iteration_count_;
    // When starting a new phase the (updated) dual objective function
    // value isn't known. Indicate this so that when the value
    // computed from scratch in rebuild() isn't checked against the
    // the updated value
    status.has_dual_objective_value = false;
    if (solve_phase == kSolvePhaseUnknown) {
      // Reset the phase 2 bounds so that true number of dual
      // infeasibilities can be determined
      ekk_instance_.initialiseBound(SimplexAlgorithm::kDual,
                                    kSolvePhaseUnknown);
      ekk_instance_.initialiseNonbasicValueAndMove();
      // Determine the number of unavoidable dual infeasibilities, and
      // hence the solve phase
      computeDualInfeasibilitiesWithFixedVariableFlips();
      dualInfeasCount = info.num_dual_infeasibilities;
      solve_phase = dualInfeasCount > 0 ? kSolvePhase1 : kSolvePhase2;
      if (info.backtracking_) {
        // Backtracking, so set the bounds and primal values
        ekk_instance_.initialiseBound(SimplexAlgorithm::kDual, solve_phase);
        ekk_instance_.initialiseNonbasicValueAndMove();
        // Can now forget that we might have been backtracking
        info.backtracking_ = false;
      }
    }
    assert(solve_phase == kSolvePhase1 || solve_phase == kSolvePhase2);
    if (solve_phase == kSolvePhase1) {
      // Phase 1
      analysis->simplexTimerStart(SimplexDualPhase1Clock);
      solvePhase1();
      analysis->simplexTimerStop(SimplexDualPhase1Clock);
      info.dual_phase1_iteration_count +=
          (ekk_instance_.iteration_count_ - it0);
    } else if (solve_phase == kSolvePhase2) {
      // Phase 2
      analysis->simplexTimerStart(SimplexDualPhase2Clock);
      solvePhase2();
      analysis->simplexTimerStop(SimplexDualPhase2Clock);
      info.dual_phase2_iteration_count +=
          (ekk_instance_.iteration_count_ - it0);
    } else {
      // Should only be kSolvePhase1 or kSolvePhase2
      model_status = HighsModelStatus::kSolveError;
      return ekk_instance_.returnFromSolve(HighsStatus::kError);
    }
    // Return if bailing out from solve
    if (ekk_instance_.solve_bailout_)
      return ekk_instance_.returnFromSolve(HighsStatus::kWarning);
    // Can have all possible cases of solve_phase
    assert(solve_phase >= kSolvePhaseMin && solve_phase <= kSolvePhaseMax);
    // Look for scenarios when the major solving loop ends
    if (solve_phase == kSolvePhaseTabooBasis) {
      // Only basis change is taboo
      ekk_instance_.model_status_ = HighsModelStatus::kUnknown;
      return ekk_instance_.returnFromSolve(HighsStatus::kWarning);
    }
    if (solve_phase == kSolvePhaseError) {
      // Solver error so return HighsStatus::kError
      assert(model_status == HighsModelStatus::kSolveError);
      return ekk_instance_.returnFromSolve(HighsStatus::kError);
    }
    if (solve_phase == kSolvePhaseExit) {
      // LP identified as not having an optimal solution
      assert(model_status == HighsModelStatus::kUnboundedOrInfeasible ||
             model_status == HighsModelStatus::kInfeasible);
      break;
    }
    if (solve_phase == kSolvePhaseOptimalCleanup ||
        solve_phase == kSolvePhasePrimalInfeasibleCleanup) {
      // Dual infeasibilities after phase 2 which ends either
      //
      // primal feasible with dual infeasibilities, so use primal
      // simplex to clean up expecting to identify optimality
      //
      // primal infeasible with dual infeasibilities so use primal
      // simplex to clean up expecting to identify (primal)
      // infeasibility
      break;
    }
    // If solve_phase == kSolvePhaseOptimal == 0 then major solving
    // loop ends naturally since solve_phase is false
  }
  // If bailing out, should have returned already
  assert(!ekk_instance_.solve_bailout_);
  // Should only have these cases
  assert(solve_phase == kSolvePhaseExit || solve_phase == kSolvePhaseUnknown ||
         solve_phase == kSolvePhaseOptimal ||
         solve_phase == kSolvePhaseOptimalCleanup ||
         solve_phase == kSolvePhasePrimalInfeasibleCleanup);
  // Can't be solve_phase == kSolvePhase1 since this requires simplex
  // solver to have continued after identifying dual infeasibility.
  if (solve_phase == kSolvePhaseOptimalCleanup ||
      solve_phase == kSolvePhasePrimalInfeasibleCleanup) {
    ekk_instance_.dual_simplex_cleanup_level_++;
    if (solve_phase == kSolvePhasePrimalInfeasibleCleanup) {
      // Primal and dual infeasibilities aren't known when cleaning up
      // after suspected unboundedness in dual phase 2. All that's
      // known is that cost shifting was required to get dual
      // feasibility after removing cost perturbations, and dual
      // simplex iterations may also have been done. This is unlike
      // clean-up of dual infeasibilities after suspected optimality,
      // when no shifting and dual simplex iterations are done after
      // removing cost perturbations.
      //
      // Determine the primal and dual infeasibilities
      ekk_instance_.computeSimplexInfeasible();
    }
    if (ekk_instance_.dual_simplex_cleanup_level_ >
        options.max_dual_simplex_cleanup_level) {
      // No clean up. Dual simplex was optimal or unbounded with
      // unperturbed costs, so say that the scaled LP has been solved
      // optimally. Optimality or infeasibility for the unscaled LP
      // are unlikely but will still be assessed honestly, so leave it
      // to the user to decide whether the solution can be accepted.
      highsLogDev(options.log_options, HighsLogType::kWarning,
                  "HEkkDual:: Cannot use level %" HIGHSINT_FORMAT
                  " primal simplex cleanup for %" HIGHSINT_FORMAT
                  " dual infeasibilities\n",
                  ekk_instance_.dual_simplex_cleanup_level_,
                  info.num_dual_infeasibilities);
      if (solve_phase == kSolvePhaseOptimalCleanup) {
        ekk_instance_.model_status_ = HighsModelStatus::kOptimal;
      } else {
        ekk_instance_.model_status_ = HighsModelStatus::kInfeasible;
      }
    } else {
      // Use primal simplex to clean up. This usually yields
      // optimality or infeasibility (according to whether solve_phase
      // is kSolvePhaseOptimalCleanup or
      // kSolvePhasePrimalInfeasibleCleanup) but can yield
      // unboundedness. Time/iteration limit return is, of course,
      // possible, as are solver error
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "HEkkDual:: Using primal simplex to try to clean up num / "
                  "max / sum = %" HIGHSINT_FORMAT
                  " / %g / %g dual infeasibilities\n",
                  info.num_dual_infeasibilities, info.max_dual_infeasibility,
                  info.sum_dual_infeasibilities);
      HighsStatus return_status = HighsStatus::kOk;
      analysis->simplexTimerStart(SimplexPrimalPhase2Clock);
      // Switch off any bound perturbation
      double save_primal_simplex_bound_perturbation_multiplier =
          info.primal_simplex_bound_perturbation_multiplier;
      info.primal_simplex_bound_perturbation_multiplier = 0;
      HEkkPrimal primal_solver(ekk_instance_);
      HighsStatus call_status = primal_solver.solve(true);
      // Restore any bound perturbation
      info.primal_simplex_bound_perturbation_multiplier =
          save_primal_simplex_bound_perturbation_multiplier;
      analysis->simplexTimerStop(SimplexPrimalPhase2Clock);
      assert(ekk_instance_.called_return_from_solve_);
      return_status = interpretCallStatus(options.log_options, call_status,
                                          return_status, "HEkkPrimal::solve");
      // Reset called_return_from_solve_ to be false, since it's
      // called for this solve
      ekk_instance_.called_return_from_solve_ = false;
      if (return_status != HighsStatus::kOk)
        return ekk_instance_.returnFromSolve(return_status);
      if (ekk_instance_.model_status_ == HighsModelStatus::kOptimal &&
          info.num_primal_infeasibilities + info.num_dual_infeasibilities)
        highsLogDev(options.log_options, HighsLogType::kWarning,
                    "HEkkDual:: Primal simplex clean up yields optimality, "
                    "but with %" HIGHSINT_FORMAT
                    " (max %g) primal infeasibilities and %" HIGHSINT_FORMAT
                    " (max %g) dual infeasibilities\n",
                    info.num_primal_infeasibilities,
                    info.max_primal_infeasibility,
                    info.num_dual_infeasibilities, info.max_dual_infeasibility);
    }
  }
  assert(model_status == HighsModelStatus::kOptimal ||
         model_status == HighsModelStatus::kInfeasible ||
         model_status == HighsModelStatus::kUnbounded ||
         model_status == HighsModelStatus::kUnboundedOrInfeasible);
  if (ekk_instance_.debugOkForSolve(SimplexAlgorithm::kDual, solve_phase) ==
      HighsDebugStatus::kLogicalError)
    return ekk_instance_.returnFromSolve(HighsStatus::kError);
  return ekk_instance_.returnFromSolve(HighsStatus::kOk);
}

void HEkkDual::initialiseInstance() {
  // Called in constructor for HEkkDual class
  // Copy size, matrix and simplex NLA

  solver_num_col = ekk_instance_.lp_.num_col_;
  solver_num_row = ekk_instance_.lp_.num_row_;
  solver_num_tot = solver_num_col + solver_num_row;
  inv_solver_num_row = 1.0 / solver_num_row;

  a_matrix = &ekk_instance_.lp_.a_matrix_;
  simplex_nla = &ekk_instance_.simplex_nla_;
  analysis = &ekk_instance_.analysis_;

  // Copy pointers
  jMove = ekk_instance_.basis_.nonbasicMove_.data();
  workDual = ekk_instance_.info_.workDual_.data();
  workValue = ekk_instance_.info_.workValue_.data();
  workRange = ekk_instance_.info_.workRange_.data();
  baseLower = ekk_instance_.info_.baseLower_.data();
  baseUpper = ekk_instance_.info_.baseUpper_.data();
  baseValue = ekk_instance_.info_.baseValue_.data();

  // Setup local vectors
  col_DSE.setup(solver_num_row);
  col_BFRT.setup(solver_num_row);
  col_aq.setup(solver_num_row);
  row_ep.setup(solver_num_row);
  row_ap.setup(solver_num_col);

  dev_row_ep.setup(solver_num_row);
  dev_col_DSE.setup(solver_num_row);

  // Setup other buffers
  dualRow.setup();
  dualRHS.setup();
}

void HEkkDual::initialiseInstanceParallel(HEkk& simplex) {
  // No need to call this with kSimplexStrategyDualPlain
  if (ekk_instance_.info_.simplex_strategy == kSimplexStrategyDualPlain) return;

  // Identify the (current) number of HiGHS tasks to be used
  const HighsInt num_concurrency = ekk_instance_.info_.num_concurrency;

  HighsInt pass_num_slice;
  if (ekk_instance_.info_.simplex_strategy == kSimplexStrategyDualTasks) {
    // Initialize for tasks
    pass_num_slice = num_concurrency - 2;
    assert(pass_num_slice > 0);
    if (pass_num_slice <= 0) {
      highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kWarning,
                  "SIP trying to use using %" HIGHSINT_FORMAT
                  " slices due to concurrency (%" HIGHSINT_FORMAT
                  ") being too small: results unpredictable\n",
                  pass_num_slice, num_concurrency);
    }
  } else {
    // Initialize for multi
    multi_num = num_concurrency;
    if (multi_num < 1) multi_num = 1;
    if (multi_num > kSimplexConcurrencyLimit)
      multi_num = kSimplexConcurrencyLimit;
    for (HighsInt i = 0; i < multi_num; i++) {
      multi_choice[i].row_out = -1;
      multi_choice[i].row_ep.setup(solver_num_row);
      multi_choice[i].col_aq.setup(solver_num_row);
      multi_choice[i].col_BFRT.setup(solver_num_row);
    }
    pass_num_slice = max(multi_num - 1, HighsInt{1});
    assert(pass_num_slice > 0);
    if (pass_num_slice <= 0) {
      highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kWarning,
                  "PAMI trying to use using %" HIGHSINT_FORMAT
                  " slices due to concurrency (%" HIGHSINT_FORMAT
                  ") being too small: results unpredictable\n",
                  pass_num_slice, num_concurrency);
    }
  }
  // Create the multiple HEkkDualRow instances: one for each column
  // slice
  for (HighsInt i = 0; i < pass_num_slice; i++)
    slice_dualRow.push_back(HEkkDualRow(simplex));
  // Initialise the column slices
  initSlice(pass_num_slice);
  multi_iteration = 0;
}

void HEkkDual::initSlice(const HighsInt initial_num_slice) {
  // Number of slices
  slice_num = initial_num_slice;
  if (slice_num < 1) slice_num = 1;
  assert(slice_num <= kHighsSlicedLimit);
  if (slice_num > kHighsSlicedLimit) {
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kWarning,
                "WARNING: %" HIGHSINT_FORMAT
                " = slice_num > kHighsSlicedLimit = %" HIGHSINT_FORMAT
                " so truncating "
                "slice_num\n",
                slice_num, kHighsSlicedLimit);
    slice_num = kHighsSlicedLimit;
  }

  // Alias to the matrix
  const HighsInt* Astart = a_matrix->start_.data();
  const HighsInt AcountX = Astart[solver_num_col];

  // Figure out partition weight
  double sliced_countX = AcountX / (double)slice_num;
  slice_start[0] = 0;
  for (HighsInt i = 0; i < slice_num - 1; i++) {
    HighsInt endColumn = slice_start[i] + 1;  // At least one column
    HighsInt endX = Astart[endColumn];
    HighsInt stopX = (i + 1) * sliced_countX;
    while (endX < stopX) {
      endX = Astart[++endColumn];
    }
    slice_start[i + 1] = endColumn;
    if (endColumn >= solver_num_col) {
      slice_num = i;  // SHRINK
      break;
    }
  }
  slice_start[slice_num] = solver_num_col;

  // Partition the matrix, row_ap and related packet
  vector<HighsInt> sliced_Astart;
  for (HighsInt i = 0; i < slice_num; i++) {
    // The matrix
    HighsInt from_col = slice_start[i];
    HighsInt to_col = slice_start[i + 1] - 1;
    HighsInt slice_num_col = slice_start[i + 1] - from_col;
    HighsInt from_el = Astart[from_col];
    sliced_Astart.resize(slice_num_col + 1);
    for (HighsInt k = 0; k <= slice_num_col; k++)
      sliced_Astart[k] = Astart[k + from_col] - from_el;
    slice_a_matrix[i].createSlice(ekk_instance_.lp_.a_matrix_, from_col,
                                  to_col);
    slice_ar_matrix[i].createRowwise(slice_a_matrix[i]);

    // The row_ap and its packages
    slice_row_ap[i].setup(slice_num_col);
    slice_dualRow[i].setupSlice(slice_num_col);
  }
}

void HEkkDual::initialiseSolve() {
  // Copy values of simplex solver options to dual simplex options
  primal_feasibility_tolerance =
      ekk_instance_.options_->primal_feasibility_tolerance;
  dual_feasibility_tolerance =
      ekk_instance_.options_->dual_feasibility_tolerance;
  objective_bound = ekk_instance_.options_->objective_bound;

  // Copy tolerances
  // ToDo: Eliminate these horribly-named unnecessary copies!
  Tp = primal_feasibility_tolerance;
  Td = dual_feasibility_tolerance;

  initial_basis_is_logical_ = true;
  for (HighsInt iRow = 0; iRow < solver_num_row; iRow++) {
    if (ekk_instance_.basis_.basicIndex_[iRow] < solver_num_col) {
      initial_basis_is_logical_ = false;
      break;
    }
  }

  interpretDualEdgeWeightStrategy(
      ekk_instance_.info_.dual_edge_weight_strategy);

  // Initialise model and run status values
  ekk_instance_.model_status_ = HighsModelStatus::kNotset;
  ekk_instance_.solve_bailout_ = false;
  ekk_instance_.called_return_from_solve_ = false;
  ekk_instance_.exit_algorithm_ = SimplexAlgorithm::kDual;

  rebuild_reason = kRebuildReasonNo;
}

void HEkkDual::solvePhase1() {
  // Performs dual phase 1 iterations. Returns solve_phase with value
  //
  // kSolvePhaseError => Solver error
  //
  // kSolvePhaseTabooBasis => Only basis change is taboo
  //
  // kSolvePhaseExit => LP identified as dual infeasible
  //
  // kSolvePhaseUnknown => Back-tracking due to singularity
  //
  // kSolvePhase1 => Dual infeasibility suspected, but have to go out
  // and back in to solvePhase1 to perform fresh rebuild. Also if
  // bailing out due to reaching time/iteration limit.
  //
  // kSolvePhase2 => Continue with dual phase 2 iterations

  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsSimplexStatus& status = ekk_instance_.status_;
  HighsModelStatus& model_status = ekk_instance_.model_status_;
  // When starting a new phase the (updated) dual objective function
  // value isn't known. Indicate this so that when the value computed
  // from scratch in build() isn't checked against the updated
  // value
  status.has_primal_objective_value = false;
  status.has_dual_objective_value = false;
  // Set rebuild_reason so that it's assigned when first tested
  rebuild_reason = kRebuildReasonNo;
  // Use to set solve_phase = kSolvePhase1 and ekk_instance_.solve_bailout_ =
  // false so they are set if solvePhase1() is called directly - but it never is
  assert(solve_phase == kSolvePhase1);
  assert(!ekk_instance_.solve_bailout_);
  if (ekk_instance_.bailout()) return;
  // Report the phase start
  highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kDetailed,
              "dual-phase-1-start\n");
  // Switch to dual phase 1 bounds
  ekk_instance_.initialiseBound(SimplexAlgorithm::kDual, solve_phase);
  ekk_instance_.initialiseNonbasicValueAndMove();
  // Free variables are boxed in phase 1, so the list of free variables
  // created at the start of any earlier phase 2 must not be used (#2300)
  dualRow.freeList.clear();

  // If there's no backtracking basis, save the initial basis in case of
  // backtracking
  if (!info.valid_backtracking_basis_) ekk_instance_.putBacktrackingBasis();

  // Main solving structure
  analysis->simplexTimerStart(IterateClock);
  for (;;) {
    analysis->simplexTimerStart(IterateDualRebuildClock);
    rebuild();
    analysis->simplexTimerStop(IterateDualRebuildClock);
    if (solve_phase == kSolvePhaseError) {
      model_status = HighsModelStatus::kSolveError;
      return;
    }
    if (solve_phase == kSolvePhaseUnknown) {
      // If backtracking, may change phase, so drop out
      analysis->simplexTimerStop(IterateClock);
      return;
    }
    if (ekk_instance_.bailout()) break;
    for (;;) {
      if (debugDualSimplex("Before iteration") ==
          HighsDebugStatus::kLogicalError) {
        solve_phase = kSolvePhaseError;
        return;
      }
      switch (info.simplex_strategy) {
        default:
        case kSimplexStrategyDualPlain:
          iterate();
          break;
        case kSimplexStrategyDualTasks:
          iterateTasks();
          break;
        case kSimplexStrategyDualMulti:
          iterateMulti();
          break;
      }
      if (ekk_instance_.bailout()) break;
      assert(solve_phase != kSolvePhaseTabooBasis);
      if (rebuild_reason) break;
    }
    if (ekk_instance_.solve_bailout_) break;
    // If the data are fresh from rebuild(), possibly break out of the
    // outer loop to see what's occurred
    //
    // Deciding whether to rebuild is now more complicated if
    // refactorization is being avoided, since
    // status.has_fresh_rebuild being true may not imply that there is
    // a fresh factorization
    //
    bool finished = status.has_fresh_rebuild &&
                    !ekk_instance_.rebuildRefactor(rebuild_reason);
    if (finished && ekk_instance_.tabooBadBasisChange()) {
      // A bad basis change has had to be made taboo without any other
      // basis changes having been performed from a fresh rebuild. In
      // other words, the only basis change that could be made is not
      // permitted, so no definitive statement about the LP can be
      // made.
      solve_phase = kSolvePhaseTabooBasis;
      return;
    }
    if (finished) break;
  }
  analysis->simplexTimerStop(IterateClock);
  // Possibly return due to bailing out, having now stopped
  // IterateClock
  if (ekk_instance_.solve_bailout_) return;

  // If bailing out, should have done so already
  assert(!ekk_instance_.solve_bailout_);
  // Assess outcome of dual phase 1
  if (row_out == kNoRowChosen) {
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kDetailed,
                "dual-phase-1-optimal\n");
    // Optimal in phase 1
    if (info.dual_objective_value == 0) {
      // Zero phase 1 objective so go to phase 2
      //
      // This is the usual way to exit phase 1. Although the test
      // looks ambitious, the dual objective is the sum of products of
      // primal and dual values for nonbasic variables. For dual
      // simplex phase 1, the primal bounds are set so that when the
      // dual value is feasible, the primal value is set to
      // zero. Otherwise the value is +1/-1 according to the required
      // sign of the dual, except for free variables, where the bounds
      // are [-1000, 1000].
      //
      // OK if costs are perturbed, since they remain perturbed in phase 2 until
      // the final clean-up
      solve_phase = kSolvePhase2;
    } else {
      // A nonzero dual objective value at an optimal solution of
      // phase 1 means that there may be dual infeasibilities. If the
      // objective value is very negative, then it's clear
      // enough. However, if it's small in magnitude, it could be the
      // sum of values, all of which are smaller than the dual
      // feasibility tolerance. Plus there may be cost perturbations
      // to remove before reliable conclusions on dual infeasibility
      // of the (scaled) LP being solved can be drawn.
      //
      // If assessPhase1Optimality() is happy that there are no dual
      // infeasibilities, it will set solve_phase = kSolvePhase2;
      assessPhase1Optimality();
    }
  } else if (rebuild_reason == kRebuildReasonChooseColumnFail ||
             rebuild_reason == kRebuildReasonExcessivePrimalValue) {
    // chooseColumn has failed or excessive primal values have been
    // created Behave as "Report strange issues" below
    solve_phase = kSolvePhaseError;
    // Solve error is opaque to users, so put in some logging
    if (rebuild_reason == kRebuildReasonChooseColumnFail) {
      highsLogUser(
          ekk_instance_.options_->log_options, HighsLogType::kError,
          "Dual simplex ratio test failed due to excessive dual values: "
          "consider scaling down the LP objective coefficients\n");
    } else {
      highsLogUser(ekk_instance_.options_->log_options, HighsLogType::kError,
                   "Dual simplex detected excessive primal values: consider "
                   "scaling down the LP bounds\n");
    }
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "dual-phase-1-not-solved\n");
    model_status = HighsModelStatus::kSolveError;
  } else if (variable_in == -1) {
    // We got dual phase 1 unbounded - strange
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "dual-phase-1-unbounded\n");
    // This use of costs_alt_perturbed is (unsurprisingly) not
    // executed by ctest
    //
    // assert(99==1);
    if (ekk_instance_.info_.costs_perturbed) {
      // Clean up perturbation
      cleanup();
      highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kWarning,
                  "Cleaning up cost perturbation when unbounded in phase 1\n");
      if (dualInfeasCount == 0) {
        // No dual infeasibilities and (since unbounded) at least zero
        // phase 1 objective so go to phase 2
        solve_phase = kSolvePhase2;
      }
    } else {
      // Report strange issues
      solve_phase = kSolvePhaseError;
      highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                  "dual-phase-1-not-solved\n");
      model_status = HighsModelStatus::kSolveError;
    }
  }

  // Debug here since not all simplex data will be correct until after
  // rebuild() when switching to Phase 2
  //
  // Also have to avoid debug when the model status is not set and
  // there are dual infeasibilities, since this happens legitimately
  // when the LP is dual infeasible. However, the model status can't
  // be set to dual infeasible until perturbations have been removed.
  //
  const bool no_debug = ekk_instance_.info_.num_dual_infeasibilities > 0 &&
                        model_status == HighsModelStatus::kNotset;
  if (!no_debug) {
    if (debugDualSimplex("End of solvePhase1") ==
        HighsDebugStatus::kLogicalError) {
      solve_phase = kSolvePhaseError;
      return;
    }
  }

  // Can also get here with solve_phase = kSolvePhaseError
  const bool solve_phase_ok =
      solve_phase == kSolvePhase1 || solve_phase == kSolvePhase2 ||
      solve_phase == kSolvePhaseExit || solve_phase == kSolvePhaseError;
  if (!solve_phase_ok)
    highsLogDev(
        ekk_instance_.options_->log_options, HighsLogType::kInfo,
        "HEkkDual::solvePhase1 solve_phase == %d (solve call %d; iter %d)\n",
        (int)solve_phase, (int)ekk_instance_.debug_solve_call_num_,
        (int)ekk_instance_.iteration_count_);
  assert(solve_phase_ok);
  if (solve_phase == kSolvePhase2 || solve_phase == kSolvePhaseExit ||
      solve_phase == kSolvePhaseError) {
    // Moving to phase 2 or exiting, so make sure that the simplex
    // bounds and nonbasic value/move correspond to the LP
    ekk_instance_.initialiseBound(SimplexAlgorithm::kDual, kSolvePhase2);
    ekk_instance_.initialiseNonbasicValueAndMove();
    if (solve_phase == kSolvePhase2) {
      // Moving to phase 2 so possibly reinstate cost perturbation
      if (ekk_instance_.dual_simplex_phase1_cleanup_level_ <
          ekk_instance_.options_->max_dual_simplex_phase1_cleanup_level) {
        // Allow cost perturbation for now, but may have to prevent it
        // to avoid cleanup-perturbation loops
        info.allow_cost_shifting = true;
        info.allow_cost_perturbation = true;
      }
      // Comment if cost perturbation is not permitted
      if (!info.allow_cost_perturbation)
        highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kWarning,
                    "Moving to phase 2, but not allowing cost perturbation\n");
    }
  }
  return;
}

void HEkkDual::solvePhase2() {
  // Performs dual phase 2 iterations. Returns solve_phase with value
  //
  // kSolvePhaseError => Solver error
  //
  // kSolvePhaseTabooBasis => Only basis change is taboo
  //
  // kSolvePhaseExit => LP identified as not having an optimal solution
  //
  // kSolvePhaseUnknown => Back-tracking due to singularity
  //
  // kSolvePhaseOptimal => Primal feasible and no dual infeasibilities =>
  // Optimal
  //
  // kSolvePhase1 => Primal feasible and dual infeasibilities.
  //
  // kSolvePhase2 => Dual unboundedness suspected, but have to go out
  // and back in to solvePhase2 to perform fresh rebuild. Also if
  // bailing out due to reaching time/iteration limit or dual
  // objective
  //
  // kSolvePhaseOptimalCleanup => Continue with primal phase 2 iterations to
  // clean up dual infeasibilities
  //
  // Have to set multi_chooseAgain = 1 since it may come in set to 0
  // from the end of Phase 1, implying that there is a set of
  // attractive candidates to choose from in minorChooseRow() so
  // majorChooseRow() is unnecessary. If there are no such candidates
  // optimality may be declared in error. Previously, the forced
  // reinversion in the rebuild at the start of phase 2 led to
  // update_count == 0 causing multi_chooseAgain to be set to 1 in
  // majorChooseRow!
  multi_chooseAgain = 1;
  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsSimplexStatus& status = ekk_instance_.status_;
  HighsModelStatus& model_status = ekk_instance_.model_status_;
  // When starting a new phase the (updated) dual objective function
  // value isn't known. Indicate this so that when the value computed
  // from scratch in build() isn't checked against the updated
  // value
  status.has_primal_objective_value = false;
  status.has_dual_objective_value = false;
  // Set rebuild_reason so that it's assigned when first tested
  rebuild_reason = kRebuildReasonNo;
  // Set solve_phase = kSolvePhase2 and ekk_instance_.solve_bailout_ = false so
  // they are set if solvePhase2() is called directly
  solve_phase = kSolvePhase2;
  ekk_instance_.solve_bailout_ = false;
  if (ekk_instance_.bailout()) return;
  // Report the phase start
  highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kDetailed,
              "dual-phase-2-start\n");
  // Collect free variables
  dualRow.createFreelist();

  // If there's no backtracking basis, save the initial basis in case
  // of backtracking
  if (!info.valid_backtracking_basis_) ekk_instance_.putBacktrackingBasis();

  // Main solving structure
  analysis->simplexTimerStart(IterateClock);
  for (;;) {
    // Outer loop of solvePhase2()
    // Rebuild all values, reinverting B if updates have been performed
    analysis->simplexTimerStart(IterateDualRebuildClock);
    rebuild();
    analysis->simplexTimerStop(IterateDualRebuildClock);
    if (solve_phase == kSolvePhaseError) {
      model_status = HighsModelStatus::kSolveError;
      return;
    }
    if (solve_phase == kSolvePhaseUnknown) {
      // If backtracking, may change phase, so drop out
      analysis->simplexTimerStop(IterateClock);
      return;
    }
    if (ekk_instance_.bailout()) break;
    if (bailoutOnDualObjective()) break;
    if (dualInfeasCount > 0) break;
    for (;;) {
      // Inner loop of solvePhase2()
      // Performs one iteration in case kSimplexStrategyDualPlain:
      if (debugDualSimplex("Before iteration") ==
          HighsDebugStatus::kLogicalError) {
        solve_phase = kSolvePhaseError;
        return;
      }
      switch (info.simplex_strategy) {
        default:
        case kSimplexStrategyDualPlain:
          iterate();
          break;
        case kSimplexStrategyDualTasks:
          iterateTasks();
          break;
        case kSimplexStrategyDualMulti:
          iterateMulti();
          break;
      }
      if (ekk_instance_.bailout()) break;
      if (bailoutOnDualObjective()) break;
      assert(solve_phase != kSolvePhaseTabooBasis);

      // If possibly dual unbounded, assess whether this implies
      // primal infeasibility.
      if (rebuild_reason == kRebuildReasonPossiblyDualUnbounded)
        assessPossiblyDualUnbounded();

      if (rebuild_reason) break;
    }
    if (ekk_instance_.solve_bailout_) break;
    // If the data are fresh from rebuild(), possibly break out of the
    // outer loop to see what's occurred
    bool finished = status.has_fresh_rebuild &&
                    !ekk_instance_.rebuildRefactor(rebuild_reason);
    // ToDo: Handle the following more elegantly as the first case of
    // "Assess outcome of dual phase 2"
    if (finished && ekk_instance_.tabooBadBasisChange()) {
      // A bad basis change has had to be made taboo without any other
      // basis changes having been performed from a fresh rebuild. In
      // other words, the only basis change that could be made is not
      // permitted, so no definitive statement about the LP can be
      // made.
      solve_phase = kSolvePhaseTabooBasis;
      return;
    }
    if (finished) break;
  }
  analysis->simplexTimerStop(IterateClock);
  // Possibly return due to bailing out, having now stopped
  // IterateClock
  if (ekk_instance_.solve_bailout_) return;

  // If bailing out, should have done so already
  assert(!ekk_instance_.solve_bailout_);
  // Assess outcome of dual phase 2
  if (dualInfeasCount > 0) {
    // There are dual infeasibilities so possibly switch to Phase 1 and
    // return. "Possibly" because, if dual infeasibility has already
    // been shown, primal simplex is used to distinguish primal
    // unboundedness from primal infeasibility
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kDetailed,
                "dual-phase-2-found-free\n");
    solve_phase = kSolvePhase1;
  } else if (row_out == kNoRowChosen) {
    // There is no candidate in CHUZR, even after rebuild so probably optimal
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kDetailed,
                "dual-phase-2-optimal\n");
    // Remove any cost perturbations and see if basis is still dual feasible
    cleanup();
    if (dualInfeasCount > 0) {
      // There are dual infeasibilities, so consider performing primal
      // simplex iterations to get dual feasibility
      solve_phase = kSolvePhaseOptimalCleanup;
    } else {
      // There are no dual infeasibilities so optimal!
      solve_phase = kSolvePhaseOptimal;
      highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kDetailed,
                  "problem-optimal\n");
      model_status = HighsModelStatus::kOptimal;
    }
  } else if (rebuild_reason == kRebuildReasonChooseColumnFail ||
             rebuild_reason == kRebuildReasonExcessivePrimalValue) {
    // chooseColumn has failed or excessive primal values have been
    // created Behave as "Report strange issues" below
    solve_phase = kSolvePhaseError;
    // Solve error is opaque to users, so put in some logging
    if (rebuild_reason == kRebuildReasonChooseColumnFail) {
      highsLogUser(
          ekk_instance_.options_->log_options, HighsLogType::kError,
          "Dual simplex ratio test failed due to excessive dual values: "
          "consider scaling down the LP objective coefficients\n");
    } else {
      highsLogUser(ekk_instance_.options_->log_options, HighsLogType::kError,
                   "Dual simplex detected excessive primal values: consider "
                   "scaling down the LP bounds\n");
    }
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "dual-phase-2-not-solved\n");
    model_status = HighsModelStatus::kSolveError;
  } else {
    // Can only be that primal infeasibility has been detected
    assert(model_status == HighsModelStatus::kInfeasible);
    assert(solve_phase == kSolvePhaseExit);
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "problem-primal-infeasible\n");
  }
  // Possibly debug unless before primal simplex clean-up (in which
  // case there will be dual infeasibilities).
  if (solve_phase != kSolvePhaseOptimalCleanup) {
    if (debugDualSimplex("End of solvePhase2") ==
        HighsDebugStatus::kLogicalError) {
      solve_phase = kSolvePhaseError;
      return;
    }
  }
  return;
}

void HEkkDual::rebuild() {
  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsSimplexStatus& status = ekk_instance_.status_;
  // Clear taboo flag from any bad basis changes
  ekk_instance_.clearBadBasisChangeTabooFlag();

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

  HighsInt alt_debug_level = -1;
  //  if (ekk_instance_.debug_solve_report_) alt_debug_level =
  //  kHighsDebugLevelExpensive;
  ekk_instance_.debugNlaCheckInvert("HEkkDual::rebuild", alt_debug_level);

  if (!ekk_instance_.status_.has_ar_matrix) {
    // Don't have the row-wise matrix, so reinitialise it
    //
    // Should only happen when backtracking
    assert(info.backtracking_);
    ekk_instance_.initialisePartitionedRowwiseMatrix();
  }
  assert(ekk_instance_.ar_matrix_.debugPartitionOk(
      ekk_instance_.basis_.nonbasicFlag_.data()));

  // Record whether the update objective value should be tested. If
  // the objective value is known, then the updated objective value
  // should be correct - once the correction due to recomputing the
  // dual values has been applied.
  //
  // Note that computePrimalObjectiveValue sets
  // has_primal_objective_value
  const bool check_updated_objective_value = status.has_dual_objective_value;
  double previous_dual_objective_value = -kHighsInf;
  if (check_updated_objective_value) {
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase,
    //    "Before computeDual");
    previous_dual_objective_value = info.updated_dual_objective_value;
  } else {
    // Reset the knowledge of previous objective values
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, -1, "");
  }
  // Recompute dual solution
  ekk_instance_.computeDual();

  if (info.backtracking_) {
    // If backtracking, may change phase, so drop out
    solve_phase = kSolvePhaseUnknown;
    return;
  }
  analysis->simplexTimerStart(CorrectDualClock);
  correctDualInfeasibilities(dualInfeasCount);
  analysis->simplexTimerStop(CorrectDualClock);

  // Recompute primal solution
  ekk_instance_.computePrimal();

  // Collect primal infeasible as a list
  analysis->simplexTimerStart(CollectPrIfsClock);
  dualRHS.createArrayOfPrimalInfeasibilities();
  dualRHS.createInfeasList(ekk_instance_.info_.col_aq_density);
  analysis->simplexTimerStop(CollectPrIfsClock);

  // Dual objective section
  //
  ekk_instance_.computeDualObjectiveValue(solve_phase);

  if (check_updated_objective_value) {
    // Apply the objective value correction due to computing duals
    // from scratch.
    const double dual_objective_value_correction =
        info.dual_objective_value - previous_dual_objective_value;
    info.updated_dual_objective_value += dual_objective_value_correction;
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm);
  }
  // Now that there's a new dual_objective_value, reset the updated
  // value
  info.updated_dual_objective_value = info.dual_objective_value;

  if (!info.run_quiet) {
    ekk_instance_.computeInfeasibilitiesForReporting(SimplexAlgorithm::kDual,
                                                     solve_phase);
    reportRebuild(local_rebuild_reason);
  }

  // Record the synthetic clock for INVERT, and zero it for UPDATE
  ekk_instance_.resetSyntheticClock();

  // Dual simplex doesn't maintain the number of primal
  // infeasibilities, so set it to an illegal value now
  ekk_instance_.invalidatePrimalInfeasibilityRecord();
  // Although dual simplex should always be dual feasible,
  // infeasibilities are only corrected in rebuild
  ekk_instance_.invalidateDualInfeasibilityRecord();

  // Data are fresh from rebuild
  status.has_fresh_rebuild = true;
}

void HEkkDual::cleanup() {
  if (solve_phase == kSolvePhase1) {
    // Take action to prevent infinite loop. Shouldn't get to this
    // stage, since cost perturbation isn't permitted unless the dual
    // simplex phase 1 cleanup level is less than the limit
    ekk_instance_.dual_simplex_phase1_cleanup_level_++;
    const bool excessive_cleanup_calls =
        ekk_instance_.dual_simplex_phase1_cleanup_level_ >
        ekk_instance_.options_->max_dual_simplex_phase1_cleanup_level;
    if (excessive_cleanup_calls) {
      highsLogDev(
          ekk_instance_.options_->log_options, HighsLogType::kError,
          "Dual simplex cleanup level has exceeded limit of %d\n",
          (int)ekk_instance_.options_->max_dual_simplex_phase1_cleanup_level);
      assert(!excessive_cleanup_calls);
    }
  }
  highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kDetailed,
              "dual-cleanup-shift\n");
  HighsSimplexInfo& info = ekk_instance_.info_;
  // Remove perturbation and don't permit further perturbation
  ekk_instance_.initialiseCost(SimplexAlgorithm::kDual, kSolvePhaseUnknown);
  info.allow_cost_perturbation = false;
  // No solve_phase term in initialiseBound is surely an omission -
  // when cleanup called in phase 1
  ekk_instance_.initialiseBound(SimplexAlgorithm::kDual, solve_phase);
  // Possibly take a copy of the original duals before recomputing them
  vector<double> original_workDual;
  if (ekk_instance_.options_->highs_debug_level > kHighsDebugLevelCheap)
    original_workDual = info.workDual_;
  // Compute the dual values
  ekk_instance_.computeDual();
  // Possibly analyse the change in duals
  //  debugCleanup(ekk_instance_, original_workDual);
  // Compute the dual infeasibilities
  ekk_instance_.computeSimplexDualInfeasible();
  dualInfeasCount = ekk_instance_.info_.num_dual_infeasibilities;

  // Compute the dual objective value
  ekk_instance_.computeDualObjectiveValue(solve_phase);
  // Now that there's a new dual_objective_value, reset the updated
  // value
  info.updated_dual_objective_value = info.dual_objective_value;

  if (!info.run_quiet) {
    // Report the primal infeasibilities
    ekk_instance_.computeSimplexPrimalInfeasible();
    // In phase 1, report the simplex LP dual infeasibilities
    // In phase 2, report the simplex dual infeasibilities (known)
    if (solve_phase == kSolvePhase1)
      ekk_instance_.computeSimplexLpDualInfeasible();
    reportRebuild(kRebuildReasonCleanup);
  }
}

void HEkkDual::iterate() {
  // This is the main iteration loop for dual revised simplex. All the
  // methods have as their first line if (rebuild_reason) return;, where
  // rebuild_reason is, for example, set to 1 when CHUZR finds no
  // candidate. This causes a break from the inner loop of
  // solvePhase% and, hence, a call to rebuild()

  // Reporting:
  // Row-wise matrix after update in updateMatrix(variable_in, variable_out);

  const HighsInt from_check_iter = 0;
  const HighsInt to_check_iter = from_check_iter + 100;
  if (ekk_instance_.debug_solve_report_) {
    ekk_instance_.debug_iteration_report_ =
        ekk_instance_.iteration_count_ >= from_check_iter &&
        ekk_instance_.iteration_count_ <= to_check_iter;
    if (ekk_instance_.debug_iteration_report_) {
      printf("HEkkDual::iterate Debug iteration %d\n",
             (int)ekk_instance_.iteration_count_);
    }
  }

  analysis->simplexTimerStart(IterateChuzrClock);
  chooseRow();
  analysis->simplexTimerStop(IterateChuzrClock);

  analysis->simplexTimerStart(IterateChuzcClock);
  chooseColumn(&row_ep);
  analysis->simplexTimerStop(IterateChuzcClock);

  if (isBadBasisChange()) return;

  analysis->simplexTimerStart(IterateFtranClock);
  updateFtranBFRT();

  // updateFtran(); computes the pivotal column in the data structure "column"
  updateFtran();

  // updateFtranDSE performs the DSE FTRAN on pi_p
  if (edge_weight_mode == EdgeWeightMode::kSteepestEdge)
    updateFtranDSE(&row_ep);
  analysis->simplexTimerStop(IterateFtranClock);

  // updateVerify() Checks row-wise pivot against column-wise pivot for
  // numerical trouble
  analysis->simplexTimerStart(IterateVerifyClock);
  updateVerify();
  analysis->simplexTimerStop(IterateVerifyClock);

  // updateDual() Updates the dual values
  analysis->simplexTimerStart(IterateDualClock);
  updateDual();
  analysis->simplexTimerStop(IterateDualClock);

  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase, "Before
  //  updatePrimal");
  // updatePrimal(&row_ep); Updates the primal values and the edge weights
  analysis->simplexTimerStart(IteratePrimalClock);
  updatePrimal(&row_ep);
  analysis->simplexTimerStop(IteratePrimalClock);
  // After primal update in dual simplex the primal objective value is not known
  ekk_instance_.status_.has_primal_objective_value = false;
  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase, "After
  //  updatePrimal");

  // Update the records of chosen rows and pivots
  //  ekk_instance_.info_.pivot_.push_back(alpha_row);
  //  ekk_instance_.info_.index_chosen_.push_back(row_out);

  // Update the basis representation
  analysis->simplexTimerStart(IteratePivotsClock);
  updatePivots();
  analysis->simplexTimerStop(IteratePivotsClock);

  if (new_devex_framework) {
    // Initialise new Devex framework
    analysis->simplexTimerStart(IterateDevexIzClock);
    initialiseDevexFramework();
    analysis->simplexTimerStop(IterateDevexIzClock);
  }

  // Analyse the iteration: possibly report; possibly switch strategy
  iterationAnalysis();
}

void HEkkDual::iterateTasks() {
  slice_PRICE = 1;

  // Group 1
  chooseRow();

  // Disable slice when too sparse
  if (1.0 * row_ep.count * inv_solver_num_row < 0.01) slice_PRICE = 0;

  analysis->simplexTimerStart(Group1Clock);
  // #pragma omp parallel
  // #pragma omp single
  {
    // #pragma omp task
    highs::parallel::spawn([&]() {
      col_DSE.copy(&row_ep);
      updateFtranDSE(&col_DSE);
    });
    // #pragma omp task
    {
      if (slice_PRICE)
        chooseColumnSlice(&row_ep);
      else
        chooseColumn(&row_ep);
      // #pragma omp task
      highs::parallel::spawn([&]() { updateFtranBFRT(); });
      // #pragma omp task
      updateFtran();
      // #pragma omp taskwait
      highs::parallel::sync();
    }

    highs::parallel::sync();
  }
  analysis->simplexTimerStop(Group1Clock);

  updateVerify();
  updateDual();
  updatePrimal(&col_DSE);
  updatePivots();
}

void HEkkDual::iterationAnalysisData() {
  double cost_scale_factor =
      pow(2.0, -ekk_instance_.options_->cost_scale_factor);
  HighsSimplexInfo& info = ekk_instance_.info_;
  analysis->simplex_strategy = info.simplex_strategy;
  analysis->edge_weight_mode = edge_weight_mode;
  analysis->solve_phase = solve_phase;
  analysis->simplex_iteration_count = ekk_instance_.iteration_count_;
  analysis->devex_iteration_count = num_devex_iterations;
  analysis->pivotal_row_index = row_out;
  analysis->leaving_variable = variable_out;
  analysis->entering_variable = variable_in;
  analysis->rebuild_reason = rebuild_reason;
  analysis->reduced_rhs_value = 0;
  analysis->reduced_cost_value = 0;
  analysis->edge_weight = 0;
  analysis->primal_delta = delta_primal;
  analysis->primal_step = theta_primal;
  analysis->dual_step = theta_dual * cost_scale_factor;
  analysis->pivot_value_from_column = alpha_col;
  analysis->pivot_value_from_row = alpha_row;
  analysis->factor_pivot_threshold = info.factor_pivot_threshold;
  analysis->numerical_trouble = numericalTrouble;
  analysis->edge_weight_error = ekk_instance_.edge_weight_error_;
  analysis->objective_value = info.updated_dual_objective_value;
  // Since maximization is achieved by minimizing the LP with negated
  // costs, in phase 2 the dual objective value is negated, so flip
  // its sign according to the LP sense
  if (solve_phase == kSolvePhase2)
    analysis->objective_value *= (HighsInt)ekk_instance_.lp_.sense_;
  analysis->num_primal_infeasibility = info.num_primal_infeasibilities;
  analysis->sum_primal_infeasibility = info.sum_primal_infeasibilities;
  if (solve_phase == kSolvePhase1) {
    analysis->num_dual_infeasibility =
        analysis->num_dual_phase_1_lp_dual_infeasibility;
    analysis->sum_dual_infeasibility =
        analysis->sum_dual_phase_1_lp_dual_infeasibility;
  } else {
    analysis->num_dual_infeasibility = info.num_dual_infeasibilities;
    analysis->sum_dual_infeasibility = info.sum_dual_infeasibilities;
  }
  if ((edge_weight_mode == EdgeWeightMode::kDevex) &&
      (num_devex_iterations == 0))
    analysis->num_devex_framework++;
  analysis->col_aq_density = info.col_aq_density;
  analysis->row_ep_density = info.row_ep_density;
  analysis->row_ap_density = info.row_ap_density;
  analysis->row_DSE_density = info.row_DSE_density;
  analysis->col_basic_feasibility_change_density =
      info.col_basic_feasibility_change_density;
  analysis->row_basic_feasibility_change_density =
      info.row_basic_feasibility_change_density;
  analysis->col_BFRT_density = info.col_BFRT_density;
  analysis->primal_col_density = info.primal_col_density;
  analysis->dual_col_density = info.dual_col_density;
  analysis->num_costly_DSE_iteration = info.num_costly_DSE_iteration;
  analysis->costly_DSE_measure = info.costly_DSE_measure;
  //  analysis-> = info.;
}

void HEkkDual::iterationAnalysis() {
  // Compute the infeasibility data (expensive) if analysing run-time
  // data and the log level is at least kIterationReportLogType
  // (Verbose)
  const bool make_iteration_report = analysis->analyse_simplex_runtime_data &&
                                     ekk_instance_.options_->log_dev_level >=
                                         (HighsInt)kIterationReportLogType;
  if (make_iteration_report)
    ekk_instance_.computeInfeasibilitiesForReporting(SimplexAlgorithm::kDual,
                                                     solve_phase);
  // Possibly report on the iteration
  iterationAnalysisData();
  analysis->iterationReport();

  // Possibly switch from DSE to Devex
  if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
    const bool switch_to_devex = ekk_instance_.switchToDevex();
    if (switch_to_devex) {
      edge_weight_mode = EdgeWeightMode::kDevex;
      initialiseDevexFramework();
    }
  }
  if (analysis->analyse_simplex_summary_data) analysis->iterationRecord();
}

void HEkkDual::reportRebuild(const HighsInt reason_for_rebuild) {
  analysis->simplexTimerStart(ReportRebuildClock);
  iterationAnalysisData();
  analysis->rebuild_reason = reason_for_rebuild;
  analysis->rebuild_reason_string =
      ekk_instance_.rebuildReason(reason_for_rebuild);
  if (ekk_instance_.options_->output_flag) analysis->invertReport();
  analysis->simplexTimerStop(ReportRebuildClock);
}

void HEkkDual::chooseRow() {
  // Choose the index of a row to leave the basis (CHUZR)
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;
  //  if (solve_phase == kSolvePhase2) dualRHS.assessOptimality();
  // Zero the infeasibility of any taboo rows
  ekk_instance_.applyTabooRowOut(dualRHS.work_infeasibility, 0);
  // Choose candidates repeatedly until candidate is OK or optimality is
  // detected
  if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
    HighsDebugStatus return_status =
        ekk_instance_.devDebugDualSteepestEdgeWeights("chooseRow");
    assert(return_status == HighsDebugStatus::kNotChecked ||
           return_status == HighsDebugStatus::kOk);
  }
  std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
  for (;;) {
    // Choose the index of a good row to leave the basis
    dualRHS.chooseNormal(&row_out);
    if (row_out == kNoRowChosen) {
      // No index found so may be dual optimal.
      rebuild_reason = kRebuildReasonPossiblyOptimal;
      return;
    }
    // Compute pi_p = B^{-T}e_p in row_ep
    analysis->simplexTimerStart(BtranClock);
    // Set up RHS for BTRAN
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = row_out;
    row_ep.array[row_out] = 1;
    row_ep.packFlag = true;
    if (analysis->analyse_simplex_summary_data)
      analysis->operationRecordBefore(kSimplexNlaBtranEp, row_ep,
                                      ekk_instance_.info_.row_ep_density);
    // Perform BTRAN
    simplex_nla->btran(row_ep, ekk_instance_.info_.row_ep_density,
                       analysis->pointer_serial_factor_clocks);
    if (analysis->analyse_simplex_summary_data)
      analysis->operationRecordAfter(kSimplexNlaBtranEp, row_ep);
    analysis->simplexTimerStop(BtranClock);
    // Verify DSE weight
    if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
      // For DSE, see how accurate the updated weight is
      // Save the updated weight
      double updated_edge_weight = edge_weight[row_out];
      // Compute the weight from row_ep and over-write the updated weight
      if (ekk_instance_.simplex_in_scaled_space_) {
        computed_edge_weight = edge_weight[row_out] = row_ep.norm2();
      } else {
        computed_edge_weight = edge_weight[row_out] =
            simplex_nla->rowEp2NormInScaledSpace(row_out, row_ep);
      }
      // If the weight error is acceptable then break out of the
      // loop. All we worry about is accepting rows with weights
      // which are not too small, since this can make the row look
      // unreasonably attractive
      if (acceptDualSteepestEdgeWeight(updated_edge_weight)) break;
      // Weight error is unacceptable so look for another
      // candidate. Of course, it's possible that the same
      // candidate is chosen, but the weight will be correct (so
      // no infinite loop).
    } else {
      // If not using DSE then accept the row by breaking out of
      // the loop
      break;
    }
  }
  // Recover the infeasibility of any taboo rows
  ekk_instance_.unapplyTabooRowOut(dualRHS.work_infeasibility);

  // Index of row to leave the basis has been found
  //
  // Assign basic info:
  //
  // Record the column (variable) associated with the leaving row
  variable_out = ekk_instance_.basis_.basicIndex_[row_out];
  // Record the change in primal variable associated with the move to the bound
  // being violated
  if (baseValue[row_out] < baseLower[row_out]) {
    // Below the lower bound so set delta_primal = value - LB < 0
    delta_primal = baseValue[row_out] - baseLower[row_out];
  } else {
    // Above the upper bound so set delta_primal = value - UB > 0
    delta_primal = baseValue[row_out] - baseUpper[row_out];
  }
  // Set move_out to be -1 if delta_primal<0, otherwise +1 (since
  // delta_primal>0)
  move_out = delta_primal < 0 ? -1 : 1;
  // Update the record of average row_ep (pi_p) density. This ignores
  // any BTRANs done for skipped candidates
  const double local_row_ep_density = (double)row_ep.count * inv_solver_num_row;
  ekk_instance_.updateOperationResultDensity(
      local_row_ep_density, ekk_instance_.info_.row_ep_density);
}

bool HEkkDual::acceptDualSteepestEdgeWeight(const double updated_edge_weight) {
  // Accept the updated weight if it is at least a quarter of the
  // computed weight. Excessively large updated weights don't matter!
  const bool accept_weight =
      updated_edge_weight >= kAcceptDseWeightThreshold * computed_edge_weight;
  //  if (analysis->analyse_simplex_summary_data)
  ekk_instance_.assessDSEWeightError(computed_edge_weight, updated_edge_weight);
  analysis->dualSteepestEdgeWeightError(computed_edge_weight,
                                        updated_edge_weight);
  return accept_weight;
}

bool HEkkDual::newDevexFramework(const double updated_edge_weight) {
  // Analyse the Devex weight to determine whether a new framework
  // should be set up
  //
  // There is a new Devex framework if either
  //
  // 1) The weight inaccuracy ratio exceeds kMaxAllowedDevexWeightRatio
  //
  // 2) There have been max(kMinAbsNumberDevexIterations,
  // numRow/kMinRlvNumberDevexIterations) Devex iterations
  //
  const HighsInt kMinAbsNumberDevexIterations = 25;
  const double kMinRlvNumberDevexIterations = 1e-2;
  const double kMaxAllowedDevexWeightRatio = 3.0;

  double devex_ratio = max(updated_edge_weight / computed_edge_weight,
                           computed_edge_weight / updated_edge_weight);
  HighsInt i_te = solver_num_row / kMinRlvNumberDevexIterations;
  i_te = max(kMinAbsNumberDevexIterations, i_te);
  // Square kMaxAllowedDevexWeightRatio due to keeping squared
  // weights
  const double accept_ratio_threshold =
      kMaxAllowedDevexWeightRatio * kMaxAllowedDevexWeightRatio;
  const bool accept_ratio = devex_ratio <= accept_ratio_threshold;
  const bool accept_it = num_devex_iterations <= i_te;
  bool return_new_devex_framework;
  return_new_devex_framework = !accept_ratio || !accept_it;
  return return_new_devex_framework;
}

void HEkkDual::chooseColumn(HVector* row_ep) {
  // Compute pivot row (PRICE) and choose the index of a column to enter the
  // basis (CHUZC)
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;
  HighsOptions* options = ekk_instance_.options_;
  HighsLp& lp = ekk_instance_.lp_;

  HighsInt debug_price_report = kDebugReportOff;
  const bool debug_price_report_on = false;
  const bool debug_small_pivot_issue_report_on = false;
  bool debug_small_pivot_issue_report = false;
  bool debug_rows_report = false;
  if (ekk_instance_.debug_iteration_report_) {
    if (debug_price_report_on) debug_price_report = kDebugReportAll;
    debug_rows_report = debug_price_report_on;
    debug_small_pivot_issue_report = debug_small_pivot_issue_report_on;
    if (debug_price_report != kDebugReportOff || debug_rows_report ||
        debug_small_pivot_issue_report)
      printf("HEkkDual::chooseColumn Check iter = %d\n",
             (int)ekk_instance_.iteration_count_);
  }
  //
  // PRICE
  //
  const bool quad_precision = false;
  ekk_instance_.tableauRowPrice(quad_precision, *row_ep, row_ap,
                                debug_price_report);
  if (debug_rows_report) {
    ekk_instance_.simplex_nla_.reportArray("Row a_p", 0, &row_ap, true);
    ekk_instance_.simplex_nla_.reportArray("Row e_p", lp.num_col_, row_ep,
                                           true);
  }
  //
  // CHUZC
  //
  // Section 0: Clear data and call createFreemove to set a value of
  // nonbasicMove for all free columns to prevent their dual values
  // from being changed.
  analysis->simplexTimerStart(Chuzc0Clock);
  dualRow.clear();
  dualRow.workDelta = delta_primal;
  dualRow.createFreemove(row_ep);
  analysis->simplexTimerStop(Chuzc0Clock);
  //
  // Section 1: Pack row_ap and row_ep
  analysis->simplexTimerStart(Chuzc1Clock);
  // Pack row_ap into the packIndex/Value of HEkkDualRow
  dualRow.chooseMakepack(&row_ap, 0);
  // Pack row_ep into the packIndex/Value of HEkkDualRow
  dualRow.chooseMakepack(row_ep, solver_num_col);
  const double row_ep_scale =
      ekk_instance_.getValueScale(dualRow.packCount, dualRow.packValue);
  analysis->simplexTimerStop(Chuzc1Clock);
  // Loop until an acceptable pivot is found. Each pass either finds a
  // pivot, identifies possible unboundedness, or reduced the number
  // of nonzeros in dualRow.pack_value
  HighsInt chuzc_pass = 0;
  for (;;) {
    //
    // Section 2: Determine the possible variables - candidates for CHUZC
    analysis->simplexTimerStart(Chuzc2Clock);
    dualRow.choosePossible();
    analysis->simplexTimerStop(Chuzc2Clock);
    //
    // Take action if the step to an expanded bound is not positive, or
    // there are no candidates for CHUZC
    variable_in = -1;
    if (dualRow.workTheta <= 0 || dualRow.workCount == 0) {
      if (chuzc_pass > 0 && debug_small_pivot_issue_report) {
        printf(
            "                                                       "
            "Negative step or no candidates after %2d CHUZC passes\n",
            (int)chuzc_pass);
      }
      if (debug_rows_report) {
        ekk_instance_.simplex_nla_.reportVector(
            "dualRow.packValue/Index", dualRow.packCount, dualRow.packValue,
            dualRow.packIndex, true);
      }
      rebuild_reason = kRebuildReasonPossiblyDualUnbounded;
      return;
    }
    //
    // Sections 3 and 4: Perform (bound-flipping) ratio test. This can
    // fail if the dual values are excessively large
    bool chooseColumnFail = (dualRow.chooseFinal() != 0);
    if (chooseColumnFail) {
      rebuild_reason = kRebuildReasonChooseColumnFail;
      return;
    }
    if (dualRow.workPivot >= 0) {
      const double growth_tolerance =
          options->dual_simplex_pivot_growth_tolerance;  // kHighsMacheps; //
      // A pivot has been chosen
      double alpha_row = dualRow.workAlpha;
      assert(alpha_row);
      const double scaled_value = row_ep_scale * alpha_row;
      if (std::abs(scaled_value) <= growth_tolerance) {
        if (chuzc_pass == 0 && debug_small_pivot_issue_report)
          printf(
              "CHUZC: Solve %6d; Iter %4d; ||e_p|| = %11.4g: Variable %6d "
              "Pivot %11.4g (dual "
              "%11.4g; ratio = %11.4g) is small",
              (int)ekk_instance_.debug_solve_call_num_,
              (int)ekk_instance_.iteration_count_, 1.0 / row_ep_scale,
              (int)dualRow.workPivot, dualRow.workAlpha,
              workDual[dualRow.workPivot],
              workDual[dualRow.workPivot] / dualRow.workAlpha);
        // On the first pass, try to make the pivotal row more accurate
        if (chuzc_pass == 0) {
          if (debug_small_pivot_issue_report) printf(": improve row\n");
          ekk_instance_.analysis_.num_improve_choose_column_row_call++;
          improveChooseColumnRow(row_ep);
        } else {
          // Remove the pivot
          ekk_instance_.analysis_.num_remove_pivot_from_pack++;
          for (HighsInt i = 0; i < dualRow.packCount; i++) {
            if (dualRow.packIndex[i] == dualRow.workPivot) {
              dualRow.packIndex[i] = dualRow.packIndex[dualRow.packCount - 1];
              dualRow.packValue[i] = dualRow.packValue[dualRow.packCount - 1];
              dualRow.packCount--;
              if (chuzc_pass == 0 && debug_small_pivot_issue_report)
                printf(": removing pivot gives pack count = %6d",
                       (int)dualRow.packCount);
              break;
            }
          }
          if (chuzc_pass == 0 && debug_small_pivot_issue_report)
            printf(" so %s\n", dualRow.packCount > 0 ? "repeat CHUZC"
                                                     : "consider unbounded");
        }
        // Indicate that no pivot has been chosen
        dualRow.workPivot = -1;
      } else if (chuzc_pass > 0 && debug_small_pivot_issue_report) {
        printf(
            "                                                       Variable "
            "%6d "
            "Pivot %11.4g (dual "
            "%11.4g; ratio = %11.4g) is OK after %2d CHUZC passes\n",
            (int)dualRow.workPivot, dualRow.workAlpha,
            workDual[dualRow.workPivot],
            workDual[dualRow.workPivot] / dualRow.workAlpha, (int)chuzc_pass);
      }
    } else {
      // No pivot has been chosen
      assert(dualRow.workPivot == -1);
      if (chuzc_pass > 0 && debug_small_pivot_issue_report) {
        printf(
            "                                                       No pivot "
            "after %2d CHUZC passes\n",
            (int)chuzc_pass);
      }
      break;
    }
    // If a pivot has been chosen, or there are no more packed values
    // then end CHUZC
    if (dualRow.workPivot >= 0 || dualRow.packCount <= 0) break;
    chuzc_pass++;
  }
  //
  // Section 5: Reset the nonbasicMove values for free columns
  analysis->simplexTimerStart(Chuzc5Clock);
  dualRow.deleteFreemove();
  analysis->simplexTimerStop(Chuzc5Clock);
  // Record values for basis change, checking for numerical problems and update
  // of dual variables
  variable_in = dualRow.workPivot;  // Index of the column entering the basis
  alpha_row = dualRow.workAlpha;    // Pivot value computed row-wise - used for
                                    // numerical checking
  theta_dual = dualRow.workTheta;   // Dual step length

  if (edge_weight_mode == EdgeWeightMode::kDevex && !new_devex_framework) {
    // When using Devex, unless a new framework is to be used, get the
    // exact weight for the pivotal row and, based on its accuracy,
    // determine that a new framework is to be used. In serial
    // new_devex_framework should only ever be false at this point in
    // this method, but in PAMI, this method may be called multiple
    // times in minor iterations and the new framework is set up in
    // majorUpdate.
    analysis->simplexTimerStart(DevexWtClock);
    // Determine the exact Devex weight
    dualRow.computeDevexWeight();
    computed_edge_weight = dualRow.computed_edge_weight;
    computed_edge_weight = max(1.0, computed_edge_weight);
    analysis->simplexTimerStop(DevexWtClock);
  }
  return;
}

void HEkkDual::improveChooseColumnRow(HVector* row_ep) {
  HighsLp& lp = ekk_instance_.lp_;
  const bool debug_price_report_on = false;
  HighsInt debug_price_report = kDebugReportOff;
  bool debug_rows_report = false;
  if (ekk_instance_.debug_iteration_report_) {
    if (debug_price_report_on) debug_price_report = kDebugReportAll;
    debug_rows_report = debug_price_report_on;
    if (debug_price_report != kDebugReportOff)
      printf("HEkkDual::chooseColumn Check iter = %d\n",
             (int)ekk_instance_.iteration_count_);
  }

  analysis->simplexTimerStart(Chuzc5Clock);
  dualRow.deleteFreemove();
  analysis->simplexTimerStop(Chuzc5Clock);

  if (debug_rows_report)
    ekk_instance_.simplex_nla_.reportArray("Row e_p.0", lp.num_col_, row_ep,
                                           true);
  ekk_instance_.unitBtranIterativeRefinement(row_out, *row_ep);
  if (debug_rows_report)
    ekk_instance_.simplex_nla_.reportArray("Row e_p.1", lp.num_col_, row_ep,
                                           true);

  const bool quad_precision = true;
  ekk_instance_.tableauRowPrice(quad_precision, *row_ep, row_ap,
                                debug_price_report);
  if (debug_rows_report) {
    ekk_instance_.simplex_nla_.reportArray("Row a_p", 0, &row_ap, true);
    ekk_instance_.simplex_nla_.reportArray("Row e_p", lp.num_col_, row_ep,
                                           true);
  }
  analysis->simplexTimerStart(Chuzc0Clock);
  dualRow.clear();
  dualRow.workDelta = delta_primal;
  dualRow.createFreemove(row_ep);
  analysis->simplexTimerStop(Chuzc0Clock);
  //
  // Section 1: Pack row_ap and row_ep
  analysis->simplexTimerStart(Chuzc1Clock);
  // Pack row_ap into the packIndex/Value of HEkkDualRow
  dualRow.chooseMakepack(&row_ap, 0);
  // Pack row_ep into the packIndex/Value of HEkkDualRow
  dualRow.chooseMakepack(row_ep, solver_num_col);
  analysis->simplexTimerStop(Chuzc1Clock);
}

void HEkkDual::chooseColumnSlice(HVector* row_ep) {
  // Choose the index of a column to enter the basis (CHUZC) by
  // exploiting slices of the pivotal row - for SIP and PAMI
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;

  analysis->simplexTimerStart(Chuzc0Clock);
  dualRow.clear();
  dualRow.workDelta = delta_primal;
  dualRow.createFreemove(row_ep);
  analysis->simplexTimerStop(Chuzc0Clock);

  //  const HighsInt solver_num_row = ekk_instance_.lp_.num_row_;
  const double local_density = 1.0 * row_ep->count * inv_solver_num_row;
  bool use_col_price;
  bool use_row_price_w_switch;
  HighsSimplexInfo& info = ekk_instance_.info_;
  ekk_instance_.choosePriceTechnique(info.price_strategy, local_density,
                                     use_col_price, use_row_price_w_switch);

  if (analysis->analyse_simplex_summary_data) {
    const HighsInt row_ep_count = row_ep->count;
    if (use_col_price) {
      analysis->operationRecordBefore(kSimplexNlaPriceAp, row_ep_count, 0.0);
      analysis->num_col_price++;
    } else if (use_row_price_w_switch) {
      analysis->operationRecordBefore(kSimplexNlaPriceAp, row_ep_count,
                                      ekk_instance_.info_.row_ep_density);
      analysis->num_row_price_with_switch++;
    } else {
      analysis->operationRecordBefore(kSimplexNlaPriceAp, row_ep_count,
                                      ekk_instance_.info_.row_ep_density);
      analysis->num_row_price++;
    }
  }
  analysis->simplexTimerStart(PriceChuzc1Clock);
  // Row_ep:         PACK + CC1

  highs::parallel::spawn([&]() {
    dualRow.chooseMakepack(row_ep, solver_num_col);
    dualRow.choosePossible();
  });

  // Row_ap: PRICE + PACK + CC1
  highs::parallel::for_each(0, slice_num, [&](HighsInt start, HighsInt end) {
    const bool quad_precision = false;
    for (HighsInt i = start; i < end; i++) {
      slice_row_ap[i].clear();

      if (use_col_price) {
        // Perform column-wise PRICE
        slice_a_matrix[i].priceByColumn(quad_precision, slice_row_ap[i],
                                        *row_ep);
      } else if (use_row_price_w_switch) {
        // Perform hyper-sparse row-wise PRICE, but switch if the density of
        // row_ap becomes extreme
        slice_ar_matrix[i].priceByRowWithSwitch(
            quad_precision, slice_row_ap[i], *row_ep,
            ekk_instance_.info_.row_ap_density, 0, kHyperPriceDensity);
      } else {
        // Perform hyper-sparse row-wise PRICE
        slice_ar_matrix[i].priceByRow(quad_precision, slice_row_ap[i], *row_ep);
      }

      slice_dualRow[i].clear();
      slice_dualRow[i].workDelta = delta_primal;
      slice_dualRow[i].chooseMakepack(&slice_row_ap[i], slice_start[i]);
      slice_dualRow[i].choosePossible();
    }
  });

  highs::parallel::sync();

  if (analysis->analyse_simplex_summary_data) {
    // Determine the nonzero count of the whole row
    HighsInt row_ap_count = 0;
    for (HighsInt i = 0; i < slice_num; i++)
      row_ap_count += slice_row_ap[i].count;
    analysis->operationRecordAfter(kSimplexNlaPriceAp, row_ap_count);
  }

  // Join CC1 results here
  for (HighsInt i = 0; i < slice_num; i++) {
    dualRow.chooseJoinpack(&slice_dualRow[i]);
  }

  analysis->simplexTimerStop(PriceChuzc1Clock);

  // Infeasible we created before
  variable_in = -1;
  if (dualRow.workTheta <= 0 || dualRow.workCount == 0) {
    rebuild_reason = kRebuildReasonPossiblyDualUnbounded;
    return;
  }

  // Choose column 2, This only happens if didn't go out
  HighsInt return_code = dualRow.chooseFinal();
  if (return_code) {
    // Only returns -1, if not zero
    assert(return_code == -1);
    if (return_code < 0) {
      rebuild_reason = kRebuildReasonChooseColumnFail;
    } else {
      rebuild_reason = kRebuildReasonPossiblyDualUnbounded;
    }
    return;
  }

  if (slice_num == 0) {
    // This check is only done for serial code - since packIndex/Value
    // is distributed
    HighsInt num_infeasibility = dualRow.debugChooseColumnInfeasibilities();
    if (num_infeasibility) {
      highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kError,
                  "chooseFinal would create %d dual infeasibilities\n",
                  (int)num_infeasibility);
      analysis->simplexTimerStop(Chuzc4dClock);
      rebuild_reason = kRebuildReasonChooseColumnFail;
      return;
    }
  }

  analysis->simplexTimerStart(Chuzc5Clock);
  dualRow.deleteFreemove();
  analysis->simplexTimerStop(Chuzc5Clock);

  variable_in = dualRow.workPivot;
  alpha_row = dualRow.workAlpha;
  theta_dual = dualRow.workTheta;

  if (edge_weight_mode == EdgeWeightMode::kDevex && !new_devex_framework) {
    // When using Devex, unless a new framework is to be used, get the
    // exact weight for the pivotal row and, based on its accuracy,
    // determine that a new framework is to be used. In serial
    // new_devex_framework should only ever be false at this point in
    // this method, but in PAMI, this method may be called multiple
    // times in minor iterations and the new framework is set up in
    // majorUpdate.
    analysis->simplexTimerStart(DevexWtClock);
    // Determine the partial sums of the exact Devex weight
    // First the partial sum for row_ep
    dualRow.computeDevexWeight();
    // Second the partial sums for the slices of row_ap
    for (HighsInt i = 0; i < slice_num; i++)
      slice_dualRow[i].computeDevexWeight(i);
    // Accumulate the partial sums
    // Initialise with the partial sum for row_ep
    computed_edge_weight = dualRow.computed_edge_weight;
    // Update with the partial sum for row_ep
    for (HighsInt i = 0; i < slice_num; i++)
      computed_edge_weight += slice_dualRow[i].computed_edge_weight;
    computed_edge_weight = max(1.0, computed_edge_weight);
    analysis->simplexTimerStop(DevexWtClock);
  }
}

void HEkkDual::updateFtran() {
  // Compute the pivotal column (FTRAN)
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;
  analysis->simplexTimerStart(FtranClock);
  // Clear the pivotal column and indicate that its values should be packed
  col_aq.clear();
  col_aq.packFlag = true;
  // Get the constraint matrix column by combining just one column
  // with unit multiplier
  a_matrix->collectAj(col_aq, variable_in, 1);
  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordBefore(kSimplexNlaFtran, col_aq,
                                    ekk_instance_.info_.col_aq_density);
  // Perform FTRAN
  simplex_nla->ftran(col_aq, ekk_instance_.info_.col_aq_density,
                     analysis->pointer_serial_factor_clocks);
  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordAfter(kSimplexNlaFtran, col_aq);
  const double local_col_aq_density = (double)col_aq.count * inv_solver_num_row;
  ekk_instance_.updateOperationResultDensity(
      local_col_aq_density, ekk_instance_.info_.col_aq_density);
  // Save the pivot value computed column-wise - used for numerical checking
  alpha_col = col_aq.array[row_out];
  analysis->simplexTimerStop(FtranClock);
}

void HEkkDual::updateFtranBFRT() {
  // Compute the RHS changes corresponding to the BFRT (FTRAN-BFRT)
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;

  // Only time updateFtranBFRT if dualRow.workCount > 0;
  // If dualRow.workCount = 0 then dualRow.updateFlip(&col_BFRT)
  // merely clears col_BFRT so no FTRAN is performed
  bool time_updateFtranBFRT = dualRow.workCount > 0;

  if (time_updateFtranBFRT) {
    analysis->simplexTimerStart(FtranBfrtClock);
  }

  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase, "Before
  //  update_flip");
  dualRow.updateFlip(&col_BFRT);
  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase, "After
  //  update_flip");

  if (col_BFRT.count) {
    if (analysis->analyse_simplex_summary_data)
      analysis->operationRecordBefore(kSimplexNlaFtranBfrt, col_BFRT,
                                      ekk_instance_.info_.col_BFRT_density);
    simplex_nla->ftran(col_BFRT, ekk_instance_.info_.col_BFRT_density,
                       analysis->pointer_serial_factor_clocks);
    if (analysis->analyse_simplex_summary_data)
      analysis->operationRecordAfter(kSimplexNlaFtranBfrt, col_BFRT);
  }
  if (time_updateFtranBFRT) {
    analysis->simplexTimerStop(FtranBfrtClock);
  }
  const double local_col_BFRT_density =
      (double)col_BFRT.count * inv_solver_num_row;
  ekk_instance_.updateOperationResultDensity(
      local_col_BFRT_density, ekk_instance_.info_.col_BFRT_density);
}

void HEkkDual::updateFtranDSE(HVector* DSE_Vector) {
  // Compute the vector required to update DSE weights - being FTRAN
  // applied to row_ep (FTRAN-DSE)
  //
  // When solving the unscaled LP with scaled NLA, have computed
  //
  // row_ep = R\bar{B}^{-T}C_B.e_p = cb.R\bar{B}^{-T}e_p
  //
  // to get row_ep in the unscaled space.
  //
  // To update DSE weights requires \bar{B}^{-1}\bar{B}^{-T}e_p, being
  //
  // (1/cp).\bar{B}^{-1}R^{-1}(cb.R\bar{B}^{-T}e_p)
  //
  // where cb.R\bar{B}^{-T}e_p is row_ep in the unscaled space.
  //
  // Operation R^{-1} is performed here. Operation (1/cp) is performed
  // when updating weights
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;
  analysis->simplexTimerStart(FtranDseClock);
  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordBefore(kSimplexNlaFtranDse, *DSE_Vector,
                                    ekk_instance_.info_.row_DSE_density);
  // Apply R{-1}
  simplex_nla->unapplyBasisMatrixRowScale(*DSE_Vector);

  // Perform FTRAN DSE
  simplex_nla->ftranInScaledSpace(*DSE_Vector,
                                  ekk_instance_.info_.row_DSE_density,
                                  analysis->pointer_serial_factor_clocks);
  if (analysis->analyse_simplex_summary_data)
    analysis->operationRecordAfter(kSimplexNlaFtranDse, *DSE_Vector);
  analysis->simplexTimerStop(FtranDseClock);
  const double local_row_DSE_density =
      (double)DSE_Vector->count * inv_solver_num_row;
  ekk_instance_.updateOperationResultDensity(
      local_row_DSE_density, ekk_instance_.info_.row_DSE_density);
}

void HEkkDual::updateVerify() {
  // Compare the pivot value computed row-wise and column-wise and
  // determine whether reinversion is advisable
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;

  // Use the two pivot values to identify numerical trouble
  if (ekk_instance_.reinvertOnNumericalTrouble(
          "HEkkDual::updateVerify", numericalTrouble, alpha_col, alpha_row,
          kNumericalTroubleTolerance)) {
    rebuild_reason = kRebuildReasonPossiblySingularBasis;
  }
}

void HEkkDual::updateDual() {
  // Update the dual values
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;

  // Update - dual (shift and back)
  if (theta_dual == 0) {
    // Little to do if theta_dual is zero
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase,
    //    "Before shift_cost");
    shiftCost(variable_in, -workDual[variable_in]);
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase,
    //    "After shift_cost");
  } else {
    // Update the whole vector of dual values
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase,
    //    "Before calling dualRow.updateDual");
    dualRow.updateDual(theta_dual);
    if (ekk_instance_.info_.simplex_strategy != kSimplexStrategyDualPlain &&
        slice_PRICE) {
      // Update the slice-by-slice copy of dual variables
      for (HighsInt i = 0; i < slice_num; i++)
        slice_dualRow[i].updateDual(theta_dual);
    }
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase,
    //    "After calling dualRow.updateDual");
  }
  // Identify the changes in the dual objective
  double dual_objective_value_change;
  const double variable_in_delta_dual = workDual[variable_in];
  const double variable_in_value = workValue[variable_in];
  const HighsInt variable_in_nonbasicFlag =
      ekk_instance_.basis_.nonbasicFlag_[variable_in];
  dual_objective_value_change =
      variable_in_nonbasicFlag * (-variable_in_value * variable_in_delta_dual);
  dual_objective_value_change *= ekk_instance_.cost_scale_;
  ekk_instance_.info_.updated_dual_objective_value +=
      dual_objective_value_change;
  // Surely variable_out_nonbasicFlag is always 0 since it's basic - so there's
  // no dual objective change
  const HighsInt variable_out_nonbasicFlag =
      ekk_instance_.basis_.nonbasicFlag_[variable_out];
  assert(variable_out_nonbasicFlag == 0);
  if (variable_out_nonbasicFlag) {
    const double variable_out_delta_dual = workDual[variable_out] - theta_dual;
    const double variable_out_value = workValue[variable_out];
    dual_objective_value_change =
        variable_out_nonbasicFlag *
        (-variable_out_value * variable_out_delta_dual);
    dual_objective_value_change *= ekk_instance_.cost_scale_;
    ekk_instance_.info_.updated_dual_objective_value +=
        dual_objective_value_change;
  }
  workDual[variable_in] = 0;
  workDual[variable_out] = -theta_dual;

  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase, "Before
  //  shift_back");
  shiftBack(variable_out);
  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase, "After
  //  shift_back");
}

void HEkkDual::updatePrimal(HVector* DSE_Vector) {
  // Update the primal values and any edge weights
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;
  std::vector<double>& edge_weight = ekk_instance_.dual_edge_weight_;
  if (edge_weight_mode == EdgeWeightMode::kDevex) {
    const double updated_edge_weight = edge_weight[row_out];
    edge_weight[row_out] = computed_edge_weight;
    new_devex_framework = newDevexFramework(updated_edge_weight);
  }
  // DSE_Vector is either col_DSE = B^{-1}B^{-T}e_p (if using dual
  // steepest edge weights) or row_ep = B^{-T}e_p.
  //
  // Update - primal and weight
  dualRHS.updatePrimal(&col_BFRT, 1);
  dualRHS.updateInfeasList(&col_BFRT);
  double x_out = baseValue[row_out];
  double l_out = baseLower[row_out];
  double u_out = baseUpper[row_out];
  theta_primal = (x_out - (delta_primal < 0 ? l_out : u_out)) / alpha_col;
  const bool ok_update_primal = dualRHS.updatePrimal(&col_aq, theta_primal);
  if (!ok_update_primal) {
    rebuild_reason = kRebuildReasonExcessivePrimalValue;
    return;
  }
  ekk_instance_.updateBadBasisChange(col_aq, theta_primal);
  if (edge_weight_mode == EdgeWeightMode::kSteepestEdge) {
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
                                                DSE_Vector->array.data());
    edge_weight[row_out] = new_pivotal_edge_weight;
  } else if (edge_weight_mode == EdgeWeightMode::kDevex) {
    // Pivotal row is for the current basis: weights are required for
    // the next basis so have to divide the current (exact) weight by
    // the pivotal value
    double new_pivotal_edge_weight =
        edge_weight[row_out] / (alpha_col * alpha_col);
    new_pivotal_edge_weight = max(1.0, new_pivotal_edge_weight);
    // nw_wt is max(use_edge_weight_[iRow], NewExactWeight*columnArray[iRow]^2);
    //
    // But NewExactWeight is new_pivotal_edge_weight = max(1.0,
    // edge_weight[row_out] / (alpha * alpha))
    //
    // so nw_wt = max(use_edge_weight_[iRow],
    // new_pivotal_edge_weight*columnArray[iRow]^2);
    //
    // Update rest of weights
    ekk_instance_.updateDualDevexWeights(&col_aq, new_pivotal_edge_weight);
    edge_weight[row_out] = new_pivotal_edge_weight;
    num_devex_iterations++;
  }
  dualRHS.updateInfeasList(&col_aq);

  // Whether or not dual steepest edge weights are being used, have to
  // add in DSE_Vector->synthetic_tick_ since this contains the
  // contribution from forming row_ep = B^{-T}e_p.
  ekk_instance_.total_synthetic_tick_ += col_aq.synthetic_tick;
  ekk_instance_.total_synthetic_tick_ += DSE_Vector->synthetic_tick;
}

// Record the shift in the cost of a particular column
void HEkkDual::shiftCost(const HighsInt iCol, const double amount) {
  HighsSimplexInfo& info = ekk_instance_.info_;
  info.costs_shifted = true;
  assert(info.workShift_[iCol] == 0);
  if (!amount) return;
  double use_amount = amount;
  info.workShift_[iCol] = use_amount;
  // Analysis
  const double shift = fabs(use_amount);
  analysis->net_num_single_cost_shift++;
  analysis->num_single_cost_shift++;
  analysis->sum_single_cost_shift += shift;
  analysis->max_single_cost_shift = max(shift, analysis->max_single_cost_shift);
}

// Undo the shift in the cost of a particular column
void HEkkDual::shiftBack(const HighsInt iCol) {
  HighsSimplexInfo& info = ekk_instance_.info_;
  if (!info.workShift_[iCol]) return;
  info.workDual_[iCol] -= info.workShift_[iCol];
  info.workShift_[iCol] = 0;
  // Analysis
  analysis->net_num_single_cost_shift--;
}

void HEkkDual::updatePivots() {
  // UPDATE
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;
  // Transform the vectors used in updateFactor if the simplex NLA involves
  // scaling
  ekk_instance_.transformForUpdate(&col_aq, &row_ep, variable_in, &row_out);
  //
  // Update the sets of indices of basic and nonbasic variables
  //
  // debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase, "Before
  // update_pivots");
  ekk_instance_.updatePivots(variable_in, row_out, move_out);
  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solve_phase, "After
  //  update_pivots");
  //
  ekk_instance_.iteration_count_++;
  //
  // Update the invertible representation of the basis matrix
  ekk_instance_.updateFactor(&col_aq, &row_ep, &row_out, &rebuild_reason);
  //
  // Update the row-wise representation of the nonbasic columns
  ekk_instance_.updateMatrix(variable_in, variable_out);
  //
  // Delete Freelist entry for variable_in
  dualRow.deleteFreelist(variable_in);
  //
  // Update the primal value for the row where the basis change has
  // occurred, and set the corresponding primal infeasibility value in
  // dualRHS.work_infeasibility
  dualRHS.updatePivots(
      row_out, ekk_instance_.info_.workValue_[variable_in] + theta_primal);

  /*
  // Determine whether to reinvert based on the synthetic clock
  bool reinvert_syntheticClock = total_synthetic_tick >= build_synthetic_tick;
  const bool performed_min_updates =
      ekk_instance_.info_.update_count >=
      kSyntheticTickReinversionMinUpdateCount;
  if (reinvert_syntheticClock && performed_min_updates)
    rebuild_reason = kRebuildReasonSyntheticClockSaysInvert;
  */
}

void HEkkDual::initialiseDevexFramework() {
  HighsSimplexInfo& info = ekk_instance_.info_;
  // Initialise the Devex framework: reference set is all basic
  // variables
  analysis->simplexTimerStart(DevexIzClock);
  // Resize in case this is the first call
  ekk_instance_.info_.devex_index_.resize(solver_num_tot);
  const vector<int8_t>& nonbasicFlag = ekk_instance_.basis_.nonbasicFlag_;
  // Initialise the devex framework. The devex reference set is
  // initialise to be the current set of basic variables - and never
  // changes until a new framework is set up. In a simplex iteration,
  // to compute the exact Devex weight for the pivotal row requires
  // summing the squares of the its entries over the indices in the
  // reference set. This is achieved by summing over all indices, but
  // multiplying the entry by the value in devex_index before
  // equating. Thus devex_index contains 1 for indices in the
  // reference set, and 0 otherwise. This is achieved by setting the
  // values of devex_index to be 1-nonbasicFlag^2, ASSUMING
  // |nonbasicFlag|=1 iff the corresponding variable is nonbasic
  for (HighsInt vr_n = 0; vr_n < solver_num_tot; vr_n++)
    info.devex_index_[vr_n] = 1 - nonbasicFlag[vr_n] * nonbasicFlag[vr_n];
  // Set all initial weights to 1, zero the count of iterations with
  // this Devex framework, increment the number of Devex frameworks
  // and indicate that there's no need for a new Devex framework
  ekk_instance_.dual_edge_weight_.assign(solver_num_row, 1.0);
  num_devex_iterations = 0;
  new_devex_framework = false;
  minor_new_devex_framework = false;
  analysis->simplexTimerStop(DevexIzClock);
}

void HEkkDual::interpretDualEdgeWeightStrategy(
    const HighsInt dual_edge_weight_strategy) {
  if (dual_edge_weight_strategy == kSimplexEdgeWeightStrategyChoose) {
    edge_weight_mode = EdgeWeightMode::kSteepestEdge;
    allow_dual_steepest_edge_to_devex_switch = true;
  } else if (dual_edge_weight_strategy == kSimplexEdgeWeightStrategyDantzig) {
    edge_weight_mode = EdgeWeightMode::kDantzig;
  } else if (dual_edge_weight_strategy == kSimplexEdgeWeightStrategyDevex) {
    edge_weight_mode = EdgeWeightMode::kDevex;
  } else if (dual_edge_weight_strategy ==
             kSimplexEdgeWeightStrategySteepestEdge) {
    edge_weight_mode = EdgeWeightMode::kSteepestEdge;
    allow_dual_steepest_edge_to_devex_switch = false;
  } else {
    assert(1 == 0);
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "HEkkDual::interpretDualEdgeWeightStrategy: "
                "unrecognised dual_edge_weight_strategy = %" HIGHSINT_FORMAT
                " - using "
                "dual steepest edge with possible switch to Devex\n",
                dual_edge_weight_strategy);
    edge_weight_mode = EdgeWeightMode::kSteepestEdge;
    allow_dual_steepest_edge_to_devex_switch = true;
  }
}

void HEkkDual::possiblyUseLiDualSteepestEdge() {
  // Decide whether to use LiDSE by not storing squared primal infeasibilities
  HighsOptions& options = *ekk_instance_.options_;
  HighsSimplexInfo& info = ekk_instance_.info_;
  info.store_squared_primal_infeasibility = true;
  if (options.less_infeasible_DSE_check) {
    if (isLessInfeasibleDSECandidate(options.log_options, ekk_instance_.lp_)) {
      // LP is a candidate for LiDSE
      if (options.less_infeasible_DSE_choose_row)
        // Use LiDSE
        info.store_squared_primal_infeasibility = false;
    }
  }
}

void HEkkDual::computeDualInfeasibilitiesWithFixedVariableFlips() {
  // Computes num/max/sum of dual infeasibilities, ignoring fixed
  // variables whose infeasibilities can be corrected by flipping at
  // the fixed value, so that decisions on the dual simplex phase can
  // be taken. It is driven by the use of nonbasicMove to identify
  // dual infeasibilities, using the bounds only to identify free
  // variables. Fixed variables are assumed to have nonbasicMove=0 so
  // that no dual infeasibility is counted for them. Indeed, when
  // called from cleanup() at the end of dual phase 1, nonbasicMove
  // relates to the phase 1 bounds, but workLower and workUpper will
  // have been set to phase 2 values! Note that there can be no free
  // variables in dual phase 1.
  HighsInt num_dual_infeasibility = 0;
  double max_dual_infeasibility = 0;
  double sum_dual_infeasibility = 0;
  HighsLp& lp = ekk_instance_.lp_;
  SimplexBasis& basis = ekk_instance_.basis_;
  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsOptions* options = ekk_instance_.options_;
  const HighsInt num_tot = lp.num_col_ + lp.num_row_;
  // Possibly verify that nonbasicMove is correct for fixed variables
  // debugFixedNonbasicMove(ekk_instance_);
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (!basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double lower = info.workLower_[iVar];
    const double upper = info.workUpper_[iVar];
    const double dual = info.workDual_[iVar];
    double dual_infeasibility = 0;
    if (lower == -kHighsInf && upper == kHighsInf) {
      dual_infeasibility = fabs(dual);
    } else {
      dual_infeasibility = -basis.nonbasicMove_[iVar] * dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= options->dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
  info.num_dual_infeasibilities = num_dual_infeasibility;
  info.max_dual_infeasibility = max_dual_infeasibility;
  info.sum_dual_infeasibilities = sum_dual_infeasibility;
}

void HEkkDual::correctDualInfeasibilities(HighsInt& free_infeasibility_count) {
  // Removes dual infeasibilities for all but free variables. For
  // fixed variables, dual infeasibilities are removed by flipping at
  // the bound. Otherwise, dual infeasibilities are removed by
  // shifting costs.
  HighsLp& lp = ekk_instance_.lp_;
  SimplexBasis& basis = ekk_instance_.basis_;
  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsSimplexAnalysis& analysis = ekk_instance_.analysis_;
  HighsRandom& random = ekk_instance_.random_;
  HighsOptions* options = ekk_instance_.options_;

  free_infeasibility_count = 0;
  const double dual_feasibility_tolerance = options->dual_feasibility_tolerance;
  double flip_dual_objective_value_change = 0;
  double shift_dual_objective_value_change = 0;
  HighsInt num_flip = 0;
  HighsInt num_shift = 0;
  double sum_flip = 0;
  double sum_shift = 0;
  double max_flip = 0;
  double max_shift = 0;
  double min_dual_infeasibility_for_flip = kHighsInf;
  double max_dual_infeasibility_for_flip = 0;
  HighsInt num_dual_infeasibilities_for_flip = 0;
  double sum_dual_infeasibilities_for_flip = 0;
  HighsInt num_dual_infeasibilities_for_shift = 0;
  double max_dual_infeasibility_for_shift = 0;
  double sum_dual_infeasibilities_for_shift = 0;
  const HighsInt num_tot = lp.num_col_ + lp.num_row_;
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (!basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double lower = info.workLower_[iVar];
    const double upper = info.workUpper_[iVar];
    const double current_dual = info.workDual_[iVar];
    const HighsInt move = basis.nonbasicMove_[iVar];
    const bool fixed = lower == upper;
    const bool boxed = lower > -kHighsInf && upper < kHighsInf;
    const bool free = lower == -kHighsInf && upper == kHighsInf;
    double dual_infeasibility = 0;
    if (free) {
      dual_infeasibility = fabs(current_dual);
      if (dual_infeasibility >= dual_feasibility_tolerance)
        free_infeasibility_count++;
      continue;
    }
    dual_infeasibility = -move * current_dual;
    if (dual_infeasibility < dual_feasibility_tolerance) continue;
    // There is a dual infeasibility to remove
    //
    // force_phase2 is set true to prevent flipping of non-fixed
    // (boxed) variables when correcting infeasibilities in the first
    // set of duals computed after cost perturbation
    if (fixed || (boxed && !force_phase2)) {
      // Flip for fixed variables and boxed variables when not forcing phase 2
      ekk_instance_.flipBound(iVar);
      // Negative dual at lower bound (move=1): flip to upper
      // bound so objective contribution is change in value (flip)
      // times dual, being move*flip*dual
      //
      // Positive dual at upper bound (move=-1): flip to lower
      // bound so objective contribution is change in value
      // (-flip) times dual, being move*flip*dual
      const double flip = upper - lower;
      double local_dual_objective_change = move * flip * current_dual;
      local_dual_objective_change *= ekk_instance_.cost_scale_;
      flip_dual_objective_value_change += local_dual_objective_change;
      num_flip++;
      max_flip = max(fabs(flip), max_flip);
      sum_flip += fabs(flip);
      // Flipping fixed variables is trivial, so only track the
      // infeasibilities involved when flipping boxed variables
      if (!fixed) {
        min_dual_infeasibility_for_flip =
            std::min(dual_infeasibility, min_dual_infeasibility_for_flip);
        if (dual_infeasibility >= dual_feasibility_tolerance)
          num_dual_infeasibilities_for_flip++;
        sum_dual_infeasibilities_for_flip += dual_infeasibility;
        max_dual_infeasibility_for_flip =
            std::max(dual_infeasibility, max_dual_infeasibility_for_flip);
      }
      continue;
    }
    // Either boxed but not fixed, of one-sided, so shift
    //
    // Cost shifting must always be possible
    assert(info.allow_cost_shifting);
    // Other variable = shift
    if (dual_infeasibility >= dual_feasibility_tolerance)
      num_dual_infeasibilities_for_shift++;
    sum_dual_infeasibilities_for_shift += dual_infeasibility;
    max_dual_infeasibility_for_shift =
        std::max(dual_infeasibility, max_dual_infeasibility_for_shift);
    info.costs_shifted = true;
    double shift;
    if (move == kNonbasicMoveUp) {
      double new_dual = (1 + random.fraction()) * dual_feasibility_tolerance;
      shift = new_dual - current_dual;
      info.workDual_[iVar] = new_dual;
      info.workCost_[iVar] = info.workCost_[iVar] + shift;
    } else {
      double new_dual = -(1 + random.fraction()) * dual_feasibility_tolerance;
      shift = new_dual - current_dual;
      info.workDual_[iVar] = new_dual;
      info.workCost_[iVar] = info.workCost_[iVar] + shift;
    }
    double local_dual_objective_change = shift * info.workValue_[iVar];
    local_dual_objective_change *= ekk_instance_.cost_scale_;
    shift_dual_objective_value_change += local_dual_objective_change;
    num_shift++;
    max_shift = max(fabs(shift), max_shift);
    sum_shift += fabs(shift);
    const std::string direction = move == kNonbasicMoveUp ? "  up" : "down";
    highsLogDev(options->log_options, HighsLogType::kVerbose,
                "Move %s: cost shift = %g; objective change = %g\n",
                direction.c_str(), shift, local_dual_objective_change);
  }
  analysis.num_correct_dual_primal_flip += num_flip;
  analysis.max_correct_dual_primal_flip =
      max(max_flip, analysis.max_correct_dual_primal_flip);
  analysis.min_correct_dual_primal_flip_dual_infeasibility =
      std::min(min_dual_infeasibility_for_flip,
               analysis.min_correct_dual_primal_flip_dual_infeasibility);
  if (num_flip && force_phase2) {
    highsLogDev(
        options->log_options, HighsLogType::kDetailed,
        "Performed num / max / sum = %" HIGHSINT_FORMAT
        " / %g / %g flip(s) for num / min / max / sum dual infeasibility of "
        "%" HIGHSINT_FORMAT " / %g / %g / %g; objective change = %g\n",
        num_flip, max_flip, sum_flip, num_dual_infeasibilities_for_flip,
        min_dual_infeasibility_for_flip, max_dual_infeasibility_for_flip,
        sum_dual_infeasibilities_for_flip, flip_dual_objective_value_change);
  }
  analysis.num_correct_dual_cost_shift += num_shift;
  analysis.max_correct_dual_cost_shift =
      max(max_shift, analysis.max_correct_dual_cost_shift);
  analysis.max_correct_dual_cost_shift_dual_infeasibility =
      max(max_dual_infeasibility_for_shift,
          analysis.max_correct_dual_cost_shift_dual_infeasibility);
  if (num_shift) {
    highsLogDev(
        options->log_options, HighsLogType::kDetailed,
        "Performed num / max / sum = %" HIGHSINT_FORMAT
        " / %g / %g shift(s) for num / max / sum dual infeasibility of "
        "%" HIGHSINT_FORMAT " / %g / %g; objective change = %g\n",
        num_shift, max_shift, sum_shift, num_dual_infeasibilities_for_shift,
        max_dual_infeasibility_for_shift, sum_dual_infeasibilities_for_shift,
        shift_dual_objective_value_change);
  }
  force_phase2 = false;
}

bool HEkkDual::proofOfPrimalInfeasibility() {
  return ekk_instance_.proofOfPrimalInfeasibility(row_ep, move_out, row_out);
}

void HEkkDual::saveDualRay() {
  assert(row_out >= 0);
  assert(move_out != kNoRaySign);
  ekk_instance_.dual_ray_record_.clear();
  ekk_instance_.dual_ray_record_.index = row_out;
  ekk_instance_.dual_ray_record_.sign = move_out;
}

void HEkkDual::assessPhase1Optimality() {
  // Should only be called when optimal in phase 1 (row_out == kNoRowChosen)
  // with nonzero dual activity, and after a fresh rebuild - so
  // "final" decisions can be made.
  assert(solve_phase == kSolvePhase1);
  assert(row_out == kNoRowChosen);
  assert(ekk_instance_.info_.dual_objective_value);
  assert(ekk_instance_.status_.has_fresh_rebuild);

  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsModelStatus& model_status = ekk_instance_.model_status_;
  double& dual_objective_value = info.dual_objective_value;
  // There are (possibly insignificant) LP dual infeasibilities that
  // can't be removed by dual Phase 1, so clean up any perturbations
  // before concluding dual infeasibility
  //
  // Interesting for Devs to know if this method is called at all
  highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
              "Optimal in phase 1 but not jumping to phase 2 since "
              "dual objective is %10.4g: Costs perturbed = %" HIGHSINT_FORMAT
              "\n",
              dual_objective_value, info.costs_perturbed);
  if (info.costs_perturbed) {
    // Clean up perturbation
    cleanup();
    assessPhase1OptimalityUnperturbed();
  } else {
    assert(dualInfeasCount == 0);
    assert(dual_objective_value != 0);
    assessPhase1OptimalityUnperturbed();
  }
  if (dualInfeasCount > 0) {
    // Must still be solve_phase = kSolvePhase1 since dual
    // infeasibilities with respect to phase 1 bounds mean that primal
    // values must change, so primal feasibility is unknown
    assert(solve_phase == kSolvePhase1);
  } else {
    // Optimal in dual phase 1, so either dual feasible wrt Phase 2
    // bounds and going to phase 2, or identified dual infeasibility and exiting
    assert(solve_phase == kSolvePhase2 ||
           (solve_phase == kSolvePhaseExit &&
            model_status == HighsModelStatus::kUnboundedOrInfeasible));
    if (solve_phase == kSolvePhase2) {
      // Reset the duals, if necessary shifting costs of free variables
      // so that their duals are zero
      exitPhase1ResetDuals();
    }
  }
}

void HEkkDual::assessPhase1OptimalityUnperturbed() {
  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsModelStatus& model_status = ekk_instance_.model_status_;
  double& dual_objective_value = info.dual_objective_value;
  assert(!info.costs_perturbed);
  if (dualInfeasCount == 0) {
    // No dual infeasibilities with respect to phase 1 bounds.
    if (dual_objective_value == 0) {
      // No dual infeasibilities with respect to phase 2 bounds so
      // go to phase 2
      highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                  "LP is dual feasible wrt Phase 2 bounds after removing cost "
                  "perturbations so go to phase 2\n");
      solve_phase = kSolvePhase2;
    } else {
      // Nonzero dual objective value: could be insignificant dual
      // infeasibilities
      highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                  "LP is dual feasible wrt Phase 1 bounds after removing cost "
                  "perturbations: "
                  "dual objective is %10.4g\n",
                  dual_objective_value);
      ekk_instance_.computeSimplexLpDualInfeasible();
      const HighsInt num_lp_dual_infeasibilities =
          ekk_instance_.analysis_.num_dual_phase_1_lp_dual_infeasibility;
      if (num_lp_dual_infeasibilities == 0) {
        highsLogDev(
            ekk_instance_.options_->log_options, HighsLogType::kInfo,
            "LP is dual feasible wrt Phase 2 bounds after removing cost "
            "perturbations so go to phase 2\n");
        solve_phase = kSolvePhase2;
      } else {
        // LP is dual infeasible if the dual objective is sufficiently
        // negative, so no conclusions on the primal LP can be deduced
        // - could be primal unbounded or primal infeasible.
        //
        // Indicate the conclusion of dual infeasibility by setting
        // the scaled model status
        reportOnPossibleLpDualInfeasibility();
        model_status = HighsModelStatus::kUnboundedOrInfeasible;
        solve_phase = kSolvePhaseExit;
      }
    }
  } else {
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "LP has %d dual feasibilities wrt Phase 1 bounds after "
                "removing cost perturbations "
                "so return to phase 1\n",
                dualInfeasCount);
    assert(solve_phase == kSolvePhase1);
  }
}

void HEkkDual::exitPhase1ResetDuals() {
  const HighsLp& lp = ekk_instance_.lp_;
  const SimplexBasis& basis = ekk_instance_.basis_;
  HighsSimplexInfo& info = ekk_instance_.info_;
  // This use of costs_alt_perturbed is not executed by ctest
  //
  //  assert(99==2);
  if (info.costs_perturbed) {
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "Costs are already perturbed in exitPhase1ResetDuals\n");
  } else {
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kDetailed,
                "Re-perturbing costs when optimal in phase 1\n");
    ekk_instance_.initialiseCost(SimplexAlgorithm::kDual, kSolvePhase2, true);
    ekk_instance_.computeDual();
  }

  const HighsInt numTot = lp.num_col_ + lp.num_row_;
  HighsInt num_shift = 0;
  double sum_shift = 0;
  for (HighsInt iVar = 0; iVar < numTot; iVar++) {
    if (basis.nonbasicFlag_[iVar]) {
      double lp_lower;
      double lp_upper;
      if (iVar < lp.num_col_) {
        lp_lower = lp.col_lower_[iVar];
        lp_upper = lp.col_upper_[iVar];
      } else {
        HighsInt iRow = iVar - lp.num_col_;
        lp_lower = lp.row_lower_[iRow];
        lp_upper = lp.row_upper_[iRow];
      }
      if (lp_lower <= -kHighsInf && lp_upper >= kHighsInf) {
        const double shift = -info.workDual_[iVar];
        info.workDual_[iVar] = 0;
        info.workCost_[iVar] = info.workCost_[iVar] + shift;
        num_shift++;
        sum_shift += fabs(shift);
        highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kVerbose,
                    "Variable %" HIGHSINT_FORMAT
                    " is free: shift cost to zero dual of %g\n",
                    iVar, shift);
      }
    }
  }
  if (num_shift) {
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kDetailed,
                "Performed %" HIGHSINT_FORMAT
                " cost shift(s) for free variables to zero "
                "dual values: total = %g\n",
                num_shift, sum_shift);
    info.costs_shifted = true;
  }
}

void HEkkDual::reportOnPossibleLpDualInfeasibility() {
  HighsSimplexInfo& info = ekk_instance_.info_;
  HighsSimplexAnalysis& analysis = ekk_instance_.analysis_;
  assert(solve_phase == kSolvePhase1);
  assert(row_out == kNoRowChosen);
  //  assert(info.dual_objective_value < 0);
  assert(!info.costs_perturbed);
  std::string lp_dual_status;
  if (analysis.num_dual_phase_1_lp_dual_infeasibility) {
    lp_dual_status = "infeasible";
  } else {
    lp_dual_status = "feasible";
  }
  highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
              "LP is dual %s with dual phase 1 objective %10.4g and num / "
              "max / sum dual infeasibilities = %" HIGHSINT_FORMAT
              " / %9.4g / %9.4g\n",
              lp_dual_status.c_str(), info.dual_objective_value,
              analysis.num_dual_phase_1_lp_dual_infeasibility,
              analysis.max_dual_phase_1_lp_dual_infeasibility,
              analysis.sum_dual_phase_1_lp_dual_infeasibility);
}

bool HEkkDual::dualInfoOk(const HighsLp& lp) const {
  HighsInt lp_num_col = lp.num_col_;
  HighsInt lp_num_row = lp.num_row_;
  bool dimensions_ok;
  dimensions_ok = lp_num_col == solver_num_col && lp_num_row == solver_num_row;
  assert(dimensions_ok);
  if (!dimensions_ok) {
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kError,
                "LP-Solver dimension incompatibility (%" HIGHSINT_FORMAT
                ", %" HIGHSINT_FORMAT ") != (%" HIGHSINT_FORMAT
                ", %" HIGHSINT_FORMAT ")\n",
                lp_num_col, solver_num_col, lp_num_row, solver_num_row);
    return false;
  }
  dimensions_ok = lp_num_col == simplex_nla->lp_->num_col_ &&
                  lp_num_row == simplex_nla->lp_->num_row_;
  assert(dimensions_ok);
  if (!dimensions_ok) {
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kError,
                "LP-Factor dimension incompatibility (%" HIGHSINT_FORMAT
                ", %" HIGHSINT_FORMAT ") != (%" HIGHSINT_FORMAT
                ", %" HIGHSINT_FORMAT ")\n",
                lp_num_col, simplex_nla->lp_->num_col_, lp_num_row,
                simplex_nla->lp_->num_row_);
    return false;
  }
  return true;
}

bool HEkkDual::bailoutOnDualObjective() {
  if (ekk_instance_.solve_bailout_) {
    // Bailout has already been decided: check that it's for one of these
    // reasons
    assert(ekk_instance_.model_status_ == HighsModelStatus::kTimeLimit ||
           ekk_instance_.model_status_ == HighsModelStatus::kIterationLimit ||
           ekk_instance_.model_status_ == HighsModelStatus::kObjectiveBound);
  } else if (ekk_instance_.lp_.sense_ == ObjSense::kMinimize &&
             solve_phase == kSolvePhase2) {
    if (ekk_instance_.info_.updated_dual_objective_value >
        ekk_instance_.options_->objective_bound)
      ekk_instance_.solve_bailout_ = reachedExactObjectiveBound();
  }
  return ekk_instance_.solve_bailout_;
}

bool HEkkDual::reachedExactObjectiveBound() {
  // Solving a minimization in dual simplex phase 2, and dual
  // objective exceeds the prescribed upper bound. However, costs
  // will be perturbed, so need to check whether exact dual
  // objective value exceeds the prescribed upper bound. This can be
  // a relatively expensive calculation, so determine whether to do
  // it according to the sparsity of the pivotal row
  bool reached_exact_objective_bound = false;
  double use_row_ap_density =
      std::min(std::max(ekk_instance_.info_.row_ap_density, 0.01), 1.0);
  HighsInt check_frequency = 1.0 / use_row_ap_density;
  assert(check_frequency > 0);

  bool check_exact_dual_objective_value =
      ekk_instance_.info_.update_count % check_frequency == 0;

  if (check_exact_dual_objective_value) {
    const double objective_bound = ekk_instance_.options_->objective_bound;
    const double perturbed_dual_objective_value =
        ekk_instance_.info_.updated_dual_objective_value;
    const double perturbed_value_residual =
        perturbed_dual_objective_value - objective_bound;
    HVector dual_col;
    HVector dual_row;
    const double exact_dual_objective_value =
        computeExactDualObjectiveValue(dual_col, dual_row);
    const double exact_value_residual =
        exact_dual_objective_value - objective_bound;
    std::string action;
    if (exact_dual_objective_value > objective_bound) {
      highsLogDev(
          ekk_instance_.options_->log_options, HighsLogType::kDetailed,
          "HEkkDual::solvePhase2: %12g = Objective > ObjectiveUB = %12g\n",
          ekk_instance_.info_.updated_dual_objective_value, objective_bound);
      action = "Have DualUB bailout";
      if (ekk_instance_.info_.costs_perturbed ||
          ekk_instance_.info_.costs_shifted) {
        // Remove cost perturbation/shifting
        ekk_instance_.initialiseCost(SimplexAlgorithm::kDual, kSolvePhase2);
      }

      // Set the duals as computed in the computeExactDualObjective call
      for (HighsInt i = 0; i < solver_num_col; i++)
        ekk_instance_.info_.workDual_[i] =
            ekk_instance_.info_.workCost_[i] - dual_row.array[i];
      for (HighsInt i = solver_num_col; i < solver_num_tot; i++)
        ekk_instance_.info_.workDual_[i] = -dual_col.array[i - solver_num_col];

      // Since the computeExactDualObjectiveValue() call succeeded, if there are
      // any dual infeasibilities they can be removed by a bound flip
      force_phase2 = false;
      correctDualInfeasibilities(dualInfeasCount);

      // no shifts should have occurred
      assert(!ekk_instance_.info_.costs_shifted);
      reached_exact_objective_bound = true;
      ekk_instance_.model_status_ = HighsModelStatus::kObjectiveBound;
    } else {
      action = "No   DualUB bailout";
    }
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "%s on iteration %" HIGHSINT_FORMAT
                ": Density %11.4g; Frequency %" HIGHSINT_FORMAT
                ": "
                "Residual(Perturbed = %g; Exact = %g)\n",
                action.c_str(), ekk_instance_.iteration_count_,
                use_row_ap_density, check_frequency, perturbed_value_residual,
                exact_value_residual);
  }
  return reached_exact_objective_bound;
}

double HEkkDual::computeExactDualObjectiveValue(HVector& dual_col,
                                                HVector& dual_row) {
  const HighsLp& lp = ekk_instance_.lp_;
  const SimplexBasis& basis = ekk_instance_.basis_;
  const HighsSimplexInfo& info = ekk_instance_.info_;
  // Create a local buffer for the pi vector
  dual_col.setup(lp.num_row_);
  dual_col.clear();
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    HighsInt iVar = basis.basicIndex_[iRow];
    if (iVar < lp.num_col_) {
      const double value = lp.col_cost_[iVar];
      if (value) {
        dual_col.array[iRow] = value;
        dual_col.index[dual_col.count++] = iRow;
      }
    }
  }
  // Create a local buffer for the dual vector
  const HighsInt numTot = lp.num_col_ + lp.num_row_;
  dual_row.setup(lp.num_col_);
  dual_row.clear();
  if (dual_col.count) {
    const bool quad_precision = false;
    const double expected_density = 1;
    simplex_nla->btran(dual_col, expected_density);
    lp.a_matrix_.priceByColumn(quad_precision, dual_row, dual_col);
  }
  // Compute dual infeasibilities
  ekk_instance_.computeSimplexDualInfeasible();
  if (info.num_dual_infeasibilities > 0)
    highsLogDev(ekk_instance_.options_->log_options, HighsLogType::kInfo,
                "When computing exact dual objective, the unperturbed costs "
                "yield num / max / sum dual "
                "infeasibilities = %d / %g / %g\n",
                (int)info.num_dual_infeasibilities, info.max_dual_infeasibility,
                info.sum_dual_infeasibilities);
  HighsCDouble dual_objective = lp.offset_;
  double norm_dual = 0;
  double norm_delta_dual = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (!basis.nonbasicFlag_[iCol]) continue;
    double exact_dual = lp.col_cost_[iCol] - dual_row.array[iCol];
    // The active value must be decided based on the exact dual. For a nonbasic
    // column the bound that must be used may flip due to cost perturbation
    // flipping the sign of its dual and for a basic variable we may need to add
    // to the dual objective using one of the bounds when its dual is not zero.
    double active_value;
    if (exact_dual > ekk_instance_.options_->small_matrix_value)
      active_value = lp.col_lower_[iCol];
    else if (exact_dual < -ekk_instance_.options_->small_matrix_value)
      active_value = lp.col_upper_[iCol];
    else
      active_value = info.workValue_[iCol];

    // when the active value is infinite the dual objective lower bound is
    // -infinity
    if (highs_isInfinity(fabs(active_value))) return -kHighsInf;

    double residual = fabs(exact_dual - info.workDual_[iCol]);
    norm_dual += fabs(exact_dual);
    norm_delta_dual += residual;
    if (residual > 1e10)
      highsLogDev(
          ekk_instance_.options_->log_options, HighsLogType::kWarning,
          "Col %4" HIGHSINT_FORMAT
          ": ExactDual = %11.4g; WorkDual = %11.4g; Residual = %11.4g\n",
          iCol, exact_dual, info.workDual_[iCol], residual);
    dual_objective += active_value * exact_dual;
  }

  for (HighsInt iVar = lp.num_col_; iVar < numTot; iVar++) {
    if (!basis.nonbasicFlag_[iVar]) continue;
    HighsInt iRow = iVar - lp.num_col_;
    double exact_dual = dual_col.array[iRow];
    // Similarly to the column case above the active value must be decided based
    // on the exact dual. For a nonbasic row the bound that must be used
    // may flip due to cost perturbation flipping the sign of its dual and for a
    // basic variable we may need to add to the dual objective using one of the
    // bounds when its dual is not zero.
    double active_value;
    if (exact_dual > ekk_instance_.options_->small_matrix_value)
      active_value = lp.row_lower_[iRow];
    else if (exact_dual < -ekk_instance_.options_->small_matrix_value)
      active_value = lp.row_upper_[iRow];
    else
      active_value = -info.workValue_[iVar];

    // when the active value is infinite the dual objective lower bound is
    // -infinity
    if (highs_isInfinity(fabs(active_value))) return -kHighsInf;

    double residual = fabs(exact_dual + info.workDual_[iVar]);
    norm_dual += fabs(exact_dual);
    norm_delta_dual += residual;
    if (residual > 1e10)
      highsLogDev(
          ekk_instance_.options_->log_options, HighsLogType::kWarning,
          "Row %4" HIGHSINT_FORMAT
          ": ExactDual = %11.4g; WorkDual = %11.4g; Residual = %11.4g\n",
          iRow, exact_dual, info.workDual_[iVar], residual);
    dual_objective += active_value * exact_dual;
  }
  double relative_delta = norm_delta_dual / std::max(norm_dual, 1.0);
  if (relative_delta > 1e-3)
    highsLogDev(
        ekk_instance_.options_->log_options, HighsLogType::kWarning,
        "||exact dual vector|| = %g; ||delta dual vector|| = %g: ratio = %g\n",
        norm_dual, norm_delta_dual, relative_delta);
  return double(dual_objective);
}

HighsDebugStatus HEkkDual::debugDualSimplex(const std::string& message,
                                            const bool initialise) {
  HighsDebugStatus return_status =
      ekk_instance_.debugSimplex(message, algorithm, solve_phase, initialise);
  if (return_status == HighsDebugStatus::kLogicalError) return return_status;
  if (initialise) return return_status;
  return HighsDebugStatus::kOk;
}

bool HEkkDual::isBadBasisChange() {
  return ekk_instance_.isBadBasisChange(SimplexAlgorithm::kDual, variable_in,
                                        row_out, rebuild_reason);
}

void HEkkDual::assessPossiblyDualUnbounded() {
  assert(rebuild_reason == kRebuildReasonPossiblyDualUnbounded);
  if (solve_phase != kSolvePhase2) return;
  if (!ekk_instance_.status_.has_fresh_rebuild) return;
  // Appears to be dual unbounded in phase 2 after fresh
  // rebuild. Normally this implies primal infeasibility, but only
  // allow this to be claimed if the proof of primal infeasibility
  // is true.
  //
  const bool proof_of_infeasibility = proofOfPrimalInfeasibility();
  if (proof_of_infeasibility) {
    // There is a proof of primal infeasibility
    solve_phase = kSolvePhaseExit;
    // Save dual ray information
    saveDualRay();
    // Model status should be unset?
    assert(ekk_instance_.model_status_ == HighsModelStatus::kNotset);
    ekk_instance_.model_status_ = HighsModelStatus::kInfeasible;
  } else {
    // No proof of primal infeasibility, so assume dual unbounded
    // claim is spurious. Make row_out taboo, and prevent rebuild
    ekk_instance_.addBadBasisChange(
        row_out, variable_out, variable_in,
        BadBasisChangeReason::kFailedInfeasibilityProof, true);
    rebuild_reason = kRebuildReasonNo;
  }
}
