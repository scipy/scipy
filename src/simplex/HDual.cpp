/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HDual.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HDual.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <set>
#include <stdexcept>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HCrash.h"
#include "simplex/HPrimal.h"
#include "simplex/HQPrimal.h"
#include "simplex/HSimplex.h"
#include "simplex/SimplexTimer.h"
#include "util/HighsTimer.h"

#ifdef OPENMP
#include "omp.h"
#endif

using std::cout;
using std::endl;
using std::fabs;
using std::flush;
using std::runtime_error;

HighsStatus HDual::solve() {
  HighsOptions& options = workHMO.options_;
  HighsSolutionParams& scaled_solution_params = workHMO.scaled_solution_params_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  workHMO.scaled_model_status_ = HighsModelStatus::NOTSET;
  bool simplex_info_ok =
      simplexInfoOk(workHMO.lp_, workHMO.simplex_lp_, simplex_info);
  if (!simplex_info_ok) {
    HighsLogMessage(workHMO.options_.logfile, HighsMessageType::ERROR,
                    "HPrimalDual::solve has error in simplex information");
    return HighsStatus::Error;
  }
  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLpSimplex
  bool positive_num_row = workHMO.simplex_lp_.numRow_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    HighsLogMessage(workHMO.options_.logfile, HighsMessageType::ERROR,
                    "HPrimal::solve called for LP with non-positive (%d) "
                    "number of constraints",
                    workHMO.simplex_lp_.numRow_);
    return HighsStatus::Error;
  }

  invertHint = INVERT_HINT_NO;

  // Set solve_bailout to be true if control is to be returned immediately to
  // calling function
  solve_bailout = false;

  // Initialise working environment. Does LOTS, including
  // initialisation of edge weights to 1s. Should only be called if
  // model dimension changes
  init();
  initParallel();

  bool dual_info_ok = dualInfoOk(workHMO.lp_);
  if (!dual_info_ok) {
    HighsLogMessage(workHMO.options_.logfile, HighsMessageType::ERROR,
                    "HPrimalDual::solve has error in dual information");
    return HighsStatus::Error;
  }

  // Decide whether to use LiDSE by not storing squared primal infeasibilities
  simplex_info.store_squared_primal_infeasibility = true;
  if (options.less_infeasible_DSE_check) {
    if (isLessInfeasibleDSECandidate(options, workHMO.simplex_lp_)) {
      // LP is a candidate for LiDSE
      if (options.less_infeasible_DSE_choose_row)
        // Use LiDSE
        simplex_info.store_squared_primal_infeasibility = false;
    }
  }

  initialise_cost(workHMO, 1);
  assert(simplex_lp_status.has_fresh_invert);
  if (!simplex_lp_status.has_fresh_invert) {
    printf(
        "ERROR: Should enter with fresh INVERT - unless no_invert_on_optimal "
        "is set\n");
  }
  // Consider initialising edge weights
  //
  // NB workEdWt is assigned and initialised to 1s in
  // dualRHS.setup(workHMO) so that CHUZR is well defined, even for
  // Dantzig pricing
  //
  if (!simplex_lp_status.has_dual_steepest_edge_weights) {
    // Edge weights are not known
    // Set up edge weights according to dual_edge_weight_mode and
    // initialise_dual_steepest_edge_weights
    if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) {
      // Using dual Devex edge weights, so set up the first framework
      simplex_info.devex_index_.assign(solver_num_tot, 0);
      initialiseDevexFramework();
    } else if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
      // Using dual steepest edge (DSE) weights
      int num_basic_structurals =
          solver_num_row - simplex_info.num_basic_logicals;
      bool computeExactDseWeights =
          num_basic_structurals > 0 && initialise_dual_steepest_edge_weights;
#ifdef HiGHSDEV
      if (computeExactDseWeights) {
        printf(
            "If (0<num_basic_structurals = %d) && %d = "
            "initialise_dual_steepest_edge_weights: Compute exact "
            "DSE weights\n",
            num_basic_structurals, initialise_dual_steepest_edge_weights);
      }
#endif
      if (computeExactDseWeights) {
        // Basis is not logical and DSE weights are to be initialised
#ifdef HiGHSDEV
        printf("Compute exact DSE weights\n");
        analysis->simplexTimerStart(SimplexIzDseWtClock);
        analysis->simplexTimerStart(DseIzClock);
#endif
        for (int i = 0; i < solver_num_row; i++) {
          row_ep.clear();
          row_ep.count = 1;
          row_ep.index[0] = i;
          row_ep.array[i] = 1;
          row_ep.packFlag = false;
          factor->btran(row_ep, analysis->row_ep_density,
                        analysis->pointer_serial_factor_clocks);
          dualRHS.workEdWt[i] = row_ep.norm2();
          const double local_row_ep_density =
              (double)row_ep.count / solver_num_row;
          analysis->updateOperationResultDensity(local_row_ep_density,
                                                 analysis->row_ep_density);
        }
#ifdef HiGHSDEV
        analysis->simplexTimerStop(SimplexIzDseWtClock);
        analysis->simplexTimerStop(DseIzClock);
        double IzDseWtTT = analysis->simplexTimerRead(SimplexIzDseWtClock);
        HighsPrintMessage(options.output, options.message_level, ML_DETAILED,
                          "Computed %d initial DSE weights in %gs\n",
                          solver_num_row, IzDseWtTT);
#endif
      }
#ifdef HiGHSDEV
      else {
        HighsPrintMessage(options.output, options.message_level, ML_DETAILED,
                          "solve:: %d basic structurals: starting from B=I so "
                          "unit initial DSE weights\n",
                          num_basic_structurals);
      }
#endif
    }
    // Indicate that edge weights are known
    simplex_lp_status.has_dual_steepest_edge_weights = true;
  }
  // Compute the dual values
  computeDual(workHMO);
  // Determine the number of dual infeasibilities, and hence the solve phase
  const bool analyse_dual_infeasibilities_with_and_without_flips = false;
  if (analyse_dual_infeasibilities_with_and_without_flips) {
    computeDualInfeasible(workHMO);
    int num_dual_infeasibilities_without_flips =
        scaled_solution_params.num_dual_infeasibilities;
    computeDualInfeasibleWithFlips(workHMO);
    int num_dual_infeasibilities_with_flips =
        scaled_solution_params.num_dual_infeasibilities;
    printf(
        "Dual infeasibilities with / without flips is %d / %d: Difference = "
        "%d\n",
        scaled_solution_params.num_dual_infeasibilities,
        num_dual_infeasibilities_without_flips,
        num_dual_infeasibilities_without_flips -
            num_dual_infeasibilities_with_flips);
  }
  if (simplex_info.allow_primal_flips_for_dual_feasibility) {
    computeDualInfeasibleWithFlips(workHMO);
  } else {
    computeDualInfeasible(workHMO);
  }
  dualInfeasCount = scaled_solution_params.num_dual_infeasibilities;
  solvePhase = dualInfeasCount > 0 ? 1 : 2;
  //
  // Check that the model is OK to solve:
  //
  // Level 0 just checks the flags
  //
  // Level 1 also checks that the basis is OK and that the necessary
  // data in work* is populated.
  //
  bool ok = ok_to_solve(workHMO, 1, solvePhase);
  if (!ok) {
    printf("NOT OK TO SOLVE???\n");
    cout << flush;
  }
  assert(ok);
#ifdef HiGHSDEV
  // reportSimplexLpStatus(simplex_lp_status, "Before HDual major solving
  // loop");
#endif
  //
  // The major solving loop
  //
  while (solvePhase) {
    int it0 = scaled_solution_params.simplex_iteration_count;
    // When starting a new phase the (updated) dual objective function
    // value isn't known. Indicate this so that when the value
    // computed from scratch in build() isn't checked against the the
    // updated value
    simplex_lp_status.has_dual_objective_value = false;
    switch (solvePhase) {
      case 1:
        analysis->simplexTimerStart(SimplexDualPhase1Clock);
        solvePhase1();
        analysis->simplexTimerStop(SimplexDualPhase1Clock);
        simplex_info.dual_phase1_iteration_count +=
            (scaled_solution_params.simplex_iteration_count - it0);
        break;
      case 2:
        analysis->simplexTimerStart(SimplexDualPhase2Clock);
        solvePhase2();
        analysis->simplexTimerStop(SimplexDualPhase2Clock);
        simplex_info.dual_phase2_iteration_count +=
            (scaled_solution_params.simplex_iteration_count - it0);
        break;
      case 4:
        break;
      default:
        solvePhase = 0;
        break;
    }
    // Possibly bail out
    if (bailout()) return HighsStatus::Warning;
    // Jump for primal
    if (solvePhase == 4) break;
  }
  // If bailing out, should have returned already
  assert(!solve_bailout);

  if (solvePhase != 0 && solvePhase != 4) {
    printf("What ho! HDual::solve() returning with solvePhase = %d\n",
           solvePhase);
  }
  if (solvePhase == 4) {
    computePrimalObjectiveValue(workHMO);
#ifdef HiGHSDEV
    vector<double> primal_value_before_cleanup;
    getPrimalValue(workHMO, primal_value_before_cleanup);
    //    printf("\nAnalyse primal objective evaluation before cleanup\n");
    //    analysePrimalObjectiveValue(workHMO);
    const double objective_before = simplex_info.primal_objective_value;
#endif
    if (options.dual_simplex_cleanup_strategy ==
        DUAL_SIMPLEX_CLEANUP_STRATEGY_NONE) {
      // No clean up. Dual simplex was optimal with perturbed costs,
      // so say that the scaled LP has been solved
      // optimally. Optimality (unlikely) for the unscaled LP will
      // still be assessed honestly, so leave it to the user to
      // deceide whether the solution can be accepted.
      workHMO.scaled_model_status_ = HighsModelStatus::OPTIMAL;
    } else {
      // Use primal to clean up
#ifdef HiGHSDEV
      initialiseValueDistribution("Cleanup primal step summary", "", 1e-16,
                                  1e16, 10.0,
                                  analysis->cleanup_primal_step_distribution);
      initialiseValueDistribution("Cleanup dual step summary", "", 1e-16, 1e16,
                                  10.0,
                                  analysis->cleanup_dual_step_distribution);
#endif
      int it0 = scaled_solution_params.simplex_iteration_count;
      const bool full_logging = false;  // true;//
      if (full_logging)
        analysis->messaging(options.logfile, options.output, ML_ALWAYS);
      analysis->simplexTimerStart(SimplexPrimalPhase2Clock);
      if (options.dual_simplex_cleanup_strategy ==
          DUAL_SIMPLEX_CLEANUP_STRATEGY_HPRIMAL) {
        // Cleanup with original primal phase 2 code
        HPrimal hPrimal(workHMO);
        hPrimal.solvePhase2();
      } else {
        // Cleanup with phase 2 for new primal code
        HQPrimal hPrimal(workHMO);
        hPrimal.solvePhase2();
      }
      analysis->simplexTimerStop(SimplexPrimalPhase2Clock);
#ifdef HiGHSDEV
      vector<double> primal_value_after_cleanup;
      getPrimalValue(workHMO, primal_value_after_cleanup);
      for (int var = 0; var < solver_num_tot; var++) {
        const double primal_change = fabs(primal_value_after_cleanup[var] -
                                          primal_value_before_cleanup[var]);
        updateValueDistribution(primal_change,
                                analysis->cleanup_primal_change_distribution);
      }
      //      printf("\nAnalyse primal objective evaluation after cleanup\n");
      //      analysePrimalObjectiveValue(workHMO);
      const double objective_after = simplex_info.primal_objective_value;
      const double abs_objective_change =
          fabs(objective_before - objective_after);
      const double rel_objective_change =
          abs_objective_change / max(1.0, fabs(objective_after));
      printf(
          "\nDuring cleanup, (abs: rel) primal objective changes is (%10.4g: "
          "%10.4g) \nfrom %20.10g\nto   %20.10g\n",
          abs_objective_change, rel_objective_change, objective_before,
          objective_after);
#endif
      simplex_info.primal_phase2_iteration_count +=
          (scaled_solution_params.simplex_iteration_count - it0);
    }
  }
  ok = ok_to_solve(workHMO, 1, solvePhase);
#ifdef HiGHSDEV
  if (!ok) printf("NOT OK After Solve???\n");
#endif
  assert(ok);
  computePrimalObjectiveValue(workHMO);
  return HighsStatus::OK;
}

void HDual::options() {
  // Set solver options from simplex options

  const HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  const HighsSolutionParams& scaled_solution_params =
      workHMO.scaled_solution_params_;

  interpretDualEdgeWeightStrategy(simplex_info.dual_edge_weight_strategy);
  //  interpretPriceStrategy(simplex_info.price_strategy);

  // Copy values of simplex solver options to dual simplex options
  primal_feasibility_tolerance =
      scaled_solution_params.primal_feasibility_tolerance;
  dual_feasibility_tolerance =
      scaled_solution_params.dual_feasibility_tolerance;
  dual_objective_value_upper_bound =
      workHMO.options_.dual_objective_value_upper_bound;
  //  perturb_costs = simplex_info.perturb_costs;
  //  iterationLimit = simplex_info.iterationLimit;

  // Set values of internal options
}

void HDual::init() {
  // Copy size, matrix and factor

  solver_num_col = workHMO.simplex_lp_.numCol_;
  solver_num_row = workHMO.simplex_lp_.numRow_;
  solver_num_tot = solver_num_col + solver_num_row;

  matrix = &workHMO.matrix_;
  factor = &workHMO.factor_;
  analysis = &workHMO.simplex_analysis_;

  // Copy pointers
  jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  workDual = &workHMO.simplex_info_.workDual_[0];
  workValue = &workHMO.simplex_info_.workValue_[0];
  workRange = &workHMO.simplex_info_.workRange_[0];
  baseLower = &workHMO.simplex_info_.baseLower_[0];
  baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  baseValue = &workHMO.simplex_info_.baseValue_[0];

  // Copy tolerances
  Tp = primal_feasibility_tolerance;
  Td = dual_feasibility_tolerance;

  // Setup local vectors
  col_DSE.setup(solver_num_row);
  col_BFRT.setup(solver_num_row);
  col_aq.setup(solver_num_row);
  row_ep.setup(solver_num_row);
  row_ap.setup(solver_num_col);
  // Setup other buffers
  dualRow.setup();
  dualRHS.setup();
}

void HDual::initParallel() {
  // Identify the (current) number of HiGHS tasks to be used
  const int num_threads = workHMO.simplex_info_.num_threads;

  // Initialize for tasks
  if (workHMO.simplex_info_.simplex_strategy == SIMPLEX_STRATEGY_DUAL_TASKS) {
    const int pass_num_slice = num_threads - 2;
    assert(pass_num_slice > 0);
    if (pass_num_slice <= 0) {
      HighsLogMessage(workHMO.options_.logfile, HighsMessageType::WARNING,
                      "SIP trying to use using %d slices due to number of "
                      "threads (%d) being too small: results unpredictable",
                      pass_num_slice, num_threads);
    }
    initSlice(pass_num_slice);
  }

  // Initialize for multi
  if (workHMO.simplex_info_.simplex_strategy == SIMPLEX_STRATEGY_DUAL_MULTI) {
    multi_num = num_threads;
    if (multi_num < 1) multi_num = 1;
    if (multi_num > HIGHS_THREAD_LIMIT) multi_num = HIGHS_THREAD_LIMIT;
    for (int i = 0; i < multi_num; i++) {
      multi_choice[i].row_ep.setup(solver_num_row);
      multi_choice[i].col_aq.setup(solver_num_row);
      multi_choice[i].col_BFRT.setup(solver_num_row);
    }
    const int pass_num_slice = max(multi_num - 1, 1);
    assert(pass_num_slice > 0);
    if (pass_num_slice <= 0) {
      HighsLogMessage(workHMO.options_.logfile, HighsMessageType::WARNING,
                      "PAMI trying to use using %d slices due to number of "
                      "threads (%d) being too small: results unpredictable",
                      pass_num_slice, num_threads);
    }
    initSlice(pass_num_slice);
  }
  multi_iteration = 0;
  //  string partitionFile = model->strOption[STROPT_PARTITION_FILE];
  //  if (partitionFile.size())
  //  {
  //    dualRHS.setup_partition(partitionFile.c_str());
  //  }
}

void HDual::initSlice(const int initial_num_slice) {
  // Number of slices
  slice_num = initial_num_slice;
  if (slice_num < 1) slice_num = 1;
  assert(slice_num <= HIGHS_SLICED_LIMIT);
  if (slice_num > HIGHS_SLICED_LIMIT) {
#ifdef HiGHSDEV
    printf(
        "WARNING: %d = slice_num > HIGHS_SLICED_LIMIT = %d so truncating "
        "slice_num\n",
        slice_num, HIGHS_SLICED_LIMIT);
#endif
    slice_num = HIGHS_SLICED_LIMIT;
  }

  // Alias to the matrix
  const int* Astart = matrix->getAstart();
  const int* Aindex = matrix->getAindex();
  const double* Avalue = matrix->getAvalue();
  const int AcountX = Astart[solver_num_col];

  // Figure out partition weight
  double sliced_countX = AcountX / slice_num;
  slice_start[0] = 0;
  for (int i = 0; i < slice_num - 1; i++) {
    int endColumn = slice_start[i] + 1;  // At least one column
    int endX = Astart[endColumn];
    int stopX = (i + 1) * sliced_countX;
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
  vector<int> sliced_Astart;
  for (int i = 0; i < slice_num; i++) {
    // The matrix
    int mystart = slice_start[i];
    int mycount = slice_start[i + 1] - mystart;
    int mystartX = Astart[mystart];
    sliced_Astart.resize(mycount + 1);
    for (int k = 0; k <= mycount; k++)
      sliced_Astart[k] = Astart[k + mystart] - mystartX;
    slice_matrix[i].setup_lgBs(mycount, solver_num_row, &sliced_Astart[0],
                               Aindex + mystartX, Avalue + mystartX);

    // The row_ap and its packages
    slice_row_ap[i].setup(mycount);
    slice_dualRow[i].setupSlice(mycount);
  }
}

void HDual::solvePhase1() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  // When starting a new phase the (updated) dual objective function
  // value isn't known. Indicate this so that when the value computed
  // from scratch in build() isn't checked against the the updated
  // value
  simplex_lp_status.has_dual_objective_value = 0;
  // Set invertHint so that it's assigned when first tested
  invertHint = INVERT_HINT_NO;
  // Set solvePhase=1 so it's set if solvePhase1() is called directly
  solvePhase = 1;
  // Report the phase start
  HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                    ML_DETAILED, "dual-phase-1-start\n");
  // Switch to dual phase 1 bounds
  initialise_bound(workHMO, 1);
  initialise_value(workHMO);
  // Main solving structure
  analysis->simplexTimerStart(IterateClock);
  for (;;) {
    analysis->simplexTimerStart(IterateDualRebuildClock);
    rebuild();
    analysis->simplexTimerStop(IterateDualRebuildClock);
    for (;;) {
      switch (simplex_info.simplex_strategy) {
        default:
        case SIMPLEX_STRATEGY_DUAL_PLAIN:
          iterate();
          break;
        case SIMPLEX_STRATEGY_DUAL_TASKS:
          iterateTasks();
          break;
        case SIMPLEX_STRATEGY_DUAL_MULTI:
          iterateMulti();
          break;
      }
      if (bailout()) break;
      if (invertHint) break;
    }
    if (bailout()) break;
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (simplex_info.update_count == 0) break;
    if (simplex_lp_status.has_fresh_rebuild) break;
  }

  analysis->simplexTimerStop(IterateClock);
  // Possibly return due to bailing out, having now stopped
  // IterateClock
  if (bailout()) return;

  // If bailing out, should have done so already
  assert(!solve_bailout);
  if (rowOut == -1) {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                      ML_DETAILED, "dual-phase-1-optimal\n");
    // Go to phase 2
    if (simplex_info.dual_objective_value == 0) {
      solvePhase = 2;
    } else {
      // We still have dual infeasible
      if (workHMO.simplex_info_.costs_perturbed) {
        // Clean up perturbation and go on
        cleanup();
        if (dualInfeasCount == 0) solvePhase = 2;
      } else {
        // Report dual infeasible
        solvePhase = -1;
        HighsPrintMessage(workHMO.options_.output,
                          workHMO.options_.message_level, ML_MINIMAL,
                          "dual-infeasible\n");
        workHMO.scaled_model_status_ = HighsModelStatus::PRIMAL_UNBOUNDED;
      }
    }
  } else if (invertHint == INVERT_HINT_CHOOSE_COLUMN_FAIL) {
    // chooseColumn has failed
    // Behave as "Report strange issues" below
    solvePhase = -1;
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                      ML_MINIMAL, "dual-phase-1-not-solved\n");
    workHMO.scaled_model_status_ = HighsModelStatus::SOLVE_ERROR;
  } else if (columnIn == -1) {
    // We got dual phase 1 unbounded - strange
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                      ML_MINIMAL, "dual-phase-1-unbounded\n");
    if (workHMO.simplex_info_.costs_perturbed) {
      // Clean up perturbation and go on
      cleanup();
      if (dualInfeasCount == 0) solvePhase = 2;
    } else {
      // Report strange issues
      solvePhase = -1;
      HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                        ML_MINIMAL, "dual-phase-1-not-solved\n");
      workHMO.scaled_model_status_ = HighsModelStatus::SOLVE_ERROR;
    }
  }

  if (solvePhase == 2) {
    initialise_bound(workHMO);
    initialise_value(workHMO);
  }
  return;
}

void HDual::solvePhase2() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  // When starting a new phase the (updated) dual objective function
  // value isn't known. Indicate this so that when the value computed
  // from scratch in build() isn't checked against the the updated
  // value
  simplex_lp_status.has_dual_objective_value = 0;
  // Set invertHint so that it's assigned when first tested
  invertHint = INVERT_HINT_NO;
  // Set solvePhase=2 so it's set if solvePhase2() is called directly
  solvePhase = 2;
  // Report the phase start
  HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                    ML_DETAILED, "dual-phase-2-start\n");
  // Collect free variables
  dualRow.createFreelist();
  // Main solving structure
  analysis->simplexTimerStart(IterateClock);
  for (;;) {
    // Outer loop of solvePhase2()
    // Rebuild all values, reinverting B if updates have been performed
    analysis->simplexTimerStart(IterateDualRebuildClock);
    rebuild();
    analysis->simplexTimerStop(IterateDualRebuildClock);
    if (dualInfeasCount > 0) break;
    for (;;) {
      // Inner loop of solvePhase2()
      // Performs one iteration in case SIMPLEX_STRATEGY_DUAL_PLAIN:
      switch (simplex_info.simplex_strategy) {
        default:
        case SIMPLEX_STRATEGY_DUAL_PLAIN:
          iterate();
          break;
        case SIMPLEX_STRATEGY_DUAL_TASKS:
          iterateTasks();
          break;
        case SIMPLEX_STRATEGY_DUAL_MULTI:
          iterateMulti();
          break;
      }
      if (bailout()) break;
      if (invertHint) break;
    }
    if (bailout()) break;
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (simplex_info.update_count == 0) break;
    if (simplex_lp_status.has_fresh_rebuild) break;
  }
  analysis->simplexTimerStop(IterateClock);
  // Possibly return due to bailing out, having now stopped
  // IterateClock
  if (bailout()) return;

  // If bailing out, should have done so already
  assert(!solve_bailout);
  if (dualInfeasCount > 0) {
    // There are dual infeasiblities so switch to Phase 1 and return
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                      ML_DETAILED, "dual-phase-2-found-free\n");
    solvePhase = 1;
  } else if (rowOut == -1) {
    // There is no candidate in CHUZR, even after rebuild so probably optimal
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                      ML_DETAILED, "dual-phase-2-optimal\n");
    // Remove any cost perturbations and see if basis is still dual feasible
    cleanup();
    if (dualInfeasCount > 0) {
      // There are dual infeasiblities, so consider performing primal
      // simplex iterations to get dual feasibility
      solvePhase = 4;
    } else {
      // There are no dual infeasiblities so optimal!
      solvePhase = 0;
      HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                        ML_DETAILED, "problem-optimal\n");
      workHMO.scaled_model_status_ = HighsModelStatus::OPTIMAL;
    }
  } else if (invertHint == INVERT_HINT_CHOOSE_COLUMN_FAIL) {
    // chooseColumn has failed
    // Behave as "Report strange issues" below
    solvePhase = -1;
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                      ML_MINIMAL, "dual-phase-2-not-solved\n");
    workHMO.scaled_model_status_ = HighsModelStatus::SOLVE_ERROR;
  } else if (columnIn == -1) {
    // There is no candidate in CHUZC, so probably dual unbounded
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                      ML_MINIMAL, "dual-phase-2-unbounded\n");
    if (workHMO.simplex_info_.costs_perturbed) {
      // If the costs have been perturbed, clean up and return
      cleanup();
    } else {
      // If the costs have not been perturbed, so dual unbounded---and hence
      // primal infeasible
      solvePhase = -1;
      HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                        ML_MINIMAL, "problem-infeasible\n");
      workHMO.scaled_model_status_ = HighsModelStatus::PRIMAL_INFEASIBLE;
    }
  }
  return;
}

void HDual::rebuild() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  // Save history information
  // Move this to Simplex class once it's created
  //  record_pivots(-1, -1, 0);  // Indicate REINVERT

  const int rebuild_invert_hint = invertHint;
  invertHint = INVERT_HINT_NO;
  // Possibly Rebuild workHMO.factor_
  bool reInvert = simplex_info.update_count > 0;
  if (!invert_if_row_out_negative) {
    // Don't reinvert if rowOut is negative [equivalently, if
    // rebuild_invert_hint == INVERT_HINT_POSSIBLY_OPTIMAL]
    if (rebuild_invert_hint == INVERT_HINT_POSSIBLY_OPTIMAL) {
      assert(rowOut == -1);
      reInvert = false;
    }
  }
  if (reInvert) {
    const int* baseIndex = &workHMO.simplex_basis_.basicIndex_[0];
    // Scatter the edge weights so that, after INVERT,
    // they can be gathered according to the new

    // permutation of baseIndex
    analysis->simplexTimerStart(PermWtClock);
    for (int i = 0; i < solver_num_row; i++)
      dualRHS.workEdWtFull[baseIndex[i]] = dualRHS.workEdWt[i];
    analysis->simplexTimerStop(PermWtClock);

    analysis->simplexTimerStart(InvertClock);

    // Call computeFactor to perform INVERT
    int rankDeficiency = computeFactor(workHMO);
    analysis->simplexTimerStop(InvertClock);

    if (rankDeficiency)
      throw runtime_error("Dual reInvert: singular-basis-matrix");
    // Gather the edge weights according to the
    // permutation of baseIndex after INVERT
    analysis->simplexTimerStart(PermWtClock);
    for (int i = 0; i < solver_num_row; i++)
      dualRHS.workEdWt[i] = dualRHS.workEdWtFull[baseIndex[i]];
    analysis->simplexTimerStop(PermWtClock);
  }

  // Recompute dual solution
  analysis->simplexTimerStart(ComputeDualClock);
  computeDual(workHMO);
  analysis->simplexTimerStop(ComputeDualClock);

  analysis->simplexTimerStart(CorrectDualClock);
  correctDual(workHMO, &dualInfeasCount);
  analysis->simplexTimerStop(CorrectDualClock);

  // Recompute primal solution
  analysis->simplexTimerStart(ComputePrimalClock);
  computePrimal(workHMO);
  analysis->simplexTimerStop(ComputePrimalClock);

  // Collect primal infeasible as a list
  analysis->simplexTimerStart(CollectPrIfsClock);
  dualRHS.createArrayOfPrimalInfeasibilities();
  dualRHS.createInfeasList(analysis->col_aq_density);
  analysis->simplexTimerStop(CollectPrIfsClock);

  analysis->simplexTimerStart(ComputePrIfsClock);
  computePrimalInfeasible(workHMO);
  analysis->simplexTimerStop(ComputePrIfsClock);

  analysis->simplexTimerStart(ComputeDuIfsClock);
  computeDualInfeasible(workHMO);
  analysis->simplexTimerStop(ComputeDuIfsClock);

  // Compute the objective value
  analysis->simplexTimerStart(ComputeDuObjClock);
  computeDualObjectiveValue(workHMO, solvePhase);
  analysis->simplexTimerStop(ComputeDuObjClock);

  double dual_objective_value = simplex_info.dual_objective_value;
#ifdef HiGHSDEV
  // Check the objective value maintained by updating against the
  // value when computed exactly - so long as there is a value to
  // check against
  /*
  if (simplex_lp_status.has_dual_objective_value) {
    double absDualObjectiveError =
        fabs(simplex_info.updated_dual_objective_value - dual_objective_value);
    double rlvDualObjectiveError =
        absDualObjectiveError / max(1.0, fabs(dual_objective_value));
    // TODO Investigate these Dual objective value errors
    if (rlvDualObjectiveError >= 1e-8) {
      HighsLogMessage(workHMO.options_.logfile, HighsMessageType::WARNING,
    "Dual objective value error |rel| = %12g (%12g)", absDualObjectiveError,
  rlvDualObjectiveError);
    }
  }
  */
#endif
  simplex_info.updated_dual_objective_value = dual_objective_value;

#ifdef HiGHSDEV
  //  checkDualObjectiveValue("After computing dual objective value");
  //  printf("Checking INVERT in rebuild()\n"); workHMO.factor_.checkInvert();
#endif

  analysis->simplexTimerStart(ReportRebuildClock);
  reportRebuild(rebuild_invert_hint);
  analysis->simplexTimerStop(ReportRebuildClock);

  build_syntheticTick = factor->build_syntheticTick;
  total_syntheticTick = 0;

#ifdef HiGHSDEV
  if (simplex_info.analyse_rebuild_time) {
    int total_rebuilds = analysis->simplexTimerNumCall(IterateDualRebuildClock);
    double total_rebuild_time =
        analysis->simplexTimerRead(IterateDualRebuildClock);
    printf(
        "Dual  Ph%-2d rebuild %4d (%1d) on iteration %9d: Total rebuild time = "
        "%11.4g\n",
        solvePhase, total_rebuilds, rebuild_invert_hint,
        workHMO.scaled_solution_params_.simplex_iteration_count,
        total_rebuild_time);
  }
#endif
  // Data are fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = true;
}

void HDual::cleanup() {
  HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                    ML_DETAILED, "dual-cleanup-shift\n");
  // Remove perturbation
  initialise_cost(workHMO);
  initialise_bound(workHMO);
  // Compute the dual values
#ifdef HiGHSDEV
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  vector<double> original_workDual = simplex_info.workDual_;
#endif
  analysis->simplexTimerStart(ComputeDualClock);
  computeDual(workHMO);
  analysis->simplexTimerStop(ComputeDualClock);
#ifdef HiGHSDEV
  int num_dual_sign_change = 0;
  for (int iCol = 0; iCol < workHMO.simplex_lp_.numCol_; iCol++) {
    const double max_dual =
        max(fabs(simplex_info.workDual_[iCol]), fabs(original_workDual[iCol]));
    if (max_dual > workHMO.scaled_solution_params_.dual_feasibility_tolerance) {
      if (simplex_info.workDual_[iCol] * original_workDual[iCol] < 0)
        num_dual_sign_change++;
    }
    const double dual_change =
        fabs(simplex_info.workDual_[iCol] - original_workDual[iCol]);
    updateValueDistribution(dual_change,
                            analysis->cleanup_dual_change_distribution);
  }
  printf(
      "grep_DuPtrb: dualCleanup for %s has %d meaningful dual sign change(s)\n",
      workHMO.simplex_lp_.model_name_.c_str(), num_dual_sign_change);
#endif

  // Compute the dual infeasibilities
  analysis->simplexTimerStart(ComputeDuIfsClock);
  computeDualInfeasible(workHMO);
  analysis->simplexTimerStop(ComputeDuIfsClock);

  // Compute the dual objective value
  analysis->simplexTimerStart(ComputeDuObjClock);
  computeDualObjectiveValue(workHMO, solvePhase);
  analysis->simplexTimerStop(ComputeDuObjClock);

  analysis->simplexTimerStart(ReportRebuildClock);
  reportRebuild();
  analysis->simplexTimerStop(ReportRebuildClock);

  computeDualInfeasible(workHMO);
  dualInfeasCount = workHMO.scaled_solution_params_.num_dual_infeasibilities;
}

void HDual::iterate() {
  // This is the main teration loop for dual revised simplex. All the
  // methods have as their first line if (invertHint) return;, where
  // invertHint is, for example, set to 1 when CHUZR finds no
  // candidate. This causes a break from the inner loop of
  // solve_phase% and, hence, a call to rebuild()

  // Reporting:
  // Row-wise matrix after update in updateMatrix(columnIn, columnOut);
  analysis->simplexTimerStart(IterateChuzrClock);
  chooseRow();
  analysis->simplexTimerStop(IterateChuzrClock);

  analysis->simplexTimerStart(IterateChuzcClock);
  chooseColumn(&row_ep);
  analysis->simplexTimerStop(IterateChuzcClock);

#ifdef HiGHSDEV
  if (rp_iter_da && rowOut >= 0) {
    // for (int row=0; row < workHMO.lp_.numRow_; row++) printf("Row %2d: Devex
    // Weight = %11.4g\n", row, dualRHS.workEdWt[row]);
    printf(
        "Iter %4d: rowOut %4d; colOut %4d; colIn %4d; Wt = %11.4g; thetaDual = "
        "%11.4g; alpha = %11.4g; Dvx = %d\n",
        workHMO.scaled_solution_params_.simplex_iteration_count, rowOut,
        columnOut, columnIn, computed_edge_weight, thetaDual, alphaRow,
        num_devex_iterations);
  }
#endif

  analysis->simplexTimerStart(IterateFtranClock);
  updateFtranBFRT();

  // updateFtran(); computes the pivotal column in the data structure "column"
  updateFtran();

  // updateFtranDSE performs the DSE FTRAN on pi_p
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE)
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

  // updatePrimal(&row_ep); Updates the primal values and the edge weights
  analysis->simplexTimerStart(IteratePrimalClock);
  updatePrimal(&row_ep);
  analysis->simplexTimerStop(IteratePrimalClock);
  // After primal update in dual simplex the primal objective value is not known
  workHMO.simplex_lp_status_.has_primal_objective_value = false;

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

void HDual::iterateTasks() {
  slice_PRICE = 1;

  // Group 1
  chooseRow();

  // Disable slice when too sparse
  if (1.0 * row_ep.count / solver_num_row < 0.01) slice_PRICE = 0;

  analysis->simplexTimerStart(Group1Clock);
#pragma omp parallel
#pragma omp single
  {
#pragma omp task
    {
      col_DSE.copy(&row_ep);
      updateFtranDSE(&col_DSE);
    }
#pragma omp task
    {
      if (slice_PRICE)
        chooseColumnSlice(&row_ep);
      else
        chooseColumn(&row_ep);
#pragma omp task
      updateFtranBFRT();
#pragma omp task
      updateFtran();
#pragma omp taskwait
    }
  }
  analysis->simplexTimerStop(Group1Clock);

  updateVerify();
  updateDual();
  updatePrimal(&col_DSE);
  updatePivots();
}

void HDual::iterationAnalysisData() {
  HighsSolutionParams& scaled_solution_params = workHMO.scaled_solution_params_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  analysis->simplex_strategy = simplex_info.simplex_strategy;
  analysis->edge_weight_mode = dual_edge_weight_mode;
  analysis->solve_phase = solvePhase;
  analysis->simplex_iteration_count =
      scaled_solution_params.simplex_iteration_count;
  analysis->devex_iteration_count = num_devex_iterations;
  analysis->pivotal_row_index = rowOut;
  analysis->leaving_variable = columnOut;
  analysis->entering_variable = columnIn;
  analysis->invert_hint = invertHint;
  analysis->freelist_size = dualRow.freeListSize;
  analysis->reduced_rhs_value = 0;
  analysis->reduced_cost_value = 0;
  analysis->edge_weight = 0;
  analysis->primal_delta = deltaPrimal;
  analysis->primal_step = thetaPrimal;
  analysis->dual_step = thetaDual;
  analysis->pivot_value_from_column = alpha;
  analysis->pivot_value_from_row = alphaRow;
  analysis->numerical_trouble = numericalTrouble;
  analysis->objective_value = simplex_info.updated_dual_objective_value;
  analysis->num_primal_infeasibilities =
      scaled_solution_params.num_primal_infeasibilities;
  analysis->num_dual_infeasibilities =
      scaled_solution_params.num_dual_infeasibilities;
  analysis->sum_primal_infeasibilities =
      scaled_solution_params.sum_primal_infeasibilities;
  analysis->sum_dual_infeasibilities =
      scaled_solution_params.sum_dual_infeasibilities;
#ifdef HiGHSDEV
  analysis->basis_condition = simplex_info.invert_condition;
#endif
  if ((dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) &&
      (num_devex_iterations == 0))
    analysis->num_devex_framework++;
}

void HDual::iterationAnalysis() {
  // Possibly report on the iteration
  iterationAnalysisData();
  analysis->iterationReport();

  // Possibly switch from DSE to Devex
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    bool switch_to_devex = false;
    switch_to_devex = analysis->switchToDevex();
    if (switch_to_devex) {
      dual_edge_weight_mode = DualEdgeWeightMode::DEVEX;
      // Using dual Devex edge weights, so set up the first framework
      workHMO.simplex_info_.devex_index_.assign(solver_num_tot, 0);
      initialiseDevexFramework();
    }
  }

#ifdef HiGHSDEV
  analysis->iterationRecord();
#endif
}

void HDual::reportRebuild(const int rebuild_invert_hint) {
  iterationAnalysisData();
  analysis->invert_hint = rebuild_invert_hint;
  analysis->invertReport();
}

void HDual::chooseRow() {
  // Choose the index of a row to leave the basis (CHUZR)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  // Choose candidates repeatedly until candidate is OK or optimality is
  // detected
  for (;;) {
    // Choose the index of a good row to leave the basis
    dualRHS.chooseNormal(&rowOut);
    if (rowOut == -1) {
      // No index found so may be dual optimal. By setting
      // invertHint>0 all subsequent methods in the iteration will
      // be skipped until reinversion and rebuild have taken place
      invertHint = INVERT_HINT_POSSIBLY_OPTIMAL;
      return;
    }
    // Compute pi_p = B^{-T}e_p in row_ep
    analysis->simplexTimerStart(BtranClock);
    // Set up RHS for BTRAN
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = rowOut;
    row_ep.array[rowOut] = 1;
    row_ep.packFlag = true;
#ifdef HiGHSDEV
    HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
    if (simplex_info.analyse_iterations)
      analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep,
                                      analysis->row_ep_density);
#endif
    // Perform BTRAN
    factor->btran(row_ep, analysis->row_ep_density,
                  analysis->pointer_serial_factor_clocks);
#ifdef HiGHSDEV
    if (simplex_info.analyse_iterations)
      analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep);
#endif
    analysis->simplexTimerStop(BtranClock);
    // Verify DSE weight
    if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
      // For DSE, see how accurate the updated weight is
      // Save the updated weight
      double updated_edge_weight = dualRHS.workEdWt[rowOut];
      // Compute the weight from row_ep and over-write the updated weight
      computed_edge_weight = dualRHS.workEdWt[rowOut] = row_ep.norm2();
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
  // Index of row to leave the basis has been found
  //
  // Assign basic info:
  //
  // Record the column (variable) associated with the leaving row
  columnOut = workHMO.simplex_basis_.basicIndex_[rowOut];
  // Record the change in primal variable associated with the move to the bound
  // being violated
  if (baseValue[rowOut] < baseLower[rowOut]) {
    // Below the lower bound so set deltaPrimal = value - LB < 0
    deltaPrimal = baseValue[rowOut] - baseLower[rowOut];
  } else {
    // Above the upper bound so set deltaPrimal = value - UB > 0
    deltaPrimal = baseValue[rowOut] - baseUpper[rowOut];
  }
  // Set sourceOut to be -1 if deltaPrimal<0, otherwise +1 (since deltaPrimal>0)
  sourceOut = deltaPrimal < 0 ? -1 : 1;
  // Update the record of average row_ep (pi_p) density. This ignores
  // any BTRANs done for skipped candidates
  const double local_row_ep_density = (double)row_ep.count / solver_num_row;
  analysis->updateOperationResultDensity(local_row_ep_density,
                                         analysis->row_ep_density);
}

bool HDual::acceptDualSteepestEdgeWeight(const double updated_edge_weight) {
  // Accept the updated weight if it is at least a quarter of the
  // computed weight. Excessively large updated weights don't matter!
  const double accept_weight_threshhold = 0.25;
  const bool accept_weight =
      updated_edge_weight >= accept_weight_threshhold * computed_edge_weight;
  analysis->dualSteepestEdgeWeightError(computed_edge_weight,
                                        updated_edge_weight);
  return accept_weight;
}

bool HDual::newDevexFramework(const double updated_edge_weight) {
  // Analyse the Devex weight to determine whether a new framework
  // should be set up
  double devex_ratio = max(updated_edge_weight / computed_edge_weight,
                           computed_edge_weight / updated_edge_weight);
  int i_te = solver_num_row / minRlvNumberDevexIterations;
  i_te = max(minAbsNumberDevexIterations, i_te);
  // Square maxAllowedDevexWeightRatio due to keeping squared
  // weights
  const double accept_ratio_threshhold =
      maxAllowedDevexWeightRatio * maxAllowedDevexWeightRatio;
  const bool accept_ratio = devex_ratio <= accept_ratio_threshhold;
  const bool accept_it = num_devex_iterations <= i_te;
  bool return_new_devex_framework;
  return_new_devex_framework = !accept_ratio || !accept_it;
  /*
  if (return_new_devex_framework) {
    printf("New Devex framework: (Iter %d) updated weight = %11.4g; computed
  weight = %11.4g; Devex ratio = %11.4g\n",
           workHMO.scaled_solution_params_.simplex_iteration_count,
           updated_edge_weight, computed_edge_weight, devex_ratio);
    return true;
  }
  */
  return return_new_devex_framework;
}

void HDual::chooseColumn(HVector* row_ep) {
  // Compute pivot row (PRICE) and choose the index of a column to enter the
  // basis (CHUZC)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  //
  // PRICE
  //
  computeTableauRowFromPiP(workHMO, *row_ep, row_ap);
  //
  // CHUZC
  //
  // Section 0: Clear data and call createFreemove to set a value of
  // nonbasicMove for all free columns to prevent their dual values
  // from being changed.
  analysis->simplexTimerStart(Chuzc0Clock);
  dualRow.clear();
  dualRow.workDelta = deltaPrimal;
  dualRow.createFreemove(row_ep);
  analysis->simplexTimerStop(Chuzc0Clock);
  //
  // Section 1: Pack row_ap and row_ep, then determine the possible
  // variables - candidates for CHUZC
  analysis->simplexTimerStart(Chuzc1Clock);
  // Pack row_ap into the packIndex/Value of HDualRow
  dualRow.chooseMakepack(&row_ap, 0);
  // Pack row_ep into the packIndex/Value of HDualRow
  dualRow.chooseMakepack(row_ep, solver_num_col);
  // Determine the possible variables - candidates for CHUZC
  dualRow.choosePossible();
  analysis->simplexTimerStop(Chuzc1Clock);
  //
  // Take action if the step to an expanded bound is not positive, or
  // there are no candidates for CHUZC
  columnIn = -1;
  if (dualRow.workTheta <= 0 || dualRow.workCount == 0) {
    invertHint = INVERT_HINT_POSSIBLY_DUAL_UNBOUNDED;
    return;
  }
  //
  // Sections 2 and 3: Perform (bound-flipping) ratio test. This can
  // fail if the dual values are excessively large
  bool chooseColumnFail = dualRow.chooseFinal();
  if (chooseColumnFail) {
    invertHint = INVERT_HINT_CHOOSE_COLUMN_FAIL;
    return;
  }
  //
  // Section 4: Reset the nonbasicMove values for free columns
  analysis->simplexTimerStart(Chuzc4Clock);
  dualRow.deleteFreemove();
  analysis->simplexTimerStop(Chuzc4Clock);
  // Record values for basis change, checking for numerical problems and update
  // of dual variables
  columnIn = dualRow.workPivot;   // Index of the column entering the basis
  alphaRow = dualRow.workAlpha;   // Pivot value computed row-wise - used for
                                  // numerical checking
  thetaDual = dualRow.workTheta;  // Dual step length

  if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX &&
      !new_devex_framework) {
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

void HDual::chooseColumnSlice(HVector* row_ep) {
  // Choose the index of a column to enter the basis (CHUZC) by
  // exploiting slices of the pivotal row - for SIP and PAMI
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;

  analysis->simplexTimerStart(Chuzc0Clock);
  dualRow.clear();
  dualRow.workDelta = deltaPrimal;
  dualRow.createFreemove(row_ep);
  analysis->simplexTimerStop(Chuzc0Clock);

  //  const int solver_num_row = highs_model_object.simplex_lp_.numRow_;
  const double local_density = 1.0 * row_ep->count / solver_num_row;
  bool use_col_price;
  bool use_row_price_w_switch;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  choosePriceTechnique(simplex_info.price_strategy, local_density,
                       use_col_price, use_row_price_w_switch);

#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations) {
    const int row_ep_count = row_ep->count;
    if (use_col_price) {
      analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP,
                                      row_ep_count, 0.0);
      analysis->num_col_price++;
    } else if (use_row_price_w_switch) {
      analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP,
                                      row_ep_count, analysis->row_ep_density);
      analysis->num_row_price_with_switch++;
    } else {
      analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP,
                                      row_ep_count, analysis->row_ep_density);
      analysis->num_row_price++;
    }
  }
#endif
  analysis->simplexTimerStart(PriceChuzc1Clock);
  // Row_ep:         PACK + CC1

  /*
  int row_ep_thread_id = 0;
  vector<int> row_ap_thread_id;
  row_ap_thread_id.resize(slice_num);
  */

#pragma omp task
  {
    dualRow.chooseMakepack(row_ep, solver_num_col);
    dualRow.choosePossible();
#ifdef OPENMP
    //    int row_ep_thread_id = omp_get_thread_num();
    //    printf("Hello world from Row_ep:         PACK + CC1 thread %d\n",
    //    row_ep_thread_id);
#endif
  }

  // Row_ap: PRICE + PACK + CC1
  for (int i = 0; i < slice_num; i++) {
#pragma omp task
    {
#ifdef OPENMP
      //      int row_ap_thread_id = omp_get_thread_num();
      //      printf("Hello world from omp Row_ap: PRICE + PACK + CC1 [%1d]
      //      thread %d\n", i, row_ap_thread_id);
#endif
      slice_row_ap[i].clear();

      //      slice_matrix[i].priceByRowSparseResult(slice_row_ap[i], *row_ep);

      if (use_col_price) {
        // Perform column-wise PRICE
        slice_matrix[i].priceByColumn(slice_row_ap[i], *row_ep);
      } else if (use_row_price_w_switch) {
        // Perform hyper-sparse row-wise PRICE, but switch if the density of
        // row_ap becomes extreme
        slice_matrix[i].priceByRowSparseResultWithSwitch(
            slice_row_ap[i], *row_ep, analysis->row_ap_density, 0,
            slice_matrix[i].hyperPRICE);
      } else {
        // Perform hyper-sparse row-wise PRICE
        slice_matrix[i].priceByRowSparseResult(slice_row_ap[i], *row_ep);
      }

      slice_dualRow[i].clear();
      slice_dualRow[i].workDelta = deltaPrimal;
      slice_dualRow[i].chooseMakepack(&slice_row_ap[i], slice_start[i]);
      slice_dualRow[i].choosePossible();
    }
  }
#pragma omp taskwait

#ifdef HiGHSDEV
  // Determine the nonzero count of the whole row
  if (simplex_info.analyse_iterations) {
    int row_ap_count = 0;
    for (int i = 0; i < slice_num; i++) row_ap_count += slice_row_ap[i].count;
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_AP,
                                   row_ap_count);
  }
#endif

  // Join CC1 results here
  for (int i = 0; i < slice_num; i++) {
    dualRow.chooseJoinpack(&slice_dualRow[i]);
  }

  analysis->simplexTimerStop(PriceChuzc1Clock);

  // Infeasible we created before
  columnIn = -1;
  if (dualRow.workTheta <= 0 || dualRow.workCount == 0) {
    invertHint = INVERT_HINT_POSSIBLY_DUAL_UNBOUNDED;
    return;
  }

  // Choose column 2, This only happens if didn't go out
  bool chooseColumnFail = dualRow.chooseFinal();
  if (chooseColumnFail) {
    invertHint = INVERT_HINT_CHOOSE_COLUMN_FAIL;
    return;
  }

  analysis->simplexTimerStart(Chuzc4Clock);
  dualRow.deleteFreemove();
  analysis->simplexTimerStop(Chuzc4Clock);

  columnIn = dualRow.workPivot;
  alphaRow = dualRow.workAlpha;
  thetaDual = dualRow.workTheta;

  if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX &&
      !new_devex_framework) {
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
    for (int i = 0; i < slice_num; i++) slice_dualRow[i].computeDevexWeight(i);
    // Accumulate the partial sums
    // Initialse with the partial sum for row_ep
    computed_edge_weight = dualRow.computed_edge_weight;
    // Update with the partial sum for row_ep
    for (int i = 0; i < slice_num; i++)
      computed_edge_weight += slice_dualRow[i].computed_edge_weight;
    computed_edge_weight = max(1.0, computed_edge_weight);
    analysis->simplexTimerStop(DevexWtClock);
  }
}

void HDual::updateFtran() {
  // Compute the pivotal column (FTRAN)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  analysis->simplexTimerStart(FtranClock);
  // Clear the picotal column and indicate that its values should be packed
  col_aq.clear();
  col_aq.packFlag = true;
  // Get the constraint matrix column by combining just one column
  // with unit multiplier
  matrix->collect_aj(col_aq, columnIn, 1);
#ifdef HiGHSDEV
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq,
                                    analysis->col_aq_density);
#endif
  // Perform FTRAN
  factor->ftran(col_aq, analysis->col_aq_density,
                analysis->pointer_serial_factor_clocks);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq);
#endif
  const double local_col_aq_density = (double)col_aq.count / solver_num_row;
  analysis->updateOperationResultDensity(local_col_aq_density,
                                         analysis->col_aq_density);
  // Save the pivot value computed column-wise - used for numerical checking
  alpha = col_aq.array[rowOut];
  analysis->simplexTimerStop(FtranClock);
}

void HDual::updateFtranBFRT() {
  // Compute the RHS changes corresponding to the BFRT (FTRAN-BFRT)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;

  // Only time updateFtranBFRT if dualRow.workCount > 0;
  // If dualRow.workCount = 0 then dualRow.updateFlip(&col_BFRT)
  // merely clears col_BFRT so no FTRAN is performed
  bool time_updateFtranBFRT = dualRow.workCount > 0;

  if (time_updateFtranBFRT) {
    analysis->simplexTimerStart(FtranBfrtClock);
  }

  dualRow.updateFlip(&col_BFRT);

  if (col_BFRT.count) {
#ifdef HiGHSDEV
    HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
    if (simplex_info.analyse_iterations)
      analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_FTRAN_BFRT,
                                      col_BFRT, analysis->col_BFRT_density);
#endif
    // Perform FTRAN BFRT
    factor->ftran(col_BFRT, analysis->col_BFRT_density,
                  analysis->pointer_serial_factor_clocks);
#ifdef HiGHSDEV
    if (simplex_info.analyse_iterations)
      analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_FTRAN_BFRT,
                                     col_BFRT);
#endif
  }
  if (time_updateFtranBFRT) {
    analysis->simplexTimerStop(FtranBfrtClock);
  }
  const double local_col_BFRT_density = (double)col_BFRT.count / solver_num_row;
  analysis->updateOperationResultDensity(local_col_BFRT_density,
                                         analysis->col_BFRT_density);
}

void HDual::updateFtranDSE(HVector* DSE_Vector) {
  // Compute the vector required to update DSE weights - being FTRAN
  // applied to the pivotal column (FTRAN-DSE)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  analysis->simplexTimerStart(FtranDseClock);
#ifdef HiGHSDEV
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_FTRAN_DSE,
                                    *DSE_Vector, analysis->row_DSE_density);
#endif
  // Perform FTRAN DSE
  factor->ftran(*DSE_Vector, analysis->row_DSE_density,
                analysis->pointer_serial_factor_clocks);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_FTRAN_DSE,
                                   *DSE_Vector);
#endif
  analysis->simplexTimerStop(FtranDseClock);
  const double local_row_DSE_density =
      (double)DSE_Vector->count / solver_num_row;
  analysis->updateOperationResultDensity(local_row_DSE_density,
                                         analysis->row_DSE_density);
}

void HDual::updateVerify() {
  // Compare the pivot value computed row-wise and column-wise and
  // determine whether reinversion is advisable
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;

  // Use the two pivot values to identify numerical trouble
  if (reinvertOnNumericalTrouble("HDual::updateVerify", workHMO,
                                 numericalTrouble, alpha, alphaRow,
                                 numerical_trouble_tolerance)) {
    invertHint = INVERT_HINT_POSSIBLY_SINGULAR_BASIS;
  }
}

void HDual::updateDual() {
  // Update the dual values
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;

  // Update - dual (shift and back)
  if (thetaDual == 0)
    // Little to do if thetaDual is zero
    shift_cost(workHMO, columnIn, -workDual[columnIn]);
  else {
    // Update the whole vector of dual values
    dualRow.updateDual(thetaDual);
    if (workHMO.simplex_info_.simplex_strategy != SIMPLEX_STRATEGY_DUAL_PLAIN &&
        slice_PRICE) {
      // Update the slice-by-slice copy of dual variables
      for (int i = 0; i < slice_num; i++)
        slice_dualRow[i].updateDual(thetaDual);
    }
  }
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;
  shift_back(workHMO, columnOut);
}

void HDual::updatePrimal(HVector* DSE_Vector) {
  // Update the primal values and any edge weights
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) {
    const double updated_edge_weight = dualRHS.workEdWt[rowOut];
    dualRHS.workEdWt[rowOut] = computed_edge_weight;
    new_devex_framework = newDevexFramework(updated_edge_weight);
  }
  // DSE_Vector is either col_DSE = B^{-1}B^{-T}e_p (if using dual
  // steepest edge weights) or row_ep = B^{-T}e_p.
  //
  // Update - primal and weight
  dualRHS.updatePrimal(&col_BFRT, 1);
  dualRHS.updateInfeasList(&col_BFRT);
  double x_out = baseValue[rowOut];
  double l_out = baseLower[rowOut];
  double u_out = baseUpper[rowOut];
  thetaPrimal = (x_out - (deltaPrimal < 0 ? l_out : u_out)) / alpha;
  dualRHS.updatePrimal(&col_aq, thetaPrimal);
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    const double new_pivotal_edge_weight =
        dualRHS.workEdWt[rowOut] / (alpha * alpha);
    const double Kai = -2 / alpha;
    dualRHS.updateWeightDualSteepestEdge(&col_aq, new_pivotal_edge_weight, Kai,
                                         &DSE_Vector->array[0]);
    dualRHS.workEdWt[rowOut] = new_pivotal_edge_weight;
  } else if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) {
    // Pivotal row is for the current basis: weights are required for
    // the next basis so have to divide the current (exact) weight by
    // the pivotal value
    double new_pivotal_edge_weight = dualRHS.workEdWt[rowOut] / (alpha * alpha);
    new_pivotal_edge_weight = max(1.0, new_pivotal_edge_weight);
    // nw_wt is max(workEdWt[iRow], NewExactWeight*columnArray[iRow]^2);
    //
    // But NewExactWeight is new_pivotal_edge_weight = max(1.0,
    // dualRHS.workEdWt[rowOut] / (alpha * alpha))
    //
    // so nw_wt = max(workEdWt[iRow],
    // new_pivotal_edge_weight*columnArray[iRow]^2);
    //
    // Update rest of weights
    dualRHS.updateWeightDevex(&col_aq, new_pivotal_edge_weight);
    dualRHS.workEdWt[rowOut] = new_pivotal_edge_weight;
    num_devex_iterations++;
  }
  dualRHS.updateInfeasList(&col_aq);

  /*
  // Move these to where the FTRAN operations are actually performed
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    const double local_row_DSE_density = (double)DSE_Vector->count /
  solver_num_row; analysis->updateOperationResultDensity(local_row_DSE_density,
  analysis->row_DSE_density);
  }
  const double local_col_aq_density = (double)col_aq.count / solver_num_row;
  analysis->updateOperationResultDensity(local_col_aq_density,
  analysis->col_aq_density);
  */
  // Whether or not dual steepest edge weights are being used, have to
  // add in DSE_Vector->syntheticTick since this contains the
  // contribution from forming row_ep = B^{-T}e_p.
  total_syntheticTick += col_aq.syntheticTick;
  total_syntheticTick += DSE_Vector->syntheticTick;
}

void HDual::updatePivots() {
  // UPDATE
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  //
  // Update the sets of indices of basic and nonbasic variables
  update_pivots(workHMO, columnIn, rowOut, sourceOut);
  //  checkDualObjectiveValue("After update_pivots");
  //
  // Update the iteration count and store the basis change if HiGHSDEV
  // is defined
  // Move this to Simplex class once it's created
  // simplex_method.record_pivots(columnIn, columnOut, alpha);
  workHMO.scaled_solution_params_.simplex_iteration_count++;
  //
  // Update the invertible representation of the basis matrix
  update_factor(workHMO, &col_aq, &row_ep, &rowOut, &invertHint);
  //
  // Update the row-wise representation of the nonbasic columns
  update_matrix(workHMO, columnIn, columnOut);
  //
  // Delete Freelist entry for columnIn
  dualRow.deleteFreelist(columnIn);
  //
  // Update the primal value for the row where the basis change has
  // occurred, and set the corresponding primal infeasibility value in
  // dualRHS.work_infeasibility
  dualRHS.updatePivots(
      rowOut, workHMO.simplex_info_.workValue_[columnIn] + thetaPrimal);
  // Determine whether to reinvert based on the synthetic clock
  bool reinvert_syntheticClock = total_syntheticTick >= build_syntheticTick;
  const bool performed_min_updates =
      workHMO.simplex_info_.update_count >=
      synthetic_tick_reinversion_min_update_count;
#ifdef HiGHSDEV
  if (rp_reinvert_syntheticClock)
    printf(
        "Synth Reinversion: total_syntheticTick = %11.4g >=? %11.4g = "
        "build_syntheticTick: (%1d, %4d)\n",
        total_syntheticTick, build_syntheticTick, reinvert_syntheticClock,
        workHMO.simplex_info_.update_count);
#endif
  if (reinvert_syntheticClock && performed_min_updates)
    invertHint = INVERT_HINT_SYNTHETIC_CLOCK_SAYS_INVERT;
}

void HDual::initialiseDevexFramework(const bool parallel) {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  // Initialise the Devex framework: reference set is all basic
  // variables
  analysis->simplexTimerStart(DevexIzClock);
  const vector<int>& nonbasicFlag = workHMO.simplex_basis_.nonbasicFlag_;
  // Initialise the devex framework. The devex reference set is
  // initialise to be the current set of basic variables - and never
  // changes until a new framework is set up. In a simplex iteration,
  // to compute the exact Devex weight for the pivotal row requires
  // summing the squares of the its entries over the indices in the
  // reference set. This is achieved by summing over all indices, but
  // multiplying the entry by the value in devex_index before
  // equaring. Thus devex_index contains 1 for indices in the
  // reference set, and 0 otherwise. This is achieved by setting the
  // values of devex_index to be 1-nonbasicFlag^2, ASSUMING
  // |nonbasicFlag|=1 iff the corresponding variable is nonbasic
  for (int vr_n = 0; vr_n < solver_num_tot; vr_n++)
    simplex_info.devex_index_[vr_n] =
        1 - nonbasicFlag[vr_n] * nonbasicFlag[vr_n];
  // Set all initial weights to 1, zero the count of iterations with
  // this Devex framework, increment the number of Devex frameworks
  // and indicate that there's no need for a new Devex framework
  dualRHS.workEdWt.assign(solver_num_row, 1.0);
  num_devex_iterations = 0;
  new_devex_framework = false;
  minor_new_devex_framework = false;
  analysis->simplexTimerStop(DevexIzClock);
}

void HDual::interpretDualEdgeWeightStrategy(
    const int dual_edge_weight_strategy) {
  if (dual_edge_weight_strategy == SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG) {
    dual_edge_weight_mode = DualEdgeWeightMode::DANTZIG;
  } else if (dual_edge_weight_strategy ==
             SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX) {
    dual_edge_weight_mode = DualEdgeWeightMode::DEVEX;
  } else if (dual_edge_weight_strategy ==
             SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE) {
    dual_edge_weight_mode = DualEdgeWeightMode::STEEPEST_EDGE;
    initialise_dual_steepest_edge_weights = true;
    allow_dual_steepest_edge_to_devex_switch = false;
  } else if (dual_edge_weight_strategy ==
             SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL) {
    dual_edge_weight_mode = DualEdgeWeightMode::STEEPEST_EDGE;
    initialise_dual_steepest_edge_weights = false;
    allow_dual_steepest_edge_to_devex_switch = false;
  } else if (dual_edge_weight_strategy ==
             SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_TO_DEVEX_SWITCH) {
    dual_edge_weight_mode = DualEdgeWeightMode::STEEPEST_EDGE;
    initialise_dual_steepest_edge_weights = true;
    allow_dual_steepest_edge_to_devex_switch = true;
  } else {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                      ML_MINIMAL,
                      "HDual::interpretDualEdgeWeightStrategy: "
                      "unrecognised dual_edge_weight_strategy = %d - using "
                      "dual steepest edge with possible switch to Devex\n",
                      dual_edge_weight_strategy);
    dual_edge_weight_mode = DualEdgeWeightMode::STEEPEST_EDGE;
    initialise_dual_steepest_edge_weights = true;
    allow_dual_steepest_edge_to_devex_switch = true;
  }
}

/*
void HDual::interpretPriceStrategy(const int price_strategy) {
  allow_priceByColumn_switch = false;
  allow_priceByRowSparseResult_switch = false;
  if (price_strategy == SIMPLEX_PRICE_STRATEGY_COL) {
    price_mode = PriceMode::COL;
  } else if (price_strategy == SIMPLEX_PRICE_STRATEGY_ROW) {
    price_mode = PriceMode::ROW;
  } else if (price_strategy == SIMPLEX_PRICE_STRATEGY_ROW_SWITCH) {
    price_mode = PriceMode::ROW;
    allow_priceByRowSparseResult_switch = true;
  } else if (price_strategy == SIMPLEX_PRICE_STRATEGY_ROW_SWITCH_COL_SWITCH) {
    price_mode = PriceMode::ROW;
    allow_priceByColumn_switch = true;
    allow_priceByRowSparseResult_switch = true;
  } else {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
ML_MINIMAL, "HDual::interpretPriceStrategy: unrecognised price_strategy = %d - "
                      "using row Price with switch or colump price switch\n",
        price_strategy);
    price_mode = PriceMode::ROW;
    allow_priceByColumn_switch = true;
    allow_priceByRowSparseResult_switch = true;
  }
}
*/
bool HDual::dualInfoOk(const HighsLp& lp) {
  int lp_numCol = lp.numCol_;
  int lp_numRow = lp.numRow_;
  bool dimensions_ok;
  dimensions_ok = lp_numCol == solver_num_col && lp_numRow == solver_num_row;
  assert(dimensions_ok);
  if (!dimensions_ok) {
    printf("LP-Solver dimension incompatibility (%d, %d) != (%d, %d)\n",
           lp_numCol, solver_num_col, lp_numRow, solver_num_row);
    return false;
  }
  dimensions_ok = lp_numCol == factor->numCol && lp_numRow == factor->numRow;
  assert(dimensions_ok);
  if (!dimensions_ok) {
    printf("LP-Factor dimension incompatibility (%d, %d) != (%d, %d)\n",
           lp_numCol, factor->numCol, lp_numRow, factor->numRow);
    return false;
  }
  return true;
}

bool HDual::bailout() {
  if (solve_bailout) {
    // Bailout has already been decided: check that it's for one of these
    // reasons
    assert(workHMO.scaled_model_status_ ==
               HighsModelStatus::REACHED_TIME_LIMIT ||
           workHMO.scaled_model_status_ ==
               HighsModelStatus::REACHED_ITERATION_LIMIT ||
           workHMO.scaled_model_status_ ==
               HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND);
  } else if (workHMO.timer_.readRunHighsClock() > workHMO.options_.time_limit) {
    solve_bailout = true;
    workHMO.scaled_model_status_ = HighsModelStatus::REACHED_TIME_LIMIT;
  } else if (workHMO.scaled_solution_params_.simplex_iteration_count >=
             workHMO.options_.simplex_iteration_limit) {
    solve_bailout = true;
    workHMO.scaled_model_status_ = HighsModelStatus::REACHED_ITERATION_LIMIT;
  } else if (solvePhase == 2 &&
             (workHMO.simplex_info_.updated_dual_objective_value >
              workHMO.options_.dual_objective_value_upper_bound)) {
#ifdef SCIP_DEV
    printf("HDual::solvePhase2: %12g = Objective > ObjectiveUB\n",
           workHMO.simplex_info_.updated_dual_objective_value,
           workHMO.options_.dual_objective_value_upper_bound);
#endif
    solve_bailout = true;
    workHMO.scaled_model_status_ =
        HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND;
  }
  return solve_bailout;
}

#ifdef HiGHSDEV
double HDual::checkDualObjectiveValue(const char* message, int phase) {
  static double previous_updated_dual_objective_value = 0;
  static double previous_dual_objective_value = 0;
  computeDualObjectiveValue(workHMO, phase);
  double updated_dual_objective_value =
      workHMO.simplex_info_.updated_dual_objective_value;
  double dual_objective_value = workHMO.simplex_info_.dual_objective_value;
  double change_in_updated_dual_objective_value =
      updated_dual_objective_value - previous_updated_dual_objective_value;
  double change_in_dual_objective_value =
      dual_objective_value - previous_dual_objective_value;
  double updated_dual_objective_error =
      dual_objective_value - updated_dual_objective_value;
  double relative_updated_dual_objective_error =
      fabs(updated_dual_objective_error) / max(1.0, fabs(dual_objective_value));
  bool error_found = relative_updated_dual_objective_error > 1e-8;
  if (error_found)
    printf(
        "Phase %1d: duObjV = %11.4g (%11.4g); updated duObjV = %11.4g "
        "(%11.4g); Error(|Rel|) = %11.4g (%11.4g) |%s\n",
        phase, dual_objective_value, change_in_dual_objective_value,
        updated_dual_objective_value, change_in_updated_dual_objective_value,
        updated_dual_objective_error, relative_updated_dual_objective_error,
        message);
  previous_dual_objective_value = dual_objective_value;
  previous_updated_dual_objective_value = dual_objective_value;
  workHMO.simplex_info_.updated_dual_objective_value = dual_objective_value;
  // Now have dual objective value
  workHMO.simplex_lp_status_.has_dual_objective_value = true;
  return updated_dual_objective_error;
}
#endif
