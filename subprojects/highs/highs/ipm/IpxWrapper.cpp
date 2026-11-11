/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ipm/IpxWrapper.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova and Michael Feldmeier
 */
#include "ipm/IpxWrapper.h"

#include <cassert>

#include "HighsExternalApi.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsSolution.h"
#include "util/stringutil.h"

using std::min;

HighsStatus solveLpIpx(HighsLpSolverObject& solver_object) {
  return solveLpIpx(solver_object.options_, solver_object.timer_,
                    solver_object.lp_, solver_object.basis_,
                    solver_object.solution_, solver_object.model_status_,
                    solver_object.highs_info_, solver_object.callback_);
}

HighsStatus solveLpIpx(const HighsOptions& options, HighsTimer& timer,
                       const HighsLp& lp, HighsBasis& highs_basis,
                       HighsSolution& highs_solution,
                       HighsModelStatus& model_status, HighsInfo& highs_info,
                       HighsCallback& callback) {
  // Use IPX to try to solve the LP
  //
  // Can return HighsModelStatus (HighsStatus) values:
  //
  // 1. kSolveError (kError) if various unlikely solution errors occur
  //
  // 2. kTimeLimit (kWarning) if time limit is reached
  //
  // 3. kIterationLimit (kWarning) if iteration limit is reached
  //
  // 4. kUnknown (kWarning) if IPM makes no progress or if
  // IPM/crossover are imprecise
  //
  // 5. kInfeasible (kOk) if IPM identifies primal infeasibility
  //
  // 6. kUnboundedOrInfeasible (kOk) if IPM identifies dual
  // infeasibility
  //
  // kOptimal (kOk) if IPM/crossover identify optimality
  //
  // With a non-error return, if just IPM has been run then a
  // non-vertex primal solution is obtained; if crossover has been run
  // then a basis and primal+dual solution are obtained.
  //
  //
  // Indicate that there is no valid primal solution, dual solution or basis
  highs_basis.valid = false;
  highs_solution.value_valid = false;
  highs_solution.dual_valid = false;
  // Indicate that no imprecise solution has (yet) been found
  resetModelStatusAndHighsInfo(model_status, highs_info);
  // Create the LpSolver instance
  ipx::LpSolver lps;
  // Use the current HiGHS time as an offset for the lps.control_
  // elapsed time
  lps.setTimerOffset(timer.read());

  // Set IPX parameters
  //
  // Cannot set internal IPX parameters directly since they are
  // private, so create instance of parameters
  ipx::Parameters parameters;
  // Set IPX parameters from options
  //
  // Set display according to output
  parameters.display = 1;
  if (!options.output_flag | !options.log_to_console) parameters.display = 0;
  // Modify parameters.debug according to log_dev_level
  parameters.debug = 0;
  if (options.log_dev_level == kHighsLogDevLevelDetailed) {
    // Default options.log_dev_level setting is kHighsLogDevLevelNone, yielding
    // default setting debug = 0
    parameters.debug = 0;
  } else if (options.log_dev_level == kHighsLogDevLevelInfo) {
    parameters.debug = 2;
  } else if (options.log_dev_level == kHighsLogDevLevelVerbose) {
    parameters.debug = 4;
  }
  parameters.highs_logging = true;
  parameters.timeless_log = options.timeless_log;
  // Use this to change print_interval from default value (5s) to
  // other value for debugging
  //
  //  parameters.print_interval = 1e-4;
  parameters.log_options = &options.log_options;
  // Just test feasibility and optimality tolerances for now
  // ToDo Set more parameters
  //
  // Translate dualization option
  //
  // parameters.dualize = -2 => Possibly dualize - Filippo style
  // parameters.dualize = -1 => Possibly dualize - Lukas style
  // parameters.dualize = 0 => No dualization
  // parameters.dualize = 1 => Perform dualization
  if (options.ipx_dualize_strategy == kIpxDualizeStrategyOn) {
    parameters.dualize = 1;
  } else if (options.ipx_dualize_strategy == kIpxDualizeStrategyOff) {
    parameters.dualize = 0;
  } else if (options.ipx_dualize_strategy == kIpxDualizeStrategyLukas) {
    parameters.dualize = -1;
  } else if (options.ipx_dualize_strategy == kIpxDualizeStrategyFilippo) {
    parameters.dualize = -2;
  } else {
    assert(111 == 222);
  }

  parameters.ipm_feasibility_tol = min(options.primal_feasibility_tolerance,
                                       options.dual_feasibility_tolerance);
  parameters.ipm_optimality_tol = options.ipm_optimality_tolerance;
  parameters.start_crossover_tol = options.start_crossover_tolerance;

  if (options.kkt_tolerance != kDefaultKktTolerance) {
    parameters.ipm_feasibility_tol = options.kkt_tolerance;
    parameters.ipm_optimality_tol = 1e-1 * options.kkt_tolerance;
    parameters.start_crossover_tol = 1e-1 * options.kkt_tolerance;
  }

  parameters.analyse_basis_data =
      kHighsAnalysisLevelNlaData & options.highs_analysis_level;
  // Now that the lps.control_ elapsed time includes the HiGHS time,
  // can use the HiGHS time limit
  parameters.time_limit = options.time_limit;
  parameters.ipm_maxiter =
      options.ipm_iteration_limit - highs_info.ipm_iteration_count;
  // Determine if crossover is to be run or not
  //
  // When doing analytic centring calculations, crossover must not be
  // run
  if (options.run_centring) {
    parameters.run_crossover = 0;
  } else if (options.run_crossover == kHighsOnString) {
    parameters.run_crossover = 1;
  } else if (options.run_crossover == kHighsOffString) {
    parameters.run_crossover = 0;
  } else {
    assert(options.run_crossover == kHighsChooseString);
    parameters.run_crossover = -1;
  }
  if (!parameters.run_crossover) {
    // If crossover is sure not to be run, then set crossover_start_ to
    // -1 so that IPX can terminate according to its feasibility and
    // optimality tolerances
    parameters.start_crossover_tol = -1;
  }

  parameters.run_centring = options.run_centring ? 1 : 0;
  parameters.max_centring_steps = options.max_centring_steps;
  parameters.centring_ratio_tolerance = options.centring_ratio_tolerance;

  // Set the internal IPX parameters
  lps.SetParameters(parameters);

  // Set pointer to any callback
  lps.SetCallback(&callback);

  ipx::Int num_col, num_row;
  double offset;
  std::vector<ipx::Int> Ap, Ai;
  std::vector<double> objective, col_lb, col_ub, Av, rhs;
  std::vector<char> constraint_type;
  fillInIpxData(lp, num_col, num_row, offset, objective, col_lb, col_ub, Ap, Ai,
                Av, rhs, constraint_type);
  highsLogUser(options.log_options, HighsLogType::kInfo,
               "IPX model has %" HIGHSINT_FORMAT " rows, %" HIGHSINT_FORMAT
               " columns and %" HIGHSINT_FORMAT " nonzeros\n",
               num_row, num_col, Ap[num_col]);

  ipx::Int load_status = lps.LoadModel(
      num_col, offset, objective.data(), col_lb.data(), col_ub.data(), num_row,
      Ap.data(), Ai.data(), Av.data(), rhs.data(), constraint_type.data());

  if (load_status) {
    model_status = HighsModelStatus::kSolveError;
    return HighsStatus::kError;
  }

  // Use IPX to solve the LP!
  ipx::Int solve_status = lps.Solve();

  const bool report_solve_data =
      kHighsAnalysisLevelSolverSummaryData & options.highs_analysis_level;
  // Get solver and solution information.
  // Struct ipx_info defined in ipx/ipx_info.h
  const ipx::Info ipx_info = lps.GetInfo();
  if (report_solve_data) reportSolveData(options.log_options, ipx_info);
  highs_info.ipm_iteration_count += (HighsInt)ipx_info.iter;
  highs_info.crossover_iteration_count += (HighsInt)ipx_info.updates_crossover;
  // Record the outcome of any centring, since the model status only
  // indicates whether the LP is solved. Centring is not tried if the
  // IPM stops before reaching optimality
  highs_info.centring_status = !ipx_info.centring_tried ? kCentringStatusNone
                               : ipx_info.centring_success
                                   ? kCentringStatusSuccess
                                   : kCentringStatusFailed;

  // If not solved...
  if (solve_status != IPX_STATUS_solved) {
    const HighsStatus solve_return_status =
        reportIpxSolveStatus(options, solve_status, ipx_info.errflag);
    // Return error if IPX solve error has occurred
    if (solve_return_status == HighsStatus::kError) {
      model_status = HighsModelStatus::kSolveError;
      return HighsStatus::kError;
    }
  }
  bool ipm_status = true;
  const HighsStatus ipm_return_status =
      reportIpxIpmCrossoverStatus(options, ipx_info.status_ipm, ipm_status);
  ipm_status = false;
  const HighsStatus crossover_return_status = reportIpxIpmCrossoverStatus(
      options, ipx_info.status_crossover, ipm_status);
  // Return error if IPX IPM or crossover error has occurred
  if (ipm_return_status == HighsStatus::kError ||
      crossover_return_status == HighsStatus::kError) {
    model_status = HighsModelStatus::kSolveError;
    return HighsStatus::kError;
  }
  // Should only reach here if Solve() returned IPX_STATUS_solved or
  // IPX_STATUS_stopped
  if (ipxStatusError(
          solve_status != IPX_STATUS_solved &&
              solve_status != IPX_STATUS_stopped,
          options, "Ipx",
          "solve_status should be solved or stopped here but value is",
          (int)solve_status))
    return HighsStatus::kError;

  // Only error returns so far
  //

  if (solve_status == IPX_STATUS_stopped) {
    // IPX stopped, so there's certainly no basic solution. Get the
    // non-vertex solution, though. This needs the model status to
    // know whether to worry about dual infeasibilities.
    const HighsModelStatus local_model_status = HighsModelStatus::kUnknown;
    getHighsNonVertexSolution(options, lp, num_col, num_row, rhs,
                              constraint_type, lps, local_model_status,
                              highs_solution);
    //
    // Look at the reason why IPX stopped
    //
    // Return error if stopped status settings occur that JAJH doesn't
    // think should happen
    //
    //==============
    // For crossover
    //==============
    if (illegalIpxStoppedCrossoverStatus(ipx_info, options))
      return HighsStatus::kError;
    // Can stop and reach time limit
    if (ipx_info.status_crossover == IPX_STATUS_time_limit) {
      model_status = HighsModelStatus::kTimeLimit;
      return HighsStatus::kWarning;
    }
    //========
    // For IPM
    //========
    //
    // Note that IPX can stop with IPM optimal, imprecise,
    // primal_infeas or dual_infeas, due to crossover stopping with
    // time limit, and this is why crossover returns are tested first
    if (illegalIpxStoppedIpmStatus(ipx_info, options))
      return HighsStatus::kError;
    // Can stop with user interrupt
    // Can stop with time limit
    // Can stop with iter limit
    // Can stop with no progress
    if (ipx_info.status_ipm == IPX_STATUS_user_interrupt) {
      model_status = HighsModelStatus::kInterrupt;
      return HighsStatus::kWarning;
    } else if (ipx_info.status_ipm == IPX_STATUS_time_limit) {
      model_status = HighsModelStatus::kTimeLimit;
      return HighsStatus::kWarning;
    } else if (ipx_info.status_ipm == IPX_STATUS_iter_limit) {
      model_status = HighsModelStatus::kIterationLimit;
      return HighsStatus::kWarning;
    } else {
      assert(ipx_info.status_ipm == IPX_STATUS_no_progress);
      reportIpmNoProgress(options, ipx_info);
      model_status = HighsModelStatus::kUnknown;
      return HighsStatus::kWarning;
    }
  }
  // Should only reach here if Solve() returned IPX_STATUS_solved
  if (ipxStatusError(solve_status != IPX_STATUS_solved, options, "Ipx",
                     "solve_status should be solved here but value is",
                     (int)solve_status))
    return HighsStatus::kError;
  // Return error if solved status settings occur that JAJH doesn't
  // think should happen
  if (illegalIpxSolvedStatus(ipx_info, options)) return HighsStatus::kError;
  //==============
  // For crossover
  //==============
  // Can be not run
  // Can solve and be optimal
  // Can solve and be imprecise
  //========
  // For IPM
  //========
  // Can solve and be optimal
  // Can solve and be imprecise
  // Can solve and be primal_infeas
  // Can solve and be dual_infeas
  if (ipx_info.status_ipm == IPX_STATUS_primal_infeas ||
      ipx_info.status_ipm == IPX_STATUS_dual_infeas) {
    // IPM identified primal or dual infeasibility: crossover will not
    // have run, so get the non-vertex solution and return
    if (ipx_info.status_ipm == IPX_STATUS_primal_infeas) {
      model_status = HighsModelStatus::kInfeasible;
    } else if (ipx_info.status_ipm == IPX_STATUS_dual_infeas) {
      model_status = HighsModelStatus::kUnboundedOrInfeasible;
    }
    getHighsNonVertexSolution(options, lp, num_col, num_row, rhs,
                              constraint_type, lps, model_status,
                              highs_solution);
    return HighsStatus::kOk;
  }

  // Should only reach here if IPM is optimal or imprecise
  if (ipxStatusError(ipx_info.status_ipm != IPX_STATUS_optimal &&
                         ipx_info.status_ipm != IPX_STATUS_imprecise,
                     options, "Ipx",
                     "ipm status should be not run, optimal or imprecise "
                     "but value is",
                     (int)ipx_info.status_ipm))
    return HighsStatus::kError;

  // Should only reach here if crossover is not run, optimal or imprecise
  if (ipxStatusError(ipx_info.status_crossover != IPX_STATUS_not_run &&
                         ipx_info.status_crossover != IPX_STATUS_optimal &&
                         ipx_info.status_crossover != IPX_STATUS_imprecise,
                     options, "Ipx",
                     "crossover status should be not run, optimal or imprecise "
                     "but value is",
                     (int)ipx_info.status_crossover))
    return HighsStatus::kError;

  // Basic solution depends on crossover being run
  const bool have_basic_solution =
      ipx_info.status_crossover != IPX_STATUS_not_run;
  // Both crossover and IPM can be imprecise
  const bool imprecise_solution =
      ipx_info.status_crossover == IPX_STATUS_imprecise ||
      ipx_info.status_ipm == IPX_STATUS_imprecise;
  if (have_basic_solution) {
    IpxSolution ipx_solution;
    ipx_solution.num_col = num_col;
    ipx_solution.num_row = num_row;
    ipx_solution.ipx_col_value.resize(num_col);
    ipx_solution.ipx_row_value.resize(num_row);
    ipx_solution.ipx_col_dual.resize(num_col);
    ipx_solution.ipx_row_dual.resize(num_row);
    ipx_solution.ipx_row_status.resize(num_row);
    ipx_solution.ipx_col_status.resize(num_col);
    ipx::Int errflag = lps.GetBasicSolution(
        ipx_solution.ipx_col_value.data(), ipx_solution.ipx_row_value.data(),
        ipx_solution.ipx_row_dual.data(), ipx_solution.ipx_col_dual.data(),
        ipx_solution.ipx_row_status.data(), ipx_solution.ipx_col_status.data());
    if (errflag != 0) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "IPX crossover getting basic solution: flag = %d\n",
                   (int)errflag);
      return HighsStatus::kError;
    }
    // Convert the IPX basic solution to a HiGHS basic solution
    HighsStatus status = ipxBasicSolutionToHighsBasicSolution(
        options.log_options, lp, rhs, constraint_type, ipx_solution,
        highs_basis, highs_solution);
    if (status != HighsStatus::kOk) {
      highsLogUser(
          options.log_options, HighsLogType::kError,
          "Failed to convert IPX basic solution to Highs basic solution\n");
      return HighsStatus::kError;
    }

  } else {
    // No basic solution, so get a non-vertex HiGHS solution. This
    // needs the model status to know whether to worry about dual
    // infeasibilities.
    const HighsModelStatus local_model_status =
        imprecise_solution ? HighsModelStatus::kUnknown
                           : HighsModelStatus::kOptimal;
    getHighsNonVertexSolution(options, lp, num_col, num_row, rhs,
                              constraint_type, lps, local_model_status,
                              highs_solution);
    assert(!highs_basis.valid);
  }
  highs_info.basis_validity =
      highs_basis.valid ? kBasisValidityValid : kBasisValidityInvalid;
  HighsStatus return_status;
  if (imprecise_solution) {
    model_status = HighsModelStatus::kUnknown;
    return_status = HighsStatus::kWarning;
  } else {
    model_status = HighsModelStatus::kOptimal;
    return_status = HighsStatus::kOk;
  }
  return return_status;
}

HighsStatus solveLpHipo(HighsLpSolverObject& solver_object) {
  return solveHipo(solver_object.options_, solver_object.timer_,
                   solver_object.lp_, HighsHessian{}, solver_object.basis_,
                   solver_object.solution_, solver_object.model_status_,
                   solver_object.highs_info_, solver_object.callback_);
}

HighsStatus solveQpHipo(HighsQpSolverObject& solver_object) {
  return solveHipo(solver_object.options_, solver_object.timer_,
                   solver_object.model_.lp_, solver_object.model_.hessian_,
                   solver_object.basis_, solver_object.solution_,
                   solver_object.model_status_, solver_object.highs_info_,
                   solver_object.callback_);
}

static HighsInt prepareOpenBLAS(const HighsOptions& options) {
  // OpenBLAS threads are a process-wide setting, so the user can opt out of
  // changing them
  if (!options.hipo_serial_blas) return 0;

  // force openblas to run in serial, for determinism and better performance
  // no-op if openblas is not used
  const int threads_used = HighsExtras::blas::openblas_set_num_threads(1);
  if (is_substring_case(HighsExtras::blas::getInfo()->provider, "openblas") &&
      threads_used != 1) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "OpenBLAS failed to set the number of threads to 1\n");
    return 1;
  }
  return 0;
}

HighsStatus solveHipo(const HighsOptions& options, HighsTimer& timer,
                      const HighsLp& lp, const HighsHessian& Q,
                      HighsBasis& highs_basis, HighsSolution& highs_solution,
                      HighsModelStatus& model_status, HighsInfo& highs_info,
                      HighsCallback& callback) {
  // Use HiPO
  //
  // Can return HighsModelStatus (HighsStatus) values:
  //
  // 1. kSolveError (kError) if various unlikely solution errors occur
  //
  // 2. kTimeLimit (kWarning) if time limit is reached
  //
  // 3. kIterationLimit (kWarning) if iteration limit is reached
  //
  // 4. kUnknown (kWarning) if HiPO makes no progress or if
  // IPM/crossover are imprecise
  //
  // 5. kInfeasible (kOk) if HiPO identifies primal infeasibility
  //
  // 6. kUnboundedOrInfeasible (kOk) if HiPO identifies dual
  // infeasibility
  //
  // kOptimal (kOk) if HiPO/crossover identify optimality
  //
  // With a non-error return, if just HiPO has been run then a
  // non-vertex primal solution is obtained; if crossover has been run
  // then a basis and primal+dual solution are obtained.
  //
  //
  // Indicate that there is no valid primal solution, dual solution or basis
  highs_basis.valid = false;
  highs_solution.value_valid = false;
  highs_solution.dual_valid = false;
  // Indicate that no imprecise solution has (yet) been found
  resetModelStatusAndHighsInfo(model_status, highs_info);

  if (prepareOpenBLAS(options)) {
    model_status = HighsModelStatus::kSolveError;
    return HighsStatus::kError;
  }

  // Create solver instance
  hipo::Solver hipo{};
  // This creates ipx::LpSolver ipx_lps_, in case HiPO has to switch
  // to IPX, so use the current HiGHS time as an offset for the
  // ipx_lps.control_ elapsed time
  hipo.setIpxTimerOffset(timer.read());

  if (options.kkt_tolerance != kDefaultKktTolerance) {
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "IpxWrapper: feasibility_tol = %g; optimality_tol = %g; "
                 "crossover_tol = %g\n",
                 options.kkt_tolerance, 1e-1 * options.kkt_tolerance,
                 1e-1 * options.kkt_tolerance);
  }

  hipo.setOptions(options);
  hipo.setTimer(timer);
  hipo.setCallback(callback);

  // Load the problem
  hipo::Int load_status = hipo.load(lp, Q);
  if (load_status) {
    model_status = HighsModelStatus::kSolveError;
    return HighsStatus::kError;
  }

  // This information about the problem loaded into HiPO is needed for later
  HighsInt num_row, num_col;
  hipo.getOriginalDims(num_row, num_col);
  std::vector<double> rhs;
  std::vector<char> constraints;
  fillInRhsAndConstraints(lp, rhs, constraints);

  hipo.solve();

  // const bool report_solve_data =
  //    kHighsAnalysisLevelSolverSummaryData & options.highs_analysis_level;

  // hipo.solved(), hipo.stopped(), hipo.failed() can be used to query if the
  // status belongs to the solved, stopped or failed group.

  // Get solver and solution information.
  const hipo::Info hipo_info = hipo.getInfo();
  hipo::Status solve_status = hipo_info.status;
  highs_info.ipm_iteration_count +=
      hipo_info.ipm_iter + hipo_info.ipx_info.iter;
  highs_info.crossover_iteration_count += hipo_info.ipx_info.updates_crossover;

  // Report hipo status
  const HighsStatus solve_return_status =
      reportHipoStatus(options, solve_status, hipo);
  if (solve_return_status == HighsStatus::kError) {
    model_status = HighsModelStatus::kSolveError;
    return HighsStatus::kError;
  }

  // Report crossover status
  const HighsStatus crossover_return_status = reportHipoCrossoverStatus(
      options, hipo_info.ipx_info.status_crossover, !Q.empty());
  if (crossover_return_status == HighsStatus::kError) {
    model_status = HighsModelStatus::kSolveError;
    return HighsStatus::kError;
  }

  // Failures should have been handled. Status should be stopped or solved.
  if (ipxStatusError(!hipo.solved() && !hipo.stopped(), options, "Hipo",
                     "status should be solved or stopped but value is",
                     solve_status))
    return HighsStatus::kError;

  if (hipo.stopped()) {
    const HighsModelStatus local_model_status = HighsModelStatus::kUnknown;

    getHipoNonVertexSolution(options, lp, num_col, num_row, rhs, constraints,
                             hipo, local_model_status, highs_solution);

    // For crossover
    if (illegalIpxStoppedCrossoverStatus(hipo_info.ipx_info, options))
      return HighsStatus::kError;
    // Can stop and reach time limit
    if (hipo_info.ipx_info.status_crossover == IPX_STATUS_time_limit) {
      model_status = HighsModelStatus::kTimeLimit;
      return HighsStatus::kWarning;
    }

    // if crossover didn't time out, then solver can only stop as follows
    if (solve_status == hipo::kStatusUserInterrupt) {
      model_status = HighsModelStatus::kInterrupt;
      return HighsStatus::kWarning;
    } else if (solve_status == hipo::kStatusTimeLimit) {
      model_status = HighsModelStatus::kTimeLimit;
      return HighsStatus::kWarning;
    } else if (solve_status == hipo::kStatusMaxIter) {
      model_status = HighsModelStatus::kIterationLimit;
      return HighsStatus::kWarning;
    } else if (solve_status == hipo::kStatusNoProgress) {
      reportHipoNoProgress(options, hipo_info);
      model_status = HighsModelStatus::kUnknown;
      return HighsStatus::kWarning;
    } else {
      assert(1 == 0);
    }
  }

  // Stopped status should have been handled. Status should be solved.
  if (ipxStatusError(!hipo.solved(), options, "Hipo",
                     "status should be solved but value is", solve_status))
    return HighsStatus::kError;

  // primal/dual infeasible
  if (solve_status == hipo::kStatusPrimalInfeasible ||
      solve_status == hipo::kStatusDualInfeasible) {
    if (solve_status == hipo::kStatusPrimalInfeasible)
      model_status = HighsModelStatus::kInfeasible;
    else
      model_status = HighsModelStatus::kUnboundedOrInfeasible;

    getHipoNonVertexSolution(options, lp, num_col, num_row, rhs, constraints,
                             hipo, model_status, highs_solution);

    return HighsStatus::kOk;
  }

  // Status should be optimal or imprecise
  if (ipxStatusError(solve_status != hipo::kStatusOptimal &&
                         solve_status != hipo::kStatusImprecise,
                     options, "Hipo",
                     "status should be optimal or imprecise but value is",
                     solve_status))
    return HighsStatus::kError;

  const bool have_basic_solution =
      hipo_info.ipx_info.status_crossover != IPX_STATUS_not_run;

  const bool imprecise_solution =
      hipo_info.status == hipo::kStatusImprecise ||
      hipo_info.ipx_info.status_crossover == IPX_STATUS_imprecise;

  if (have_basic_solution) {
    IpxSolution hipo_solution;
    hipo_solution.num_col = num_col;
    hipo_solution.num_row = num_row;
    hipo_solution.ipx_col_value.resize(num_col);
    hipo_solution.ipx_row_value.resize(num_row);
    hipo_solution.ipx_col_dual.resize(num_col);
    hipo_solution.ipx_row_dual.resize(num_row);
    hipo_solution.ipx_row_status.resize(num_row);
    hipo_solution.ipx_col_status.resize(num_col);

    hipo::Int errflag = hipo.getBasicSolution(
        hipo_solution.ipx_col_value, hipo_solution.ipx_row_value,
        hipo_solution.ipx_row_dual, hipo_solution.ipx_col_dual,
        hipo_solution.ipx_row_status.data(),
        hipo_solution.ipx_col_status.data());

    if (errflag) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "IPX crossover getting basic solution: flag = %d\n",
                   (int)errflag);
      return HighsStatus::kError;
    }
    // Convert the IPX basic solution to a HiGHS basic solution
    HighsStatus status = ipxBasicSolutionToHighsBasicSolution(
        options.log_options, lp, rhs, constraints, hipo_solution, highs_basis,
        highs_solution);
    if (status != HighsStatus::kOk) {
      highsLogUser(
          options.log_options, HighsLogType::kError,
          "Failed to convert IPX basic solution to Highs basic solution\n");
      return HighsStatus::kError;
    }
  } else {
    const HighsModelStatus local_model_status =
        imprecise_solution ? HighsModelStatus::kUnknown
                           : HighsModelStatus::kOptimal;
    getHipoNonVertexSolution(options, lp, num_col, num_row, rhs, constraints,
                             hipo, local_model_status, highs_solution);
    assert(!highs_basis.valid);
  }

  highs_info.basis_validity =
      highs_basis.valid ? kBasisValidityValid : kBasisValidityInvalid;

  HighsStatus return_status;
  if (imprecise_solution) {
    model_status = HighsModelStatus::kUnknown;
    return_status = HighsStatus::kWarning;
  } else {
    model_status = HighsModelStatus::kOptimal;
    return_status = HighsStatus::kOk;
  }
  return return_status;
}

void fillInIpxData(const HighsLp& lp, ipx::Int& num_col, ipx::Int& num_row,
                   double& offset, std::vector<double>& obj,
                   std::vector<double>& col_lb, std::vector<double>& col_ub,
                   std::vector<ipx::Int>& Ap, std::vector<ipx::Int>& Ai,
                   std::vector<double>& Ax, std::vector<double>& rhs,
                   std::vector<char>& constraint_type) {
  num_col = lp.num_col_;
  num_row = lp.num_row_;

  // For each row with both a lower and an upper bound introduce one new column
  // so num_col may increase. Ignore each free row so num_row may decrease.
  // lba <= a'x <= uba becomes
  // a'x-s = 0 and lba <= s <= uba.

  // For each row with bounds on both sides introduce explicit slack and
  // transfer bounds.
  assert((HighsInt)lp.row_lower_.size() == num_row);
  assert((HighsInt)lp.row_upper_.size() == num_row);

  std::vector<HighsInt> general_bounded_rows;
  std::vector<HighsInt> free_rows;

  for (HighsInt row = 0; row < num_row; row++)
    if (lp.row_lower_[row] < lp.row_upper_[row] &&
        lp.row_lower_[row] > -kHighsInf && lp.row_upper_[row] < kHighsInf)
      general_bounded_rows.push_back(row);
    else if (lp.row_lower_[row] <= -kHighsInf &&
             lp.row_upper_[row] >= kHighsInf)
      free_rows.push_back(row);

  const HighsInt num_slack = general_bounded_rows.size();

  fillInRhsAndConstraints(lp, rhs, constraint_type);

  std::vector<HighsInt> reduced_rowmap(lp.num_row_, -1);
  if (free_rows.size() > 0) {
    HighsInt counter = 0;
    HighsInt findex = 0;
    const HighsInt free_rows_size = static_cast<HighsInt>(free_rows.size());
    for (HighsInt row = 0; row < lp.num_row_; row++) {
      if (findex < free_rows_size && free_rows[findex] == row) {
        findex++;
        continue;
      } else {
        reduced_rowmap[row] = counter;
        counter++;
      }
    }
  } else {
    for (HighsInt k = 0; k < lp.num_row_; k++) reduced_rowmap[k] = k;
  }
  num_row -= free_rows.size();
  num_col += num_slack;

  std::vector<HighsInt> sizes(num_col, 0);

  for (HighsInt col = 0; col < lp.num_col_; col++)
    for (HighsInt k = lp.a_matrix_.start_[col];
         k < lp.a_matrix_.start_[col + 1]; k++) {
      HighsInt row = lp.a_matrix_.index_[k];
      if (lp.row_lower_[row] > -kHighsInf || lp.row_upper_[row] < kHighsInf)
        sizes[col]++;
    }
  // Copy Astart and Aindex to ipx::Int array.
  HighsInt nnz = lp.a_matrix_.index_.size();
  Ap.resize(num_col + 1);
  Ai.reserve(nnz + num_slack);
  Ax.reserve(nnz + num_slack);

  // Set starting points of original and newly introduced columns.
  Ap[0] = 0;
  for (HighsInt col = 0; col < lp.num_col_; col++) {
    Ap[col + 1] = Ap[col] + sizes[col];
  }
  for (HighsInt col = lp.num_col_; col < (HighsInt)num_col; col++) {
    Ap[col + 1] = Ap[col] + 1;
  }
  for (HighsInt k = 0; k < nnz; k++) {
    HighsInt row = lp.a_matrix_.index_[k];
    if (lp.row_lower_[row] > -kHighsInf || lp.row_upper_[row] < kHighsInf) {
      Ai.push_back(reduced_rowmap[lp.a_matrix_.index_[k]]);
      Ax.push_back(lp.a_matrix_.value_[k]);
    }
  }

  for (HighsInt k = 0; k < num_slack; k++) {
    Ai.push_back((ipx::Int)reduced_rowmap[general_bounded_rows[k]]);
    Ax.push_back(-1);
  }

  // Column bound vectors.
  col_lb.resize(num_col);
  col_ub.resize(num_col);
  for (HighsInt col = 0; col < lp.num_col_; col++) {
    if (lp.col_lower_[col] <= -kHighsInf)
      col_lb[col] = -INFINITY;
    else
      col_lb[col] = lp.col_lower_[col];

    if (lp.col_upper_[col] >= kHighsInf)
      col_ub[col] = INFINITY;
    else
      col_ub[col] = lp.col_upper_[col];
  }
  for (HighsInt slack = 0; slack < num_slack; slack++) {
    const int row = general_bounded_rows[slack];
    col_lb[lp.num_col_ + slack] = lp.row_lower_[row];
    col_ub[lp.num_col_ + slack] = lp.row_upper_[row];
  }

  offset = HighsInt(lp.sense_) * lp.offset_;
  obj.resize(num_col);
  for (HighsInt col = 0; col < lp.num_col_; col++) {
    obj[col] = (HighsInt)lp.sense_ * lp.col_cost_[col];
  }
}

void fillInRhsAndConstraints(const HighsLp& lp, std::vector<double>& rhs,
                             std::vector<char>& constraint_type) {
  // For each row except free rows add entry to char array and set up rhs
  // vector

  const HighsInt num_row = lp.num_row_;

  rhs.reserve(num_row);
  constraint_type.reserve(num_row);

  for (int row = 0; row < num_row; row++) {
    if (lp.row_lower_[row] > -kHighsInf && lp.row_upper_[row] >= kHighsInf) {
      rhs.push_back(lp.row_lower_[row]);
      constraint_type.push_back('>');
    } else if (lp.row_lower_[row] <= -kHighsInf &&
               lp.row_upper_[row] < kHighsInf) {
      rhs.push_back(lp.row_upper_[row]);
      constraint_type.push_back('<');
    } else if (lp.row_lower_[row] == lp.row_upper_[row]) {
      rhs.push_back(lp.row_upper_[row]);
      constraint_type.push_back('=');
    } else if (lp.row_lower_[row] > -kHighsInf &&
               lp.row_upper_[row] < kHighsInf) {
      // general bounded
      rhs.push_back(0);
      constraint_type.push_back('=');
    }
  }
}

HighsStatus reportIpxSolveStatus(const HighsOptions& options,
                                 const ipx::Int solve_status,
                                 const ipx::Int error_flag) {
  if (solve_status == IPX_STATUS_solved) {
    highsLogUser(options.log_options, HighsLogType::kInfo, "Ipx: Solved\n");
    return HighsStatus::kOk;
  } else if (solve_status == IPX_STATUS_stopped) {
    highsLogUser(options.log_options, HighsLogType::kWarning, "Ipx: Stopped\n");
    return HighsStatus::kWarning;
    // Remaining cases are errors so drop through to return HighsStatus::kError;
  } else if (solve_status == IPX_STATUS_no_model) {
    if (error_flag == IPX_ERROR_argument_null) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Ipx: Invalid input - argument_null\n");
    } else if (error_flag == IPX_ERROR_invalid_dimension) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Ipx: Invalid input - invalid dimension\n");
    } else if (error_flag == IPX_ERROR_invalid_matrix) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Ipx: Invalid input - invalid matrix\n");
    } else if (error_flag == IPX_ERROR_invalid_vector) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Ipx: Invalid input - invalid vector\n");
    } else if (error_flag == IPX_ERROR_invalid_basis) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Ipx: Invalid input - invalid basis\n");
    } else {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Ipx: Invalid input - unrecognised error\n");
    }
  } else if (solve_status == IPX_STATUS_out_of_memory) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Ipx: Out of memory\n");
  } else if (solve_status == IPX_STATUS_internal_error) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Ipx: Internal error %" HIGHSINT_FORMAT "\n", (int)error_flag);
  } else {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Ipx: unrecognised solve status = %" HIGHSINT_FORMAT "\n",
                 (int)solve_status);
  }
  return HighsStatus::kError;
}

HighsStatus reportIpxIpmCrossoverStatus(const HighsOptions& options,
                                        const ipx::Int status,
                                        const bool ipm_status) {
  std::string method_name;
  if (ipm_status)
    method_name = "IPM      ";
  else
    method_name = "Crossover";
  if (status == IPX_STATUS_not_run) {
    if (ipm_status || options.run_crossover == kHighsOnString) {
      // Warn if method not run is IPM or method not run is crossover
      // and run_crossover option is "on"
      highsLogUser(options.log_options, HighsLogType::kWarning,
                   "Ipx: %s not run\n", method_name.c_str());
      return HighsStatus::kWarning;
    }
    // OK if method not run is crossover and run_crossover option is
    // not "on"
    return HighsStatus::kOk;
  } else if (status == IPX_STATUS_optimal) {
    highsLogUser(options.log_options, HighsLogType::kInfo, "Ipx: %s optimal\n",
                 method_name.c_str());
    return HighsStatus::kOk;
  } else if (status == IPX_STATUS_imprecise) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s imprecise\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_primal_infeas) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s primal infeasible\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_dual_infeas) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s dual infeasible\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_user_interrupt) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s user interrupt\n", method_name.c_str());
    return HighsStatus::kOk;
  } else if (status == IPX_STATUS_time_limit) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s reached time limit\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_iter_limit) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s reached iteration limit\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_no_progress) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ipx: %s no progress\n", method_name.c_str());
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_failed) {
    highsLogUser(options.log_options, HighsLogType::kError, "Ipx: %s failed\n",
                 method_name.c_str());
    return HighsStatus::kError;
  } else if (status == IPX_STATUS_debug) {
    highsLogUser(options.log_options, HighsLogType::kError, "Ipx: %s debug\n",
                 method_name.c_str());
    return HighsStatus::kError;
  } else {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Ipx: %s unrecognised status\n", method_name.c_str());
    return HighsStatus::kError;
  }
  return HighsStatus::kError;
}

bool ipxStatusError(const bool status_error, const HighsOptions& options,
                    const std::string& solver, const std::string& message,
                    const int value) {
  if (status_error) {
    if (value < 0) {
      highsLogUser(options.log_options, HighsLogType::kError, "%s: %s\n",
                   solver.c_str(), message.c_str());
    } else {
      highsLogUser(options.log_options, HighsLogType::kError, "%s: %s %d\n",
                   solver.c_str(), message.c_str(), value);
    }
    fflush(NULL);
  }
  assert(!status_error);
  return status_error;
}

bool illegalIpxSolvedStatus(const ipx::Info& ipx_info,
                            const HighsOptions& options) {
  bool found_illegal_status = false;
  //========
  // For IPX
  //========
  // Can solve and be optimal
  // Can solve and be imprecise
  // Can solve and be primal_infeas
  // Can solve and be dual_infeas
  // Cannot solve and reach time limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_time_limit, options,
                     "Ipx",
                     "solved  status_ipm should not be IPX_STATUS_time_limit");
  // Cannot solve and reach iteration limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_iter_limit, options,
                     "Ipx",
                     "solved  status_ipm should not be IPX_STATUS_iter_limit");
  // Cannot solve and make no progress
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_no_progress, options,
                     "Ipx",
                     "solved  status_ipm should not be IPX_STATUS_no_progress");
  // Cannot solve and failed
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_failed, options, "Ipx",
                     "solved  status_ipm should not be IPX_STATUS_failed");
  // Cannot solve and debug
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_debug, options, "Ipx",
                     "solved  status_ipm should not be IPX_STATUS_debug");
  //==============
  // For crossover
  //==============
  // Can solve and be optimal
  // Can solve and be imprecise
  // Cannot solve with primal infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_primal_infeas, options, "Ipx",
          "solved  status_crossover should not be IPX_STATUS_primal_infeas");
  // Cannot solve with dual infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_dual_infeas, options, "Ipx",
          "solved  status_crossover should not be IPX_STATUS_dual_infeas");
  // Cannot solve and reach time limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_time_limit, options, "Ipx",
          "solved  status_crossover should not be IPX_STATUS_time_limit");
  // Cannot solve and reach time limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_iter_limit, options, "Ipx",
          "solved  status_crossover should not be IPX_STATUS_iter_limit");
  // Cannot solve and make no progress
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_no_progress, options, "Ipx",
          "solved  status_crossover should not be IPX_STATUS_no_progress");
  // Cannot solve and failed
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_failed, options, "Ipx",
          "solved  status_crossover should not be IPX_STATUS_failed");
  // Cannot solve and debug
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_crossover == IPX_STATUS_debug, options,
                     "Ipx",
                     "solved  status_crossover should not be IPX_STATUS_debug");
  return found_illegal_status;
}

bool illegalIpxStoppedIpmStatus(const ipx::Info& ipx_info,
                                const HighsOptions& options) {
  bool found_illegal_status = false;
  // Cannot stop and be optimal
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_optimal, options, "Ipx",
                     "stopped status_ipm should not be IPX_STATUS_optimal");
  // Cannot stop and be imprecise
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_imprecise, options,
                     "Ipx",
                     "stopped status_ipm should not be IPX_STATUS_imprecise");
  // Cannot stop with primal infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_ipm == IPX_STATUS_primal_infeas, options, "Ipx",
          "stopped status_ipm should not be IPX_STATUS_primal_infeas");
  // Cannot stop with dual infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_dual_infeas, options,
                     "Ipx",
                     "stopped status_ipm should not be IPX_STATUS_dual_infeas");
  // Can stop with time limit
  // Can stop with iter limit
  // Can stop with no progress
  // Cannot stop and failed - should be error return earlier
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_failed, options, "Ipx",
                     "stopped status_ipm should not be IPX_STATUS_failed");
  // Cannot stop and debug - should be error return earlier
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_debug, options, "Ipx",
                     "stopped status_ipm should not be IPX_STATUS_debug");
  return found_illegal_status;
}

bool illegalIpxStoppedCrossoverStatus(const ipx::Info& ipx_info,
                                      const HighsOptions& options) {
  bool found_illegal_status = false;
  // Cannot stop and be optimal
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_optimal, options, "Ipx",
          "stopped status_crossover should not be IPX_STATUS_optimal");
  // Cannot stop and be imprecise
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_imprecise, options, "Ipx",
          "stopped status_crossover should not be IPX_STATUS_imprecise");
  // Cannot stop with primal infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_primal_infeas, options, "Ipx",
          "stopped status_crossover should not be IPX_STATUS_primal_infeas");
  // Cannot stop with dual infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_dual_infeas, options, "Ipx",
          "stopped status_crossover should not be IPX_STATUS_dual_infeas");
  // Cannot stop and reach iteration limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_iter_limit, options, "Ipx",
          "stopped status_crossover should not be IPX_STATUS_iter_limit");
  // Can stop and reach time limit
  // Cannot stop with no_progress
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_no_progress, options, "Ipx",
          "stopped status_crossover should not be IPX_STATUS_no_progress");
  // Cannot stop and failed - should be error return earlier
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_failed, options, "Ipx",
          "stopped status_crossover should not be IPX_STATUS_failed");
  // Cannot stop and debug - should be error return earlier
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_crossover == IPX_STATUS_debug, options,
                     "Ipx",
                     "stopped status_crossover should not be IPX_STATUS_debug");
  return found_illegal_status;
}

void reportIpmNoProgress(const HighsOptions& options,
                         const ipx::Info& ipx_info) {
  highsLogUser(options.log_options, HighsLogType::kWarning,
               "No progress: primal objective value       = %11.4g\n",
               ipx_info.pobjval);
  highsLogUser(options.log_options, HighsLogType::kWarning,
               "No progress: max absolute primal residual = %11.4g\n",
               ipx_info.abs_presidual);
  highsLogUser(options.log_options, HighsLogType::kWarning,
               "No progress: max absolute   dual residual = %11.4g\n",
               ipx_info.abs_dresidual);
}

void reportHipoNoProgress(const HighsOptions& options,
                          const hipo::Info& hipo_info) {
  highsLogUser(options.log_options, HighsLogType::kWarning,
               "No progress: primal objective value       = %11.4g\n",
               hipo_info.p_obj);
  highsLogUser(options.log_options, HighsLogType::kWarning,
               "No progress: max absolute primal residual = %11.4g\n",
               hipo_info.p_res_abs);
  highsLogUser(options.log_options, HighsLogType::kWarning,
               "No progress: max absolute   dual residual = %11.4g\n",
               hipo_info.d_res_abs);
}

void getHighsNonVertexSolution(const HighsOptions& options, const HighsLp& lp,
                               const ipx::Int num_col, const ipx::Int num_row,
                               const std::vector<double>& rhs,
                               const std::vector<char>& constraint_type,
                               const ipx::LpSolver& lps,
                               const HighsModelStatus model_status,
                               HighsSolution& highs_solution) {
  // Get the interior solution (available if IPM was started).
  // GetInteriorSolution() returns the final IPM iterate, regardless if the
  // IPM terminated successfully or not. (Only in case of out-of-memory no
  // solution exists.)
  std::vector<double> x(num_col);
  std::vector<double> xl(num_col);
  std::vector<double> xu(num_col);
  std::vector<double> zl(num_col);
  std::vector<double> zu(num_col);
  std::vector<double> slack(num_row);
  std::vector<double> y(num_row);

  lps.GetInteriorSolution(x.data(), xl.data(), xu.data(), slack.data(),
                          y.data(), zl.data(), zu.data());

  ipxSolutionToHighsSolution(options, lp, rhs, constraint_type, num_col,
                             num_row, x, slack, y, zl, zu, highs_solution);
}

void getHipoNonVertexSolution(const HighsOptions& options, const HighsLp& lp,
                              const hipo::Int num_col, const hipo::Int num_row,
                              const std::vector<double>& rhs,
                              const std::vector<char>& constraint_type,
                              const hipo::Solver& hipo,
                              const HighsModelStatus model_status,
                              HighsSolution& highs_solution) {
  std::vector<double> x(num_col);
  std::vector<double> xl(num_col);
  std::vector<double> xu(num_col);
  std::vector<double> zl(num_col);
  std::vector<double> zu(num_col);
  std::vector<double> slack(num_row);
  std::vector<double> y(num_row);

  hipo.getInteriorSolution(x, xl, xu, slack, y, zl, zu);
  ipxSolutionToHighsSolution(options, lp, rhs, constraint_type, num_col,
                             num_row, x, slack, y, zl, zu, highs_solution);
}

void reportSolveData(const HighsLogOptions& log_options,
                     const ipx::Info& ipx_info) {
  highsLogDev(log_options, HighsLogType::kInfo, "\nIPX Solve data\n");
  highsLogDev(log_options, HighsLogType::kInfo, "    IPX       status = %4d\n",
              (int)ipx_info.status);
  highsLogDev(log_options, HighsLogType::kInfo, "    IPM       status = %4d\n",
              (int)ipx_info.status_ipm);
  highsLogDev(log_options, HighsLogType::kInfo, "    Crossover status = %4d\n",
              (int)ipx_info.status_crossover);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    IPX errflag      = %4d\n\n", (int)ipx_info.errflag);

  highsLogDev(log_options, HighsLogType::kInfo, "    LP variables   = %8d\n",
              (int)ipx_info.num_var);
  highsLogDev(log_options, HighsLogType::kInfo, "    LP constraints = %8d\n",
              (int)ipx_info.num_constr);
  highsLogDev(log_options, HighsLogType::kInfo, "    LP entries     = %8d\n\n",
              (int)ipx_info.num_entries);

  highsLogDev(log_options, HighsLogType::kInfo, "    Solver columns = %8d\n",
              (int)ipx_info.num_cols_solver);
  highsLogDev(log_options, HighsLogType::kInfo, "    Solver rows    = %8d\n",
              (int)ipx_info.num_rows_solver);
  highsLogDev(log_options, HighsLogType::kInfo, "    Solver entries = %8d\n\n",
              (int)ipx_info.num_entries_solver);

  highsLogDev(log_options, HighsLogType::kInfo, "    Dualized = %d\n",
              (int)ipx_info.dualized);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Number of dense columns detected = %d\n\n",
              (int)ipx_info.dense_cols);

  highsLogDev(log_options, HighsLogType::kInfo, "    Dependent rows    = %d\n",
              (int)ipx_info.dependent_rows);
  highsLogDev(log_options, HighsLogType::kInfo, "    Dependent cols    = %d\n",
              (int)ipx_info.dependent_cols);
  highsLogDev(log_options, HighsLogType::kInfo, "    Inconsistent rows = %d\n",
              (int)ipx_info.rows_inconsistent);
  highsLogDev(log_options, HighsLogType::kInfo, "    Inconsistent cols = %d\n",
              (int)ipx_info.cols_inconsistent);
  highsLogDev(log_options, HighsLogType::kInfo, "    Primal dropped    = %d\n",
              (int)ipx_info.primal_dropped);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Dual   dropped    = %d\n\n", (int)ipx_info.dual_dropped);

  highsLogDev(log_options, HighsLogType::kInfo,
              "    |Absolute primal residual| = %11.4g\n",
              ipx_info.abs_presidual);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    |Absolute   dual residual| = %11.4g\n",
              ipx_info.abs_dresidual);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    |Relative primal residual| = %11.4g\n",
              ipx_info.rel_presidual);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    |Relative   dual residual| = %11.4g\n\n",
              ipx_info.rel_dresidual);

  highsLogDev(log_options, HighsLogType::kInfo,
              "    Primal objective value     = %11.4g\n", ipx_info.pobjval);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Dual   objective value     = %11.4g\n", ipx_info.dobjval);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Relative objective gap     = %11.4g\n", ipx_info.rel_objgap);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Complementarity            = %11.4g\n\n",
              ipx_info.complementarity);

  highsLogDev(log_options, HighsLogType::kInfo, "    |x| = %11.4g\n",
              ipx_info.normx);
  highsLogDev(log_options, HighsLogType::kInfo, "    |y| = %11.4g\n",
              ipx_info.normy);
  highsLogDev(log_options, HighsLogType::kInfo, "    |z| = %11.4g\n\n",
              ipx_info.normz);

  highsLogDev(log_options, HighsLogType::kInfo,
              "    Objective value       = %11.4g\n", ipx_info.objval);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Primal infeasibility = %11.4g\n", ipx_info.primal_infeas);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Dual infeasibility   = %11.4g\n\n", ipx_info.dual_infeas);

  highsLogDev(log_options, HighsLogType::kInfo, "    IPM iter   = %d\n",
              (int)ipx_info.iter);
  highsLogDev(log_options, HighsLogType::kInfo, "    KKT iter 1 = %d\n",
              (int)ipx_info.kktiter1);
  highsLogDev(log_options, HighsLogType::kInfo, "    KKT iter 2 = %d\n",
              (int)ipx_info.kktiter2);
  highsLogDev(log_options, HighsLogType::kInfo, "    Basis repairs = %d\n",
              (int)ipx_info.basis_repairs);
  highsLogDev(log_options, HighsLogType::kInfo, "    Updates start     = %d\n",
              (int)ipx_info.updates_start);
  highsLogDev(log_options, HighsLogType::kInfo, "    Updates ipm       = %d\n",
              (int)ipx_info.updates_ipm);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Updates crossover = %d\n\n",
              (int)ipx_info.updates_crossover);

  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time total          = %8.2f\n\n", ipx_info.time_total);
  double sum_time = 0;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time IPM 1          = %8.2f\n", ipx_info.time_ipm1);
  sum_time += ipx_info.time_ipm1;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time IPM 2          = %8.2f\n", ipx_info.time_ipm2);
  sum_time += ipx_info.time_ipm2;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time starting basis = %8.2f\n",
              ipx_info.time_starting_basis);
  sum_time += ipx_info.time_starting_basis;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time crossover      = %8.2f\n", ipx_info.time_crossover);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Sum                 = %8.2f\n\n", sum_time);

  sum_time = 0;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time kkt_factorize  = %8.2f\n", ipx_info.time_kkt_factorize);
  sum_time += ipx_info.time_kkt_factorize;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time kkt_solve      = %8.2f\n", ipx_info.time_kkt_solve);
  sum_time += ipx_info.time_kkt_solve;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Sum                 = %8.2f\n\n", sum_time);

  sum_time = 0;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time maxvol         = %8.2f\n", ipx_info.time_maxvol);
  sum_time += ipx_info.time_maxvol;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time cr1            = %8.2f\n", ipx_info.time_cr1);
  sum_time += ipx_info.time_cr1;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time cr2            = %8.2f\n", ipx_info.time_cr2);
  sum_time += ipx_info.time_cr2;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Sum                 = %8.2f\n\n", sum_time);

  sum_time = 0;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time cr1_AAt        = %8.2f\n", ipx_info.time_cr1_AAt);
  sum_time += ipx_info.time_cr1_AAt;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time cr1_pre        = %8.2f\n", ipx_info.time_cr1_pre);
  sum_time += ipx_info.time_cr1_pre;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Sum  cr1            = %8.2f\n\n", sum_time);

  sum_time = 0;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time cr2_NNt        = %8.2f\n", ipx_info.time_cr2_NNt);
  sum_time += ipx_info.time_cr2_NNt;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time cr2_B          = %8.2f\n", ipx_info.time_cr2_B);
  sum_time += ipx_info.time_cr2_B;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time cr2_Bt         = %8.2f\n", ipx_info.time_cr2_Bt);
  sum_time += ipx_info.time_cr2_Bt;
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Sum  cr2            = %8.2f\n\n", sum_time);

  highsLogDev(log_options, HighsLogType::kInfo,
              "    Proportion of sparse FTRAN = %11.4g\n",
              ipx_info.ftran_sparse);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Proportion of sparse BTRAN = %11.4g\n\n",
              ipx_info.btran_sparse);

  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time FTRAN       = %8.2f\n", ipx_info.time_ftran);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time BTRAN       = %8.2f\n", ipx_info.time_btran);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time LU INVERT   = %8.2f\n", ipx_info.time_lu_invert);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time LU UPDATE   = %8.2f\n", ipx_info.time_lu_update);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Mean fill-in     = %11.4g\n", ipx_info.mean_fill);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Max fill-in      = %11.4g\n", ipx_info.max_fill);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Time symb INVERT = %11.4g\n\n", ipx_info.time_symb_invert);

  highsLogDev(log_options, HighsLogType::kInfo,
              "    Maxvol updates       = %d\n", (int)ipx_info.maxvol_updates);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Maxvol skipped       = %d\n", (int)ipx_info.maxvol_skipped);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Maxvol passes        = %d\n", (int)ipx_info.maxvol_passes);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Tableau num nonzeros = %d\n", (int)ipx_info.tbl_nnz);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Tbl max?             = %11.4g\n", ipx_info.tbl_max);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Frobnorm squared     = %11.4g\n", ipx_info.frobnorm_squared);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Lambda max           = %11.4g\n", ipx_info.lambdamax);
  highsLogDev(log_options, HighsLogType::kInfo,
              "    Volume increase      = %11.4g\n\n",
              ipx_info.volume_increase);
}

HighsStatus reportHipoStatus(const HighsOptions& options,
                             const hipo::Int status, const hipo::Solver& hipo) {
  if (hipo.solved()) {
    highsLogUser(options.log_options, HighsLogType::kInfo, "Hipo: Solved\n");
    return HighsStatus::kOk;
  }

  else if (hipo.stopped()) {
    highsLogUser(options.log_options, HighsLogType::kWarning, "Hipo: %s\n",
                 hipo::statusString((hipo::Status)status).c_str());
    return HighsStatus::kWarning;
  }

  else if (hipo.failed()) {
    highsLogUser(options.log_options, HighsLogType::kError, "Hipo: %s\n",
                 hipo::statusString((hipo::Status)status).c_str());
  }
  return HighsStatus::kError;
}

HighsStatus reportHipoCrossoverStatus(const HighsOptions& options,
                                      const ipx::Int status, bool is_qp) {
  if (is_qp && options.run_crossover == kHighsOnString) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Hipo: Crossover not available for QP\n");
    return HighsStatus::kWarning;
  }

  if (status == IPX_STATUS_not_run) {
    if (options.run_crossover == kHighsOnString) {
      // Warn if crossover not run and run_crossover option is "on"
      highsLogUser(options.log_options, HighsLogType::kWarning,
                   "Hipo: Crossover not run\n");
      return HighsStatus::kWarning;
    }
    return HighsStatus::kOk;
  } else if (status == IPX_STATUS_optimal) {
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "Hipo: Crossover optimal\n");
    return HighsStatus::kOk;
  } else if (status == IPX_STATUS_imprecise) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Hipo: Crossover imprecise\n");
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_primal_infeas) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Hipo: Crossover primal infeasible\n");
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_dual_infeas) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Hipo: Crossover dual infeasible\n");
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_user_interrupt) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Hipo: Crossover user interrupt\n");
    return HighsStatus::kOk;
  } else if (status == IPX_STATUS_time_limit) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Hipo: Crossover reached time limit\n");
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_iter_limit) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Hipo: Crossover reached iteration limit\n");
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_no_progress) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Hipo: Crossover no progress\n");
    return HighsStatus::kWarning;
  } else if (status == IPX_STATUS_failed) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Hipo: Crossover failed\n");
    return HighsStatus::kError;
  } else if (status == IPX_STATUS_debug) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Hipo: Crossover debug\n");
    return HighsStatus::kError;
  } else {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Hipo: Crossover unrecognised status\n");
    return HighsStatus::kError;
  }
}
