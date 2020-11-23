/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ipm/IpxWrapper.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef IPM_IPX_WRAPPER_H_
#define IPM_IPX_WRAPPER_H_

#include <algorithm>

#include "ipm/IpxSolution.h"
#include "ipm/IpxStatus.h"
#include "ipm/ipx/include/ipx_status.h"
#include "ipm/ipx/src/lp_solver.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsSolution.h"

IpxStatus fillInIpxData(const HighsLp& lp, ipx::Int& num_col,
                        std::vector<double>& obj, std::vector<double>& col_lb,
                        std::vector<double>& col_ub, ipx::Int& num_row,
                        std::vector<ipx::Int>& Ap, std::vector<ipx::Int>& Ai,
                        std::vector<double>& Ax, std::vector<double>& rhs,
                        std::vector<char>& constraint_type) {
  num_col = lp.numCol_;
  num_row = lp.numRow_;

  // For each row with both a lower and an upper bound introduce one new column
  // so num_col may increase. Ignore each free row so num_row may decrease.
  // lba <= a'x <= uba becomes
  // a'x-s = 0 and lba <= s <= uba.

  // For each row with bounds on both sides introduce explicit slack and
  // transfer bounds.
  assert(lp.rowLower_.size() == (unsigned int)num_row);
  assert(lp.rowUpper_.size() == (unsigned int)num_row);

  std::vector<int> general_bounded_rows;
  std::vector<int> free_rows;

  for (int row = 0; row < num_row; row++)
    if (lp.rowLower_[row] < lp.rowUpper_[row] &&
        lp.rowLower_[row] > -HIGHS_CONST_INF &&
        lp.rowUpper_[row] < HIGHS_CONST_INF)
      general_bounded_rows.push_back(row);
    else if (lp.rowLower_[row] <= -HIGHS_CONST_INF &&
             lp.rowUpper_[row] >= HIGHS_CONST_INF)
      free_rows.push_back(row);

  const int num_slack = general_bounded_rows.size();

  // For each row except free rows add entry to char array and set up rhs
  // vector
  rhs.reserve(num_row);
  constraint_type.reserve(num_row);

  for (int row = 0; row < num_row; row++) {
    if (lp.rowLower_[row] > -HIGHS_CONST_INF &&
        lp.rowUpper_[row] >= HIGHS_CONST_INF) {
      rhs.push_back(lp.rowLower_[row]);
      constraint_type.push_back('>');
    } else if (lp.rowLower_[row] <= -HIGHS_CONST_INF &&
               lp.rowUpper_[row] < HIGHS_CONST_INF) {
      rhs.push_back(lp.rowUpper_[row]);
      constraint_type.push_back('<');
    } else if (lp.rowLower_[row] == lp.rowUpper_[row]) {
      rhs.push_back(lp.rowUpper_[row]);
      constraint_type.push_back('=');
    } else if (lp.rowLower_[row] > -HIGHS_CONST_INF &&
               lp.rowUpper_[row] < HIGHS_CONST_INF) {
      // general bounded
      rhs.push_back(0);
      constraint_type.push_back('=');
    }
  }

  std::vector<int> reduced_rowmap(lp.numRow_, -1);
  if (free_rows.size() > 0) {
    int counter = 0;
    int findex = 0;
    for (int row = 0; row < lp.numRow_; row++) {
      if (free_rows[findex] == row) {
        findex++;
        continue;
      } else {
        reduced_rowmap[row] = counter;
        counter++;
      }
    }
  } else {
    for (int k = 0; k < lp.numRow_; k++) reduced_rowmap[k] = k;
  }
  num_row -= free_rows.size();
  num_col += num_slack;

  std::vector<int> sizes(num_col, 0);

  for (int col = 0; col < lp.numCol_; col++)
    for (int k = lp.Astart_[col]; k < lp.Astart_[col + 1]; k++) {
      int row = lp.Aindex_[k];
      if (lp.rowLower_[row] > -HIGHS_CONST_INF ||
          lp.rowUpper_[row] < HIGHS_CONST_INF)
        sizes[col]++;
    }
  // Copy Astart and Aindex to ipx::Int array.
  int nnz = lp.Aindex_.size();
  Ap.resize(num_col + 1);
  Ai.reserve(nnz + num_slack);
  Ax.reserve(nnz + num_slack);

  // Set starting points of original and newly introduced columns.
  Ap[0] = 0;
  for (int col = 0; col < lp.numCol_; col++) {
    Ap[col + 1] = Ap[col] + sizes[col];
    //    printf("Struc Ap[%2d] = %2d; Al[%2d] = %2d\n", col, (int)Ap[col], col,
    //    (int)sizes[col]);
  }
  for (int col = lp.numCol_; col < (int)num_col; col++) {
    Ap[col + 1] = Ap[col] + 1;
    //    printf("Slack Ap[%2d] = %2d\n", col, (int)Ap[col]);
  }
  //  printf("Fictn Ap[%2d] = %2d\n", (int)num_col, (int)Ap[num_col]);
  for (int k = 0; k < nnz; k++) {
    int row = lp.Aindex_[k];
    if (lp.rowLower_[row] > -HIGHS_CONST_INF ||
        lp.rowUpper_[row] < HIGHS_CONST_INF) {
      Ai.push_back(reduced_rowmap[lp.Aindex_[k]]);
      Ax.push_back(lp.Avalue_[k]);
    }
  }

  for (int k = 0; k < num_slack; k++) {
    Ai.push_back((ipx::Int)general_bounded_rows[k]);
    Ax.push_back(-1);
  }

  // Column bound vectors.
  col_lb.resize(num_col);
  col_ub.resize(num_col);
  for (int col = 0; col < lp.numCol_; col++) {
    if (lp.colLower_[col] <= -HIGHS_CONST_INF)
      col_lb[col] = -INFINITY;
    else
      col_lb[col] = lp.colLower_[col];

    if (lp.colUpper_[col] >= HIGHS_CONST_INF)
      col_ub[col] = INFINITY;
    else
      col_ub[col] = lp.colUpper_[col];
  }
  for (int slack = 0; slack < num_slack; slack++) {
    const int row = general_bounded_rows[slack];
    col_lb[lp.numCol_ + slack] = lp.rowLower_[row];
    col_ub[lp.numCol_ + slack] = lp.rowUpper_[row];
  }

  obj.resize(num_col);
  for (int col = 0; col < lp.numCol_; col++) {
    obj[col] = (int)lp.sense_ * lp.colCost_[col];
  }
  obj.insert(obj.end(), num_slack, 0);
  /*
  for (int col = 0; col < num_col; col++)
    printf("Col %2d: [%11.4g, %11.4g] Cost = %11.4g; Start = %d\n", col,
  col_lb[col], col_ub[col], obj[col], (int)Ap[col]); for (int row = 0; row <
  num_row; row++) printf("Row %2d: RHS = %11.4g; Type = %d\n", row, rhs[row],
  constraint_type[row]); for (int col = 0; col < num_col; col++) { for (int el =
  Ap[col]; el < Ap[col+1]; el++) { printf("El %2d: [%2d, %11.4g]\n", el,
  (int)Ai[el], Ax[el]);
    }
  }
  */

  return IpxStatus::OK;
}

HighsStatus reportIpxSolveStatus(const HighsOptions& options,
                                 const ipx::Int solve_status,
                                 const ipx::Int error_flag) {
  if (solve_status == IPX_STATUS_solved) {
    HighsLogMessage(options.logfile, HighsMessageType::INFO, "Ipx: Solved");
    return HighsStatus::OK;
  } else if (solve_status == IPX_STATUS_stopped) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING, "Ipx: Stopped");
    return HighsStatus::Warning;
  } else if (solve_status == IPX_STATUS_invalid_input) {
    if (error_flag == IPX_ERROR_argument_null) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Ipx: Invalid input - argument_null");
      return HighsStatus::Error;
    } else if (error_flag == IPX_ERROR_invalid_dimension) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Ipx: Invalid input - invalid dimension");
      return HighsStatus::Error;
    } else if (error_flag == IPX_ERROR_invalid_matrix) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Ipx: Invalid input - invalid matrix");
      return HighsStatus::Error;
    } else if (error_flag == IPX_ERROR_invalid_vector) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Ipx: Invalid input - invalid vector");
      return HighsStatus::Error;
    } else if (error_flag == IPX_ERROR_invalid_basis) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Ipx: Invalid input - invalid basis");
      return HighsStatus::Error;
    } else {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                      "Ipx: Invalid input - unrecognised error");
      return HighsStatus::Error;
    }
  } else if (solve_status == IPX_STATUS_out_of_memory) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Ipx: Out of memory");
    return HighsStatus::Error;
  } else if (solve_status == IPX_STATUS_internal_error) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Ipx: Internal error %d", (int)error_flag);
    return HighsStatus::Error;
  } else {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Ipx: unrecognised solve status = %d", (int)solve_status);
    return HighsStatus::Error;
  }
  return HighsStatus::Error;
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
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Ipx: %s not run", method_name.c_str());
    return HighsStatus::Warning;
  } else if (status == IPX_STATUS_optimal) {
    HighsLogMessage(options.logfile, HighsMessageType::INFO, "Ipx: %s optimal",
                    method_name.c_str());
    return HighsStatus::OK;
  } else if (status == IPX_STATUS_imprecise) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Ipx: %s imprecise", method_name.c_str());
    return HighsStatus::Warning;
  } else if (status == IPX_STATUS_primal_infeas) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Ipx: %s primal infeasible", method_name.c_str());
    return HighsStatus::Warning;
  } else if (status == IPX_STATUS_dual_infeas) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Ipx: %s dual infeasible", method_name.c_str());
    return HighsStatus::Warning;
  } else if (status == IPX_STATUS_time_limit) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Ipx: %s reached time limit", method_name.c_str());
    return HighsStatus::Warning;
  } else if (status == IPX_STATUS_iter_limit) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Ipx: %s reached iteration limit", method_name.c_str());
    return HighsStatus::Warning;
  } else if (status == IPX_STATUS_no_progress) {
    HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                    "Ipx: %s no progress", method_name.c_str());
    return HighsStatus::Warning;
  } else if (status == IPX_STATUS_failed) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR, "Ipx: %s failed",
                    method_name.c_str());
    return HighsStatus::Error;
  } else if (status == IPX_STATUS_debug) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR, "Ipx: %s debug",
                    method_name.c_str());
    return HighsStatus::Error;
  } else {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "Ipx: %s unrecognised status", method_name.c_str());
    return HighsStatus::Error;
  }
  return HighsStatus::Error;
}

bool ipxStatusError(const bool status_error, const HighsOptions& options,
                    std::string message, const int value = -1) {
  if (status_error) {
    if (value < 0) {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR, "Ipx: %s",
                      message.c_str());
    } else {
      HighsLogMessage(options.logfile, HighsMessageType::ERROR, "Ipx: %s %d",
                      message.c_str(), value);
    }
    fflush(NULL);
  }
  assert(!status_error);
  return status_error;
}

bool illegalIpxSolvedStatus(ipx::Info& ipx_info, const HighsOptions& options) {
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
                     "solved  status_ipm should not be IPX_STATUS_time_limit");
  // Cannot solve and reach iteration limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_iter_limit, options,
                     "solved  status_ipm should not be IPX_STATUS_iter_limit");
  // Cannot solve and make no progress
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_no_progress, options,
                     "solved  status_ipm should not be IPX_STATUS_no_progress");
  // Cannot solve and failed
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_failed, options,
                     "solved  status_ipm should not be IPX_STATUS_failed");
  // Cannot solve and debug
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_debug, options,
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
          ipx_info.status_crossover == IPX_STATUS_primal_infeas, options,
          "solved  status_crossover should not be IPX_STATUS_primal_infeas");
  // Cannot solve with dual infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_dual_infeas, options,
          "solved  status_crossover should not be IPX_STATUS_dual_infeas");
  // Cannot solve and reach time limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_time_limit, options,
          "solved  status_crossover should not be IPX_STATUS_time_limit");
  // Cannot solve and reach time limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_iter_limit, options,
          "solved  status_crossover should not be IPX_STATUS_iter_limit");
  // Cannot solve and make no progress
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_no_progress, options,
          "solved  status_crossover should not be IPX_STATUS_no_progress");
  // Cannot solve and failed
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_failed, options,
          "solved  status_crossover should not be IPX_STATUS_failed");
  // Cannot solve and debug
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_crossover == IPX_STATUS_debug, options,
                     "solved  status_crossover should not be IPX_STATUS_debug");
  return found_illegal_status;
}

bool illegalIpxStoppedIpmStatus(ipx::Info& ipx_info,
                                const HighsOptions& options) {
  bool found_illegal_status = false;
  // Cannot stop and be optimal
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_optimal, options,
                     "stopped status_ipm should not be IPX_STATUS_optimal");
  // Cannot stop and be imprecise
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_imprecise, options,
                     "stopped status_ipm should not be IPX_STATUS_imprecise");
  // Cannot stop with primal infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_ipm == IPX_STATUS_primal_infeas, options,
          "stopped status_ipm should not be IPX_STATUS_primal_infeas");
  // Cannot stop with dual infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_dual_infeas, options,
                     "stopped status_ipm should not be IPX_STATUS_dual_infeas");
  // Can stop with time limit
  // Can stop with iter limit
  // Can stop with no progress
  // Cannot stop and failed - should be error return earlier
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_failed, options,
                     "stopped status_ipm should not be IPX_STATUS_failed");
  // Cannot stop and debug - should be error return earlier
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_ipm == IPX_STATUS_debug, options,
                     "stopped status_ipm should not be IPX_STATUS_debug");
  return found_illegal_status;
}

bool illegalIpxStoppedCrossoverStatus(ipx::Info& ipx_info,
                                      const HighsOptions& options) {
  bool found_illegal_status = false;
  // Cannot stop and be optimal
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_optimal, options,
          "stopped status_crossover should not be IPX_STATUS_optimal");
  // Cannot stop and be imprecise
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_imprecise, options,
          "stopped status_crossover should not be IPX_STATUS_imprecise");
  // Cannot stop with primal infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_primal_infeas, options,
          "stopped status_crossover should not be IPX_STATUS_primal_infeas");
  // Cannot stop with dual infeasibility
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_dual_infeas, options,
          "stopped status_crossover should not be IPX_STATUS_dual_infeas");
  // Cannot stop and reach iteration limit
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_iter_limit, options,
          "stopped status_crossover should not be IPX_STATUS_iter_limit");
  // Can stop and reach time limit
  // Cannot stop with no_progress
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_no_progress, options,
          "stopped status_crossover should not be IPX_STATUS_no_progress");
  // Cannot stop and failed - should be error return earlier
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(
          ipx_info.status_crossover == IPX_STATUS_failed, options,
          "stopped status_crossover should not be IPX_STATUS_failed");
  // Cannot stop and debug - should be error return earlier
  found_illegal_status =
      found_illegal_status ||
      ipxStatusError(ipx_info.status_crossover == IPX_STATUS_debug, options,
                     "stopped status_crossover should not be IPX_STATUS_debug");
  return found_illegal_status;
}

void reportIpmNoProgress(const HighsOptions& options,
                         const ipx::Info& ipx_info) {
  HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                  "No progress: primal objective value       = %11.4g",
                  ipx_info.pobjval);
  HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                  "No progress: max absolute primal residual = %11.4g",
                  ipx_info.abs_presidual);
  HighsLogMessage(options.logfile, HighsMessageType::WARNING,
                  "No progress: max absolute   dual residual = %11.4g",
                  ipx_info.abs_dresidual);
}

HighsStatus analyseIpmNoProgress(const ipx::Info& ipx_info,
                                 const ipx::Parameters& parameters,
                                 HighsModelStatus& unscaled_model_status) {
  if (ipx_info.abs_presidual > parameters.ipm_feasibility_tol) {
    // Looks like the LP is infeasible
    unscaled_model_status = HighsModelStatus::PRIMAL_INFEASIBLE;
    return HighsStatus::OK;
  } else if (ipx_info.abs_dresidual > parameters.ipm_optimality_tol) {
    // Looks like the LP is unbounded
    unscaled_model_status = HighsModelStatus::PRIMAL_UNBOUNDED;
    return HighsStatus::OK;
  } else if (ipx_info.pobjval < -HIGHS_CONST_INF) {
    // Looks like the LP is unbounded
    unscaled_model_status = HighsModelStatus::PRIMAL_UNBOUNDED;
    return HighsStatus::OK;
  } else {
    // Don't know
    unscaled_model_status = HighsModelStatus::SOLVE_ERROR;
    return HighsStatus::Error;
  }
  return HighsStatus::Warning;
}
HighsStatus solveLpIpx(const HighsOptions& options, HighsTimer& timer,
                       const HighsLp& lp, bool& imprecise_solution,
                       HighsBasis& highs_basis, HighsSolution& highs_solution,
                       HighsIterationCounts& iteration_counts,
                       HighsModelStatus& unscaled_model_status,
                       HighsSolutionParams& unscaled_solution_params) {
  imprecise_solution = false;
  resetModelStatusAndSolutionParams(unscaled_model_status,
                                    unscaled_solution_params, options);
  // Create the LpSolver instance
  ipx::LpSolver lps;
  // Set IPX parameters
  //
  // Cannot set internal IPX parameters directly since they are
  // private, so create instance of parameters
  ipx::Parameters parameters;
  // Set IPX parameters from options
  //
  // Set display according to output
  parameters.display = 1;
  if (options.output == NULL) parameters.display = 0;
  // Set debug according to message_level
  parameters.debug = 0;
  if (options.message_level % HighsPrintMessageLevel::ML_MINIMAL == 0) {
    // Default options.message_level is
    // HighsPrintMessageLevel::ML_MINIMAL, yielding default setting
    // debug = 0
    parameters.debug = 0;
  } else if (options.message_level % HighsPrintMessageLevel::ML_DETAILED == 0) {
    parameters.debug = 3;
  } else if (options.message_level % HighsPrintMessageLevel::ML_VERBOSE == 0) {
    parameters.debug = 4;
  }
  // Just test feasibility and optimality tolerances for now
  // ToDo Set more parameters
  parameters.ipm_feasibility_tol =
      min(unscaled_solution_params.primal_feasibility_tolerance,
          unscaled_solution_params.dual_feasibility_tolerance);

  parameters.ipm_optimality_tol = options.ipm_optimality_tolerance;
  // Determine the run time allowed for IPX
  parameters.time_limit = options.time_limit - timer.readRunHighsClock();
  parameters.ipm_maxiter = options.ipm_iteration_limit - iteration_counts.ipm;
  // Determine if crossover is to be run or not
  parameters.crossover = options.run_crossover;
  if (!parameters.crossover) {
    // If crossover is not run, then set crossover_start to -1 so that
    // IPX can terminate according to its feasibility and optimality
    // tolerances
    parameters.crossover_start = -1;
  }

  // Set the internal IPX parameters
  lps.SetParameters(parameters);

  ipx::Int num_col, num_row;
  std::vector<ipx::Int> Ap, Ai;
  std::vector<double> objective, col_lb, col_ub, Av, rhs;
  std::vector<char> constraint_type;
  IpxStatus result = fillInIpxData(lp, num_col, objective, col_lb, col_ub,
                                   num_row, Ap, Ai, Av, rhs, constraint_type);
  if (result != IpxStatus::OK) return HighsStatus::Error;
  HighsLogMessage(options.logfile, HighsMessageType::INFO,
                  "IPX model has %d rows, %d columns and %d nonzeros",
                  (int)num_row, (int)num_col, (int)Ap[num_col]);

  ipx::Int solve_status =
      lps.Solve(num_col, &objective[0], &col_lb[0], &col_ub[0], num_row, &Ap[0],
                &Ai[0], &Av[0], &rhs[0], &constraint_type[0]);

  // Get solver and solution information.
  // Struct ipx_info defined in ipx/include/ipx_info.h
  ipx::Info ipx_info = lps.GetInfo();
  iteration_counts.ipm += (int)ipx_info.iter;
  //  iteration_counts.crossover += (int)ipx_info.updates_crossover;
  iteration_counts.crossover += (int)ipx_info.pushes_crossover;

  // If not solved...
  if (solve_status != IPX_STATUS_solved) {
    const HighsStatus solve_return_status =
        reportIpxSolveStatus(options, solve_status, ipx_info.errflag);
    // Return error if IPX solve error has occurred
    if (solve_return_status == HighsStatus::Error) {
      unscaled_model_status = HighsModelStatus::SOLVE_ERROR;
      return HighsStatus::Error;
    }
  }
  bool ipm_status = true;
  const HighsStatus ipm_return_status =
      reportIpxIpmCrossoverStatus(options, ipx_info.status_ipm, ipm_status);
  ipm_status = false;
  const HighsStatus crossover_return_status = reportIpxIpmCrossoverStatus(
      options, ipx_info.status_crossover, ipm_status);
  // Return error if IPX IPM or crossover error has occurred
  if (ipm_return_status == HighsStatus::Error ||
      crossover_return_status == HighsStatus::Error) {
    unscaled_model_status = HighsModelStatus::SOLVE_ERROR;
    return HighsStatus::Error;
  }
  // Should only reach here if Solve() returned IPX_STATUS_solved or
  // IPX_STATUS_stopped
  if (ipxStatusError(
          solve_status != IPX_STATUS_solved &&
              solve_status != IPX_STATUS_stopped,
          options, "solve_status should be solved or stopped here but value is",
          (int)solve_status))
    return HighsStatus::Error;

  if (solve_status == IPX_STATUS_stopped) {
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
      return HighsStatus::Error;
    // Can stop and reach time limit
    if (ipx_info.status_crossover == IPX_STATUS_time_limit) {
      unscaled_model_status = HighsModelStatus::REACHED_TIME_LIMIT;
      return HighsStatus::Warning;
    }
    //========
    // For IPM
    //========
    //
    // Note that IPX can stop with IPM optimal, imprecise,
    // primal_infeas or dual_infeas, due to crossover stopping with
    // time limit, and this is why crossover returns are tested first
    if (illegalIpxStoppedIpmStatus(ipx_info, options))
      return HighsStatus::Error;
    // Can stop with time limit
    // Can stop with iter limit
    // Can stop with no progress
    if (ipx_info.status_ipm == IPX_STATUS_time_limit) {
      unscaled_model_status = HighsModelStatus::REACHED_TIME_LIMIT;
      return HighsStatus::Warning;
    } else if (ipx_info.status_ipm == IPX_STATUS_iter_limit) {
      unscaled_model_status = HighsModelStatus::REACHED_ITERATION_LIMIT;
      return HighsStatus::Warning;
    } else if (ipx_info.status_ipm == IPX_STATUS_no_progress) {
      reportIpmNoProgress(options, ipx_info);
      return analyseIpmNoProgress(ipx_info, lps.GetParameters(),
                                  unscaled_model_status);
    }
  }
  // Should only reach here if Solve() returned IPX_STATUS_solved
  if (ipxStatusError(solve_status != IPX_STATUS_solved, options,
                     "solve_status should be solved here but value is",
                     (int)solve_status))
    return HighsStatus::Error;
  // Return error if solved status settings occur that JAJH doesn't
  // think should happen
  if (illegalIpxSolvedStatus(ipx_info, options)) return HighsStatus::Error;
  //==============
  // For crossover
  //==============
  // Can be not run
  // Can solve and be optimal
  // Can solve and be imprecise
  //========
  // For IPX
  //========
  // Can solve and be optimal
  // Can solve and be imprecise
  // Can solve and be primal_infeas
  // Can solve and be dual_infeas
  if (ipx_info.status_ipm == IPX_STATUS_primal_infeas) {
    unscaled_model_status = HighsModelStatus::PRIMAL_INFEASIBLE;
    return HighsStatus::OK;
  } else if (ipx_info.status_ipm == IPX_STATUS_dual_infeas) {
    unscaled_model_status = HighsModelStatus::PRIMAL_UNBOUNDED;
    return HighsStatus::OK;
  }

  // Should only reach here if crossover is not run, optimal or imprecise
  if (ipxStatusError(ipx_info.status_crossover != IPX_STATUS_not_run &&
                         ipx_info.status_crossover != IPX_STATUS_optimal &&
                         ipx_info.status_crossover != IPX_STATUS_imprecise,
                     options,
                     "crossover status should be not run, optimal or imprecise "
                     "but value is",
                     (int)ipx_info.status_crossover))
    return HighsStatus::Error;

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

  lps.GetInteriorSolution(&x[0], &xl[0], &xu[0], &slack[0], &y[0], &zl[0],
                          &zu[0]);

  // Basic solution depends on crossover being run
  const bool have_basic_solution =
      ipx_info.status_crossover != IPX_STATUS_not_run;
  imprecise_solution = ipx_info.status_crossover == IPX_STATUS_imprecise;
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
    lps.GetBasicSolution(
        &ipx_solution.ipx_col_value[0], &ipx_solution.ipx_row_value[0],
        &ipx_solution.ipx_row_dual[0], &ipx_solution.ipx_col_dual[0],
        &ipx_solution.ipx_row_status[0], &ipx_solution.ipx_col_status[0]);

    // Convert the IPX basic solution to a HiGHS basic solution
    ipxBasicSolutionToHighsBasicSolution(options.logfile, lp, rhs,
                                         constraint_type, ipx_solution,
                                         highs_basis, highs_solution);
  } else {
    ipxSolutionToHighsSolution(options.logfile, lp, rhs, constraint_type,
                               num_col, num_row, x, slack, highs_solution);
    highs_basis.valid_ = false;
  }
  HighsStatus return_status;
  if (imprecise_solution) {
    unscaled_model_status = HighsModelStatus::NOTSET;
    return_status = HighsStatus::Warning;
  } else {
    unscaled_model_status = HighsModelStatus::OPTIMAL;
    unscaled_solution_params.primal_status =
        PrimalDualStatus::STATUS_FEASIBLE_POINT;
    // Currently only have a dual solution if there is a basic solution
    if (have_basic_solution)
      unscaled_solution_params.dual_status =
          PrimalDualStatus::STATUS_FEASIBLE_POINT;
    return_status = HighsStatus::OK;
  }
  double objective_function_value = lp.offset_;
  for (int iCol = 0; iCol < lp.numCol_; iCol++)
    objective_function_value +=
        highs_solution.col_value[iCol] * lp.colCost_[iCol];
  unscaled_solution_params.objective_function_value = objective_function_value;
  if (highs_basis.valid_)
    getPrimalDualInfeasibilities(lp, highs_basis, highs_solution,
                                 unscaled_solution_params);
  return return_status;
}
#endif
