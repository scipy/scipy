/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/ICrashX.h"

#include <algorithm>
#include <iostream>

#include "HConfig.h"
#include "io/HighsIO.h"
#include "ipm/IpxWrapper.h"

HighsStatus callCrossover(const HighsOptions& options, const HighsLp& lp,
                          HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& highs_callback) {
  ipx::Int num_col, num_row;
  double offset;
  std::vector<ipx::Int> Ap, Ai;
  std::vector<double> objective, col_lb, col_ub, Av, rhs;
  std::vector<char> constraint_type;

  fillInIpxData(lp, num_col, num_row, offset, objective, col_lb, col_ub, Ap, Ai,
                Av, rhs, constraint_type);
  // if (res != IpxStatus::OK) return HighsStatus::kError;

  const HighsLogOptions& log_options = options.log_options;

  ipx::Parameters parameters;
  parameters.run_crossover = 1;  // 1 = "on"
  parameters.crash_basis = 1;    // 0 = slack basis; 1 = crash basis
  parameters.display = 1;
  if (!options.output_flag) parameters.display = 0;
  // Modify parameters.debug according to log_dev_level
  parameters.debug = 0;
  if (options.log_dev_level == kHighsLogDevLevelDetailed) {
    // Default options.log_dev_level setting is kHighsLogDevLevelNone, yielding
    // default setting debug = 0
  } else if (options.log_dev_level == kHighsLogDevLevelInfo) {
    parameters.debug = 2;
  } else if (options.log_dev_level == kHighsLogDevLevelVerbose) {
    parameters.debug = 4;
  }
  parameters.highs_logging = true;
  parameters.log_options = &options.log_options;
  parameters.time_limit = options.time_limit;

  ipx::LpSolver lps;
  lps.SetParameters(parameters);

  // Set pointer to any callback
  lps.SetCallback(&highs_callback);

  ipx::Int load_status = lps.LoadModel(
      num_col, offset, objective.data(), col_lb.data(), col_ub.data(), num_row,
      Ap.data(), Ai.data(), Av.data(), rhs.data(), constraint_type.data());
  if (load_status != 0) {
    highsLogUser(log_options, HighsLogType::kError,
                 "Error loading ipx model\n");
    return HighsStatus::kError;
  }

  // IPX's columns are the LP's plus one slack per row with two
  // finite, distinct bounds (fillInIpxData; coefficient -1, the row
  // becomes '='), and its rows are the LP's less any free rows.
  const bool use_duals =
      highs_solution.dual_valid &&
      (HighsInt)highs_solution.col_dual.size() == lp.num_col_ &&
      (HighsInt)highs_solution.row_dual.size() == lp.num_row_;

  // IPX accepts only an exactly complementary start: a column with z
  // > 0 must sit at its lower bound, z < 0 at its upper; a row's
  // slack must be at its bound wherever its dual is nonzero. An
  // approximate (first-order or early-stopped) point is not, and IPX
  // would refuse it (invalid vector).  Per coordinate, either move
  // the primal onto the bound or zero the dual, whichever is the
  // smaller relative change.
  double cost_scale = 1.0;
  for (const double c : objective)
    cost_scale = std::max(cost_scale, std::fabs(c));
  HighsInt num_snapped = 0;
  HighsInt num_dual_zeroed = 0;
  double sum_snapped = 0;
  double sum_dual_zeroed = 0;
  const double use_primal_tolerance =
      0;                                // options.primal_feasibility_tolerance;
  const double use_dual_tolerance = 0;  // options.dual_feasibility_tolerance;

  std::vector<double> x(num_col, 0.0);
  std::vector<double> y;
  std::vector<double> z;

  // Duals of maximizations must be negated
  const HighsInt sense = lp.sense_ == ObjSense::kMinimize ? 1 : -1;

  // Boxed rows l <= a^Tx <= u are converted using explicit slack s to
  // a^Tx - s = 0; l <= s <= u
  //
  // Hence, i \in [lp.num_col_, num_col), Ap[i + 1] == Ap[i] + 1,
  // Av[Ap[i]] = -1 and rhs[Ai[Ap[i]]] = 0
  //
  // Since IPX excludes any free rows, skip them in mapping the LP's
  // row duals onto y
  HighsInt ipx_num_row = 0;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (lp.row_lower_[iRow] <= -kHighsInf && lp.row_upper_[iRow] >= kHighsInf) {
      // Free row not in IPX: any dual it has should be small, but is ignored
      double dual = sense * highs_solution.row_dual[iRow];
      if (dual) {
        num_dual_zeroed++;
        sum_dual_zeroed += std::fabs(dual);
      }
    } else {
      if (use_duals) y.push_back(sense * highs_solution.row_dual[iRow]);
      ipx_num_row++;
    }
  }

  for (HighsInt i = 0; i < lp.num_col_; i++)
    x[i] =
        std::min(std::max(highs_solution.col_value[i], col_lb[i]), col_ub[i]);
  if (use_duals) {
    z.assign(num_col, 0.0);
    for (HighsInt i = 0; i < lp.num_col_; i++)
      z[i] = sense * highs_solution.col_dual[i];
    for (HighsInt i = lp.num_col_; i < num_col; i++) {
      assert(Av[Ap[i]] == -1);
      z[i] = y[Ai[Ap[i]]];
    }
  }

  auto complement = [&](double& value, double& dual, const double lower,
                        const double upper) {
    if (dual > use_dual_tolerance) {
      if (value == lower) return;
      if (std::isfinite(lower) &&
          (value - lower) / (1 + std::fabs(lower)) <= dual / cost_scale) {
        sum_snapped += std::fabs(value - lower);
        num_snapped++;
        value = lower;
      } else {
        sum_dual_zeroed += std::fabs(dual);
        num_dual_zeroed++;
        dual = 0;
      }
    } else if (dual < -use_dual_tolerance) {
      if (value == upper) return;
      if (std::isfinite(upper) &&
          (upper - value) / (1 + std::fabs(upper)) <= -dual / cost_scale) {
        sum_snapped += std::fabs(value - upper);
        num_snapped++;
        value = upper;
      } else {
        sum_dual_zeroed += std::fabs(dual);
        num_dual_zeroed++;
        dual = 0;
      }
    }
  };
  if (use_duals) {
    for (HighsInt i = 0; i < lp.num_col_; i++)
      complement(x[i], z[i], col_lb[i], col_ub[i]);
  }

  // Compute the row activites for the original LP
  std::vector<HighsCDouble> activity(num_row, HighsCDouble(0.0));
  for (HighsInt i = 0; i < lp.num_col_; i++) {
    for (HighsInt p = Ap[i]; p < Ap[i + 1]; ++p)
      activity[Ai[p]] += Av[p] * x[i];
  }

  if (num_col > lp.num_col_) {
    for (HighsInt i = lp.num_col_; i < num_col; i++) {
      const HighsInt row = Ai[Ap[i]];
      const double value = double(activity[row]);
      x[i] = std::min(std::max(value, col_lb[i]), col_ub[i]);
      if (use_duals) {
        complement(x[i], z[i], col_lb[i], col_ub[i]);
        assert(rhs[row] == 0);
        y[row] = z[i];
      }
    }
  }

  // Build slack variables from rhs-A*x but subject to sign conditions.
  // Complete the calculation of the activities for IPX
  for (HighsInt i = lp.num_col_; i < num_col; i++) activity[Ai[Ap[i]]] -= x[i];

  std::vector<double> slack(num_row);
  for (HighsInt i = 0; i < num_row; i++) {
    slack[i] = double(HighsCDouble(rhs[i]) - activity[i]);
    switch (constraint_type[i]) {
      case '=':
        slack[i] = 0.0;
        break;
      case '<':
        // slack in [0, inf), reduced cost -y
        slack[i] = std::max(slack[i], 0.0);
        if (use_duals) {
          double dual = -y[i];
          complement(slack[i], dual, 0.0, kHighsInf);
          y[i] = -dual;
        }
        break;
      case '>':
        // slack in (-inf, 0], reduced cost -y
        slack[i] = std::min(slack[i], 0.0);
        if (use_duals) {
          double dual = -y[i];
          complement(slack[i], dual, -kHighsInf, 0.0);
          y[i] = -dual;
        }
        break;
    }
  }

  HighsInt crossover_status;
  if (use_duals) {
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Calling IPX crossover with primal and dual values");
    if (num_snapped || num_dual_zeroed) {
      highsLogUser(log_options, HighsLogType::kInfo,
                   ": made complementary with\n");
      if (num_snapped)
        highsLogUser(log_options, HighsLogType::kInfo,
                     "   %" HIGHSINT_FORMAT
                     " primal value%s moved to a bound (total changes = %g)\n",
                     num_snapped, highsIntToPlural(num_snapped).c_str(),
                     sum_snapped);
      if (num_dual_zeroed)
        highsLogUser(log_options, HighsLogType::kInfo,
                     "   %" HIGHSINT_FORMAT
                     " dual%s zeroed (total changes = %g)\n",
                     num_dual_zeroed, highsIntToPlural(num_dual_zeroed).c_str(),
                     sum_dual_zeroed);
    } else {
      highsLogUser(log_options, HighsLogType::kInfo, "\n");
    }
    crossover_status = lps.CrossoverFromStartingPoint(x.data(), slack.data(),
                                                      y.data(), z.data());
  } else {
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Calling IPX crossover with only primal values\n");
    crossover_status = lps.CrossoverFromStartingPoint(x.data(), slack.data(),
                                                      nullptr, nullptr);
  }

  if (crossover_status != 0) {
    highsLogUser(log_options, HighsLogType::kError,
                 "IPX crossover error: flag = %d\n", (int)crossover_status);
    return HighsStatus::kError;
  }
  ipx::Info ipx_info = lps.GetInfo();
  highs_info.crossover_iteration_count += ipx_info.updates_crossover;
  const bool imprecise_solution =
      ipx_info.status_crossover == IPX_STATUS_imprecise;
  if (ipx_info.status_crossover != IPX_STATUS_optimal &&
      ipx_info.status_crossover != IPX_STATUS_imprecise &&
      ipx_info.status_crossover != IPX_STATUS_time_limit) {
    highsLogUser(log_options, HighsLogType::kError,
                 "IPX crossover failed: status = %d\n",
                 int(ipx_info.status_crossover));
    return HighsStatus::kError;
  }
  if (ipx_info.status_crossover == IPX_STATUS_time_limit) {
    model_status = HighsModelStatus::kTimeLimit;
    return HighsStatus::kWarning;
  }

  // Get basis
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
    highsLogUser(log_options, HighsLogType::kError,
                 "IPX crossover getting basic solution: flag = %d\n",
                 int(errflag));
    return HighsStatus::kError;
  }

  // Convert the IPX basic solution to a HiGHS basic solution
  HighsStatus status = ipxBasicSolutionToHighsBasicSolution(
      options.log_options, lp, rhs, constraint_type, ipx_solution, highs_basis,
      highs_solution);

  if (status != HighsStatus::kOk) {
    highsLogUser(
        log_options, HighsLogType::kError,
        "Failed to convert IPX basic solution to Highs basic solution\n");
    return HighsStatus::kError;
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
