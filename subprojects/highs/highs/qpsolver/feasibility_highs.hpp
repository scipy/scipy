/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef __SRC_LIB_FEASIBILITYHIGHS_HPP__
#define __SRC_LIB_FEASIBILITYHIGHS_HPP__

#include "Highs.h"
#include "qpsolver/a_asm.hpp"

static void computeStartingPointHighs(
    Instance& instance, Settings& settings, Statistics& stats,
    QpModelStatus& modelstatus, QpHotstartInformation& result,
    HighsModelStatus& highs_model_status, HighsBasis& highs_basis,
    HighsSolution& highs_solution, HighsTimer& timer) {
  bool have_starting_point = false;
  const bool debug_report = false;
  if (settings.allow_hot_start && highs_solution.value_valid) {
    // #1350 add primal_feasibility_tolerance to settings
    const double primal_feasibility_tolerance = settings.lambda_zero_threshold;

    HighsInt num_var_infeasibilities = 0;
    double max_var_infeasibility = 0;
    double sum_var_infeasibilities = 0;
    HighsInt num_con_infeasibilities = 0;
    double max_con_infeasibility = 0;
    double sum_con_infeasibilities = 0;
    double max_con_residual = 0;
    double sum_con_residuals = 0;

    assessQpPrimalFeasibility(
        instance, primal_feasibility_tolerance, highs_solution.col_value,
        highs_solution.row_value, num_var_infeasibilities,
        max_var_infeasibility, sum_var_infeasibilities, num_con_infeasibilities,
        max_con_infeasibility, sum_con_infeasibilities, max_con_residual,
        sum_con_residuals);

    if (debug_report)
      printf(
          "computeStartingPointHighs highs_solution has (num / max / sum) "
          "var (%d / %g / %g) and "
          "con (%d / %g / %g) infeasibilities "
          "with (max = %g; sum = %g) residuals\n",
          int(num_var_infeasibilities), max_var_infeasibility,
          sum_var_infeasibilities, int(num_con_infeasibilities),
          max_con_infeasibility, sum_con_infeasibilities, max_con_residual,
          sum_con_residuals);
    have_starting_point = num_var_infeasibilities == 0 &&
                          num_con_infeasibilities == 0 && highs_basis.valid;
  }
  // compute initial feasible point
  HighsBasis use_basis;
  HighsSolution use_solution;
  auto isfreevar = [&](const HighsInt iVar) {
    return instance.var_lo[iVar] == -kHighsInf &&
           instance.var_up[iVar] == kHighsInf;
  };
  if (have_starting_point) {
    use_basis = highs_basis;
    use_solution = highs_solution;
    // Have to assume that the supplied basis is feasible
    modelstatus = QpModelStatus::kNotset;
  } else {
    Highs highs;

    // set HiGHS to be silent
    highs.setOptionValue("output_flag", false);
    highs.setOptionValue("presolve", kHighsOnString);
    // Set the residual time limit
    const double use_time_limit =
        std::max(settings.time_limit - timer.read(), 0.001);
    highs.setOptionValue("time_limit", use_time_limit);

    HighsLp lp;
    lp.a_matrix_.index_ = instance.A.mat.index;
    lp.a_matrix_.start_ = instance.A.mat.start;
    lp.a_matrix_.value_ = instance.A.mat.value;
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.col_cost_.assign(instance.num_var, 0.0);
    // lp.col_cost_ = runtime.instance.c.value;
    lp.col_lower_ = instance.var_lo;
    lp.col_upper_ = instance.var_up;
    lp.row_lower_ = instance.con_lo;
    lp.row_upper_ = instance.con_up;
    lp.num_col_ = instance.num_var;
    lp.num_row_ = instance.num_con;

    // create artificial bounds for free variables: false by default
    assert(!settings.phase1boundfreevars);
    if (settings.phase1boundfreevars) {
      for (HighsInt i = 0; i < instance.num_var; i++) {
        if (isfreevar(i)) {
          lp.col_lower_[i] = -1E5;
          lp.col_upper_[i] = 1E5;
        }
      }
    }

    highs.passModel(lp);
    // Make free variables basic: false by default
    assert(!settings.phase1movefreevarsbasic);
    if (settings.phase1movefreevarsbasic) {
      HighsBasis basis;
      basis.alien = true;  // Set true when basis is instantiated
      for (HighsInt i = 0; i < instance.num_con; i++) {
        basis.row_status.push_back(HighsBasisStatus::kNonbasic);
      }

      for (HighsInt i = 0; i < instance.num_var; i++) {
        // make free variables basic
        if (instance.var_lo[i] == -kHighsInf &&
            instance.var_up[i] == kHighsInf) {
          // free variable
          basis.col_status.push_back(HighsBasisStatus::kBasic);
        } else {
          basis.col_status.push_back(HighsBasisStatus::kNonbasic);
        }
      }

      highs.setBasis(basis);

      highs.setOptionValue("simplex_strategy", kSimplexStrategyPrimal);
    }

    HighsStatus status = highs.run();
    if (status == HighsStatus::kError) {
      modelstatus = QpModelStatus::kError;
      return;
    }

    HighsModelStatus phase1stat = highs.getModelStatus();
    switch (phase1stat) {
      case HighsModelStatus::kOptimal:
        modelstatus = QpModelStatus::kNotset;
        break;
      case HighsModelStatus::kInfeasible:
        modelstatus = QpModelStatus::kInfeasible;
        break;
      case HighsModelStatus::kTimeLimit:
        modelstatus = QpModelStatus::kTimeLimit;
        break;
      case HighsModelStatus::kInterrupt:
        modelstatus = QpModelStatus::kInterrupt;
        break;
      default:
        modelstatus = QpModelStatus::kError;
    }

    stats.phase1_iterations = highs.getInfo().simplex_iteration_count;

    if (modelstatus != QpModelStatus::kNotset) return;

    // Should only get here if feasibility problem is solved to
    // optimality - hence there is a feasible basis
    assert(phase1stat == HighsModelStatus::kOptimal);

    use_basis = highs.getBasis();
    use_solution = highs.getSolution();
  }

  HighsInt num_small_x0 = 0;
  HighsInt num_small_ra = 0;
  const double zero_activity_tolerance = have_starting_point ? 0 : 1e-4;
  QpVector x0(instance.num_var);
  QpVector ra(instance.num_con);
  for (HighsInt i = 0; i < x0.dim; i++) {
    if (fabs(use_solution.col_value[i]) > zero_activity_tolerance) {
      x0.value[i] = use_solution.col_value[i];
      x0.index[x0.num_nz++] = i;
    } else if (fabs(use_solution.col_value[i]) > 0) {
      num_small_x0++;
    }
  }

  for (HighsInt i = 0; i < ra.dim; i++) {
    if (fabs(use_solution.row_value[i]) > zero_activity_tolerance) {
      ra.value[i] = use_solution.row_value[i];
      ra.index[ra.num_nz++] = i;
    } else if (fabs(use_solution.row_value[i]) > 0) {
      num_small_ra++;
    }
  }
  if (debug_report && num_small_x0 + num_small_ra)
    printf(
        "feasibility_highs has %d small col values and %d small row values\n",
        int(num_small_x0), int(num_small_ra));
  std::vector<HighsInt> initial_active;
  std::vector<HighsInt> initial_inactive;
  std::vector<BasisStatus> initial_status;

  const HighsInt num_highs_basis_status =
      HighsInt(HighsBasisStatus::kNonbasic) + 1;
  std::vector<HighsInt> debug_row_status_count;
  debug_row_status_count.assign(num_highs_basis_status, 0);
  for (HighsInt i = 0; i < HighsInt(use_basis.row_status.size()); i++) {
    HighsBasisStatus status = use_basis.row_status[i];
    debug_row_status_count[HighsInt(status)]++;
    if (status == HighsBasisStatus::kLower) {
      initial_active.push_back(i);
      initial_status.push_back(BasisStatus::kActiveAtLower);
    } else if (status == HighsBasisStatus::kUpper) {
      initial_active.push_back(i);
      initial_status.push_back(BasisStatus::kActiveAtUpper);
    } else if (status == HighsBasisStatus::kZero) {
      // Shouldn't happen, since free rows are basic in a logical
      // basis and remain basic, or are removed by presolve and
      // restored as basic in postsolve
      assert(111 == 222);
      // That said, a free row that is nonbasic in the Highs basis
      // must be counted as inactive in the QP basis for accounting
      // purposes
      initial_inactive.push_back(i);
    } else if (status != HighsBasisStatus::kBasic) {
      assert(status == HighsBasisStatus::kNonbasic);
      // Surely an error, but not a problem before, since simplex
      // solver cannot return a HighsBasisStatus::kNonbasic
      // variable. Does matter now, since a saved QP basis will
      // generally have such variables.
      //
      //      initial_inactive.push_back(instance.num_con + i);
      //
      // A HighsBasisStatus::kNonbasic variable corresponds one-to-one
      // with being inactive in the QP basis
      initial_inactive.push_back(i);
    } else {
      assert(status == HighsBasisStatus::kBasic);
    }
  }

  std::vector<HighsInt> debug_col_status_count;
  debug_col_status_count.assign(num_highs_basis_status, 0);
  for (HighsInt i = 0; i < HighsInt(use_basis.col_status.size()); i++) {
    HighsBasisStatus status = use_basis.col_status[i];
    debug_col_status_count[HighsInt(status)]++;
    if (status == HighsBasisStatus::kLower) {
      if (isfreevar(i)) {
        initial_inactive.push_back(instance.num_con + i);
      } else {
        initial_active.push_back(instance.num_con + i);
        initial_status.push_back(BasisStatus::kActiveAtLower);
      }

    } else if (status == HighsBasisStatus::kUpper) {
      if (isfreevar(i)) {
        initial_inactive.push_back(instance.num_con + i);
      } else {
        initial_active.push_back(instance.num_con + i);
        initial_status.push_back(BasisStatus::kActiveAtUpper);
      }

    } else if (status == HighsBasisStatus::kZero) {
      initial_inactive.push_back(instance.num_con + i);
    } else if (status != HighsBasisStatus::kBasic) {
      assert(status == HighsBasisStatus::kNonbasic);
      initial_inactive.push_back(instance.num_con + i);
    } else {
      assert(status == HighsBasisStatus::kBasic);
    }
  }

  if (debug_report) {
    printf("QP solver initial basis: (Lo / Bs / Up / Ze / Nb) for cols (");
    for (HighsInt k = 0; k < num_highs_basis_status; k++)
      printf("%s%d", k == 0 ? "" : " / ", int(debug_col_status_count[k]));
    printf(") and rows (");
    for (HighsInt k = 0; k < num_highs_basis_status; k++)
      printf("%s%d", k == 0 ? "" : " / ", int(debug_row_status_count[k]));
    printf(")\n");
  }

  // This used to be an assert
  if ((HighsInt)(initial_active.size() + initial_inactive.size()) !=
      instance.num_var) {
    modelstatus = QpModelStatus::kError;
    return;
  }

  if (!have_starting_point) {
    // When starting from a feasible basis, there will generally be
    // inactive variables in the basis that aren't free
    for (HighsInt ia : initial_inactive) {
      if (ia < instance.num_con) {
        // printf("free row %d\n", (int)ia);
        assert(instance.con_lo[ia] == -kHighsInf);
        assert(instance.con_up[ia] == kHighsInf);
      } else {
        // printf("free col %d\n", (int)ia);
        assert(instance.var_lo[ia - instance.num_con] == -kHighsInf);
        assert(instance.var_up[ia - instance.num_con] == kHighsInf);
      }
    }
  }

  result.status = initial_status;
  result.active = initial_active;
  result.inactive = initial_inactive;
  result.primal = x0;
  result.rowact = ra;
}

#endif
