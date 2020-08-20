/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "test/DevKkt.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "presolve/PresolveUtils.h"

namespace presolve {
namespace dev_kkt_check {

constexpr int dev_print = 1;
constexpr double tol = 1e-07;

KktInfo initInfo() {
  KktInfo info;
  info.rules[KktCondition::kColBounds] =
      KktConditionDetails(KktCondition::kColBounds);
  info.rules[KktCondition::kPrimalFeasibility] =
      KktConditionDetails(KktCondition::kPrimalFeasibility);
  info.rules[KktCondition::kDualFeasibility] =
      KktConditionDetails(KktCondition::kDualFeasibility);
  info.rules[KktCondition::kComplementarySlackness] =
      KktConditionDetails(KktCondition::kComplementarySlackness);
  info.rules[KktCondition::kStationarityOfLagrangian] =
      KktConditionDetails(KktCondition::kStationarityOfLagrangian);
  info.rules[KktCondition::kBasicFeasibleSolution] =
      KktConditionDetails(KktCondition::kBasicFeasibleSolution);
  return info;
}

void checkPrimalBounds(const State& state, KktConditionDetails& details) {
  details.type = KktCondition::kColBounds;
  details.checked = 0;
  details.violated = 0;
  details.max_violation = 0.0;
  details.sum_violation_2 = 0.0;

  for (int i = 0; i < state.numCol; i++) {
    if (state.flagCol[i]) {
      details.checked++;
      double infeas = 0;

      if ((state.colLower[i] - state.colValue[i] > tol) ||
          (state.colValue[i] - state.colUpper[i] > tol)) {
        if (state.colLower[i] - state.colValue[i] > tol)
          infeas = state.colLower[i] - state.colValue[i];
        else
          infeas = state.colValue[i] - state.colUpper[i];

        if (dev_print == 1)
          std::cout << "Variable " << i
                    << " infeasible: lb=" << state.colLower[i]
                    << ", vaule=" << state.colValue[i]
                    << ",  ub=" << state.colUpper[i] << std::endl;

        details.violated++;
        details.sum_violation_2 += infeas * infeas;

        if (details.max_violation < infeas) details.max_violation = infeas;
      }
    }
  }
}

void checkPrimalFeasMatrix(const State& state, KktConditionDetails& details) {
  details.type = KktCondition::kPrimalFeasibility;
  details.checked = 0;
  details.violated = 0;
  details.max_violation = 0.0;
  details.sum_violation_2 = 0.0;

  for (int i = 0; i < state.numRow; i++) {
    if (state.flagRow[i]) {
      details.checked++;
      double rowV = state.rowValue[i];

      if (state.rowLower[i] < rowV && rowV < state.rowUpper[i]) continue;
      double infeas = 0;

      if (((rowV - state.rowLower[i]) < 0) &&
          (fabs(rowV - state.rowLower[i]) > tol)) {
        infeas = state.rowLower[i] - rowV;
        if (dev_print == 1)
          std::cout << "Row " << i << " infeasible: Row value=" << rowV
                    << "  L=" << state.rowLower[i]
                    << "  U=" << state.rowUpper[i] << std::endl;
      }

      if (((rowV - state.rowUpper[i]) > 0) &&
          (fabs(rowV - state.rowUpper[i]) > tol)) {
        infeas = rowV - state.rowUpper[i];
        if (dev_print == 1)
          std::cout << "Row " << i << " infeasible: Row value=" << rowV
                    << "  L=" << state.rowLower[i]
                    << "  U=" << state.rowUpper[i] << std::endl;
      }
      if (infeas > 0) {
        details.violated++;
        details.sum_violation_2 += infeas * infeas;

        if (details.max_violation < infeas) details.max_violation = infeas;
      }
    }
  }

  if (details.violated == 0) {
    if (dev_print == 1) std::cout << "Primal feasible.\n";
  } else {
    if (dev_print == 1) std::cout << "KKT check error: Primal infeasible.\n";
  }
}

void checkDualFeasibility(const State& state, KktConditionDetails& details) {
  details.type = KktCondition::kPrimalFeasibility;
  details.checked = 0;
  details.violated = 0;
  details.max_violation = 0.0;
  details.sum_violation_2 = 0.0;

  // check values of z_j are dual feasible
  for (int i = 0; i < state.numCol; i++) {
    if (state.flagCol[i]) {
      details.checked++;
      double infeas = 0;
      // j not in L or U
      if (state.colLower[i] <= -HIGHS_CONST_INF &&
          state.colUpper[i] >= HIGHS_CONST_INF) {
        if (fabs(state.colDual[i]) > tol) {
          if (dev_print == 1)
            std::cout << "Dual feasibility fail: l=-inf, x[" << i
                      << "]=" << state.colValue[i] << ", u=inf, z[" << i
                      << "]=" << state.colDual[i] << std::endl;
          infeas = fabs(state.colDual[i]);
        }
      }
      // j in L: x=l and l<u
      else if (state.colValue[i] == state.colLower[i] &&
               state.colLower[i] < state.colUpper[i]) {
        if (state.colDual[i] < 0 && fabs(state.colDual[i]) > tol) {
          if (dev_print == 1)
            std::cout << "Dual feasibility fail: l[" << i
                      << "]=" << state.colLower[i] << " = x[" << i
                      << "]=" << state.colValue[i] << ", z[" << i
                      << "]=" << state.colDual[i] << std::endl;
          infeas = fabs(state.colDual[i]);
        }
      }
      // j in U: x=u and l<u
      else if (state.colValue[i] == state.colUpper[i] &&
               state.colLower[i] < state.colUpper[i]) {
        if (state.colDual[i] > tol) {
          if (dev_print == 1)
            std::cout << "Dual feasibility fail: x[" << i << "]=" << i << "=u["
                      << i << "], z[" << i << "]=" << state.colDual[i]
                      << std::endl;
          infeas = fabs(state.colDual[i]);
        }
      }

      if (infeas > 0) {
        details.violated++;
        details.sum_violation_2 += infeas * infeas;

        if (details.max_violation < infeas) details.max_violation = infeas;
      }
    }
  }

  // check values of y_i are dual feasible
  for (int i = 0; i < state.numRow; i++) {
    if (state.flagRow[i]) {
      details.checked++;

      double rowV = state.rowValue[i];

      // L = Ax = U can be any sign
      if (fabs(state.rowLower[i] - rowV) < tol &&
          fabs(state.rowUpper[i] - rowV) < tol)
        continue;

      double infeas = 0;
      // L = Ax < U
      if (fabs(state.rowLower[i] - rowV) < tol && rowV < state.rowUpper[i]) {
        if (state.rowDual[i] > tol) {
          if (dev_print == 1)
            std::cout << "Dual feasibility fail for row " << i
                      << ": L= " << state.rowLower[i] << ", Ax=" << rowV
                      << ", U=" << state.rowUpper[i]
                      << ", y=" << state.rowDual[i] << std::endl;
          infeas = state.rowDual[i];
        }
      }
      // L < Ax = U
      else if (state.rowLower[i] < rowV &&
               fabs(rowV - state.rowUpper[i]) < tol) {
        if (state.rowDual[i] < -tol) {
          if (dev_print == 1)
            std::cout << "Dual feasibility fail for row " << i
                      << ": L= " << state.rowLower[i] << ", Ax=" << rowV
                      << ", U=" << state.rowUpper[i]
                      << ", y=" << state.rowDual[i] << std::endl;
          infeas = -state.rowDual[i];
        }
      }
      // L < Ax < U
      else if ((state.rowLower[i] < (rowV + tol)) &&
               (rowV < (state.rowUpper[i] + tol))) {
        if (fabs(state.rowDual[i]) > tol) {
          if (dev_print == 1)
            std::cout << "Dual feasibility fail for row " << i
                      << ": L= " << state.rowLower[i] << ", Ax=" << rowV
                      << ", U=" << state.rowUpper[i]
                      << ", y=" << state.rowDual[i] << std::endl;
          infeas = state.rowDual[i];
        }
      }
      if (infeas > 0) {
        details.violated++;
        details.sum_violation_2 += infeas * infeas;

        if (details.max_violation < infeas) details.max_violation = infeas;
      }
    }
  }

  if (details.violated == 0) {
    if (dev_print == 1) std::cout << "Dual feasible.\n";
  } else {
    if (dev_print == 1)
      std::cout << "KKT check error: Dual feasibility fail.\n";
  }
}

void checkComplementarySlackness(const State& state,
                                 KktConditionDetails& details) {
  details.type = KktCondition::kComplementarySlackness;
  details.checked = 0;
  details.violated = 0;
  details.max_violation = 0.0;
  details.sum_violation_2 = 0.0;

  for (int i = 0; i < state.numCol; i++) {
    if (state.flagCol[i]) {
      double infeas = 0;
      details.checked++;
      if (state.colLower[i] > -HIGHS_CONST_INF &&
          fabs(state.colValue[i] - state.colLower[i]) > tol) {
        if (fabs(state.colDual[i]) > tol &&
            fabs(state.colValue[i] - state.colUpper[i]) > tol) {
          if (dev_print)
            std::cout << "Comp. slackness fail: "
                      << "l[" << i << "]=" << state.colLower[i] << ", x[" << i
                      << "]=" << state.colValue[i] << ", z[" << i
                      << "]=" << state.colDual[i] << std::endl;
          infeas = fabs(state.colDual[i]);
        }
      }
      if (state.colUpper[i] < HIGHS_CONST_INF &&
          fabs(state.colUpper[i] - state.colValue[i]) > tol) {
        if (fabs(state.colDual[i]) > tol &&
            fabs(state.colValue[i] - state.colLower[i]) > tol) {
          if (dev_print == 1)
            std::cout << "Comp. slackness fail: x[" << i
                      << "]=" << state.colValue[i] << ", u[" << i
                      << "]=" << state.colUpper[i] << ", z[" << i
                      << "]=" << state.colDual[i] << std::endl;
          infeas = fabs(state.colDual[i]);
        }
      }

      if (infeas > 0) {
        details.violated++;
        details.sum_violation_2 += infeas * infeas;

        if (details.max_violation < infeas) details.max_violation = infeas;
      }
    }
  }

  if (details.violated == 0) {
    if (dev_print == 1) std::cout << "Complementary Slackness.\n";
  } else {
    if (dev_print == 1) std::cout << "KKT check error: Comp slackness fail.\n";
  }
}

void checkStationarityOfLagrangian(const State& state,
                                   KktConditionDetails& details) {
  details.type = KktCondition::kStationarityOfLagrangian;
  details.checked = 0;
  details.violated = 0;
  details.max_violation = 0.0;
  details.sum_violation_2 = 0.0;

  // A'y + c - z = 0
  for (int j = 0; j < state.numCol; j++) {
    if (state.flagCol[j]) {
      details.checked++;
      double infeas = 0;

      double lagrV = state.colCost[j] - state.colDual[j];
      for (int k = state.Astart[j]; k < state.Aend[j]; k++) {
        const int row = state.Aindex[k];
        assert(row >= 0 && row < state.numRow);
        if (state.flagRow[row])
          lagrV = lagrV + state.rowDual[row] * state.Avalue[k];
      }

      if (fabs(lagrV) > tol) {
        if (dev_print == 1)
          std::cout << "Column " << j
                    << " fails stationary of Lagrangian: dL/dx" << j << " = "
                    << lagrV << ", rather than zero." << std::endl;
        infeas = fabs(lagrV);
      }

      if (infeas > 0) {
        details.violated++;
        details.sum_violation_2 += infeas * infeas;

        if (details.max_violation < infeas) details.max_violation = infeas;
      }
    }
  }

  if (details.violated == 0) {
    if (dev_print == 1) std::cout << "Stationarity of Lagrangian.\n";
  } else {
    if (dev_print == 1)
      std::cout << "KKT check error: Lagrangian is not stationary.\n";
  }
}

void checkBasicFeasibleSolution(const State& state,
                                KktConditionDetails& details) {
  // Go over cols and check that the duals of basic values are zero.
  assert((int)state.col_status.size() == state.numCol);
  assert((int)state.colDual.size() == state.numCol);
  for (int j = 0; j < state.numCol; j++) {
    if (state.flagCol[j]) {
      details.checked++;
      double infeas = 0;
      if (state.col_status[j] == HighsBasisStatus::BASIC &&
          fabs(state.colDual[j]) > tol) {
        if (dev_print == 1)
          std::cout << "Col " << j << " is basic but has nonzero dual."
                    << std::endl;
        infeas = fabs(state.colDual[j]);
      }

      if (infeas > 0) {
        details.violated++;
        details.sum_violation_2 += infeas * infeas;

        if (details.max_violation < infeas) details.max_violation = infeas;
      }
    }
  }

  // Go over rows and check that the duals of basic values are zero.
  assert((int)state.row_status.size() == state.numRow);
  assert((int)state.rowDual.size() == state.numRow);
  for (int i = 0; i < state.numRow; i++) {
    if (state.flagRow[i]) {
      details.checked++;
      double infeas = 0;
      if (state.row_status[i] == HighsBasisStatus::BASIC &&
          fabs(state.rowDual[i]) > tol) {
        if (dev_print == 1)
          std::cout << "Row " << i << " is basic but has nonzero dual."
                    << std::endl;
        infeas = fabs(state.rowDual[i]);
      }
      if (infeas > 0) {
        details.violated++;
        details.sum_violation_2 += infeas * infeas;

        if (details.max_violation < infeas) details.max_violation = infeas;
      }
    }
  }

  if (details.violated == 0) {
    if (dev_print == 1) std::cout << "BFS." << std::endl;
  } else {
    if (dev_print == 1)
      std::cout << "BFS X Violated: " << details.violated << std::endl;
  }

  // check number of basic rows during postsolve.
  int current_n_rows = 0;
  int current_n_rows_basic = 0;
  int current_n_cols_basic = 0;

  for (int i = 0; i < state.numRow; i++) {
    if (state.flagRow[i]) current_n_rows++;

    if (state.flagRow[i] && (state.row_status[i] == HighsBasisStatus::BASIC ||
                             state.row_status[i] == HighsBasisStatus::SUPER))
      current_n_rows_basic++;
  }

  for (int i = 0; i < state.numCol; i++) {
    if (state.flagCol[i] && (state.col_status[i] == HighsBasisStatus::BASIC ||
                             state.col_status[i] == HighsBasisStatus::SUPER))
      current_n_cols_basic++;
  }

  bool holds = current_n_cols_basic + current_n_rows_basic == current_n_rows;
  if (!holds)
    std::cout << "BFS X Violated WRONG basis count: "
              << current_n_cols_basic + current_n_rows_basic << " "
              << current_n_rows << std::endl;
  assert(current_n_cols_basic + current_n_rows_basic == current_n_rows);
}

bool checkKkt(const State& state, KktInfo info) {
  if (state.numCol == 0) {
    std::cout << "KKT warning: empty problem" << std::endl;
    return true;
  }

  std::cout << std::endl;

  checkPrimalBounds(state, info.rules[KktCondition::kColBounds]);
  checkPrimalFeasMatrix(state, info.rules[KktCondition::kPrimalFeasibility]);
  checkDualFeasibility(state, info.rules[KktCondition::kDualFeasibility]);
  checkComplementarySlackness(
      state, info.rules[KktCondition::kComplementarySlackness]);
  checkStationarityOfLagrangian(
      state, info.rules[KktCondition::kStationarityOfLagrangian]);
  checkBasicFeasibleSolution(state,
                             info.rules[KktCondition::kBasicFeasibleSolution]);

  assert(info.rules.size() == 6);

  if (info.rules[KktCondition::kColBounds].violated == 0)
    info.pass_col_bounds = true;
  if (info.rules[KktCondition::kPrimalFeasibility].violated == 0)
    info.pass_primal_feas_matrix = true;
  if (info.rules[KktCondition::kDualFeasibility].violated == 0)
    info.pass_dual_feas = true;
  if (info.rules[KktCondition::kComplementarySlackness].violated == 0)
    info.pass_comp_slackness = true;
  if (info.rules[KktCondition::kStationarityOfLagrangian].violated == 0)
    info.pass_st_of_L = true;
  if (info.rules[KktCondition::kBasicFeasibleSolution].violated == 0)
    info.pass_bfs = true;

  if (info.pass_primal_feas_matrix && info.pass_col_bounds &&
      info.pass_dual_feas && info.pass_comp_slackness && info.pass_st_of_L)
    return true;

  return false;
}

}  // namespace dev_kkt_check
}  // namespace presolve