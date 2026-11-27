/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolve.h"

namespace presolve {

HPresolve::Result HPresolve::zeroCostSingleton(
    HighsPostsolveStack& postsolve_stack, HighsInt col) {
  assert(this->allow_rule_[kPresolveRuleZeroCostSingleton]);
  // For a double-sided row b_0 <= a^Tx + cs <= b_1,
  // where s is a singleton continuous column with 0 cost,
  // relax s out as its value can be determined in postsolve.
  // The row may now admit additional reductions afterwards.
  // Dual fixing already handles single-sided row case with fixings.
  if (model->integrality_[col] == HighsVarType::kInteger ||
      model->col_cost_[col] != 0.0 || colsize[col] != 1) {
    return Result::kOk;
  }
  assert(!colDeleted[col]);

  HighsInt nzPos = colhead[col];
  HighsInt row = Arow[nzPos];
  assert(!rowDeleted[row]);
  // Row must be ranged, but why can't it be an equation?
  const bool was_equation = isEquation(row);
  if (!isRanged(row)) return Result::kOk;

  double coef = Avalue[nzPos];

  if (std::abs(coef) == kHighsInf) return Result::kOk;

  const bool logging_on = analysis_.logging_on_;
  if (logging_on)
    analysis_.startPresolveRuleLog(kPresolveRuleZeroCostSingleton);
  storeRow(row);

  double lb = model->col_lower_[col];
  double ub = model->col_upper_[col];
  double change_from_col_lb = coef * lb;
  double change_from_col_ub = coef * ub;

  double newRowLower =
      model->row_lower_[row] - std::max(change_from_col_lb, change_from_col_ub);
  double newRowUpper =
      model->row_upper_[row] - std::min(change_from_col_lb, change_from_col_ub);

  postsolve_stack.zeroCostSingleton(row, col, model->row_lower_[row],
                                    model->row_upper_[row], newRowLower,
                                    newRowUpper, lb, ub, coef, getStoredRow());

  changeRowLower(row, newRowLower);
  changeRowUpper(row, newRowUpper);
  if (was_equation && newRowLower != newRowUpper &&
      eqiters[row] != equations.end()) {
    equations.erase(eqiters[row]);
    eqiters[row] = equations.end();
  }

  // Delete the singleton column
  markColDeleted(col);
  unlink(nzPos);

  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleZeroCostSingleton);

  return checkLimits(postsolve_stack);
}

}  // namespace presolve
