/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolve.h"

namespace presolve {

HPresolve::Result HPresolve::presolveRuleTest(
    HighsPostsolveStack& postsolve_stack) {
  assert(options->presolve_rule_test);
  if (options->presolve_rule_test == kPresolveRuleColStuffing) {
    return presolveRuleTestColStuffing(postsolve_stack);
  } else if (options->presolve_rule_test == kPresolveRuleParallelRowsAndCols) {
    return presolveRuleTestParallelRowsAndCols(postsolve_stack);
  } else if (options->presolve_rule_test == kPresolveRuleProbing) {
    return presolveRuleTestProbing(postsolve_stack);
  } else if (options->presolve_rule_test == kPresolveRuleFourierMotzkin) {
    return presolveRuleTestFourierMotzkin(postsolve_stack);
  } else if (options->presolve_rule_test == kPresolveRuleAggregator) {
    return presolveRuleTestAggregator(postsolve_stack);
  }
  return Result::kOk;
}
HPresolve::Result HPresolve::presolveRuleTestColStuffing(
    HighsPostsolveStack& postsolve_stack) {
  assert(options->presolve_rule_test == kPresolveRuleColStuffing);
  assert(model->num_row_ == 1);
  highsLogUser(options->log_options, HighsLogType::kInfo,
               "HPresolve::presolveRuleTestColStuffing\n");
  HPresolve::Result result = Result::kOk;
  for (HighsInt col = 0; col < model->num_col_; col++) {
    if (colDeleted[col]) continue;
    result = singletonColStuffing(postsolve_stack, col);
    if (result != Result::kOk) return result;
  }
  highsLogUser(options->log_options, HighsLogType::kInfo,
               "HPresolve::presolveRuleTestColStuffing: Stuffing removed %d "
               "rows and %d columns\n",
               static_cast<int>(numDeletedRows),
               static_cast<int>(numDeletedCols));
  // Possibly remove the row
  return rowPresolve(postsolve_stack, 0);
}
HPresolve::Result HPresolve::presolveRuleTestParallelRowsAndCols(
    HighsPostsolveStack& postsolve_stack) {
  assert(options->presolve_rule_test == kPresolveRuleParallelRowsAndCols);
  highsLogUser(options->log_options, HighsLogType::kInfo,
               "HPresolve::presolveRuleTestParallelRowsAndCols\n");
  return detectParallelRowsAndCols(postsolve_stack);
}
HPresolve::Result HPresolve::presolveRuleTestProbing(
    HighsPostsolveStack& postsolve_stack) {
  assert(options->presolve_rule_test == kPresolveRuleProbing);
  highsLogUser(options->log_options, HighsLogType::kInfo,
               "HPresolve::presolveRuleTestProbing\n");
  if (mipsolver == nullptr) return Result::kStopped;
  return runProbing(postsolve_stack);
}
HPresolve::Result HPresolve::presolveRuleTestFourierMotzkin(
    HighsPostsolveStack& postsolve_stack) {
  assert(options->presolve_rule_test == kPresolveRuleFourierMotzkin);
  highsLogUser(options->log_options, HighsLogType::kInfo,
               "HPresolve::presolveRuleTestFourierMotzkin\n");
  HighsInt numColsEliminated;
  return fourierMotzkin(postsolve_stack, numColsEliminated);
}
HPresolve::Result HPresolve::presolveRuleTestAggregator(
    HighsPostsolveStack& postsolve_stack) {
  assert(options->presolve_rule_test == kPresolveRuleAggregator);
  highsLogUser(options->log_options, HighsLogType::kInfo,
               "HPresolve::presolveRuleTestAggregator\n");
  // Compute the implied row dual bounds, which may make rows dual implied
  // free, so that substitution opportunities are identified
  for (HighsInt col = 0; col < model->num_col_; col++)
    if (!colDeleted[col]) updateRowDualImpliedBounds(col);
  return aggregator(postsolve_stack);
}
}  // namespace presolve
