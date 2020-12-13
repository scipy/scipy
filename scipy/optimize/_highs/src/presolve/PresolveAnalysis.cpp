/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/PresolveAnalysis.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "presolve/PresolveAnalysis.h"

#include <limits>

namespace presolve {

void initializePresolveRuleInfo(std::vector<PresolveRuleInfo>& rules) {
  assert((int)rules.size() == 0);

  rules.push_back(PresolveRuleInfo(EMPTY_ROW, "Empty row", "EMR"));
  rules.push_back(PresolveRuleInfo(FIXED_COL, "Fixed col", "FXC"));
  rules.push_back(PresolveRuleInfo(SING_ROW, "Sing row", "SGR"));
  rules.push_back(PresolveRuleInfo(DOUBLETON_EQUATION, "Doubleton eq", "DEQ"));
  rules.push_back(
      PresolveRuleInfo(REMOVE_FORCING_CONSTRAINTS, "Rm forcing cs", "RFC"));
  rules.push_back(PresolveRuleInfo(FORCING_ROW, "Forcing row", "FRR"));
  rules.push_back(PresolveRuleInfo(REDUNDANT_ROW, "Redundant row", "RDR"));
  rules.push_back(
      PresolveRuleInfo(DOMINATED_ROW_BOUNDS, "Dom row bounds", "DRB"));
  rules.push_back(
      PresolveRuleInfo(REMOVE_COLUMN_SINGLETONS, "Remove col sing", "RCS"));
  rules.push_back(PresolveRuleInfo(FREE_SING_COL, "Free sing col", "FSC"));
  rules.push_back(
      PresolveRuleInfo(SING_COL_DOUBLETON_INEQ, "Sing col dbtn ineq", "SCD"));
  rules.push_back(
      PresolveRuleInfo(IMPLIED_FREE_SING_COL, "Impl free sing col", "IFS"));
  rules.push_back(
      PresolveRuleInfo(REMOVE_DOMINATED_COLUMNS, "Rm dom col", "RDC"));
  rules.push_back(PresolveRuleInfo(DOMINATED_COLS, "Dominated col", "DMC"));
  rules.push_back(
      PresolveRuleInfo(WEAKLY_DOMINATED_COLS, "Weakly dom col", "WDC"));
  rules.push_back(
      PresolveRuleInfo(DOMINATED_COL_BOUNDS, "Dom col bounds", "DCB"));
  rules.push_back(PresolveRuleInfo(EMPTY_COL, "Empty col", "EMC"));
  rules.push_back(PresolveRuleInfo(MATRIX_COPY, "Initialize matrix", "INM"));
  rules.push_back(PresolveRuleInfo(RESIZE_MATRIX, "Resize matrix", "RSM"));
  //
  rules.push_back(PresolveRuleInfo(RUN_PRESOLVERS, "Run Presolvers", "RPr"));
  rules.push_back(
      PresolveRuleInfo(REMOVE_ROW_SINGLETONS, "Rm row sing", "RRS"));
  rules.push_back(
      PresolveRuleInfo(REMOVE_DOUBLETON_EQUATIONS, "Rm dbleton eq", "RDE"));
  rules.push_back(PresolveRuleInfo(REMOVE_EMPTY_ROW, "Rm empty row", "RER"));
  //
  rules.push_back(
      PresolveRuleInfo(TOTAL_PRESOLVE_TIME, "Total presolve time", "TPT"));

  // Plus one for the total resize time.
  assert((int)rules.size() == PRESOLVE_RULES_COUNT);
}

void PresolveTimer::updateInfo() {
  for (PresolveRuleInfo& rule : rules_) {
    rule.total_time = timer_.read(rule.clock_id);
  }
}

}  // namespace presolve
