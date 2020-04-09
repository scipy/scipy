/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/PresolveAnalysis.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef PRESOLVE_PRESOLVE_ANALYSIS_H_
#define PRESOLVE_PRESOLVE_ANALYSIS_H_

#include <iostream>
#include <string>
#include <vector>

#include "util/HighsTimer.h"

enum PresolveRule {
  // Presolve rules.
  EMPTY_ROW,
  FIXED_COL,
  SING_ROW,
  DOUBLETON_EQUATION,
  FORCING_ROW,
  REDUNDANT_ROW,
  DOMINATED_ROW_BOUNDS,
  FREE_SING_COL,
  SING_COL_DOUBLETON_INEQ,
  IMPLIED_FREE_SING_COL,
  DOMINATED_COLS,
  WEAKLY_DOMINATED_COLS,
  DOMINATED_COL_BOUNDS,
  EMPTY_COL,
  // HTICK_PRE_DUPLICATE_ROWS,
  // HTICK_PRE_DUPLICATE_COLUMNS,

  // For timing.
  MATRIX_COPY,
  RESIZE_MATRIX,

  TOTAL_PRESOLVE_TIME,
  // Number of presolve rules.
  PRESOLVE_RULES_COUNT,

  // Items required by postsolve
  DOUBLETON_EQUATION_ROW_BOUNDS_UPDATE,
  DOUBLETON_EQUATION_X_ZERO_INITIALLY,
  DOUBLETON_EQUATION_NEW_X_NONZERO,
  DOUBLETON_EQUATION_NEW_X_ZERO_AR_UPDATE,
  DOUBLETON_EQUATION_NEW_X_ZERO_A_UPDATE,
  SING_COL_DOUBLETON_INEQ_SECOND_SING_COL,
  FORCING_ROW_VARIABLE
};

struct PresolveRuleInfo {
  PresolveRuleInfo(PresolveRule id, std::string name, std::string name_ch3)
      : rule_id(id),
        rule_name(std::move(name)),
        rule_name_ch3(std::move(name_ch3)) {}
  PresolveRule rule_id;

  std::string rule_name;
  std::string rule_name_ch3;

  int count_applied = 0;
  int rows_removed = 0;
  int cols_removed = 0;

  int clock_id = 0;
  double total_time = 0;
};

void initializePresolveRuleInfo(std::vector<PresolveRuleInfo>& rules);

class PresolveTimer {
 public:
  PresolveTimer(HighsTimer& timer) : timer_(timer) {
    initializePresolveRuleInfo(rules_);
    for (PresolveRuleInfo& rule : rules_) {
      int clock_id =
          timer_.clock_def(rule.rule_name.c_str(), rule.rule_name_ch3.c_str());
      rule.clock_id = clock_id;
    }
  }

  void recordStart(PresolveRule rule) {
    assert(rule >= 0 && rule < PRESOLVE_RULES_COUNT);
    assert((int)rules_.size() == (int)PRESOLVE_RULES_COUNT);
    timer_.start(rules_[rule].clock_id);
  }

  void recordFinish(PresolveRule rule) {
    assert(rule >= 0 && rule < PRESOLVE_RULES_COUNT);
    assert((int)rules_.size() == (int)PRESOLVE_RULES_COUNT);
    timer_.stop(rules_[rule].clock_id);

    if (rule == TOTAL_PRESOLVE_TIME)
      total_time_ = timer_.read(rules_[rule].clock_id);
  }

  void addChange(PresolveRule rule) {
    assert(rule >= 0 && rule < PRESOLVE_RULES_COUNT);
    assert((int)rules_.size() == (int)PRESOLVE_RULES_COUNT);
    rules_[rule].count_applied++;
  }

  void increaseCount(bool row_count, PresolveRule rule) {
    assert(rule >= 0 && rule < PRESOLVE_RULES_COUNT);
    assert((int)rules_.size() == (int)PRESOLVE_RULES_COUNT);
    if (row_count)
      rules_[rule].rows_removed++;
    else
      rules_[rule].cols_removed++;
  }

  void reportClocks() {
    std::vector<int> clocks(PRESOLVE_RULES_COUNT - 1);
    for (int id = 0; id < PRESOLVE_RULES_COUNT - 1; id++) {
      assert(rules_[id].rule_id == id);
      clocks[id] = rules_[id].clock_id;
    }
    std::cout << std::endl;
    timer_.report("grep-Presolve", clocks);
    std::cout << std::endl;
  }

  void updateInfo();
  double getTotalTime() { return total_time_; }

 private:
  HighsTimer& timer_;
  std::vector<PresolveRuleInfo> rules_;

  double total_time_ = 0.0;
};

#endif
