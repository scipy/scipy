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

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsTimer.h"

namespace presolve {

using std::min;

constexpr double inf = std::numeric_limits<double>::infinity();

enum PresolveRule {
  // Presolve rules.
  EMPTY_ROW,
  FIXED_COL,
  SING_ROW,
  DOUBLETON_EQUATION,
  REMOVE_FORCING_CONSTRAINTS,
  FORCING_ROW,
  REDUNDANT_ROW,
  DOMINATED_ROW_BOUNDS,
  REMOVE_COLUMN_SINGLETONS,
  FREE_SING_COL,
  SING_COL_DOUBLETON_INEQ,
  IMPLIED_FREE_SING_COL,
  REMOVE_DOMINATED_COLUMNS,
  DOMINATED_COLS,
  WEAKLY_DOMINATED_COLS,
  DOMINATED_COL_BOUNDS,
  EMPTY_COL,
  // HTICK_PRE_DUPLICATE_ROWS,
  // HTICK_PRE_DUPLICATE_COLUMNS,

  // For timing.
  MATRIX_COPY,
  RESIZE_MATRIX,

  RUN_PRESOLVERS,
  REMOVE_ROW_SINGLETONS,
  REMOVE_DOUBLETON_EQUATIONS,
  REMOVE_EMPTY_ROW,

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

enum presolveNumerics {
  INCONSISTENT_BOUNDS,
  FIXED_COLUMN,
  DOUBLETON_EQUATION_BOUND,
  DOUBLETON_INEQUALITY_BOUND,
  SMALL_MATRIX_VALUE,
  EMPTY_ROW_BOUND,
  DOMINATED_COLUMN,
  WEAKLY_DOMINATED_COLUMN,
  PRESOLVE_NUMERICS_COUNT
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

struct numericsRecord {
  std::string name;
  double tolerance;
  int num_test;
  int num_zero_true;
  int num_tol_true;
  int num_10tol_true;
  int num_clear_true;
  double min_positive_true;
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

  std::vector<numericsRecord> presolve_numerics;

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
    std::vector<int> clocks;
    for (int id = 0; id < PRESOLVE_RULES_COUNT - 1; id++) {
      assert(rules_[id].rule_id == id);
      if (id == RUN_PRESOLVERS) continue;
      if (id == REMOVE_ROW_SINGLETONS) continue;
      if (id == REMOVE_DOUBLETON_EQUATIONS) continue;
      if (id == REMOVE_EMPTY_ROW) continue;
      clocks.push_back(rules_[id].clock_id);
    }
    int ideal_time_rule;
    double ideal_time;
    ideal_time_rule = TOTAL_PRESOLVE_TIME;
    ideal_time = getRuleTime(ideal_time_rule);
    std::cout << std::endl;
    timer_.report_tl("grep-Presolve", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    clocks.push_back(rules_[RUN_PRESOLVERS].clock_id);
    clocks.push_back(rules_[RESIZE_MATRIX].clock_id);
    std::cout << std::endl;
    timer_.report_tl("grep-Presolve", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    ideal_time_rule = RUN_PRESOLVERS;
    ideal_time = getRuleTime(ideal_time_rule);
    clocks.push_back(rules_[REMOVE_ROW_SINGLETONS].clock_id);
    clocks.push_back(rules_[REMOVE_FORCING_CONSTRAINTS].clock_id);
    clocks.push_back(rules_[REMOVE_COLUMN_SINGLETONS].clock_id);
    clocks.push_back(rules_[REMOVE_DOUBLETON_EQUATIONS].clock_id);
    clocks.push_back(rules_[REMOVE_DOMINATED_COLUMNS].clock_id);
    timer_.report_tl("grep-Presolve", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    ideal_time_rule = REMOVE_FORCING_CONSTRAINTS;
    ideal_time = getRuleTime(ideal_time_rule);
    clocks.push_back(rules_[REMOVE_EMPTY_ROW].clock_id);
    clocks.push_back(rules_[FORCING_ROW].clock_id);
    clocks.push_back(rules_[REDUNDANT_ROW].clock_id);
    clocks.push_back(rules_[DOMINATED_ROW_BOUNDS].clock_id);
    timer_.report_tl("grep--RmFrcCs", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    ideal_time_rule = REMOVE_COLUMN_SINGLETONS;
    ideal_time = getRuleTime(ideal_time_rule);
    clocks.push_back(rules_[FREE_SING_COL].clock_id);
    clocks.push_back(rules_[SING_COL_DOUBLETON_INEQ].clock_id);
    clocks.push_back(rules_[IMPLIED_FREE_SING_COL].clock_id);
    timer_.report_tl("grep-RmColSng", clocks, ideal_time, 0);
    std::cout << std::endl;

    clocks.clear();
    ideal_time_rule = REMOVE_DOMINATED_COLUMNS;
    ideal_time = getRuleTime(ideal_time_rule);
    clocks.push_back(rules_[DOMINATED_COLS].clock_id);
    clocks.push_back(rules_[WEAKLY_DOMINATED_COLS].clock_id);
    timer_.report_tl("grep-RmDomCol", clocks, ideal_time, 0);
    std::cout << std::endl;
  }

  void initialiseNumericsRecord(int record, std::string name,
                                const double tolerance) {
    // Make sure that the tolerance has been set to a positive value
    assert(tolerance > 0);
    numericsRecord& numerics_record = presolve_numerics[record];
    numerics_record.name = name;
    numerics_record.tolerance = tolerance;
    numerics_record.num_test = 0;
    numerics_record.num_zero_true = 0;
    numerics_record.num_tol_true = 0;
    numerics_record.num_10tol_true = 0;
    numerics_record.num_clear_true = 0;
    numerics_record.min_positive_true = HIGHS_CONST_INF;
  }

  void updateNumericsRecord(int record, const double value) {
    numericsRecord& numerics_record = presolve_numerics[record];
    double tolerance = numerics_record.tolerance;
    numerics_record.num_test++;
    if (value < 0) return;
    if (value == 0) {
      numerics_record.num_zero_true++;
    } else if (value <= tolerance) {
      numerics_record.num_tol_true++;
    } else if (value <= 10 * tolerance) {
      numerics_record.num_10tol_true++;
    } else {
      numerics_record.num_clear_true++;
    }
    if (value > 0)
      numerics_record.min_positive_true =
          min(value, numerics_record.min_positive_true);
  }

  void reportNumericsRecord(const numericsRecord& numerics_record) {
    if (!numerics_record.num_test) return;
    printf(
        "%-26s: tolerance =%6.1g: Zero =%9d; Tol =%9d; 10Tol =%9d; Clear =%9d; "
        "MinPositive =%7.2g; Tests =%9d\n",
        numerics_record.name.c_str(), numerics_record.tolerance,
        numerics_record.num_zero_true, numerics_record.num_tol_true,
        numerics_record.num_10tol_true, numerics_record.num_clear_true,
        numerics_record.min_positive_true, numerics_record.num_test);
  }

  void reportNumericsCsvRecord(const numericsRecord& numerics_record) {
    printf(",%d,%d,%d", numerics_record.num_zero_true,
           numerics_record.num_tol_true + numerics_record.num_10tol_true,
           numerics_record.num_clear_true);
  }

  void reportNumericsRecords() {
    assert((int)presolve_numerics.size() == PRESOLVE_NUMERICS_COUNT);
    if (presolve_numerics.size() < PRESOLVE_NUMERICS_COUNT) return;
    printf("Presolve numerics analysis for %s:\n\n", model_name.c_str());
    for (int record = 0; record < PRESOLVE_NUMERICS_COUNT; record++)
      reportNumericsRecord(presolve_numerics[record]);
    printf("grep_presolveNumerics:,%s", model_name.c_str());
    for (int record = 0; record < PRESOLVE_NUMERICS_COUNT; record++)
      reportNumericsCsvRecord(presolve_numerics[record]);
    printf("\n\n");
  }

  void updateInfo();
  double getTotalTime() { return total_time_; }

  HighsTimer& timer_;

  double getRuleTime(const int rule_id) {
    return timer_.read(rules_[rule_id].clock_id);
  }

  inline double getTime() { return timer_.readRunHighsClock(); }

  inline bool reachLimit() {
    if (time_limit == inf || time_limit <= 0) return false;
    if (getTime() < time_limit) return false;
    return true;
  }

  double start_time = 0.0;
  double time_limit = 0.0;
  std::string model_name;

 private:
  std::vector<PresolveRuleInfo> rules_;

  double total_time_ = 0.0;
};

}  // namespace presolve

#endif
