/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "lp_data/HighsModelUtils.h"
#include "presolve/HPresolve.h"
#include "presolve/PresolveTimer.h"

void HPresolveAnalysis::setup(const HighsLp* model_,
                              const HighsOptions* options_,
                              const HighsInt& numDeletedRows_,
                              const HighsInt& numDeletedCols_,
                              HighsTimer* timer) {
  model = model_;
  options = options_;
  numDeletedRows = &numDeletedRows_;
  numDeletedCols = &numDeletedCols_;

  timer_ = timer;
  const bool lp_presolve = !model_->isMip() || options->solve_relaxation;
  analyse_presolve_time_ =
      kHighsAnalysisLevelPresolveTime & options->highs_analysis_level &&
      lp_presolve;
  if (analyse_presolve_time_) {
    HighsTimerClock clock;
    clock.timer_pointer_ = timer_;
    PresolveTimer presolve_timer;
    presolve_timer.initialisePresolveClocks(clock);
    presolve_clocks_ = clock;
  }

  // Allow logging if option is set and LP presolve is being used
  allow_logging_ = options_->presolve_rule_logging && lp_presolve;
  // NB logging_on_ is also used to determine whether logging has
  // started to prevent double-accounting. Specifically,
  //
  // * If allow_logging_ is false, then startPresolveRuleLog is never called
  //
  // * If allow_logging_ is true, then startPresolveRuleLog will be
  //   called at some point, and within startPresolveRuleLog,
  //   allow_logging_ is set to false. Hence, it's never called again -
  //   which would lead to double-accounting. When stopPresolveRuleLog
  //   is called - and this is done by saving the value of allow_logging_
  //   before a posible call to startPresolveRuleLog - allow_logging_
  //   is set to true
  //
  // Since all reductions are logged, logging_on_ must be false in all
  // calls to markRowDeleted and markColDeleted

  logging_on_ = allow_logging_;
  called_analyse_presolve_rule_log_ = false;
  log_rule_type_ = kPresolveRuleIllegal;
  resetNumDeleted();
  presolve_log_.clear();
  original_num_col_ = model->num_col_;
  original_num_row_ = model->num_row_;
}

void HighsPresolveLog::clear() {
  this->rule.resize(kPresolveRuleCount);
  for (HighsInt rule_type = 0; rule_type < kPresolveRuleCount; rule_type++) {
    this->rule[rule_type].call = 0;
    this->rule[rule_type].col_removed = 0;
    this->rule[rule_type].row_removed = 0;
  }
}

void HPresolveAnalysis::resetNumDeleted() {
  num_deleted_rows0_ = 0;
  num_deleted_cols0_ = 0;
}

void HPresolveAnalysis::startPresolveRuleLog(const HighsInt rule_type) {
  const bool debug_print = false;
  assert(logging_on_);
  assert(rule_type >= kPresolveRuleMin && rule_type <= kPresolveRuleMax);
  // Prevent any future calls to "start" until logging is on again
  logging_on_ = false;
  const HighsInt check_rule = kPresolveRuleIllegal;
  if (debug_print)
    printf("   startPresolveRuleLog [%6d, %6d] for (%2d) %s\n",
           int(*numDeletedRows), int(*numDeletedCols), int(rule_type),
           utilPresolveRuleTypeToString(rule_type).c_str());
  if (rule_type == check_rule) {
    printf(">> startPresolveRuleLog [%6d, %6d] for (%2d) %s\n", int(check_rule),
           int(*numDeletedRows), int(*numDeletedCols),
           utilPresolveRuleTypeToString(check_rule).c_str());
  }
  presolve_log_.rule[rule_type].call++;
  // Check that stop has been called since the last start
  assert(log_rule_type_ == kPresolveRuleIllegal);
  log_rule_type_ = rule_type;
  // Check that no un-logged reductions have been performed
  if (num_deleted_rows0_ != *numDeletedRows ||
      num_deleted_cols0_ != *numDeletedCols) {
    printf(
        "ERROR: Model %s: "
        "%d = num_deleted_rows0_ != *numDeletedRows = %d ||"
        "%d = num_deleted_cols0_ != *numDeletedCols = %d\n",
        model->model_name_.c_str(), int(num_deleted_rows0_),
        int(*numDeletedRows), int(num_deleted_cols0_), int(*numDeletedCols));
    fflush(stdout);
  }
  assert(num_deleted_rows0_ == *numDeletedRows);
  assert(num_deleted_cols0_ == *numDeletedCols);
  num_deleted_rows0_ = *numDeletedRows;
  num_deleted_cols0_ = *numDeletedCols;
}

void HPresolveAnalysis::stopPresolveRuleLog(const HighsInt rule_type) {
  if (rule_type < 0) {
    // Calling from HPRESOLVE_CHECKED_CALL, so return if no logging is
    // currently being performed (indicated by log_rule_type_ ==
    // kPresolveRuleIllegal)
    if (log_rule_type_ == kPresolveRuleIllegal) return;
    assert(!logging_on_);
    // The trick of setting logging_on_ to be true before calling
    // stopPresolveRuleLog (so that start/stop calls are paired in a
    // method) won't have been done, so do it here
    logging_on_ = true;
  }
  assert(logging_on_);
  const int use_rule_type = rule_type >= 0 ? rule_type : log_rule_type_;
  const bool debug_print = false;
  if (debug_print)
    printf("    stopPresolveRuleLog [%6d, %6d] for (%2d) %s\n",
           int(*numDeletedRows), int(*numDeletedCols), int(use_rule_type),
           utilPresolveRuleTypeToString(use_rule_type).c_str());
  const int check_rule = kPresolveRuleIllegal;
  if (use_rule_type == check_rule) {
    printf(">>  stopPresolveRuleLog [%6d, %6d] for (%2d) %s\n", int(check_rule),
           int(*numDeletedRows), int(*numDeletedCols),
           utilPresolveRuleTypeToString(check_rule).c_str());
  }
  const HighsInt num_removed_row = *numDeletedRows - num_deleted_rows0_;
  const HighsInt num_removed_col = *numDeletedCols - num_deleted_cols0_;
  assert(num_removed_row >= 0);
  assert(num_removed_col >= 0);
  presolve_log_.rule[use_rule_type].col_removed += num_removed_col;
  presolve_log_.rule[use_rule_type].row_removed += num_removed_row;

  // Set the rule type to be illegal to indicate that stop has been
  // called, and update the record of num_deleted_rows/cols
  log_rule_type_ = kPresolveRuleIllegal;
  num_deleted_rows0_ = *numDeletedRows;
  num_deleted_cols0_ = *numDeletedCols;

  const bool report = false;
  if (report)
    printf("%-25s Call %9d: (%3d, %3d) (%3d, %3d)\n",
           utilPresolveRuleTypeToString(use_rule_type).c_str(),
           (int)presolve_log_.rule[use_rule_type].call, (int)num_removed_col,
           (int)num_removed_row,
           (int)presolve_log_.rule[use_rule_type].col_removed,
           (int)presolve_log_.rule[use_rule_type].row_removed);
  const int check_num_deleted_rows0_ = -212;
  const int check_num_deleted_cols0_ = -637;
  if (num_deleted_rows0_ == check_num_deleted_rows0_ &&
      num_deleted_cols0_ == check_num_deleted_cols0_) {
    printf("num_deleted (%d, %d)\n", int(num_deleted_rows0_),
           int(num_deleted_cols0_));
  }
}

bool HPresolveAnalysis::analysePresolveRuleLog(const bool report) {
  if (!allow_logging_) return true;
  if (called_analyse_presolve_rule_log_) return true;
  if (report) called_analyse_presolve_rule_log_ = true;
  const HighsLogOptions& log_options = options->log_options;
  HighsInt sum_removed_row = 0;
  HighsInt sum_removed_col = 0;
  for (HighsInt rule_type = kPresolveRuleMin; rule_type < kPresolveRuleCount;
       rule_type++) {
    sum_removed_row += presolve_log_.rule[rule_type].row_removed;
    sum_removed_col += presolve_log_.rule[rule_type].col_removed;
  }
  if (report && sum_removed_row + sum_removed_col) {
    const std::string rule =
        "-------------------------------------------------------";
    highsLogUser(log_options, HighsLogType::kInfo, "%s\n", rule.c_str());
    highsLogUser(log_options, HighsLogType::kInfo,
                 "%-25s      Rows      Cols     Calls\n",
                 "Presolve rule removed");
    highsLogUser(log_options, HighsLogType::kInfo, "%s\n", rule.c_str());
    for (HighsInt k = kPresolveRuleMin; k < kPresolveRuleCount; k++) {
      HighsInt rule_type = k;
      // Hack so that initial logging is of initial sweep
      if (k == 0) {
        rule_type = kPresolveRuleInitialSweep;
      } else if (k <= kPresolveRuleInitialSweep) {
        rule_type = k - 1;
      } else {
        assert(k > kPresolveRuleInitialSweep);
        rule_type = k;
      }
      if (presolve_log_.rule[rule_type].call ||
          presolve_log_.rule[rule_type].row_removed ||
          presolve_log_.rule[rule_type].col_removed)
        highsLogUser(log_options, HighsLogType::kInfo, "%-25s %9d %9d %9d\n",
                     utilPresolveRuleTypeToString(rule_type).c_str(),
                     (int)presolve_log_.rule[rule_type].row_removed,
                     (int)presolve_log_.rule[rule_type].col_removed,
                     (int)presolve_log_.rule[rule_type].call);
    }
    highsLogUser(log_options, HighsLogType::kInfo, "%s\n", rule.c_str());
    highsLogUser(log_options, HighsLogType::kInfo, "%-25s %9d %9d\n",
                 "Total reductions", (int)sum_removed_row,
                 (int)sum_removed_col);
    highsLogUser(log_options, HighsLogType::kInfo, "%s\n", rule.c_str());
    highsLogUser(log_options, HighsLogType::kInfo, "%-25s %9d %9d\n",
                 "Original  model", (int)original_num_row_,
                 (int)original_num_col_);
    highsLogUser(log_options, HighsLogType::kInfo, "%-25s %9d %9d\n",
                 "Presolved model", (int)(original_num_row_ - sum_removed_row),
                 (int)(original_num_col_ - sum_removed_col));
    highsLogUser(log_options, HighsLogType::kInfo, "%s\n", rule.c_str());
  }
  if (original_num_row_ == model->num_row_ &&
      original_num_col_ == model->num_col_) {
    if (sum_removed_row != *numDeletedRows) {
      highsLogDev(log_options, HighsLogType::kError,
                  "%d = sum_removed_row != numDeletedRows = %d\n",
                  (int)sum_removed_row, (int)*numDeletedRows);
      fflush(stdout);
      assert(sum_removed_row == *numDeletedRows);
      return false;
    }
    if (sum_removed_col != *numDeletedCols) {
      highsLogDev(log_options, HighsLogType::kError,
                  "%d = sum_removed_col != numDeletedCols = %d\n",
                  (int)sum_removed_col, (int)*numDeletedCols);
      fflush(stdout);
      assert(sum_removed_col == *numDeletedCols);
      return false;
    }
  }
  return true;
}

void HPresolveAnalysis::presolveTimerStart(
    const HighsInt presolve_clock) const {
  if (!analyse_presolve_time_) return;
  HighsInt highs_timer_clock = presolve_clocks_.clock_[presolve_clock];
  presolve_clocks_.timer_pointer_->start(highs_timer_clock);
}

void HPresolveAnalysis::presolveTimerStop(const HighsInt presolve_clock) const {
  if (!analyse_presolve_time_) return;
  HighsInt highs_timer_clock = presolve_clocks_.clock_[presolve_clock];
  presolve_clocks_.timer_pointer_->stop(highs_timer_clock);
}

void HPresolveAnalysis::reportPresolveTimer() {
  if (!analyse_presolve_time_) return;
  PresolveTimer presolve_timer;
  presolve_timer.reportPresolveCoreClock(model->model_name_, presolve_clocks_);
  presolve_timer.reportPresolveInitialColPresolveClock(model->model_name_,
                                                       presolve_clocks_);
  presolve_timer.reportPresolveSingletonColPresolveClock(model->model_name_,
                                                         presolve_clocks_);
}
