/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/HPresolveAnalysis.h
 * @brief
 */
#ifndef PRESOLVE_HPRESOLVEANALYSIS_H_
#define PRESOLVE_HPRESOLVEANALYSIS_H_

#include "util/HighsTimer.h"

class HPresolveAnalysis {
 public:
  HPresolveAnalysis() : timer_(nullptr), analyse_presolve_time_(false) {}

  HighsTimer* timer_;
  const HighsLp* model;
  const HighsOptions* options;
  const HighsInt* numDeletedRows;
  const HighsInt* numDeletedCols;

  // store original problem sizes for reference
  HighsInt original_num_col_;
  HighsInt original_num_row_;

 public:
  bool allow_logging_;
  bool logging_on_;
  bool called_analyse_presolve_rule_log_;

  int log_rule_type_;
  HighsInt num_deleted_rows0_;
  HighsInt num_deleted_cols0_;
  HighsPresolveLog presolve_log_;

  HighsTimerClock presolve_clocks_;
  bool analyse_presolve_time_;

  void setup(const HighsLp* model_, const HighsOptions* options_,
             const HighsInt& numDeletedRows_, const HighsInt& numDeletedCols_,
             HighsTimer* timer);
  void setupPresolveTime(const HighsOptions& options);
  void resetNumDeleted();

  std::string presolveReductionTypeToString(const HighsInt reduction_type);
  void startPresolveRuleLog(const HighsInt rule_type);
  void stopPresolveRuleLog(const HighsInt rule_type = -1);
  bool analysePresolveRuleLog(const bool report = false);
  void presolveTimerStart(const HighsInt presolve_clock = 0) const;
  void presolveTimerStop(const HighsInt presolve_clock = 0) const;
  void reportPresolveTimer();
  friend class HPresolve;

 private:
  // Clock interface for HighsClockScope
  void start(const HighsInt presolve_clock) const {
    presolveTimerStart(presolve_clock);
  }
  void stop(const HighsInt presolve_clock) const {
    presolveTimerStop(presolve_clock);
  }
  friend class HighsClockScope<const HPresolveAnalysis>;
};

// Starts a presolve timer and stops it when the scope is left, on every exit
using PresolveTimerScope = HighsClockScope<const HPresolveAnalysis>;

#endif /* PRESOLVE_HPRESOLVEANALYSIS_H_ */
