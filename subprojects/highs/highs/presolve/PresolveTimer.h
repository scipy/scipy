/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/PresolveTimer.h
 * @brief Indices of presolve iClocks
 */
#ifndef PRESOLVE_PRESOLVETIMER_H_
#define PRESOLVE_PRESOLVETIMER_H_

// Clocks for profiling presolve
enum iClockPresolve {
  kPresolveClockTotal = 0,
  kPresolveClockPresolve,
  kPresolveClockInitialSweep,
  kPresolveClockSetupResize,
  kPresolveClockSetupToCsc,
  kPresolveClockSetupSubstitutionOpportunities,
  kPresolveClockInitial,
  kPresolveClockInitialRow,
  kPresolveClockInitialCol,
  kPresolveClockInitialColIsFixed,
  kPresolveClockInitialColIsEmpty,
  kPresolveClockInitialColIsSingleton,
  kPresolveClockInitialColDominated,
  kPresolveClockInitialColImpliedInteger,
  kPresolveClockInitialColDualFixing,
  kPresolveClockInitialColSingletonStuffing,
  kPresolveClockSingletonColSingletonRow,
  kPresolveClockSingletonColDominated,
  kPresolveClockSingletonColDualFixing,
  kPresolveClockSingletonColStuffing,
  kPresolveClockSingletonColImpliedBounds,
  kPresolveClockSingletonColRowDualImpliedBounds,
  kPresolveClockSingletonColDualImpliedFree,
  kPresolveClockFastLoop,
  kPresolveClockFastLoopRowSingletons,
  kPresolveClockFastLoopColSingletons,
  kPresolveClockFastLoopDoubletonEquations,
  kPresolveClockFastLoopChangedRows,
  kPresolveClockFastLoopChangedCols,
  kPresolveClockAggregator,
  kPresolveClockSparsify,
  kPresolveClockParallelRowsAndCols,
  kPresolveClockDependentEquations,
  kPresolveClockDependentFreeCol,
  kPresolveClockShrinkProblem,
  //  kPresolveClock@,
  kNumPresolveClock  //!< Number of PRESOLVE clocks
};

static const double kPresolveClockTolerancePercentReport = 0.1;

class PresolveTimer {
 public:
  void initialisePresolveClocks(HighsTimerClock& presolve_timer_clock) {
    HighsTimer* timer_pointer = presolve_timer_clock.timer_pointer_;
    std::vector<HighsInt>& clock = presolve_timer_clock.clock_;

    clock.resize(kNumPresolveClock);
    clock[kPresolveClockTotal] = 0;
    clock[kPresolveClockPresolve] = timer_pointer->clock_def("Presolve");
    clock[kPresolveClockInitialSweep] =
        timer_pointer->clock_def("Initial sweep");
    clock[kPresolveClockSetupResize] =
        timer_pointer->clock_def("Setup: resize");
    clock[kPresolveClockSetupToCsc] = timer_pointer->clock_def("Setup: to CSC");
    clock[kPresolveClockSetupSubstitutionOpportunities] =
        timer_pointer->clock_def("Setup substitution opportunities");
    clock[kPresolveClockInitial] = timer_pointer->clock_def("Initial");
    clock[kPresolveClockInitialRow] = timer_pointer->clock_def("Initial row");
    clock[kPresolveClockInitialCol] = timer_pointer->clock_def("Initial col");
    clock[kPresolveClockInitialColIsFixed] =
        timer_pointer->clock_def("I-col: is fixed");
    clock[kPresolveClockInitialColIsEmpty] =
        timer_pointer->clock_def("I-col: is empty");
    clock[kPresolveClockInitialColIsSingleton] =
        timer_pointer->clock_def("I-col: is singleton");
    clock[kPresolveClockInitialColDominated] =
        timer_pointer->clock_def("I-col: dominated");
    clock[kPresolveClockInitialColImpliedInteger] =
        timer_pointer->clock_def("I-col: implied integer");
    clock[kPresolveClockInitialColDualFixing] =
        timer_pointer->clock_def("I-col: dual fixing");
    clock[kPresolveClockInitialColSingletonStuffing] =
        timer_pointer->clock_def("I-col: singleton stuffing");
    clock[kPresolveClockSingletonColSingletonRow] =
        timer_pointer->clock_def("S-col: singleton row");
    clock[kPresolveClockSingletonColDominated] =
        timer_pointer->clock_def("S-col: dominated");
    clock[kPresolveClockSingletonColDualFixing] =
        timer_pointer->clock_def("S-col: dual fixing");
    clock[kPresolveClockSingletonColStuffing] =
        timer_pointer->clock_def("S-col: singleton stuffing");
    clock[kPresolveClockSingletonColImpliedBounds] =
        timer_pointer->clock_def("S-col: impl bounds");
    clock[kPresolveClockSingletonColRowDualImpliedBounds] =
        timer_pointer->clock_def("S-col: row dual impl bounds");
    clock[kPresolveClockSingletonColDualImpliedFree] =
        timer_pointer->clock_def("S-col: dual impl free");
    clock[kPresolveClockFastLoop] = timer_pointer->clock_def("Fast loop");
    clock[kPresolveClockFastLoopRowSingletons] =
        timer_pointer->clock_def("Fast loop: row singletons");
    clock[kPresolveClockFastLoopColSingletons] =
        timer_pointer->clock_def("Fast loop: col singletons");
    clock[kPresolveClockFastLoopDoubletonEquations] =
        timer_pointer->clock_def("Fast loop: doubleton equations");
    clock[kPresolveClockFastLoopChangedRows] =
        timer_pointer->clock_def("Fast loop: changed rows");
    clock[kPresolveClockFastLoopChangedCols] =
        timer_pointer->clock_def("Fast loop: changed cols");
    clock[kPresolveClockAggregator] = timer_pointer->clock_def("Aggregator");
    clock[kPresolveClockSparsify] = timer_pointer->clock_def("Sparsify");
    clock[kPresolveClockParallelRowsAndCols] =
        timer_pointer->clock_def("Parallel rows and cols");
    clock[kPresolveClockDependentEquations] =
        timer_pointer->clock_def("Dependent equations");
    clock[kPresolveClockDependentFreeCol] =
        timer_pointer->clock_def("Dependent free columns");
    clock[kPresolveClockShrinkProblem] =
        timer_pointer->clock_def("Shrink problem");
    //    clock[kPresolveClock@] = timer_pointer->clock_def("@");
  };

  bool reportPresolveClockList(
      const char* grepStamp, const std::vector<HighsInt>& presolve_clock_list,
      const HighsTimerClock& presolve_timer_clock,
      const HighsInt kPresolveClockIdeal = kPresolveClockPresolve,
      const double tolerance_percent_report_ = -1) {
    HighsTimer* timer_pointer = presolve_timer_clock.timer_pointer_;
    if (!timer_pointer->printf_flag) return false;
    const std::vector<HighsInt>& clock = presolve_timer_clock.clock_;
    HighsInt presolve_clock_list_size = presolve_clock_list.size();
    std::vector<HighsInt> clockList;
    clockList.resize(presolve_clock_list_size);
    for (HighsInt en = 0; en < presolve_clock_list_size; en++) {
      clockList[en] = clock[presolve_clock_list[en]];
    }
    const double ideal_sum_time =
        timer_pointer->clock_time[clock[kPresolveClockIdeal]];
    const double tolerance_percent_report =
        tolerance_percent_report_ >= 0 ? tolerance_percent_report_ : 1e-8;
    return timer_pointer->reportOnTolerance(
        grepStamp, clockList, ideal_sum_time, tolerance_percent_report);
  };

  void csvPresolveClockList(const std::string& grep_query,
                            const std::string& model_name,
                            const std::vector<HighsInt>& presolve_clock_list,
                            const HighsTimerClock& presolve_timer_clock,
                            const HighsInt kPresolveClockIdeal,
                            const bool header, const bool end_line) {
    HighsTimer* timer_pointer = presolve_timer_clock.timer_pointer_;
    if (!timer_pointer->printf_flag) return;
    const std::vector<HighsInt>& clock = presolve_timer_clock.clock_;
    const double ideal_sum_time =
        timer_pointer->clock_time[clock[kPresolveClockIdeal]];
    if (ideal_sum_time < 1e-2) return;
    const HighsInt num_clock = presolve_clock_list.size();
    if (header) {
      printf("grep_%s,model,ideal", grep_query.c_str());
      for (HighsInt iX = 0; iX < num_clock; iX++) {
        HighsInt iclock = clock[presolve_clock_list[iX]];
        printf(",%s", timer_pointer->clock_names[iclock].c_str());
      }
      printf(",Unaccounted");
      if (end_line) printf("\n");
      return;
    }
    double sum_time = 0;
    printf("grep_%s,%s,%11.4g", grep_query.c_str(), model_name.c_str(),
           ideal_sum_time);
    for (HighsInt iX = 0; iX < num_clock; iX++) {
      HighsInt iclock = clock[presolve_clock_list[iX]];
      double time = timer_pointer->read(iclock);
      sum_time += time;
      printf(",%11.4g", time);
    }
    printf(",%11.4g", ideal_sum_time - sum_time);
    if (end_line) printf("\n");
  }

  void reportPresolveCoreClock(const std::string& model_name,
                               const HighsTimerClock& presolve_timer_clock) {
    const std::vector<HighsInt> presolve_clock_list{
        kPresolveClockInitialSweep, kPresolveClockSetupResize,
        kPresolveClockSetupToCsc, kPresolveClockSetupSubstitutionOpportunities,
        //      kPresolveClockInitial,
        kPresolveClockInitialRow, kPresolveClockInitialCol,
        //	kPresolveClockFastLoop,
        kPresolveClockFastLoopRowSingletons,
        kPresolveClockFastLoopColSingletons,
        kPresolveClockFastLoopDoubletonEquations,
        kPresolveClockFastLoopChangedRows, kPresolveClockFastLoopChangedCols,
        kPresolveClockAggregator, kPresolveClockSparsify,
        kPresolveClockParallelRowsAndCols, kPresolveClockDependentEquations,
        kPresolveClockDependentFreeCol, kPresolveClockShrinkProblem
        //	kPresolveClock@
    };
    reportPresolveClockList("PresolveCore_", presolve_clock_list,
                            presolve_timer_clock, kPresolveClockPresolve, 0.1);
    const bool csv_output = false;
    if (csv_output) {
      csvPresolveClockList("GrepPresolveCore_", model_name, presolve_clock_list,
                           presolve_timer_clock, kPresolveClockPresolve, true,
                           true);
      csvPresolveClockList("GrepPresolveCore_", model_name, presolve_clock_list,
                           presolve_timer_clock, kPresolveClockPresolve, false,
                           true);
    }
  };

  void reportPresolveInitialColPresolveClock(
      const std::string& model_name,
      const HighsTimerClock& presolve_timer_clock) {
    const std::vector<HighsInt> presolve_clock_list{
        kPresolveClockInitialColIsFixed,
        kPresolveClockInitialColIsEmpty,
        kPresolveClockInitialColIsSingleton,
        kPresolveClockInitialColDominated,
        kPresolveClockInitialColImpliedInteger,
        kPresolveClockInitialColDualFixing,
        kPresolveClockInitialColSingletonStuffing};
    reportPresolveClockList("PresolveInitialCol_", presolve_clock_list,
                            presolve_timer_clock, kPresolveClockInitialCol,
                            0.1);
  };

  void reportPresolveSingletonColPresolveClock(
      const std::string& model_name,
      const HighsTimerClock& presolve_timer_clock) {
    const std::vector<HighsInt> presolve_clock_list{
        kPresolveClockSingletonColSingletonRow,
        kPresolveClockSingletonColDominated,
        kPresolveClockSingletonColDualFixing,
        kPresolveClockSingletonColStuffing,
        kPresolveClockSingletonColImpliedBounds,
        kPresolveClockSingletonColRowDualImpliedBounds,
        kPresolveClockSingletonColDualImpliedFree};
    reportPresolveClockList("PresolveSingletonCol_", presolve_clock_list,
                            presolve_timer_clock,
                            kPresolveClockInitialColIsSingleton, 0.1);
  };
};

#endif /* PRESOLVE_PRESOLVETIMER_H_ */
