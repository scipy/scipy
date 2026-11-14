/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/HiPdlpTimer.h
 * @brief Indices of HiPDLP iClocks
 */
#ifndef PDLP_HIPDLPTIMER_H_
#define PDLP_HIPDLPTIMER_H_

// Clocks for profiling HiPDLP
enum iClockHipdlp {
  kHipdlpClockTotal = 0,
  kHipdlpClockPreprocess,
  kHipdlpClockSolve,
  kHipdlpClockPostprocess,
  // Solve clocks
  kHipdlpClockIterateUpdate,
  kHipdlpClockConvergenceCheck,
  kHipdlpClockRestartCheck,
  kHipdlpClockAverageIterate,
  // IterateUpdate clocks
  kHipdlpClockMatrixMultiply,
  kHipdlpClockMatrixTransposeMultiply,
  kHipdlpClockProjectX,
  kHipdlpClockProjectY,
  kHipdlpClockStepSizeAdjustment,
  // AverageIterate clocks
  kHipdlpClockAverageIterateUpdateX,
  kHipdlpClockAverageIterateUpdateY,
  kHipdlpClockAverageIterateComputeX,
  kHipdlpClockAverageIterateComputeY,
  kHipdlpClockAverageIterateMatrixMultiply,
  kHipdlpClockAverageIterateMatrixTransposeMultiply,

  kNumHipdlpClock  //!< Number of HIPDLP clocks
};

const double tolerance_percent_report = 0.1;

class HipdlpTimer {
 public:
  void initialiseHipdlpClocks(HighsTimerClock& hipdlp_timer_clock) {
    HighsTimer* timer_pointer = hipdlp_timer_clock.timer_pointer_;
    std::vector<HighsInt>& clock = hipdlp_timer_clock.clock_;

    clock.resize(kNumHipdlpClock);
    clock[kHipdlpClockTotal] = 0;
    clock[kHipdlpClockPreprocess] =
        timer_pointer->clock_def("HiPDLP preprocess");
    clock[kHipdlpClockSolve] = timer_pointer->clock_def("HiPDLP solve");
    clock[kHipdlpClockPostprocess] =
        timer_pointer->clock_def("HiPDLP postprocess");
    clock[kHipdlpClockIterateUpdate] =
        timer_pointer->clock_def("Iterate update");
    clock[kHipdlpClockMatrixMultiply] = timer_pointer->clock_def("Ax");
    clock[kHipdlpClockMatrixTransposeMultiply] =
        timer_pointer->clock_def("Aty");
    clock[kHipdlpClockProjectX] = timer_pointer->clock_def("Project X");
    clock[kHipdlpClockProjectY] = timer_pointer->clock_def("Project Y");
    clock[kHipdlpClockStepSizeAdjustment] =
        timer_pointer->clock_def("Step size adjustment");
    clock[kHipdlpClockConvergenceCheck] =
        timer_pointer->clock_def("Convergence check");
    clock[kHipdlpClockRestartCheck] = timer_pointer->clock_def("Restart check");
    clock[kHipdlpClockAverageIterate] =
        timer_pointer->clock_def("Average iterate");
    clock[kHipdlpClockAverageIterateUpdateX] =
        timer_pointer->clock_def("Update X");
    clock[kHipdlpClockAverageIterateUpdateY] =
        timer_pointer->clock_def("Update Y");
    clock[kHipdlpClockAverageIterateComputeX] =
        timer_pointer->clock_def("Compute X");
    clock[kHipdlpClockAverageIterateComputeY] =
        timer_pointer->clock_def("Compute Y");
    clock[kHipdlpClockAverageIterateMatrixMultiply] =
        timer_pointer->clock_def("Ax  [A]");
    clock[kHipdlpClockAverageIterateMatrixTransposeMultiply] =
        timer_pointer->clock_def("Aty [A]");
  };

  bool reportHipdlpClockList(
      const char* grepStamp, const std::vector<HighsInt> hipdlp_clock_list,
      const HighsTimerClock& hipdlp_timer_clock,
      const HighsInt kHipdlpClockIdeal = kHipdlpClockTotal,
      const double tolerance_percent_report_ = -1) {
    HighsTimer* timer_pointer = hipdlp_timer_clock.timer_pointer_;
    if (!timer_pointer->printf_flag) return false;
    const std::vector<HighsInt>& clock = hipdlp_timer_clock.clock_;
    HighsInt hipdlp_clock_list_size = hipdlp_clock_list.size();
    std::vector<HighsInt> clockList;
    clockList.resize(hipdlp_clock_list_size);
    for (HighsInt en = 0; en < hipdlp_clock_list_size; en++) {
      clockList[en] = clock[hipdlp_clock_list[en]];
    }
    const double ideal_sum_time =
        timer_pointer->clock_time[clock[kHipdlpClockIdeal]];
    const double tolerance_percent_report =
        tolerance_percent_report_ >= 0 ? tolerance_percent_report_ : 1e-8;
    return timer_pointer->reportOnTolerance(
        grepStamp, clockList, ideal_sum_time, tolerance_percent_report);
  };

  void reportHipdlpCoreClock(const HighsTimerClock& hipdlp_timer_clock) {
    const std::vector<HighsInt> hipdlp_clock_list{
        kHipdlpClockPreprocess, kHipdlpClockSolve, kHipdlpClockPostprocess};
    reportHipdlpClockList("HipdlpCore___", hipdlp_clock_list,
                          hipdlp_timer_clock, kHipdlpClockTotal);
  };

  void reportHipdlpSolveClock(const HighsTimerClock& hipdlp_timer_clock) {
    const std::vector<HighsInt> hipdlp_clock_list{
        kHipdlpClockIterateUpdate, kHipdlpClockConvergenceCheck,
        kHipdlpClockRestartCheck, kHipdlpClockAverageIterate};
    reportHipdlpClockList("HipdlpSolve__", hipdlp_clock_list,
                          hipdlp_timer_clock, kHipdlpClockSolve);
  };

  void reportHipdlpIterateUpdateClock(
      const HighsTimerClock& hipdlp_timer_clock) {
    const std::vector<HighsInt> hipdlp_clock_list{
        kHipdlpClockMatrixMultiply, kHipdlpClockMatrixTransposeMultiply,
        kHipdlpClockProjectX, kHipdlpClockProjectY};
    reportHipdlpClockList("HipdlpIterUpd", hipdlp_clock_list,
                          hipdlp_timer_clock, kHipdlpClockIterateUpdate);
  };

  void reportHipdlpAverageIterateClock(
      const HighsTimerClock& hipdlp_timer_clock) {
    const std::vector<HighsInt> hipdlp_clock_list{
        kHipdlpClockAverageIterateUpdateX,
        kHipdlpClockAverageIterateUpdateY,
        kHipdlpClockAverageIterateComputeX,
        kHipdlpClockAverageIterateComputeY,
        kHipdlpClockAverageIterateMatrixMultiply,
        kHipdlpClockAverageIterateMatrixTransposeMultiply};
    reportHipdlpClockList("HipdlpAvgIter", hipdlp_clock_list,
                          hipdlp_timer_clock, kHipdlpClockAverageIterate);
  };

  void reportHipdlpMatrixMultiplyClock(
      const HighsTimerClock& hipdlp_timer_clock) {
    const std::vector<HighsInt> hipdlp_clock_list{
        kHipdlpClockMatrixMultiply, kHipdlpClockMatrixTransposeMultiply,
        kHipdlpClockAverageIterateMatrixMultiply,
        kHipdlpClockAverageIterateMatrixTransposeMultiply};
    reportHipdlpClockList("HipdlpMtxMult", hipdlp_clock_list,
                          hipdlp_timer_clock, kHipdlpClockSolve);
  };
};

#endif /* PDLP_HIPDLPTIMER_H_ */
