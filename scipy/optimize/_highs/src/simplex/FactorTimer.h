/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/FactorTimer.h
 * @brief Indices of factor iClocks
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_FACTORTIMER_H_
#define SIMPLEX_FACTORTIMER_H_

#include "HConfig.h"
#include "lp_data/HighsAnalysis.h"
#include "util/HighsTimer.h"

// Clocks for profiling the dual simplex solver
enum iClockFactor {
  FactorInvert = 0,        //!< INVERT
  FactorInvertSimple,      //!< INVERT simple
  FactorInvertKernel,      //!< INVERT kernel
  FactorInvertDeficient,   //!< INVERT deficient
  FactorInvertFinish,      //!< INVERT finish
  FactorFtran,             //!< FTRAN
  FactorFtranLower,        //!< FTRAN Lower part
  FactorFtranLowerAPF,     //!< FTRAN Lower part APF
  FactorFtranLowerSps,     //!< FTRAN Lower part sparse
  FactorFtranLowerHyper,   //!< FTRAN Lower part hyper-sparse
  FactorFtranUpper,        //!< FTRAN Upper part
  FactorFtranUpperFT,      //!< FTRAN Upper part FT
  FactorFtranUpperMPF,     //!< FTRAN Upper part MPF
  FactorFtranUpperSps0,    //!< FTRAN Upper part sparse
  FactorFtranUpperSps1,    //!< FTRAN Upper part sparse
  FactorFtranUpperSps2,    //!< FTRAN Upper part sparse
  FactorFtranUpperHyper0,  //!< FTRAN Upper part hyper-sparse
  FactorFtranUpperHyper1,  //!< FTRAN Upper part hyper-sparse
  FactorFtranUpperHyper2,  //!< FTRAN Upper part hyper-sparse
  FactorFtranUpperHyper3,  //!< FTRAN Upper part hyper-sparse
  FactorFtranUpperHyper4,  //!< FTRAN Upper part hyper-sparse
  FactorFtranUpperHyper5,  //!< FTRAN Upper part hyper-sparse
  FactorFtranUpperPF,      //!< FTRAN Upper part PF
  FactorBtran,             //!< BTRAN
  FactorBtranLower,        //!< BTRAN Lower part
  FactorBtranLowerSps,     //!< BTRAN Lower part sparse
  FactorBtranLowerHyper,   //!< BTRAN Lower part hyper-sparse
  FactorBtranLowerAPF,     //!< BTRAN Lower part APF
  FactorBtranUpper,        //!< BTRAN Upper part
  FactorBtranUpperPF,      //!< BTRAN Upper part PF
  FactorBtranUpperSps,     //!< BTRAN Upper part sparse
  FactorBtranUpperHyper,   //!< BTRAN Upper part hyper-sparse
  FactorBtranUpperFT,      //!< BTRAN Upper part FT
  FactorBtranUpperMPF,     //!< BTRAN Upper part MPF
  FactorNumClock           //!< Number of factor clocks
};

class FactorTimer {
 public:
  void start(const int factor_clock,
             HighsTimerClock* factor_timer_clock_pointer) {
    if (factor_timer_clock_pointer != NULL)
      factor_timer_clock_pointer->timer_.start(
          factor_timer_clock_pointer->clock_[factor_clock]);
  };

  void stop(const int factor_clock,
            HighsTimerClock* factor_timer_clock_pointer) {
    if (factor_timer_clock_pointer != NULL)
      factor_timer_clock_pointer->timer_.stop(
          factor_timer_clock_pointer->clock_[factor_clock]);
  };

  double read(const int factor_clock,
              HighsTimerClock* factor_timer_clock_pointer) {
    double argument = 0;
    if (factor_timer_clock_pointer != NULL)
      argument = factor_timer_clock_pointer->timer_.read(
          factor_timer_clock_pointer->clock_[factor_clock]);
    return argument;
  };
#ifdef HiGHSDEV
  void initialiseFactorClocks(HighsTimerClock& factor_timer_clock) {
    HighsTimer& timer = factor_timer_clock.timer_;
    std::vector<int>& clock = factor_timer_clock.clock_;
    clock.resize(FactorNumClock);
    clock[FactorInvert] = timer.clock_def("INVERT", "INV");
    clock[FactorInvertSimple] = timer.clock_def("INVERT Simple", "IVS");
    clock[FactorInvertKernel] = timer.clock_def("INVERT Kernel", "IVK");
    clock[FactorInvertDeficient] = timer.clock_def("INVERT Deficient", "IVD");
    clock[FactorInvertFinish] = timer.clock_def("INVERT Finish", "IVF");
    clock[FactorFtran] = timer.clock_def("FTRAN", "FTR");
    clock[FactorFtranLower] = timer.clock_def("FTRAN Lower", "FTL");
    clock[FactorFtranLowerAPF] = timer.clock_def("FTRAN Lower APF", "FLA");
    clock[FactorFtranLowerSps] = timer.clock_def("FTRAN Lower Sps", "FLS");
    clock[FactorFtranLowerHyper] = timer.clock_def("FTRAN Lower Hyper", "FLH");
    clock[FactorFtranUpper] = timer.clock_def("FTRAN Upper", "FTU");
    clock[FactorFtranUpperFT] = timer.clock_def("FTRAN Upper FT", "FUF");
    clock[FactorFtranUpperMPF] = timer.clock_def("FTRAN Upper MPF", "FUM");
    clock[FactorFtranUpperSps0] = timer.clock_def("FTRAN Upper Sps0", "FUS");
    clock[FactorFtranUpperSps1] = timer.clock_def("FTRAN Upper Sps1", "FUS");
    clock[FactorFtranUpperSps2] = timer.clock_def("FTRAN Upper Sps2", "FUS");
    clock[FactorFtranUpperHyper0] =
        timer.clock_def("FTRAN Upper Hyper0", "FUH");
    clock[FactorFtranUpperHyper1] =
        timer.clock_def("FTRAN Upper Hyper1", "FUH");
    clock[FactorFtranUpperHyper2] =
        timer.clock_def("FTRAN Upper Hyper2", "FUH");
    clock[FactorFtranUpperHyper3] =
        timer.clock_def("FTRAN Upper Hyper3", "FUH");
    clock[FactorFtranUpperHyper4] =
        timer.clock_def("FTRAN Upper Hyper4", "FUH");
    clock[FactorFtranUpperHyper5] =
        timer.clock_def("FTRAN Upper Hyper5", "FUH");
    clock[FactorFtranUpperPF] = timer.clock_def("FTRAN Upper PF", "FUP");
    clock[FactorBtran] = timer.clock_def("BTRAN", "BTR");
    clock[FactorBtranLower] = timer.clock_def("BTRAN Lower", "BTL");
    clock[FactorBtranLowerSps] = timer.clock_def("BTRAN Lower Sps", "BLS");
    clock[FactorBtranLowerHyper] = timer.clock_def("BTRAN Lower Hyper", "BLH");
    clock[FactorBtranLowerAPF] = timer.clock_def("BTRAN Lower APF", "BLA");
    clock[FactorBtranUpper] = timer.clock_def("BTRAN Upper", "BTU");
    clock[FactorBtranUpperPF] = timer.clock_def("BTRAN Upper PF", "BUP");
    clock[FactorBtranUpperSps] = timer.clock_def("BTRAN Upper Sps", "BUS");
    clock[FactorBtranUpperHyper] = timer.clock_def("BTRAN Upper Hyper", "BUH");
    clock[FactorBtranUpperFT] = timer.clock_def("BTRAN Upper FT", "BUF");
    clock[FactorBtranUpperMPF] = timer.clock_def("BTRAN Upper MPS", "BUM");
  };

  void reportFactorClockList(const char* grepStamp,
                             HighsTimerClock& factor_timer_clock,
                             std::vector<int> factor_clock_list) {
    HighsTimer& timer = factor_timer_clock.timer_;
    std::vector<int>& clock = factor_timer_clock.clock_;
    int factor_clock_list_size = factor_clock_list.size();
    std::vector<int> clockList;
    clockList.resize(factor_clock_list_size);
    for (int en = 0; en < factor_clock_list_size; en++) {
      clockList[en] = clock[factor_clock_list[en]];
    }
    double ideal_sum_time = 0;
    ideal_sum_time += timer.read(clock[FactorInvert]);
    ideal_sum_time += timer.read(clock[FactorFtran]);
    ideal_sum_time += timer.read(clock[FactorBtran]);
    timer.report_tl(grepStamp, clockList, ideal_sum_time, 1e-8);
  };

  void reportFactorLevel0Clock(HighsTimerClock& factor_timer_clock) {
    std::vector<int> factor_clock_list{FactorInvert, FactorFtran, FactorBtran};
    reportFactorClockList("FactorLevel0", factor_timer_clock,
                          factor_clock_list);
  };

  void reportFactorLevel1Clock(HighsTimerClock& factor_timer_clock) {
    std::vector<int> factor_clock_list{
        FactorInvertSimple, FactorInvertKernel, FactorInvertDeficient,
        FactorInvertFinish, FactorFtranLower,   FactorFtranUpper,
        FactorBtranLower,   FactorBtranUpper};
    reportFactorClockList("FactorLevel1", factor_timer_clock,
                          factor_clock_list);
  };

  void reportFactorLevel2Clock(HighsTimerClock& factor_timer_clock) {
    std::vector<int> factor_clock_list{
        FactorInvertSimple,     FactorInvertKernel,     FactorInvertDeficient,
        FactorInvertFinish,     FactorFtranLowerAPF,    FactorFtranLowerSps,
        FactorFtranLowerHyper,  FactorFtranUpperFT,     FactorFtranUpperMPF,
        FactorFtranUpperSps0,   FactorFtranUpperSps1,   FactorFtranUpperSps2,
        FactorFtranUpperHyper0, FactorFtranUpperHyper1, FactorFtranUpperHyper2,
        FactorFtranUpperHyper3, FactorFtranUpperHyper4, FactorFtranUpperHyper5,
        FactorFtranUpperPF,     FactorBtranLowerSps,    FactorBtranLowerHyper,
        FactorBtranLowerAPF,    FactorBtranUpperPF,     FactorBtranUpperSps,
        FactorBtranUpperHyper,  FactorBtranUpperFT,     FactorBtranUpperMPF};
    reportFactorClockList("FactorLevel2", factor_timer_clock,
                          factor_clock_list);
  };

  void reportFactorClock(HighsTimerClock& factor_timer_clock) {
    reportFactorLevel0Clock(factor_timer_clock);
    reportFactorLevel1Clock(factor_timer_clock);
    reportFactorLevel2Clock(factor_timer_clock);
  }

#endif
};
#endif /* SIMPLEX_FACTORTIMER_H_ */
