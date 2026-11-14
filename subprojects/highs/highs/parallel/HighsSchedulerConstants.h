/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_SCHEDULER_CONSTANTS_H_
#define HIGHS_SCHEDULER_CONSTANTS_H_

struct HighsSchedulerConstants {
  static constexpr int kNumTryFac = 16;
  static constexpr int kMicroSecsBeforeSleep = 5000;
  static constexpr int kMicroSecsBeforeGlobalSync = 1000;
  static constexpr int kNumSpinTries = 10;
  static constexpr size_t kTaskArraySize = 8192;
  static constexpr size_t kMaxTaskSize = 64;
  static constexpr size_t kCacheLineSize = 64;
};

#endif
