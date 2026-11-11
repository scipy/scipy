/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HIGHS_RACE_TIMER_H_
#define HIGHS_RACE_TIMER_H_

#include <atomic>
#include <limits>

template <typename T>
class HighsRaceTimer {
  std::atomic<T> limit;

 public:
  HighsRaceTimer() : limit(std::numeric_limits<T>::max()) {}

  void decreaseLimit(T newLimit) {
    T current = limit.load(std::memory_order_relaxed);

    while (current > newLimit) {
      if (limit.compare_exchange_weak(current, newLimit,
                                      std::memory_order_relaxed,
                                      std::memory_order_relaxed))
        break;
    }
  }

  bool limitReached(const T currentTime) const {
    return currentTime > limit.load(std::memory_order_relaxed);
  }
};

#endif
