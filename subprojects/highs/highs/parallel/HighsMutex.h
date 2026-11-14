/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HIGHS_MUTEX_H_
#define HIGHS_MUTEX_H_

#include "parallel/HighsSpinMutex.h"
#include "parallel/HighsTaskExecutor.h"

class HighsMutex {
  std::atomic<unsigned int> state{0u};

 public:
  bool try_lock() {
    unsigned int uncontendedState = 0;
    return state.compare_exchange_weak(uncontendedState, 1,
                                       std::memory_order_acquire,
                                       std::memory_order_relaxed);
  }

  void lock() {
    // First try to acquire the lock directly to have a fast path for
    // uncontended access
    if (try_lock()) return;

    // Now spin a few times to check if the lock becomes available
    for (int i = 0; i < HighsSchedulerConstants::kNumSpinTries; ++i) {
      if (state.load(std::memory_order_relaxed) == 0) {
        if (try_lock()) return;
      }

      HighsSpinMutex::yieldProcessor();
    }

    // Lock is still not available, so now start a loop where we yield
    // to the scheduler after each time we observe the lock as unavailable
    // so that we can see if there are other tasks and then check
    // the lock again. We do so for some microseconds which we measure starting
    // from now.
    HighsSplitDeque* thisDeque = HighsTaskExecutor::getThisWorkerDeque();
    auto tStart = std::chrono::high_resolution_clock::now();

    while (true) {
      int numTries = HighsSchedulerConstants::kNumSpinTries;

      for (int i = 0; i < numTries; ++i) {
        if (state.load(std::memory_order_relaxed) == 0) {
          if (try_lock()) return;
        }
        thisDeque->yield();
      }

      auto numMicroSecs =
          std::chrono::duration_cast<std::chrono::microseconds>(
              std::chrono::high_resolution_clock::now() - tStart)
              .count();

      if (numMicroSecs < HighsSchedulerConstants::kMicroSecsBeforeGlobalSync)
        numTries *= 2;
      else
        break;
    }

    // The lock is still not available, now we will try to set ourselves as the
    // next worker to acquire the lock and start wait until we are
    // notified by the current worker holding the lock.
    unsigned int ownerId = (thisDeque->getOwnerId() + 1) << 1;
    unsigned int s = state.load(std::memory_order_relaxed);
    while (true) {
      // as long as we observe that the lock is available we try to lock it, so
      // that we are guaranteed to have a locked state stored within s.
      while (s == 0)
        if (state.compare_exchange_weak(s, 1, std::memory_order_acquire,
                                        std::memory_order_relaxed))
          return;

      if (s & 1u) {
        if (state.compare_exchange_weak(s, ownerId | 1,
                                        std::memory_order_release,
                                        std::memory_order_relaxed))
          break;
      } else {
        // if we observe that the semaphore is unlocked but has its state not
        // equal to zero it means another worker has just been notified to take
        // the lock, but not yet taken it. Before we can set ourselves as the
        // next worker we need to wait for its state to be in locked state.
        HighsSpinMutex::yieldProcessor();
        s = state.load(std::memory_order_relaxed);
      }
    }

    // once we are here we have successfully set ourselves as the next worker
    // and once the previous worker released the lock we will be notified. Hence
    // we wait until a notify and then are free to acquire the lock
    // unconditionally.
    thisDeque->wait();

    // set the state to have its locked bit set and store the previous next
    // worker as the new next worker. As we only get here when the compare
    // exchange in the loop above succeeds with s having its locked bit set we
    // can simply store back s.
    state.store(s, std::memory_order_relaxed);
  }

  void unlock() {
    unsigned int prevState = state.fetch_add(
        std::numeric_limits<unsigned int>::max(), std::memory_order_release);

    if (prevState != 1) {
      unsigned int notifyWorkerId = (prevState >> 1) - 1;
      HighsTaskExecutor::getThisWorkerDeque()
          ->getWorkerById(notifyWorkerId)
          ->notify();
    }
  }
};

#endif
