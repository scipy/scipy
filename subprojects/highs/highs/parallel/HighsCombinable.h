/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HIGHS_COMBINABLE_H_
#define HIGHS_COMBINABLE_H_

#include <functional>

#include "parallel/HighsCacheAlign.h"
#include "parallel/HighsTaskExecutor.h"

template <typename T, typename FConstruct_ = std::function<T(void)>>
class HighsCombinable {
  FConstruct_ construct_;
  struct PaddedData {
    alignas(highs::cache_aligned::alignment()) bool initialized_;
    T data_;
  };
  int numThreads;
  highs::cache_aligned::unique_ptr<PaddedData[]> threadCopies_;

 public:
  HighsCombinable()
      : construct_([]() { return T(); }),
        numThreads(HighsTaskExecutor::getNumWorkerThreads()),
        threadCopies_(
            highs::cache_aligned::make_unique_array<PaddedData>(numThreads)) {
    for (int i = 0; i < numThreads; ++i) threadCopies_[i].initialized_ = false;
  }

  template <typename FConstruct>
  HighsCombinable(FConstruct&& fConstruct)
      : construct_(std::forward<FConstruct>(fConstruct)),
        numThreads(HighsTaskExecutor::getNumWorkerThreads()),
        threadCopies_(
            highs::cache_aligned::make_unique_array<PaddedData>(numThreads)) {
    for (int i = 0; i < numThreads; ++i) threadCopies_[i].initialized_ = false;
  }

  HighsCombinable<T, FConstruct_>& operator=(
      HighsCombinable<T, FConstruct_>&&) = default;
  HighsCombinable(HighsCombinable<T, FConstruct_>&&) = default;

  void clear() {
    for (int i = 0; i < numThreads; ++i) {
      if (threadCopies_[i].initialized_) {
        threadCopies_[i].initialized_ = false;
        threadCopies_[i].data_.~T();
      }
    }
  }

  T& local() {
    int threadId = HighsTaskExecutor::getThisWorkerDeque()->getOwnerId();
    if (!threadCopies_[threadId].initialized_) {
      threadCopies_[threadId].initialized_ = true;
      new (&threadCopies_[threadId].data_) T(construct_());
    }

    return threadCopies_[threadId].data_;
  }

  const T& local() const {
    int threadId = HighsTaskExecutor::getThisWorkerDeque()->getOwnerId();
    if (!threadCopies_[threadId].initialized_) {
      threadCopies_[threadId].initialized_ = true;
      new (&threadCopies_[threadId].data_) T(construct_());
    }

    return threadCopies_[threadId].data_;
  }

  template <typename FCombine>
  void combine_each(FCombine&& combine) {
    for (int i = 0; i < numThreads; ++i)
      if (threadCopies_[i].initialized_) combine(threadCopies_[i].data_);
  }

  template <typename FCombine>
  T combine(FCombine&& combine) {
    T combined;
    int i;
    for (i = 0; i < numThreads; ++i) {
      if (threadCopies_[i].initialized_) {
        combined = std::move(threadCopies_[i].data_);
        break;
      }
    }

    for (++i; i < numThreads; ++i) {
      if (threadCopies_[i].initialized_) {
        combined =
            combine(std::move(combined), std::move(threadCopies_[i].data_));
        break;
      }
    }

    return combined;
  }

  ~HighsCombinable() {
    if (threadCopies_) clear();
  }
};

template <typename U, typename FConstruct>
HighsCombinable<U, FConstruct> makeHighsCombinable(FConstruct&& fconstruct) {
  return HighsCombinable<U, FConstruct>(std::forward<FConstruct>(fconstruct));
}

#endif
