/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_TASK_H_
#define HIGHS_TASK_H_

#include <atomic>
#include <cassert>
#include <cstring>
#include <type_traits>

#include "parallel/HighsSpinMutex.h"

class HighsSplitDeque;

class HighsTask {
  friend class HighsSplitDeque;

 public:
  class Interrupt {};

 private:
  static constexpr uint64_t kFinishedFlag = 1;
  static constexpr uint64_t kCancelFlag = 2;
  static constexpr uint64_t kPtrMask = ~(kFinishedFlag | kCancelFlag);
  struct Metadata {
    std::atomic<uintptr_t> stealer;
  };

  class CallableBase {
   public:
    virtual void operator()() = 0;
  };

  template <typename F>
  class Callable : public CallableBase {
    F functor;

   public:
    Callable(F&& functor) : functor(std::forward<F>(functor)) {}

    virtual void operator()() override {
      F callFunctor = std::move(functor);
      callFunctor();
    }
  };

  char taskData[HighsSchedulerConstants::kMaxTaskSize - sizeof(Metadata)];
  Metadata metadata;

  CallableBase& getCallable() {
    union {
      CallableBase* callablePtr;
      char* storagePtr;
    } u;

    u.storagePtr = this->taskData;
    return *u.callablePtr;
  }

  HighsSplitDeque* markAsFinished(HighsSplitDeque* stealer) {
    uintptr_t state =
        metadata.stealer.exchange(kFinishedFlag, std::memory_order_release);
    HighsSplitDeque* waitingOwner =
        reinterpret_cast<HighsSplitDeque*>(state & kPtrMask);
    if (waitingOwner != stealer) return waitingOwner;

    return nullptr;
  }

  /// run task as stealer and return the owner's split deque if the owner is
  /// waiting and needs to be signaled
  HighsSplitDeque* run(HighsSplitDeque* stealer) {
    uintptr_t state = metadata.stealer.fetch_or(
        reinterpret_cast<uintptr_t>(stealer), std::memory_order_acquire);
    if (state == 0) getCallable()();

    return markAsFinished(stealer);
  }

 public:
  /// initialize the task with given callable type. Task is considered
  /// unfinished after setTaskData
  template <typename F>
  void setTaskData(F&& f) {
    static_assert(sizeof(F) <= sizeof(taskData),
                  "given task type exceeds maximum size allowed for deque\n");
    static_assert(std::is_trivially_destructible<F>::value,
                  "given task type must be trivially destructible\n");
    metadata.stealer.store(0, std::memory_order_relaxed);

    new (taskData) Callable<F>(std::forward<F>(f));

    assert(static_cast<CallableBase*>(reinterpret_cast<Callable<F>*>(
               taskData)) == reinterpret_cast<CallableBase*>(taskData));
  }

  void cancel() {
    metadata.stealer.fetch_or(kCancelFlag, std::memory_order_release);
  }

  /// run task as owner, if not cancelled
  void run() {
    if (metadata.stealer.load(std::memory_order_relaxed) == 0) getCallable()();
  }

  /// request notification of the owner when the task is finished.
  /// Should be called while the owner holds its wait mutex
  /// and only after getStealerIfUnfinished() has been called.
  /// Returns true if the notification was set and false if it was not set
  /// because the task was finished in the meantime.
  bool requestNotifyWhenFinished(HighsSplitDeque* owner,
                                 HighsSplitDeque* stealer) {
    uintptr_t xormask = reinterpret_cast<uintptr_t>(owner) ^
                        reinterpret_cast<uintptr_t>(stealer);
    uintptr_t state =
        metadata.stealer.fetch_xor(xormask, std::memory_order_acquire);

    assert(stealer != nullptr);

    return (state & kFinishedFlag) == 0;
  }

  /// check if task is finished
  bool isFinished() const {
    uintptr_t state = metadata.stealer.load(std::memory_order_acquire);
    return state & kFinishedFlag;
  }

  /// check if task is cancelled
  bool isCancelled() const {
    uintptr_t state = metadata.stealer.load(std::memory_order_relaxed);
    return state & kCancelFlag;
  }

  /// get the stealer of a stolen task, or nullptr if the stealer finished
  /// executing the task. Spin waits for the stealer to have started executing
  /// the task if necessary.
  HighsSplitDeque* getStealerIfUnfinished(bool* cancelled = nullptr) {
    uintptr_t state = metadata.stealer.load(std::memory_order_acquire);
    if (state & kFinishedFlag)
      return nullptr;
    else {
      while ((state & ~kCancelFlag) == 0) {
        // the task has been stolen, but the stealer has not yet started
        // executing the task in this case, yield and check again in a spin
        // loop until the stealer executes the task and becomes visible to
        // this thread
        HighsSpinMutex::yieldProcessor();
        state = metadata.stealer.load(std::memory_order_acquire);
      }
    }

    if (state & kFinishedFlag) return nullptr;

    if (cancelled) *cancelled = state & kCancelFlag;

    return reinterpret_cast<HighsSplitDeque*>(state & kPtrMask);
  }
};

#endif
