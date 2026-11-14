#include "Control.h"

#include <cassert>

#include "Status.h"
#include "parallel/HighsParallel.h"

namespace hipo {

void Control::setCallback(HighsCallback& callback) { callback_ = &callback; }
void Control::setTimer(const HighsTimer& timer) { timer_ = &timer; }
void Control::setOptions(const Options& options) { options_ = &options; }

HighsCallback* Control::callback() const { return callback_; }

double Control::elapsed() const { return timer_ ? timer_->read() : -1.0; }

Int Control::interruptCheck(const Int ipm_iteration_count) const {
  HighsTaskExecutor::getThisWorkerDeque()->checkInterrupt();

  if (options_ && options_->time_limit >= 0 &&
      elapsed() > options_->time_limit) {
    return kStatusTimeLimit;
  }

  if (callback_) {
    if (callback_->user_callback && callback_->active[kCallbackIpmInterrupt]) {
      callback_->clearHighsCallbackOutput();
      callback_->data_out.ipm_iteration_count = ipm_iteration_count;
      if (callback_->callbackAction(kCallbackIpmInterrupt, "IPM interrupt"))
        return kStatusUserInterrupt;
    }
  }

  return 0;
}

}  // namespace hipo