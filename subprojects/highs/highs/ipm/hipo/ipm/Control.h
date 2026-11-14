#ifndef HIPO_CONTROL_H
#define HIPO_CONTROL_H

#include "Options.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "lp_data/HighsCallback.h"
#include "util/HighsTimer.h"

// Interface to Highs timer and callback

namespace hipo {

class Control {
  // disallow copy and copy construction
  Control& operator=(const Control&) = delete;
  Control(const Control&) = delete;

  // disallow move and move construction
  Control& operator=(Control&&) = delete;
  Control(const Control&&) = delete;

  HighsCallback* callback_ = nullptr;
  const HighsTimer* timer_ = nullptr;
  const Options* options_ = nullptr;

 public:
  Control() = default;

  void setCallback(HighsCallback& callback);
  void setTimer(const HighsTimer& timer);
  void setOptions(const Options& options);

  HighsCallback* callback() const;

  double elapsed() const;
  Int interruptCheck(const Int ipm_iteration_count = -1) const;
};

}  // namespace hipo

#endif