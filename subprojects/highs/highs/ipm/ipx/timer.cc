#include "ipm/ipx/timer.h"
#include <cstdio>
namespace ipx {

Timer::Timer(const double offset)
    : offset_(offset) {
    Reset(true);
}

void Timer::Reset(const bool first) {
  if (!first) offset_ -= t0_;
  t0_ = read();
  if (!first) offset_ += t0_;
}

}  // namespace ipx
