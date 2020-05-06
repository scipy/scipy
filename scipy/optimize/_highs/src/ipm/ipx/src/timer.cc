// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "timer.h"

namespace ipx {

Timer::Timer() {
    Reset();
}

double Timer::Elapsed() const {
    return toc(t0_);
}

void Timer::Reset() {
    t0_ = tic();
}

Timer::TimePoint Timer::tic() {
    return std::chrono::high_resolution_clock::now();
}

double Timer::toc(TimePoint start) {
    TimePoint end = tic();
    std::chrono::duration<double> diff = end-start;
    return diff.count();
}

}  // namespace ipx
