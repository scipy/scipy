// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#ifndef IPX_TIMER_H_
#define IPX_TIMER_H_

#include <chrono>

namespace ipx {

class Timer {
public:
    Timer();
    double Elapsed() const;
    void Reset();

private:
    typedef std::chrono::time_point<std::chrono::high_resolution_clock>
        TimePoint;
    static TimePoint tic();
    static double toc(TimePoint start);
    TimePoint t0_;
};

}  // namespace ipx

#endif  // IPX_TIMER_H_
