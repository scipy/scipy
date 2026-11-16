#ifndef IPX_TIMER_H_
#define IPX_TIMER_H_

#include <chrono>

namespace ipx {
    using namespace std::chrono;

class Timer {
public:
    Timer(const double offset=0);
    double Elapsed() const { return toc(t0_); }
    void Reset(const bool first = false);

private:
    static double toc(double start) { return read() - start; }
    static double read() { return duration_cast<duration<double>>(high_resolution_clock::now().time_since_epoch()).count(); }
    double t0_;
public:
    double offset_;
};

}  // namespace ipx

#endif  // IPX_TIMER_H_
