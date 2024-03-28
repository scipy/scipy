#pragma once

#include "specfun.h"

namespace special {

template <typename T>
void oblate_aswfa_nocv(T m, T n, T c, T x, T *s1, T *s1d);

template <>
inline void oblate_aswfa_nocv(double m, double n, double c, double x, double *s1, double *s1d) {
    *s1 = oblate_aswfa_nocv(m, n, c, x, s1d);
}

template <>
inline void oblate_aswfa_nocv(float mf, float nf, float cf, float xf, float *s1f, float *s1df) {
    double m = mf;
    double n = nf;
    double c = cf;
    double x = xf;
    double s1;
    double s1d;
    oblate_aswfa_nocv(m, n, c, x, &s1, &s1d);

    *s1f = s1;
    *s1df = s1d;
}

} // namespace special
