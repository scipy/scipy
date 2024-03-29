#pragma once

#include "specfun.h"

namespace special {

/* Integrals of bessel functions */

/* int(j0(t),t=0..x) */
/* int(y0(t),t=0..x) */

template <typename T>
void it1j0y0(T x, T *j0int, T *y0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::itjya(x, j0int, y0int);
    if (flag) {
        *j0int = -(*j0int);
        *y0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
    }
}

/* int((1-j0(t))/t,t=0..x) */
/* int(y0(t)/t,t=x..inf) */

template <typename T>
void it2j0y0(T x, T *j0int, T *y0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::ittjya(x, j0int, y0int);
    if (flag) {
        *y0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
    }
}

/* Integrals of modified bessel functions */

template <typename T>
void it1i0k0(T x, T *i0int, T *k0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::itika(x, i0int, k0int);
    if (flag) {
        *i0int = -(*i0int);
        *k0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
    }
}

template <typename T>
void it2i0k0(T x, T *i0int, T *k0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::ittika(x, i0int, k0int);
    if (flag) {
        *k0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
    }
}

} // namespace special
