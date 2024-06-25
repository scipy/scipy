#pragma once

#include "config.h"

namespace special {
namespace detail {

    inline double itth0(double x) {

        // ===========================================================
        // Purpose: Evaluate the integral H0(t)/t with respect to t
        //          from x to infinity
        // Input :  x   --- Lower limit  ( x ≥ 0 )
        // Output:  TTH --- Integration of H0(t)/t from x to infinity
        // ===========================================================

        int k;
        double f0, g0, r, s, t, tth, tty, xt;
        const double pi = 3.141592653589793;
        s = 1.0;
        r = 1.0;
        if (x < 24.5) {
            for (k = 1; k < 61; k++) {
                r = -r * x * x * (2.0 * k - 1.0) / pow(2.0 * k + 1.0, 3);
                s += r;
                if (fabs(r) < fabs(s) * 1.0e-12) {
                    break;
                }
            }
            tth = pi / 2.0 - 2.0 / pi * x * s;
        } else {
            for (k = 1; k < 11; k++) {
                r = -r * pow(2.0 * k - 1.0, 3) / ((2.0 * k + 1.0) * x * x);
                s += r;
                if (fabs(r) < fabs(s) * 1.0e-12) {
                    break;
                }
            }
            tth = 2.0 / (pi * x) * s;
            t = 8.0 / x;
            xt = x + 0.25 * pi;
            f0 = (((((0.18118e-2 * t - 0.91909e-2) * t + 0.017033) * t - 0.9394e-3) * t - 0.051445) * t - 0.11e-5) * t +
                 0.7978846;
            g0 =
                (((((-0.23731e-2 * t + 0.59842e-2) * t + 0.24437e-2) * t - 0.0233178) * t + 0.595e-4) * t + 0.1620695) *
                t;
            tty = (f0 * sin(xt) - g0 * cos(xt)) / (sqrt(x) * x);
            tth = tth + tty;
        }
        return tth;
    }

    inline double itsh0(double x) {

        // ===================================================
        // Purpose: Evaluate the integral of Struve function
        //          H0(t) with respect to t from 0 and x
        // Input :  x   --- Upper limit  ( x ≥ 0 )
        // Output:  TH0 --- Integration of H0(t) from 0 and x
        // ===================================================

        int k;
        double a[25], a0, a1, af, bf, bg, r, rd, s, s0, th0, ty, xp;
        const double pi = 3.141592653589793;
        const double el = 0.57721566490153;

        r = 1.0;
        if (x <= 30.0) {
            s = 0.5;
            for (k = 1; k < 101; k++) {
                rd = 1.0;
                if (k == 1) {
                    rd = 0.5;
                }
                r = -r * rd * k / (k + 1.0) * pow(x / (2.0 * k + 1.0), 2);
                s += r;
                if (fabs(r) < fabs(s) * 1.0e-12) {
                    break;
                }
            }
            th0 = 2.0 / pi * x * x * s;
        } else {
            s = 1.0;
            for (k = 1; k < 13; k++) {
                r = -r * k / (k + 1.0) * pow((2.0 * k + 1.0) / x, 2);
                s += r;
                if (fabs(r) < fabs(s) * 1.0e-12) {
                    break;
                }
            }
            s0 = s / (pi * x * x) + 2.0 / pi * (log(2.0 * x) + el);
            a0 = 1.0;
            a1 = 5.0 / 8.0;
            a[0] = a1;
            for (k = 1; k < 21; k++) {
                af = ((1.5 * (k + 0.5) * (k + 5.0 / 6.0) * a1 - 0.5 * (k + 0.5) * (k + 0.5) * (k - 0.5) * a0)) /
                     (k + 1.0);
                a[k] = af;
                a0 = a1;
                a1 = af;
            }
            bf = 1.0;
            r = 1.0;
            for (k = 1; k < 11; k++) {
                r = -r / (x * x);
                bf += a[2 * k - 1] * r;
            }
            bg = a[0] * x;
            r = 1.0 / x;
            for (k = 1; k < 10; k++) {
                r = -r / (x * x);
                bg += a[2 * k] * r;
            }
            xp = x + 0.25 * pi;
            ty = sqrt(2.0 / (pi * x)) * (bg * cos(xp) - bf * sin(xp));
            th0 = ty + s0;
        }
        return th0;
    }

    inline double itsl0(double x) {

        // ===========================================================
        // Purpose: Evaluate the integral of modified Struve function
        //          L0(t) with respect to t from 0 to x
        // Input :  x   --- Upper limit  ( x ≥ 0 )
        // Output:  TL0 --- Integration of L0(t) from 0 to x
        // ===========================================================

        int k;
        double a[18], a0, a1, af, r, rd, s, s0, ti, tl0;
        const double pi = 3.141592653589793;
        const double el = 0.57721566490153;
        r = 1.0;
        if (x <= 20.0) {
            s = 0.5;
            for (k = 1; k < 101; k++) {
                rd = 1.0;
                if (k == 1) {
                    rd = 0.5;
                }
                r = r * rd * k / (k + 1.0) * pow(x / (2.0 * k + 1.0), 2);
                s += r;
                if (fabs(r / s) < 1.0e-12) {
                    break;
                }
            }
            tl0 = 2.0 / pi * x * x * s;
        } else {
            s = 1.0;
            for (k = 1; k < 11; k++) {
                r = r * k / (k + 1.0) * pow((2.0 * k + 1.0) / x, 2);
                s += r;
                if (fabs(r / s) < 1.0e-12) {
                    break;
                }
            }
            s0 = -s / (pi * x * x) + 2.0 / pi * (log(2.0 * x) + el);
            a0 = 1.0;
            a1 = 5.0 / 8.0;
            a[0] = a1;
            for (k = 1; k < 11; k++) {
                af = ((1.5 * (k + .50) * (k + 5.0 / 6.0) * a1 - 0.5 * pow(k + 0.5, 2) * (k - 0.5) * a0)) / (k + 1.0);
                a[k] = af;
                a0 = a1;
                a1 = af;
            }
            ti = 1.0;
            r = 1.0;
            for (k = 1; k < 11; k++) {
                r = r / x;
                ti += a[k - 1] * r;
            }
            tl0 = ti / sqrt(2 * pi * x) * exp(x) + s0;
        }
        return tl0;
    }

} // namespace detail

template <typename T>
T itstruve0(T x) {
    if (x < 0) {
        x = -x;
    }

    T out = detail::itsh0(x);
    SPECFUN_CONVINF("itstruve0", out);
    return out;
}

template <typename T>
T it2struve0(T x) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }

    T out = detail::itth0(x);
    SPECFUN_CONVINF("it2struve0", out);
    if (flag) {
        out = M_PI - out;
    }
    return out;
}

template <typename T>
T itmodstruve0(T x) {
    if (x < 0) {
        x = -x;
    }

    T out = detail::itsl0(x);
    SPECFUN_CONVINF("itmodstruve0", out);
    return out;
}

} // namespace special
