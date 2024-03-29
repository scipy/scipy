#pragma once

#include "config.h"

namespace special {
namespace detail {

    template <typename T>
    void itairy(T x, T &apt, T &bpt, T &ant, T &bnt) {

        // ======================================================
        // Purpose: Compute the integrals of Airy fnctions with
        //          respect to t from 0 and x ( x ≥ 0 )
        // Input  : x   --- Upper limit of the integral
        // Output : APT --- Integration of Ai(t) from 0 and x
        //          BPT --- Integration of Bi(t) from 0 and x
        //          ANT --- Integration of Ai(-t) from 0 and x
        //          BNT --- Integration of Bi(-t) from 0 and x
        // ======================================================

        int k, l;
        T fx, gx, r, su1, su2, su3, su4, su5, su6, xp6, xe, xr1, xr2;

        const T pi = 3.141592653589793;
        const T c1 = 0.355028053887817;
        const T c2 = 0.258819403792807;
        const T sr3 = 1.732050807568877;
        const T q0 = 1.0 / 3.0;
        const T q1 = 2.0 / 3.0;
        const T q2 = 1.4142135623730951;
        const T eps = 1e-5;
        static const T a[16] = {
            0.569444444444444,     0.891300154320988,     0.226624344493027e+01, 0.798950124766861e+01,
            0.360688546785343e+02, 0.198670292131169e+03, 0.129223456582211e+04, 0.969483869669600e+04,
            0.824184704952483e+05, 0.783031092490225e+06, 0.822210493622814e+07, 0.945557399360556e+08,
            0.118195595640730e+10, 0.159564653040121e+11, 0.231369166433050e+12, 0.358622522796969e+13};

        if (x == 0.0) {
            apt = 0.0;
            bpt = 0.0;
            ant = 0.0;
            bnt = 0.0;
        } else {
            if (std::abs(x) <= 9.25) {
                for (l = 0; l < 2; l++) {
                    x = x * std::pow(-1, l);
                    fx = x;
                    r = x;
                    for (k = 1; k < 41; k++) {
                        r = r * (3.0 * k - 2.0) / (3.0 * k + 1.0) * x / (3.0 * k) * x / (3.0 * k - 1.0) * x;
                        fx += r;
                        if (std::abs(r) < std::abs(fx) * eps) {
                            break;
                        }
                    }
                    gx = 0.5 * x * x;
                    r = gx;
                    for (k = 1; k < 41; k++) {
                        r = r * (3.0 * k - 1.0) / (3.0 * k + 2.0) * x / (3.0 * k) * x / (3.0 * k + 1.0) * x;
                        gx += r;
                        if (std::abs(r) < std::abs(gx) * eps) {
                            break;
                        }
                    }
                    ant = c1 * fx - c2 * gx;
                    bnt = sr3 * (c1 * fx + c2 * gx);
                    if (l == 0) {
                        apt = ant;
                        bpt = bnt;
                    } else {
                        ant = -ant;
                        bnt = -bnt;
                        x = -x;
                    }
                }
            } else {
                xe = x * std::sqrt(x) / 1.5;
                xp6 = 1.0 / std::sqrt(6.0 * pi * xe);
                su1 = 1.0;
                r = 1.0;
                xr1 = 1.0 / xe;
                for (k = 1; k < 17; k++) {
                    r = -r * xr1;
                    su1 += a[k - 1] * r;
                }
                su2 = 1.0;
                r = 1.0;
                for (k = 1; k < 17; k++) {
                    r = r * xr1;
                    su2 += a[k - 1] * r;
                }
                apt = q0 - std::exp(-xe) * xp6 * su1;
                bpt = 2.0 * std::exp(xe) * xp6 * su2;
                su3 = 1.0;
                r = 1.0;
                xr2 = 1.0 / (xe * xe);
                for (k = 1; k < 9; k++) {
                    r = -r * xr2;
                    su3 += a[2 * k - 1] * r;
                }
                su4 = a[0] * xr1;
                r = xr1;
                for (k = 1; k < 8; k++) {
                    r = -r * xr2;
                    su4 += a[2 * k] * r;
                }
                su5 = su3 + su4;
                su6 = su3 - su4;
                ant = q1 - q2 * xp6 * (su5 * std::cos(xe) - su6 * std::sin(xe));
                bnt = q2 * xp6 * (su5 * std::sin(xe) + su6 * std::cos(xe));
            }
        }
    }

} // namespace detail

template <typename T>
void itairy(T x, T *apt, T *bpt, T *ant, T *bnt) {
    bool x_signbit = std::signbit(x);
    if (x_signbit) {
        x = -x;
    }

    detail::itairy(x, *apt, *bpt, *ant, *bnt);
    if (x_signbit) { /* negative limit -- switch signs and roles */
        T tmp = *apt;
        *apt = -*ant;
        *ant = -tmp;

        tmp = *bpt;
        *bpt = -*bnt;
        *bnt = -tmp;
    }
}

} // namespace special
