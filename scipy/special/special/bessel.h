#pragma once

#include "specfun.h"

namespace special {
namespace detail {

    template <typename T>
    void itika(T x, T *ti, T *tk) {

        // =======================================================
        // Purpose: Integrate modified Bessel functions I0(t) and
        //          K0(t) with respect to t from 0 to x
        // Input :  x  --- Upper limit of the integral  ( x ≥ 0 )
        // Output:  TI --- Integration of I0(t) from 0 to x
        //          TK --- Integration of K0(t) from 0 to x
        // =======================================================

        int k;
        T rc1, rc2, e0, b1, b2, rs, r, tw, x2;
        static const T a[10] = {0.625,
                                1.0078125,
                                2.5927734375,
                                9.1868591308594,
                                4.1567974090576e+1,
                                2.2919635891914e+2,
                                1.491504060477e+3,
                                1.1192354495579e+4,
                                9.515939374212e+4,
                                9.0412425769041e+5};
        const T pi = 3.141592653589793;
        const T el = 0.5772156649015329;

        if (x == 0.0) {
            *ti = 0.0;
            *tk = 0.0;
            return;
        } else if (x < 20.0) {
            x2 = x * x;
            *ti = 1.0;
            r = 1.0;
            for (k = 1; k <= 50; k++) {
                r = 0.25 * r * (2 * k - 1.0) / (2 * k + 1.0) / (k * k) * x2;
                *ti += r;
                if (fabs(r / (*ti)) < 1.0e-12) {
                    break;
                }
            }
            *ti *= x;
        } else {
            x2 = 0.0;
            *ti = 1.0;
            r = 1.0;
            for (k = 1; k <= 10; k++) {
                r = r / x;
                *ti += a[k - 1] * r;
            }
            rc1 = 1.0 / sqrt(2.0 * pi * x);
            *ti = rc1 * exp(x) * (*ti);
        }

        if (x < 12.0) {
            e0 = el + log(x / 2.0);
            b1 = 1.0 - e0;
            b2 = 0.0;
            rs = 0.0;
            r = 1.0;
            tw = 0.0;
            for (k = 1; k <= 50; k++) {
                r = 0.25 * r * (2 * k - 1.0) / (2 * k + 1.0) / (k * k) * x2;
                b1 += r * (1.0 / (2 * k + 1) - e0);
                rs += 1.0 / k;
                b2 += r * rs;
                *tk = b1 + b2;
                if (fabs((*tk - tw) / (*tk)) < 1.0e-12) {
                    break;
                }
                tw = *tk;
            }
            *tk *= x;
        } else {
            *tk = 1.0;
            r = 1.0;
            for (k = 1; k <= 10; k++) {
                r = -r / x;
                *tk = *tk + a[k - 1] * r;
            }
            rc2 = sqrt(pi / (2.0 * x));
            *tk = pi / 2.0 - rc2 * (*tk) * exp(-x);
        }
        return;
    }

    template <typename T>
    void itjya(T x, T *tj, T *ty) {
        int k;
        T a[18], a0, a1, af, bf, bg, r, r2, rc, rs, ty1, ty2, x2, xp;
        const T pi = 3.141592653589793;
        const T el = 0.5772156649015329;
        const T eps = 1.0e-12;

        if (x == 0.0) {
            *tj = 0.0;
            *ty = 0.0;
        } else if (x <= 20.0) {
            x2 = x * x;
            *tj = x;
            r = x;
            for (k = 1; k < 61; k++) {
                r = -0.25 * r * (2.0 * k - 1.0) / (2.0 * k + 1.0) / (k * k) * x2;
                *tj += r;
                if (fabs(r) < fabs(*tj) * eps) {
                    break;
                }
            }
            ty1 = (el + log(x / 2.0)) * (*tj);
            rs = 0.0;
            ty2 = 1.0;
            r = 1.0;
            for (k = 1; k < 61; k++) {
                r = -0.25 * r * (2.0 * k - 1.0) / (2.0 * k + 1.0) / (k * k) * x2;
                rs += 1.0 / k;
                r2 = r * (rs + 1.0 / (2.0 * k + 1.0));
                ty2 += r2;
                if (fabs(r2) < fabs(ty2) * eps) {
                    break;
                }
            }
            *ty = (ty1 - x * ty2) * 2.0 / pi;
        } else {
            a0 = 1.0;
            a1 = 5.0 / 8.0;
            a[0] = a1;

            for (int k = 1; k <= 16; k++) {
                af = ((1.5 * (k + 0.5) * (k + 5.0 / 6.0) * a1 - 0.5 * (k + 0.5) * (k + 0.5) * (k - 0.5) * a0)) /
                     (k + 1.0);
                a[k] = af;
                a0 = a1;
                a1 = af;
            }
            bf = 1.0;
            r = 1.0;
            for (int k = 1; k <= 8; k++) {
                r = -r / (x * x);
                bf += a[2 * k - 1] * r;
            }
            bg = a[0] / x;
            r = 1.0 / x;
            for (int k = 1; k <= 8; k++) {
                r = -1.0 / (x * x);
                bg += a[2 * k] * r;
            }
            xp = x + 0.25 * pi;
            rc = sqrt(2.0 / (pi * x));
            *tj = 1.0 - rc * (bf * cos(xp) + bg * sin(xp));
            *ty = rc * (bg * cos(xp) - bf * sin(xp));
        }
        return;
    }

    template <typename T>
    void ittika(T x, T *tti, T *ttk) {

        // =========================================================
        // Purpose: Integrate [I0(t)-1]/t with respect to t from 0
        //          to x, and K0(t)/t with respect to t from x to ∞
        // Input :  x   --- Variable in the limits  ( x ≥ 0 )
        // Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
        //          TTK --- Integration of K0(t)/t from x to ∞
        // =========================================================

        int k;
        T b1, e0, r, r2, rc, rs;
        const T pi = 3.141592653589793;
        const T el = 0.5772156649015329;
        static const T c[8] = {1.625,           4.1328125,       1.45380859375,   6.553353881835,
                               3.6066157150269, 2.3448727161884, 1.7588273098916, 1.4950639538279};

        if (x == 0.0) {
            *tti = 0.0;
            *ttk = 1.0e+300;
            return;
        }
        if (x < 40.0) {
            *tti = 1.0;
            r = 1.0;
            for (int k = 2; k <= 50; k++) {
                r = 0.25 * r * (k - 1.0) / (k * k * k) * x * x;
                *tti += r;
                if (fabs(r / (*tti)) < 1.0e-12) {
                    break;
                }
            }
            *tti *= 0.125 * x * x;
        } else {
            *tti = 1.0;
            r = 1.0;
            for (int k = 1; k <= 8; k++) {
                r = r / x;
                *tti += c[k - 1] * r;
            }
            rc = x * sqrt(2.0 * pi * x);
            *tti = (*tti) * exp(x) / rc;
        }
        if (x <= 12.0) {
            e0 = (0.5 * log(x / 2.0) + el) * log(x / 2.0) + pi * pi / 24.0 + 0.5 * el * el;
            b1 = 1.5 - (el + log(x / 2.0));
            rs = 1.0;
            r = 1.0;
            for (k = 2; k <= 50; k++) {
                r = 0.25 * r * (k - 1.0) / (k * k * k) * x * x;
                rs += 1.0 / k;
                r2 = r * (rs + 1.0 / (2.0 * k) - (el + log(x / 2.0)));
                b1 += r2;
                if (fabs(r2 / b1) < 1.0e-12) {
                    break;
                }
            }
            *ttk = e0 - 0.125 * x * x * b1;
        } else {
            *ttk = 1.0;
            r = 1.0;
            for (k = 1; k <= 8; k++) {
                r = -r / x;
                *ttk += c[k - 1] * r;
            }
            rc = x * sqrt(2.0 / (pi * x));
            *ttk = (*ttk) * exp(-x) / rc;
        }
        return;
    }

    template <typename T>
    void ittjya(T x, T *ttj, T *tty) {

        // =========================================================
        // Purpose: Integrate [1-J0(t)]/t with respect to t from 0
        //          to x, and Y0(t)/t with respect to t from x to ∞
        // Input :  x   --- Variable in the limits  ( x ≥ 0 )
        // Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
        //          TTY --- Integration of Y0(t)/t from x to ∞
        // =========================================================

        int k, l;
        T a0, bj0, bj1, by0, by1, e0, b1, g0, g1, px, qx, r, r0, r1, r2, rs, t, vt, xk;
        const T pi = 3.141592653589793;
        const T el = 0.5772156649015329;

        if (x == 0.0) {
            *ttj = 0.0;
            *tty = -1.0e+300;
        } else if (x <= 20.0) {
            *ttj = 1.0;
            r = 1.0;
            for (k = 2; k <= 100; k++) {
                r = -0.25 * r * (k - 1.0) / (k * k * k) * x * x;
                *ttj += r;
                if (fabs(r) < fabs(*ttj) * 1.0e-12) {
                    break;
                }
            }
            *ttj *= 0.125 * x * x;
            e0 = 0.5 * (pi * pi / 6.0 - el * el) - (0.5 * log(x / 2.0) + el) * log(x / 2.0);
            b1 = el + log(x / 2.0) - 1.5;
            rs = 1.0;
            r = -1.0;
            for (k = 2; k <= 100; k++) {
                r = -0.25 * r * (k - 1.0) / (k * k * k) * x * x;
                rs = rs + 1.0 / k;
                r2 = r * (rs + 1.0 / (2.0 * k) - (el + log(x / 2.0)));
                b1 = b1 + r2;
                if (fabs(r2) < fabs(b1) * 1.0e-12) {
                    break;
                }
            }
            *tty = 2.0 / pi * (e0 + 0.125 * x * x * b1);
        } else {
            a0 = sqrt(2.0 / (pi * x));
            bj0 = 0.0;
            by0 = 0.0;
            bj1 = 0.0;
            for (l = 0; l <= 1; l++) {
                vt = 4.0 * l * l;
                px = 1.0;
                r = 1.0;
                for (k = 1; k <= 14; k++) {
                    r = -0.0078125 * r * (vt - pow((4.0 * k - 3.0), 2)) / (x * k) * (vt - pow((4.0 * k - 1.0), 2)) /
                        ((2.0 * k - 1.0) * x);
                    px += r;
                    if (fabs(r) < fabs(px) * 1.0e-12) {
                        break;
                    }
                }
                qx = 1.0;
                r = 1.0;
                for (k = 1; k <= 14; k++) {
                    r = -0.0078125 * r * (vt - pow((4.0 * k - 1.0), 2)) / (x * k) * (vt - pow((4.0 * k + 1.0), 2)) /
                        ((2.0 * k + 1.0) * x);
                    qx += r;
                    if (fabs(r) < fabs(qx) * 1.0e-12) {
                        break;
                    }
                }
                qx = 0.125 * (vt - 1.0) / x * qx;
                xk = x - (0.25 + 0.5 * l) * pi;
                bj1 = a0 * (px * cos(xk) - qx * sin(xk));
                by1 = a0 * (px * sin(xk) + qx * cos(xk));
                if (l == 0) {
                    bj0 = bj1;
                    by0 = by1;
                }
            }
            t = 2.0 / x;
            g0 = 1.0;
            r0 = 1.0;
            for (k = 1; k <= 10; k++) {
                r0 = -r0 * k * k * t * t;
                g0 += r0;
            }
            g1 = 1.0;
            r1 = 1.0;
            for (k = 1; k <= 10; k++) {
                r1 = -r1 * k * (k + 1.0) * t * t;
                g1 += r1;
            }
            *ttj = 2.0 * g1 * bj0 / (x * x) - g0 * bj1 / x + el + log(x / 2.0);
            *tty = 2.0 * g1 * by0 / (x * x) - g0 * by1 / x;
        }
        return;
    }

} // namespace detail

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
    detail::itjya(x, j0int, y0int);
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
    detail::ittjya(x, j0int, y0int);
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
    detail::itika(x, i0int, k0int);
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
    detail::ittika(x, i0int, k0int);
    if (flag) {
        *k0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
    }
}

template <typename T, typename OutputVec1, typename OutputVec2>
void rctj(T x, int *nm, OutputVec1 rj, OutputVec2 dj) {

    // ========================================================
    // Purpose: Compute Riccati-Bessel functions of the first
    //          kind and their derivatives
    // Input:   x --- Argument of Riccati-Bessel function
    //          n --- Order of jn(x)  ( n = 0,1,2,... )
    // Output:  RJ(n) --- x·jn(x)
    //          DJ(n) --- [x·jn(x)]'
    //          NM --- Highest order computed
    // Routines called:
    //          MSTA1 and MSTA2 for computing the starting
    //          point for backward recurrence
    // ========================================================

    int n = rj.extent(0) - 1;

    int k, m;
    T cs, f, f0, f1, rj0, rj1;

    *nm = n;
    if (fabs(x) < 1.0e-100) {
        for (int k = 0; k <= n; k++) {
            rj[k] = 0.0;
            dj[k] = 0.0;
        }
        dj[0] = 1.0;
        return;
    }
    rj[0] = sin(x);
    rj[1] = rj[0] / x - cos(x);
    rj0 = rj[0];
    rj1 = rj[1];
    cs = 0.0;
    f = 0.0;

    if (n >= 2) {
        m = specfun::msta1(x, 200);
        if (m < n) {
            *nm = m;
        } else {
            m = specfun::msta2(x, n, 15);
        }

        f0 = 0.0;
        f1 = 1.0e-100;

        for (k = m; k >= 0; k--) {
            f = (2.0 * k + 3.0) * f1 / x - f0;
            if (k <= *nm) {
                rj[k] = f;
            }
            f0 = f1;
            f1 = f;
        }
        cs = (fabs(rj0) > fabs(rj1) ? rj0 / f : rj1 / f0);
        for (k = 0; k <= *nm; k++) {
            rj[k] = cs * rj[k];
        }
    }
    dj[0] = cos(x);
    for (int k = 1; k <= *nm; k++) {
        dj[k] = -k * rj[k] / x + rj[k - 1];
    }
}

template <typename T, typename OutputVec1, typename OutputVec2>
void rctj(T x, OutputVec1 rj, OutputVec2 dj) {
    int nm;
    rctj(x, &nm, rj, dj);
}

template <typename T, typename OutputVec1, typename OutputVec2>
void rcty(T x, int *nm, OutputVec1 ry, OutputVec2 dy) {

    // ========================================================
    // Purpose: Compute Riccati-Bessel functions of the second
    //          kind and their derivatives
    // Input:   x --- Argument of Riccati-Bessel function
    //          n --- Order of yn(x)
    // Output:  RY(n) --- x·yn(x)
    //          DY(n) --- [x·yn(x)]'
    //          NM --- Highest order computed
    // ========================================================

    int n = ry.extent(0) - 1;

    int k;
    T rf0, rf1, rf2;
    *nm = n;
    if (x < 1.0e-60) {
        for (k = 0; k <= n; k++) {
            ry[k] = -1.0e+300;
            dy[k] = 1.0e+300;
        }
        ry[0] = -1.0;
        dy[0] = 0.0;
        return;
    }

    ry[0] = -cos(x);
    ry[1] = ry[0] / x - sin(x);
    rf0 = ry[0];
    rf1 = ry[1];

    for (k = 2; k <= n; k++) {
        rf2 = (2.0 * k - 1.0) * rf1 / x - rf0;
        if (fabs(rf2) > 1.0e+300) {
            break;
        }
        ry[k] = rf2;
        rf0 = rf1;
        rf1 = rf2;
    }

    *nm = k - 1;
    dy[0] = sin(x);
    for (k = 1; k <= *nm; k++) {
        dy[k] = -k * ry[k] / x + ry[k - 1];
    }
    return;
}

template <typename T, typename OutputVec1, typename OutputVec2>
void rcty(T x, OutputVec1 ry, OutputVec2 dy) {
    int nm;
    rcty(x, &nm, ry, dy);
}

} // namespace special
