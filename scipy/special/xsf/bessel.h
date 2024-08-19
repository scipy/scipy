#pragma once

#include "amos.h"
#include "cephes/jv.h"
#include "cephes/scipy_iv.h"
#include "cephes/yv.h"
#include "error.h"
#include "specfun.h"
#include "trig.h"

extern "C" double cephes_iv(double v, double x);

namespace xsf {
namespace detail {

    template <typename T>
    std::complex<T> rotate(std::complex<T> z, T v) {
        T c = cospi(v);
        T s = sinpi(v);

        return {c * std::real(z) - s * std::imag(z), s * std::real(z) + c * std::imag(z)};
    }

    template <typename T>
    std::complex<T> rotate_jy(std::complex<T> j, std::complex<T> y, T v) {
        T c = cospi(v);
        T s = sinpi(v);

        return c * j - s * y;
    }

    template <typename T>
    int reflect_jy(std::complex<T> *jy, T v) {
        /* NB: Y_v may be huge near negative integers -- so handle exact
         *     integers carefully
         */
        int i;
        if (v != floor(v))
            return 0;

        i = v - 16384.0 * floor(v / 16384.0);
        if (i & 1) {
            *jy = -(*jy);
        }
        return 1;
    }

    template <typename T>
    int reflect_i(std::complex<T> *ik, T v) {
        if (v != floor(v)) {
            return 0;
        }

        return 1; /* I is symmetric for integer v */
    }

    template <typename T>
    std::complex<T> rotate_i(std::complex<T> i, std::complex<T> k, T v) {
        T s = std::sin(v * M_PI) * (2 / M_PI);

        return i + s * k;
    }

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
        static const T a[10] = {
            0.625,
            1.0078125,
            2.5927734375,
            9.1868591308594,
            4.1567974090576e+1,
            2.2919635891914e+2,
            1.491504060477e+3,
            1.1192354495579e+4,
            9.515939374212e+4,
            9.0412425769041e+5
        };
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
void it1j0y0(T x, T &j0int, T &y0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    detail::itjya(x, &j0int, &y0int);
    if (flag) {
        j0int = -j0int;
        y0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
    }
}

/* int((1-j0(t))/t,t=0..x) */
/* int(y0(t)/t,t=x..inf) */

template <typename T>
void it2j0y0(T x, T &j0int, T &y0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    detail::ittjya(x, &j0int, &y0int);
    if (flag) {
        y0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
    }
}

/* Integrals of modified bessel functions */

template <typename T>
void it1i0k0(T x, T &i0int, T &k0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    detail::itika(x, &i0int, &k0int);
    if (flag) {
        i0int = -i0int;
        k0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
    }
}

template <typename T>
void it2i0k0(T x, T &i0int, T &k0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    detail::ittika(x, &i0int, &k0int);
    if (flag) {
        k0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
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

inline std::complex<double> cyl_bessel_je(double v, std::complex<double> z) {
    int n = 1;
    int kode = 2;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy_j, cy_y;

    cy_j.real(NAN);
    cy_j.imag(NAN);
    cy_y.real(NAN);
    cy_y.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy_j;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besj(z, v, kode, n, &cy_j, &ierr);
    set_error_and_nan("jve:", ierr_to_sferr(nz, ierr), cy_j);
    if (sign == -1) {
        if (!detail::reflect_jy(&cy_j, v)) {
            nz = amos::besy(z, v, kode, n, &cy_y, &ierr);
            set_error_and_nan("jve(yve):", ierr_to_sferr(nz, ierr), cy_y);
            cy_j = detail::rotate_jy(cy_j, cy_y, v);
        }
    }
    return cy_j;
}

inline std::complex<float> cyl_bessel_je(float v, std::complex<float> x) {
    return static_cast<std::complex<float>>(cyl_bessel_je(static_cast<double>(v), static_cast<std::complex<double>>(x))
    );
}

template <typename T>
T cyl_bessel_je(T v, T x) {
    if (v != floor(v) && x < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }

    return std::real(cyl_bessel_je(v, std::complex(x)));
}

inline std::complex<double> cyl_bessel_ye(double v, std::complex<double> z) {
    int n = 1;
    int kode = 2;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy_y, cy_j;

    cy_j.real(NAN);
    cy_j.imag(NAN);
    cy_y.real(NAN);
    cy_y.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy_y;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besy(z, v, kode, n, &cy_y, &ierr);
    set_error_and_nan("yve:", ierr_to_sferr(nz, ierr), cy_y);
    if (ierr == 2) {
        if (z.real() >= 0 && z.imag() == 0) {
            /* overflow */
            cy_y.real(INFINITY);
            cy_y.imag(0);
        }
    }

    if (sign == -1) {
        if (!detail::reflect_jy(&cy_y, v)) {
            nz = amos::besj(z, v, kode, n, &cy_j, &ierr);
            set_error_and_nan("yv(jv):", ierr_to_sferr(nz, ierr), cy_j);
            cy_y = detail::rotate_jy(cy_y, cy_j, -v);
        }
    }
    return cy_y;
}

inline std::complex<float> cyl_bessel_ye(float v, std::complex<float> x) {
    return static_cast<std::complex<float>>(cyl_bessel_ye(static_cast<double>(v), static_cast<std::complex<double>>(x))
    );
}

template <typename T>
T cyl_bessel_ye(T v, T x) {
    if (x < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }

    return std::real(cyl_bessel_ye(v, std::complex(x)));
}

inline std::complex<double> cyl_bessel_ie(double v, std::complex<double> z) {
    int n = 1;
    int kode = 2;
    int sign = 1;
    int nz, ierr;
    std::complex<double> cy, cy_k;

    cy.real(NAN);
    cy.imag(NAN);
    cy_k.real(NAN);
    cy_k.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besi(z, v, kode, n, &cy, &ierr);
    set_error_and_nan("ive:", ierr_to_sferr(nz, ierr), cy);

    if (sign == -1) {
        if (!detail::reflect_i(&cy, v)) {
            nz = amos::besk(z, v, kode, n, &cy_k, &ierr);
            set_error_and_nan("ive(kv):", ierr_to_sferr(nz, ierr), cy_k);
            /* adjust scaling to match zbesi */
            cy_k = detail::rotate(cy_k, -z.imag() / M_PI);
            if (z.real() > 0) {
                cy_k.real(cy_k.real() * exp(-2 * z.real()));
                cy_k.imag(cy_k.imag() * exp(-2 * z.real()));
            }
            /* v -> -v */
            cy = detail::rotate_i(cy, cy_k, v);
        }
    }

    return cy;
}

inline std::complex<float> cyl_bessel_ie(float v, std::complex<float> x) {
    return static_cast<std::complex<float>>(cyl_bessel_ie(static_cast<double>(v), static_cast<std::complex<double>>(x))
    );
}

template <typename T>
T cyl_bessel_ie(T v, T x) {
    if (v != floor(v) && x < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }

    return std::real(cyl_bessel_ie(v, std::complex(x)));
}

inline std::complex<double> cyl_bessel_ke(double v, std::complex<double> z) {
    std::complex<double> cy{NAN, NAN};
    if (isnan(v) || isnan(std::real(z)) || isnan(std::imag(z))) {
        return cy;
    }

    if (v < 0) {
        /* K_v == K_{-v} even for non-integer v */
        v = -v;
    }

    int n = 1;
    int kode = 2;
    int ierr;
    int nz = amos::besk(z, v, kode, n, &cy, &ierr);
    set_error_and_nan("kve:", ierr_to_sferr(nz, ierr), cy);
    if (ierr == 2) {
        if (std::real(z) >= 0 && std::imag(z) == 0) {
            /* overflow */
            cy = INFINITY;
        }
    }

    return cy;
}

inline std::complex<float> cyl_bessel_ke(float v, std::complex<float> x) {
    return static_cast<std::complex<float>>(cyl_bessel_ke(static_cast<double>(v), static_cast<std::complex<double>>(x))
    );
}

template <typename T>
T cyl_bessel_ke(T v, T x) {
    if (x < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (x == 0) {
        return std::numeric_limits<T>::infinity();
    }

    return std::real(cyl_bessel_ke(v, std::complex(x)));
}

inline std::complex<double> cyl_hankel_1e(double v, std::complex<double> z) {
    int n = 1;
    int kode = 2;
    int m = 1;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy;

    cy.real(NAN);
    cy.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besh(z, v, kode, m, n, &cy, &ierr);
    set_error_and_nan("hankel1e:", ierr_to_sferr(nz, ierr), cy);
    if (sign == -1) {
        cy = detail::rotate(cy, v);
    }
    return cy;
}

inline std::complex<float> cyl_hankel_1e(float v, std::complex<float> z) {
    return static_cast<std::complex<float>>(cyl_hankel_1e(static_cast<double>(v), static_cast<std::complex<double>>(z))
    );
}

inline std::complex<double> cyl_hankel_2e(double v, std::complex<double> z) {
    int n = 1;
    int kode = 2;
    int m = 2;
    int nz, ierr;
    int sign = 1;

    std::complex<double> cy{NAN, NAN};
    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }

    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besh(z, v, kode, m, n, &cy, &ierr);
    set_error_and_nan("hankel2e:", ierr_to_sferr(nz, ierr), cy);
    if (sign == -1) {
        cy = detail::rotate(cy, -v);
    }
    return cy;
}

inline std::complex<float> cyl_hankel_2e(float v, std::complex<float> z) {
    return static_cast<std::complex<float>>(cyl_hankel_2e(static_cast<double>(v), static_cast<std::complex<double>>(z))
    );
}

inline std::complex<double> cyl_bessel_j(double v, std::complex<double> z) {
    int n = 1;
    int kode = 1;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy_j, cy_y;

    cy_j.real(NAN);
    cy_j.imag(NAN);
    cy_y.real(NAN);
    cy_y.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy_j;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besj(z, v, kode, n, &cy_j, &ierr);
    set_error_and_nan("jv:", ierr_to_sferr(nz, ierr), cy_j);
    if (ierr == 2) {
        /* overflow */
        cy_j = cyl_bessel_je(v, z);
        cy_j.real(cy_j.real() * INFINITY);
        cy_j.imag(cy_j.imag() * INFINITY);
    }

    if (sign == -1) {
        if (!detail::reflect_jy(&cy_j, v)) {
            nz = amos::besy(z, v, kode, n, &cy_y, &ierr);
            set_error_and_nan("jv(yv):", ierr_to_sferr(nz, ierr), cy_y);
            cy_j = detail::rotate_jy(cy_j, cy_y, v);
        }
    }
    return cy_j;
}

inline std::complex<float> cyl_bessel_j(float v, std::complex<float> x) {
    return static_cast<std::complex<float>>(cyl_bessel_j(static_cast<double>(v), static_cast<std::complex<double>>(x)));
}

template <typename T>
T cyl_bessel_j(T v, T x) {
    if (v != static_cast<int>(v) && x < 0) {
        set_error("jv", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<T>::quiet_NaN();
    }

    std::complex<T> res = cyl_bessel_j(v, std::complex(x));
    if (std::real(res) != std::real(res)) {
        /* AMOS returned NaN, possibly due to overflow */
        return cephes::jv(v, x);
    }

    return std::real(res);
}

inline std::complex<double> cyl_bessel_y(double v, std::complex<double> z) {
    int n = 1;
    int kode = 1;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy_y, cy_j;

    cy_j.real(NAN);
    cy_j.imag(NAN);
    cy_y.real(NAN);
    cy_y.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy_y;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }

    if (z.real() == 0 && z.imag() == 0) {
        /* overflow */
        cy_y.real(-INFINITY);
        cy_y.imag(0);
        set_error("yv", SF_ERROR_OVERFLOW, NULL);
    } else {
        nz = amos::besy(z, v, kode, n, &cy_y, &ierr);
        set_error_and_nan("yv:", ierr_to_sferr(nz, ierr), cy_y);
        if (ierr == 2) {
            if (z.real() >= 0 && z.imag() == 0) {
                /* overflow */
                cy_y.real(-INFINITY);
                cy_y.imag(0);
            }
        }
    }

    if (sign == -1) {
        if (!detail::reflect_jy(&cy_y, v)) {
            nz = amos::besj(z, v, kode, n, &cy_j, &ierr);
            // F_FUNC(zbesj,ZBESJ)(CADDR(z), &v,  &kode, &n, CADDR(cy_j), &nz, &ierr);
            set_error_and_nan("yv(jv):", ierr_to_sferr(nz, ierr), cy_j);
            cy_y = detail::rotate_jy(cy_y, cy_j, -v);
        }
    }
    return cy_y;
}

inline std::complex<float> cyl_bessel_y(float v, std::complex<float> x) {
    return static_cast<std::complex<float>>(cyl_bessel_y(static_cast<double>(v), static_cast<std::complex<double>>(x)));
}

template <typename T>
T cyl_bessel_y(T v, T x) {
    if (x < 0) {
        set_error("yv", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    std::complex<T> res = cyl_bessel_y(v, std::complex(x));
    if (std::real(res) != std::real(res)) {
        return cephes::yv(v, x);
    }

    return std::real(res);
}

inline double cyl_bessel_i(double v, double x) { return cephes::iv(v, x); }

inline float cyl_bessel_i(float v, float x) { return cyl_bessel_i(static_cast<double>(v), static_cast<double>(x)); }

inline std::complex<double> cyl_bessel_i(double v, std::complex<double> z) {
    int n = 1;
    int kode = 1;
    int sign = 1;
    int nz, ierr;
    std::complex<double> cy, cy_k;

    cy.real(NAN);
    cy.imag(NAN);
    cy_k.real(NAN);
    cy_k.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besi(z, v, kode, n, &cy, &ierr);
    set_error_and_nan("iv:", ierr_to_sferr(nz, ierr), cy);
    if (ierr == 2) {
        /* overflow */
        if (z.imag() == 0 && (z.real() >= 0 || v == floor(v))) {
            if (z.real() < 0 && v / 2 != floor(v / 2))
                cy.real(-INFINITY);
            else
                cy.real(INFINITY);
            cy.imag(0);
        } else {
            cy = cyl_bessel_ie(v * sign, z);
            cy.real(cy.real() * INFINITY);
            cy.imag(cy.imag() * INFINITY);
        }
    }

    if (sign == -1) {
        if (!detail::reflect_i(&cy, v)) {
            nz = amos::besk(z, v, kode, n, &cy_k, &ierr);
            set_error_and_nan("iv(kv):", ierr_to_sferr(nz, ierr), cy_k);
            cy = detail::rotate_i(cy, cy_k, v);
        }
    }

    return cy;
}

inline std::complex<float> cyl_bessel_i(float v, std::complex<float> x) {
    return static_cast<std::complex<float>>(cyl_bessel_i(static_cast<double>(v), static_cast<std::complex<double>>(x)));
}

inline std::complex<double> cyl_bessel_k(double v, std::complex<double> z) {
    std::complex<double> cy(NAN, NAN);
    if (std::isnan(v) || std::isnan(std::real(z)) || isnan(std::imag(z))) {
        return cy;
    }

    if (v < 0) {
        /* K_v == K_{-v} even for non-integer v */
        v = -v;
    }

    int n = 1;
    int kode = 1;
    int ierr;
    int nz = amos::besk(z, v, kode, n, &cy, &ierr);
    set_error_and_nan("kv:", ierr_to_sferr(nz, ierr), cy);
    if (ierr == 2) {
        if (std::real(z) >= 0 && std::imag(z) == 0) {
            /* overflow */
            cy = INFINITY;
        }
    }

    return cy;
}

inline std::complex<float> cyl_bessel_k(float v, std::complex<float> x) {
    return static_cast<std::complex<float>>(cyl_bessel_k(static_cast<double>(v), static_cast<std::complex<double>>(x)));
}

template <typename T>
T cyl_bessel_k(T v, T z) {
    if (z < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (z == 0) {
        return std::numeric_limits<T>::infinity();
    }

    if (z > 710 * (1 + std::abs(v))) {
        /* Underflow. See uniform expansion https://dlmf.nist.gov/10.41
         * This condition is not a strict bound (it can underflow earlier),
         * rather, we are here working around a restriction in AMOS.
         */
        return 0;
    }

    return std::real(cyl_bessel_k(v, std::complex(z)));
}

inline std::complex<double> cyl_hankel_1(double v, std::complex<double> z) {
    int n = 1;
    int kode = 1;
    int m = 1;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy;

    cy.real(NAN);
    cy.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besh(z, v, kode, m, n, &cy, &ierr);
    set_error_and_nan("hankel1:", ierr_to_sferr(nz, ierr), cy);
    if (sign == -1) {
        cy = detail::rotate(cy, v);
    }
    return cy;
}

inline std::complex<float> cyl_hankel_1(float v, std::complex<float> z) {
    return static_cast<std::complex<float>>(cyl_hankel_1(static_cast<double>(v), static_cast<std::complex<double>>(z)));
}

inline std::complex<double> cyl_hankel_2(double v, std::complex<double> z) {
    int n = 1;
    int kode = 1;
    int m = 2;
    int nz, ierr;
    int sign = 1;
    std::complex<double> cy;

    cy.real(NAN);
    cy.imag(NAN);

    if (isnan(v) || isnan(z.real()) || isnan(z.imag())) {
        return cy;
    }
    if (v < 0) {
        v = -v;
        sign = -1;
    }
    nz = amos::besh(z, v, kode, m, n, &cy, &ierr);
    set_error_and_nan("hankel2:", ierr_to_sferr(nz, ierr), cy);
    if (sign == -1) {
        cy = detail::rotate(cy, -v);
    }
    return cy;
}

inline std::complex<float> cyl_hankel_2(float v, std::complex<float> z) {
    return static_cast<std::complex<float>>(cyl_hankel_2(static_cast<double>(v), static_cast<std::complex<double>>(z)));
}

} // namespace xsf
