#pragma once

#include "specfun/specfun.h"

namespace xsf {
namespace detail {

    template <typename T>
    T vvla(T x, T va);

    template <typename T>
    T dvsa(T x, T va) {

        // ===================================================
        // Purpose: Compute parabolic cylinder function Dv(x)
        //          for small argument
        // Input:   x  --- Argument
        //          va --- Order
        // Output:  PD --- Dv(x)
        // Routine called: GAMMA2 for computing Г(x)
        // ===================================================

        int m;
        T ep, a0, va0, pd, ga0, g1, g0, r, vm, gm, r1, vt;

        const T pi = 3.141592653589793;
        const T eps = 1.0e-15;
        const T sq2 = sqrt(2);

        ep = exp(-0.25 * x * x);
        va0 = 0.5 * (1.0 - va);
        if (va == 0.0) {
            pd = ep;
        } else {
            if (x == 0.0) {
                if ((va0 <= 0.0) && (va0 == (int) va0)) {
                    pd = 0.0;
                } else {
                    ga0 = specfun::gamma2(va0);
                    pd = sqrt(pi) / (pow(2.0, -0.5 * va) * ga0);
                }
            } else {
                g1 = specfun::gamma2(-va);
                a0 = pow(2.0, -0.5 * va - 1.0) * ep / g1;
                vt = -0.5 * va;
                g0 = specfun::gamma2(vt);
                pd = g0;
                r = 1.0;
                for (m = 1; m <= 250; m++) {
                    vm = 0.5 * (m - va);
                    gm = specfun::gamma2(vm);
                    r = -r * sq2 * x / m;
                    r1 = gm * r;
                    pd += r1;
                    if (fabs(r1) < fabs(pd) * eps) {
                        break;
                    }
                }
                pd *= a0;
            }
        }
        return pd;
    }

    template <typename T>
    T dvla(T x, T va) {

        // ====================================================
        // Purpose: Compute parabolic cylinder functions Dv(x)
        //          for large argument
        // Input:   x  --- Argument
        //          va --- Order
        // Output:  PD --- Dv(x)
        // Routines called:
        //       (1) VVLA for computing Vv(x) for large |x|
        //       (2) GAMMA2 for computing Г(x)
        // ====================================================

        int k;
        T ep, a0, gl, pd, r, vl, x1;

        const T pi = 3.141592653589793;
        const T eps = 1.0e-12;
        ep = exp(-.25 * x * x);
        a0 = pow(fabs(x), va) * ep;
        r = 1.0;
        pd = 1.0;
        for (k = 1; k <= 16; k++) {
            r = -0.5 * r * (2.0 * k - va - 1.0) * (2.0 * k - va - 2.0) / (k * x * x);
            pd += r;
            if (fabs(r / pd) < eps) {
                break;
            }
        }
        pd *= a0;
        if (x < 0.0) {
            x1 = -x;
            vl = vvla(x1, va);
            gl = specfun::gamma2(-va);
            pd = pi * vl / gl + cos(pi * va) * pd;
        }
        return pd;
    }

    template <typename T>
    T vvla(T x, T va) {

        // ===================================================
        // Purpose: Compute parabolic cylinder function Vv(x)
        //          for large argument
        // Input:   x  --- Argument
        //          va --- Order
        // Output:  PV --- Vv(x)
        // Routines called:
        //       (1) DVLA for computing Dv(x) for large |x|
        //       (2) GAMMA2 for computing Г(x)
        // ===================================================

        int k;
        T pv, qe, a0, r, x1, gl, dsl, pdl;

        const T pi = 3.141592653589793;
        const T eps = 1e-12;

        qe = exp(0.25 * x * x);
        a0 = pow(fabs(x), -va - 1) * sqrt(2.0 / pi) * qe;
        r = 1.0;
        pv = 1.0;
        for (k = 1; k <= 18; k++) {
            r = 0.5 * r * (2.0 * k + va - 1.0) * (2.0 * k + va) / (k * x * x);
            pv += r;
            if (fabs(r / pv) < eps) {
                break;
            }
        }
        pv *= a0;
        if (x < 0.0) {
            x1 = -x;
            pdl = dvla(x1, va);
            gl = specfun::gamma2(-va);
            dsl = sin(pi * va) * sin(pi * va);
            pv = dsl * gl / pi * pdl - cos(pi * va) * pv;
        }
        return pv;
    }

    template <typename T>
    T vvsa(T x, T va) {

        // ===================================================
        // Purpose: Compute parabolic cylinder function Vv(x)
        //          for small argument
        // Input:   x  --- Argument
        //          va --- Order
        // Output:  PV --- Vv(x)
        // Routine called : GAMMA2 for computing Г(x)
        // ===================================================

        T a0, fac, g1, ga0, gm, gw, r, r1, sq2, sv, sv0, v1, vb0, vm, pv;

        const T eps = 1.0e-15;
        const T pi = 3.141592653589793;
        const T ep = exp(-0.25 * x * x);
        T va0 = 1.0 + 0.5 * va;

        if (x == 0.0) {
            if (((va0 <= 0.0) && (va0 == (int) va0)) || va == 0.0) {
                pv = 0.0;
            } else {
                vb0 = -0.5 * va;
                sv0 = sin(va0 * pi);
                ga0 = specfun::gamma2(va0);
                pv = pow(2.0, vb0) * sv0 / ga0;
            }
        } else {
            sq2 = sqrt(2.0);
            a0 = pow(2.0, -0.5 * va) * ep / (2.0 * pi);
            sv = sin(-(va + 0.5) * pi);
            v1 = -0.5 * va;
            g1 = specfun::gamma2(v1);
            pv = (sv + 1.0) * g1;
            r = 1.0;
            fac = 1.0;

            for (int m = 1; m <= 250; m++) {
                vm = 0.5 * (m - va);
                gm = specfun::gamma2(vm);
                r = r * sq2 * x / m;
                fac = -fac;
                gw = fac * sv + 1.0;
                r1 = gw * r * gm;
                pv += r1;
                if ((fabs(r1 / pv) < eps) && (gw != 0.0)) {
                    break;
                }
            }
            pv *= a0;
        }
        return pv;
    }

    template <typename T>
    void pbdv(T x, T v, T *dv, T *dp, T *pdf, T *pdd) {

        // ====================================================
        // Purpose: Compute parabolic cylinder functions Dv(x)
        //          and their derivatives
        // Input:   x --- Argument of Dv(x)
        //          v --- Order of Dv(x)
        // Output:  DV(na) --- Dn+v0(x)
        //          DP(na) --- Dn+v0'(x)
        //          ( na = |n|, v0 = v-n, |v0| < 1,
        //            n = 0,±1,±2,… )
        //          PDF --- Dv(x)
        //          PDD --- Dv'(x)
        // Routines called:
        //       (1) DVSA for computing Dv(x) for small |x|
        //       (2) DVLA for computing Dv(x) for large |x|
        // ====================================================

        int ja, k, l, m, nk, nv, na;
        T xa, vh, ep, f, f0, f1, v0, v1, v2, pd, pd0, pd1, s0;

        xa = fabs(x);
        vh = v;
        v += copysign(1.0, v);
        nv = (int) v;
        v0 = v - nv;
        na = abs(nv);
        ep = exp(-0.25 * x * x);
        ja = 0;
        if (na >= 1) {
            ja = 1;
        }
        if (v >= 0.0) {
            if (v0 == 0.0) {
                pd0 = ep;
                pd1 = x * ep;
            } else {
                for (l = 0; l <= ja; l++) {
                    v1 = v0 + l;
                    if (xa <= 5.8) {
                        pd1 = dvsa(x, v1);
                    } else {
                        pd1 = dvla(x, v1);
                    }
                    if (l == 0) {
                        pd0 = pd1;
                    }
                }
            }
            dv[0] = pd0;
            dv[1] = pd1;
            for (k = 2; k <= na; k++) {
                *pdf = x * pd1 - (k + v0 - 1.0) * pd0;
                dv[k] = *pdf;
                pd0 = pd1;
                pd1 = *pdf;
            }
        } else {
            if (x <= 0.0) {
                if (xa <= 5.8) {
                    pd0 = dvsa(x, v0);
                    v1 = v0 - 1.0;
                    pd1 = dvsa(x, v1);
                } else {
                    pd0 = dvla(x, v0);
                    v1 = v0 - 1.0;
                    pd1 = dvla(x, v1);
                }
                dv[0] = pd0;
                dv[1] = pd1;
                for (k = 2; k <= na; k++) {
                    pd = (-x * pd1 + pd0) / (k - 1.0 - v0);
                    dv[k] = pd;
                    pd0 = pd1;
                    pd1 = pd;
                }
            } else if (x <= 2.0) {
                v2 = nv + v0;
                if (nv == 0) {
                    v2 -= 1.0;
                }
                nk = (int) (-v2);
                f1 = dvsa(x, v2);
                v1 = v2 + 1.0;
                f0 = dvsa(x, v1);
                dv[nk] = f1;
                dv[nk - 1] = f0;
                for (k = nk - 2; k >= 0; k--) {
                    f = x * f0 + (k - v0 + 1.0) * f1;
                    dv[k] = f;
                    f1 = f0;
                    f0 = f;
                }
            } else {
                if (xa <= 5.8) {
                    pd0 = dvsa(x, v0);
                } else {
                    pd0 = dvla(x, v0);
                }
                dv[0] = pd0;
                m = 100 + na;
                f1 = 0.0;
                f0 = 1e-30;
                f = 0.0;
                for (k = m; k >= 0; k--) {
                    f = x * f0 + (k - v0 + 1.0) * f1;
                    if (k <= na) {
                        dv[k] = f;
                    }
                    f1 = f0;
                    f0 = f;
                }
                s0 = pd0 / f;
                for (k = 0; k <= na; k++) {
                    dv[k] *= s0;
                }
            }
        }
        for (k = 0; k < na; k++) {
            v1 = fabs(v0) + k;
            if (v >= 0.0) {
                dp[k] = 0.5 * x * dv[k] - dv[k + 1];
            } else {
                dp[k] = -0.5 * x * dv[k] - v1 * dv[k + 1];
            }
        }
        *pdf = dv[na - 1];
        *pdd = dp[na - 1];
        v = vh;
        return;
    }

    template <typename T>
    void pbvv(T x, T v, T *vv, T *vp, T *pvf, T *pvd) {

        // ===================================================
        // Purpose: Compute parabolic cylinder functions Vv(x)
        //          and their derivatives
        // Input:   x --- Argument of Vv(x)
        //          v --- Order of Vv(x)
        // Output:  VV(na) --- Vv(x)
        //          VP(na) --- Vv'(x)
        //          ( na = |n|, v = n+v0, |v0| < 1
        //            n = 0,±1,±2,… )
        //          PVF --- Vv(x)
        //          PVD --- Vv'(x)
        // Routines called:
        //       (1) VVSA for computing Vv(x) for small |x|
        //       (2) VVLA for computing Vv(x) for large |x|
        // ===================================================

        int ja, k, kv, l, m, na, nv;
        T f, f0, f1, pv0, q2p, qe, s0, v0, v1, v2, vh, xa;

        const T pi = 3.141592653589793;

        xa = fabs(x);
        vh = v;
        v += copysign(1.0, v);
        nv = (int) v;
        v0 = v - nv;
        na = abs(nv);
        qe = exp(0.25 * x * x);
        q2p = sqrt(2.0 / pi);
        ja = 0;
        if (na >= 1) {
            ja = 1;
        }
        f = 0.0;
        if (v <= 0.0) {
            if (v0 == 0.0) {
                if (xa <= 7.5) {
                    pv0 = vvsa(x, v0);
                } else {
                    pv0 = vvla(x, v0);
                }
                f0 = q2p * qe;
                f1 = x * f0;
                vv[0] = pv0;
                vv[1] = f0;
                vv[2] = f1;
            } else {
                for (l = 0; l <= ja; l++) {
                    v1 = v0 - l;
                    if (xa <= 7.5) {
                        f1 = vvsa(x, v1);
                    } else {
                        f1 = vvla(x, v1);
                    }
                    if (l == 0) {
                        f0 = f1;
                    }
                }
                vv[0] = f0;
                vv[1] = f1;
            }
            kv = 2;
            if (v0 == 0.0) {
                kv = 3;
            }
            for (k = kv; k <= na; k++) {
                f = x * f1 + (k - v0 - 2.0) * f0;
                vv[k] = f;
                f0 = f1;
                f1 = f;
            }
        } else {
            if ((x >= 0.0) && (x <= 7.5)) {
                v2 = v;
                if (v2 < 1.0) {
                    v2 = v2 + 1.0;
                }
                f1 = vvsa(x, v2);
                v1 = v2 - 1.0;
                kv = (int) v2;
                f0 = vvsa(x, v1);
                vv[kv] = f1;
                vv[kv - 1] = f0;
                for (k = kv - 2; k >= 0; k--) {
                    f = x * f0 - (k + v0 + 2.0) * f1;
                    if (k <= na) {
                        vv[k] = f;
                    }
                    f1 = f0;
                    f0 = f;
                }
            } else if (x > 7.5) {
                pv0 = vvla(x, v0);
                m = 100 + abs(na);
                vv[1] = pv0;
                f1 = 0.0;
                f0 = 1.0e-40;
                for (k = m; k >= 0; k--) {
                    f = x * f0 - (k + v0 + 2.0) * f1;
                    if (k <= na) {
                        vv[k] = f;
                    }
                    f1 = f0;
                    f0 = f;
                }
                s0 = pv0 / f;
                for (k = 0; k <= na; k++) {
                    vv[k] *= s0;
                }
            } else {
                if (xa <= 7.5) {
                    f0 = vvsa(x, v0);
                    v1 = v0 + 1.0;
                    f1 = vvsa(x, v1);
                } else {
                    f0 = vvla(x, v0);
                    v1 = v0 + 1.0;
                    f1 = vvla(x, v1);
                }
                vv[0] = f0;
                vv[1] = f1;
                for (k = 2; k <= na; k++) {
                    f = (x * f1 - f0) / (k + v0);
                    vv[k] = f;
                    f0 = f1;
                    f1 = f;
                }
            }
        }
        for (k = 0; k < na; k++) {
            v1 = v0 + k;
            if (v >= 0.0) {
                vp[k] = 0.5 * x * vv[k] - (v1 + 1.0) * vv[k + 1];
            } else {
                vp[k] = -0.5 * x * vv[k] + vv[k + 1];
            }
        }
        *pvf = vv[na - 1];
        *pvd = vp[na - 1];
        v = vh;
        return;
    }

    template <typename T>
    void pbwa(T a, T x, T *w1f, T *w1d, T *w2f, T *w2d) {

        // ======================================================
        // Purpose: Compute parabolic cylinder functions W(a,±x)
        //          and their derivatives
        // Input  : a --- Parameter  ( 0 ≤ |a| ≤ 5 )
        //          x --- Argument of W(a,±x)  ( 0 ≤ |x| ≤ 5 )
        // Output : W1F --- W(a,x)
        //          W1D --- W'(a,x)
        //          W2F --- W(a,-x)
        //          W2D --- W'(a,-x)
        // Routine called:
        //         CGAMA for computing complex gamma function
        // ======================================================

        int k, L1, L2;
        T d[80], d1, d2, dl, f1, f2, g1, g2, h[100], h0, h1, hl, r, r1, y1d, y2d, y1f, y2f;
        std::complex<T> ug, vg;
        const T eps = 1e-15;
        const T p0 = 0.59460355750136;

        if (a == 0.0) {
            g1 = 3.625609908222;
            g2 = 1.225416702465;
        } else {
            ug = specfun::cgama(std::complex<T>(0.25, 0.5 * a), 1);
            g1 = std::abs(ug);
            vg = specfun::cgama(std::complex<T>(0.75, 0.5 * a), 1);
            g2 = std::abs(vg);
        }
        f1 = sqrt(g1 / g2);
        f2 = sqrt(2.0 * g2 / g1);
        h0 = 1.0;
        h1 = a;
        h[0] = a;
        for (L1 = 2; L1 <= 100; L1++) {
            hl = a * h1 - 0.25 * (2 * L1 - 2.0) * (2 * L1 - 3.0) * h0;
            h[L1 - 1] = hl;
            h0 = h1;
            h1 = hl;
        }
        y1f = 1.0;
        r = 1.0;
        for (k = 1; k <= 100; k++) {
            r = 0.5 * r * x * x / (k * (2.0 * k - 1.0));
            r1 = h[k - 1] * r;
            y1f += r1;
            if ((fabs(r1) <= eps * fabs(y1f)) && (k > 30)) {
                break;
            }
        }
        y1d = a;
        r = 1.0;
        for (k = 1; k < 100; k++) {
            r = 0.5 * r * x * x / (k * (2.0 * k + 1.0));
            r1 = h[k] * r;
            y1d += r1;
            if ((fabs(r1) <= eps * fabs(y1d)) && (k > 30)) {
                break;
            }
        }
        y1d *= x;
        d1 = 1.0;
        d2 = a;
        d[0] = 1.0;
        d[1] = a;
        for (L2 = 3; L2 <= 80; L2++) {
            dl = a * d2 - 0.25 * ((2 * L2 - 1) - 2.0) * ((2 * L2 - 1) - 3.0) * d1;
            d[L2 - 1] = dl;
            d1 = d2;
            d2 = dl;
        }
        y2f = 1.0;
        r = 1.0;
        for (k = 1; k < 80; k++) {
            r = 0.5 * r * x * x / (k * (2.0 * k + 1.0));
            r1 = d[k] * r;
            y2f += r1;
            if ((fabs(r1) <= eps * fabs(y2f)) && (k > 30)) {
                break;
            }
        }
        y2f *= x;
        y2d = 1.0;
        r = 1.0;
        for (k = 1; k < 80; k++) {
            r = 0.5 * r * x * x / (k * (2.0 * k - 1.0));
            r1 = d[k] * r;
            y2d += r1;
            if ((fabs(r1) <= eps * fabs(y2d)) && (k > 30)) {
                break;
            }
        }
        *w1f = p0 * (f1 * y1f - f2 * y2f);
        *w2f = p0 * (f1 * y1f + f2 * y2f);
        *w1d = p0 * (f1 * y1d - f2 * y2d);
        *w2d = p0 * (f1 * y1d + f2 * y2d);
        return;
    }

} // namespace detail

/*
 * If x > 0 return w1f and w1d. Otherwise set x = abs(x) and return
 * w2f and -w2d.
 */
template <typename T>
void pbwa(T a, T x, T &wf, T &wd) {
    int flag = 0;
    T w1f = 0.0, w1d = 0.0, w2f = 0.0, w2d = 0.0;

    if (x < -5 || x > 5 || a < -5 || a > 5) {
        /*
         * The Zhang and Jin implementation only uses Taylor series;
         * return NaN outside of the range which they are accurate.
         */
        wf = std::numeric_limits<T>::quiet_NaN();
        wd = std::numeric_limits<T>::quiet_NaN();
        set_error("pbwa", SF_ERROR_LOSS, NULL);
    } else {
        if (x < 0) {
            x = -x;
            flag = 1;
        }
        detail::pbwa(a, x, &w1f, &w1d, &w2f, &w2d);
        if (flag) {
            wf = w2f;
            wd = -w2d;
        } else {
            wf = w1f;
            wd = w1d;
        }
    }
}

template <typename T>
void pbdv(T v, T x, T &pdf, T &pdd) {
    T *dv;
    T *dp;
    int num;

    if (isnan(v) || isnan(x)) {
        pdf = std::numeric_limits<T>::quiet_NaN();
        pdd = std::numeric_limits<T>::quiet_NaN();
    } else {
        /* NB. Indexing of DV/DP in specfun.f:PBDV starts from 0, hence +2 */
        num = std::abs((int) v) + 2;
        dv = (T *) malloc(sizeof(T) * 2 * num);
        if (dv == NULL) {
            set_error("pbdv", SF_ERROR_OTHER, "memory allocation error");
            pdf = std::numeric_limits<T>::quiet_NaN();
            pdd = std::numeric_limits<T>::quiet_NaN();
        } else {
            dp = dv + num;
            detail::pbdv(x, v, dv, dp, &pdf, &pdd);
            free(dv);
        }
    }
}

template <typename T>
void pbvv(T v, T x, T &pvf, T &pvd) {
    T *vv;
    T *vp;
    int num;

    if (isnan(v) || isnan(x)) {
        pvf = std::numeric_limits<T>::quiet_NaN();
        pvd = std::numeric_limits<T>::quiet_NaN();
    } else {
        /* NB. Indexing of DV/DP in specfun.f:PBVV starts from 0, hence +2 */
        num = std::abs((int) v) + 2;
        vv = (T *) malloc(sizeof(T) * 2 * num);
        if (vv == NULL) {
            set_error("pbvv", SF_ERROR_OTHER, "memory allocation error");
            pvf = std::numeric_limits<T>::quiet_NaN();
            pvd = std::numeric_limits<T>::quiet_NaN();
        } else {
            vp = vv + num;
            detail::pbvv(x, v, vv, vp, &pvf, &pvd);
            free(vv);
        }
    }
}

} // namespace xsf
