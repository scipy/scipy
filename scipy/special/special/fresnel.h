#pragma once

#include "config.h"

namespace special {
namespace detail {

    inline void cfc(std::complex<double> z, std::complex<double> *zf, std::complex<double> *zd) {

        // =========================================================
        // Purpose: Compute complex Fresnel integral C(z) and C'(z)
        // Input :  z --- Argument of C(z)
        // Output:  ZF --- C(z)
        //          ZD --- C'(z)
        // =========================================================

        int k, m;
        double wa0, wa;
        std::complex<double> c, cr, cf, cf0, cf1, cg, d;
        const double eps = 1.0e-14;
        const double pi = 3.141592653589793;

        double w0 = std::abs(z);
        std::complex<double> zp = 0.5 * pi * z * z;
        std::complex<double> zp2 = zp * zp;
        std::complex<double> z0 = 0.0;

        if (z == z0) {
            c = z0;
        } else if (w0 <= 2.5) {
            cr = z;
            c = cr;
            wa0 = 0.0;
            for (k = 1; k <= 80; k++) {
                cr = -0.5 * cr * (4.0 * k - 3.0) / static_cast<double>(k) / (2.0 * k - 1.0) / (4.0 * k + 1.0) * zp2;
                c += cr;
                wa = std::abs(c);
                if ((fabs((wa - wa0) / wa) < eps) && (k > 10)) {
                    *zf = c;
                    *zd = std::cos(0.5 * pi * z * z);
                    return;
                }
                wa0 = wa;
            }
        } else if ((w0 > 2.5) && (w0 < 4.5)) {
            m = 85;
            c = z0;
            cf1 = z0;
            cf0 = 1.0e-100;
            for (k = m; k >= 0; k--) {
                cf = (2.0 * k + 3.0) * cf0 / zp - cf1;
                if (k % 2 == 0) {
                    c += cf;
                }
                cf1 = cf0;
                cf0 = cf;
            }
            c *= 2.0 / (pi * z) * std::sin(zp) / cf;
        } else {
            // See comment at CFS(), use C(z) = iC(-iz)
            if ((z.imag() > -z.real()) && (z.imag() <= z.real())) {
                // right quadrant
                d = 0.5;
            } else if ((z.imag() > z.real()) && (z.imag() >= -z.real())) {
                // upper quadrant
                d = std::complex<double>(0, 0.5);
            } else if ((z.imag() < -z.real()) && (z.imag() >= z.real())) {
                // left quadrant
                d = -0.5;
            } else {
                d = std::complex<double>(0, -0.5);
            }
            cr = 1.0;
            cf = 1.0;
            for (k = 1; k <= 20; k++) {
                cr = -0.25 * cr * (4.0 * k - 1.0) * (4.0 * k - 3.0) / zp2;
                cf += cr;
            }
            cr = 1.0 / (pi * z * z);
            cg = cr;
            for (k = 1; k <= 12; k++) {
                cr = -0.25 * cr * (4.0 * k + 1.0) * (4.0 * k - 1.0) / zp2;
                cg += cr;
            }
            c = d + (cf * std::sin(zp) - cg * std::cos(zp)) / (pi * z);
        }
        *zf = c;
        *zd = std::cos(0.5 * pi * z * z);
        return;
    }

    inline void cfs(std::complex<double> z, std::complex<double> *zf, std::complex<double> *zd) {

        // =========================================================
        // Purpose: Compute complex Fresnel Integral S(z) and S'(z)
        // Input :  z  --- Argument of S(z)
        // Output:  ZF --- S(z)
        //          ZD --- S'(z)
        // =========================================================

        int k, m;
        double wb0, wb;
        std::complex<double> s, cr, cf, cf0, cf1, cg, d;
        const double eps = 1.0e-14;
        const double pi = 3.141592653589793;

        double w0 = std::abs(z);
        std::complex<double> zp = 0.5 * pi * z * z;
        std::complex<double> zp2 = zp * zp;
        std::complex<double> z0 = 0.0;

        if (z == z0) {
            s = z0;
        } else if (w0 <= 2.5) {
            s = z * zp / 3.0;
            cr = s;
            wb0 = 0.0;
            for (k = 1; k <= 80; k++) {
                cr = -0.5 * cr * (4.0 * k - 1.0) / static_cast<double>(k) / (2.0 * k + 1.0) / (4.0 * k + 3.0) * zp2;
                s += cr;
                wb = std::abs(s);
                if ((fabs(wb - wb0) < eps) && (k > 10)) {
                    *zf = s;
                    *zd = std::sin(0.5 * pi * z * z);
                    return;
                }
                wb0 = wb;
            }
        } else if ((w0 > 2.5) && (w0 < 4.5)) {
            m = 85;
            s = z0;
            cf1 = z0;
            cf0 = 1.0e-100;
            for (k = m; k >= 0; k--) {
                cf = (2.0 * k + 3.0) * cf0 / zp - cf1;
                if (k % 2 == 1) {
                    s += cf;
                }
                cf1 = cf0;
                cf0 = cf;
            }
            s = 2.0 / (pi * z) * std::sin(zp) / cf * s;
        } else {
            // Auxiliary functions f(z) and g(z) can be computed using an
            // asymptotic expansion in the right quadrant |arg(z)| <= pi/4, not pi/2
            // as sometimes suggested. Use the symmetry S(z) = -iS(-iz).
            // Interestingly, most of the expansion code is the same across
            // the quadrants. (The forth power in Z is the equalizer here.)
            // Only one constant has to be adapted.
            if ((z.imag() > -z.real()) && (z.imag() <= z.real())) {
                // right quadrant
                d = 0.5;
            } else if ((z.imag() > z.real()) && (z.imag() >= -z.real())) {
                // upper quadrant
                d = std::complex<double>(0, -0.5);
            } else if ((z.imag() < -z.real()) && (z.imag() >= z.real())) {
                // left quadrant
                d = -0.5;
            } else {
                d = std::complex<double>(0, 0.5);
            }
            cr = 1.0;
            cf = 1.0;
            for (k = 1; k <= 20; k++) {
                cr = -0.25 * cr * (4.0 * k - 1.0) * (4.0 * k - 3.0) / zp2;
                cf += cr;
            }
            cr = 1.0;
            cg = 1.0;
            for (k = 1; k <= 12; k++) {
                cr = -0.25 * cr * (4.0 * k + 1.0) * (4.0 * k - 1.0) / zp2;
                cg += cr;
            }
            cg = cg / (pi * z * z);
            s = d - (cf * std::cos(zp) + cg * std::sin(zp)) / (pi * z);
        }
        *zf = s;
        *zd = std::sin(0.5 * pi * z * z);
        return;
    }

    template <typename T>
    void ffk(int ks, T x, std::complex<T> &f, std::complex<T> &g) {

        // =======================================================
        // Purpose: Compute modified Fresnel integrals F±(x)
        //          and K±(x)
        // Input :  x   --- Argument of F±(x) and K±(x)
        //          KS  --- Sign code
        //                  KS=0 for calculating F+(x) and K+(x)
        //                  KS=1 for calculating F_(x) and K_(x)
        // Output:  FR  --- Re[F±(x)]
        //          FI  --- Im[F±(x)]
        //          GR  --- Re[K±(x)]
        //          GI  --- Im[K±(x)]
        // ======================================================

        const T eps = 1.0e-15;
        const T pi = 3.141592653589793;
        const T pp2 = 1.2533141373155;
        const T p2p = 0.7978845608028654;

        T fi0, c1, s1, cs, ss, xa, x2, x4, xc, xf, xf0, xf1, xg, xp, xq, xq2, xr, xs, xsu, xw;

        xa = fabs(x);
        x2 = x * x;
        x4 = x2 * x2;

        if (x == 0.0) {
            f.real(0.5 * sqrt(0.5 * pi));
            f.imag(pow(-1, ks) * f.real());
            g = 0.5;
        } else {
            if (xa <= 2.5) {
                xr = p2p * xa;
                c1 = xr;

                for (int k = 1; k <= 50; ++k) {
                    xr = -0.5 * xr * (4.0 * k - 3.0) / k / (2.0 * k - 1.0) / (4.0 * k + 1.0) * x4;
                    c1 += xr;
                    if (fabs(xr / c1) < eps)
                        break;
                }

                s1 = p2p * xa * xa * xa / 3.0;
                xr = s1;

                for (int k = 1; k <= 50; ++k) {
                    xr = -0.5 * xr * (4.0 * k - 1.0) / k / (2.0 * k + 1.0) / (4.0 * k + 3.0) * x4;
                    s1 += xr;
                    if (fabs(xr / s1) < eps)
                        break;
                }

                f.real(pp2 * (0.5 - c1));
                fi0 = pp2 * (0.5 - s1);
                f.imag(pow(-1, ks) * fi0);
            } else if (xa < 5.5) {
                int m = (int) (42 + 1.75 * x2);
                xsu = 0.0;
                xc = 0.0;
                xs = 0.0;
                xf1 = 0.0;
                xf0 = 1.0e-100;

                for (int k = m; k >= 0; --k) {
                    xf = (2.0 * k + 3.0) * xf0 / x2 - xf1;
                    if (k % 2 == 0) {
                        xc += xf;
                    } else {
                        xs += xf;
                    }
                    xsu += (2.0 * k + 1.0) * xf * xf;
                    xf1 = xf0;
                    xf0 = xf;
                }

                xq = sqrt(xsu);
                xw = p2p * xa / xq;
                c1 = xc * xw;
                s1 = xs * xw;
            } else {
                xr = 1.0;
                xf = 1.0;

                for (int k = 1; k <= 12; ++k) {
                    xr = -0.25 * xr * (4.0 * k - 1.0) * (4.0 * k - 3.0) / x4;
                    xf += xr;
                }

                xr = 1.0 / (2.0 * xa * xa);
                xg = xr;

                for (int k = 1; k <= 12; ++k) {
                    xr = -0.25 * xr * (4.0 * k + 1.0) * (4.0 * k - 1.0) / x4;
                    xg += xr;
                }

                c1 = 0.5 + (xf * sin(x2) - xg * cos(x2)) / sqrt(2.0 * pi) / xa;
                s1 = 0.5 - (xf * cos(x2) + xg * sin(x2)) / sqrt(2.0 * pi) / xa;
            }

            f.real(pp2 * (0.5 - c1));
            fi0 = pp2 * (0.5 - s1);
            f.imag(pow(-1, ks) * fi0);

            xp = x2 + pi / 4.0;
            cs = cos(xp);
            ss = sin(xp);
            xq2 = 1.0 / sqrt(pi);

            g.real(xq2 * (f.real() * cs + fi0 * ss));
            g.imag(pow(-1, ks) * xq2 * (fi0 * cs - f.real() * ss));

            if (x < 0.0) {
                f.real(pp2 - f.real());
                f.imag(pow(-1, ks) * pp2 - f.real());
                g.real(cos(x2) - g.real());
                g.imag(-pow(-1, ks) * sin(x2) - g.imag());
            }
        }
    }

} // namespace detail

/* Fresnel integrals of complex numbers */

inline void cfresnl(std::complex<double> z, std::complex<double> *zfs, std::complex<double> *zfc) {
    std::complex<double> zfd;

    detail::cfs(z, zfs, &zfd);
    detail::cfc(z, zfc, &zfd);
}

template <typename T>
void modified_fresnel_plus(T x, std::complex<T> &Fplus, std::complex<T> &Kplus) {
    detail::ffk(0, x, Fplus, Kplus);
}

template <typename T>
void modified_fresnel_minus(T x, std::complex<T> &Fminus, std::complex<T> &Kminus) {
    detail::ffk(1, x, Fminus, Kminus);
}

inline void fcszo(int kf, int nt, std::complex<double> *zo) {

    // ===============================================================
    // Purpose: Compute the complex zeros of Fresnel integral C(z)
    //          or S(z) using modified Newton's iteration method
    // Input :  KF  --- Function code
    //                  KF=1 for C(z) or KF=2 for S(z)
    //          NT  --- Total number of zeros
    // Output:  ZO(L) --- L-th zero of C(z) or S(z)
    // Routines called:
    //      (1) CFC for computing Fresnel integral C(z)
    //      (2) CFS for computing Fresnel integral S(z)
    // ==============================================================

    int it;
    double psq, px, py, w, w0;
    std::complex<double> z, zp, zf, zd, zfd, zgd, zq, zw;
    const double pi = 3.141592653589793;
    psq = 0.0;
    w = 0.0;

    for (int nr = 1; nr <= nt; ++nr) {
        if (kf == 1)
            psq = sqrt(4.0 * nr - 1.0);
        if (kf == 2)
            psq = 2.0 * sqrt(nr);

        px = psq - log(pi * psq) / (pi * pi * psq * psq * psq);
        py = log(pi * psq) / (pi * psq);
        z = std::complex<double>(px, py);

        if (kf == 2) {
            if (nr == 2) {
                z = std::complex<double>(2.8334, 0.2443);
            }
            if (nr == 3) {
                z = std::complex<double>(3.4674, 0.2185);
            }
            if (nr == 4) {
                z = std::complex<double>(4.0025, 0.2008);
            }
        }

        it = 0;
        do {
            it++;
            if (kf == 1) {
                detail::cfc(z, &zf, &zd);
            }
            if (kf == 2) {
                detail::cfs(z, &zf, &zd);
            }

            zp = 1.0;
            for (int i = 1; i < nr; i++)
                zp *= (z - zo[i - 1]);

            zfd = zf / zp;
            zq = 0.0;
            for (int i = 1; i < nr; i++) {
                zw = 1.0;
                for (int j = 1; j < nr; j++) {
                    if (j == i) {
                        continue;
                    }
                    zw *= (z - zo[j - 1]);
                }
                zq += zw;
            }
            zgd = (zd - zq * zfd) / zp;
            z -= zfd / zgd;
            w0 = w;
            w = std::abs(z);
        } while ((it <= 50) && (fabs((w - w0) / w) > 1.0e-12));
        zo[nr - 1] = z;
    }
    return;
}

} // namespace special
