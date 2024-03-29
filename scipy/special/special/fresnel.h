#pragma once

#include "config.h"

namespace special {
namespace detail {

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

    specfun::cfs(z, zfs, &zfd);
    specfun::cfc(z, zfc, &zfd);
}

template <typename T>
void modified_fresnel_plus(T x, std::complex<T> *Fplus, std::complex<T> *Kplus) {
    detail::ffk(0, x, *Fplus, *Kplus);
}

template <typename T>
void modified_fresnel_minus(T x, std::complex<T> *Fminus, std::complex<T> *Kminus) {
    detail::ffk(1, x, *Fminus, *Kminus);
}

} // namespace special
