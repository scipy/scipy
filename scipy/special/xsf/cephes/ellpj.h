/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     ellpj.c
 *
 *     Jacobian Elliptic Functions
 *
 *
 *
 * SYNOPSIS:
 *
 * double u, m, sn, cn, dn, phi;
 * int ellpj();
 *
 * ellpj( u, m, _&sn, _&cn, _&dn, _&phi );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 * Evaluates the Jacobian elliptic functions sn(u|m), cn(u|m),
 * and dn(u|m) of parameter m between 0 and 1, and real
 * argument u.
 *
 * These functions are periodic, with quarter-period on the
 * real axis equal to the complete elliptic integral
 * ellpk(m).
 *
 * Relation to incomplete elliptic integral:
 * If u = ellik(phi,m), then sn(u|m) = sin(phi),
 * and cn(u|m) = cos(phi).  Phi is called the amplitude of u.
 *
 * Computation is by means of the arithmetic-geometric mean
 * algorithm, except when m is within 1e-9 of 0 or 1.  In the
 * latter case with m close to 1, the approximation applies
 * only for phi < pi/2.
 *
 * ACCURACY:
 *
 * Tested at random points with u between 0 and 10, m between
 * 0 and 1.
 *
 *            Absolute error (* = relative error):
 * arithmetic   function   # trials      peak         rms
 *    IEEE      phi         10000       9.2e-16*    1.4e-16*
 *    IEEE      sn          50000       4.1e-15     4.6e-16
 *    IEEE      cn          40000       3.6e-15     4.4e-16
 *    IEEE      dn          10000       1.3e-12     1.8e-14
 *
 *  Peak error observed in consistency check using addition
 * theorem for sn(u+v) was 4e-16 (absolute).  Also tested by
 * the above relation to the incomplete elliptic integral.
 * Accuracy deteriorates when u is large.
 *
 */

/*                                                     ellpj.c         */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

/* Scipy changes:
 * - 07-18-2016: improve evaluation of dn near quarter periods
 */
#pragma once

#include "../config.h"
#include "../error.h"
#include "const.h"

namespace xsf {
namespace cephes {

    XSF_HOST_DEVICE inline int ellpj(double u, double m, double *sn, double *cn, double *dn, double *ph) {
        double ai, b, phi, t, twon, dnfac;
        double a[9], c[9];
        int i;

        /* Check for special cases */
        if (m < 0.0 || m > 1.0 || std::isnan(m)) {
            set_error("ellpj", SF_ERROR_DOMAIN, NULL);
            *sn = std::numeric_limits<double>::quiet_NaN();
            *cn = std::numeric_limits<double>::quiet_NaN();
            *ph = std::numeric_limits<double>::quiet_NaN();
            *dn = std::numeric_limits<double>::quiet_NaN();
            return (-1);
        }
        if (m < 1.0e-9) {
            t = std::sin(u);
            b = std::cos(u);
            ai = 0.25 * m * (u - t * b);
            *sn = t - ai * b;
            *cn = b + ai * t;
            *ph = u - ai;
            *dn = 1.0 - 0.5 * m * t * t;
            return (0);
        }
        if (m >= 0.9999999999) {
            ai = 0.25 * (1.0 - m);
            b = std::cosh(u);
            t = std::tanh(u);
            phi = 1.0 / b;
            twon = b * std::sinh(u);
            *sn = t + ai * (twon - u) / (b * b);
            *ph = 2.0 * std::atan(exp(u)) - M_PI_2 + ai * (twon - u) / b;
            ai *= t * phi;
            *cn = phi - ai * (twon - u);
            *dn = phi + ai * (twon + u);
            return (0);
        }

        /* A. G. M. scale. See DLMF 22.20(ii) */
        a[0] = 1.0;
        b = std::sqrt(1.0 - m);
        c[0] = std::sqrt(m);
        twon = 1.0;
        i = 0;

        while (std::abs(c[i] / a[i]) > detail::MACHEP) {
            if (i > 7) {
                set_error("ellpj", SF_ERROR_OVERFLOW, NULL);
                goto done;
            }
            ai = a[i];
            ++i;
            c[i] = (ai - b) / 2.0;
            t = std::sqrt(ai * b);
            a[i] = (ai + b) / 2.0;
            b = t;
            twon *= 2.0;
        }

    done:
        /* backward recurrence */
        phi = twon * a[i] * u;
        do {
            t = c[i] * std::sin(phi) / a[i];
            b = phi;
            phi = (std::asin(t) + phi) / 2.0;
        } while (--i);

        *sn = std::sin(phi);
        t = std::cos(phi);
        *cn = t;
        dnfac = std::cos(phi - b);
        /* See discussion after DLMF 22.20.5 */
        if (std::abs(dnfac) < 0.1) {
            *dn = std::sqrt(1 - m * (*sn) * (*sn));
        } else {
            *dn = t / dnfac;
        }
        *ph = phi;
        return (0);
    }

} // namespace cephes
} // namespace xsf
