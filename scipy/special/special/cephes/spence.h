/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     spence.c
 *
 *     Dilogarithm
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, spence();
 *
 * y = spence( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the integral
 *
 *                    x
 *                    -
 *                   | | log t
 * spence(x)  =  -   |   ----- dt
 *                 | |   t - 1
 *                  -
 *                  1
 *
 * for x >= 0.  A rational approximation gives the integral in
 * the interval (0.5, 1.5).  Transformation formulas for 1/x
 * and 1-x are employed outside the basic expansion range.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,4         30000       3.9e-15     5.4e-16
 *
 *
 */

/*                                                     spence.c */

/*
 * Cephes Math Library Release 2.1:  January, 1989
 * Copyright 1985, 1987, 1989 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
#pragma once

#include "../config.h"
#include "../error.h"

#include "const.h"
#include "polevl.h"

namespace special {
namespace cephes {

    namespace detail {

        constexpr double spence_A[8] = {
            4.65128586073990045278E-5, 7.31589045238094711071E-3, 1.33847639578309018650E-1, 8.79691311754530315341E-1,
            2.71149851196553469920E0,  4.25697156008121755724E0,  3.29771340985225106936E0,  1.00000000000000000126E0,
        };

        constexpr double spence_B[8] = {
            6.90990488912553276999E-4, 2.54043763932544379113E-2, 2.82974860602568089943E-1, 1.41172597751831069617E0,
            3.63800533345137075418E0,  5.03278880143316990390E0,  3.54771340985225096217E0,  9.99999999999999998740E-1,
        };

    } // namespace detail

    SPECFUN_HOST_DEVICE inline double spence(double x) {
        double w, y, z;
        int flag;

        if (x < 0.0) {
            set_error("spence", SF_ERROR_DOMAIN, NULL);
            return (std::numeric_limits<double>::quiet_NaN());
        }

        if (x == 1.0) {
            return (0.0);
        }

        if (x == 0.0) {
            return (M_PI * M_PI / 6.0);
        }

        flag = 0;

        if (x > 2.0) {
            x = 1.0 / x;
            flag |= 2;
        }

        if (x > 1.5) {
            w = (1.0 / x) - 1.0;
            flag |= 2;
        } else if (x < 0.5) {
            w = -x;
            flag |= 1;
        } else {
            w = x - 1.0;
        }

        y = -w * polevl(w, detail::spence_A, 7) / polevl(w, detail::spence_B, 7);

        if (flag & 1) {
            y = (M_PI * M_PI) / 6.0 - std::log(x) * std::log(1.0 - x) - y;
        }

        if (flag & 2) {
            z = std::log(x);
            y = -0.5 * z * z - y;
        }

        return (y);
    }

} // namespace cephes
} // namespace special
