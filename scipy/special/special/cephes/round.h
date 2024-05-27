/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     round.c
 *
 *     Round double to nearest or even integer valued double
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, round();
 *
 * y = round(x);
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the nearest integer to x as a double precision
 * floating point result.  If x ends in 0.5 exactly, the
 * nearest even integer is chosen.
 *
 *
 *
 * ACCURACY:
 *
 * If x is greater than 1/(2*MACHEP), its closest machine
 * representation is already an integer, so rounding does
 * not change it.
 */

/*
 * Cephes Math Library Release 2.1:  January, 1989
 * Copyright 1984, 1987, 1989 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
#pragma once

#include "../config.h"

namespace special {
namespace cephes {

    double round(double x) {
        double y, r;

        /* Largest integer <= x */
        y = std::floor(x);

        /* Fractional part */
        r = x - y;

        /* Round up to nearest. */
        if (r > 0.5) {
            goto rndup;
        }

        /* Round to even */
        if (r == 0.5) {
            r = y - 2.0 * std::floor(0.5 * y);
            if (r == 1.0) {
            rndup:
                y += 1.0;
            }
        }

        /* Else round down. */
        return (y);
    }

} // namespace cephes
} // namespace special
