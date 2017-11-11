/* MIT License
 *
 * Copyright (c) 2016--2017 Felix Lenders
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include "trlib.h"
#include "trlib_private.h"

#include "_c99compat.h"

trlib_int_t trlib_quadratic_zero(trlib_flt_t c_abs, trlib_flt_t c_lin, trlib_flt_t tol,
        trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
        trlib_flt_t *t1, trlib_flt_t *t2) {
    trlib_int_t n  = 0;   // number of roots
    trlib_flt_t q = 0.0;
    trlib_flt_t dq = 0.0;
    trlib_flt_t lin_sq = c_lin*c_lin;
    *t1 = 0.0;    // first root
    *t2 = 0.0;    // second root

    if (fabs(c_abs) > tol*lin_sq) {
        // well behaved non-degenerate quadratic
        // compute discriminant
        q = lin_sq - 4.0 * c_abs;
        if ( fabs(q) <= (TRLIB_EPS*c_lin)*(TRLIB_EPS*c_lin) ) {
            // two distinct zeros, but discrimant tiny --> numeric double zero
            // initialize on same root obtained by standard formula with zero discrement, let newton refinement do the rest
            n = 2;
            *t1 = -.5*c_lin; *t2 = *t1;
        }
        else if ( q < 0.0 ) {
            n = 2;
            *t1 = 0.0; *t2 = 0.0;
            return n;
        }
        else {
            // discriminant large enough, two distinc zeros
            n = 2;
            // start with root according to plus sign to avoid cancellation
            *t1 = -.5 * ( c_lin + copysign( sqrt(q), c_lin ) );
            *t2 = c_abs/(*t1);
            if (*t2 < *t1) { q = *t2; *t2 = *t1; *t1 = q; }
        }
    }
    else {
        n = 2;
        if (c_lin < 0.0) { *t1 = 0.0; *t2 = - c_lin; }
        else { *t1 = - c_lin; *t2 = 0.0; }
    }

    // newton correction
    q = (*t1+c_lin)*(*t1)+c_abs; dq = 2.0*(*t1)+c_lin;
    if (dq != 0.0) { *t1 = *t1 - q/dq; }
    q = (*t2+c_lin)*(*t2)+c_abs; dq = 2.0*(*t2)+c_lin;
    if (dq != 0.0) { *t2 = *t2 - q/dq; }
    return n;
}
