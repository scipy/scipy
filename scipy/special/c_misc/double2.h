/*
 * 2-double floating point numbers
 *
 * Almost 2 times the precision of double, by using 2 doubles.
 */

/*
   Portions from:

   Quad Double computation package
   Copyright (C) 2003-2012
   ================================================

   Revision date:  13 Mar 2012

   Authors:
   Yozo Hida               U.C. Berkeley                   yozo@cs.berkeley.edu
   Xiaoye S. Li            Lawrence Berkeley Natl Lab      xiaoye@nersc.gov
   David H. Bailey         Lawrence Berkeley Natl Lab      dhbailey@lbl.gov
  
   C++ usage guide:
   Alex Kaiser             Lawrence Berkeley Natl Lab      adkaiser@lbl.gov
   

   Berkeley Software Distribution Agreement
   
   This License Agreement is entered into by The Regents of the University
   of California, Department of Energy contract-operators of the Lawrence
   Berkeley National Laboratory, 1 Cyclotron Road, Berkeley, CA 94720
   (“Berkeley Lab”), and the entity listed below (“you” or
   "Licensee") having its place of business at the address below:
   
   The parties now agree as follows:
   
   1. Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
   
   (1) Redistributions of source code must retain the copyright notice,
   this list of conditions and the following disclaimer.
   
   (2) Redistributions in binary form must reproduce the copyright notice,
   this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   
   (3) Neither the name of the University of California, Lawrence Berkeley
   National Laboratory, U.S. Dept. of Energy nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.
   
   2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
   PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   
   3. You are under no obligation whatsoever to provide any bug fixes,
   patches, or upgrades to the features, functionality or performance of
   the source code ("Enhancements") to anyone; however, if you choose to
   make your Enhancements available either publicly, or directly to
   Lawrence Berkeley National Laboratory, without imposing a separate
   written license agreement for such Enhancements, then you hereby grant
   the following license: a non-exclusive, royalty-free perpetual license
   to install, use, modify, prepare derivative works, incorporate into
   other computer software, distribute, and sublicense such enhancements or
   derivative works thereof, in binary and source code form.
*/

#ifndef DOUBLEN_H_
#define DOUBLEN_H_

#include <math.h>


typedef struct
{
    double x[2];
} double2_t;


static void double2_init(double2_t *a, double y)
{
    a->x[0] = y;
    a->x[1] = 0;
}


static void double2_init2(double2_t *a, double hi, double lo)
{
    a->x[0] = hi;
    a->x[1] = lo;
}


static double double_sum_err(double a, double b, double *err)
{
    double tmp;
    volatile double c, d, e, g, h, f, x;

    if (fabs(a) < fabs(b)) {
        tmp = b;
        b = a;
        a = tmp;
    }

    c = a + b;
    e = c - a;
    g = c - e;
    h = g - a;
    f = b - h;
    d = f - e;
    x = d + e;
    if (x != f) {
        c = a;
        d = b;
    }
    *err = d;
    return c;
}


#define _QD_SPLIT_THRESH 6.69692879491417e+299   /* = 2^996 */
#define _QD_SPLITTER 134217729.0                 /* = 2^27 + 1 */

static void double_split(double a, double *hi, double *lo)
{
    volatile double b, c;
    if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH) {
        a *= 3.7252902984619140625e-09;  /* 2^-28 */
        b = _QD_SPLITTER * a;
        c = b - a;
        *hi = b - c;
        *lo = a - *hi;
        *hi *= 268435456.0;          /* 2^28 */
        *lo *= 268435456.0;          /* 2^28 */
    } else {
        b = _QD_SPLITTER * a;
        c = b - a;
        *hi = b - c;
        *lo = a - *hi;
    }
}


static double double_mul_err(double a, double b, double *err)
{
    double a_hi, a_lo, b_hi, b_lo;
    double p = a * b;
    volatile double c, d;
    double_split(a, &a_hi, &a_lo);
    double_split(b, &b_hi, &b_lo);
    c = a_hi * b_hi - p;
    d = c + a_hi * b_lo + a_lo * b_hi;
    *err = d + a_lo * b_lo;
    return p;
}


static void double2_add(double2_t *a, double2_t *b, double2_t *c)
{
    double s1, s2, t1, t2;

    s1 = double_sum_err(a->x[0], b->x[0], &s2);
    t1 = double_sum_err(a->x[1], b->x[1], &t2);

    s2 += t1;
    s1 = double_sum_err(s1, s2, &s2);
    s2 += t2;
    s1 = double_sum_err(s1, s2, &s2);

    double2_init2(c, s1, s2);
}


static void double2_neg(double2_t *a, double2_t *b)
{
    b->x[0] = -a->x[0];
    b->x[1] = -a->x[1];
}


static void double2_sub(double2_t *a, double2_t *b, double2_t *c)
{
    double2_t d;
    double2_neg(b, &d);
    double2_add(a, &d, c);
}


static void double2_mul(double2_t *a, double2_t *b, double2_t *c)
{
    double p1, p2;
    
    p1 = double_mul_err(a->x[0], b->x[0], &p2);
    p2 += (a->x[0] * b->x[1] + a->x[1] * b->x[0]);
    p1 = double_sum_err(p1, p2, &p2);

    double2_init2(c, p1, p2);
}


static void double2_div(double2_t *a, double2_t *b, double2_t *c)
{
    double q1, q2, q3;
    double2_t r, e, f;

    q1 = a->x[0] / b->x[0];  /* approximate quotient */

    double2_init(&f, q1);
    double2_mul(&f, b, &e);
    double2_sub(a, &e, &r);

    q2 = r.x[0] / b->x[0];

    double2_init(&f, q2);
    double2_mul(&f, b, &e);
    double2_sub(&r, &e, &r);

    q3 = r.x[0] / b->x[0];

    q1 = double_sum_err(q1, q2, &q2);

    double2_init2(&e, q1, q2);
    double2_init(&f, q3);
    double2_add(&e, &f, c);
}


static double double2_double(double2_t *a)
{
    return a->x[0] + a->x[1];
}

#endif /* DOUBLE2_H_ */
