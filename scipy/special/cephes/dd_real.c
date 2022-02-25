/*
 * src/double2.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract numbers DE-AC03-76SF00098 and
 * DE-AC02-05CH11231.
 * 
 * Copyright (c) 2003-2009, The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory (subject to receipt of
 * any required approvals from U.S. Dept. of Energy) All rights reserved.
 *
 * By downloading or using this software you are agreeing to the modified
 * BSD license "BSD-LBNL-License.doc" (see LICENSE.txt).
 */
/*
 * Contains implementation of non-inlined functions of double-double
 * package.  Inlined functions are found in dd_real_inline.h.
 */

/*
 * This code taken from v2.3.18 of the qd package.
*/

#include "mconf.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

#include "dd_real.h"

#define _DD_REAL_INIT(A, B)  {{A, B}}

const double DD_C_EPS = 4.93038065763132e-32; // 2^-104
const double DD_C_MIN_NORMALIZED = 2.0041683600089728e-292; // = 2^(-1022 + 53)

/*  Compile-time initialization of const double2 structs */

const double2 DD_C_MAX =
    _DD_REAL_INIT(1.79769313486231570815e+308, 9.97920154767359795037e+291);
const double2 DD_C_SAFE_MAX =
    _DD_REAL_INIT(1.7976931080746007281e+308, 9.97920154767359795037e+291);
const int _DD_C_NDIGITS = 31;

const double2 DD_C_ZERO = _DD_REAL_INIT(0.0, 0.0);
const double2 DD_C_ONE = _DD_REAL_INIT(1.0, 0.0);

const double2 DD_C_2PI =
    _DD_REAL_INIT(6.283185307179586232e+00, 2.449293598294706414e-16);
const double2 DD_C_PI =
    _DD_REAL_INIT(3.141592653589793116e+00, 1.224646799147353207e-16);
const double2 DD_C_PI2 =
    _DD_REAL_INIT(1.570796326794896558e+00, 6.123233995736766036e-17);
const double2 DD_C_PI4 =
    _DD_REAL_INIT(7.853981633974482790e-01, 3.061616997868383018e-17);
const double2 DD_C_PI16 =
    _DD_REAL_INIT(1.963495408493620697e-01, 7.654042494670957545e-18);
const double2 DD_C_3PI4 =
    _DD_REAL_INIT(2.356194490192344837e+00, 9.1848509936051484375e-17);

const double2 DD_C_E =
    _DD_REAL_INIT(2.718281828459045091e+00, 1.445646891729250158e-16);
const double2 DD_C_LOG2 =
    _DD_REAL_INIT(6.931471805599452862e-01, 2.319046813846299558e-17);
const double2 DD_C_LOG10 =
    _DD_REAL_INIT(2.302585092994045901e+00, -2.170756223382249351e-16);

#ifdef DD_C_NAN_IS_CONST
const double2 DD_C_NAN = _DD_REAL_INIT(NAN, NAN);
const double2 DD_C_INF = _DD_REAL_INIT(INFINITY, INFINITY);
const double2 DD_C_NEGINF = _DD_REAL_INIT(-INFINITY, -INFINITY);
#endif /* NAN */


/* This routine is called whenever a fatal error occurs. */
static volatile int errCount = 0;
void
dd_error(const char *msg)
{
    errCount++;
    /* if (msg) { */
    /*     fprintf(stderr, "ERROR %s\n", msg); */
    /* } */
}


int
get_double_expn(double x)
{
    int i = 0;
    double y;
    if (x == 0.0) {
        return INT_MIN;
    }
    if (DD_ISINF(x) || DD_ISNAN(x)) {
        return INT_MAX;
    }

    y = fabs(x);
    if (y < 1.0) {
        while (y < 1.0) {
            y *= 2.0;
            i++;
        }
        return -i;
    } else if (y >= 2.0) {
        while (y >= 2.0) {
            y *= 0.5;
            i++;
        }
        return i;
    }
    return 0;
}

/* ######################################################################## */
/* # Exponentiation */
/* ######################################################################## */

/* Computes the square root of the double-double number dd.
   NOTE: dd must be a non-negative number.                   */

double2
dd_sqrt(const double2 a)
{
    /* Strategy:  Use Karp's trick:  if x is an approximation
       to sqrt(a), then

          sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)

       The approximation is accurate to twice the accuracy of x.
       Also, the multiplication (a*x) and [-]*x can be done with
       only half the precision.
    */
    double x, ax;

    if (dd_is_zero(a))
        return DD_C_ZERO;

    if (dd_is_negative(a)) {
        dd_error("(dd_sqrt): Negative argument.");
        return DD_C_NAN;
    }

    x = 1.0 / sqrt(a.x[0]);
    ax = a.x[0] * x;
    return dd_add_d_d(ax, dd_sub(a, dd_sqr_d(ax)).x[0] * (x * 0.5));
}

/* Computes the square root of a double in double-double precision.
   NOTE: d must not be negative.                                   */

double2
dd_sqrt_d(double d)
{
    return dd_sqrt(dd_create_d(d));
}

/* Computes the n-th root of the double-double number a.
   NOTE: n must be a positive integer.
   NOTE: If n is even, then a must not be negative.       */

double2
dd_nroot(const double2 a, int n)
{
    /* Strategy:  Use Newton iteration for the function

            f(x) = x^(-n) - a

       to find its root a^{-1/n}.  The iteration is thus

            x' = x + x * (1 - a * x^n) / n

       which converges quadratically.  We can then find
      a^{1/n} by taking the reciprocal.
    */
    double2 r, x;

    if (n <= 0) {
        dd_error("(dd_nroot): N must be positive.");
        return DD_C_NAN;
    }

    if (n % 2 == 0 && dd_is_negative(a)) {
        dd_error("(dd_nroot): Negative argument.");
        return DD_C_NAN;
    }

    if (n == 1) {
        return a;
    }
    if (n == 2) {
        return dd_sqrt(a);
    }

    if (dd_is_zero(a))
        return DD_C_ZERO;

    /* Note  a^{-1/n} = exp(-log(a)/n) */
    r = dd_abs(a);
    x = dd_create_d(exp(-log(r.x[0]) / n));

    /* Perform Newton's iteration. */
    x = dd_add(
        x, dd_mul(x, dd_sub_d_dd(1.0, dd_div_dd_d(dd_mul(r, dd_npwr(x, n)),
                                                  DD_STATIC_CAST(double, n)))));
    if (a.x[0] < 0.0) {
        x = dd_neg(x);
    }
    return dd_inv(x);
}

/* Computes the n-th power of a double-double number.
   NOTE:  0^0 causes an error.                         */

double2
dd_npwr(const double2 a, int n)
{
    double2 r = a;
    double2 s = DD_C_ONE;
    int N = abs(n);
    if (N == 0) {
        if (dd_is_zero(a)) {
            dd_error("(dd_npwr): Invalid argument.");
            return DD_C_NAN;
        }
        return DD_C_ONE;
    }

    if (N > 1) {
        /* Use binary exponentiation */
        while (N > 0) {
            if (N % 2 == 1) {
                s = dd_mul(s, r);
            }
            N /= 2;
            if (N > 0) {
                r = dd_sqr(r);
            }
        }
    }
    else {
        s = r;
    }

    /* Compute the reciprocal if n is negative. */
    if (n < 0) {
        return dd_inv(s);
    }

    return s;
}

double2
dd_npow(const double2 a, int n)
{
    return dd_npwr(a, n);
}

double2
dd_pow(const double2 a, const double2 b)
{
    return dd_exp(dd_mul(b, dd_log(a)));
}

/* ######################################################################## */
/* # Exp/Log functions */
/* ######################################################################## */

static const double2 inv_fact[] = {
    {{1.66666666666666657e-01,  9.25185853854297066e-18}},
    {{4.16666666666666644e-02,  2.31296463463574266e-18}},
    {{8.33333333333333322e-03,  1.15648231731787138e-19}},
    {{1.38888888888888894e-03, -5.30054395437357706e-20}},
    {{1.98412698412698413e-04,  1.72095582934207053e-22}},
    {{2.48015873015873016e-05,  2.15119478667758816e-23}},
    {{2.75573192239858925e-06, -1.85839327404647208e-22}},
    {{2.75573192239858883e-07,  2.37677146222502973e-23}},
    {{2.50521083854417202e-08, -1.44881407093591197e-24}},
    {{2.08767569878681002e-09, -1.20734505911325997e-25}},
    {{1.60590438368216133e-10,  1.25852945887520981e-26}},
    {{1.14707455977297245e-11,  2.06555127528307454e-28}},
    {{7.64716373181981641e-13,  7.03872877733453001e-30}},
    {{4.77947733238738525e-14,  4.39920548583408126e-31}},
    {{2.81145725434552060e-15,  1.65088427308614326e-31}}
};
//static const int n_inv_fact = sizeof(inv_fact) / sizeof(inv_fact[0]);

/* Exponential.  Computes exp(x) in double-double precision. */

double2
dd_exp(const double2 a)
{
    /* Strategy:  We first reduce the size of x by noting that

            exp(kr + m * log(2)) = 2^m * exp(r)^k

       where m and k are integers.  By choosing m appropriately
       we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
       evaluated using the familiar Taylor series.  Reducing the
       argument substantially speeds up the convergence.       */

    const double k = 512.0;
    const double inv_k = 1.0 / k;
    double m;
    double2 r, s, t, p;
    int i = 0;

    if (a.x[0] <= -709.0) {
        return DD_C_ZERO;
    }

    if (a.x[0] >= 709.0) {
        return DD_C_INF;
    }

    if (dd_is_zero(a)) {
        return DD_C_INF;
    }

    if (dd_is_one(a)) {
        return DD_C_E;
    }

    m = floor(a.x[0] / DD_C_LOG2.x[0] + 0.5);
    r = dd_mul_pwr2(dd_sub(a, dd_mul_dd_d(DD_C_LOG2, m)), inv_k);

    p = dd_sqr(r);
    s = dd_add(r, dd_mul_pwr2(p, 0.5));
    p = dd_mul(p, r);
    t = dd_mul(p, inv_fact[0]);
    do {
        s = dd_add(s, t);
        p = dd_mul(p, r);
        ++i;
        t = dd_mul(p, inv_fact[i]);
    } while (fabs(dd_to_double(t)) > inv_k * DD_C_EPS && i < 5);

    s = dd_add(s, t);

    s = dd_add(dd_mul_pwr2(s, 2.0), dd_sqr(s));
    s = dd_add(dd_mul_pwr2(s, 2.0), dd_sqr(s));
    s = dd_add(dd_mul_pwr2(s, 2.0), dd_sqr(s));
    s = dd_add(dd_mul_pwr2(s, 2.0), dd_sqr(s));
    s = dd_add(dd_mul_pwr2(s, 2.0), dd_sqr(s));
    s = dd_add(dd_mul_pwr2(s, 2.0), dd_sqr(s));
    s = dd_add(dd_mul_pwr2(s, 2.0), dd_sqr(s));
    s = dd_add(dd_mul_pwr2(s, 2.0), dd_sqr(s));
    s = dd_add(dd_mul_pwr2(s, 2.0), dd_sqr(s));
    s = dd_add(s, DD_C_ONE);

    return dd_ldexp(s, DD_STATIC_CAST(int, m));
}

double2
dd_exp_d(const double a)
{
    return dd_exp(dd_create(a, 0));
}


/* Logarithm.  Computes log(x) in double-double precision.
   This is a natural logarithm (i.e., base e).            */
double2
dd_log(const double2 a)
{
    /* Strategy.  The Taylor series for log converges much more
       slowly than that of exp, due to the lack of the factorial
       term in the denominator.  Hence this routine instead tries
       to determine the root of the function

           f(x) = exp(x) - a

       using Newton iteration.  The iteration is given by

           x' = x - f(x)/f'(x)
              = x - (1 - a * exp(-x))
              = x + a * exp(-x) - 1.

       Only one iteration is needed, since Newton's iteration
       approximately doubles the number of digits per iteration. */
    double2 x;

    if (dd_is_one(a)) {
        return DD_C_ZERO;
    }

    if (a.x[0] <= 0.0) {
        dd_error("(dd_log): Non-positive argument.");
        return DD_C_NAN;
    }

    x = dd_create_d(log(a.x[0])); /* Initial approximation */

    /* x = x + a * exp(-x) - 1.0; */
    x = dd_add(x, dd_sub(dd_mul(a, dd_exp(dd_neg(x))), DD_C_ONE));
    return x;
}


double2
dd_log1p(const double2 a)
{
    double2 ans;
    double la, elam1, ll;
    if (a.x[0] <= -1.0) {
        return DD_C_NEGINF;
    }
    la = log1p(a.x[0]);
    elam1 = expm1(la);
    ll = log1p(a.x[1] / (1 + a.x[0]));
    if (a.x[0] > 0) {
        ll -= (elam1 - a.x[0])/(elam1+1);
    }
    ans = dd_add_d_d(la, ll);
    return ans;
}

double2
dd_expm1(const double2 a)
{
   double eam1 = expm1(a.x[0]);
   return dd_add_d_d(eam1, expm1(a.x[1])*(eam1+1));
}

double2
dd_log10(const double2 a)
{
    return dd_div(dd_log(a), DD_C_LOG10);
}

double2
dd_log_d(double a)
{
    return dd_log(dd_create(a, 0));
}


double2
dd_rand(void)
{
    static const double m_const = 4.6566128730773926e-10; /* = 2^{-31} */
    double m = m_const;
    double2 r = DD_C_ZERO;
    double d;
    int i;

    /* Strategy:  Generate 31 bits at a time, using lrand48
       random number generator.  Shift the bits, and reapeat
       4 times. */

    for (i = 0; i < 4; i++, m *= m_const) {
        //    d = lrand48() * m;
        d = rand() * m;
        r = dd_add_dd_d(r, d);
    }

    return r;
}

/* polyeval(c, n, x)
   Evaluates the given n-th degree polynomial at x.
   The polynomial is given by the array of (n+1) coefficients. */

double2
polyeval(const double2 *c, int n, const double2 x)
{
    /* Just use Horner's method of polynomial evaluation. */
    double2 r = c[n];
    int i;

    for (i = n - 1; i >= 0; i--) {
        r = dd_mul(r, x);
        r = dd_add(r, c[i]);
    }

    return r;
}

/* polyroot(c, n, x0)
   Given an n-th degree polynomial, finds a root close to
   the given guess x0.  Note that this uses simple Newton
   iteration scheme, and does not work for multiple roots.  */

double2
polyroot(const double2 *c, int n, const double2 x0, int max_iter,
         double thresh)
{
    double2 x = x0;
    double2 f;
    double2 *d = DD_STATIC_CAST(double2 *, calloc(sizeof(double2), n));
    int conv = 0;
    int i;
    double max_c = fabs(dd_to_double(c[0]));
    double v;

    if (thresh == 0.0) {
        thresh = DD_C_EPS;
    }

    /* Compute the coefficients of the derivatives. */
    for (i = 1; i <= n; i++) {
        v = fabs(dd_to_double(c[i]));
        if (v > max_c) {
            max_c = v;
        }
        d[i - 1] = dd_mul_dd_d(c[i], DD_STATIC_CAST(double, i));
    }
    thresh *= max_c;

    /* Newton iteration. */
    for (i = 0; i < max_iter; i++) {
        f = polyeval(c, n, x);

        if (fabs(dd_to_double(f)) < thresh) {
            conv = 1;
            break;
        }
        x = dd_sub(x, (dd_div(f, polyeval(d, n - 1, x))));
    }
    free(d);

    if (!conv) {
        dd_error("(dd_polyroot): Failed to converge.");
        return DD_C_NAN;
    }

    return x;
}
