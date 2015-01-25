/*                                                     ellie.c
 *
 *     Incomplete elliptic integral of the second kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double phi, m, y, ellie();
 *
 * y = ellie( phi, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *                phi
 *                 -
 *                | |
 *                |                   2
 * E(phi_\m)  =    |    sqrt( 1 - m sin t ) dt
 *                |
 *              | |    
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
 *
 *
 * ACCURACY:
 *
 * Tested at random arguments with phi in [-10, 10] and m in
 * [0, 1].
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC        0,2         2000       1.9e-16     3.4e-17
 *    IEEE     -10,10      150000       3.3e-15     1.4e-16
 *
 *
 */


/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987, 1993 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */
/* Copyright 2014, Eric W. Moore */

/*     Incomplete elliptic integral of second kind     */

#include "mconf.h"

extern double MACHEP;

static double ellie_neg_m(double phi, double m);

double ellie(double phi, double m)
{
    double a, b, c, e, temp;
    double lphi, t, E, denom, npio2;
    int d, mod, sign;

    if (isnan(phi) || isnan(m))
        return NPY_NAN;
    if (m > 1.0)
        return NPY_NAN;
    if (isinf(phi))
        return phi;
    if (isinf(m))
        return -m;
    if (m == 0.0)
	return (phi);
    lphi = phi;
    npio2 = floor(lphi / NPY_PI_2);
    if (fmod(fabs(npio2), 2.0) == 1.0)
	npio2 += 1;
    lphi = lphi - npio2 * NPY_PI_2;
    if (lphi < 0.0) {
	lphi = -lphi;
	sign = -1;
    }
    else {
	sign = 1;
    }
    a = 1.0 - m;
    E = ellpe(m);
    if (a == 0.0) {
	temp = sin(lphi);
	goto done;
    }
    if (a > 1.0) {
        temp = ellie_neg_m(lphi, m);
        goto done;
    }
 
    if (lphi < 0.135) {
        double m11= (((((-7.0/2816.0)*m + (5.0/1056.0))*m - (7.0/2640.0))*m
                    + (17.0/41580.0))*m - (1.0/155925.0))*m;
        double m9 = ((((-5.0/1152.0)*m + (1.0/144.0))*m - (1.0/360.0))*m
                    + (1.0/5670.0))*m;
        double m7 = ((-m/112.0 + (1.0/84.0))*m - (1.0/315.0))*m;
        double m5 = (-m/40.0 + (1.0/30))*m;
        double m3 = -m/6.0;
        double p2 = lphi * lphi;

        temp = ((((m11*p2 + m9)*p2 + m7)*p2 + m5)*p2 + m3)*p2*lphi + lphi;
        goto done;
    }
    t = tan(lphi);
    b = sqrt(a);
    /* Thanks to Brian Fitzgerald <fitzgb@mml0.meche.rpi.edu>
     * for pointing out an instability near odd multiples of pi/2.  */
    if (fabs(t) > 10.0) {
	/* Transform the amplitude */
	e = 1.0 / (b * t);
	/* ... but avoid multiple recursions.  */
	if (fabs(e) < 10.0) {
	    e = atan(e);
	    temp = E + m * sin(lphi) * sin(e) - ellie(e, m);
	    goto done;
	}
    }
    c = sqrt(m);
    a = 1.0;
    d = 1;
    e = 0.0;
    mod = 0;

    while (fabs(c / a) > MACHEP) {
	temp = b / a;
	lphi = lphi + atan(t * temp) + mod * NPY_PI;
        denom = 1 - temp * t * t;
        if (fabs(denom) > 10*MACHEP) {
            t = t * (1.0 + temp) / denom;
            mod = (lphi + NPY_PI_2) / NPY_PI;
        }
        else {
            t = tan(lphi);
            mod = (int)floor((lphi - atan(t))/NPY_PI);
        }
	c = (a - b) / 2.0;
	temp = sqrt(a * b);
	a = (a + b) / 2.0;
	b = temp;
	d += d;
	e += c * sin(lphi);
    }

    temp = E / ellpk(1.0 - m);
    temp *= (atan(t) + mod * NPY_PI) / (d * a);
    temp += e;

  done:

    if (sign < 0)
	temp = -temp;
    temp += npio2 * E;
    return (temp);
}

/* N.B. This will evaluate its arguments multiple times. */
#define MAX3(a, b, c) (a > b ? (a > c ? a : c) : (b > c ? b : c))

/* To calculate legendre's incomplete elliptical integral of the second kind for
 * negative m, we use a power series in phi for small m*phi*phi, an asymptotic
 * series in m for large m*phi*phi* and the relation to Carlson's symmetric
 * integrals, R_F(x,y,z) and R_D(x,y,z).
 * 
 * E(phi, m) = sin(phi) * R_F(cos(phi)^2, 1 - m * sin(phi)^2, 1.0)
 *             - m * sin(phi)^3 * R_D(cos(phi)^2, 1 - m * sin(phi)^2, 1.0) / 3
 *             
 *           = R_F(c-1, c-m, c) - m * R_D(c-1, c-m, c) / 3
 *
 * where c = csc(phi)^2. We use the second form of this for (approximately)
 * phi > 1/(sqrt(DBL_MAX) ~ 1e-154, where csc(phi)^2 overflows. Elsewhere we
 * use the first form, accounting for the smallness of phi.
 * 
 * The algorithm used is described in Carlson, B. C. Numerical computation of
 * real or complex elliptic integrals. (1994) http://arxiv.org/abs/math/9409227
 * Most variable names reflect Carlson's usage.
 *
 * In this routine, we assume m < 0 and  0 > phi > pi/2.
 */
double ellie_neg_m(double phi, double m)
{
    double x, y, z, x1, y1, z1, ret, Q;
    double A0f, Af, Xf, Yf, Zf, E2f, E3f, scalef;
    double A0d, Ad, seriesn, seriesd, Xd, Yd, Zd, E2d, E3d, E4d, E5d, scaled;
    int n = 0;
    double mpp = (m*phi)*phi;
   
    if (-mpp < 1e-6 && phi < -m) {
        return phi + (mpp*phi*phi/30.0 - mpp*mpp/40.0 - mpp/6.0)*phi;
    }

    if (-mpp > 1e6) {
        double sm = sqrt(-m);
        double sp = sin(phi);
        double cp = cos(phi);

        double a = -cosm1(phi);
        double b1 = log(4*sp*sm/(1+cp));
        double b = -(0.5 + b1) / 2.0 / m;
        double c = (0.75 + cp/sp/sp - b1) / 16.0 / m / m;
        return (a + b + c) * sm;
    }

    if (phi > 1e-153 && m > -1e200) {
        double s = sin(phi);
        double csc2 = 1.0 / s / s;
        scalef = 1.0;
        scaled = m / 3.0;
        x = 1.0 / tan(phi) / tan(phi);
        y = csc2 - m;
        z = csc2;
    }
    else {
        scalef = phi;
        scaled = mpp * phi / 3.0;
        x = 1.0;
        y = 1 - mpp;
        z = 1.0;
    }
    
    if (x == y && x == z) {
        return (scalef + scaled/x)/sqrt(x);
    }

    A0f = (x + y + z) / 3.0;
    Af = A0f;
    A0d = (x + y + 3.0*z) / 5.0;
    Ad = A0d;
    x1 = x; y1 = y; z1 = z; seriesd = 0.0; seriesn = 1.0;
    /* Carlson gives 1/pow(3*r, 1.0/6.0) for this constant. if r == eps,
     * it is ~338.38. */
    Q = 400.0 * MAX3(fabs(A0f-x), fabs(A0f-y), fabs(A0f-z));
    
    while (Q > fabs(Af) && Q > fabs(Ad) && n <= 100) {
        double sx = sqrt(x1);
        double sy = sqrt(y1);
        double sz = sqrt(z1);
        double lam = sx*sy + sx*sz + sy*sz;
        seriesd += seriesn / (sz * (z1 + lam));
        x1 = (x1 + lam) / 4.0;
        y1 = (y1 + lam) / 4.0;
        z1 = (z1 + lam) / 4.0;
        Af = (x1 + y1 + z1) / 3.0;
        Ad = (Ad + lam) / 4.0;
        n += 1;
        Q /= 4.0;
        seriesn /= 4.0;
    }

    Xf = (A0f - x) / Af / (1 << 2*n);
    Yf = (A0f - y) / Af / (1 << 2*n);
    Zf = -(Xf + Yf);

    E2f = Xf*Yf - Zf*Zf;
    E3f = Xf*Yf*Zf;

    ret = scalef * (1.0 - E2f/10.0 + E3f/14.0 + E2f*E2f/24.0
                    - 3.0*E2f*E3f/44.0) / sqrt(Af);

    Xd = (A0d - x) / Ad / (1 << 2*n);
    Yd = (A0d - y) / Ad / (1 << 2*n);
    Zd = -(Xd + Yd)/3.0;

    E2d = Xd*Yd - 6.0*Zd*Zd;
    E3d = (3*Xd*Yd - 8.0*Zd*Zd)*Zd;
    E4d = 3.0*(Xd*Yd - Zd*Zd)*Zd*Zd;
    E5d = Xd*Yd*Zd*Zd*Zd;

    ret -= scaled * (1.0 - 3.0*E2d/14.0 + E3d/6.0 + 9.0*E2d*E2d/88.0
                     - 3.0*E4d/22.0 - 9.0*E2d*E3d/52.0 + 3.0*E5d/26.0)
                     /(1 << 2*n) / Ad / sqrt(Ad);
    ret -= 3.0 * scaled * seriesd;
    return ret;
}

