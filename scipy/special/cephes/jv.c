/*                                                     jv.c
 *
 *     Bessel function of noninteger order
 *
 *
 *
 * SYNOPSIS:
 *
 * double v, x, y, jv();
 *
 * y = jv( v, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order v of the argument,
 * where v is real.  Negative x is allowed if v is an integer.
 *
 * Several expansions are included: the ascending power
 * series, the Hankel expansion, and two transitional
 * expansions for large v.  If v is not too large, it
 * is reduced by recurrence to a region of best accuracy.
 * The transitional expansions give 12D accuracy for v > 500.
 *
 *
 *
 * ACCURACY:
 * Results for integer v are indicated by *, where x and v
 * both vary from -125 to +125.  Otherwise,
 * x ranges from 0 to 125, v ranges as indicated by "domain."
 * Error criterion is absolute, except relative when |jv()| > 1.
 *
 * arithmetic  v domain  x domain    # trials      peak       rms
 *    IEEE      0,125     0,125      100000      4.6e-15    2.2e-16
 *    IEEE   -125,0       0,125       40000      5.4e-11    3.7e-13
 *    IEEE      0,500     0,500       20000      4.4e-15    4.0e-16
 * Integer v:
 *    IEEE   -125,125   -125,125      50000      3.5e-15*   1.9e-16*
 *
 */


/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
 */


#include "mconf.h"
#define CEPHES_DEBUG 0

#if CEPHES_DEBUG
#include <stdio.h>
#endif

#define MAXGAM 171.624376956302725

extern double MACHEP, MINLOG, MAXLOG;

#define BIG  1.44115188075855872E+17

static double jvs(double n, double x);
static double hankel(double n, double x);
static double recur(double *n, double x, double *newn, int cancel);
static double jnx(double n, double x);
static double jnt(double n, double x);

double jv(double n, double x)
{
    double k, q, t, y, an;
    int i, sign, nint;

    nint = 0;			/* Flag for integer n */
    sign = 1;			/* Flag for sign inversion */
    an = fabs(n);
    y = floor(an);
    if (y == an) {
	nint = 1;
	i = an - 16384.0 * floor(an / 16384.0);
	if (n < 0.0) {
	    if (i & 1)
		sign = -sign;
	    n = an;
	}
	if (x < 0.0) {
	    if (i & 1)
		sign = -sign;
	    x = -x;
	}
	if (n == 0.0)
	    return (j0(x));
	if (n == 1.0)
	    return (sign * j1(x));
    }

    if ((x < 0.0) && (y != an)) {
	mtherr("Jv", DOMAIN);
	y = NPY_NAN;
	goto done;
    }

    if (x == 0 && n < 0 && !nint) {
        mtherr("Jv", OVERFLOW);
        return NPY_INFINITY / gamma(n + 1);
    }

    y = fabs(x);

    if (y * y < fabs(n + 1) * MACHEP) {
	return pow(0.5 * x, n) / gamma(n + 1);
    }

    k = 3.6 * sqrt(y);
    t = 3.6 * sqrt(an);
    if ((y < t) && (an > 21.0))
	return (sign * jvs(n, x));
    if ((an < k) && (y > 21.0))
	return (sign * hankel(n, x));

    if (an < 500.0) {
	/* Note: if x is too large, the continued fraction will fail; but then the
	 * Hankel expansion can be used. */
	if (nint != 0) {
	    k = 0.0;
	    q = recur(&n, x, &k, 1);
	    if (k == 0.0) {
		y = j0(x) / q;
		goto done;
	    }
	    if (k == 1.0) {
		y = j1(x) / q;
		goto done;
	    }
	}

	if (an > 2.0 * y)
	    goto rlarger;

	if ((n >= 0.0) && (n < 20.0)
	    && (y > 6.0) && (y < 20.0)) {
	    /* Recur backwards from a larger value of n */
	  rlarger:
	    k = n;

	    y = y + an + 1.0;
	    if (y < 30.0)
		y = 30.0;
	    y = n + floor(y - n);
	    q = recur(&y, x, &k, 0);
	    y = jvs(y, x) * q;
	    goto done;
	}

	if (k <= 30.0) {
	    k = 2.0;
	}
	else if (k < 90.0) {
	    k = (3 * k) / 4;
	}
	if (an > (k + 3.0)) {
	    if (n < 0.0)
		k = -k;
	    q = n - floor(n);
	    k = floor(k) + q;
	    if (n > 0.0)
		q = recur(&n, x, &k, 1);
	    else {
		t = k;
		k = n;
		q = recur(&t, x, &k, 1);
		k = t;
	    }
	    if (q == 0.0) {
		y = 0.0;
		goto done;
	    }
	}
	else {
	    k = n;
	    q = 1.0;
	}

	/* boundary between convergence of
	 * power series and Hankel expansion
	 */
	y = fabs(k);
	if (y < 26.0)
	    t = (0.0083 * y + 0.09) * y + 12.9;
	else
	    t = 0.9 * y;

	if (x > t)
	    y = hankel(k, x);
	else
	    y = jvs(k, x);
#if CEPHES_DEBUG
	printf("y = %.16e, recur q = %.16e\n", y, q);
#endif
	if (n > 0.0)
	    y /= q;
	else
	    y *= q;
    }

    else {
	/* For large n, use the uniform expansion or the transitional expansion.
	 * But if x is of the order of n**2, these may blow up, whereas the
	 * Hankel expansion will then work.
	 */
	if (n < 0.0) {
	    mtherr("Jv", TLOSS);
	    y = NPY_NAN;
	    goto done;
	}
	t = x / n;
	t /= n;
	if (t > 0.3)
	    y = hankel(n, x);
	else
	    y = jnx(n, x);
    }

  done:return (sign * y);
}

/* Reduce the order by backward recurrence.
 * AMS55 #9.1.27 and 9.1.73.
 */

static double recur(double *n, double x, double *newn, int cancel)
{
    double pkm2, pkm1, pk, qkm2, qkm1;

    /* double pkp1; */
    double k, ans, qk, xk, yk, r, t, kf;
    static double big = BIG;
    int nflag, ctr;
    int miniter, maxiter;

    /* Continued fraction for Jn(x)/Jn-1(x)
     * AMS 9.1.73
     *
     *    x       -x^2      -x^2
     * ------  ---------  ---------   ...
     * 2 n +   2(n+1) +   2(n+2) +
     *
     * Compute it with the simplest possible algorithm.
     *
     * This continued fraction starts to converge when (|n| + m) > |x|.
     * Hence, at least |x|-|n| iterations are necessary before convergence is
     * achieved. There is a hard limit set below, m <= 30000, which is chosen
     * so that no branch in `jv` requires more iterations to converge.
     * The exact maximum number is (500/3.6)^2 - 500 ~ 19000
     */

    maxiter = 22000;
    miniter = fabs(x) - fabs(*n);
    if (miniter < 1)
	miniter = 1;

    if (*n < 0.0)
	nflag = 1;
    else
	nflag = 0;

  fstart:

#if CEPHES_DEBUG
    printf("recur: n = %.6e, newn = %.6e, cfrac = ", *n, *newn);
#endif

    pkm2 = 0.0;
    qkm2 = 1.0;
    pkm1 = x;
    qkm1 = *n + *n;
    xk = -x * x;
    yk = qkm1;
    ans = 0.0;			/* ans=0.0 ensures that t=1.0 in the first iteration */
    ctr = 0;
    do {
	yk += 2.0;
	pk = pkm1 * yk + pkm2 * xk;
	qk = qkm1 * yk + qkm2 * xk;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;

	/* check convergence */
	if (qk != 0 && ctr > miniter)
	    r = pk / qk;
	else
	    r = 0.0;

	if (r != 0) {
	    t = fabs((ans - r) / r);
	    ans = r;
	}
	else {
	    t = 1.0;
	}

	if (++ctr > maxiter) {
	    mtherr("jv", UNDERFLOW);
	    goto done;
	}
	if (t < MACHEP)
	    goto done;

	/* renormalize coefficients */
	if (fabs(pk) > big) {
	    pkm2 /= big;
	    pkm1 /= big;
	    qkm2 /= big;
	    qkm1 /= big;
	}
    }
    while (t > MACHEP);

  done:
    if (ans == 0)
	ans = 1.0;

#if CEPHES_DEBUG
    printf("%.6e\n", ans);
#endif

    /* Change n to n-1 if n < 0 and the continued fraction is small */
    if (nflag > 0) {
	if (fabs(ans) < 0.125) {
	    nflag = -1;
	    *n = *n - 1.0;
	    goto fstart;
	}
    }


    kf = *newn;

    /* backward recurrence
     *              2k
     *  J   (x)  =  --- J (x)  -  J   (x)
     *   k-1         x   k         k+1
     */

    pk = 1.0;
    pkm1 = 1.0 / ans;
    k = *n - 1.0;
    r = 2 * k;
    do {
	pkm2 = (pkm1 * r - pk * x) / x;
	/*      pkp1 = pk; */
	pk = pkm1;
	pkm1 = pkm2;
	r -= 2.0;
	/*
	 * t = fabs(pkp1) + fabs(pk);
	 * if( (k > (kf + 2.5)) && (fabs(pkm1) < 0.25*t) )
	 * {
	 * k -= 1.0;
	 * t = x*x;
	 * pkm2 = ( (r*(r+2.0)-t)*pk - r*x*pkp1 )/t;
	 * pkp1 = pk;
	 * pk = pkm1;
	 * pkm1 = pkm2;
	 * r -= 2.0;
	 * }
	 */
	k -= 1.0;
    }
    while (k > (kf + 0.5));

    /* Take the larger of the last two iterates
     * on the theory that it may have less cancellation error.
     */

    if (cancel) {
	if ((kf >= 0.0) && (fabs(pk) > fabs(pkm1))) {
	    k += 1.0;
	    pkm2 = pk;
	}
    }
    *newn = k;
#if CEPHES_DEBUG
    printf("newn %.6e rans %.6e\n", k, pkm2);
#endif
    return (pkm2);
}



/* Ascending power series for Jv(x).
 * AMS55 #9.1.10.
 */

static double jvs(double n, double x)
{
    double t, u, y, z, k;
    int ex, sgngam;

    z = -x * x / 4.0;
    u = 1.0;
    y = u;
    k = 1.0;
    t = 1.0;

    while (t > MACHEP) {
	u *= z / (k * (n + k));
	y += u;
	k += 1.0;
	if (y != 0)
	    t = fabs(u / y);
    }
#if CEPHES_DEBUG
    printf("power series=%.5e ", y);
#endif
    t = frexp(0.5 * x, &ex);
    ex = ex * n;
    if ((ex > -1023)
	&& (ex < 1023)
	&& (n > 0.0)
	&& (n < (MAXGAM - 1.0))) {
	t = pow(0.5 * x, n) / gamma(n + 1.0);
#if CEPHES_DEBUG
	printf("pow(.5*x, %.4e)/gamma(n+1)=%.5e\n", n, t);
#endif
	y *= t;
    }
    else {
#if CEPHES_DEBUG
	z = n * log(0.5 * x);
	k = lgam(n + 1.0);
	t = z - k;
	printf("log pow=%.5e, lgam(%.4e)=%.5e\n", z, n + 1.0, k);
#else
	t = n * log(0.5 * x) - lgam_sgn(n + 1.0, &sgngam);
#endif
	if (y < 0) {
	    sgngam = -sgngam;
	    y = -y;
	}
	t += log(y);
#if CEPHES_DEBUG
	printf("log y=%.5e\n", log(y));
#endif
	if (t < -MAXLOG) {
	    return (0.0);
	}
	if (t > MAXLOG) {
	    mtherr("Jv", OVERFLOW);
	    return (NPY_INFINITY);
	}
	y = sgngam * exp(t);
    }
    return (y);
}

/* Hankel's asymptotic expansion
 * for large x.
 * AMS55 #9.2.5.
 */

static double hankel(double n, double x)
{
    double t, u, z, k, sign, conv;
    double p, q, j, m, pp, qq;
    int flag;

    m = 4.0 * n * n;
    j = 1.0;
    z = 8.0 * x;
    k = 1.0;
    p = 1.0;
    u = (m - 1.0) / z;
    q = u;
    sign = 1.0;
    conv = 1.0;
    flag = 0;
    t = 1.0;
    pp = 1.0e38;
    qq = 1.0e38;

    while (t > MACHEP) {
	k += 2.0;
	j += 1.0;
	sign = -sign;
	u *= (m - k * k) / (j * z);
	p += sign * u;
	k += 2.0;
	j += 1.0;
	u *= (m - k * k) / (j * z);
	q += sign * u;
	t = fabs(u / p);
	if (t < conv) {
	    conv = t;
	    qq = q;
	    pp = p;
	    flag = 1;
	}
	/* stop if the terms start getting larger */
	if ((flag != 0) && (t > conv)) {
#if CEPHES_DEBUG
	    printf("Hankel: convergence to %.4E\n", conv);
#endif
	    goto hank1;
	}
    }

  hank1:
    u = x - (0.5 * n + 0.25) * NPY_PI;
    t = sqrt(2.0 / (NPY_PI * x)) * (pp * cos(u) - qq * sin(u));
#if CEPHES_DEBUG
    printf("hank: %.6e\n", t);
#endif
    return (t);
}


/* Asymptotic expansion for large n.
 * AMS55 #9.3.35.
 */

static double lambda[] = {
    1.0,
    1.041666666666666666666667E-1,
    8.355034722222222222222222E-2,
    1.282265745563271604938272E-1,
    2.918490264641404642489712E-1,
    8.816272674437576524187671E-1,
    3.321408281862767544702647E+0,
    1.499576298686255465867237E+1,
    7.892301301158651813848139E+1,
    4.744515388682643231611949E+2,
    3.207490090890661934704328E+3
};

static double mu[] = {
    1.0,
    -1.458333333333333333333333E-1,
    -9.874131944444444444444444E-2,
    -1.433120539158950617283951E-1,
    -3.172272026784135480967078E-1,
    -9.424291479571202491373028E-1,
    -3.511203040826354261542798E+0,
    -1.572726362036804512982712E+1,
    -8.228143909718594444224656E+1,
    -4.923553705236705240352022E+2,
    -3.316218568547972508762102E+3
};

static double P1[] = {
    -2.083333333333333333333333E-1,
    1.250000000000000000000000E-1
};

static double P2[] = {
    3.342013888888888888888889E-1,
    -4.010416666666666666666667E-1,
    7.031250000000000000000000E-2
};

static double P3[] = {
    -1.025812596450617283950617E+0,
    1.846462673611111111111111E+0,
    -8.912109375000000000000000E-1,
    7.324218750000000000000000E-2
};

static double P4[] = {
    4.669584423426247427983539E+0,
    -1.120700261622299382716049E+1,
    8.789123535156250000000000E+0,
    -2.364086914062500000000000E+0,
    1.121520996093750000000000E-1
};

static double P5[] = {
    -2.8212072558200244877E1,
    8.4636217674600734632E1,
    -9.1818241543240017361E1,
    4.2534998745388454861E1,
    -7.3687943594796316964E0,
    2.27108001708984375E-1
};

static double P6[] = {
    2.1257013003921712286E2,
    -7.6525246814118164230E2,
    1.0599904525279998779E3,
    -6.9957962737613254123E2,
    2.1819051174421159048E2,
    -2.6491430486951555525E1,
    5.7250142097473144531E-1
};

static double P7[] = {
    -1.9194576623184069963E3,
    8.0617221817373093845E3,
    -1.3586550006434137439E4,
    1.1655393336864533248E4,
    -5.3056469786134031084E3,
    1.2009029132163524628E3,
    -1.0809091978839465550E2,
    1.7277275025844573975E0
};


static double jnx(double n, double x)
{
    double zeta, sqz, zz, zp, np;
    double cbn, n23, t, z, sz;
    double pp, qq, z32i, zzi;
    double ak, bk, akl, bkl;
    int sign, doa, dob, nflg, k, s, tk, tkp1, m;
    static double u[8];
    static double ai, aip, bi, bip;

    /* Test for x very close to n. Use expansion for transition region if so. */
    cbn = cbrt(n);
    z = (x - n) / cbn;
    if (fabs(z) <= 0.7)
	return (jnt(n, x));

    z = x / n;
    zz = 1.0 - z * z;
    if (zz == 0.0)
	return (0.0);

    if (zz > 0.0) {
	sz = sqrt(zz);
	t = 1.5 * (log((1.0 + sz) / z) - sz);	/* zeta ** 3/2          */
	zeta = cbrt(t * t);
	nflg = 1;
    }
    else {
	sz = sqrt(-zz);
	t = 1.5 * (sz - acos(1.0 / z));
	zeta = -cbrt(t * t);
	nflg = -1;
    }
    z32i = fabs(1.0 / t);
    sqz = cbrt(t);

    /* Airy function */
    n23 = cbrt(n * n);
    t = n23 * zeta;

#if CEPHES_DEBUG
    printf("zeta %.5E, Airy(%.5E)\n", zeta, t);
#endif
    airy(t, &ai, &aip, &bi, &bip);

    /* polynomials in expansion */
    u[0] = 1.0;
    zzi = 1.0 / zz;
    u[1] = polevl(zzi, P1, 1) / sz;
    u[2] = polevl(zzi, P2, 2) / zz;
    u[3] = polevl(zzi, P3, 3) / (sz * zz);
    pp = zz * zz;
    u[4] = polevl(zzi, P4, 4) / pp;
    u[5] = polevl(zzi, P5, 5) / (pp * sz);
    pp *= zz;
    u[6] = polevl(zzi, P6, 6) / pp;
    u[7] = polevl(zzi, P7, 7) / (pp * sz);

#if CEPHES_DEBUG
    for (k = 0; k <= 7; k++)
	printf("u[%d] = %.5E\n", k, u[k]);
#endif

    pp = 0.0;
    qq = 0.0;
    np = 1.0;
    /* flags to stop when terms get larger */
    doa = 1;
    dob = 1;
    akl = NPY_INFINITY;
    bkl = NPY_INFINITY;

    for (k = 0; k <= 3; k++) {
	tk = 2 * k;
	tkp1 = tk + 1;
	zp = 1.0;
	ak = 0.0;
	bk = 0.0;
	for (s = 0; s <= tk; s++) {
	    if (doa) {
		if ((s & 3) > 1)
		    sign = nflg;
		else
		    sign = 1;
		ak += sign * mu[s] * zp * u[tk - s];
	    }

	    if (dob) {
		m = tkp1 - s;
		if (((m + 1) & 3) > 1)
		    sign = nflg;
		else
		    sign = 1;
		bk += sign * lambda[s] * zp * u[m];
	    }
	    zp *= z32i;
	}

	if (doa) {
	    ak *= np;
	    t = fabs(ak);
	    if (t < akl) {
		akl = t;
		pp += ak;
	    }
	    else
		doa = 0;
	}

	if (dob) {
	    bk += lambda[tkp1] * zp * u[0];
	    bk *= -np / sqz;
	    t = fabs(bk);
	    if (t < bkl) {
		bkl = t;
		qq += bk;
	    }
	    else
		dob = 0;
	}
#if CEPHES_DEBUG
	printf("a[%d] %.5E, b[%d] %.5E\n", k, ak, k, bk);
#endif
	if (np < MACHEP)
	    break;
	np /= n * n;
    }

    /* normalizing factor ( 4*zeta/(1 - z**2) )**1/4    */
    t = 4.0 * zeta / zz;
    t = sqrt(sqrt(t));

    t *= ai * pp / cbrt(n) + aip * qq / (n23 * n);
    return (t);
}

/* Asymptotic expansion for transition region,
 * n large and x close to n.
 * AMS55 #9.3.23.
 */

static double PF2[] = {
    -9.0000000000000000000e-2,
    8.5714285714285714286e-2
};

static double PF3[] = {
    1.3671428571428571429e-1,
    -5.4920634920634920635e-2,
    -4.4444444444444444444e-3
};

static double PF4[] = {
    1.3500000000000000000e-3,
    -1.6036054421768707483e-1,
    4.2590187590187590188e-2,
    2.7330447330447330447e-3
};

static double PG1[] = {
    -2.4285714285714285714e-1,
    1.4285714285714285714e-2
};

static double PG2[] = {
    -9.0000000000000000000e-3,
    1.9396825396825396825e-1,
    -1.1746031746031746032e-2
};

static double PG3[] = {
    1.9607142857142857143e-2,
    -1.5983694083694083694e-1,
    6.3838383838383838384e-3
};


static double jnt(double n, double x)
{
    double z, zz, z3;
    double cbn, n23, cbtwo;
    double ai, aip, bi, bip;	/* Airy functions */
    double nk, fk, gk, pp, qq;
    double F[5], G[4];
    int k;

    cbn = cbrt(n);
    z = (x - n) / cbn;
    cbtwo = cbrt(2.0);

    /* Airy function */
    zz = -cbtwo * z;
    airy(zz, &ai, &aip, &bi, &bip);

    /* polynomials in expansion */
    zz = z * z;
    z3 = zz * z;
    F[0] = 1.0;
    F[1] = -z / 5.0;
    F[2] = polevl(z3, PF2, 1) * zz;
    F[3] = polevl(z3, PF3, 2);
    F[4] = polevl(z3, PF4, 3) * z;
    G[0] = 0.3 * zz;
    G[1] = polevl(z3, PG1, 1);
    G[2] = polevl(z3, PG2, 2) * z;
    G[3] = polevl(z3, PG3, 2) * zz;
#if CEPHES_DEBUG
    for (k = 0; k <= 4; k++)
	printf("F[%d] = %.5E\n", k, F[k]);
    for (k = 0; k <= 3; k++)
	printf("G[%d] = %.5E\n", k, G[k]);
#endif
    pp = 0.0;
    qq = 0.0;
    nk = 1.0;
    n23 = cbrt(n * n);

    for (k = 0; k <= 4; k++) {
	fk = F[k] * nk;
	pp += fk;
	if (k != 4) {
	    gk = G[k] * nk;
	    qq += gk;
	}
#if CEPHES_DEBUG
	printf("fk[%d] %.5E, gk[%d] %.5E\n", k, fk, k, gk);
#endif
	nk /= n23;
    }

    fk = cbtwo * ai * pp / cbn + cbrt(4.0) * aip * qq / n;
    return (fk);
}
