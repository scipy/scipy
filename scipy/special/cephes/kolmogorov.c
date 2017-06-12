/* File altered for inclusion in cephes module for Python:
 * Main loop commented out.... */
/*  Travis Oliphant Nov. 1998 */


/* Re Kolmogorov statistics, here is Birnbaum and Tingey's formula for the
 * distribution of D+, the maximum of all positive deviations between a
 * theoretical distribution function P(x) and an empirical one Sn(x)
 * from n samples.
 * 
 *     +
 *    D  =         sup     [P(x) - S (x)]
 *     n     -inf < x < inf         n
 *
 *
 *                  [n(1-e)]
 *        +            -                    v-1              n-v
 *    Pr{D   > e} =    >    C    e (e + v/n)    (1 - e - v/n)
 *        n            -   n v
 *                    v=0
 * 
 * [n(1-e)] is the largest integer not exceeding n(1-e).
 * nCv is the number of combinations of n things taken v at a time.  */


#include "mconf.h"
#include "float.h"
extern double MAXLOG;

/* Exact Smirnov statistic, for one-sided test.  */
double smirnov(int n, double e)
{
    int v, nn;
    double evn, omevn, p, t, c, lgamnp1;

    /* This comparison should assure returning NaN whenever
     * e is NaN itself.  In original || form it would proceed */
    if (!(n > 0 && e >= 0.0 && e <= 1.0))
	return (NPY_NAN);
    if (e == 0.0)
	return 1.0;
    nn = (int) (floor((double) n * (1.0 - e)));
    p = 0.0;
    if (n < 1013) {
	c = 1.0;
	for (v = 0; v <= nn; v++) {
	    evn = e + ((double) v) / n;
	    p += c * pow(evn, (double) (v - 1))
		* pow(1.0 - evn, (double) (n - v));
	    /* Next combinatorial term; worst case error = 4e-15.  */
	    c *= ((double) (n - v)) / (v + 1);
	}
    }
    else {
	lgamnp1 = lgam((double) (n + 1));
	for (v = 0; v <= nn; v++) {
	    evn = e + ((double) v) / n;
	    omevn = 1.0 - evn;
	    if (fabs(omevn) > 0.0) {
		t = lgamnp1 - lgam((double) (v + 1))
		    - lgam((double) (n - v + 1))
		    + (v - 1) * log(evn)
		    + (n - v) * log(omevn);
		if (t > -MAXLOG)
		    p += exp(t);
	    }
	}
    }
    return (p * e);
}


/* Kolmogorov's limiting distribution of two-sided test, returns
 * probability that sqrt(n) * max deviation > x,
 * or that max deviation > x/sqrt(n).
 * The approximation is useful for the tail of the distribution
 * when n is large.  */

 /* Two series for kolmogorov(x), a Jacobi theta function
  *  sum (-1)^k exp(-2k^2 x^2) (over all integer k); or
  *  sqrt(2pi)/x * sum exp((2k-1)^2pi^2/(8x^2)) (over positive integer k)
 *  The first is good for x not close to 0, the second for x close to 0
 */
#define X_MIN_USE_ORIGINAL 0.5
#define KOLMOG_RTOL (DBL_EPSILON)

double kolmogorov(double x)
{
    double p, t;

    if (x <= 0)
	return 1.0;

    if (x >= X_MIN_USE_ORIGINAL) {
	double alpha = -2.0 * x * x;
	double sign = 1.0;
	double r = 1.0;
	p = 0.0;
	do {
	    t = exp(alpha * r * r);
	    p += sign * t;
	    if (t == 0.0)
		break;
	    r += 1.0;
	    sign = -sign;
	}
	while ((t / p) > KOLMOG_RTOL);
	p = 2*p;
    }  else {
	double alpha = - NPY_PI * NPY_PI / (8 * x * x);
	double r = 1;
	p = 0.0;
	do {
	    t = exp(alpha * r * r);
	    p += t;
	    if (fabs(t) == 0.0)
		break;
	    r +=  2;
	} while ((t / p) >= KOLMOG_RTOL);
	p *= sqrt(2 * NPY_PI) / x;
	p = 1 - p;
    }
    return p;
}


double kolmogorovp(double x)
{
    double pp, t;

    if (x <= 0)
	return 0.0;

    if (x >= X_MIN_USE_ORIGINAL) {
	double alpha = -2.0 * x * x;
	double sign = 1.0;
	double r = 1.0;
	pp = 0.0;
	do {
	    double r2 = r*r;
	    t = exp(alpha * r2);
	    if (t == 0.0)
		break;
	    pp += sign * t * r2;
	    r += 1.0;
	    sign = -sign;
	}  while ((t / pp) > KOLMOG_RTOL);
	pp = -8 * pp;
    } else {
	double alpha = - NPY_PI * NPY_PI / (8 * x * x);
	double r = 1;
	double pp1 = 0.0;
	double sqrt2pi = sqrt(2 * NPY_PI);
	pp = 0.0;
	do {
	    double r2 = r*r;
	    double q2n = exp(alpha * r2);
	    t = r2 * q2n;
	    pp += t;
	    pp1 += q2n;
	    if (t == 0.0)
		break;
	    r +=  2;
	} while ((t / pp) >= KOLMOG_RTOL);
	pp1 *= sqrt2pi/x/x;
	pp *= pow(NPY_PI, 2) * sqrt2pi / pow(x, 4) / 4;
	pp = -pp + pp1;
    }
    return pp;
}

/* Functional inverse of Smirnov distribution
 * finds e such that smirnov(n,e) = p.  */
double smirnovi(int n, double p)
{
    double e, t, dpde;
    int iterations;

    if (!(p > 0.0 && p <= 1.0)) {
	mtherr("smirnovi", DOMAIN);
	return (NPY_NAN);
    }
    /* Start with approximation p = exp(-2 n e^2).  */
    e = sqrt(-log(p) / (2.0 * n));
    iterations = 0;
    do {
	/* Use approximate derivative in Newton iteration. */
	t = -2.0 * n * e;
	dpde = 2.0 * t * exp(t * e);
	if (fabs(dpde) > 0.0)
	    t = (p - smirnov(n, e)) / dpde;
	else {
	    mtherr("smirnovi", UNDERFLOW);
	    return 0.0;
	}
	e = e + t;
	if (e >= 1.0 || e <= 0.0) {
	    mtherr("smirnovi", OVERFLOW);
	    return 0.0;
	}
	if (++iterations > MAXITER) {
	    mtherr("smirnovi", TOOMANY);
	    return (e);
	}
    }
    while (fabs(t / e) > 1e-10);
    return (e);
}


/* Functional inverse of Kolmogorov statistic for two-sided test.
 * Finds x such that kolmogorov(x) = p.
 * If x = smirnovi (n, p), then kolmogi(2 * p) / sqrt(n) should
 * be close to x.  */
double kolmogi(double p)
{
    double x, t;
    int iterations;

    if (!(p >= 0.0 && p <= 1.0)) {
	mtherr("kolmogi", DOMAIN);
	return (NPY_NAN);
    }
    if ((1.0 - p) < 1e-16)
	return 0.0;

    /* For x between 0.5 and 1, kolmogorov(x) is close to the straight line
     connecting (0.5, 1) to (1.0, 0.25). I.e. p ~ (-6x+7)/4.
     Otherwise use the approximation p ~ 2 exp(-2x^2) */
    if (p > 0.25) {
	x = (7-4*p)/6.0;
    } else {
	x = sqrt(-0.5 * log(0.5 * p));
    }

    iterations = 0;
    do {
	double x0 = x;
	double val = kolmogorov(x0);
	double df = val - p;
	double dpdy;
	if (fabs(df) == 0) {
	    break;
	}
	dpdy = kolmogorovp(x0);
	if (fabs(dpdy) <= 0.0) {
	    mtherr("kolmogi", UNDERFLOW);
	    return 0.0;
	}
	t = df/dpdy;
	x = x0 - t;

	if (fabs(t/x) < KOLMOG_RTOL) {
	break;
	}

	if (++iterations > MAXITER) {
	    mtherr("kolmogi", TOOMANY);
	    break;
	}
    } while(1);
    // while (fabs(t / x) > 1.0e-10);
    return (x);
}





/* Type in a number.  */
/* void
 * getnum (s, px)
 * char *s;
 * double *px;
 * {
 * char str[30];
 * 
 * printf (" %s (%.15e) ? ", s, *px);
 * gets (str);
 * if (str[0] == '\0' || str[0] == '\n')
 * return;
 * sscanf (str, "%lf", px);
 * printf ("%.15e\n", *px);
 * }
 */
/* Type in values, get answers.  */
/*
 * void
 * main ()
 * {
 * int n;
 * double e, p, ps, pk, ek, y;
 * 
 * n = 5;
 * e = 0.0;
 * p = 0.1;
 * loop:
 * ps = n;
 * getnum ("n", &ps);
 * n = ps;
 * if (n <= 0)
 * {
 * printf ("? Operator error.\n");
 * goto loop;
 * }
 */
  /*
   * getnum ("e", &e);
   * ps = smirnov (n, e);
   * y = sqrt ((double) n) * e;
   * printf ("y = %.4e\n", y);
   * pk = kolmogorov (y);
   * printf ("Smirnov = %.15e, Kolmogorov/2 = %.15e\n", ps, pk / 2.0);
   */
/*
 * getnum ("p", &p);
 * e = smirnovi (n, p);
 * printf ("Smirnov e = %.15e\n", e);
 * y = kolmogi (2.0 * p);
 * ek = y / sqrt ((double) n);
 * printf ("Kolmogorov e = %.15e\n", ek);
 * goto loop;
 * }
 */
