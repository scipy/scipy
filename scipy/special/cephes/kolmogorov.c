/* File altered for inclusion in cephes module for Python:
 * Main loop commented out.... */
/*  Travis Oliphant Nov. 1998 */


/* Re Kolmogorov statistics, here is Birnbaum and Tingey's formula for the
 * distribution of D+, the maximum of all positive deviations between a
 * theoretical distribution function P(x) and an empirical one Sn(x)
 * from n samples.
 * 
 * +
 * D  =         sup        [ P(x) - Sn(x) ]
 * n     -inf < x < inf
 * 
 * 
 * [n(1-e)]
 * +            -                    v-1              n-v
 * Pr{D   > e} =    >    C    e (e + v/n)    (1 - e - v/n)
 * n            -   n v
 * v=0
 * 
 * [n(1-e)] is the largest integer not exceeding n(1-e).
 * nCv is the number of combinations of n things taken v at a time.  */


#include "mconf.h"
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
 * probability that sqrt(n) * max deviation > y,
 * or that max deviation > y/sqrt(n).
 * The approximation is useful for the tail of the distribution
 * when n is large.  */
double kolmogorov(double y)
{
    double p, t, r, sign, x;

    if (y < 1.1e-16)
	return 1.0;
    x = -2.0 * y * y;
    sign = 1.0;
    p = 0.0;
    r = 1.0;
    do {
	t = exp(x * r * r);
	p += sign * t;
	if (t == 0.0)
	    break;
	r += 1.0;
	sign = -sign;
    }
    while ((t / p) > 1.1e-16);
    return (p + p);
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
 * Finds y such that kolmogorov(y) = p.
 * If e = smirnovi (n,p), then kolmogi(2 * p) / sqrt(n) should
 * be close to e.  */
double kolmogi(double p)
{
    double y, t, dpdy;
    int iterations;

    if (!(p > 0.0 && p <= 1.0)) {
	mtherr("kolmogi", DOMAIN);
	return (NPY_NAN);
    }
    if ((1.0 - p) < 1e-16)
	return 0.0;
    /* Start with approximation p = 2 exp(-2 y^2).  */
    y = sqrt(-0.5 * log(0.5 * p));
    iterations = 0;
    do {
	/* Use approximate derivative in Newton iteration. */
	t = -2.0 * y;
	dpdy = 4.0 * t * exp(t * y);
	if (fabs(dpdy) > 0.0)
	    t = (p - kolmogorov(y)) / dpdy;
	else {
	    mtherr("kolmogi", UNDERFLOW);
	    return 0.0;
	}
	y = y + t;
	if (++iterations > MAXITER) {
	    mtherr("kolmogi", TOOMANY);
	    return (y);
	}
    }
    while (fabs(t / y) > 1.0e-10);
    return (y);
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
