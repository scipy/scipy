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
 * nCv is the number of combinations of n things taken v at a time.
 * "e" here is being used as an abbreviation of "epsilon", not of e=2.71828...
 */


#include "mconf.h"
#include "float.h"
extern double MAXLOG;

#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

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

/* Now the Smirnov functions for one-sided testing. */

/* For n bigger than 1028, nCv is no longer always representable in
   64bit floating point. Back off a bit further to account for
   the other factors in each term.
   */
#define MAX_N_NOLOG 1013

/* Exact Smirnov statistic, for one-sided test.  */
double smirnov(int n, double e)
{
    int v, nn;
    double p;

    /* This comparison should assure returning NaN whenever
     * e is NaN itself.  In original || form it would proceed */
    if (!(n > 0 && e >= 0.0 && e <= 1.0))
	return (NPY_NAN);
    if (e == 0.0)
	return 1.0;
    if (e == 1.0)
	return 0.0;

    nn = (int) (floor((double) n * (1.0 - e)));
    nn = MIN(nn, n-1);
    p = 0.0;
    if (n < MAX_N_NOLOG) {
	double c = 1.0;
	double p0 = pow(1-e, (double)(n));
	c *= n;
	for (v = 1; v <= nn; v++) {
	    double evn = e + ((double) v) / n;
	    double omevn = 1.0-evn;
	    double t;
	    if (omevn <= 0) {
		continue;
	    }
	    t = c * pow(evn, (double)(v - 1)) * pow(omevn, (double)(n - v));
	    if (t == 0) {
		/* Can happen if one of the two pow()'s above is very small so 0 is returned.
		  Even though the full product would be non-zero, the order matters.
		  Redo but include the pow()s half at a time.*/
		double v2 = ((double)(v-1))/2;
		int nmv2 = ((double)(n-v))/2;
		double pwhalf = pow(evn, (double)v2) * pow(omevn, (double)nmv2);
		if (pwhalf != 0) {
		    t = (pwhalf * c) * pwhalf;
		}
	    }
	    p += t;
	    /* Next combinatorial term; worst case error = 4e-15.  */
	    c *= ((double) (n - v)) / (v + 1);
	}
	p = p*e + p0;
	return p;
    }
    else {
	double lgamnp1 = lgam((double) (n + 1));
	double lt, lc = 0.0;
	for (v = 0; v <= nn; v++) {
	    double evn = e + ((double) v) / n;
	    if (1-evn <= 0) {
		continue;
	    }
	    lc = lgamnp1 - lgam((double)(v + 1)) - lgam((double)(n - v + 1));
	    lt = lc + (v - 1) * log(evn)  + (n - v) * log1p(-evn);
	    if (lt >= -MAXLOG) {
		p += exp(lt);
	    }
	}
	p *= e;
    }
    return p;
}


/* Derivative of smirnov(n, x)

    smirnov(n, x) = sum f_v(x)  over {v : 0 <= v <= (n-1)*x}
    f_v  = nCv * x * pow(xvn, v-1) * pow(omxvn, n-v)     [but special case v=0,1,n below]
    f_v' = nCv *     pow(xvn, v-2) * pow(omxvn, n-v-1) * (xvn*omxvn + (v-1)*x*omxvn -(n-v)*x*xvn)
      where xvn = 1+v/n, omxvn = 1-xvn = 1-x-v/n
    Specializations for v=0,1,n
    f_0 = pow(omxvn, n)
    f_1 = n * x * pow(omxvn, n-1)
    f_n = x * pow(xvn, n-1)
    f_v(x) has a zero of order (n-v) at x=(1-v/n) for v=0,1,...n-1;
      and f_n(x) has a simple zero at x=0
    As the zero of f_(n-1) at x=1/n is only a simple zero, there is a
      discontinuity in smirnovp at x=1/n, and a fixable one at x=0.
    smirnovp is continuous at all other values of x.

    Note that the terms in these products can have wildly different magnitudes.
    If one of the pow()s in the product is "zero", recalculate the product
    in pieces and different order to avoid underflow.
*/

/* Use for small values of n, where nCv is representable as a double */
static double smirnovp_pow(int n, double x)
{
    int v, nn;
    double xvn, omxvn, pp, t;
    double c = 1.0; /* Holds the nCv combinations value */

    assert(n > 1);
    assert(x > 0);
    assert(x < 1);

    nn = (int)floor(n * (1.0 - x));
    /* The f_n term is not needed, and would introduce an inconvenient discontuity at x = 0.
    Should only appear if x == 0 (or x is less than DBL_EPSILON, so that 1-x==1) */
    nn = MIN(nn, n-1);
    pp = 0.0;

    /* Pull out v=0 */
    pp = c * pow(1-x, n-1) * (-n);
    c *= n;
    for (v = 1; v <= nn; v++) {
	xvn = x + ((double) v) / n;
	omxvn = 1.0 - xvn;
	if (omxvn <= 0.0) {  /* rounding issue */
	    continue;
	}
	if (v == 1) {
	    double coeff = (omxvn - (n-1)*x);
	    if (coeff != 0) {
		t = c * pow(omxvn, n-2) * coeff;
		if (t == 0) {
		    double nmv2 = (n-v-1)/2.0;
		    double powhf = pow(omxvn, nmv2);
		    if (powhf != 0) {
			t = (powhf * c) * coeff * powhf;
		    }
		}
		pp += t;
	    }
	}
	else {
	    double coeff = (xvn * omxvn + (v-1)*x*omxvn - (n-v) * x * xvn);
	    if (coeff != 0) {
		t = c * pow(omxvn, n - v -1) * pow(xvn, v - 2) * coeff;
		if (t == 0) {
		    double v2 = (v-2)/2.0;
		    double nmv2 = (n-v-1)/2.0;
		    double pwhalf = pow(xvn, v2) * pow(omxvn, nmv2);
		    if (pwhalf != 0) {
			t = (pwhalf * c) * coeff * pwhalf;
		    }
		}
		pp += t;
	    }
	}
	/* Next combinatorial term; worst case error = 4e-15.  */
	c *= ((double) (n - v)) / (v + 1);
    }
    return pp;
}

/* Needed for larger values of n.
   nCv is not always representable as a double if N >= 1028.
   Use log(gamma(m+1)) and exponentiate later */
static double smirnovp_gamma(int n, double x)
{
    int v, nn;
    double xvn, omxvn, pp, t;

    assert(n > 1);
    assert(x > 0);
    assert(x < 1);

    nn = (int)floor(n * (1.0 - x));
    /* The f_n term is not needed, and would introduce an additional
        discontinuity at x = 0, so drop it.
       nn=n  <==>  x == 0 (or x is less than DBL_EPSILON, so that 1-x=1) */
    nn = MIN(nn, n-1);
    pp = 0.0;
    {
	double logn = log((double) (n));
	double lgamnp1 = lgam((double) (n + 1));
	/* Pull out v=0 */
	t = (n-1) * log1p(-x) + logn;
	if (t > -MAXLOG) {
	    pp = -exp(t);
	}

	for (v = 1; v <= nn; v++) {
	    int signcoeff = 1;
	    xvn = x + ((double) v) / n;
	    omxvn = 1.0 - xvn;
	    if (fabs(omxvn) <= 0.0) {
		continue;
	    }
	    if (v == 1) {
		/* f_1' = nC1 * pow(omxvn, n-2) * (omxvn-(n-1)x)*/
		double coeff = (omxvn - (n-1)*x);
		if (coeff == 0) {
		    continue;
		}
		if (coeff < 0) {
		    signcoeff = -1;
		    coeff = -coeff;
		}
		t = (n-2) * log1p(-xvn) + log(coeff)+ logn;
	    }
	    else {
		/*  f_v' = nCv * pow(xvn, v-2) * pow(omxvn, n-v-1) *
			(xvn*omxvn + (v-1)*x*omxvn -(n-v)*x*xvn)   */
		double lc;
		double coeff = (xvn * omxvn) + (v-1)*x*omxvn - (n-v) * x * xvn;
		if (coeff == 0) {
		    continue;
		}
		if (coeff < 0) {
		    signcoeff = -1;
		    coeff = -coeff;
		}
		lc = (lgamnp1 - lgam((double)(v+1)) - lgam((double)(n-v+1)));
		t = (v-2) * log(xvn) + (n-v-1) * log1p(-xvn);
		t += lc + log(coeff);
	    }
	    if (t > -MAXLOG) {
		pp += signcoeff * exp(t);
	    }
	}
    }
    return pp;
}


/* Derivative of smirnov(n, x)
   One point of discontinuity at x=1/n
*/
double smirnovp(int n, double x)
{
    /* This comparison should assure returning NaN whenever
     * x is NaN itself.  */
    if (!(n > 0 && x >= 0.0 && x <= 1.0)) {
	return (NPY_NAN);
    }
    if (n == 1) {
	 /* Slope is always -1 for n ==1, even at x = 1.0 */
	 return -1.0;
    }
    if (x == 1.0) {
	return -0.0;
    }
    /* If x is 0, the derivative is discontinuous, but approaching
       from the right the limit is -1 */
    if (x == 0.0) {
	return -1.0;
    }
    /* If n is too big, calculate using logs */
    if (n < MAX_N_NOLOG) {
	return smirnovp_pow(n, x);
    }
    return smirnovp_gamma(n, x);
}


static double _xtol = 2e-14;
static double _rtol = 4*DBL_EPSILON;

static int _within_tol(double x, double y, double atol, double rtol)
{
    double diff = fabs(x-y);
    int result = (diff <=  (atol + rtol * fabs(y)));
    return result;
}

#define MAX_BISECTION_STEPS 10
#define MIN_N_ADJUST_ESTIMATE 5

/* Solve smirnov(n, x) == p for p in [0, 1] with constraint 0<=x<=1.

    Discussion:
    smirnov(n, _) is easily invertible except for
    values of p close to the endpoints, where it can be quite brutal.
    Useful Asymptotic formula:
	As n -> infinity,  smirnov(n, x*sqrt(n)) ~ exp(-2 * x**2).
    Using this approximation to generate a starting point, and then
    approximating the derivative with the derivative of the limit
    sometimes (often?) succeeds, but also often results in
    slow convergence or non-convergence (indicated by a return of NPY_NAN.)
    Either the original estimate isn't a good starting point
    for Newton-Raphson(NR), or the derivative of the asymptote isn't close
    enough to the actual derivative.

    Current Algorithm:
    1. First handle the two endpoints: p=0, or p=1.
    2. Handle n==1: smirnov(1,x) = 1-x
    3. Exactly handle the case of extremely small p with a formula.
       (corresponds to (n-1)/n < x < 1.0).
    4a. Generate an initial estimate for x using the asymptotic limit.
    4b. If indications are that the initial estimate is not
        going to be good enough, use a few bisection steps
        to get a better initial estimate, which, importantly, is
        on the "correct" side of the true value.
    4c. Use Newton-Raphson to find the root, ensuring estimates stay within the interval [0,1].
        [Instrumentation suggest it converges in ~6-10 iterations for n<=100.]

    Note on 3.
    Small p (extremely small!):
     smirnov(n, x) = (1-x)**n,     if n*(1-x) <= 1, so
     smirnovi(n, p) = 1-p**(1/n),  if p is small enough, log(p) < -n*log(n).
    Asymptotically this cutoff handles fewer and fewer values of p,
    but smirnov(n, x) is not smooth at (1-1/n), so it is still useful.

    Note on 4b.  The estimate is always > correct value of x.
	 This can present a big problem if the true value lies on
	 the "sloped" part of the graph, but the initial estimate lies on
	  the "flat" part of the graph, (or if the estimate is bigger than 1.)
     For x about 1/sqrt(n), should be OK.
     For x >> 1/sqrt(n), then the problems appear.

    Alternative approaches considered.
    1. Only use bisection.
       Pros: Convergence guaranteed.
       Cons: Many iterations required to be within same tolerance.
    2. Use a smarter bracketing algorithm, such as a modified false position or brentq.
       Pros: Faster convergence than bisection.
       Cons: Exists elsewhere in scipy (scipy.optimize) but not easily
	     callable from here, so an implementation needed.
       [This was tested from the Python side, and brentq had fewer correct digits with
       the same number of iterations than did this combination of bisection & NR.]

    Using the derivative results in needing fewer iterations
    to achieve same (or better) accuracy.
*/

/* Functional inverse of Smirnov distribution
 * finds x such that smirnov(n, x) = p.  */
double smirnovi(int n, double p)
{

    double x, logp;
    int iterations = 0;
    int function_calls = 0;
    double sqrtn;

    if (!(n > 0 && p >= 0.0 && p <= 1.0)) {
	mtherr("smirnovi", DOMAIN);
	return (NPY_NAN);
    }
    if (p == 0.0) {
	return 1.0;
    } else if (p == 1.0) {
	return 0.0;
    } else if (n == 1) {
	return 1-p;
    }

    /* Handle p *very* close to 0.  Correspond to (n-1)/n < x < 1  */
    logp = log(p);
    if (logp < -n * log(n)) {
	x = - expm1(logp / n);
	return x;
    }

    sqrtn = sqrt(n);
    /* Start with approximation of x from smirnov(n, sqrt(n)*x) ~ exp(-2 n x^2) */
    x = sqrt(-logp / (2.0 * n));
    x = MIN(x, (n-1.0)/n);

    if (n >= MIN_N_ADJUST_ESTIMATE && x > 1/sqrtn) {
	/* Do some bisections to get an estimate close enough to use NR.
	 If final is [low, high]:
	   a) If high < sqrt(n), use x = high/n.
	   b) If low > sqrt(n), use x = low/n
	   c) Otherwise use x = (low+high)/2/n

	 On the interval [j/n, (j+1)/n], smirnov(n, x) is a poly of degree n
	 and NR will converge quickly after it gets close enough.
	 */
	int high = (int)ceil(n * x);
	int low = 0;
	high = MIN(high, n-1);
	if (x >= 2/sqrtn) {
	    low = (int)floor(sqrtn);
	}
	assert(smirnov(n, low*1.0/n) >= p);
	function_calls++;
	do {
	    int idx = (int)((low+high)/2);
	    double xidx = idx*1.0/n;
	    double fidx = smirnov(n, xidx);
	    ++function_calls;
	    /* If fidx is NAN, that is actually ok.
	       Just means we need to start further to the left. */
	    if (npy_isnan(fidx) || fidx < p) {
		high = idx;
	    } else if (fidx == p) {
		return xidx;
	    } else if (fidx > p) {
		assert (fidx > p);
		low = idx;
	    }
	    if (++iterations >= MAX_BISECTION_STEPS) {
		break;
	    }
	} while (low < high - 1);
	x = (low > sqrtn ? low : (high < sqrtn ? high : (low+high)/2.0));
	x /= n;
    }
    assert (x < 1);

    /* smirnov should be well-enough behaved for NR starting at this location */
    do {
	double deriv, diff, x0 = x, deltax, val;
	assert(x < 1);
	assert(x > 0);

	val = smirnov(n, x0);
	++function_calls;
	diff = val - p;
	if (diff == 0) {
	    return x;
	}
	deriv = smirnovp(n, x0);
	++function_calls;
	if (deriv == 0) {
	    /* x was not within tolerance, but now we hit a 0 derivative.
	    This implies that x >> 1/sqrt(n), and even then |smirnovp| >= |smirnov|
	    so this condition is unexpected. */
	    mtherr("smirnovi", UNDERFLOW);
	    return (NPY_NAN);
	}

	deltax = diff / deriv;
	x = x0 - deltax;
	/* Check out-of-bounds.
	Not expecting this to happen --- smirnov is convex near x=1 and
	  concave near x=0, and we should be approaching from
	  the "correct" side.
	  If out-of-bounds, replace x with a value halfway
	  to the endpoint. */
	if (x <= 0) {
	    x = x0 / 2;
	    if (x <= 0.0) {
		return 0;
	    }
	} else if (x >= 1) {
	    x = (x0+1)/2;
	    if (x >= 1.0) {
		/* Could happen if true value of x lies in interval [1-DBL_EPSILON, 1]
		Requires p**(1/n) < DBL_EPSILON, which is uncommon. */
		return 1.0;
	    }
	}
	/* Note that if p is close to 1, f(x) -> 1, f'(x) -> -1.
	   => abs difference |x-x0| is approx |f(x)-p| >= DBL_EPSILON,
	   => |x-x0|/x >= DBL_EPSILON/x.
	   => cannot use a purely relative criteria as it will fail for x close to 0.
	*/
	if (_within_tol(x, x0, _xtol, _rtol)) {
	    break;
	}
	if (++iterations > MAXITER) {
	    mtherr("smirnovi", TOOMANY);
	    return (x);
	}
    } while (1);
    return x;
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
