/*                                                     ndtr.c
 *
 *     Normal distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, ndtr();
 *
 * y = ndtr( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the area under the Gaussian probability density
 * function, integrated from minus infinity to x:
 *
 *                            x
 *                             -
 *                   1        | |          2
 *    ndtr(x)  = ---------    |    exp( - t /2 ) dt
 *               sqrt(2pi)  | |
 *                           -
 *                          -inf.
 *
 *             =  ( 1 + erf(z) ) / 2
 *             =  erfc(z) / 2
 *
 * where z = x/sqrt(2). Computation is via the functions
 * erf and erfc.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -13,0        30000       3.4e-14     6.7e-15
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition         value returned
 * erfc underflow    x > 37.519379347       0.0
 *
 */
/*							erf.c
 *
 *	Error function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, erf();
 *
 * y = erf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * The integral is
 *
 *                           x
 *                            -
 *                 2         | |          2
 *   erf(x)  =  --------     |    exp( - t  ) dt.
 *              sqrt(pi)   | |
 *                          -
 *                           0
 *
 * For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise
 * erf(x) = 1 - erfc(x).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,1         30000       3.7e-16     1.0e-16
 *
 */
/*							erfc.c
 *
 *	Complementary error function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, erfc();
 *
 * y = erfc( x );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *  1 - erf(x) =
 *
 *                           inf.
 *                             -
 *                  2         | |          2
 *   erfc(x)  =  --------     |    exp( - t  ) dt
 *               sqrt(pi)   | |
 *                           -
 *                            x
 *
 *
 * For small x, erfc(x) = 1 - erf(x); otherwise rational
 * approximations are computed.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,26.6417   30000       5.7e-14     1.5e-14
 */


/*
 * Cephes Math Library Release 2.2:  June, 1992
 * Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

#include <float.h>		/* DBL_EPSILON */
#include "mconf.h"

extern double MAXLOG;

static double P[] = {
    2.46196981473530512524E-10,
    5.64189564831068821977E-1,
    7.46321056442269912687E0,
    4.86371970985681366614E1,
    1.96520832956077098242E2,
    5.26445194995477358631E2,
    9.34528527171957607540E2,
    1.02755188689515710272E3,
    5.57535335369399327526E2
};

static double Q[] = {
    /* 1.00000000000000000000E0, */
    1.32281951154744992508E1,
    8.67072140885989742329E1,
    3.54937778887819891062E2,
    9.75708501743205489753E2,
    1.82390916687909736289E3,
    2.24633760818710981792E3,
    1.65666309194161350182E3,
    5.57535340817727675546E2
};

static double R[] = {
    5.64189583547755073984E-1,
    1.27536670759978104416E0,
    5.01905042251180477414E0,
    6.16021097993053585195E0,
    7.40974269950448939160E0,
    2.97886665372100240670E0
};

static double S[] = {
    /* 1.00000000000000000000E0, */
    2.26052863220117276590E0,
    9.39603524938001434673E0,
    1.20489539808096656605E1,
    1.70814450747565897222E1,
    9.60896809063285878198E0,
    3.36907645100081516050E0
};

static double T[] = {
    9.60497373987051638749E0,
    9.00260197203842689217E1,
    2.23200534594684319226E3,
    7.00332514112805075473E3,
    5.55923013010394962768E4
};

static double U[] = {
    /* 1.00000000000000000000E0, */
    3.35617141647503099647E1,
    5.21357949780152679795E2,
    4.59432382970980127987E3,
    2.26290000613890934246E4,
    4.92673942608635921086E4
};

#define UTHRESH 37.519379347


double ndtr(double a)
{
    double x, y, z;

    if (cephes_isnan(a)) {
	sf_error("ndtr", SF_ERROR_DOMAIN, NULL);
	return (NPY_NAN);
    }

    x = a * NPY_SQRT1_2;
    z = fabs(x);

    if (z < NPY_SQRT1_2)
	y = 0.5 + 0.5 * erf(x);

    else {
	y = 0.5 * erfc(z);

	if (x > 0)
	    y = 1.0 - y;
    }

    return (y);
}

//
// Gauss-Legendre points.
//

// Gauss-Legendre n = 7 x coordinates in [-1, 1]
static double glx7[] = {
    -0.9491079123427584,
    -0.7415311855993945,
    -0.4058451513773972,
     0.0,
     0.4058451513773972,
     0.7415311855993945,
     0.9491079123427584
};

// Gauss-Legendre n = 7 weights
static double glw7[] = {
    0.12948496616886998,
    0.27970539148927664,
    0.38183005050511876,
    0.4179591836734691,
    0.38183005050511876,
    0.27970539148927664,
    0.12948496616886998
};

// Gauss-Legendre n = 21 x coordinates in [-1, 1]
static double glx21[] = {
    -0.9937521706203896,
    -0.9672268385663063,
    -0.9200993341504009,
    -0.8533633645833172,
    -0.768439963475678,
    -0.6671388041974123,
    -0.5516188358872198,
    -0.4243421202074388,
    -0.2880213168024011,
    -0.14556185416089512,
     0.0,
     0.14556185416089512,
     0.2880213168024011,
     0.4243421202074388,
     0.5516188358872198,
     0.6671388041974123,
     0.768439963475678,
     0.8533633645833172,
     0.9200993341504009,
     0.9672268385663063,
     0.9937521706203896
};

// Gauss-Legendre n = 21 weights
static double glw21[] = {
    0.016017228257774067,
    0.03695378977085202,
    0.057134425426857,
    0.07610011362837958,
    0.09344442345603363,
    0.10879729916714861,
    0.12183141605372863,
    0.13226893863333766,
    0.13988739479107337,
    0.1445244039899701,
    0.1460811336496905,
    0.1445244039899701,
    0.13988739479107337,
    0.13226893863333766,
    0.12183141605372863,
    0.10879729916714861,
    0.09344442345603363,
    0.07610011362837958,
    0.057134425426857,
    0.03695378977085202,
    0.016017228257774067
};


#define NDTR_DELTA_TAIL_LIMIT  38.5
#define SQRT_2PI  2.5066282746310007

//
// Estimate integral of the normal PDF from a to b using
// Gauss-Legendre quadrature.
//
// n is the length of the arrays x and weight.
// x and weight are the Gauss-Legendre points and weights.
//
double norm_quad_gl(double a, double b, int n, double x[], double weight[])
{
    double hw = (b - a)/2;  // half width of the interval
    double m = (a + b)/2;   // midpoint of the interval
    double delta = 0.0;
    for (int i = 0; i < n; ++i) {
        double z = hw*x[i] + m;
        double pdf = exp(-z*z/2);
        delta += pdf * weight[i];
    }
    return delta * hw / SQRT_2PI;
}

//
// Compute ndtr(b) - ndtr(a) robustly.
//
double ndtr_delta(double a, double b)
{
    int sgn = 1;
    double w;

    if (cephes_isnan(a) || cephes_isnan(b)) {
        sf_error("ndtr_delta", SF_ERROR_DOMAIN, NULL);
        return NPY_NAN;
    }
    if (a == b) {
        return 0.0;
    }
    if (a > b) {
        double tmp = a;
        a = b;
        b = tmp;
        sgn = -1;
    }
    if ((a > NDTR_DELTA_TAIL_LIMIT) || (b < -NDTR_DELTA_TAIL_LIMIT)) {
        // The endpoints of the interval are so far out in the tail
        // that the exact result is smaller than can be represented
        // with double precision.
        return 0.0;
    }

    w = b - a;
    if (w < 1e-10) {
        // Degree 3 Taylor series approximation at the midpoint
        // of [a, b].
        double w = b - a;
        double m = (a + b)/2;
        double pdf = exp(-m*m/2)/SQRT_2PI;
        return sgn * w * pdf * (1 + w*w*(4*m*m - 2)/6);
    }
    else if (w < 1e-7) {
        // 7 point Gauss-Legendre quadrature.
        return sgn * norm_quad_gl(a, b, 7, glx7, glw7);
    }
    else if (w < 1.0) {
        // 21 point Gauss-Legendre quadrature.
        return sgn * norm_quad_gl(a, b, 21, glx21, glw21);
    }
    else {
        // The interval length isn't "small", so compute
        // the result using the difference in ndtr.
        double delta;

        if (a > 0) {
            delta = ndtr(-a) - ndtr(-b);
        }
        else {
            delta = ndtr(b) - ndtr(a);
        }
        if (delta < 0) {
            delta = 0.0;
        }
        return sgn*delta;
    }
}

//
// Compute the CDF for the truncated normal distribution.
//
double trunc_ndtr(double a, double b, double x)
{
    if (cephes_isnan(x) || cephes_isnan(a) || cephes_isnan(b) || (a >= b)) {
        sf_error("trunc_ndtr", SF_ERROR_DOMAIN, NULL);
        return NPY_NAN;
    }
    if (x <= a) {
        return 0.0;
    }
    else if (x >= b) {
        return 1.0;
    }
    else {
        double delta = ndtr_delta(a, b);
        // printf("delta = %16.12e\n", delta);
        if (delta > 0) {
            // printf("using ndtr_delta\n");
            return ndtr_delta(a, x) / delta;
        }
        else {
            // Let F(x) be the CDF of the normal distribution.
            // The function ndtr_delta(a, b), which computes F(b) - F(a),
            // returned 0, so the straightforward calculation
            //   ndtr_delta(a, x) / ndtr(a, b) = (F(x) - F(a))/(F(b) - F(a))
            // won't work.  Instead, we'll express this formula in terms
            // of logarithms of F.  For brevity, define Fx = F(x), Fa = F(a),
            // etc., and logFx = log(F(x)), etc.
            //
            // log((Fx - Fa)/(Fb - Fa))
            //   = log(Fx - Fa) - log(Fb - Fa)
            //   = log(Fx) + log(1 - Fa/Fx) - log(Fb) - log(1 + Fa/Fb)
            //   = logFx + log(1 - exp(logFa - logFx))
            //     - logFb - log(1 - exp(logFa - logFb))
            //   = logFx + log1p(-exp(logFa - logFx))
            //     - logFb - log1p(-exp(logFa - logFb))
            //
            // The truncated normal CDF can also be expressed in terms of
            // the survivial functon S(x) as (S(a) - S(x))/(S(a) - S(b)).
            // Expressing the log of that in terms of the logarithms of the
            // survival function gives
            //
            // log(S(a) - S(x))/(S(a) - S(b))
            //   = log(Sa - Sx) - log(Sa - Sb)
            //   = logSa + log(1 - Sx/Sa) - logSa - log(1 - Sb/Sa)
            //   = log1p(-exp(logSx - logSa)) - log1p(-exp(logSb - logSa))
            //

            double p;
            // printf("using log_ndtr\n");
            if (a < 0) {
                double logFa = log_ndtr(a);
                double logFb = log_ndtr(b);
                double logFx = log_ndtr(x);
                double tab = log1p(-exp(logFa - logFb));
                double tax = log1p(-exp(logFa - logFx));
                double logp = logFx + tax - (logFb + tab);
                //p = expm1(logFx - logFa) / expm1(logFb - logFa);
                p = exp(logp);
            }
            else {
                double logSa = log_ndtr(-a);
                double logSb = log_ndtr(-b);
                double logSx = log_ndtr(-x);
                /*
                logp = log1p(-exp(logSx - logSa))
                       - log1p(-exp(logSb - logSa));
                */
                p = expm1(logSx - logSa) / expm1(logSb - logSa);
            }
            return p;
        }
    }
}


double erfc(double a)
{
    double p, q, x, y, z;

    if (cephes_isnan(a)) {
        sf_error("erfc", SF_ERROR_DOMAIN, NULL);
	return (NPY_NAN);
    }

    if (a < 0.0)
	x = -a;
    else
	x = a;

    if (x < 1.0)
	return (1.0 - erf(a));

    z = -a * a;

    if (z < -MAXLOG) {
      under:
	sf_error("erfc", SF_ERROR_UNDERFLOW, NULL);
	if (a < 0)
	    return (2.0);
	else
	    return (0.0);
    }

    z = exp(z);

    if (x < 8.0) {
	p = polevl(x, P, 8);
	q = p1evl(x, Q, 8);
    }
    else {
	p = polevl(x, R, 5);
	q = p1evl(x, S, 6);
    }
    y = (z * p) / q;

    if (a < 0)
	y = 2.0 - y;

    if (y == 0.0)
	goto under;

    return (y);
}



double erf(double x)
{
    double y, z;

    if (cephes_isnan(x)) {
	sf_error("erf", SF_ERROR_DOMAIN, NULL);
	return (NPY_NAN);
    }

    if (x < 0.0) {
	return -erf(-x);
    }

    if (fabs(x) > 1.0)
	return (1.0 - erfc(x));
    z = x * x;

    y = x * polevl(z, T, 4) / p1evl(z, U, 5);
    return (y);

}

/*
 * double log_ndtr(double a)
 *
 * For a > -20, use the existing ndtr technique and take a log.
 * for a <= -20, we use the Taylor series approximation of erf to compute
 * the log CDF directly. The Taylor series consists of two parts which we will name "left"
 * and "right" accordingly.  The right part involves a summation which we compute until the
 * difference in terms falls below the machine-specific EPSILON.
 *
 * \Phi(z) &=&
 *   \frac{e^{-z^2/2}}{-z\sqrt{2\pi}}  * [1 +  \sum_{n=1}^{N-1}  (-1)^n \frac{(2n-1)!!}{(z^2)^n}]
 *   + O(z^{-2N+2})
 *   = [\mbox{LHS}] * [\mbox{RHS}] + \mbox{error}.
 *
 */

double log_ndtr(double a)
{

    double log_LHS,		/* we compute the left hand side of the approx (LHS) in one shot */
     last_total = 0,		/* variable used to check for convergence */
	right_hand_side = 1,	/* includes first term from the RHS summation */
	numerator = 1,		/* numerator for RHS summand */
	denom_factor = 1,	/* use reciprocal for denominator to avoid division */
	denom_cons = 1.0 / (a * a);	/* the precomputed division we use to adjust the denominator */
    long sign = 1, i = 0;

    if (a > 6) {
	return -ndtr(-a);     /* log(1+x) \approx x */
    }
    if (a > -20) {
	return log(ndtr(a));
    }
    log_LHS = -0.5 * a * a - log(-a) - 0.5 * log(2 * M_PI);

    while (fabs(last_total - right_hand_side) > DBL_EPSILON) {
	i += 1;
	last_total = right_hand_side;
	sign = -sign;
	denom_factor *= denom_cons;
	numerator *= 2 * i - 1;
	right_hand_side += sign * numerator * denom_factor;

    }
    return log_LHS + log(right_hand_side);
}
