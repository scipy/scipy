/*                                                     unity.c
 *
 * Relative error approximations for function arguments near
 * unity.
 *
 *    log1p(x) = log(1+x)
 *    expm1(x) = exp(x) - 1
 *    cosm1(x) = cos(x) - 1
 *
 */

#include "mconf.h"
/* log1p(x) = log(1 + x)  */

/* Coefficients for log(1+x) = x - x**2/2 + x**3 P(x)/Q(x)
 * 1/sqrt(2) <= x < sqrt(2)
 * Theoretical peak relative error = 2.32e-20
 */
static double LP[] = {
    4.5270000862445199635215E-5,
    4.9854102823193375972212E-1,
    6.5787325942061044846969E0,
    2.9911919328553073277375E1,
    6.0949667980987787057556E1,
    5.7112963590585538103336E1,
    2.0039553499201281259648E1,
};

static double LQ[] = {
    /* 1.0000000000000000000000E0, */
    1.5062909083469192043167E1,
    8.3047565967967209469434E1,
    2.2176239823732856465394E2,
    3.0909872225312059774938E2,
    2.1642788614495947685003E2,
    6.0118660497603843919306E1,
};

#define SQRTH 0.70710678118654752440
#define SQRT2 1.41421356237309504880

double log1p(double x)
{
    double z;

    z = 1.0 + x;
    if ((z < SQRTH) || (z > SQRT2))
	return (log(z));
    z = x * x;
    z = -0.5 * z + x * (z * polevl(x, LP, 6) / p1evl(x, LQ, 6));
    return (x + z);
}



/* expm1(x) = exp(x) - 1  */

/*  e^x =  1 + 2x P(x^2)/( Q(x^2) - P(x^2) )
 * -0.5 <= x <= 0.5
 */

static double EP[3] = {
    1.2617719307481059087798E-4,
    3.0299440770744196129956E-2,
    9.9999999999999999991025E-1,
};

static double EQ[4] = {
    3.0019850513866445504159E-6,
    2.5244834034968410419224E-3,
    2.2726554820815502876593E-1,
    2.0000000000000000000897E0,
};

double expm1(double x)
{
    double r, xx;

    if (!cephes_isfinite(x)) {
	if (cephes_isnan(x)) {
	    return x;
	}
	else if (x > 0) {
	    return x;
	}
	else {
	    return -1.0;
	}

    }
    if ((x < -0.5) || (x > 0.5))
	return (exp(x) - 1.0);
    xx = x * x;
    r = x * polevl(xx, EP, 2);
    r = r / (polevl(xx, EQ, 3) - r);
    return (r + r);
}



/* cosm1(x) = cos(x) - 1  */

static double coscof[7] = {
    4.7377507964246204691685E-14,
    -1.1470284843425359765671E-11,
    2.0876754287081521758361E-9,
    -2.7557319214999787979814E-7,
    2.4801587301570552304991E-5,
    -1.3888888888888872993737E-3,
    4.1666666666666666609054E-2,
};

double cosm1(double x)
{
    double xx;

    if ((x < -NPY_PI_4) || (x > NPY_PI_4))
	return (cos(x) - 1.0);
    xx = x * x;
    xx = -0.5 * xx + xx * xx * polevl(xx, coscof, 6);
    return xx;
}
