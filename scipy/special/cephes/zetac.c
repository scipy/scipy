/*                                                     zetac.c
 *
 *     Riemann zeta function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, zetac();
 *
 * y = zetac( x );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *
 *                inf.
 *                 -    -x
 *   zetac(x)  =   >   k   ,   x > 1,
 *                 -
 *                k=2
 *
 * is related to the Riemann zeta function by
 *
 *     Riemann zeta(x) = zetac(x) + 1.
 *
 * Extension of the function definition for x < 1 is implemented.
 * Zero is returned for x > log2(NPY_INFINITY).
 *
 * An overflow error may occur for large negative x, due to the
 * Gamma function in the reflection formula.
 *
 * ACCURACY:
 *
 * Tabulated values have full machine accuracy.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      1,50        10000       9.8e-16            1.3e-16
 *
 *
 */

/*
 * Cephes Math Library Release 2.1:  January, 1989
 * Copyright 1984, 1987, 1989 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

#include "mconf.h"

/* Riemann zeta(x) - 1
 * for integer arguments between 0 and 30.
 */
static double azetac[] = {
    -1.50000000000000000000E0,
    1.70141183460469231730E38,	/* infinity. */
    6.44934066848226436472E-1,
    2.02056903159594285400E-1,
    8.23232337111381915160E-2,
    3.69277551433699263314E-2,
    1.73430619844491397145E-2,
    8.34927738192282683980E-3,
    4.07735619794433937869E-3,
    2.00839282608221441785E-3,
    9.94575127818085337146E-4,
    4.94188604119464558702E-4,
    2.46086553308048298638E-4,
    1.22713347578489146752E-4,
    6.12481350587048292585E-5,
    3.05882363070204935517E-5,
    1.52822594086518717326E-5,
    7.63719763789976227360E-6,
    3.81729326499983985646E-6,
    1.90821271655393892566E-6,
    9.53962033872796113152E-7,
    4.76932986787806463117E-7,
    2.38450502727732990004E-7,
    1.19219925965311073068E-7,
    5.96081890512594796124E-8,
    2.98035035146522801861E-8,
    1.49015548283650412347E-8,
    7.45071178983542949198E-9,
    3.72533402478845705482E-9,
    1.86265972351304900640E-9,
    9.31327432419668182872E-10
};

/* 2**x (1 - 1/x) (zeta(x) - 1) = P(1/x)/Q(1/x), 1 <= x <= 10 */
static double P[9] = {
    5.85746514569725319540E11,
    2.57534127756102572888E11,
    4.87781159567948256438E10,
    5.15399538023885770696E9,
    3.41646073514754094281E8,
    1.60837006880656492731E7,
    5.92785467342109522998E5,
    1.51129169964938823117E4,
    2.01822444485997955865E2,
};

static double Q[8] = {
    /*  1.00000000000000000000E0, */
    3.90497676373371157516E11,
    5.22858235368272161797E10,
    5.64451517271280543351E9,
    3.39006746015350418834E8,
    1.79410371500126453702E7,
    5.66666825131384797029E5,
    1.60382976810944131506E4,
    1.96436237223387314144E2,
};

/* log(zeta(x) - 1 - 2**-x), 10 <= x <= 50 */
static double A[11] = {
    8.70728567484590192539E6,
    1.76506865670346462757E8,
    2.60889506707483264896E10,
    5.29806374009894791647E11,
    2.26888156119238241487E13,
    3.31884402932705083599E14,
    5.13778997975868230192E15,
    -1.98123688133907171455E15,
    -9.92763810039983572356E16,
    7.82905376180870586444E16,
    9.26786275768927717187E16,
};

static double B[10] = {
    /* 1.00000000000000000000E0, */
    -7.92625410563741062861E6,
    -1.60529969932920229676E8,
    -2.37669260975543221788E10,
    -4.80319584350455169857E11,
    -2.07820961754173320170E13,
    -2.96075404507272223680E14,
    -4.86299103694609136686E15,
    5.34589509675789930199E15,
    5.71464111092297631292E16,
    -1.79915597658676556828E16,
};

/* (1-x) (zeta(x) - 1), 0 <= x <= 1 */
static double R[6] = {
    -3.28717474506562731748E-1,
    1.55162528742623950834E1,
    -2.48762831680821954401E2,
    1.01050368053237678329E3,
    1.26726061410235149405E4,
    -1.11578094770515181334E5,
};

static double S[5] = {
    /* 1.00000000000000000000E0, */
    1.95107674914060531512E1,
    3.17710311750646984099E2,
    3.03835500874445748734E3,
    2.03665876435770579345E4,
    7.43853965136767874343E4,
};

#define MAXL2 127

/*
 * Riemann zeta function, minus one
 */
extern double MACHEP;

double zetac(x)
double x;
{
    int i;
    double a, b, s, w;

    if (x < 0.0) {
	if (x < -30.8148) {
	    mtherr("zetac", OVERFLOW);
	    return (0.0);
	}
	s = 1.0 - x;
	w = zetac(s);
	b = sin(0.5 * NPY_PI * x) * pow(2.0 * NPY_PI,
				    x) * Gamma(s) * (1.0 + w) / NPY_PI;
	return (b - 1.0);
    }

    if (x >= MAXL2)
	return (0.0);		/* because first term is 2**-x */

    /* Tabulated values for integer argument */
    w = floor(x);
    if (w == x) {
	i = x;
	if (i < 31) {
#ifdef UNK
	    return (azetac[i]);
#else
	    return (*(double *) &azetac[4 * i]);
#endif
	}
    }


    if (x < 1.0) {
	w = 1.0 - x;
	a = polevl(x, R, 5) / (w * p1evl(x, S, 5));
	return (a);
    }

    if (x == 1.0) {
	mtherr("zetac", SING);
	return (NPY_INFINITY);
    }

    if (x <= 10.0) {
	b = pow(2.0, x) * (x - 1.0);
	w = 1.0 / x;
	s = (x * polevl(w, P, 8)) / (b * p1evl(w, Q, 8));
	return (s);
    }

    if (x <= 50.0) {
	b = pow(2.0, -x);
	w = polevl(x, A, 10) / p1evl(x, B, 10);
	w = exp(w) + b;
	return (w);
    }


    /* Basic sum of inverse powers */


    s = 0.0;
    a = 1.0;
    do {
	a += 2.0;
	b = pow(a, -x);
	s += b;
    }
    while (b / s > MACHEP);

    b = pow(2.0, -x);
    s = (s + b) / (1.0 - b);
    return (s);
}
