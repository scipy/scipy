/*                                                     fresnl.c
 *
 *     Fresnel integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, S, C;
 * void fresnl();
 *
 * fresnl( x, _&S, _&C );
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the Fresnel integrals
 *
 *           x
 *           -
 *          | |
 * C(x) =   |   cos(pi/2 t**2) dt,
 *        | |
 *         -
 *          0
 *
 *           x
 *           -
 *          | |
 * S(x) =   |   sin(pi/2 t**2) dt.
 *        | |
 *         -
 *          0
 *
 *
 * The integrals are evaluated by a power series for x < 1.
 * For x >= 1 auxiliary functions f(x) and g(x) are employed
 * such that
 *
 * C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
 * S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )
 *
 *
 *
 * ACCURACY:
 *
 *  Relative error.
 *
 * Arithmetic  function   domain     # trials      peak         rms
 *   IEEE       S(x)      0, 10       10000       2.0e-15     3.2e-16
 *   IEEE       C(x)      0, 10       10000       1.8e-15     3.3e-16
 */

/*
 * Cephes Math Library Release 2.1:  January, 1989
 * Copyright 1984, 1987, 1989 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

#include "mconf.h"

/* S(x) for small x */
static double sn[6] = {
    -2.99181919401019853726E3,
    7.08840045257738576863E5,
    -6.29741486205862506537E7,
    2.54890880573376359104E9,
    -4.42979518059697779103E10,
    3.18016297876567817986E11,
};

static double sd[6] = {
    /* 1.00000000000000000000E0, */
    2.81376268889994315696E2,
    4.55847810806532581675E4,
    5.17343888770096400730E6,
    4.19320245898111231129E8,
    2.24411795645340920940E10,
    6.07366389490084639049E11,
};

/* C(x) for small x */
static double cn[6] = {
    -4.98843114573573548651E-8,
    9.50428062829859605134E-6,
    -6.45191435683965050962E-4,
    1.88843319396703850064E-2,
    -2.05525900955013891793E-1,
    9.99999999999999998822E-1,
};

static double cd[7] = {
    3.99982968972495980367E-12,
    9.15439215774657478799E-10,
    1.25001862479598821474E-7,
    1.22262789024179030997E-5,
    8.68029542941784300606E-4,
    4.12142090722199792936E-2,
    1.00000000000000000118E0,
};

/* Auxiliary function f(x) */
static double fn[10] = {
    4.21543555043677546506E-1,
    1.43407919780758885261E-1,
    1.15220955073585758835E-2,
    3.45017939782574027900E-4,
    4.63613749287867322088E-6,
    3.05568983790257605827E-8,
    1.02304514164907233465E-10,
    1.72010743268161828879E-13,
    1.34283276233062758925E-16,
    3.76329711269987889006E-20,
};

static double fd[10] = {
    /*  1.00000000000000000000E0, */
    7.51586398353378947175E-1,
    1.16888925859191382142E-1,
    6.44051526508858611005E-3,
    1.55934409164153020873E-4,
    1.84627567348930545870E-6,
    1.12699224763999035261E-8,
    3.60140029589371370404E-11,
    5.88754533621578410010E-14,
    4.52001434074129701496E-17,
    1.25443237090011264384E-20,
};

/* Auxiliary function g(x) */
static double gn[11] = {
    5.04442073643383265887E-1,
    1.97102833525523411709E-1,
    1.87648584092575249293E-2,
    6.84079380915393090172E-4,
    1.15138826111884280931E-5,
    9.82852443688422223854E-8,
    4.45344415861750144738E-10,
    1.08268041139020870318E-12,
    1.37555460633261799868E-15,
    8.36354435630677421531E-19,
    1.86958710162783235106E-22,
};

static double gd[11] = {
    /*  1.00000000000000000000E0, */
    1.47495759925128324529E0,
    3.37748989120019970451E-1,
    2.53603741420338795122E-2,
    8.14679107184306179049E-4,
    1.27545075667729118702E-5,
    1.04314589657571990585E-7,
    4.60680728146520428211E-10,
    1.10273215066240270757E-12,
    1.38796531259578871258E-15,
    8.39158816283118707363E-19,
    1.86958710162783236342E-22,
};

extern double MACHEP;

int fresnl(xxa, ssa, cca)
double xxa, *ssa, *cca;
{
    double f, g, cc, ss, c, s, t, u;
    double x, x2;

    if (cephes_isinf(xxa)) {
    cc = 0.5;
    ss = 0.5;
    goto done;
    }

    x = fabs(xxa);
    x2 = x * x;
    if (x2 < 2.5625) {
	t = x2 * x2;
	ss = x * x2 * polevl(t, sn, 5) / p1evl(t, sd, 6);
	cc = x * polevl(t, cn, 5) / polevl(t, cd, 6);
	goto done;
    }

    if (x > 36974.0) {
        /*
         * http://functions.wolfram.com/GammaBetaErf/FresnelC/06/02/
         * http://functions.wolfram.com/GammaBetaErf/FresnelS/06/02/
         */
	cc = 0.5 + 1/(M_PI*x) * sin(M_PI*x*x/2);
	ss = 0.5 - 1/(M_PI*x) * cos(M_PI*x*x/2);
	goto done;
    }


    /*             Asymptotic power series auxiliary functions
     *             for large argument
     */
    x2 = x * x;
    t = M_PI * x2;
    u = 1.0 / (t * t);
    t = 1.0 / t;
    f = 1.0 - u * polevl(u, fn, 9) / p1evl(u, fd, 10);
    g = t * polevl(u, gn, 10) / p1evl(u, gd, 11);

    t = M_PI_2 * x2;
    c = cos(t);
    s = sin(t);
    t = M_PI * x;
    cc = 0.5 + (f * s - g * c) / t;
    ss = 0.5 - (f * c + g * s) / t;

  done:
    if (xxa < 0.0) {
	cc = -cc;
	ss = -ss;
    }

    *cca = cc;
    *ssa = ss;
    return (0);
}
