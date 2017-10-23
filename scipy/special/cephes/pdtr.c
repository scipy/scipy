/*                                                     pdtr.c
 *
 *     Poisson distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k;
 * double m, y, pdtr();
 *
 * y = pdtr( k, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the sum of the first k terms of the Poisson
 * distribution:
 *
 *   k         j
 *   --   -m  m
 *   >   e    --
 *   --       j!
 *  j=0
 *
 * The terms are not summed directly; instead the incomplete
 * Gamma integral is employed, according to the relation
 *
 * y = pdtr( k, m ) = igamc( k+1, m ).
 *
 * The arguments must both be nonnegative.
 *
 *
 *
 * ACCURACY:
 *
 * See igamc().
 *
 */
/*  pdtrc()
 *
 *  Complemented poisson distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k;
 * double m, y, pdtrc();
 *
 * y = pdtrc( k, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms k+1 to infinity of the Poisson
 * distribution:
 *
 *  inf.       j
 *   --   -m  m
 *   >   e    --
 *   --       j!
 *  j=k+1
 *
 * The terms are not summed directly; instead the incomplete
 * Gamma integral is employed, according to the formula
 *
 * y = pdtrc( k, m ) = igam( k+1, m ).
 *
 * The arguments must both be nonnegative.
 *
 *
 *
 * ACCURACY:
 *
 * See igam.c.
 *
 */
/*  pdtri()
 *
 *  Inverse Poisson distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k;
 * double m, y, pdtr();
 *
 * m = pdtri( k, y );
 *
 *
 *
 *
 * DESCRIPTION:
 *
 * Finds the Poisson variable x such that the integral
 * from 0 to x of the Poisson density is equal to the
 * given probability y.
 *
 * This is accomplished using the inverse Gamma integral
 * function and the relation
 *
 *    m = igamci( k+1, y ).
 *
 *
 *
 *
 * ACCURACY:
 *
 * See igami.c.
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * pdtri domain    y < 0 or y >= 1       0.0
 *                     k < 0
 *
 */

/*
 * Cephes Math Library Release 2.3:  March, 1995
 * Copyright 1984, 1987, 1995 by Stephen L. Moshier
 */

#include "mconf.h"

double pdtrc(k, m)
int k;
double m;
{
    double v;

    if ((k < 0) || (m < 0.0)) {
        mtherr("pdtrc", DOMAIN);
        return (NPY_NAN);
    }
    if (m == 0.0) {
        return 0.0;
    }
    v = k + 1;
    return (igam(v, m));
}


double pdtr(k, m)
int k;
double m;
{
    double v;

    if ((k < 0) || (m < 0.0)) {
        mtherr("pdtr", DOMAIN);
        return (NPY_NAN);
    }
    if (m == 0.0) {
        return 1.0;
    }
    v = k + 1;
    return (igamc(v, m));
}


double pdtri(k, y)
int k;
double y;
{
    double v;

    if ((k < 0) || (y < 0.0) || (y >= 1.0)) {
        mtherr("pdtri", DOMAIN);
        return (NPY_NAN);
    }
    v = k + 1;
    v = igamci(v, y);
    return (v);
}
