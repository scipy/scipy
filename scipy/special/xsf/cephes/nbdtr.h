/* Translated into C++ by SciPy developers in 2024.
 * Original header with Copyright information appears below.
 */

/*                                                     nbdtr.c
 *
 *     Negative binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * double p, y, nbdtr();
 *
 * y = nbdtr( k, n, p );
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms 0 through k of the negative
 * binomial distribution:
 *
 *   k
 *   --  ( n+j-1 )   n      j
 *   >   (       )  p  (1-p)
 *   --  (   j   )
 *  j=0
 *
 * In a sequence of Bernoulli trials, this is the probability
 * that k or fewer failures precede the nth success.
 *
 * The terms are not computed individually; instead the incomplete
 * beta integral is employed, according to the formula
 *
 * y = nbdtr( k, n, p ) = incbet( n, k+1, p ).
 *
 * The arguments must be positive, with p ranging from 0 to 1.
 *
 * ACCURACY:
 *
 * Tested at random points (a,b,p), with p between 0 and 1.
 *
 *               a,b                     Relative error:
 * arithmetic  domain     # trials      peak         rms
 *    IEEE     0,100       100000      1.7e-13     8.8e-15
 * See also incbet.c.
 *
 */
/*							nbdtrc.c
 *
 *	Complemented negative binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * double p, y, nbdtrc();
 *
 * y = nbdtrc( k, n, p );
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms k+1 to infinity of the negative
 * binomial distribution:
 *
 *   inf
 *   --  ( n+j-1 )   n      j
 *   >   (       )  p  (1-p)
 *   --  (   j   )
 *  j=k+1
 *
 * The terms are not computed individually; instead the incomplete
 * beta integral is employed, according to the formula
 *
 * y = nbdtrc( k, n, p ) = incbet( k+1, n, 1-p ).
 *
 * The arguments must be positive, with p ranging from 0 to 1.
 *
 * ACCURACY:
 *
 * Tested at random points (a,b,p), with p between 0 and 1.
 *
 *               a,b                     Relative error:
 * arithmetic  domain     # trials      peak         rms
 *    IEEE     0,100       100000      1.7e-13     8.8e-15
 * See also incbet.c.
 */

/*                                                     nbdtrc
 *
 *     Complemented negative binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * double p, y, nbdtrc();
 *
 * y = nbdtrc( k, n, p );
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms k+1 to infinity of the negative
 * binomial distribution:
 *
 *   inf
 *   --  ( n+j-1 )   n      j
 *   >   (       )  p  (1-p)
 *   --  (   j   )
 *  j=k+1
 *
 * The terms are not computed individually; instead the incomplete
 * beta integral is employed, according to the formula
 *
 * y = nbdtrc( k, n, p ) = incbet( k+1, n, 1-p ).
 *
 * The arguments must be positive, with p ranging from 0 to 1.
 *
 * ACCURACY:
 *
 * See incbet.c.
 */
/*							nbdtri
 *
 *	Functional inverse of negative binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * double p, y, nbdtri();
 *
 * p = nbdtri( k, n, y );
 *
 * DESCRIPTION:
 *
 * Finds the argument p such that nbdtr(k,n,p) is equal to y.
 *
 * ACCURACY:
 *
 * Tested at random points (a,b,y), with y between 0 and 1.
 *
 *               a,b                     Relative error:
 * arithmetic  domain     # trials      peak         rms
 *    IEEE     0,100       100000      1.5e-14     8.5e-16
 * See also incbi.c.
 */

/*
 * Cephes Math Library Release 2.3:  March, 1995
 * Copyright 1984, 1987, 1995 by Stephen L. Moshier
 */
#pragma once

#include "../config.h"
#include "../error.h"

#include "incbet.h"
#include "incbi.h"

namespace xsf {
namespace cephes {

    XSF_HOST_DEVICE inline double nbdtrc(int k, int n, double p) {
        double dk, dn;

        if ((p < 0.0) || (p > 1.0)) {
            goto domerr;
        }
        if (k < 0) {
        domerr:
            set_error("nbdtr", SF_ERROR_DOMAIN, NULL);
            return (std::numeric_limits<double>::quiet_NaN());
        }

        dk = k + 1;
        dn = n;
        return (incbet(dk, dn, 1.0 - p));
    }

    XSF_HOST_DEVICE inline double nbdtr(int k, int n, double p) {
        double dk, dn;

        if ((p < 0.0) || (p > 1.0)) {
            goto domerr;
        }
        if (k < 0) {
        domerr:
            set_error("nbdtr", SF_ERROR_DOMAIN, NULL);
            return (std::numeric_limits<double>::quiet_NaN());
        }
        dk = k + 1;
        dn = n;
        return (incbet(dn, dk, p));
    }

    XSF_HOST_DEVICE inline double nbdtri(int k, int n, double p) {
        double dk, dn, w;

        if ((p < 0.0) || (p > 1.0)) {
            goto domerr;
        }
        if (k < 0) {
        domerr:
            set_error("nbdtri", SF_ERROR_DOMAIN, NULL);
            return (std::numeric_limits<double>::quiet_NaN());
        }
        dk = k + 1;
        dn = n;
        w = incbi(dn, dk, p);
        return (w);
    }

} // namespace cephes
} // namespace xsf
