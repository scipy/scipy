#ifndef _ELLINT_RJ_GENERIC_GUARD
#define _ELLINT_RJ_GENERIC_GUARD

#include <stdlib.h>
#include "ellint_common.h"
#include "ellint_poly.h"
#include "ellint_argcheck.h"

#define ELLINT_RJ_1over3	(0.3333333333333333333333333)
#define ELLINT_RJ_1over5	(0.2000000000000000000000000)


#define ELLINT_RJ_CALC(VOID)	\
do { 				\
    xrm = sqrt(xm);		\
    yrm = sqrt(ym);		\
    zrm = sqrt(zm);		\
    prm = sqrt(pm);		\
    lam = xrm * (yrm + zrm) + yrm * zrm;		\
    dm = (prm + xrm) * (prm + yrm) * (prm + zrm);	\
} while ( 0 )

#define ELLINT_RJ_UPDT(VOID)	\
do {				\
    Am = (Am + lam) * 0.25;	\
    xm = (xm + lam) * 0.25;	\
    ym = (ym + lam) * 0.25;	\
    zm = (zm + lam) * 0.25;	\
    pm = (pm + lam) * 0.25;	\
    d4m *= 0.25;		\
    xxm *= 0.25;		\
    yym *= 0.25;		\
    zzm *= 0.25;		\
    fterm *= 0.25;		\
} while ( 0 )


/* Comparison function based on the ordering of real part */
static int fcmp(const void * a, const void * b)
{
    double fdiff;
    fdiff = (creal( *((const EllInt_Num_t *)a) ) -
             creal( *((const EllInt_Num_t *)b) ));
    if ( fdiff < 0.0 )
    {
	return -1;
    } else if ( fdiff == 0.0 ) {
	return 0;
    }
    return 1;
}


/* Prevent division by zero due to underflow in atan(sqrt(z)) / sqrt(z) and
 * square root of negative number in the real context. */
static inline EllInt_Num_t safe_atan_sqrt_div(EllInt_Num_t z)
{
    if ( too_small(creal(z)) && too_small(cimag(z)) )
    {
	return (EllInt_Num_t)1.0;
    }
#ifdef ELLINT_POLY_COMPLEX
    EllInt_Num_t s;
    s = sqrt(z);
    return atan(s) / s;
#else
    double s;
    if ( z < 0.0 )
    {
	s = sqrt(-z);
	return atanh(s) / s;
    }
    s = sqrt(z);
    return atan(s) / s;
#endif
}


/* Check whether the input arguments for RJ are out of domain.
 * If "retry" is set, p is real negative and x, y, z are real and non-negative,
 * and we may try to evaluate the Cauchy principal value.
 *
 * NOTE: x, y, z must be in-order by real parts.
 */
static inline bool RJ_good_args(EllInt_Num_t x, EllInt_Num_t y, EllInt_Num_t z,
                                EllInt_Num_t p, bool * restrict retry)
{
    double xr, xi, yr, yi, zr, zi, pr, pi;
    bool allreal_nonneg_atmost1z;

    xr = creal(x);
    xi = cimag(x);

    yr = creal(y);
    yi = cimag(y);

    zr = creal(z);
    zi = cimag(z);

    pr = creal(p);
    pi = cimag(p);

    allreal_nonneg_atmost1z = ( too_small(xi) && too_small(yi) &&
                                too_small(zi) && (xr >= 0.0) && (yr > 0.0) );
    /* "If x, y, z are real and nonnegative, at most one of them is 0, and the
     * fourth variable of RJ is negative, the Cauchy principal value ..." */
    if ( allreal_nonneg_atmost1z && too_small(pi) && (pr < 0.0) )
    {
	*retry = true;
	return false;
    }

    /* "Let x, y, z have nonnegative real part and at most one of them [1] be
     * 0, while Re p > 0."
     *     [1] By "them", Carlson seems to have meant the numbers x, y, z
     *         themselves, rather than their "real parts". */
    if ( ( pr > 0.0 ) && ( xr >= 0.0 ) && C_ATMOST_1Z(x, y, z) )
    {
	return true;
    }

    /* "Alternatively, if p != 0 ... " */
    if ( !( too_small(pr) && too_small(pi) ) )
    {
	/* "... and |ph p| < pi ..." */
	/* "... either let x, y, z be real and nonnegative and at most one of
	 * them be 0, ..." */
	unsigned char flag = 0u;
	if (  allreal_nonneg_atmost1z )
	{
	    flag = 1u;
	/* "... or else let two of the variables x, y, z be nonzero and
	 * conjugate complex with phase less in magnitude than pi and the
	 * third variable be real and nonnegative." */
	}
	if ( ELLINT_1R_2CONJ(xr, xi, y, z) ||
	     ELLINT_1R_2CONJ(yr, yi, z, x) ||
	     ELLINT_1R_2CONJ(zr, zi, x, y) )
	{
	    flag ^= 1u;
	}
	return (bool)flag;
    }
    return false;
}


int
ELLINT_POLY_FCN(ellint_RJ) (EllInt_Num_t x, EllInt_Num_t y, EllInt_Num_t z,
                            EllInt_Num_t p,
	                    double rerr, EllInt_Num_t * restrict res)
{
    int status;
    bool retry;

    if ( ELLINT_BAD_RERR(rerr, 1.0e-4) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_RERR);
    }

    /* Put the symmetric arguments in the order of real parts. */
    EllInt_Num_t na[3] = {x, y, z};
    qsort(na, 3, sizeof (EllInt_Num_t), fcmp);
    x = na[0];
    y = na[1];
    z = na[2];

    retry = false;
    if ( !RJ_good_args(x, y, z, p, &retry) )
    {
	if ( retry )
	{
	    /* Retry with principal value evaluation, valid for reals. */
	    double r, pp;

	    r = rerr / 3.0;

	    double xx, yy, zz, q, pn, rjv, rfv, rcv, pq, xy, xypq, tmpres;
	    pp = creal(p);
	    xx = creal(x);
	    yy = creal(y);
	    zz = creal(z);

	    q = -pp;
	    xy = xx * yy;
	    pn = (zz * (xx + yy + q) - xy) / (zz + q);

	    status = fellint_RJ(xx, yy, zz, pn, r, &rjv);
	    if ( status != ELLINT_STATUS_SUCCESS )
	    {
		ELLINT_FAIL_WITH(status);
	    }

	    status = fellint_RF(xx, yy, zz, r, &rfv);
	    if ( status != ELLINT_STATUS_SUCCESS )
	    {
		ELLINT_FAIL_WITH(status);
	    }

	    pq = pn * q;
	    xypq = xy + pq;
	    status = fellint_RC(xypq, pq, r, &rcv);
	    if ( status != ELLINT_STATUS_SUCCESS )
	    {
		ELLINT_FAIL_WITH(status);
	    }

	    tmpres = (pn - zz) * rjv - 3.0 * (rfv -
					      sqrt(xy * zz / xypq) * rcv);
	    tmpres /= (q + zz);

	    *res = (EllInt_Num_t)tmpres;
	    return ELLINT_STATUS_SUCCESS;
	} else {
	    ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_ARGS);
	}
    }

    unsigned int m;
    EllInt_Num_t A0, delta, A0mx, A0my, A0mz;
    EllInt_Num_t Am, xm, ym, zm, pm;
    EllInt_Num_t xrm, yrm, zrm, prm, lam, dm;
    EllInt_Num_t sm, rm;
    EllInt_Num_t d4m, xxm, yym, zzm;
    EllInt_Num_t sprev;
    EllInt_Num_t pp, pp2, xyz, e2, e3, e4, e5;
    EllInt_Num_t tmp, t;
    double Q, fterm;

    A0 = (x + y + z + 2.0 * p) / 5.0;
    delta = (p - x) * (p - y) * (p - z);
    A0mx = A0 - x;
    A0my = A0 - y;
    A0mz = A0 - z;
    Q = (ELLINT_FMAX4(fabs(A0mx), fabs(A0my), fabs(A0mz), fabs(A0 - p)) /
         ELLINT_SXRT(rerr * 0.25));

    /* m = 0; */
    fterm = Q;
    Am = A0;
    d4m = 1.0;
    xm = x;
    ym = y;
    zm = z;
    pm = p;
    xxm = A0mx;
    yym = A0my;
    zzm = A0mz;

    ELLINT_RJ_CALC();
    sprev = dm * 0.5;

    /* next */
    ELLINT_RJ_UPDT();
    m = 1;

    while ( 1 )
    {
        if ( fterm < fabs(Am) )
	{
            break;
	}

	if ( m > ELLINT_MAXITER )
	{
	    ELLINT_FAIL_WITH(ELLINT_STATUS_NITER);
	}

	rm = sprev * (1.0 + sqrt(1.0 + d4m * delta / (sprev * sprev)));
	ELLINT_RJ_CALC();
	sm = 0.5 * (rm * dm - d4m * d4m * delta) / (dm + rm * d4m);

	/* next */
	ELLINT_RJ_UPDT();
	sprev = sm;
	m += 1;
    }

    xxm /= Am;
    yym /= Am;
    zzm /= Am;
    pp = -0.5 * (xxm + yym + zzm);
    pp2 = pp * pp;
    xyz = xxm * yym * zzm;
    e2 = xxm * yym + zzm * xxm + yym * zzm - 3.0 * pp2;
    e3 = xyz + 2.0 * pp * (e2 + 2.0 * pp2);
    e4 = (2.0 * xyz + (e2 + 3.0 * pp2) * pp) * pp;
    e5 = xyz * pp2;
    tmp = d4m * pow(sqrt(Am), -3);
    tmp *= (1.0 + e2 * (ELLINT_RDJ_c02 + e2 * ELLINT_RDJ_c22 +
                        e3 * ELLINT_RDJ_c23) + e3 * ELLINT_RDJ_c03 +
	    e4 * ELLINT_RDJ_c04 + e5 * ELLINT_RDJ_c05);
    t = (d4m * delta) / (sprev * sprev);
    tmp += 3.0 * safe_atan_sqrt_div(t) / sprev;

    *res = tmp;
    status = ELLINT_STATUS_SUCCESS;
    return status;
}


#endif  /* _ELLINT_RJ_GENERIC_GUARD */
