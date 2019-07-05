#ifndef _ELLINT_RJ_GENERIC_GUARD
#define _ELLINT_RJ_GENERIC_GUARD

#include <stdlib.h>
#include <string.h>
#include "ellint_common.h"


#define ELLINT_RJ_CALC(VOID)	\
do { 	\
    prm = SQRT(pm);	\
/*    lam = ADD(ADD(MULcc(xrm, yrm), MULcc(zrm, xrm)), MULcc(yrm, zrm)); */\
    cct1[0] = cct2[1] = SQRT(xm);	\
    cct1[1] = cct2[2] = SQRT(ym);	\
    cct1[2] = cct2[0] = SQRT(zm);	\
    lam = ELLINT_POLY_FCN(dot2)(cct1, cct2, 3);	\
    dm = MULcc(MULcc(ADD(prm, cct1[0]), ADD(prm, cct1[1])), ADD(prm, cct1[2]));\
} while ( 0 )

#define ELLINT_RJ_UPDT(VOID)	\
do {	\
    Am = MULcr(ADD(Am, lam), 0.25);	\
    xm = MULcr(ADD(xm, lam), 0.25);	\
    ym = MULcr(ADD(ym, lam), 0.25);	\
    zm = MULcr(ADD(zm, lam), 0.25);	\
    pm = MULcr(ADD(pm, lam), 0.25);	\
    xxm = MULcr(xxm, 0.25);	\
    yym = MULcr(yym, 0.25);	\
    zzm = MULcr(zzm, 0.25);	\
    d4m *= 0.25;	\
    fterm *= 0.25;	\
} while ( 0 )


#define ASYMP_ZERO(r)		( ( 0.0 < (r) ) && ( (r) <= 5e-14 ) )
#define ABS_CLOSE_ZERO(r)	( ( 0.0 < (r) ) && ( (r) <= 1e-10 ) )


typedef enum rj_asymp_flag
{
    asymp_hugep,	/* x, y, z << p */
    asymp_tinyp,	/* p << geom. mean of x, y, z */
    asymp_tinyy,	/* max(x, y) == y << min(z, p) */
    asymp_hugey,	/* max(x, p) << min(y, z) == y */
    asymp_tinyx,	/* x << min(y, z, p) == min(y, p) */
    asymp_hugez,	/* max(x, y, p) == max(y, p) << z */
    asymp_nothing
} rj_asymp_t;


typedef struct rj_arg_cases
{
    bool retry_caupv;	/* should use Cauchy principal value */
    bool hit_pole;	/* singular */
    bool good_infinity;	/* the "good" kind of directed infinity */
    bool maybe_asymp;	/* might be good for an asympt. case */
} rj_arg_cases_t;


/* Comparison function based on the ordering of real part */
static int ELLINT_POLY_FCN(rcmp)(const void * a, const void * b)
{
    if ( (CREAL( *((const EllInt_Num_t *)a) ) <
          CREAL( *((const EllInt_Num_t *)b) )) )
    {
	return -1;
    }
    return 1;
}


/* Check whether the input arguments for RJ are out of domain, while set
 * corresponding flags in the class(ification) variable.
 *
 * NOTE: x, y, z must be in-order by real parts.
 */
static bool RJ_good_args(EllInt_Num_t x, EllInt_Num_t y, EllInt_Num_t z,
                         EllInt_Num_t p,
			 rj_arg_cases_t * restrict classify)
{
    double xr, xi, yr, yi, zr, zi, pr, pi;
    bool xyzreal_nonneg_atmost1z;

    xr = CREAL(x);
    xi = CIMAG(x);

    yr = CREAL(y);
    yi = CIMAG(y);

    zr = CREAL(z);
    zi = CIMAG(z);

    pr = CREAL(p);
    pi = CIMAG(p);

    if ( (classify->hit_pole = ( too_small(FABS(x)) &&
                                 too_small(FABS(y)) &&
				 (ph_is_not_pm_pi(z)) &&
				 !too_small(FABS(p)) )) )
    {
	return false;
    }
    if ( (classify->good_infinity = ( (Z_INFTY(x) || Z_INFTY(y) ||
                                       Z_INFTY(z) || Z_INFTY(p)) &&
                                      (ph_is_not_pm_pi(x) &&
				       ph_is_not_pm_pi(y) &&
				       ph_is_not_pm_pi(z)) )) )
    {
	return false;
    }
    /* "If x, y, z are real and nonnegative, at most one of them is 0, and the
     * fourth variable of RJ is negative, the Cauchy principal value ..." */
    xyzreal_nonneg_atmost1z = ( too_small(xi) && too_small(yi) &&
				too_small(zi) && (xr >= 0.0) && (yr > 0.0) );
    if ( too_small(pi) )
    {
	if ( (classify->retry_caupv = ( xyzreal_nonneg_atmost1z &&
	                                (pr < 0.0) )) )
	{
	    return false;
	}
	/* "Assume x, y, and z are real and nonnegative, at most one of them is
	 * 0, and p > 0" */
	if ( (classify->maybe_asymp = ( xyzreal_nonneg_atmost1z &&
	                                (pr > 0.0) )) )
	{
	    return true;
	}
    }

    /* "Let x, y, z have nonnegative real part and at most one of them [1] be
     * 0, while Re p > 0."
     *     [1] By "them", Carlson seems to have meant the numbers x, y, z
     *         themselves, rather than their "real parts". */
    if ( ( pr > 0.0 ) && ( xr >= 0.0 ) && C_ATMOST_1Z(x, y, z) )
    {
	return true;
    }

    /* "Alternatively, if p != 0 and |ph p| < pi ..." */
    if ( !( too_small(pr) && too_small(pi) ) && ph_is_not_pm_pi(p) )
    {
	/* "... either let x, y, z be real and nonnegative and at most one of
	 * them be 0, ..." */
	unsigned char flag = 0u;
	if (  xyzreal_nonneg_atmost1z )
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


/* Cauchy principal value dispatcher */
static EllInt_Status_t ellint_RJ_cpv_dispatch(EllInt_Num_t x, EllInt_Num_t y,
                                              EllInt_Num_t z, EllInt_Num_t p,
                                              double rerr,
					      double * restrict res)
{
    /* Retry with principal value evaluation, valid for reals. */
    EllInt_Status_t status, status_tmp;
    double r, pn, pq, xy, xypq, tmpres;
    double xct1[4];
    double xct2[4];

    r = rerr / 3.0;
    /*
    xx = CREAL(x);
    yy = CREAL(y);
    zz = CREAL(z);
    pp = CREAL(p);
    */
    xct1[0] = CREAL(x);
    xct1[1] = CREAL(y);
    xct1[2] = -CREAL(p);
    xct1[3] = CREAL(z);
    /* q = -pp; */
    /* xy = xx * yy; */
    xy = xct1[0] * xct1[1];
    xct2[3] = xct1[2] / xct1[3] + 1.0;
    pn = (fsum2(xct1, 3) - xy / xct1[3]) / xct2[3];

    status = ELLINT_STATUS_SUCCESS;

    status_tmp = fellint_RJ(xct1[0], xct1[1], xct1[3], pn, r, xct2);
    if ( HORRIBLE_STATUS(status_tmp) )
    {
	return status_tmp;
    } else if ( TROUBLESOME_STATUS(status_tmp) ) {
	status = status_tmp;
    }

    status_tmp = fellint_RF(xct1[0], xct1[1], xct1[3], r, xct2 + 1);
    if ( HORRIBLE_STATUS(status_tmp) )
    {
	return status_tmp;
    } else if ( TROUBLESOME_STATUS(status_tmp) ) {
	status = status_tmp;
    }

    pq = pn * xct1[2];
    xypq = xy + pq;
    status_tmp = fellint_RC(xypq, pq, r, xct2 + 2);
    if ( HORRIBLE_STATUS(status_tmp) )
    {
	return status_tmp;
    } else if ( TROUBLESOME_STATUS(status_tmp) ) {
	status = status_tmp;
    }
    xct1[0] = pn / xct1[3] - 1.0;
    xct1[1] = -3.0 / xct1[3];
    xct1[2] = 3.0 * sqrt(xy / (xypq * xct1[3]));

    /* tmpres = (pn - zz) * rjv - 3.0 * (rfv - sqrt(xy * zz / xypq) * rcv); */
    tmpres = fdot2(xct1, xct2, 3);
    /* tmpres /= q + zz */
    tmpres /= xct2[3];
    *res = tmpres;

    return status;
}


static inline rj_asymp_t RJ_asymp_conf( const double * restrict x,
	const double * restrict y, const double * restrict z,
	const double * restrict p,
	double * restrict a, double * restrict b, double * restrict c,
	double * restrict f, double * restrict g, double * restrict h)
{
    double t;

    /*
    t = (*z) / (*p);
    */
    /* this bound is neither sharp enough nor useful */
    /*
    if ( ASYMP_ZERO(t) )
    {
	*c = ((*x) + (*y) + (*z)) / 3.0;
	return asymp_hugep;
    }
    */

    /* this bound is sharp. RJ in this case behaves with logarithmic
     * singularity as p -> +0 */
    if ( ABS_CLOSE_ZERO(*p) ||
	 ( (fpclassify(*x) != FP_ZERO) && ASYMP_ZERO((*p) / (*x)) ) )
    {
	*f = sqrt((*x) * (*y) * (*z));
	return asymp_tinyp;
    }

    /* XXX until RJ(0, y, z, p) is implemented, this is useless */
    /*
    t = (*x) / fmin((*y), (*p));
    if ( ( 0.0 < (*x) && (*x) <= 1e-26 ) || ASYMP_ZERO(t) )
    {
	*h = sqrt((*y) * (*z));
	if ( (*h) / (*p) + 0.5 * ((*y) + (*z)) / (*h) <= sqrt((*h) / (*x)) )
	{
	    return asymp_tinyx;
	}
    }
    */

    t = (*y) / fmin((*z), (*p));
    if ( ( 0 < (*y) && (*y) <= 1e-26 ) || ASYMP_ZERO(t) )
    {
	/* bound fairly sharp even if p is large */
	*a = 0.5 * ((*x) + (*y));
	*g = sqrt((*x) * (*y));
	if ( ( (*a) / (*z) + (*a) / (*p) ) * fabs( log((*p) / (*a)) ) <= 1.0 )
	{
	    return asymp_tinyy;
	}
	return asymp_tinyy;
    }

    if ( (fpclassify(*x) != FP_ZERO) && ASYMP_ZERO(fmax((*z), (*p)) / (*x)) )
    {
	/* bound might not be sharp if x + 2p much larger than (yz) ** 2,
	 * but this is unlikely to be true anyway. */
	return asymp_hugey;
    }

    if ( (fpclassify(*z) != FP_ZERO) && ASYMP_ZERO(fmax((*y), (*p)) / (*z)) )
    {
	*b = 0.5 * ((*x) + (*y));
	*h = sqrt((*x) * (*y));
	/* when bounds are sharp */
	if ( fabs(log((*z) / ((*b) + (*h)))) <= sqrt((*z)) )
	{
	    return asymp_hugez;
	}
    }

    return asymp_nothing;
}


/* Prevent division by zero due to underflow in atan(sqrt(z)) / sqrt(z) and
 * square root of negative number in the real context. */
static inline EllInt_Num_t
safe_atan_sqrt_div(EllInt_Num_t z)
{
    EllInt_Num_t s;

    if ( too_small(CREAL(z)) && too_small(CIMAG(z)) )
    {
	return MKCR(1.0);
    }
#ifdef ELLINT_POLY_COMPLEX
    s = csqrt(z);
    return DIVcc(catan(s), s);
#else
    if ( z < 0.0 )
    {
	s = sqrt(-z);
	return atanh(s) / s;
    }
    s = sqrt(z);
    return atan(s) / s;
#endif
}


EllInt_Status_t
ELLINT_POLY_FCN(ellint_RJ) (EllInt_Num_t x, EllInt_Num_t y, EllInt_Num_t z,
                            EllInt_Num_t p,
	                    double rerr, EllInt_Num_t * restrict res)
{
    EllInt_Status_t status;
    rj_arg_cases_t classify;
    rj_asymp_t cres;
    unsigned int m;
    EllInt_Num_t A0, delta;
    EllInt_Num_t Am, xm, ym, zm, pm;
    EllInt_Num_t prm, lam, dm;
    EllInt_Num_t sm, rm;
    EllInt_Num_t d4m, xxm, yym, zzm;
    EllInt_Num_t pp, pp2, xyz, e2, e3, e4, e5;
    EllInt_Num_t tmp, t;
    double fterm, aAm;
    EllInt_Num_t cct1[6];
    EllInt_Num_t cct2[6];

    status = ELLINT_STATUS_SUCCESS;
    if ( ELLINT_BAD_RERR(rerr, 1.0e-4) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_RERR);
    }

    /* Put the symmetric arguments in the order of real parts. */
    cct1[0] = x;
    cct1[1] = y;
    cct1[2] = z;
    qsort(cct1, 3, sizeof (EllInt_Num_t), ELLINT_POLY_FCN(rcmp));
    x = cct1[0];
    y = cct1[1];
    z = cct1[2];

    memset(&classify, 0, sizeof (classify));
    if ( !RJ_good_args(x, y, z, p, &classify) )
    {
	if ( classify.good_infinity )
	{
	    ELLINT_RETURN_WITH(ELLINT_STATUS_SUCCESS, CZERO);
	}

	if ( classify.retry_caupv )
	{
	    /* Retry with principal value evaluation, valid for reals. */
	    double tmpres;

	    status = ellint_RJ_cpv_dispatch(x, y, z, p, rerr, &tmpres);
	    if ( HORRIBLE_STATUS(status) )
	    {
		*res = ELLINT_NAN;
	    } else {
		*res = MKCR(tmpres);
	    }
	    return status;
	} else if ( classify.hit_pole ) {
	    ELLINT_RETURN_WITH(ELLINT_STATUS_SINGULAR, CPINF);
	} else {
	    ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_ARGS);
	}
    }

    if ( classify.maybe_asymp )
    {
	/* might be dealt with by asymptotic expansion of real-arg RJ */
	double xx, yy, zz, ppr, a, b, c, f, g, h;
	double tmpres;

	ppr = CREAL(p);
	xx = CREAL(x);
	yy = CREAL(y);
	zz = CREAL(z);
	cres = RJ_asymp_conf(&xx, &yy, &zz, &ppr, &a, &b, &c, &f, &g, &h);
	switch ( cres )
	{
	    case asymp_hugep :
		{
		    status = fellint_RF(xx, yy, zz, rerr, &tmpres);
		    tmpres = 3.0 * (tmpres - 0.5 * M_PI / sqrt(ppr)) / ppr;
		}
		break;
	    case asymp_tinyp :
		{
		    EllInt_Status_t status_tmp;
		    double lamt, alpha, beta;
		    double xct1[6];
		    double xct2[6];
		    rerr *= 0.5;
		    xct1[1] = xct2[0] = sqrt(xx);
		    xct1[2] = xct2[1] = sqrt(yy);
		    xct1[0] = xct2[2] = sqrt(zz);
		    lamt = fdot2(xct1, xct2, 3);
		    xct2[3] = xct2[4] = xct2[5] = ppr;
		    alpha = fdot2(xct1, xct2 + 3, 3) + f;
		    alpha = alpha * alpha;
		    beta = ppr + lamt;
		    beta = beta * beta * ppr;
		    status_tmp = fellint_RC(alpha, beta, rerr, xct2);
		    status = fellint_RJ(xx + lamt, yy + lamt, zz + lamt,
					ppr + lamt, rerr, xct2 + 1);
		    if ( status_tmp != ELLINT_STATUS_SUCCESS )
		    {
			status = status_tmp;
		    }
		    xct1[0] = 3.0;
		    xct1[1] = 2.0;
		    tmpres = fdot2(xct1, xct2, 2);
		}
		break;
	    case asymp_tinyy :
		{
		    double r_est_h, r_est_l, tx;
		    status = fellint_RC(1.0, ppr / zz, rerr, &tx);
		    tmpres = (log(8.0 * zz / (a + g)) - 2.0 * tx);
		    tx = log(2.0 * ppr / (a + g)) / (tmpres * ppr);
		    r_est_l = tx * g / (1.0 - g / ppr);
		    r_est_h = tx * a * (1.0 + 0.5 * ppr / zz) / (1.0 - a / ppr);
		    /* if asymptotic expansion found to violate error bound
		     * after the fact */
		    if ( r_est_h - r_est_l >= 2.0 * rerr )
		    {
			cres = asymp_nothing;
			break;
		    } else {
			tmpres += r_est_l;
			tmpres *= 1.5 / (sqrt(zz) * ppr);
		    }
		}
		break;
	    case asymp_hugey :
		{
		    double rr, t1, t2;
		    EllInt_Status_t status_tmp;
		    rr = rerr / 3.0;
		    tmpres = 1.0 / sqrt(yy * zz);
		    status_tmp = fellint_RC(xx, ppr, rr, &t1);
		    status = fellint_RG(0.0, yy, zz, rr, &t2);
		    if ( status_tmp != ELLINT_STATUS_SUCCESS )
		    {
			status = status_tmp;
		    }
		    tmpres *= (3.0 * t1 - 2.0 * t2 * tmpres);
		}
		break;
	    case asymp_tinyx :
		{
		    status = fellint_RJ(0.0, yy, zz, ppr, rerr, &tmpres);
		    tmpres -= 3.0 * sqrt(xx) / (h * ppr);
		}
		break;
	    case asymp_hugez :
		{
		    double tt;
		    double r_est;
#if 0
		    if ( ASYMP_ZERO(xx / fmin(yy, ppr)) )
		    {
			/* XXX */
			status = fellint_RC(ppr, yy, rerr, &tmpres);
			tmpres *= 3.0 / sqrt(ppr * zz);
		    } else {
#endif
			tt = h + ppr;
			tt *= tt;
			status = fellint_RC(tt, 2.0 * (b + h) * ppr,
			                    rerr, &tmpres);
			r_est = 0.25 * (0.5 +
				log1p(2.0 * zz / sqrt(h * ppr))) / (tmpres *
				zz);
			/* if asymptotic expansion found to violate error bound
			 * after the fact */
			if ( r_est >= rerr )
			{
			    cres = asymp_nothing;
			    break;
			} else {
			    tmpres *= 3.0 / sqrt(zz);
			}
#if 0
		    }
#endif
		}
		break;
	    case asymp_nothing:
		break;
	}	/* end of switch (argument classification for asympt.) */
	if ( cres != asymp_nothing )
	{
	    *res = tmpres;
	    return status;
	}
    } 

    /* A0 = DIVcr(ADD(ADD(ADD(x, y), z), MULcr(p, 2.0)), 5.0); */
    cct1[3] = p;
    cct1[4] = p;
    A0 = DIVcr(ELLINT_POLY_FCN(sum2)(cct1, 5), 5.0);
    delta = MULcc(MULcc(SUB(p, x), SUB(p, y)), SUB(p, z));
    rerr = ELLINT_OCRT(rerr / 5.0);
    xxm = SUB(A0, x);
    yym = SUB(A0, y);
    zzm = SUB(A0, z);
    fterm = ELLINT_FMAX4(FABS(xxm), FABS(yym), FABS(zzm),
			 FABS(SUB(A0, p))) / rerr;

    /* m = 0; */
    Am = A0;
    d4m = 1.0;
    xm = x;
    ym = y;
    zm = z;
    pm = p;

    ELLINT_RJ_CALC();
    sm = MULcr(dm, 0.5);

    /* next */
    ELLINT_RJ_UPDT();

    m = 1;
    while ( (aAm = FABS(Am)) <= fterm ||
	    aAm <= ELLINT_FMAX4(FABS(Am - xm), FABS(Am - ym),
                                FABS(Am - zm), FABS(Am - pm)) ) 
    {
	if ( m > ELLINT_MAXITER )
	{
	    status = ELLINT_STATUS_NITER;
	    break;
	}

	rm = MULcc(sm, ADDcr(SQRT(ADDcr(DIVcc(MULcr(delta, d4m),
				              MULcc(sm, sm)), 1.0)),
	                        1.0));
	ELLINT_RJ_CALC();
	sm = DIVcc(MULcr(SUB(MULcc(rm, dm), MULcr(delta, d4m * d4m)), 0.5),
	           ADD(dm, MULcr(rm, d4m)));

	/* next */
	ELLINT_RJ_UPDT();
	m += 1;
    }

    /* Burn some extra cycles re-balancing Am as the "true" centroid */
    cct1[0] = xm;
    cct1[1] = ym;
    cct1[2] = zm;
    cct1[3] = pm;
    cct1[4] = pm;
    Am = DIVcr(ELLINT_POLY_FCN(sum2)(cct1, 5), 5.0);
    /* {xx,yy,zz}m /= Am */
    xxm = DIVcc(xxm, Am);
    yym = DIVcc(yym, Am);
    zzm = DIVcc(zzm, Am);
    /* pp = -0.5 * (xxm + yym + zzm) */
    cct1[0] = cct2[2] = xxm;
    cct1[1] = cct2[0] = yym;
    cct1[2] = cct2[1] = zzm;
    pp = MULcr(ELLINT_POLY_FCN(sum2)(cct1, 3), -0.5);
    cct1[3] = cct2[3] = pp;
    cct1[3] = MULcr(cct1[3], -3.0);
    /* pp2 = pp * pp */
    pp2 = MULcc(pp, pp);
    /* xyz = xxm * yym * zzm */
    xyz = MULcc(MULcc(xxm, yym), zzm);
    /* e2 = xxm * yym + zzm * xxm + yym * zzm - pp2 * 3.0 */
    e2 = ELLINT_POLY_FCN(dot2)(cct1, cct2, 4);
    /* e3 = xyz + 2.0 * pp * (e2 + 2.0 * pp2) */
    e3 = ADD(xyz,
             MULcc(MULcr(pp, 2.0),
                   ADD(e2, MULcr(pp2, 2.0))));
    /* e4 = (2.0 * xyz + (e2 + 3.0 * pp2) * pp) * pp */
    e4 = MULcc(ADD(MULcr(xyz, 2.0),
		   MULcc(ADD(e2,
                             MULcr(pp2, 3.0)),
                         pp)),
               pp);
    /* e5 = xyz * pp2 */
    e5 = MULcc(xyz, pp2);
    /* tmp = d4m * pow(sqrt(Am), -3) */
    t = SQRT(Am);
    tmp = MULcc(MULcc(t, t), t);
    tmp = DIVrc(d4m, tmp);
    cct1[0] = HORNER(e2, ELLINT_RDJ_C1, 3);
    cct1[1] = HORNER(e3, ELLINT_RDJ_C2, 2);
    cct1[2] = HORNER(e2, ELLINT_RDJ_C3, 2);
    cct1[3] = HORNER(e2, ELLINT_RDJ_C4, 1);
    cct1[4] = HORNER(e2, ELLINT_RDJ_C5, 1);
    cct1[5] = MULcr(e3, ELLINT_RDJ_C5[1]);

    cct2[0] = MKCR(1.0);
    cct2[1] = MKCR(1.0);
    cct2[2] = e3;
    cct2[3] = e4;
    cct2[4] = e5;
    cct2[5] = e4;
    t = ADDcr(DIVcr(ELLINT_POLY_FCN(dot2)(cct1, cct2, 6), ELLINT_RDJ_DENOM),
              1.0);
    tmp = MULcc(tmp, t);
    t = DIVcc(MULcr(delta, d4m), MULcc(sm, sm));
    tmp = ADD(tmp, DIVcc(MULcr(safe_atan_sqrt_div(t), 3.0), sm));

    *res = tmp;
    return status;
}


#endif  /* _ELLINT_RJ_GENERIC_GUARD */
