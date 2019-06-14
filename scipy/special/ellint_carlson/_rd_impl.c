#ifndef _ELLINT_RD_GENERIC_GUARD
#define _ELLINT_RD_GENERIC_GUARD

#include "ellint_common.h"
#include "ellint_poly.h"
#include "ellint_argcheck.h"


int
ELLINT_POLY_FCN(ellint_RD) (EllInt_Num_t x, EllInt_Num_t y, EllInt_Num_t z,
                            double rerr, EllInt_Num_t * restrict res)
{
    int status;

    if ( ELLINT_BAD_RERR(rerr, 1.0e-4) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_RERR);
    }

    if ( ( ( too_small(fabs(x)) ) && ( too_small(fabs(y)) ) ) ||
         ( too_small(fabs(z)) ) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_ARGS);
    }

    unsigned int m;
    EllInt_Num_t A0, A0mx, A0my, A0mz;
    EllInt_Num_t xm, ym, zm, Am;
    EllInt_Num_t adt, adtc, tmp;
    EllInt_Num_t e2, e3, e4, e5, xxm, yym, zzm, xy, zz2;
    double fterm, Q, d4m;

    A0 = (x + y + 3.0 * z) / 5.0;
    A0mx = A0 - x;
    A0my = A0 - y;
    A0mz = A0 - z;
    xxm = A0mx;
    yym = A0my;
    Q = ELLINT_FMAX3(fabs(A0mx), fabs(A0my), fabs(A0mz)) / ELLINT_SXRT(0.25 *
	                                                               rerr);
    Am = A0;
    d4m = 1.0;
    adt = 0.0;
    adtc = 0.0;
    xm = x;
    ym = y;
    zm = z;
    fterm = Q;
    m = 0;

    while ( 1 )
    {
	EllInt_Num_t xrm, yrm, zrm, lam, adttmp;

        if ( fterm < fabs(Am) )
	{
            break;
	}

	if ( m > ELLINT_MAXITER )
	{
	    ELLINT_FAIL_WITH(ELLINT_STATUS_NITER);
	}

        xrm = sqrt(xm);
        yrm = sqrt(ym);
        zrm = sqrt(zm);
        lam = xrm * (yrm + zrm) + yrm * zrm;

        /* Neumaier summation -- begin */
        /* https://doi.org/10.1002%2Fzamm.19740540106 */
        tmp = d4m / zrm / (zm + lam);
        adttmp = adt + tmp;
#ifdef ELLINT_POLY_REAL
        if ( fabs(adt) >= fabs(tmp) )
	{
            adtc = (adt - adttmp) + tmp;
	} else {
            adtc += (tmp - adttmp) + adt;
	}
#else	/* complex case, split into real and imaginary parts */
	{
	    double adttmp_r, adttmp_i, adt_r, adt_i, adtc_r, adtc_i;

	    adttmp_r = creal(adttmp);
	    adttmp_i = cimag(adttmp);
	    adt_r = creal(adt);
	    adt_i = cimag(adt);
	    adtc_r = creal(adtc);
	    adtc_i = cimag(adtc);

	    if ( fabs(adt_r) >= fabs(creal(tmp)) )
	    {
		adtc_r = (adt_r - adttmp_r) + creal(tmp);
	    } else {
		adtc_r += (creal(tmp) - adttmp_r) + adt_r;
	    }

	    if ( fabs(adt_i) >= fabs(cimag(tmp)) )
	    {
		adtc_i = (adt_i - adttmp_i) + cimag(tmp);
	    } else {
		adtc_i += (cimag(tmp) - adttmp_i) + adt_i;
	    }
	    adttmp = adttmp_r + adttmp_i * I;
	    adtc = adtc_r + adtc_i * I;
	}
#endif
        adt = adttmp;
        /* Neumaier summation -- end */

        Am = (Am + lam) * 0.25;
        xm = (xm + lam) * 0.25;
        ym = (ym + lam) * 0.25;
        zm = (zm + lam) * 0.25;
        fterm *= 0.25;
        d4m *= 0.25;
        xxm *= 0.25;
        yym *= 0.25;

        m += 1;
    }
    /* Neumaier summation -- post-processing */
    adt += adtc;

    xxm /= Am;
    yym /= Am;
    zzm = -(xxm + yym) / 3.0;
    xy = xxm * yym;
    zz2 = zzm * zzm;
    e2 = xy - 6.0 * zz2;
    e3 = (3.0 * xy - 8.0 * zz2) * zzm;
    e4 = 3.0 * (xy - zz2) * zz2;
    e5 = xy * zz2 * zzm;
    tmp = d4m * pow(sqrt(Am), -3);
    tmp *= (1.0 + e2 * (ELLINT_RDJ_c02 + e2 * ELLINT_RDJ_c22 +
                        e3 * ELLINT_RDJ_c23) + e3 * ELLINT_RDJ_c03 +
	    e4 * ELLINT_RDJ_c04 + e5 * ELLINT_RDJ_c05);
    tmp += 3.0 * adt;

    *res = tmp;
    status = ELLINT_STATUS_SUCCESS;
    return status;
}


#endif  /* _ELLINT_RD_GENERIC_GUARD */
