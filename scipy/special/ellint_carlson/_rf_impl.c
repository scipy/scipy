#ifndef _ELLINT_RF_GENERIC_GUARD
#define _ELLINT_RF_GENERIC_GUARD

#include "ellint_common.h"
#include "ellint_poly.h"
#include "ellint_argcheck.h"

#define ELLINT_RF_c22	(+0.04166666666666666666666667)	/*  1 / 24 */
#define ELLINT_RF_c23	(-0.06818181818181818181818182)	/* -3 / 44 */
#define ELLINT_RF_c03	(+0.07142857142857142857142857)	/*  1 / 14 */
#define ELLINT_RF_c222	(-0.02403846153846153846153846)	/*  5 / 208 */
#define ELLINT_RF_c33	(+0.02884615384615384615384615)	/*  3 / 104 */
#define ELLINT_RF_c223	(+0.06250000000000000000000000)	/*  1 / 16 */

int
ELLINT_POLY_FCN(ellint_RF) (EllInt_Num_t x, EllInt_Num_t y, EllInt_Num_t z,
	                    double rerr, EllInt_Num_t * restrict res)
{
    int status;

    if ( ELLINT_BAD_RERR((rerr), (3.0e-4)) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_RERR);
    }

    if ( !(C_ATMOST_1Z(x, y, z)) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_ARGS);
    }

    unsigned int m;
    EllInt_Num_t A0, A0mx, A0my, A0mz;
    EllInt_Num_t xm, ym, zm, Am;
    EllInt_Num_t e2, e3, xxm, yym, zzm, s;
    double fterm, Q;

    A0 = (x + y + z) / 3.0;
    A0mx = A0 - x;
    A0my = A0 - y;
    A0mz = A0 - z;
    xxm = A0mx;
    yym = A0my;
    Q = ELLINT_FMAX3((fabs(A0mx)), \
                     (fabs(A0my)), (fabs(A0mz))) / ELLINT_OCRT(3.0 * rerr);
    Am = A0;
    xm = x;
    ym = y;
    zm = z;
    fterm = Q;
    m = 0;

    while ( 1 )
    {
	EllInt_Num_t xrm, yrm, zrm, lam;

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

        Am = (Am + lam) * 0.25;
        xm = (xm + lam) * 0.25;
        ym = (ym + lam) * 0.25;
        zm = (zm + lam) * 0.25;
        xxm *= 0.25;
        yym *= 0.25;
        fterm *= 0.25;

        m += 1;
    }
    xxm /= Am;
    yym /= Am;
    zzm = -(xxm + yym);
    e2 = xxm * yym - zzm * zzm;
    e3 = xxm * yym * zzm;
    s = (1.0 + e2 * (e3 * ELLINT_RF_c23 - 0.1 +
		     e2 * (ELLINT_RF_c22 + e3 * ELLINT_RF_c223 +
			   e2 * ELLINT_RF_c222)) +
	 e3 * (ELLINT_RF_c03 + e3 * ELLINT_RF_c33)) / sqrt(Am);

    *res = s;
    status = ELLINT_STATUS_SUCCESS;
    return status;
}


#endif  /* _ELLINT_RF_GENERIC_GUARD */
