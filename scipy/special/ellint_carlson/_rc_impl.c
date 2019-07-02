#ifndef _ELLINT_RC_GENERIC_GUARD
#define _ELLINT_RC_GENERIC_GUARD

#include "ellint_common.h"
#include "_arithmetic.c"


static const double ELLINT_RC_C[8] = {
    80080,
    0,
    24024,	/* 3 / 10 */
    11440,	/* 1 / 7 */
    30030,	/* 3 / 8 */
    32760,	/* 9 / 22 */
    61215,	/* 159 / 208 */
    90090	/* 9 / 8 */
};


EllInt_Status_t
ELLINT_POLY_FCN(ellint_RC) (EllInt_Num_t x, EllInt_Num_t y,
	                    double rerr, EllInt_Num_t * restrict res)
{
    EllInt_Status_t status;
    unsigned int m;
    EllInt_Num_t A0;
    EllInt_Num_t xm, ym, sm, Am;
    double fterm;

    status = ELLINT_STATUS_SUCCESS;
    if ( ELLINT_BAD_RERR(rerr, 2.0e-4) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_RERR);
    }

    if ( ( too_small(CIMAG(y)) ) && ELLINT_RNEG(y) )
    {
	/* Cauchy principal value with negative real y */
	EllInt_Num_t tmpres;

	status = ELLINT_POLY_FCN(ellint_RC)(SUB(x, y), NEG(y), rerr, &tmpres);

	if ( status == ELLINT_STATUS_SUCCESS )
	{
	    *res = MULcc(tmpres, SQRT(DIVcc(x, SUB(x, y))));
	}
	return status;
    } else if ( Z_INFTY(x) && ph_is_not_pm_pi(x) && ( !too_small(FABS(y)) ) ) {
	ELLINT_RETURN_WITH(ELLINT_STATUS_SUCCESS, CZERO);
    } else if ( ELLINT_RNEG(x) || ( too_small(FABS(y)) ) ) {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_ARGS);
    }

    A0 = DIVcr(ADD(x, MULcr(y, 2.0)), 3.0);
    fterm = FABS(SUB(A0, x)) / ELLINT_OCRT(3.0 * rerr);

    Am = A0;
    xm = x;
    ym = y;
    sm = SUB(y, A0);

    m = 0;
    while ( fmax(5.0 * FABS(SUB(xm, ym)), fterm) >= FABS(Am) )
    {
	EllInt_Num_t lam;

	if ( m > ELLINT_MAXITER )
	{
	    status = ELLINT_STATUS_NITER;
	    break;
	}

	lam = ADD(MULcr(MULcc(SQRT(xm), SQRT(ym)), 2.0), ym);

        Am = MULcr(ADD(Am, lam), 0.25);
        xm = MULcr(ADD(xm, lam), 0.25);
        ym = MULcr(ADD(ym, lam), 0.25);
	sm = MULcr(sm, 0.25);
	fterm *= 0.25;

	m += 1;
    }
    Am = DIVcr(ADD(xm, ADD(ym, ym)),  3.0);
    sm = DIVcc(sm, Am);
    *res = DIVcc(HORNER(sm, ELLINT_RC_C, 7), MULcr(SQRT(Am), ELLINT_RC_C[0]));
    return status;
}


#endif  /* _ELLINT_RC_GENERIC_GUARD */
