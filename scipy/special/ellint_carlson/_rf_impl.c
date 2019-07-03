#ifndef _ELLINT_RF_GENERIC_GUARD
#define _ELLINT_RF_GENERIC_GUARD

#include <stdlib.h>
#include "ellint_common.h"


static const double ELLINT_RF_C1[4] = {0, -24024, 10010, -5775};
static const double ELLINT_RF_C2[3] = {17160, -16380, 15015};
#define ELLINT_RF_c33	(+6930)		/*  3 / 104 */
#define ELLINT_RF_DENOM	(240240.0)


static int ELLINT_POLY_FCN(rabscmp)(const void * a, const void * b)
{
    if ( (FABS( *((const EllInt_Num_t *)a) ) <
          FABS( *((const EllInt_Num_t *)b) )) )
    {
	return -1;
    }
    return 1;
}


static EllInt_Status_t ELLINT_POLY_FCN(ellint_RF0)(EllInt_Num_t x,
						   EllInt_Num_t y,
						   double rerr,
						   EllInt_Num_t * restrict res)
{
    EllInt_Status_t status;
    unsigned int m;
    EllInt_Num_t xm, ym;
    double rsq;

    status = ELLINT_STATUS_SUCCESS;

    if ( ELLINT_BAD_RERR(rerr, 1.0e-4) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_RERR);
    }

    rsq = 2.0 * sqrt(rerr);

    xm = SQRT(x);
    ym = SQRT(y);
    m = 0;
    while ( FABS(SUB(xm, ym)) >= rsq * fmin(FABS(xm), FABS(ym)) )
    {
	EllInt_Num_t xnext, ynext;

	if ( m > ELLINT_MAXITER )
	{
	    status = ELLINT_STATUS_NITER;
	    break;
	}

	xnext = MULcr(ADD(xm, ym), 0.5);
	ynext = SQRT(MULcc(xm, ym));

	xm = xnext;
	ym = ynext;
	m += 1;
    }

    *res = DIVrc(M_PI, ADD(xm, ym));
    return status;
}


EllInt_Status_t
ELLINT_POLY_FCN(ellint_RF) (EllInt_Num_t x, EllInt_Num_t y, EllInt_Num_t z,
	                    double rerr, EllInt_Num_t * restrict res)
{
    EllInt_Status_t status;
    unsigned int m;
    EllInt_Num_t A0;
    EllInt_Num_t xm, ym, zm, Am;
    EllInt_Num_t e2, e3, xxm, yym, zzm, s;
    double fterm;
    EllInt_Num_t cct1[3];
    EllInt_Num_t cct2[3];

    status = ELLINT_STATUS_SUCCESS;
    if ( ELLINT_BAD_RERR(rerr, 3.0e-4) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_RERR);
    }

    if ( ( Z_INFTY(x) || Z_INFTY(y) || Z_INFTY(z) ) &&
         ( ph_is_not_pm_pi(x) && ph_is_not_pm_pi(y) && ph_is_not_pm_pi(z) ) )
    {
	ELLINT_RETURN_WITH(ELLINT_STATUS_SUCCESS, CZERO);
    }

    if ( PH_IS_PMPI_Z(x) || PH_IS_PMPI_Z(y) || PH_IS_PMPI_Z(z) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_ARGS);
    }

    if ( !(C_ATMOST_1Z(x, y, z)) )
    {
	ELLINT_RETURN_WITH(ELLINT_STATUS_SINGULAR, CPINF);
    }

    /* A0 = DIVcr(ADD(ADD(x, y), z), 3.0); */
    cct1[0] = x;
    cct1[1] = y;
    cct1[2] = z;
    qsort(cct1, 3, sizeof (EllInt_Num_t), ELLINT_POLY_FCN(rabscmp));
    x = cct1[0];
    y = cct1[1];
    z = cct1[2];
    if ( too_small(FABS(x)) )
    {
	EllInt_Num_t tmpres;
	status = ELLINT_POLY_FCN(ellint_RF0)(y, z, rerr, &tmpres);
	*res = SUB(tmpres, DIVcc(SQRT(x), SQRT(MULcc(y, z))));
	return status;
    }
    A0 = DIVcr(ELLINT_POLY_FCN(sum2)(cct1, 3), 3.0);
    xxm = SUB(A0, x);
    yym = SUB(A0, y);
    rerr = ELLINT_OCRT(3.0 * rerr);
    fterm = ELLINT_FMAX3(FABS(xxm), FABS(yym), FABS(SUB(A0, z))) / rerr;
    Am = A0;
    xm = x;
    ym = y;
    zm = z;

    m = 0;
    while ( fmax(fterm,
                 ELLINT_FMAX3(FABS(Am - xm), FABS(Am - ym),
			      FABS(Am - zm))) >= FABS(Am) )
    {
	EllInt_Num_t lam;

	if ( m > ELLINT_MAXITER )
	{
	    status = ELLINT_STATUS_NITER;
	    break;
	}

        /* lam = ADD(ADD(MULcc(xrm, yrm), MULcc(xrm, zrm)), MULcc(yrm, zrm)); */
	cct1[0] = cct2[2] = SQRT(xm);
	cct1[1] = cct2[0] = SQRT(ym);
	cct1[2] = cct2[1] = SQRT(zm);
#ifdef ELLINT_POLY_REAL
	lam = fdot2(cct1, cct2, 3);
#else
	lam = cdot2(cct1, cct2, 6);
#endif
        Am = MULcr(ADD(Am, lam), 0.25);
        xm = MULcr(ADD(xm, lam), 0.25);
        ym = MULcr(ADD(ym, lam), 0.25);
        zm = MULcr(ADD(zm, lam), 0.25);
        xxm = MULcr(xxm, 0.25);
        yym = MULcr(yym, 0.25);
        fterm *= 0.25;

        m += 1;
    }

    /* Burn some extra cycles re-balancing Am as the "true" centroid */
    cct1[0] = xm;
    cct1[1] = ym;
    cct1[2] = zm;
    Am = DIVcr(ELLINT_POLY_FCN(sum2)(cct1, 3), 3.0);
    xxm = DIVcc(xxm, Am);
    yym = DIVcc(yym, Am);
    zzm = NEG(ADD(xxm, yym));
    e2 = SUB(MULcc(xxm, yym), MULcc(zzm, zzm));
    e3 = MULcc(xxm, MULcc(yym, zzm));
    s = HORNER(e2, ELLINT_RF_C1, 3);
    s = ADD(s,
	    MULcc(e3, ADD(HORNER(e2, ELLINT_RF_C2, 2),
                          MULcr(e3, ELLINT_RF_c33))));
    s = ADDcr(DIVcr(s, ELLINT_RF_DENOM), 1.0);

    *res = DIVcc(s, SQRT(Am));
    return status;
}


#endif  /* _ELLINT_RF_GENERIC_GUARD */
