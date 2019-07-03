#ifndef _ELLINT_RD_GENERIC_GUARD
#define _ELLINT_RD_GENERIC_GUARD


EllInt_Status_t
ELLINT_POLY_FCN(ellint_RD) (EllInt_Num_t x, EllInt_Num_t y, EllInt_Num_t z,
                            double rerr, EllInt_Num_t * restrict res)
{
    EllInt_Status_t status;
    unsigned int m;
    EllInt_Num_t A0;
    EllInt_Num_t xm, ym, zm, Am;
    EllInt_Num_t adt, tmp, t;
    EllInt_Num_t e2, e3, e4, e5, xxm, yym, zzm, xy, zz2;
    double fterm, d4m, aAm;
    double adt_r, adt_i, ade_r, ade_i;
    /* Workspaces */
    EllInt_Num_t cct1[6];
    EllInt_Num_t cct2[6];

    status = ELLINT_STATUS_SUCCESS;
    if ( ELLINT_BAD_RERR(rerr, 1.0e-4) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_RERR);
    }

    if ( ( Z_INFTY(x) || Z_INFTY(y) || Z_INFTY(z) ) &&
         ( ph_is_not_pm_pi(x) && ph_is_not_pm_pi(y) && ph_is_not_pm_pi(z) ) )
    {
	ELLINT_RETURN_WITH(ELLINT_STATUS_SUCCESS, CZERO);
    }

    if ( too_small(FABS(z)) || PH_IS_PMPI_Z(z) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_ARGS);
    }

    if ( !(ph_is_not_pm_pi(x) && ph_is_not_pm_pi(y)) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_ARGS);
    }

    if ( too_small(FABS(x)) && too_small(FABS(y)) )
    {
	ELLINT_RETURN_WITH(ELLINT_STATUS_SINGULAR, CPINF);
    }

    ade_r = ade_i = 0.0;
    /* A0 = DIVcr(ADD(ADD(x, y), MULcr(z, 3.0)), 5.0); */
    cct1[0] = x;
    cct1[1] = y;
    cct1[2] = z;
    cct1[3] = z;
    cct1[4] = z;
    A0 = DIVcr(ELLINT_POLY_FCN(sum2)(cct1, 5), 5.0);
    xxm = SUB(A0, x);
    yym = SUB(A0, y);
    rerr = ELLINT_OCRT(rerr / 5.0);
    fterm = ELLINT_FMAX3(FABS(xxm), FABS(yym), FABS(SUB(A0, z))) / rerr;
    Am = A0;
    d4m = 1.0;
    adt_r = 0.0;
    adt_i = 0.0;
    xm = x;
    ym = y;
    zm = z;
    m = 0;

    while ( (aAm = FABS(Am)) <= fterm ||
	    aAm <= ELLINT_FMAX3(FABS(Am - xm), FABS(Am - ym), FABS(Am - zm)) ) 
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
	{
	    double ts, te;

	    tmp = DIVrc(d4m, MULcc(cct1[2], ADD(zm, lam)));
	    feft_sum(CREAL(tmp), adt_r, &ts, &te);
	    adt_r = ts;
	    ade_r += te;
	    feft_sum(CIMAG(tmp), adt_i, &ts, &te);
	    adt_i = ts;
	    ade_i += te;
	}

        Am = MULcr(ADD(Am, lam), 0.25);
        xm = MULcr(ADD(xm, lam), 0.25);
        ym = MULcr(ADD(ym, lam), 0.25);
        zm = MULcr(ADD(zm, lam), 0.25);
        fterm *= 0.25;
        d4m *= 0.25;
        xxm = MULcr(xxm, 0.25);
        yym = MULcr(yym, 0.25);

        m += 1;
    }

    adt = MKCMPLX(adt_r + ade_r, adt_i + ade_i);

    /* Burn some extra cycles re-balancing Am as the "true" centroid */
    cct1[0] = xm;
    cct1[1] = ym;
    cct1[2] = zm;
    cct1[3] = zm;
    cct1[4] = zm;
    Am = DIVcr(ELLINT_POLY_FCN(sum2)(cct1, 5), 5.0);
    xxm = DIVcc(xxm, Am);
    yym = DIVcc(yym, Am);
    zzm = DIVcr(ADD(xxm, yym), -3.0);
    xy = MULcc(xxm, yym);
    zz2 = MULcc(zzm, zzm);
    e2 = SUB(xy, MULcr(zz2, 6.0));
    e3 = MULcc(SUB(MULcr(xy, 3.0), MULcr(zz2, 8.0)), zzm);
    e4 = MULcr(MULcc(SUB(xy, zz2), zz2), 3.0);
    e5 = MULcc(MULcc(xy, zz2), zzm);
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
#ifdef ELLINT_POLY_REAL
    t = 1.0 + fdot2(cct1, cct2, 6) / ELLINT_RDJ_DENOM;
#else
    t = ADDcr(DIVcr(cdot2(cct1, cct2, 12), ELLINT_RDJ_DENOM), 1.0);
#endif
    tmp = MULcc(tmp, t);
    tmp = ADD(tmp, MULcc(adt, 3.0));

    *res = tmp;
    return status;
}


#endif  /* _ELLINT_RD_GENERIC_GUARD */
