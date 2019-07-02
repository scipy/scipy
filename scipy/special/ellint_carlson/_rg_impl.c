#ifndef _ELLINT_RG_GENERIC_GUARD
#define _ELLINT_RG_GENERIC_GUARD

#include "ellint_common.h"


EllInt_Status_t
ELLINT_POLY_FCN(ellint_RG) (EllInt_Num_t x, EllInt_Num_t y, EllInt_Num_t z,
	                    double rerr, EllInt_Num_t * restrict res)
{
    EllInt_Status_t status;
    EllInt_Num_t rfv, rdv, tmp;

    if ( ELLINT_BAD_RERR(rerr, 1.0e-4) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_RERR);
    }

    /* Select pivot as the last argument to the symmetric functions */
    {
	double ta1, ta2;

	ta1 = FABS(x);
	if ( (ta2 = FABS(y)) < ta1 )
	{
	    if ( FABS(z) < ta1 )
	    {
		ELLINT_SWAP(x, z);
	    }
	} else if ( FABS(z) < ta2 ) {
	    ELLINT_SWAP(y, z);
	}
    }
    /* Special case -- also covers the case of z ~ zero. */
    if ( too_small(FABS(x)) && too_small(FABS(y)) )
    {
	*res = MULcr(SQRT(z), 0.5);
	return ELLINT_STATUS_SUCCESS;
    }

    status = ELLINT_POLY_FCN(ellint_RF)(x, y, z, rerr * 0.5, &rfv);
    if ( status != ELLINT_STATUS_SUCCESS )
    {
	ELLINT_FAIL_WITH(status);
    }

    status = ELLINT_POLY_FCN(ellint_RD)(x, y, z, rerr * 0.5, &rdv);
    if ( status != ELLINT_STATUS_SUCCESS )
    {
	ELLINT_FAIL_WITH(status);
    }

    tmp = MULcc(z, rfv);
    tmp = ADD(tmp,
              DIVcr(MULcc(MULcc(SUB(z, x), SUB(z, y)), rdv),
		    -3.0));
    tmp = ADD(tmp, SQRT(DIVcc(MULcc(x, y), z)));
    tmp = MULcr(tmp, 0.5);

    *res = tmp;
    status = ELLINT_STATUS_SUCCESS;
    return status;
}


#endif  /* _ELLINT_RG_GENERIC_GUARD */
