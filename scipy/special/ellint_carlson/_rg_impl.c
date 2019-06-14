#ifndef _ELLINT_RG_GENERIC_GUARD
#define _ELLINT_RG_GENERIC_GUARD

#include "ellint_common.h"
#include "ellint_poly.h"
#include "ellint_argcheck.h"

#define ELLINT_RG_m1over3	(-0.3333333333333333333333333)


int
ELLINT_POLY_FCN(ellint_RG) (EllInt_Num_t x, EllInt_Num_t y, EllInt_Num_t z,
	                    double rerr, EllInt_Num_t * restrict res)
{
    int status;

    if ( ELLINT_BAD_RERR(rerr, 1.0e-4) )
    {
	ELLINT_FAIL_WITH(ELLINT_STATUS_BAD_RERR);
    }

    /* Select pivot as the last argument to the symmetric functions */
    {
	double ta1, ta2;
	ta1 = fabs(x);
	if ( (ta2 = fabs(y)) < ta1 )
	{
	    if ( fabs(z) < ta1 )
	    {
		ELLINT_SWAP(x, z);
	    }
	} else if ( fabs(z) < ta2 ) {
	    ELLINT_SWAP(y, z);
	}
    }
    /* Special case -- also covers the case of z ~ zero. */
    if ( too_small(fabs(x)) && too_small(fabs(y)) )
    {
	*res = (0.5 * sqrt(z));
	return ELLINT_STATUS_SUCCESS;
    }

    EllInt_Num_t rfv, rdv, tmp;

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
    tmp = z * rfv;
    tmp += ELLINT_RG_m1over3 * (z - x) * (z - y) * rdv;
    tmp += sqrt(x * y / z);
    tmp *= 0.5;

    *res = tmp;
    status = ELLINT_STATUS_SUCCESS;
    return status;
}


#endif  /* _ELLINT_RG_GENERIC_GUARD */
