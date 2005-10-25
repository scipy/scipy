/* ========================================================================== */
/* === UMFPACK_report_vector ================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Prints a real or complex vector.
    See umfpack_report_vector.h for details.
*/

#include "umf_internal.h"
#include "umf_report_vector.h"

GLOBAL Int UMFPACK_report_vector
(
    Int n,
    const double Xx [ ],
#ifdef COMPLEX
    const double Xz [ ],
#endif
    const double Control [UMFPACK_CONTROL]
)
{
    Int prl ;

#ifndef COMPLEX
    double *Xz = (double *) NULL ;
#endif

    prl = GET_CONTROL (UMFPACK_PRL, UMFPACK_DEFAULT_PRL) ;

    if (prl <= 2)
    {
	return (UMFPACK_OK) ;
    }

    return (UMF_report_vector (n, Xx, Xz, prl, TRUE, FALSE)) ;
}
