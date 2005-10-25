/* ========================================================================== */
/* === AMD_control ========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* AMD Version 1.0 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A. Davis,   */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README for License.          */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/amd                           */
/* -------------------------------------------------------------------------- */

/* User-callable.  Prints the control parameters for AMD.  See amd.h
 * for details.  If the Control array is not present, the defaults are
 * printed instead.
 */

#include "amd_internal.h"

GLOBAL void AMD_control
(
    double Control [ ]
)
{
    double alpha ;
    Int aggressive ;

    if (Control != (double *) NULL)
    {
	alpha = Control [AMD_DENSE] ;
	aggressive = Control [AMD_AGGRESSIVE] != 0 ;
    }
    else
    {
	alpha = AMD_DEFAULT_DENSE ;
	aggressive = AMD_DEFAULT_AGGRESSIVE ;
    }

    PRINTF (("\namd:  approximate minimum degree ordering, parameters:\n"
	"    dense row parameter: %g\n", alpha)) ;

    if (alpha < 0)
    {
	PRINTF (("    no rows treated as dense\n")) ;
    }
    else
    {
	PRINTF ((
	"    (rows with more than max (%g * sqrt (n), 16) entries are\n"
	"    considered \"dense\", and placed last in output permutation)\n",
	alpha)) ;
    }

    if (aggressive)
    {
	PRINTF (("    aggressive absorption:  yes\n\n")) ;
    }
    else
    {
	PRINTF (("    aggressive absorption:  no\n\n")) ;
    }
}
