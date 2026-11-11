//------------------------------------------------------------------------------
// AMD/Source/amd_control: print control parameters for AMD
//------------------------------------------------------------------------------

// AMD, Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and
// Iain S. Duff.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

/* User-callable.  Prints the control parameters for AMD.  See amd.h
 * for details.  If the Control array is not present, the defaults are
 * printed instead.
 */

#include "amd_internal.h"

#ifndef NDEBUG
#define HIGHS_ORDERING_PRINT(...) printf(__VA_ARGS__)
#else
#define HIGHS_ORDERING_PRINT(...)
#endif

void Highs_amd_control
(
    double Control [ ]
)
{
    double alpha ;
    amd_int aggressive ;

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

    HIGHS_ORDERING_PRINT (
        "\nAMD version %d.%d.%d, %s: approximate minimum degree ordering\n"
	"    dense row parameter: %g\n", AMD_MAIN_VERSION, AMD_SUB_VERSION,
	AMD_SUBSUB_VERSION, AMD_DATE, alpha);

    if (alpha < 0)
    {
	HIGHS_ORDERING_PRINT ("    no rows treated as dense\n") ;
    }
    else
    {
	HIGHS_ORDERING_PRINT (
	"    (rows with more than max (%g * sqrt (n), 16) entries are\n"
	"    considered \"dense\", and placed last in output permutation)\n",
	alpha) ;
    }

    if (aggressive)
    {
	HIGHS_ORDERING_PRINT ("    aggressive absorption:  yes\n") ;
    }
    else
    {
	HIGHS_ORDERING_PRINT ("    aggressive absorption:  no\n") ;
    }

    HIGHS_ORDERING_PRINT ("    size of AMD integer: %zu\n\n", sizeof (amd_int)) ;
}
