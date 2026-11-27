//------------------------------------------------------------------------------
// AMD/Source/amd_defaults: set defaults for AMD
//------------------------------------------------------------------------------

// AMD, Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and
// Iain S. Duff.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

/* User-callable.  Sets default control parameters for AMD.  See amd.h
 * for details.
 */

#include "amd_internal.h"

/* ========================================================================= */
/* === AMD defaults ======================================================== */
/* ========================================================================= */

void Highs_amd_defaults
(
    double Control [ ]
)
{
    amd_int i ;

    if (Control != (double *) NULL)
    {
	for (i = 0 ; i < AMD_CONTROL ; i++)
	{
	    Control [i] = 0 ;
	}
	Control [AMD_DENSE] = AMD_DEFAULT_DENSE ;
	Control [AMD_AGGRESSIVE] = AMD_DEFAULT_AGGRESSIVE ;
    }
}
