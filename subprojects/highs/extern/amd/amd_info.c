//------------------------------------------------------------------------------
// AMD/Source/amd_info: print output statistics for AMD
//------------------------------------------------------------------------------

// AMD, Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and
// Iain S. Duff.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

/* User-callable.  Prints the output statistics for AMD.  See amd.h
 * for details.  If the Info array is not present, nothing is printed.
 */

#include "amd_internal.h"

#ifndef NDEBUG
#define HIGHS_ORDERING_PRINT(...) printf(__VA_ARGS__)
#else
#define HIGHS_ORDERING_PRINT(...)
#endif

#define PRI(format,x) { if (x >= 0) { HIGHS_ORDERING_PRINT (format, x) ; }}

void Highs_amd_info
(
    double Info [ ]
)
{
    double n, ndiv, nmultsubs_ldl, nmultsubs_lu, lnz, lnzd ;

    HIGHS_ORDERING_PRINT ("\nAMD version %d.%d.%d, %s, results:\n",
	AMD_MAIN_VERSION, AMD_SUB_VERSION, AMD_SUBSUB_VERSION, AMD_DATE) ;

    if (!Info)
    {
	return ;
    }

    n = Info [AMD_N] ;
    ndiv = Info [AMD_NDIV] ;
    nmultsubs_ldl = Info [AMD_NMULTSUBS_LDL] ;
    nmultsubs_lu = Info [AMD_NMULTSUBS_LU] ;
    lnz = Info [AMD_LNZ] ;
    lnzd = (n >= 0 && lnz >= 0) ? (n + lnz) : (-1) ;

    /* AMD return status */
    HIGHS_ORDERING_PRINT ("    status: ") ;
    if (Info [AMD_STATUS] == AMD_OK)
    {
	HIGHS_ORDERING_PRINT ("OK\n") ;
    }
    else if (Info [AMD_STATUS] == AMD_OUT_OF_MEMORY)
    {
	HIGHS_ORDERING_PRINT ("out of memory\n") ;
    }
    else if (Info [AMD_STATUS] == AMD_INVALID)
    {
	HIGHS_ORDERING_PRINT ("invalid matrix\n") ;
    }
    else if (Info [AMD_STATUS] == AMD_OK_BUT_JUMBLED)
    {
	HIGHS_ORDERING_PRINT ("OK, but jumbled\n") ;
    }
    else
    {
	HIGHS_ORDERING_PRINT ("unknown\n") ;
    }

    /* statistics about the input matrix */
    PRI ("    n, dimension of A:                                  %.20g\n", n);
    PRI ("    nz, number of nonzeros in A:                        %.20g\n",
	Info [AMD_NZ]) ;
    PRI ("    symmetry of A:                                      %.4f\n",
	Info [AMD_SYMMETRY]) ;
    PRI ("    number of nonzeros on diagonal:                     %.20g\n",
	Info [AMD_NZDIAG]) ;
    PRI ("    nonzeros in pattern of A+A' (excl. diagonal):       %.20g\n",
	Info [AMD_NZ_A_PLUS_AT]) ;
    PRI ("    # dense rows/columns of A+A':                       %.20g\n",
	Info [AMD_NDENSE]) ;

    /* statistics about AMD's behavior  */
    PRI ("    memory used, in bytes:                              %.20g\n",
	Info [AMD_MEMORY]) ;
    PRI ("    # of memory compactions:                            %.20g\n",
	Info [AMD_NCMPA]) ;

    /* statistics about the ordering quality */
    HIGHS_ORDERING_PRINT ("\n"
	"    The following approximate statistics are for a subsequent\n"
	"    factorization of A(P,P) + A(P,P)'.  They are slight upper\n"
	"    bounds if there are no dense rows/columns in A+A', and become\n"
	"    looser if dense rows/columns exist.\n\n") ;

    PRI ("    nonzeros in L (excluding diagonal):                 %.20g\n",
	lnz) ;
    PRI ("    nonzeros in L (including diagonal):                 %.20g\n",
	lnzd) ;
    PRI ("    # divide operations for LDL' or LU:                 %.20g\n",
	ndiv) ;
    PRI ("    # multiply-subtract operations for LDL':            %.20g\n",
	nmultsubs_ldl) ;
    PRI ("    # multiply-subtract operations for LU:              %.20g\n",
	nmultsubs_lu) ;
    PRI ("    max nz. in any column of L (incl. diagonal):        %.20g\n",
	Info [AMD_DMAX]) ;

    /* total flop counts for various factorizations */

    if (n >= 0 && ndiv >= 0 && nmultsubs_ldl >= 0 && nmultsubs_lu >= 0)
    {
	HIGHS_ORDERING_PRINT ("\n"
	"    chol flop count for real A, sqrt counted as 1 flop: %.20g\n"
	"    LDL' flop count for real A:                         %.20g\n"
	"    LDL' flop count for complex A:                      %.20g\n"
	"    LU flop count for real A (with no pivoting):        %.20g\n"
	"    LU flop count for complex A (with no pivoting):     %.20g\n\n",
	n + ndiv + 2*nmultsubs_ldl,
	    ndiv + 2*nmultsubs_ldl,
	  9*ndiv + 8*nmultsubs_ldl,
	    ndiv + 2*nmultsubs_lu,
	  9*ndiv + 8*nmultsubs_lu) ;
    }
}
