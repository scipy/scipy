/* ========================================================================== */
/* === UMF_blas3_update ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

#include "umf_internal.h"

GLOBAL void UMF_blas3_update
(
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Entry *L, *U, *C, *LU ;
    Int k, m, n, d, nb, dc ;

    DEBUG5 (("In UMF_blas3_update "ID" "ID" "ID"\n",
	Work->fnpiv, Work->fnrows, Work->fncols)) ;

    k = Work->fnpiv ;
    if (k == 0)
    {
	/* no work to do */
	return ;
    }

    m = Work->fnrows ;
    n = Work->fncols ;

    d = Work->fnr_curr ;
    dc = Work->fnc_curr ;
    nb = Work->nb ;
    ASSERT (d >= 0 && (d % 2) == 1) ;
    C = Work->Fcblock ;	    /* ldc is fnr_curr */
    L =	Work->Flblock ;	    /* ldl is fnr_curr */
    U = Work->Fublock ;	    /* ldu is fnc_curr, stored by rows */
    LU = Work->Flublock ;   /* nb-by-nb */

#ifndef NDEBUG
    DEBUG5 (("DO RANK-NB UPDATE of frontal:\n")) ;
    DEBUG5 (("DGEMM : "ID" "ID" "ID"\n", k, m, n)) ;
    DEBUG7 (("C  block: ")) ; UMF_dump_dense (C,  d, m, n) ;
    DEBUG7 (("A  block: ")) ; UMF_dump_dense (L,  d, m, k) ;
    DEBUG7 (("B' block: ")) ; UMF_dump_dense (U, dc, n, k) ;
    DEBUG7 (("LU block: ")) ; UMF_dump_dense (LU, nb, k, k) ;
#endif

    if (k == 1)
    {

#ifdef USE_NO_BLAS

	/* no BLAS available - use plain C code instead */
	Int i, j ;

	/* rank-1 outer product to update the C block */
	for (j = 0 ; j < n ; j++)
	{
	    Entry u_j = U [j] ;
	    if (IS_NONZERO (u_j))
	    {
		Entry *c_ij, *l_is ;
		c_ij = & C [j*d] ;
		l_is = & L [0] ;
#pragma ivdep
		for (i = 0 ; i < m ; i++)
		{
		    /* C [i+j*d]-= L [i] * U [j] */
		    MULT_SUB (*c_ij, *l_is, u_j) ;
		    c_ij++ ;
		    l_is++ ;
		}
	    }
	}

#else
	BLAS_GER (m, n, L, U, C, d) ;

#endif /* USE_NO_BLAS */

    }
    else
    {

#ifdef USE_NO_BLAS

	/* no BLAS available - use plain C code instead */
	Int i, j, s ;

	/* triangular solve to update the U block */
	for (s = 0 ; s < k ; s++)
	{
	    for (i = s+1 ; i < k ; i++)
	    {
		Entry l_is = LU [i+s*nb] ;
		if (IS_NONZERO (l_is))
		{
		    Entry *u_ij, *u_sj ;
		    u_ij = & U [i*dc] ;
		    u_sj = & U [s*dc] ;
#pragma ivdep
		    for (j = 0 ; j < n ; j++)
		    {
			/* U [i*dc+j] -= LU [i+s*nb] * U [s*dc+j] ; */
			MULT_SUB (*u_ij, l_is, *u_sj) ;
			u_ij++ ;
			u_sj++ ;
		    }
		}
	    }
	}

	/* rank-k outer product to update the C block */
	/* C = C - L*U' (U is stored by rows, not columns) */
	for (s = 0 ; s < k ; s++)
	{
	    for (j = 0 ; j < n ; j++)
	    {
		Entry u_sj = U [j+s*dc] ;
		if (IS_NONZERO (u_sj))
		{
		    Entry *c_ij, *l_is ;
		    c_ij = & C [j*d] ;
		    l_is = & L [s*d] ;
#pragma ivdep
		    for (i = 0 ; i < m ; i++)
		    {
			/* C [i+j*d]-= L [i+s*d] * U [s*dc+j] */
			MULT_SUB (*c_ij, *l_is, u_sj) ;
			c_ij++ ;
			l_is++ ;
		    }
		}
	    }
	}

#else

	BLAS_TRSM_RIGHT (n, k, LU, nb, U, dc) ;
	BLAS_GEMM (m, n, k, L, U, dc, C, d) ;

#endif /* USE_NO_BLAS */

    }

#ifndef NDEBUG
    DEBUG5 (("RANK-NB UPDATE of frontal done:\n")) ;
    DEBUG5 (("DGEMM : "ID" "ID" "ID"\n", k, m, n)) ;
    DEBUG7 (("C  block: ")) ; UMF_dump_dense (C,  d, m, n) ;
    DEBUG7 (("A  block: ")) ; UMF_dump_dense (L,  d, m, k) ;
    DEBUG7 (("B' block: ")) ; UMF_dump_dense (U, dc, n, k) ;
    DEBUG7 (("LU block: ")) ; UMF_dump_dense (LU, nb, k, k) ;
#endif

    DEBUG2 (("blas3 "ID" "ID" "ID"\n", k, Work->fnrows, Work->fncols)) ;
}
