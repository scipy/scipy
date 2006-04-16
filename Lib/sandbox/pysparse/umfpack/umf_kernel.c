/* ========================================================================== */
/* === UMF_kernel =========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Primary factorization routine.   Called by UMFPACK_numeric.
    Returns:
	UMFPACK_OK if successful,
	UMFPACK_ERROR_out_of_memory if out of memory, or
	UMFPACK_ERROR_different_pattern if pattern of matrix (Ap and/or Ai)
	   has changed since the call to UMFPACK_*symbolic.
*/

#include "umf_internal.h"
#include "umf_kernel_init.h"
#include "umf_init_front.h"
#include "umf_start_front.h"
#include "umf_assemble.h"
#include "umf_scale_column.h"
#include "umf_local_search.h"
#include "umf_create_element.h"
#include "umf_extend_front.h"
#include "umf_blas3_update.h"
#include "umf_store_lu.h"
#include "umf_kernel_wrapup.h"

/* perform an action, and return if out of memory */
#define DO(action) { if (! (action)) { return (UMFPACK_ERROR_out_of_memory) ; }}

GLOBAL Int UMF_kernel
(
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    NumericType *Numeric,
    WorkType *Work,
    SymbolicType *Symbolic
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int j, f1, f2, chain, nchains, *Chain_start, status, fixQ, evaporate,
	*Front_npivcol, jmax, nb ;

    /* ---------------------------------------------------------------------- */
    /* initialize memory space and load the matrix. Optionally scale. */
    /* ---------------------------------------------------------------------- */

    if (!UMF_kernel_init (Ap, Ai, Ax,
#ifdef COMPLEX
	Az,
#endif
	Numeric, Work, Symbolic))
    {
	/* UMF_kernel_init is guaranteed to succeed, since UMFPACK_numeric */
	/* either allocates enough space or if not, UMF_kernel does not get */
	/* called.  So running out of memory here is a fatal error, and means */
	/* that the user changed Ap and/or Ai since the call to */
	/* UMFPACK_*symbolic. */
	DEBUGm4 (("kernel init failed\n")) ;
	return (UMFPACK_ERROR_different_pattern) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get the symbolic factorization */
    /* ---------------------------------------------------------------------- */

    nchains = Symbolic->nchains ;
    Chain_start = Symbolic->Chain_start ;
    Front_npivcol = Symbolic->Front_npivcol ;
    nb = Symbolic->nb ;
    fixQ = Symbolic->fixQ ;

#ifndef NDEBUG
    for (chain = 0 ; chain < nchains ; chain++)
    {
	Int i ;
	f1 = Chain_start [chain] ;
	f2 = Chain_start [chain+1] - 1 ;
	DEBUG1 (("\nCHain: "ID" start "ID" end "ID"\n", chain, f1, f2)) ;
	for (i = f1 ; i <= f2 ; i++)
	{
	    DEBUG1 (("Front "ID", npivcol "ID"\n", i, Front_npivcol [i])) ;
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* factorize each chain of frontal matrices */
    /* ---------------------------------------------------------------------- */

    for (chain = 0 ; chain < nchains ; chain++)
    {
	f1 = Chain_start [chain] ;
	f2 = Chain_start [chain+1] - 1 ;

	/* ------------------------------------------------------------------ */
	/* get the initial frontal matrix size for this chain */
	/* ------------------------------------------------------------------ */

	DO (UMF_start_front (chain, Numeric, Work, Symbolic)) ;

	/* ------------------------------------------------------------------ */
	/* factorize each front in the chain */
	/* ------------------------------------------------------------------ */

	for (Work->frontid = f1 ; Work->frontid <= f2 ; Work->frontid++)
	{

	    /* -------------------------------------------------------------- */
	    /* Initialize the pivot column candidate set  */
	    /* -------------------------------------------------------------- */

	    Work->ncand = Front_npivcol [Work->frontid] ;
	    Work->lo = Work->nextcand ;
	    Work->hi = Work->nextcand + Work->ncand - 1 ;
	    jmax = MIN (MAX_CANDIDATES, Work->ncand) ;
	    DEBUGm1 ((">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Starting front "
		ID", npivcol: "ID"\n", Work->frontid, Work->ncand)) ;
	    if (fixQ)
	    {
		/* do not modify the column order */
		jmax = 1 ;
	    }
	    DEBUGm1 (("Initial candidates: ")) ;
	    for (j = 0 ; j < jmax ; j++)
	    {
		DEBUGm1 ((" "ID, Work->nextcand)) ;
		ASSERT (Work->nextcand <= Work->hi) ;
		Work->Candidates [j] = Work->nextcand++ ;
	    }
	    Work->nCandidates = jmax ;
	    DEBUGm1 (("\n")) ;

	    /* -------------------------------------------------------------- */
	    /* Assemble and factorize the current frontal matrix */
	    /* -------------------------------------------------------------- */

	    while (Work->ncand > 0)
	    {

		/* ---------------------------------------------------------- */
		/* get the pivot row and column */
		/* ---------------------------------------------------------- */

		status = UMF_local_search (Numeric, Work, Symbolic) ;
		if (status == UMFPACK_ERROR_different_pattern)
		{
		    /* :: pattern change detected in umf_local_search :: */
		    /* input matrix has changed since umfpack_*symbolic */
		    DEBUGm4 (("local search failed\n")) ;
		    return (UMFPACK_ERROR_different_pattern) ;
		}
		if (status == UMFPACK_WARNING_singular_matrix)
		{
		    /* no pivot found, discard and try again */
		    continue ;
		}

		/* ---------------------------------------------------------- */
		/* update if front not extended or too many zeros in L,U */
		/* ---------------------------------------------------------- */

		if (Work->do_update)
		{
		    UMF_blas3_update (Work) ;
		    DO (UMF_store_lu (Numeric, Work)) ;
		}

		/* ---------------------------------------------------------- */
		/* extend the frontal matrix, or start a new one */
		/* ---------------------------------------------------------- */

		if (Work->do_extend)
		{
		    /* extend the current front */
		    DO (UMF_extend_front (Numeric, Work)) ;
		}
		else
		{
		    /* finish the current front (if any) and start a new one */
		    DO (UMF_create_element (Numeric, Work, Symbolic)) ;
		    DO (UMF_init_front (Numeric, Work)) ;
		}

		/* ---------------------------------------------------------- */
		/* Numerical & symbolic assembly into current frontal matrix */
		/* ---------------------------------------------------------- */

		if (fixQ)
		{
		    UMF_assemble_fixq (Numeric, Work) ;
		}
		else
		{
		    UMF_assemble (Numeric, Work) ;
		}

		/* ---------------------------------------------------------- */
		/* scale the pivot column */
		/* ---------------------------------------------------------- */

		UMF_scale_column (Numeric, Work) ;

		/* ---------------------------------------------------------- */
		/* Numerical update if enough pivots accumulated */
		/* ---------------------------------------------------------- */

		evaporate = Work->fnrows == 0 || Work->fncols == 0 ;
		if (Work->fnpiv >= nb || evaporate)
		{
		    UMF_blas3_update (Work) ;
		    DO (UMF_store_lu (Numeric, Work)) ;
		}

		Work->pivrow_in_front = FALSE ;
		Work->pivcol_in_front = FALSE ;

		/* ---------------------------------------------------------- */
		/* If front is empty, evaporate it */
		/* ---------------------------------------------------------- */

		if (evaporate)
		{
		    /* This does not create an element, just evaporates it.
		     * It ensures that a front is not 0-by-c or r-by-0.  No
		     * memory is allocated, so it is guaranteed to succeed. */
		    (void) UMF_create_element (Numeric, Work, Symbolic) ;
		    Work->fnrows = 0 ;
		    Work->fncols = 0 ;
		}
	    }
	}

	/* ------------------------------------------------------------------
	 * Wrapup the current frontal matrix.  This is the last in a chain
	 * in the column elimination tree.  The next frontal matrix
	 * cannot overlap with the current one, which will be its sibling
	 * in the column etree.
	 * ------------------------------------------------------------------ */

	UMF_blas3_update (Work) ;
	DO (UMF_store_lu (Numeric, Work)) ;
	Work->fnrows_new = Work->fnrows ;
	Work->fncols_new = Work->fncols ;
	DO (UMF_create_element (Numeric, Work, Symbolic)) ;

	/* ------------------------------------------------------------------ */
	/* current front is now empty */
	/* ------------------------------------------------------------------ */

	Work->fnrows = 0 ;
	Work->fncols = 0 ;
    }

    /* ---------------------------------------------------------------------- */
    /* end the last Lchain and Uchain and finalize the LU factors */
    /* ---------------------------------------------------------------------- */

    UMF_kernel_wrapup (Numeric, Symbolic, Work) ;

    /* note that the matrix may be singular (this is OK) */
    return (UMFPACK_OK) ;
}
