/* ========================================================================== */
/* === UMF_malloc =========================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Allocate a block of n objects, each of a given size.  This routine does not
    handle the case when the size is 1 (allocating char's) because of potential
    integer overflow.  UMFPACK never does that.
    Also maintains the UMFPACK malloc count.
*/

#include "umf_internal.h"

#if defined (UMF_MALLOC_COUNT) || !defined (NDEBUG)

/*
    UMF_malloc_count is a count of the objects malloc'd by UMFPACK.  If you
    suspect a memory leak in your program (caused by not properly destroying
    the Symbolic and Numeric objects) then compile with -DUMF_MALLOC_COUNT and
    check value of UMF_malloc_count.  By default, UMF_MALLOC_COUNT is not
    defined, and thus UMFPACK has no global variables.
*/

GLOBAL Int UMF_malloc_count = 0 ;

#endif

#ifdef UMF_TCOV_TEST
/* For exhaustive statement coverage testing only! */
GLOBAL Int umf_fail, umf_fail_lo, umf_fail_hi ;
GLOBAL Int umf_realloc_fail, umf_realloc_lo, umf_realloc_hi ;
#endif

GLOBAL void *UMF_malloc
(
    Int n_objects,
    size_t size_of_object
)
{
    size_t size ;
    void *p ;

#ifdef UMF_TCOV_TEST
    /* For exhaustive statement coverage testing only! */
    /* Pretend to fail, to test out-of-memory conditions. */
    umf_fail-- ;
    if (umf_fail <= umf_fail_hi && umf_fail >= umf_fail_lo)
    {
	DEBUG0 (("umf_malloc: Pretend to fail %d %d %d\n",
	    umf_fail, umf_fail_hi, umf_fail_lo)) ;
	return ((void *) NULL) ;
    }
#endif

    DEBUG0 (("UMF_malloc: ")) ;

    /* make sure that we allocate something */
    n_objects = MAX (1, n_objects) ;

    size = (size_t) n_objects ;
    ASSERT (size_of_object > 1) ;
    if (size > Int_MAX / size_of_object)
    {
	/* object is too big for integer pointer arithmetic */
	return ((void *) NULL) ;
    }
    size *= size_of_object ;

    /* see umf_config.h for the memory allocator selection */
    p = ALLOCATE (size) ;

    DEBUG0 ((ID"\n", (Int) p)) ;

#if defined (UMF_MALLOC_COUNT) || !defined (NDEBUG)
    if (p)
    {
	/* One more object has been malloc'ed.  Keep track of the count. */
	/* (purely for sanity checks). */
	UMF_malloc_count++ ;
	DEBUG0 (("  successful, new malloc count: "ID"\n", UMF_malloc_count)) ;
    }
#endif

    return (p) ;
}
