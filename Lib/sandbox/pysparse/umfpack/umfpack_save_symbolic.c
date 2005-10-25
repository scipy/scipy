/* ========================================================================== */
/* === UMFPACK_save_symbolic ================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Saves a Symbolic object to a file.  It can later be read
    back in via a call to umfpack_*_load_symbolic.
*/

#include "umf_internal.h"
#include "umf_valid_symbolic.h"

#define WRITE(object,type,n) \
{ \
    ASSERT (object != (type *) NULL) ; \
    if (fwrite (object, sizeof (type), n, f) != n) \
    { \
	fclose (f) ; \
	return (UMFPACK_ERROR_file_IO) ; \
    } \
}

/* ========================================================================== */
/* === UMFPACK_save_symbolic ================================================ */
/* ========================================================================== */

GLOBAL Int UMFPACK_save_symbolic
(
    void *SymbolicHandle,
    char *user_filename
)
{
    SymbolicType *Symbolic ;
    char *filename ;
    FILE *f ;

    /* get the Symbolic object */
    Symbolic = (SymbolicType *) SymbolicHandle ;

    /* make sure the Symbolic object is valid */
    if (!UMF_valid_symbolic (Symbolic))
    {
	return (UMFPACK_ERROR_invalid_Symbolic_object) ;
    }

    /* get the filename, or use the default name if filename is NULL */
    if (user_filename == (char *) NULL)
    {
	filename = "symbolic.umf" ;
    }
    else
    {
	filename = user_filename ;
    }
    f = fopen (filename, "wb") ;
    if (!f)
    {
	return (UMFPACK_ERROR_file_IO) ;
    }

    /* write the Symbolic object to the file, in binary */
    WRITE (Symbolic,                     SymbolicType, 1) ;
    WRITE (Symbolic->Cperm_init,         Int, Symbolic->n_col+1) ;
    WRITE (Symbolic->Rperm_init,         Int, Symbolic->n_row+1) ;
    WRITE (Symbolic->Front_npivcol,      Int, Symbolic->nfr+1) ;
    WRITE (Symbolic->Front_parent,       Int, Symbolic->nfr+1) ;
    WRITE (Symbolic->Front_1strow,       Int, Symbolic->nfr+1) ;
    WRITE (Symbolic->Front_leftmostdesc, Int, Symbolic->nfr+1) ;
    WRITE (Symbolic->Chain_start,        Int, Symbolic->nchains+1) ;
    WRITE (Symbolic->Chain_maxrows,      Int, Symbolic->nchains+1) ;
    WRITE (Symbolic->Chain_maxcols,      Int, Symbolic->nchains+1) ;
    WRITE (Symbolic->Cdeg,               Int, Symbolic->n_col+1) ;
    WRITE (Symbolic->Rdeg,               Int, Symbolic->n_row+1) ;
    if (Symbolic->esize > 0)
    {
	/* only when dense rows are present */
	WRITE (Symbolic->Esize, Int, Symbolic->esize) ;
    }
    if (Symbolic->prefer_diagonal)
    {
	/* only when diagonal pivoting is prefered */
	WRITE (Symbolic->Diagonal_map, Int, Symbolic->n_col+1) ;
    }

    /* close the file */
    fclose (f) ;

    return (UMFPACK_OK) ;
}
