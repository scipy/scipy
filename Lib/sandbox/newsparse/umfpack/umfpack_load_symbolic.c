/* ========================================================================== */
/* === UMFPACK_load_symbolic ================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Loads a Symbolic object from a file created by
    umfpack_*_save_symbolic.
*/

#include "umf_internal.h"
#include "umf_valid_symbolic.h"
#include "umf_malloc.h"
#include "umf_free.h"

#define READ(object,type,n) \
{ \
    object = (type *) UMF_malloc (n, sizeof (type)) ; \
    if (object == (type *) NULL) \
    { \
	UMFPACK_free_symbolic ((void **) &Symbolic) ; \
	fclose (f) ; \
	return (UMFPACK_ERROR_out_of_memory) ; \
    } \
    if (fread (object, sizeof (type), n, f) != n) \
    { \
	UMFPACK_free_symbolic ((void **) &Symbolic) ; \
	fclose (f) ; \
	return (UMFPACK_ERROR_file_IO) ; \
    } \
    if (ferror (f)) \
    { \
	UMFPACK_free_symbolic ((void **) &Symbolic) ; \
	fclose (f) ; \
	return (UMFPACK_ERROR_file_IO) ; \
    } \
}

/* ========================================================================== */
/* === UMFPACK_load_symbolic ================================================ */
/* ========================================================================== */

GLOBAL Int UMFPACK_load_symbolic
(
    void **SymbolicHandle,
    char *user_filename
)
{
    SymbolicType *Symbolic ;
    char *filename ;
    FILE *f ;

    *SymbolicHandle = (void *) NULL ;

    /* ---------------------------------------------------------------------- */
    /* get the filename, or use the default name if filename is NULL */
    /* ---------------------------------------------------------------------- */

    if (user_filename == (char *) NULL)
    {
	filename = "symbolic.umf" ;
    }
    else
    {
	filename = user_filename ;
    }
    f = fopen (filename, "rb") ;
    if (!f)
    {
	return (UMFPACK_ERROR_file_IO) ;
    }

    /* ---------------------------------------------------------------------- */
    /* read the Symbolic header from the file, in binary */
    /* ---------------------------------------------------------------------- */

    Symbolic = (SymbolicType *) UMF_malloc (1, sizeof (SymbolicType)) ;
    if (Symbolic == (SymbolicType *) NULL)
    {
	fclose (f) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }
    if (fread (Symbolic, sizeof (SymbolicType), 1, f) != 1)
    {
	(void) UMF_free ((void *) Symbolic) ;
	fclose (f) ;
	return (UMFPACK_ERROR_file_IO) ;
    }
    if (ferror (f))
    {
	(void) UMF_free ((void *) Symbolic) ;
	fclose (f) ;
	return (UMFPACK_ERROR_file_IO) ;
    }

    if (Symbolic->valid != SYMBOLIC_VALID || Symbolic->n_row <= 0 ||
	Symbolic->n_col <= 0 || Symbolic->nfr < 0 || Symbolic->nchains < 0 ||
	Symbolic->esize < 0)
    {
	/* Symbolic does not point to a Symbolic object */
	(void) UMF_free ((void *) Symbolic) ;
	fclose (f) ;
	return (UMFPACK_ERROR_invalid_Symbolic_object) ;
    }

    Symbolic->Cperm_init         = (Int *) NULL ;
    Symbolic->Rperm_init         = (Int *) NULL ;
    Symbolic->Front_npivcol      = (Int *) NULL ;
    Symbolic->Front_parent       = (Int *) NULL ;
    Symbolic->Front_1strow       = (Int *) NULL ;
    Symbolic->Front_leftmostdesc = (Int *) NULL ;
    Symbolic->Chain_start        = (Int *) NULL ;
    Symbolic->Chain_maxrows      = (Int *) NULL ;
    Symbolic->Chain_maxcols      = (Int *) NULL ;
    Symbolic->Cdeg               = (Int *) NULL ;
    Symbolic->Rdeg               = (Int *) NULL ;
    Symbolic->Esize              = (Int *) NULL ;
    Symbolic->Diagonal_map       = (Int *) NULL ;

    /* umfpack_free_symbolic can now be safely called if an error occurs */

    /* ---------------------------------------------------------------------- */
    /* read the rest of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    READ (Symbolic->Cperm_init,         Int, Symbolic->n_col+1) ;
    READ (Symbolic->Rperm_init,         Int, Symbolic->n_row+1) ;
    READ (Symbolic->Front_npivcol,      Int, Symbolic->nfr+1) ;
    READ (Symbolic->Front_parent,       Int, Symbolic->nfr+1) ;
    READ (Symbolic->Front_1strow,       Int, Symbolic->nfr+1) ;
    READ (Symbolic->Front_leftmostdesc, Int, Symbolic->nfr+1) ;
    READ (Symbolic->Chain_start,        Int, Symbolic->nchains+1) ;
    READ (Symbolic->Chain_maxrows,      Int, Symbolic->nchains+1) ;
    READ (Symbolic->Chain_maxcols,      Int, Symbolic->nchains+1) ;
    READ (Symbolic->Cdeg,               Int, Symbolic->n_col+1) ;
    READ (Symbolic->Rdeg,               Int, Symbolic->n_row+1) ;
    if (Symbolic->esize > 0)
    {
	/* only when dense rows are present */
	READ (Symbolic->Esize, Int, Symbolic->esize) ;
    }
    if (Symbolic->prefer_diagonal)
    {
	/* only when diagonal pivoting is prefered */
	READ (Symbolic->Diagonal_map, Int, Symbolic->n_col+1) ;
    }

    /* close the file */
    fclose (f) ;

    /* make sure the Symbolic object is valid */
    if (!UMF_valid_symbolic (Symbolic))
    {
	UMFPACK_free_symbolic ((void **) &Symbolic) ;
	return (UMFPACK_ERROR_invalid_Symbolic_object) ;
    }

    *SymbolicHandle = (void *) Symbolic ;
    return (UMFPACK_OK) ;
}
