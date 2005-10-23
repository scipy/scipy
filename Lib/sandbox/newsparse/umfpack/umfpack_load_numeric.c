/* ========================================================================== */
/* === UMFPACK_load_numeric ================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Loads a Numeric object from a file created by
    umfpack_*_save_numeric.
*/

#include "umf_internal.h"
#include "umf_valid_numeric.h"
#include "umf_malloc.h"
#include "umf_free.h"

#define READ(object,type,n) \
{ \
    object = (type *) UMF_malloc (n, sizeof (type)) ; \
    if (object == (type *) NULL) \
    { \
	UMFPACK_free_numeric ((void **) &Numeric) ; \
	fclose (f) ; \
	return (UMFPACK_ERROR_out_of_memory) ; \
    } \
    if (fread (object, sizeof (type), n, f) != n) \
    { \
	UMFPACK_free_numeric ((void **) &Numeric) ; \
	fclose (f) ; \
	return (UMFPACK_ERROR_file_IO) ; \
    } \
    if (ferror (f)) \
    { \
	UMFPACK_free_numeric ((void **) &Numeric) ; \
	fclose (f) ; \
	return (UMFPACK_ERROR_file_IO) ; \
    } \
}

/* ========================================================================== */
/* === UMFPACK_load_numeric ================================================= */
/* ========================================================================== */

GLOBAL Int UMFPACK_load_numeric
(
    void **NumericHandle,
    char *user_filename
)
{
    NumericType *Numeric ;
    char *filename ;
    FILE *f ;

    *NumericHandle = (void *) NULL ;

    /* ---------------------------------------------------------------------- */
    /* get the filename, or use the default name if filename is NULL */
    /* ---------------------------------------------------------------------- */

    if (user_filename == (char *) NULL)
    {
	filename = "numeric.umf" ;
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
    /* read the Numeric header from the file, in binary */
    /* ---------------------------------------------------------------------- */

    Numeric = (NumericType *) UMF_malloc (1, sizeof (NumericType)) ;
    if (Numeric == (NumericType *) NULL)
    {
	fclose (f) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }
    if (fread (Numeric, sizeof (NumericType), 1, f) != 1)
    {
	(void) UMF_free ((void *) Numeric) ;
	fclose (f) ;
	return (UMFPACK_ERROR_file_IO) ;
    }
    if (ferror (f))
    {
	(void) UMF_free ((void *) Numeric) ;
	fclose (f) ;
	return (UMFPACK_ERROR_file_IO) ;
    }

    if (Numeric->valid != NUMERIC_VALID || Numeric->n_row <= 0 ||
	Numeric->n_col <= 0 || Numeric->npiv < 0 || Numeric->ulen < 0 ||
	Numeric->size <= 0)
    {
	/* Numeric does not point to a NumericType object */
	(void) UMF_free ((void *) Numeric) ;
	fclose (f) ;
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }

    Numeric->D        = (Entry *) NULL ;
    Numeric->Rperm    = (Int *) NULL ;
    Numeric->Cperm    = (Int *) NULL ;
    Numeric->Lpos     = (Int *) NULL ;
    Numeric->Lilen    = (Int *) NULL ;
    Numeric->Lip      = (Int *) NULL ;
    Numeric->Upos     = (Int *) NULL ;
    Numeric->Uilen    = (Int *) NULL ;
    Numeric->Uip      = (Int *) NULL ;
    Numeric->Rs       = (double *) NULL ;
    Numeric->Memory   = (Unit *) NULL ;
    Numeric->Upattern = (Int *) NULL ;

    /* umfpack_free_numeric can now be safely called if an error occurs */

    /* ---------------------------------------------------------------------- */
    /* read the rest of the Numeric object */
    /* ---------------------------------------------------------------------- */

    READ (Numeric->D,     Entry, MIN (Numeric->n_row, Numeric->n_col)+1) ;
    READ (Numeric->Rperm, Int,   Numeric->n_row+1) ;
    READ (Numeric->Cperm, Int,   Numeric->n_col+1) ;
    READ (Numeric->Lpos,  Int,   Numeric->npiv+1) ;
    READ (Numeric->Lilen, Int,   Numeric->npiv+1) ;
    READ (Numeric->Lip,   Int,   Numeric->npiv+1) ;
    READ (Numeric->Upos,  Int,   Numeric->npiv+1) ;
    READ (Numeric->Uilen, Int,   Numeric->npiv+1) ;
    READ (Numeric->Uip,   Int,   Numeric->npiv+1) ;
    if (Numeric->scale != UMFPACK_SCALE_NONE)
    {
	READ (Numeric->Rs, double, Numeric->n_row) ;
    }
    if (Numeric->ulen > 0)
    {
	READ (Numeric->Upattern, Int, Numeric->ulen+1) ;
    }
    READ (Numeric->Memory, Unit, Numeric->size) ;

    /* close the file */
    fclose (f) ;

    /* make sure the Numeric object is valid */
    if (!UMF_valid_numeric (Numeric))
    {
	UMFPACK_free_numeric ((void **) &Numeric) ;
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }

    *NumericHandle = (void *) Numeric ;
    return (UMFPACK_OK) ;
}
