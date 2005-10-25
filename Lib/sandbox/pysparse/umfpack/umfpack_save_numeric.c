/* ========================================================================== */
/* === UMFPACK_save_numeric ================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Saves a Numeric object to a file.  It can later be read back
    in via a call to umfpack_*_load_numeric.
*/

#include "umf_internal.h"
#include "umf_valid_numeric.h"

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
/* === UMFPACK_save_numeric ================================================= */
/* ========================================================================== */

GLOBAL Int UMFPACK_save_numeric
(
    void *NumericHandle,
    char *user_filename
)
{
    NumericType *Numeric ;
    char *filename ;
    FILE *f ;

    /* get the Numeric object */
    Numeric = (NumericType *) NumericHandle ;

    /* make sure the Numeric object is valid */
    if (!UMF_valid_numeric (Numeric))
    {
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }

    /* get the filename, or use the default name if filename is NULL */
    if (user_filename == (char *) NULL)
    {
	filename = "numeric.umf" ;
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

    /* write the Numeric object to the file, in binary */
    WRITE (Numeric,        NumericType, 1) ;
    WRITE (Numeric->D,     Entry, MIN (Numeric->n_row, Numeric->n_col)+1) ;
    WRITE (Numeric->Rperm, Int,   Numeric->n_row+1) ;
    WRITE (Numeric->Cperm, Int,   Numeric->n_col+1) ;
    WRITE (Numeric->Lpos,  Int,   Numeric->npiv+1) ;
    WRITE (Numeric->Lilen, Int,   Numeric->npiv+1) ;
    WRITE (Numeric->Lip,   Int,   Numeric->npiv+1) ;
    WRITE (Numeric->Upos,  Int,   Numeric->npiv+1) ;
    WRITE (Numeric->Uilen, Int,   Numeric->npiv+1) ;
    WRITE (Numeric->Uip,   Int,   Numeric->npiv+1) ;
    if (Numeric->scale != UMFPACK_SCALE_NONE)
    {
	WRITE (Numeric->Rs, double, Numeric->n_row) ;
    }
    if (Numeric->ulen > 0)
    {
	WRITE (Numeric->Upattern, Int, Numeric->ulen+1) ;
    }
    WRITE (Numeric->Memory, Unit, Numeric->size) ;

    /* close the file */
    fclose (f) ;

    return (UMFPACK_OK) ;
}
