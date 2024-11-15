/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/* ========================================================================== */
/* === colamd/symamd prototypes and definitions ============================= */
/* ========================================================================== */

/* COLAMD / SYMAMD include file

    You must include this file (colamd.h) in any routine that uses colamd,
    symamd, or the related macros and definitions.

    Authors:

	The authors of the code itself are Stefan I. Larimore and Timothy A.
	Davis (DrTimothyAldenDavis@gmail.com).  The algorithm was
	developed in collaboration with John Gilbert, Xerox PARC, and Esmond
	Ng, Oak Ridge National Laboratory.

    Acknowledgements:

	This work was supported by the National Science Foundation, under
	grants DMS-9504974 and DMS-9803599.

    Notice:

	Copyright (c) 1998-2007, Timothy A. Davis, All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

	Permission is hereby granted to use, copy, modify, and/or distribute
	this program, provided that the Copyright, this License, and the
	Availability of the original version is retained on all copies and made
	accessible to the end-user of any code or package that includes COLAMD
	or any modified version of COLAMD. 

    Availability:

	The colamd/symamd library is available at http://www.suitesparse.com
	This file is required by the colamd.c, colamdmex.c, and symamdmex.c
	files, and by any C code that calls the routines whose prototypes are
	listed below, or that uses the colamd/symamd definitions listed below.

*/

#ifndef COLAMD_H
#define COLAMD_H

#include "superlu_config.h"

/* make it easy for C++ programs to include COLAMD */
#ifdef __cplusplus
extern "C" {
#endif

#if defined ( _LONGINT )
#define DLONG
#endif


/* ========================================================================== */
/* === Include files ======================================================== */
/* ========================================================================== */

#include <stdlib.h>
    /*#include <stdint.h>*/

/* ========================================================================== */
/* === COLAMD version ======================================================= */
/* ========================================================================== */

/* COLAMD Version 2.4 and later will include the following definitions.
 * As an example, to test if the version you are using is 2.4 or later:
 *
 * #ifdef COLAMD_VERSION
 *	if (COLAMD_VERSION >= COLAMD_VERSION_CODE (2,4)) ...
 * #endif
 *
 * This also works during compile-time:
 *
 *  #if defined(COLAMD_VERSION) && (COLAMD_VERSION >= COLAMD_VERSION_CODE (2,4))
 *    printf ("This is version 2.4 or later\n") ;
 *  #else
 *    printf ("This is an early version\n") ;
 *  #endif
 *
 * Versions 2.3 and earlier of COLAMD do not include a #define'd version number.
 */

#define COLAMD_DATE "Oct 10, 2014"
#define COLAMD_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define COLAMD_MAIN_VERSION 2
#define COLAMD_SUB_VERSION 9
#define COLAMD_SUBSUB_VERSION 1
#define COLAMD_VERSION \
	COLAMD_VERSION_CODE(COLAMD_MAIN_VERSION,COLAMD_SUB_VERSION)

/* ========================================================================== */
/* === Knob and statistics definitions ====================================== */
/* ========================================================================== */

/* size of the knobs [ ] array.  Only knobs [0..1] are currently used. */
#define COLAMD_KNOBS 20

/* number of output statistics.  Only stats [0..6] are currently used. */
#define COLAMD_STATS 20

/* knobs [0] and stats [0]: dense row knob and output statistic. */
#define COLAMD_DENSE_ROW 0

/* knobs [1] and stats [1]: dense column knob and output statistic. */
#define COLAMD_DENSE_COL 1

/* knobs [2]: aggressive absorption */
#define COLAMD_AGGRESSIVE 2

/* stats [2]: memory defragmentation count output statistic */
#define COLAMD_DEFRAG_COUNT 2

/* stats [3]: colamd status:  zero OK, > 0 warning or notice, < 0 error */
#define COLAMD_STATUS 3

/* stats [4..6]: error info, or info on jumbled columns */ 
#define COLAMD_INFO1 4
#define COLAMD_INFO2 5
#define COLAMD_INFO3 6

/* error codes returned in stats [3]: */
#define COLAMD_OK				(0)
#define COLAMD_OK_BUT_JUMBLED			(1)
#define COLAMD_ERROR_A_not_present		(-1)
#define COLAMD_ERROR_p_not_present		(-2)
#define COLAMD_ERROR_nrow_negative		(-3)
#define COLAMD_ERROR_ncol_negative		(-4)
#define COLAMD_ERROR_nnz_negative		(-5)
#define COLAMD_ERROR_p0_nonzero			(-6)
#define COLAMD_ERROR_A_too_small		(-7)
#define COLAMD_ERROR_col_length_negative	(-8)
#define COLAMD_ERROR_row_index_out_of_bounds	(-9)
#define COLAMD_ERROR_out_of_memory		(-10)
#define COLAMD_ERROR_internal_error		(-999)


#define SUITESPARSE_PRINTF  printf


/* ========================================================================== */
/* === Prototypes of user-callable routines ================================= */
/* ========================================================================== */

/*#include "SuiteSparse_config.h" */
/**** The following is from SuiteSparse_config.h  ****/
#ifndef SuiteSparse_long

#ifdef _WIN64

#define SuiteSparse_long __int64
#define SuiteSparse_long_max _I64_MAX
#define SuiteSparse_long_idd "I64d"

#else

#if 0    /* commented out by Sherry */
#define SuiteSparse_long long
#define SuiteSparse_long_max LONG_MAX
#define SuiteSparse_long_idd "ld"
#endif

#if 1
#define SuiteSparse_long long long int
#else
#define SuiteSparse_long int64_t
#endif
#define SuiteSparse_long_max LONG_MAX
#define SuiteSparse_long_idd "lld"

#endif
#define SuiteSparse_long_id "%" SuiteSparse_long_idd
#endif
/**** end from SuiteSparse_config.h ********/


/* ========================================================================== */
/* === int or SuiteSparse_long ============================================== */
/* ========================================================================== */

#ifdef DLONG

#define Int SuiteSparse_long
#define ID  SuiteSparse_long_id
#define Int_MAX SuiteSparse_long_max

#define COLAMD_recommended colamd_l_recommended
#define COLAMD_set_defaults colamd_l_set_defaults
#define COLAMD_MAIN colamd_l
#define SYMAMD_MAIN symamd_l
#define COLAMD_report colamd_l_report
#define SYMAMD_report symamd_l_report

#else

#define Int int
#define ID "%d"
#define Int_MAX INT_MAX

#define COLAMD_recommended colamd_recommended
#define COLAMD_set_defaults colamd_set_defaults
#define COLAMD_MAIN colamd
#define SYMAMD_MAIN symamd
#define COLAMD_report colamd_report
#define SYMAMD_report symamd_report

#endif


size_t colamd_recommended	/* returns recommended value of Alen, */
				/* or 0 if input arguments are erroneous */
(
    int nnz,			/* nonzeros in A */
    int n_row,			/* number of rows in A */
    int n_col			/* number of columns in A */
) ;

size_t colamd_l_recommended	/* returns recommended value of Alen, */
				/* or 0 if input arguments are erroneous */
(
    SuiteSparse_long nnz,       /* nonzeros in A */
    SuiteSparse_long n_row,     /* number of rows in A */
    SuiteSparse_long n_col      /* number of columns in A */
) ;

void colamd_set_defaults	/* sets default parameters */
(				/* knobs argument is modified on output */
    double knobs [COLAMD_KNOBS]	/* parameter settings for colamd */
) ;

void colamd_l_set_defaults	/* sets default parameters */
(				/* knobs argument is modified on output */
    double knobs [COLAMD_KNOBS]	/* parameter settings for colamd */
) ;

int colamd			/* returns (1) if successful, (0) otherwise*/
(				/* A and p arguments are modified on output */
    int n_row,			/* number of rows in A */
    int n_col,			/* number of columns in A */
    int Alen,			/* size of the array A */
    int A [],			/* row indices of A, of size Alen */
    int p [],			/* column pointers of A, of size n_col+1 */
    double knobs [COLAMD_KNOBS],/* parameter settings for colamd */
    int stats [COLAMD_STATS]	/* colamd output statistics and error codes */
) ;

SuiteSparse_long colamd_l       /* returns (1) if successful, (0) otherwise*/
(				/* A and p arguments are modified on output */
    SuiteSparse_long n_row,     /* number of rows in A */
    SuiteSparse_long n_col,     /* number of columns in A */
    SuiteSparse_long Alen,      /* size of the array A */
    SuiteSparse_long A [],      /* row indices of A, of size Alen */
    SuiteSparse_long p [],      /* column pointers of A, of size n_col+1 */
    double knobs [COLAMD_KNOBS],/* parameter settings for colamd */
    SuiteSparse_long stats [COLAMD_STATS]   /* colamd output statistics
                                             * and error codes */
) ;

int symamd				/* return (1) if OK, (0) otherwise */
(
    int n,				/* number of rows and columns of A */
    int A [],				/* row indices of A */
    int p [],				/* column pointers of A */
    int perm [],			/* output permutation, size n_col+1 */
    double knobs [COLAMD_KNOBS],	/* parameters (uses defaults if NULL) */
    int stats [COLAMD_STATS],		/* output statistics and error codes */
    void * (*allocate) (size_t, size_t),
    					/* pointer to calloc (ANSI C) or */
					/* mxCalloc (for MATLAB mexFunction) */
    void (*release) (void *)
    					/* pointer to free (ANSI C) or */
    					/* mxFree (for MATLAB mexFunction) */
) ;

SuiteSparse_long symamd_l               /* return (1) if OK, (0) otherwise */
(
    SuiteSparse_long n,                 /* number of rows and columns of A */
    SuiteSparse_long A [],              /* row indices of A */
    SuiteSparse_long p [],              /* column pointers of A */
    SuiteSparse_long perm [],           /* output permutation, size n_col+1 */
    double knobs [COLAMD_KNOBS],	/* parameters (uses defaults if NULL) */
    SuiteSparse_long stats [COLAMD_STATS],  /* output stats and error codes */
    void * (*allocate) (size_t, size_t),
    					/* pointer to calloc (ANSI C) or */
					/* mxCalloc (for MATLAB mexFunction) */
    void (*release) (void *)
    					/* pointer to free (ANSI C) or */
    					/* mxFree (for MATLAB mexFunction) */
) ;

void colamd_report
(
    int stats [COLAMD_STATS]
) ;

void colamd_l_report
(
    SuiteSparse_long stats [COLAMD_STATS]
) ;

void symamd_report
(
    int stats [COLAMD_STATS]
) ;

void symamd_l_report
(
    SuiteSparse_long stats [COLAMD_STATS]
) ;

#ifdef __cplusplus
}
#endif

#endif /* COLAMD_H */
