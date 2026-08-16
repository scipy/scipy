/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file slu_sdefs.h
 * \brief Header file for real operations
 * 
 * <pre> 
 * -- SuperLU routine (version 7.0.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
 * August 2024
 * 
 * Global data structures used in LU factorization -
 * 
 *   nsuper: \#supernodes = nsuper + 1, numbered [0, nsuper].
 *   (xsup,supno): supno[i] is the supernode no to which i belongs;
 *	xsup(s) points to the beginning of the s-th supernode.
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new 
 *	      row pivoting ordering
 *
 *   (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub 
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript interchange
 *	performed, similar interchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numerical LU factorization phase.
 *
 *   (xlusup,lusup): lusup[*] contains the numerical values of the
 *	rectangular supernodes; xlusup[j] points to the starting
 *	location of the j-th column in storage vector lusup[*]
 *	Note: xlusup is indexed by column.
 *	Each rectangular supernode is stored by column-major
 *	scheme, consistent with Fortran 2-dim array storage.
 *
 *   (xusub,ucol,usub): ucol[*] stores the numerical values of
 *	U-columns outside the rectangular supernodes. The row
 *	subscript of nonzero ucol[k] is stored in usub[k].
 *	xusub[i] points to the starting location of column i in ucol.
 *	Storage: new row subscripts; that is subscripts of PA.
 * </pre>
 */
#ifndef __SUPERLU_sSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_sSP_DEFS

/*
 * File name:		ssp_defs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */

#ifdef _CRAY
#include <fortran.h>
#endif

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "slu_Cnames.h"
#include "superlu_config.h"
#include "supermatrix.h"
#include "slu_util.h"


/* -------- Prototypes -------- */

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Driver routines */
extern void
sgssv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int_t *info);
extern void
sgssvx(superlu_options_t *, SuperMatrix *, int *, int *, int *,
       char *, float *, float *, SuperMatrix *, SuperMatrix *,
       void *, int_t lwork, SuperMatrix *, SuperMatrix *,
       float *, float *, float *, float *,
       GlobalLU_t *, mem_usage_t *, SuperLUStat_t *, int_t *info);
    /* ILU */
extern void
sgsisv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
extern void
sgsisx(superlu_options_t *options, SuperMatrix *A, int *perm_c, int *perm_r,
       int *etree, char *equed, float *R, float *C,
       SuperMatrix *L, SuperMatrix *U, void *work, int_t lwork,
       SuperMatrix *B, SuperMatrix *X, float *recip_pivot_growth, float *rcond,
       GlobalLU_t *Glu, mem_usage_t *mem_usage, SuperLUStat_t *stat, int_t *info);


/*! \brief Supernodal LU factor related */
extern void
sCreate_CompCol_Matrix(SuperMatrix *, int, int, int_t, float *,
		       int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
extern void
sCreate_CompRow_Matrix(SuperMatrix *, int, int, int_t, float *,
		       int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
extern void sCompRow_to_CompCol(int, int, int_t, float*, int_t*, int_t*,
		                   float **, int_t **, int_t **);
extern void
sCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
sCreate_Dense_Matrix(SuperMatrix *, int, int, float *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
sCreate_SuperNode_Matrix(SuperMatrix *, int, int, int_t, float *, 
		         int_t *, int_t *, int_t *, int *, int *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
sCopy_Dense_Matrix(int, int, float *, int, float *, int);

extern void    sallocateA (int, int_t, float **, int_t **, int_t **);
extern void    sgstrf (superlu_options_t*, SuperMatrix*,
                       int, int, int*, void *, int_t, int *, int *, 
                       SuperMatrix *, SuperMatrix *, GlobalLU_t *,
		       SuperLUStat_t*, int_t *info);
extern int_t   ssnode_dfs (const int, const int, const int_t *, const int_t *,
			     const int_t *, int_t *, int *, GlobalLU_t *);
extern int     ssnode_bmod (const int, const int, const int, float *,
                              float *, GlobalLU_t *, SuperLUStat_t*);
extern void    spanel_dfs (const int, const int, const int, SuperMatrix *,
			   int *, int *, float *, int *, int *, int *,
			   int_t *, int *, int *, int_t *, GlobalLU_t *);
extern void    spanel_bmod (const int, const int, const int, const int,
                           float *, float *, int *, int *,
			   GlobalLU_t *, SuperLUStat_t*);
extern int     scolumn_dfs (const int, const int, int *, int *, int *, int *,
			   int *, int_t *, int *, int *, int_t *, GlobalLU_t *);
extern int     scolumn_bmod (const int, const int, float *,
			   float *, int *, int *, int,
                           GlobalLU_t *, SuperLUStat_t*);
extern int     scopy_to_ucol (int, int, int *, int *, int *,
                              float *, GlobalLU_t *);         
extern int     spivotL (const int, const double, int *, int *, 
                         int *, int *, int *, GlobalLU_t *, SuperLUStat_t*);
extern void    spruneL (const int, const int *, const int, const int,
			  const int *, const int *, int_t *, GlobalLU_t *);
extern void    sreadmt (int *, int *, int_t *, float **, int_t **, int_t **);
extern void    sGenXtrue (int, int, float *, int);
extern void    sFillRHS (trans_t, int, float *, int, SuperMatrix *,
			  SuperMatrix *);
extern void    sgstrs (trans_t, SuperMatrix *, SuperMatrix *, const int *, const int *,
                        SuperMatrix *, SuperLUStat_t*, int *);
/* ILU */
extern void    sgsitrf (superlu_options_t*, SuperMatrix*, int, int, int*,
		        void *, int_t, int *, int *, SuperMatrix *, SuperMatrix *,
                        GlobalLU_t *, SuperLUStat_t*, int_t *info);
extern int     sldperm(int, int, int_t, int_t [], int_t [], float [],
                        int [],	float [], float []);
extern int     ilu_ssnode_dfs (const int, const int, const int_t *, const int_t *,
			       const int_t *, int *, GlobalLU_t *);
extern void    ilu_spanel_dfs (const int, const int, const int, SuperMatrix *,
			       int *, int *, float *, float *, int *, int *,
			       int *, int *, int *, int_t *, GlobalLU_t *);
extern int     ilu_scolumn_dfs (const int, const int, int *, int *, int *,
				int *, int *, int *, int *, int_t *, GlobalLU_t *);
extern int     ilu_scopy_to_ucol (int, int, int *, int *, int *,
                                  float *, int, milu_t, double, int,
                                  float *, int *, GlobalLU_t *, float *);
extern int     ilu_spivotL (const int, const double, int *, int *, int, int *,
			    int *, int *, int *, double, milu_t,
                            float, GlobalLU_t *, SuperLUStat_t*);
extern int     ilu_sdrop_row (superlu_options_t *, int, int, double,
                              int, int *, double *, GlobalLU_t *, 
                              float *, float *, int);


/*! \brief Driver related */

extern void    sgsequ (SuperMatrix *, float *, float *, float *,
			float *, float *, int *);
extern void    slaqgs (SuperMatrix *, float *, float *, float,
                        float, float, char *);
extern void    sgscon (char *, SuperMatrix *, SuperMatrix *, 
		         float, float *, SuperLUStat_t*, int *);
extern float   sPivotGrowth(int, SuperMatrix *, int *, 
                            SuperMatrix *, SuperMatrix *);
extern void    sgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, int *, int *, char *, float *, 
                       float *, SuperMatrix *, SuperMatrix *,
                       float *, float *, SuperLUStat_t*, int *);

extern int     sp_strsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, float *, SuperLUStat_t*, int *);
extern int     sp_sgemv (char *, float, SuperMatrix *, float *,
			int, float, float *, int);

extern int     sp_sgemm (char *, char *, int, int, int, float,
			SuperMatrix *, float *, int, float, 
			float *, int);
extern         float smach(char *);   /* from C99 standard, in float.h */

/*! \brief Memory-related */
extern int_t   sLUMemInit (fact_t, void *, int_t, int, int, int_t, int,
                            float, SuperMatrix *, SuperMatrix *,
                            GlobalLU_t *, int **, float **);
extern void    sSetRWork (int, int, float *, float **, float **);
extern void    sLUWorkFree (int *, float *, GlobalLU_t *);
extern int_t   sLUMemXpand (int, int_t, MemType, int_t *, GlobalLU_t *);

extern float  *floatMalloc(size_t);
extern float  *floatCalloc(size_t);
extern int_t   smemory_usage(const int_t, const int_t, const int_t, const int);
extern int     sQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
extern int     ilu_sQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

/*! \brief Auxiliary routines */
extern void    sreadhb(FILE *, int *, int *, int_t *, float **, int_t **, int_t **);
extern void    sreadrb(int *, int *, int_t *, float **, int_t **, int_t **);
extern void    sreadtriple(int *, int *, int_t *, float **, int_t **, int_t **);
extern void    sreadtriple_noheader(int *, int *, int_t *, float **, int_t **, int_t **);
extern void    sreadMM(FILE *, int *, int *, int_t *, float **, int_t **, int_t **);
extern void    sfill (float *, int, float);
extern void    sinf_norm_error (int, SuperMatrix *, float *);
extern float  sqselect(int, float *, int);


/*! \brief Routines for debugging */
extern void    sPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    sPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    sPrint_Dense_Matrix(char *, SuperMatrix *);
extern void    sprint_lu_col(char *, int, int, int_t *, GlobalLU_t *);
extern int     print_float_vec(const char *, int, const float *);
extern void    scheck_tempv(int, float *);

/*! \brief BLAS */

extern void scopy_(int *, float *, int *, float *, int *);
extern void saxpy_(int *, float *, float *, int *, float *, int *);
extern void sgemm_(const char*, const char*, const int*, const int*, const int*,
                  const float*, const float*, const int*, const float*,
		  const int*, const float*, float*, const int*);
extern void strsv_(char*, char*, char*, int*, float*, int*,
                  float*, int*);
extern void strsm_(char*, char*, char*, char*, int*, int*,
                  float*, float*, int*, float*, int*);
extern void sgemv_(char *, int *, int *, float *, float *a, int *,
                  float *, int *, float *, float *, int *);

extern void susolve(int, int, float*, float*);
extern void slsolve(int, int, float*, float*);
extern void smatvec(int, int, int, float*, float*, float*);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_sSP_DEFS */

