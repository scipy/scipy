/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/** @file slu_util.h
 * \brief Utility header file 
 *
 * -- SuperLU routine (version 4.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
 *
 */

#ifndef __SUPERLU_UTIL /* allow multiple inclusions */
#define __SUPERLU_UTIL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
#ifndef __STDC__
#include <malloc.h>
#endif
*/
#include <assert.h>
#include "superlu_enum_consts.h"

#include "scipy_slu_config.h"

/***********************************************************************
 * Macros
 ***********************************************************************/
/*                                                                                           
 * You can support older version of SuperLU.                                              
 * At compile-time, you can catch the new release as:                                          
 *   #ifdef SUPERLU_MAJOR_VERSION == 5
 *       use the new interface                                                                 
 *   #else                                                                                     
 *       use the old interface                                                                 
 *   #endif                                                                                    
 * Versions 4.x and earlier do not include a #define'd version numbers.                        
 */
#define SUPERLU_MAJOR_VERSION     6
#define SUPERLU_MINOR_VERSION     0
#define SUPERLU_PATCH_VERSION     1


#define FIRSTCOL_OF_SNODE(i)	(xsup[i])
/* No of marker arrays used in the symbolic factorization,
   each of size n */
#define NO_MARKER     3
#define NUM_TEMPV(m,w,t,b)  ( SUPERLU_MAX(m, (t + b)*w) )

#ifndef USER_ABORT
#define USER_ABORT(msg) superlu_abort_and_exit(msg)
#endif

#define ABORT(err_msg) \
 { char msg[256];\
   sprintf(msg,"%s at line %d in file %s\n",err_msg,__LINE__, __FILE__);\
   USER_ABORT(msg); }


#ifndef USER_MALLOC
#if 1
#define USER_MALLOC(size) superlu_malloc(size)
#else
/* The following may check out some uninitialized data */
#define USER_MALLOC(size) memset (superlu_malloc(size), '\x0F', size)
#endif
#endif

#define SUPERLU_MALLOC(size) USER_MALLOC(size)

#ifndef USER_FREE
#define USER_FREE(addr) superlu_free(addr)
#endif

#define SUPERLU_FREE(addr) USER_FREE(addr)

#define CHECK_MALLOC(where) {                 \
    extern int superlu_malloc_total;        \
    printf("%s: malloc_total %d Bytes\n",     \
	   where, superlu_malloc_total); \
}

#define SUPERLU_MAX(x, y) 	( (x) > (y) ? (x) : (y) )
#define SUPERLU_MIN(x, y) 	( (x) < (y) ? (x) : (y) )

/*********************************************************
 * Macros used for easy access of sparse matrix entries. *
 *********************************************************/
#define L_SUB_START(col)     ( Lstore->rowind_colptr[col] )
#define L_SUB(ptr)           ( Lstore->rowind[ptr] )
#define L_NZ_START(col)      ( Lstore->nzval_colptr[col] )
#define L_FST_SUPC(superno)  ( Lstore->sup_to_col[superno] )
#define U_NZ_START(col)      ( Ustore->colptr[col] )
#define U_SUB(ptr)           ( Ustore->rowind[ptr] )


/***********************************************************************
 * Constants 
 ***********************************************************************/
#define EMPTY	(-1)
/*#define NO	(-1)*/
#define FALSE	0
#define TRUE	1

#if 0 // this was old; new one is 6, defined in superlu_enum_consts.h 
#define NO_MEMTYPE  4      /* 0: lusup;
			      1: ucol;
			      2: lsub;
			      3: usub */
#endif

#define GluIntArray(n)   (5 * (n) + 5)

/* Dropping rules */
#define  NODROP	        ( 0x0000 )
#define	 DROP_BASIC	( 0x0001 )  /* ILU(tau) */
#define  DROP_PROWS	( 0x0002 )  /* ILUTP: keep p maximum rows */
#define  DROP_COLUMN	( 0x0004 )  /* ILUTP: for j-th column, 
				              p = gamma * nnz(A(:,j)) */
#define  DROP_AREA 	( 0x0008 )  /* ILUTP: for j-th column, use
 		 			      nnz(F(:,1:j)) / nnz(A(:,1:j))
					      to limit memory growth  */
#define  DROP_SECONDARY	( 0x000E )  /* PROWS | COLUMN | AREA */
#define  DROP_DYNAMIC	( 0x0010 )  /* adaptive tau */
#define  DROP_INTERP	( 0x0100 )  /* use interpolation */


#if 1
#define MILU_ALPHA (1.0e-2) /* multiple of drop_sum to be added to diagonal */
#else
#define MILU_ALPHA  1.0 /* multiple of drop_sum to be added to diagonal */
#endif


/***********************************************************************
 * Type definitions
 ***********************************************************************/
typedef float    flops_t;
typedef unsigned char Logical;

/* 
 *-- This contains the options used to control the solution process.
 *
 * Fact   (fact_t)
 *        Specifies whether or not the factored form of the matrix
 *        A is supplied on entry, and if not, how the matrix A should
 *        be factorizaed.
 *        = DOFACT: The matrix A will be factorized from scratch, and the
 *             factors will be stored in L and U.
 *        = SamePattern: The matrix A will be factorized assuming
 *             that a factorization of a matrix with the same sparsity
 *             pattern was performed prior to this one. Therefore, this
 *             factorization will reuse column permutation vector 
 *             ScalePermstruct->perm_c and the column elimination tree
 *             LUstruct->etree.
 *        = SamePattern_SameRowPerm: The matrix A will be factorized
 *             assuming that a factorization of a matrix with the same
 *             sparsity	pattern and similar numerical values was performed
 *             prior to this one. Therefore, this factorization will reuse
 *             both row and column scaling factors R and C, both row and
 *             column permutation vectors perm_r and perm_c, and the
 *             L & U data structures set up from the previous factorization.
 *        = FACTORED: On entry, L, U, perm_r and perm_c contain the 
 *              factored form of A. If DiagScale is not NOEQUIL, the matrix
 *              A has been equilibrated with scaling factors R and C.
 *
 * Equil  (yes_no_t)
 *        Specifies whether to equilibrate the system (scale A's row and
 *        columns to have unit norm).
 *
 * ColPerm (colperm_t)
 *        Specifies what type of column permutation to use to reduce fill.
 *        = NATURAL: use the natural ordering 
 *        = MMD_ATA: use minimum degree ordering on structure of A'*A
 *        = MMD_AT_PLUS_A: use minimum degree ordering on structure of A'+A
 *        = COLAMD: use approximate minimum degree column ordering
 *        = MY_PERMC: use the ordering specified by the user
 *         
 * Trans  (trans_t)
 *        Specifies the form of the system of equations:
 *        = NOTRANS: A * X = B        (No transpose)
 *        = TRANS:   A**T * X = B     (Transpose)
 *        = CONJ:    A**H * X = B     (Transpose)
 *
 * IterRefine (IterRefine_t)
 *        Specifies whether to perform iterative refinement.
 *        = NO: no iterative refinement
 *        = SLU_SINGLE: perform iterative refinement in single precision
 *        = SLU_DOUBLE: perform iterative refinement in double precision
 *        = SLU_EXTRA: perform iterative refinement in extra precision
 *
 * DiagPivotThresh (double, in [0.0, 1.0]) (only for sequential SuperLU)
 *        Specifies the threshold used for a diagonal entry to be an
 *        acceptable pivot.
 *
 * SymmetricMode (yest_no_t)
 *        Specifies whether to use symmetric mode. Symmetric mode gives 
 *        preference to diagonal pivots, and uses an (A'+A)-based column
 *        permutation algorithm.
 *
 * PivotGrowth (yes_no_t)
 *        Specifies whether to compute the reciprocal pivot growth.
 *
 * ConditionNumber (ues_no_t)
 *        Specifies whether to compute the reciprocal condition number.
 *
 * RowPerm (rowperm_t) (only for SuperLU_DIST or ILU)
 *        Specifies whether to permute rows of the original matrix.
 *        = NO: not to permute the rows
 *        = LargeDiag: make the diagonal large relative to the off-diagonal
 *        = MY_PERMR: use the permutation given by the user
 *
 * ILU_DropRule (int)
 *        Specifies the dropping rule:
 *	  = DROP_BASIC:   Basic dropping rule, supernodal based ILUTP(tau).
 *	  = DROP_PROWS:   Supernodal based ILUTP(p,tau), p = gamma * nnz(A)/n.
 *	  = DROP_COLUMN:  Variant of ILUTP(p,tau), for j-th column,
 *			      p = gamma * nnz(A(:,j)).
 *	  = DROP_AREA:    Variation of ILUTP, for j-th column, use
 *			      nnz(F(:,1:j)) / nnz(A(:,1:j)) to control memory.
 *	  = DROP_DYNAMIC: Modify the threshold tau during factorizaion:
 *			  If nnz(L(:,1:j)) / nnz(A(:,1:j)) > gamma
 *				  tau_L(j) := MIN(tau_0, tau_L(j-1) * 2);
 *			  Otherwise
 *				  tau_L(j) := MAX(tau_0, tau_L(j-1) / 2);
 *			  tau_U(j) uses the similar rule.
 *			  NOTE: the thresholds used by L and U are separate.
 *	  = DROP_INTERP:  Compute the second dropping threshold by
 *	                  interpolation instead of sorting (default).
 *  		          In this case, the actual fill ratio is not
 *			  guaranteed to be smaller than gamma.
 *   	  Note: DROP_PROWS, DROP_COLUMN and DROP_AREA are mutually exclusive.
 *	  ( Default: DROP_BASIC | DROP_AREA )
 *
 * ILU_DropTol (double)
 *        numerical threshold for dropping.
 *
 * ILU_FillFactor (double) 
 *        Gamma in the secondary dropping.
 *
 * ILU_Norm (norm_t)
 *        Specify which norm to use to measure the row size in a
 *        supernode: infinity-norm, 1-norm, or 2-norm.
 *
 * ILU_FillTol (double)
 *        numerical threshold for zero pivot perturbation.
 *
 * ILU_MILU (milu_t)
 *        Specifies which version of MILU to use.
 *
 * ILU_MILU_Dim (double) 
 *        Dimension of the PDE if available.
 *
 * ReplaceTinyPivot (yes_no_t) (only for SuperLU_DIST)
 *        Specifies whether to replace the tiny diagonals by
 *        sqrt(epsilon)*||A|| during LU factorization.
 *
 * SolveInitialized (yes_no_t) (only for SuperLU_DIST)
 *        Specifies whether the initialization has been performed to the
 *        triangular solve.
 *
 * RefineInitialized (yes_no_t) (only for SuperLU_DIST)
 *        Specifies whether the initialization has been performed to the
 *        sparse matrix-vector multiplication routine needed in iterative
 *        refinement.
 *
 * PrintStat (yes_no_t)
 *        Specifies whether to print the solver's statistics.
 */
typedef struct {
    fact_t        Fact;
    yes_no_t      Equil;
    colperm_t     ColPerm;
    trans_t       Trans;
    IterRefine_t  IterRefine;
    double        DiagPivotThresh;
    yes_no_t      SymmetricMode;
    yes_no_t      PivotGrowth;
    yes_no_t      ConditionNumber;
    rowperm_t     RowPerm;
    int 	  ILU_DropRule;
    double	  ILU_DropTol;    /* threshold for dropping */
    double	  ILU_FillFactor; /* gamma in the secondary dropping */
    norm_t	  ILU_Norm;       /* infinity-norm, 1-norm, or 2-norm */
    double	  ILU_FillTol;    /* threshold for zero pivot perturbation */
    milu_t	  ILU_MILU;
    double	  ILU_MILU_Dim;   /* Dimension of PDE (if available) */
    yes_no_t      ParSymbFact;
    yes_no_t      ReplaceTinyPivot; /* used in SuperLU_DIST */
    yes_no_t      SolveInitialized;
    yes_no_t      RefineInitialized;
    yes_no_t      PrintStat;
    int           nnzL, nnzU;      /* used to store nnzs for now       */
    int           num_lookaheads;  /* num of levels in look-ahead      */
    yes_no_t      lookahead_etree; /* use etree computed from the
				      serial symbolic factorization */
    yes_no_t      SymPattern;      /* symmetric factorization          */
} superlu_options_t;

/*! \brief Headers for 4 types of dynamatically managed memory */
typedef struct e_node {
    int_t size;    /* length of the memory mem[] that has been used */
    void *mem;     /* pointer to the new malloc'd store */
} ExpHeader;

typedef struct {
    int_t  size;
    int_t  used;
    int_t  top1;  /* grow upward, relative to &array[0] */
    int_t  top2;  /* grow downward */
    void *array;
} LU_stack_t;

typedef struct {
    int     *panel_histo; /* histogram of panel size distribution */
    double  *utime;       /* running time at various phases */
    flops_t *ops;         /* operation count at various phases */
    int     TinyPivots;   /* number of tiny pivots */
    int     RefineSteps;  /* number of iterative refinement steps */
    int     expansions;   /* number of memory expansions */
} SuperLUStat_t;

typedef struct {
    float for_lu;
    float total_needed;
} mem_usage_t;


typedef struct {
    int     *xsup;    /* supernode and column mapping */
    int     *supno;   
    int_t   *lsub;    /* compressed L subscripts */
    int_t   *xlsub;
    void    *lusup;   /* L supernodes */
    int_t   *xlusup;
    void    *ucol;    /* U columns */
    int_t   *usub;
    int_t   *xusub;
    int_t   nzlmax;   /* current max size of lsub */
    int_t   nzumax;   /*    "    "    "      ucol */
    int_t   nzlumax;  /*    "    "    "     lusup */
    int     n;        /* number of columns in the matrix */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
    int     num_expansions;
    ExpHeader *expanders; /* Array of pointers to 4 types of memory */
    LU_stack_t stack;     /* use user supplied memory */
} GlobalLU_t;


/***********************************************************************
 * Prototypes
 ***********************************************************************/
#ifdef __cplusplus
extern "C" {
#endif

extern int     input_error(char *, int *);

extern void    Destroy_SuperMatrix_Store(SuperMatrix *);
extern void    Destroy_CompCol_Matrix(SuperMatrix *);
extern void    Destroy_CompRow_Matrix(SuperMatrix *);
extern void    Destroy_SuperNode_Matrix(SuperMatrix *);
extern void    Destroy_CompCol_Permuted(SuperMatrix *);
extern void    Destroy_Dense_Matrix(SuperMatrix *);
extern void    get_perm_c(int, SuperMatrix *, int *);
extern void    set_default_options(superlu_options_t *options);
extern void    ilu_set_default_options(superlu_options_t *options);
extern void    sp_preorder (superlu_options_t *, SuperMatrix*, int*, int*,
			    SuperMatrix*);
extern void    superlu_abort_and_exit(char*);
extern void    *superlu_malloc (size_t);
extern int     *int32Malloc (int);
extern int     *int32Calloc (int);
extern int_t   *intMalloc (int_t);
extern int_t   *intCalloc (int_t);
extern void    superlu_free (void*);
extern void    SetIWork (int, int, int, int *, int **, int **, int_t **xplore,
                         int **, int **, int_t **xprune, int **);
extern int     sp_coletree (int_t *, int_t *, int_t *, int, int, int *);
extern void    relax_snode (const int, int *, const int, int *, int *);
extern void    heap_relax_snode (const int, int *, const int, int *, int *);
extern int     mark_relax(int, int *, int *, int_t *, int_t *, int_t *, int *);
extern void    countnz(const int n, int_t *xprune, int_t *nnzL, int_t *nnzU, GlobalLU_t *);
extern void    ilu_countnz (const int, int_t *, int_t *, GlobalLU_t *);
extern void    fixupL (const int, const int *, GlobalLU_t *);
extern void    ilu_relax_snode (const int, int *, const int, int *,
				int *, int *);
extern void    ilu_heap_relax_snode (const int, int *, const int, int *,
				     int *, int*);
extern void    resetrep_col (const int, const int *, int *);
extern int     spcoletree (int *, int *, int *, int, int, int *);
extern int     *TreePostorder (int, int *);
extern double  SuperLU_timer_ (void);
extern int     sp_ienv (int);
extern int     xerbla_ (char *, int *);
extern void    ifill (int *, int, int);
extern void    snode_profile (int, int *);
extern void    super_stats (int, int *);
extern void    check_repfnz(int, int, int, int *);
extern void    PrintSumm (char *, int, int, int);
extern void    StatInit(SuperLUStat_t *);
extern void    StatPrint (SuperLUStat_t *);
extern void    StatFree(SuperLUStat_t *);
extern void    print_panel_seg(int, int, int, int, int *, int *);
extern void    print_int_vec(char *what, int n, int *vec);
extern void    slu_PrintInt10(char *name, int len, int *x);
extern void    check_perm(char *what, int n, int *perm);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_UTIL */
