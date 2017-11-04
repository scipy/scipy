/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/** @file superlu_enum_consts.h
 * \brief enum constants header file 
 *
 * -- SuperLU routine (version 4.1) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley, 
 * October 1, 2010
 *
 */

#ifndef __SUPERLU_ENUM_CONSTS /* allow multiple inclusions */
#define __SUPERLU_ENUM_CONSTS

/***********************************************************************
 * Enumerate types
 ***********************************************************************/
typedef enum {NO, YES}                                          yes_no_t;
typedef enum {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED} fact_t;
typedef enum {NOROWPERM, LargeDiag, MY_PERMR}                   rowperm_t;
typedef enum {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD,
	      METIS_AT_PLUS_A, PARMETIS, ZOLTAN, MY_PERMC}      colperm_t;
typedef enum {NOTRANS, TRANS, CONJ}                             trans_t;
typedef enum {NOEQUIL, ROW, COL, BOTH}                          DiagScale_t;
typedef enum {NOREFINE, SLU_SINGLE=1, SLU_DOUBLE, SLU_EXTRA}    IterRefine_t;
typedef enum {LUSUP, UCOL, LSUB, USUB, LLVL, ULVL}              MemType;
typedef enum {HEAD, TAIL}                                       stack_end_t;
typedef enum {SYSTEM, USER}                                     LU_space_t;
typedef enum {ONE_NORM, TWO_NORM, INF_NORM}			norm_t;
typedef enum {SILU, SMILU_1, SMILU_2, SMILU_3}			milu_t;
#if 0
typedef enum {NODROP		= 0x0000,
	      DROP_BASIC	= 0x0001, /* ILU(tau) */
	      DROP_PROWS	= 0x0002, /* ILUTP: keep p maximum rows */
	      DROP_COLUMN	= 0x0004, /* ILUTP: for j-th column, 
					     p = gamma * nnz(A(:,j)) */
	      DROP_AREA 	= 0x0008, /* ILUTP: for j-th column, use
					     nnz(F(:,1:j)) / nnz(A(:,1:j))
					     to limit memory growth  */
	      DROP_SECONDARY	= 0x000E, /* PROWS | COLUMN | AREA */
	      DROP_DYNAMIC	= 0x0010,
	      DROP_INTERP	= 0x0100}			rule_t;
#endif


/* 
 * The following enumerate type is used by the statistics variable 
 * to keep track of flop count and time spent at various stages.
 *
 * Note that not all of the fields are disjoint.
 */
typedef enum {
    COLPERM, /* find a column ordering that minimizes fills */
    ROWPERM, /* find a row ordering maximizes diagonal. */
    RELAX,   /* find artificial supernodes */
    ETREE,   /* compute column etree */
    EQUIL,   /* equilibrate the original matrix */
    SYMBFAC, /* symbolic factorization. */
    DIST,    /* distribute matrix. */
    FACT,    /* perform LU factorization */
    COMM,    /* communication for factorization */
    SOL_COMM,/* communication for solve */
    RCOND,   /* estimate reciprocal condition number */
    SOLVE,   /* forward and back solves */
    REFINE,  /* perform iterative refinement */
    TRSV,    /* fraction of FACT spent in xTRSV */
    GEMV,    /* fraction of FACT spent in xGEMV */
    FERR,    /* estimate error bounds after iterative refinement */
    NPHASES  /* total number of phases */
} PhaseType;


#endif /* __SUPERLU_ENUM_CONSTS */
