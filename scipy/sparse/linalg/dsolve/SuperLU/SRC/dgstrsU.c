/*! @file dgstrsU.c
 * \brief Performs the U-solve using the LU factorization computed by DGSTRF
 *
 * <pre>
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 * 
 * Copyright (c) 1994 by Xerox Corporation.  All rights reserved.
 *
 * THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
 * EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
 *
 * Permission is hereby granted to use or copy this program for any
 * purpose, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is
 * granted, provided the above notices are retained, and a notice that
 * the code was modified is included with the above copyright notice.
 * </pre>
 */


#include "slu_ddefs.h"


/* 
 * Function prototypes 
 */
void dusolve(int, int, double*, double*);
void dlsolve(int, int, double*, double*);
void dmatvec(int, int, int, double*, double*, double*);

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * dgstrsU only performs the U-solve using the LU factorization computed
 * by DGSTRF.
 *
 * See supermatrix.h for the definition of 'SuperMatrix' structure.
 *
 * Arguments
 * =========
 *
 * trans   (input) trans_t
 *          Specifies the form of the system of equations:
 *          = NOTRANS: A * X = B  (No transpose)
 *          = TRANS:   A'* X = B  (Transpose)
 *          = CONJ:    A**H * X = B  (Conjugate transpose)
 *
 * L       (input) SuperMatrix*
 *         The factor L from the factorization Pr*A*Pc=L*U as computed by
 *         dgstrf(). Use compressed row subscripts storage for supernodes,
 *         i.e., L has types: Stype = SLU_SC, Dtype = SLU_D, Mtype = SLU_TRLU.
 *
 * U       (input) SuperMatrix*
 *         The factor U from the factorization Pr*A*Pc=L*U as computed by
 *         dgstrf(). Use column-wise storage scheme, i.e., U has types:
 *         Stype = SLU_NC, Dtype = SLU_D, Mtype = SLU_TRU.
 *
 * perm_c  (input) int*, dimension (L->ncol)
 *	   Column permutation vector, which defines the 
 *         permutation matrix Pc; perm_c[i] = j means column i of A is 
 *         in position j in A*Pc.
 *
 * perm_r  (input) int*, dimension (L->nrow)
 *         Row permutation vector, which defines the permutation matrix Pr; 
 *         perm_r[i] = j means row i of A is in position j in Pr*A.
 *
 * B       (input/output) SuperMatrix*
 *         B has types: Stype = SLU_DN, Dtype = SLU_D, Mtype = SLU_GE.
 *         On entry, the right hand side matrix.
 *         On exit, the solution matrix if info = 0;
 *
 * stat     (output) SuperLUStat_t*
 *          Record the statistics on runtime and floating-point operation count.
 *          See util.h for the definition of 'SuperLUStat_t'.
 *
 * info    (output) int*
 * 	   = 0: successful exit
 *	   < 0: if info = -i, the i-th argument had an illegal value
 * </pre>
 */
void
dgstrsU(trans_t trans, SuperMatrix *L, SuperMatrix *U,
        int *perm_c, int *perm_r, SuperMatrix *B,
        SuperLUStat_t *stat, int *info)
{
#ifdef _CRAY
    _fcd ftcs1, ftcs2, ftcs3, ftcs4;
#endif
    int      incx = 1, incy = 1;
#ifdef USE_VENDOR_BLAS
    double   alpha = 1.0, beta = 1.0;
    double   *work_col;
#endif
    DNformat *Bstore;
    double   *Bmat;
    SCformat *Lstore;
    NCformat *Ustore;
    double   *Lval, *Uval;
    int      fsupc, nrow, nsupr, nsupc, luptr, istart, irow;
    int      i, j, k, iptr, jcol, n, ldb, nrhs;
    double   *rhs_work, *soln;
    flops_t  solve_ops;
    void dprint_soln();

    /* Test input parameters ... */
    *info = 0;
    Bstore = B->Store;
    ldb = Bstore->lda;
    nrhs = B->ncol;
    if ( trans != NOTRANS && trans != TRANS && trans != CONJ ) *info = -1;
    else if ( L->nrow != L->ncol || L->nrow < 0 ||
	      L->Stype != SLU_SC || L->Dtype != SLU_D || L->Mtype != SLU_TRLU )
	*info = -2;
    else if ( U->nrow != U->ncol || U->nrow < 0 ||
	      U->Stype != SLU_NC || U->Dtype != SLU_D || U->Mtype != SLU_TRU )
	*info = -3;
    else if ( ldb < SUPERLU_MAX(0, L->nrow) ||
	      B->Stype != SLU_DN || B->Dtype != SLU_D || B->Mtype != SLU_GE )
	*info = -6;
    if ( *info ) {
	i = -(*info);
	xerbla_("dgstrs", &i);
	return;
    }

    n = L->nrow;
    soln = doubleMalloc(n);
    if ( !soln ) ABORT("Malloc fails for local soln[].");

    Bmat = Bstore->nzval;
    Lstore = L->Store;
    Lval = Lstore->nzval;
    Ustore = U->Store;
    Uval = Ustore->nzval;
    solve_ops = 0;
    
    if ( trans == NOTRANS ) {
	/*
	 * Back solve Ux=y.
	 */
	for (k = Lstore->nsuper; k >= 0; k--) {
	    fsupc = L_FST_SUPC(k);
	    istart = L_SUB_START(fsupc);
	    nsupr = L_SUB_START(fsupc+1) - istart;
	    nsupc = L_FST_SUPC(k+1) - fsupc;
	    luptr = L_NZ_START(fsupc);

	    solve_ops += nsupc * (nsupc + 1) * nrhs;

	    if ( nsupc == 1 ) {
		rhs_work = &Bmat[0];
		for (j = 0; j < nrhs; j++) {
		    rhs_work[fsupc] /= Lval[luptr];
		    rhs_work += ldb;
		}
	    } else {
#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
		ftcs1 = _cptofcd("L", strlen("L"));
		ftcs2 = _cptofcd("U", strlen("U"));
		ftcs3 = _cptofcd("N", strlen("N"));
		STRSM( ftcs1, ftcs2, ftcs3, ftcs3, &nsupc, &nrhs, &alpha,
		       &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
#else
		dtrsm_("L", "U", "N", "N", &nsupc, &nrhs, &alpha,
		       &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
#endif
#else		
		for (j = 0; j < nrhs; j++)
		    dusolve ( nsupr, nsupc, &Lval[luptr], &Bmat[fsupc+j*ldb] );
#endif		
	    }

	    for (j = 0; j < nrhs; ++j) {
		rhs_work = &Bmat[j*ldb];
		for (jcol = fsupc; jcol < fsupc + nsupc; jcol++) {
		    solve_ops += 2*(U_NZ_START(jcol+1) - U_NZ_START(jcol));
		    for (i = U_NZ_START(jcol); i < U_NZ_START(jcol+1); i++ ){
			irow = U_SUB(i);
			rhs_work[irow] -= rhs_work[jcol] * Uval[i];
		    }
		}
	    }
	    
	} /* for U-solve */

#ifdef DEBUG
  	printf("After U-solve: x=\n");
	dprint_soln(n, nrhs, Bmat);
#endif

	/* Compute the final solution X := Pc*X. */
	for (i = 0; i < nrhs; i++) {
	    rhs_work = &Bmat[i*ldb];
	    for (k = 0; k < n; k++) soln[k] = rhs_work[perm_c[k]];
	    for (k = 0; k < n; k++) rhs_work[k] = soln[k];
	}
	
        stat->ops[SOLVE] = solve_ops;

    } else { /* Solve U'x = b */
	/* Permute right hand sides to form Pc'*B. */
	for (i = 0; i < nrhs; i++) {
	    rhs_work = &Bmat[i*ldb];
	    for (k = 0; k < n; k++) soln[perm_c[k]] = rhs_work[k];
	    for (k = 0; k < n; k++) rhs_work[k] = soln[k];
	}

	for (k = 0; k < nrhs; ++k) {
	    /* Multiply by inv(U'). */
	    sp_dtrsv("U", "T", "N", L, U, &Bmat[k*ldb], stat, info);
	}

    }

    SUPERLU_FREE(soln);
}

