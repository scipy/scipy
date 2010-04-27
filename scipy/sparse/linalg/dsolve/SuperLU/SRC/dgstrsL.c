/*! @file dgstrsL.c
 * \brief Performs the L-solve using the LU factorization computed by DGSTRF
 *
 * <pre>
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * September 15, 2003
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
#include "slu_util.h"


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
 * dgstrsL only performs the L-solve using the LU factorization computed
 * by DGSTRF.
 *
 * See supermatrix.h for the definition of 'SuperMatrix' structure.
 *
 * Arguments
 * =========
 *
 * trans   (input) char*
 *          Specifies the form of the system of equations:
 *          = 'N':  A * X = B  (No transpose)
 *          = 'T':  A'* X = B  (Transpose)
 *          = 'C':  A**H * X = B  (Conjugate transpose)
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
 * perm_r  (input) int*, dimension (L->nrow)
 *         Row permutation vector, which defines the permutation matrix Pr; 
 *         perm_r[i] = j means row i of A is in position j in Pr*A.
 *
 * B       (input/output) SuperMatrix*
 *         B has types: Stype = SLU_DN, Dtype = SLU_D, Mtype = SLU_GE.
 *         On entry, the right hand side matrix.
 *         On exit, the solution matrix if info = 0;
 *
 * info    (output) int*
 * 	   = 0: successful exit
 *	   < 0: if info = -i, the i-th argument had an illegal value
 * </pre>
 */
void
dgstrsL(char *trans, SuperMatrix *L, int *perm_r, SuperMatrix *B, int *info)
{
#ifdef _CRAY
    _fcd ftcs1, ftcs2, ftcs3, ftcs4;
#endif
    int      incx = 1, incy = 1;
    double   alpha = 1.0, beta = 1.0;
    DNformat *Bstore;
    double   *Bmat;
    SCformat *Lstore;
    double   *Lval, *Uval;
    int      nrow, notran;
    int      fsupc, nsupr, nsupc, luptr, istart, irow;
    int      i, j, k, iptr, jcol, n, ldb, nrhs;
    double   *work, *work_col, *rhs_work, *soln;
    flops_t  solve_ops;
    extern SuperLUStat_t SuperLUStat;
    void dprint_soln();

    /* Test input parameters ... */
    *info = 0;
    Bstore = B->Store;
    ldb = Bstore->lda;
    nrhs = B->ncol;
    notran = lsame_(trans, "N");
    if ( !notran && !lsame_(trans, "T") && !lsame_(trans, "C") ) *info = -1;
    else if ( L->nrow != L->ncol || L->nrow < 0 ||
	      L->Stype != SLU_SC || L->Dtype != SLU_D || L->Mtype != SLU_TRLU )
	*info = -2;
    else if ( ldb < SUPERLU_MAX(0, L->nrow) ||
	      B->Stype != SLU_DN || B->Dtype != SLU_D || B->Mtype != SLU_GE )
	*info = -4;
    if ( *info ) {
	i = -(*info);
	xerbla_("dgstrsL", &i);
	return;
    }

    n = L->nrow;
    work = doubleCalloc(n * nrhs);
    if ( !work ) ABORT("Malloc fails for local work[].");
    soln = doubleMalloc(n);
    if ( !soln ) ABORT("Malloc fails for local soln[].");

    Bmat = Bstore->nzval;
    Lstore = L->Store;
    Lval = Lstore->nzval;
    solve_ops = 0;
    
    if ( notran ) {
	/* Permute right hand sides to form Pr*B */
	for (i = 0; i < nrhs; i++) {
	    rhs_work = &Bmat[i*ldb];
	    for (k = 0; k < n; k++) soln[perm_r[k]] = rhs_work[k];
	    for (k = 0; k < n; k++) rhs_work[k] = soln[k];
	}
	
	/* Forward solve PLy=Pb. */
	for (k = 0; k <= Lstore->nsuper; k++) {
	    fsupc = L_FST_SUPC(k);
	    istart = L_SUB_START(fsupc);
	    nsupr = L_SUB_START(fsupc+1) - istart;
	    nsupc = L_FST_SUPC(k+1) - fsupc;
	    nrow = nsupr - nsupc;

	    solve_ops += nsupc * (nsupc - 1) * nrhs;
	    solve_ops += 2 * nrow * nsupc * nrhs;
	    
	    if ( nsupc == 1 ) {
		for (j = 0; j < nrhs; j++) {
		    rhs_work = &Bmat[j*ldb];
	    	    luptr = L_NZ_START(fsupc);
		    for (iptr=istart+1; iptr < L_SUB_START(fsupc+1); iptr++){
			irow = L_SUB(iptr);
			++luptr;
			rhs_work[irow] -= rhs_work[fsupc] * Lval[luptr];
		    }
		}
	    } else {
	    	luptr = L_NZ_START(fsupc);
#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
		ftcs1 = _cptofcd("L", strlen("L"));
		ftcs2 = _cptofcd("N", strlen("N"));
		ftcs3 = _cptofcd("U", strlen("U"));
		STRSM( ftcs1, ftcs1, ftcs2, ftcs3, &nsupc, &nrhs, &alpha,
		       &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
		
		SGEMM( ftcs2, ftcs2, &nrow, &nrhs, &nsupc, &alpha, 
			&Lval[luptr+nsupc], &nsupr, &Bmat[fsupc], &ldb, 
			&beta, &work[0], &n );
#else
		dtrsm_("L", "L", "N", "U", &nsupc, &nrhs, &alpha,
		       &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
		
		dgemm_( "N", "N", &nrow, &nrhs, &nsupc, &alpha, 
			&Lval[luptr+nsupc], &nsupr, &Bmat[fsupc], &ldb, 
			&beta, &work[0], &n );
#endif
		for (j = 0; j < nrhs; j++) {
		    rhs_work = &Bmat[j*ldb];
		    work_col = &work[j*n];
		    iptr = istart + nsupc;
		    for (i = 0; i < nrow; i++) {
			irow = L_SUB(iptr);
			rhs_work[irow] -= work_col[i]; /* Scatter */
			work_col[i] = 0.0;
			iptr++;
		    }
		}
#else		
		for (j = 0; j < nrhs; j++) {
		    rhs_work = &Bmat[j*ldb];
		    dlsolve (nsupr, nsupc, &Lval[luptr], &rhs_work[fsupc]);
		    dmatvec (nsupr, nrow, nsupc, &Lval[luptr+nsupc],
			    &rhs_work[fsupc], &work[0] );

		    iptr = istart + nsupc;
		    for (i = 0; i < nrow; i++) {
			irow = L_SUB(iptr);
			rhs_work[irow] -= work[i];
			work[i] = 0.0;
			iptr++;
		    }
		}
#endif		    
	    } /* else ... */
	} /* for L-solve */

#ifdef DEBUG
  	printf("After L-solve: y=\n");
	dprint_soln(n, nrhs, Bmat);
#endif
	
        SuperLUStat.ops[SOLVE] = solve_ops;

    } else { 
      printf("Transposed solve not implemented.\n");
      exit(0);
    }

    SUPERLU_FREE(work);
    SUPERLU_FREE(soln);
}

/*
 * Diagnostic print of the solution vector 
 */
void
dprint_soln(int n, int nrhs, double *soln)
{
    int i;

    for (i = 0; i < n; i++) 
  	printf("\t%d: %.4f\n", i, soln[i]);
}
