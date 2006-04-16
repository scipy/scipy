

/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#include "dsp_defs.h"
#include "util.h"

void
dgssv(SuperMatrix *A, int *perm_c, int *perm_r, SuperMatrix *L,
      SuperMatrix *U, SuperMatrix *B, int *info )
{
/*
 * Purpose
 * =======
 *
 * DGSSV solves the system of linear equations A*X=B, using the
 * LU factorization from DGSTRF. It performs the following steps:
 *
 *   1. If A is stored column-wise (A->Stype = NC):
 *
 *      1.1. Permute the columns of A, forming A*Pc, where Pc
 *           is a permutation matrix. For more details of this step, 
 *           see sp_preorder.c.
 *
 *      1.2. Factor A as Pr*A*Pc=L*U with the permutation Pr determined
 *           by Gaussian elimination with partial pivoting.
 *           L is unit lower triangular with offdiagonal entries
 *           bounded by 1 in magnitude, and U is upper triangular.
 *
 *      1.3. Solve the system of equations A*X=B using the factored
 *           form of A.
 *
 *   2. If A is stored row-wise (A->Stype = NR), apply the
 *      above algorithm to the transpose of A:
 *
 *      2.1. Permute columns of transpose(A) (rows of A),
 *           forming transpose(A)*Pc, where Pc is a permutation matrix. 
 *           For more details of this step, see sp_preorder.c.
 *
 *      2.2. Factor A as Pr*transpose(A)*Pc=L*U with the permutation Pr
 *           determined by Gaussian elimination with partial pivoting.
 *           L is unit lower triangular with offdiagonal entries
 *           bounded by 1 in magnitude, and U is upper triangular.
 *
 *      2.3. Solve the system of equations A*X=B using the factored
 *           form of A.
 *
 *   See supermatrix.h for the definition of 'SuperMatrix' structure.
 * 
 * Arguments
 * =========
 *
 * A       (input) SuperMatrix*
 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *         of linear equations is A->nrow. Currently, the type of A can be:
 *         Stype = NC or NR; Dtype = D_; Mtype = GE. In the future, more
 *         general A will be handled.
 *
 * perm_c  (input/output) int*
 *         If A->Stype = NC, column permutation vector of size A->ncol
 *         which defines the permutation matrix Pc; perm_c[i] = j means 
 *         column i of A is in position j in A*Pc.
 *         On exit, perm_c may be overwritten by the product of the input
 *         perm_c and a permutation that postorders the elimination tree
 *         of Pc'*A'*A*Pc; perm_c is not changed if the elimination tree
 *         is already in postorder.
 *
 *         If A->Stype = NR, column permutation vector of size A->nrow
 *         which describes permutation of columns of transpose(A) 
 *         (rows of A) as described above.
 * 
 * perm_r  (output) int*
 *         If A->Stype = NC, row permutation vector of size A->nrow, 
 *         which defines the permutation matrix Pr, and is determined 
 *         by partial pivoting.  perm_r[i] = j means row i of A is in 
 *         position j in Pr*A.
 *
 *         If A->Stype = NR, permutation vector of size A->ncol, which
 *         determines permutation of rows of transpose(A)
 *         (columns of A) as described above.
 *
 * L       (output) SuperMatrix*
 *         The factor L from the factorization 
 *             Pr*A*Pc=L*U              (if A->Stype = NC) or
 *             Pr*transpose(A)*Pc=L*U   (if A->Stype = NR).
 *         Uses compressed row subscripts storage for supernodes, i.e.,
 *         L has types: Stype = SC, Dtype = D_, Mtype = TRLU.
 *         
 * U       (output) SuperMatrix*
 *	   The factor U from the factorization 
 *             Pr*A*Pc=L*U              (if A->Stype = NC) or
 *             Pr*transpose(A)*Pc=L*U   (if A->Stype = NR).
 *         Uses column-wise storage scheme, i.e., U has types:
 *         Stype = NC, Dtype = D_, Mtype = TRU.
 *
 * B       (input/output) SuperMatrix*
 *         B has types: Stype = DN, Dtype = D_, Mtype = GE.
 *         On entry, the right hand side matrix.
 *         On exit, the solution matrix if info = 0;
 *
 * info    (output) int*
 *	   = 0: successful exit
 *         > 0: if info = i, and i is
 *             <= A->ncol: U(i,i) is exactly zero. The factorization has
 *                been completed, but the factor U is exactly singular,
 *                so the solution could not be computed.
 *             > A->ncol: number of bytes allocated when memory allocation
 *                failure occurred, plus A->ncol.
 *   
 */
    double   t1;	/* Temporary time */
    char     refact[1], trans[1];
    DNformat *Bstore;
    SuperMatrix *AA; /* A in NC format used by the factorization routine.*/
    SuperMatrix AC; /* Matrix postmultiplied by Pc */
    int      lwork = 0, *etree, i;
    
    /* Set default values for some parameters */
    double   diag_pivot_thresh = 1.0;
    double   drop_tol = 0;
    int      panel_size;     /* panel size */
    int      relax;          /* no of columns in a relaxed snodes */
    double   *utime;
    extern SuperLUStat_t SuperLUStat;

    /* Test the input parameters ... */
    *info = 0;
    Bstore = B->Store;
    if ( A->nrow != A->ncol || A->nrow < 0 ||
	 (A->Stype != NC && A->Stype != NR) ||
	 A->Dtype != D_ || A->Mtype != GE )
	*info = -1;
    else if ( B->ncol < 0 || Bstore->lda < SUPERLU_MAX(0, A->nrow) ||
	B->Stype != DN || B->Dtype != D_ || B->Mtype != GE )
	*info = -6;
    if ( *info != 0 ) {
	i = -(*info);
	xerbla_("dgssv", &i);
	return;
    }
    
    *refact = 'N';
    *trans = 'N';
    panel_size = sp_ienv(1);
    relax = sp_ienv(2);

    StatInit(panel_size, relax);
    utime = SuperLUStat.utime;
 
    /* Convert A to NC format when necessary. */
    if ( A->Stype == NR ) {
	NRformat *Astore = A->Store;
	AA = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	dCreate_CompCol_Matrix(AA, A->ncol, A->nrow, Astore->nnz, 
			       Astore->nzval, Astore->colind, Astore->rowptr,
			       NC, A->Dtype, A->Mtype);
	*trans = 'T';
    } else if ( A->Stype == NC ) AA = A;

    etree = intMalloc(A->ncol);

    t1 = SuperLU_timer_();
    sp_preorder(refact, AA, perm_c, etree, &AC);
    utime[ETREE] = SuperLU_timer_() - t1;

    /*printf("Factor PA = LU ... relax %d\tw %d\tmaxsuper %d\trowblk %d\n", 
	  relax, panel_size, sp_ienv(3), sp_ienv(4));*/
    t1 = SuperLU_timer_(); 
    /* Compute the LU factorization of A. */
    dgstrf(refact, &AC, diag_pivot_thresh, drop_tol, relax, panel_size,
	   etree, NULL, lwork, perm_r, perm_c, L, U, info);
    utime[FACT] = SuperLU_timer_() - t1;

    t1 = SuperLU_timer_();
    if ( *info == 0 ) {
        /* Solve the system A*X=B, overwriting B with X. */
        dgstrs (trans, L, U, perm_r, perm_c, B, info);
    }
    utime[SOLVE] = SuperLU_timer_() - t1;

    SUPERLU_FREE (etree);
    Destroy_CompCol_Permuted(&AC);
    if ( A->Stype == NR ) {
	Destroy_SuperMatrix_Store(AA);
	SUPERLU_FREE(AA);
    }

/*     PrintStat( &SuperLUStat ); */
    StatFree();

}
