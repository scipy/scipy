

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
dgssvx(char *fact, char *trans, char *refact,
       SuperMatrix *A, factor_param_t *factor_params, int *perm_c,
       int *perm_r, int *etree, char *equed, double *R, double *C,
       SuperMatrix *L, SuperMatrix *U, void *work, int lwork,
       SuperMatrix *B, SuperMatrix *X, double *recip_pivot_growth, 
       double *rcond, double *ferr, double *berr, 
       mem_usage_t *mem_usage, int *info )
{
/*
 * Purpose
 * =======
 *
 * DGSSVX solves the system of linear equations A*X=B or A'*X=B, using
 * the LU factorization from dgstrf(). Error bounds on the solution and
 * a condition estimate are also provided. It performs the following steps:
 *
 *   1. If A is stored column-wise (A->Stype = NC):
 *  
 *      1.1. If fact = 'E', scaling factors are computed to equilibrate the
 *           system:
 *             trans = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
 *             trans = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
 *             trans = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
 *           Whether or not the system will be equilibrated depends on the
 *           scaling of the matrix A, but if equilibration is used, A is
 *           overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if trans='N')
 *           or diag(C)*B (if trans = 'T' or 'C').
 *
 *      1.2. Permute columns of A, forming A*Pc, where Pc is a permutation
 *           matrix that usually preserves sparsity.
 *           For more details of this step, see sp_preorder.c.
 *
 *      1.3. If fact = 'N' or 'E', the LU decomposition is used to factor the
 *           matrix A (after equilibration if fact = 'E') as Pr*A*Pc = L*U,
 *           with Pr determined by partial pivoting.
 *
 *      1.4. Compute the reciprocal pivot growth factor.
 *
 *      1.5. If some U(i,i) = 0, so that U is exactly singular, then the
 *           routine returns with info = i. Otherwise, the factored form of 
 *           A is used to estimate the condition number of the matrix A. If
 *           the reciprocal of the condition number is less than machine
 *           precision, info = A->ncol+1 is returned as a warning, but the
 *           routine still goes on to solve for X and computes error bounds
 *           as described below.
 *
 *      1.6. The system of equations is solved for X using the factored form
 *           of A.
 *
 *      1.7. Iterative refinement is applied to improve the computed solution
 *           matrix and calculate error bounds and backward error estimates
 *           for it.
 *
 *      1.8. If equilibration was used, the matrix X is premultiplied by
 *           diag(C) (if trans = 'N') or diag(R) (if trans = 'T' or 'C') so
 *           that it solves the original system before equilibration.
 *
 *   2. If A is stored row-wise (A->Stype = NR), apply the above algorithm
 *      to the transpose of A:
 *
 *      2.1. If fact = 'E', scaling factors are computed to equilibrate the
 *           system:
 *             trans = 'N':  diag(R)*A'*diag(C)     *inv(diag(C))*X = diag(R)*B
 *             trans = 'T': (diag(R)*A'*diag(C))**T *inv(diag(R))*X = diag(C)*B
 *             trans = 'C': (diag(R)*A'*diag(C))**H *inv(diag(R))*X = diag(C)*B
 *           Whether or not the system will be equilibrated depends on the
 *           scaling of the matrix A, but if equilibration is used, A' is
 *           overwritten by diag(R)*A'*diag(C) and B by diag(R)*B 
 *           (if trans='N') or diag(C)*B (if trans = 'T' or 'C').
 *
 *      2.2. Permute columns of transpose(A) (rows of A), 
 *           forming transpose(A)*Pc, where Pc is a permutation matrix that 
 *           usually preserves sparsity.
 *           For more details of this step, see sp_preorder.c.
 *
 *      2.3. If fact = 'N' or 'E', the LU decomposition is used to factor the
 *           transpose(A) (after equilibration if fact = 'E') as 
 *           Pr*transpose(A)*Pc = L*U with the permutation Pr determined by
 *           partial pivoting.
 *
 *      2.4. Compute the reciprocal pivot growth factor.
 *
 *      2.5. If some U(i,i) = 0, so that U is exactly singular, then the
 *           routine returns with info = i. Otherwise, the factored form 
 *           of transpose(A) is used to estimate the condition number of the
 *           matrix A. If the reciprocal of the condition number
 *           is less than machine precision, info = A->nrow+1 is returned as
 *           a warning, but the routine still goes on to solve for X and
 *           computes error bounds as described below.
 *
 *      2.6. The system of equations is solved for X using the factored form
 *           of transpose(A).
 *
 *      2.7. Iterative refinement is applied to improve the computed solution
 *           matrix and calculate error bounds and backward error estimates
 *           for it.
 *
 *      2.8. If equilibration was used, the matrix X is premultiplied by
 *           diag(C) (if trans = 'N') or diag(R) (if trans = 'T' or 'C') so
 *           that it solves the original system before equilibration.
 *
 *   See supermatrix.h for the definition of 'SuperMatrix' structure.
 *
 * Arguments
 * =========
 *
 * fact    (input) char*
 *         Specifies whether or not the factored form of the matrix
 *         A is supplied on entry, and if not, whether the matrix A should
 *         be equilibrated before it is factored.
 *         = 'F': On entry, L, U, perm_r and perm_c contain the factored
 *                form of A. If equed is not 'N', the matrix A has been
 *                equilibrated with scaling factors R and C.
 *                A, L, U, perm_r are not modified.
 *         = 'N': The matrix A will be factored, and the factors will be
 *                stored in L and U.
 *         = 'E': The matrix A will be equilibrated if necessary, then
 *                factored into L and U.
 *
 * trans   (input) char*
 *         Specifies the form of the system of equations:
 *         = 'N': A * X = B        (No transpose)
 *         = 'T': A**T * X = B     (Transpose)
 *         = 'C': A**H * X = B     (Transpose)
 *
 * refact  (input) char*
 *         Specifies whether we want to re-factor the matrix.
 *         = 'N': Factor the matrix A.
 *         = 'Y': Matrix A was factored before, now we want to re-factor
 *                matrix A with perm_r and etree as inputs. Use
 *                the same storage for the L\U factors previously allocated,
 *                expand it if necessary. User should insure to use the same
 *                memory model.  In this case, perm_r may be modified due to
 *                different pivoting determined by diagonal threshold.
 *         If fact = 'F', then refact is not accessed.
 *
 * A       (input/output) SuperMatrix*
 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *         of the linear equations is A->nrow. Currently, the type of A can be:
 *         Stype = NC or NR, Dtype = D_, Mtype = GE. In the future,
 *         more general A can be handled.
 *
 *         On entry, If fact = 'F' and equed is not 'N', then A must have
 *         been equilibrated by the scaling factors in R and/or C.  
 *         A is not modified if fact = 'F' or 'N', or if fact = 'E' and 
 *         equed = 'N' on exit.
 *
 *         On exit, if fact = 'E' and equed is not 'N', A is scaled as follows:
 *         If A->Stype = NC:
 *           equed = 'R':  A := diag(R) * A
 *           equed = 'C':  A := A * diag(C)
 *           equed = 'B':  A := diag(R) * A * diag(C).
 *         If A->Stype = NR:
 *           equed = 'R':  transpose(A) := diag(R) * transpose(A)
 *           equed = 'C':  transpose(A) := transpose(A) * diag(C)
 *           equed = 'B':  transpose(A) := diag(R) * transpose(A) * diag(C).
 *
 * factor_params (input) factor_param_t*
 *         The structure defines the input scalar parameters, consisting of
 *         the following fields. If factor_params = NULL, the default
 *         values are used for all the fields; otherwise, the values
 *         are given by the user.
 *         - panel_size (int): Panel size. A panel consists of at most
 *             panel_size consecutive columns. If panel_size = -1, use 
 *             default value 8.
 *         - relax (int): To control degree of relaxing supernodes. If the
 *             number of nodes (columns) in a subtree of the elimination
 *             tree is less than relax, this subtree is considered as one
 *             supernode, regardless of the row structures of those columns.
 *             If relax = -1, use default value 8.
 *         - diag_pivot_thresh (double): Diagonal pivoting threshold.
 *             At step j of the Gaussian elimination, if
 *                 abs(A_jj) >= diag_pivot_thresh * (max_(i>=j) abs(A_ij)),
 *             then use A_jj as pivot. 0 <= diag_pivot_thresh <= 1.
 *             If diag_pivot_thresh = -1, use default value 1.0,
 *             which corresponds to standard partial pivoting.
 *         - drop_tol (double): Drop tolerance threshold. (NOT IMPLEMENTED)
 *             At step j of the Gaussian elimination, if
 *                 abs(A_ij)/(max_i abs(A_ij)) < drop_tol,
 *             then drop entry A_ij. 0 <= drop_tol <= 1.
 *             If drop_tol = -1, use default value 0.0, which corresponds to
 *             standard Gaussian elimination.
 *
 * perm_c  (input/output) int*
 *	   If A->Stype = NC, Column permutation vector of size A->ncol,
 *         which defines the permutation matrix Pc; perm_c[i] = j means
 *         column i of A is in position j in A*Pc.
 *         On exit, perm_c may be overwritten by the product of the input
 *         perm_c and a permutation that postorders the elimination tree
 *         of Pc'*A'*A*Pc; perm_c is not changed if the elimination tree
 *         is already in postorder.
 *
 *         If A->Stype = NR, column permutation vector of size A->nrow,
 *         which describes permutation of columns of transpose(A) 
 *         (rows of A) as described above.
 * 
 * perm_r  (input/output) int*
 *         If A->Stype = NC, row permutation vector of size A->nrow, 
 *         which defines the permutation matrix Pr, and is determined
 *         by partial pivoting.  perm_r[i] = j means row i of A is in 
 *         position j in Pr*A.
 *
 *         If A->Stype = NR, permutation vector of size A->ncol, which
 *         determines permutation of rows of transpose(A)
 *         (columns of A) as described above.
 *
 *         If refact is not 'Y', perm_r is output argument;
 *         If refact = 'Y', the pivoting routine will try to use the input
 *         perm_r, unless a certain threshold criterion is violated.
 *         In that case, perm_r is overwritten by a new permutation
 *         determined by partial pivoting or diagonal threshold pivoting.
 * 
 * etree   (input/output) int*,  dimension (A->ncol)
 *         Elimination tree of Pc'*A'*A*Pc.
 *         If fact is not 'F' and refact = 'Y', etree is an input argument,
 *         otherwise it is an output argument.
 *         Note: etree is a vector of parent pointers for a forest whose
 *         vertices are the integers 0 to A->ncol-1; etree[root]==A->ncol.
 *
 * equed   (input/output) char*
 *         Specifies the form of equilibration that was done.
 *         = 'N': No equilibration.
 *         = 'R': Row equilibration, i.e., A was premultiplied by diag(R).
 *         = 'C': Column equilibration, i.e., A was postmultiplied by diag(C).
 *         = 'B': Both row and column equilibration, i.e., A was replaced 
 *                by diag(R)*A*diag(C).
 *         If fact = 'F', equed is an input argument, otherwise it is
 *         an output argument.
 *
 * R       (input/output) double*, dimension (A->nrow)
 *         The row scale factors for A or transpose(A).
 *         If equed = 'R' or 'B', A (if A->Stype = NC) or transpose(A) (if
 *             A->Stype = NR) is multiplied on the left by diag(R).
 *         If equed = 'N' or 'C', R is not accessed.
 *         If fact = 'F', R is an input argument; otherwise, R is output.
 *         If fact = 'F' and equed = 'R' or 'B', each element of R must
 *            be positive.
 * 
 * C       (input/output) double*, dimension (A->ncol)
 *         The column scale factors for A or transpose(A).
 *         If equed = 'C' or 'B', A (if A->Stype = NC) or transpose(A) (if 
 *             A->Stype = NR) is multiplied on the right by diag(C).
 *         If equed = 'N' or 'R', C is not accessed.
 *         If fact = 'F', C is an input argument; otherwise, C is output.
 *         If fact = 'F' and equed = 'C' or 'B', each element of C must
 *            be positive.
 *         
 * L       (output) SuperMatrix*
 *	   The factor L from the factorization
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
 * work    (workspace/output) void*, size (lwork) (in bytes)
 *         User supplied workspace, should be large enough
 *         to hold data structures for factors L and U.
 *         On exit, if fact is not 'F', L and U point to this array.
 *
 * lwork   (input) int
 *         Specifies the size of work array in bytes.
 *         = 0:  allocate space internally by system malloc;
 *         > 0:  use user-supplied work array of length lwork in bytes,
 *               returns error if space runs out.
 *         = -1: the routine guesses the amount of space needed without
 *               performing the factorization, and returns it in
 *               mem_usage->total_needed; no other side effects.
 *
 *         See argument 'mem_usage' for memory usage statistics.
 *
 * B       (input/output) SuperMatrix*
 *         B has types: Stype = DN, Dtype = D_, Mtype = GE.
 *         On entry, the right hand side matrix.
 *         On exit,
 *            if equed = 'N', B is not modified; otherwise
 *            if A->Stype = NC:
 *               if trans = 'N' and equed = 'R' or 'B', B is overwritten by
 *                  diag(R)*B;
 *               if trans = 'T' or 'C' and equed = 'C' of 'B', B is
 *                  overwritten by diag(C)*B;
 *            if A->Stype = NR:
 *               if trans = 'N' and equed = 'C' or 'B', B is overwritten by
 *                  diag(C)*B;
 *               if trans = 'T' or 'C' and equed = 'R' of 'B', B is
 *                  overwritten by diag(R)*B.
 *
 * X       (output) SuperMatrix*
 *         X has types: Stype = DN, Dtype = D_, Mtype = GE. 
 *         If info = 0 or info = A->ncol+1, X contains the solution matrix
 *         to the original system of equations. Note that A and B are modified
 *         on exit if equed is not 'N', and the solution to the equilibrated
 *         system is inv(diag(C))*X if trans = 'N' and equed = 'C' or 'B',
 *         or inv(diag(R))*X if trans = 'T' or 'C' and equed = 'R' or 'B'.
 *
 * recip_pivot_growth (output) double*
 *         The reciprocal pivot growth factor max_j( norm(A_j)/norm(U_j) ).
 *         The infinity norm is used. If recip_pivot_growth is much less
 *         than 1, the stability of the LU factorization could be poor.
 *
 * rcond   (output) double*
 *         The estimate of the reciprocal condition number of the matrix A
 *         after equilibration (if done). If rcond is less than the machine
 *         precision (in particular, if rcond = 0), the matrix is singular
 *         to working precision. This condition is indicated by a return
 *         code of info > 0.
 *
 * FERR    (output) double*, dimension (B->ncol)   
 *         The estimated forward error bound for each solution vector   
 *         X(j) (the j-th column of the solution matrix X).   
 *         If XTRUE is the true solution corresponding to X(j), FERR(j) 
 *         is an estimated upper bound for the magnitude of the largest 
 *         element in (X(j) - XTRUE) divided by the magnitude of the   
 *         largest element in X(j).  The estimate is as reliable as   
 *         the estimate for RCOND, and is almost always a slight   
 *         overestimate of the true error.
 *
 * BERR    (output) double*, dimension (B->ncol)
 *         The componentwise relative backward error of each solution   
 *         vector X(j) (i.e., the smallest relative change in   
 *         any element of A or B that makes X(j) an exact solution).
 *
 * mem_usage (output) mem_usage_t*
 *         Record the memory usage statistics, consisting of following fields:
 *         - for_lu (float)
 *           The amount of space used in bytes for L\U data structures.
 *         - total_needed (float)
 *           The amount of space needed in bytes to perform factorization.
 *         - expansions (int)
 *           The number of memory expansions during the LU factorization.
 *
 * info    (output) int*
 *         = 0: successful exit   
 *         < 0: if info = -i, the i-th argument had an illegal value   
 *         > 0: if info = i, and i is   
 *              <= A->ncol: U(i,i) is exactly zero. The factorization has   
 *                    been completed, but the factor U is exactly   
 *                    singular, so the solution and error bounds   
 *                    could not be computed.   
 *              = A->ncol+1: U is nonsingular, but RCOND is less than machine
 *                    precision, meaning that the matrix is singular to
 *                    working precision. Nevertheless, the solution and
 *                    error bounds are computed because there are a number
 *                    of situations where the computed solution can be more
 *                    accurate than the value of RCOND would suggest.   
 *              > A->ncol+1: number of bytes allocated when memory allocation
 *                    failure occurred, plus A->ncol.
 *
 */

    DNformat  *Bstore, *Xstore;
    double    *Bmat, *Xmat;
    int       ldb, ldx, nrhs;
    SuperMatrix *AA; /* A in NC format used by the factorization routine.*/
    SuperMatrix AC; /* Matrix postmultiplied by Pc */
    int       colequ, equil, nofact, notran, rowequ;
    char      trant[1], norm[1];
    int       i, j, info1;
    double    amax, anorm, bignum, smlnum, colcnd, rowcnd, rcmax, rcmin;
    int       relax, panel_size;
    double    diag_pivot_thresh, drop_tol;
    double    t0;      /* temporary time */
    double    *utime;
    extern SuperLUStat_t SuperLUStat;

    /* External functions */
    extern double dlangs(char *, SuperMatrix *);
    extern double dlamch_(char *);

    Bstore = B->Store;
    Xstore = X->Store;
    Bmat   = Bstore->nzval;
    Xmat   = Xstore->nzval;
    ldb    = Bstore->lda;
    ldx    = Xstore->lda;
    nrhs   = B->ncol;

#if 0
printf("dgssvx: fact=%c, trans=%c, refact=%c, equed=%c\n",
       *fact, *trans, *refact, *equed);
#endif
    
    *info = 0;
    nofact = lsame_(fact, "N");
    equil = lsame_(fact, "E");
    notran = lsame_(trans, "N");
    if (nofact || equil) {
	*(unsigned char *)equed = 'N';
	rowequ = FALSE;
	colequ = FALSE;
    } else {
	rowequ = lsame_(equed, "R") || lsame_(equed, "B");
	colequ = lsame_(equed, "C") || lsame_(equed, "B");
	smlnum = dlamch_("Safe minimum");
	bignum = 1. / smlnum;
    }

    /* Test the input parameters */
    if (!nofact && !equil && !lsame_(fact, "F")) *info = -1;
    else if (!notran && !lsame_(trans, "T") && !lsame_(trans, "C")) *info = -2;
    else if ( !(lsame_(refact,"Y") || lsame_(refact, "N")) ) *info = -3;
    else if ( A->nrow != A->ncol || A->nrow < 0 ||
	      (A->Stype != NC && A->Stype != NR) ||
	      A->Dtype != D_ || A->Mtype != GE )
	*info = -4;
    else if (lsame_(fact, "F") && !(rowequ || colequ || lsame_(equed, "N")))
	*info = -9;
    else {
	if (rowequ) {
	    rcmin = bignum;
	    rcmax = 0.;
	    for (j = 0; j < A->nrow; ++j) {
		rcmin = SUPERLU_MIN(rcmin, R[j]);
		rcmax = SUPERLU_MAX(rcmax, R[j]);
	    }
	    if (rcmin <= 0.) *info = -10;
	    else if ( A->nrow > 0)
		rowcnd = SUPERLU_MAX(rcmin,smlnum) / SUPERLU_MIN(rcmax,bignum);
	    else rowcnd = 1.;
	}
	if (colequ && *info == 0) {
	    rcmin = bignum;
	    rcmax = 0.;
	    for (j = 0; j < A->nrow; ++j) {
		rcmin = SUPERLU_MIN(rcmin, C[j]);
		rcmax = SUPERLU_MAX(rcmax, C[j]);
	    }
	    if (rcmin <= 0.) *info = -11;
	    else if (A->nrow > 0)
		colcnd = SUPERLU_MAX(rcmin,smlnum) / SUPERLU_MIN(rcmax,bignum);
	    else colcnd = 1.;
	}
	if (*info == 0) {
	    if ( lwork < -1 ) *info = -15;
	    else if ( B->ncol < 0 || Bstore->lda < SUPERLU_MAX(0, A->nrow) ||
		      B->Stype != DN || B->Dtype != D_ || 
		      B->Mtype != GE )
		*info = -16;
	    else if ( X->ncol < 0 || Xstore->lda < SUPERLU_MAX(0, A->nrow) ||
		      B->ncol != X->ncol || X->Stype != DN ||
		      X->Dtype != D_ || X->Mtype != GE )
		*info = -17;
	}
    }
    if (*info != 0) {
	i = -(*info);
	xerbla_("dgssvx", &i);
	return;
    }
    
    /* Default values for factor_params */
    panel_size = sp_ienv(1);
    relax      = sp_ienv(2);
    diag_pivot_thresh = 1.0;
    drop_tol   = 0.0;
    if ( factor_params != NULL ) {
	if ( factor_params->panel_size != -1 )
	    panel_size = factor_params->panel_size;
	if ( factor_params->relax != -1 ) relax = factor_params->relax;
	if ( factor_params->diag_pivot_thresh != -1 )
	    diag_pivot_thresh = factor_params->diag_pivot_thresh;
	if ( factor_params->drop_tol != -1 )
	    drop_tol = factor_params->drop_tol;
    }

    StatInit(panel_size, relax);
    utime = SuperLUStat.utime;
    
    /* Convert A to NC format when necessary. */
    if ( A->Stype == NR ) {
	NRformat *Astore = A->Store;
	AA = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	dCreate_CompCol_Matrix(AA, A->ncol, A->nrow, Astore->nnz, 
			       Astore->nzval, Astore->colind, Astore->rowptr,
			       NC, A->Dtype, A->Mtype);
	if ( notran ) { /* Reverse the transpose argument. */
	    *trant = 'T';
	    notran = 0;
	} else {
	    *trant = 'N';
	    notran = 1;
	}
    } else { /* A->Stype == NC */
	*trant = *trans;
	AA = A;
    }

    if ( equil ) {
	t0 = SuperLU_timer_();
	/* Compute row and column scalings to equilibrate the matrix A. */
	dgsequ(AA, R, C, &rowcnd, &colcnd, &amax, &info1);
	
	if ( info1 == 0 ) {
	    /* Equilibrate matrix A. */
	    dlaqgs(AA, R, C, rowcnd, colcnd, amax, equed);
	    rowequ = lsame_(equed, "R") || lsame_(equed, "B");
	    colequ = lsame_(equed, "C") || lsame_(equed, "B");
	}
	utime[EQUIL] = SuperLU_timer_() - t0;
    }

    /* Scale the right hand side if equilibration was performed. */
    if ( notran ) {
	if ( rowequ ) {
	    for (j = 0; j < nrhs; ++j)
		for (i = 0; i < A->nrow; ++i) {
		  Bmat[i + j*ldb] *= R[i];
	        }
	}
    } else if ( colequ ) {
	for (j = 0; j < nrhs; ++j)
	    for (i = 0; i < A->nrow; ++i) {
	      Bmat[i + j*ldb] *= C[i];
	    }
    }

    if ( nofact || equil ) {
	
	t0 = SuperLU_timer_();
	sp_preorder(refact, AA, perm_c, etree, &AC);
	utime[ETREE] = SuperLU_timer_() - t0;
    
/*	printf("Factor PA = LU ... relax %d\tw %d\tmaxsuper %d\trowblk %d\n", 
	       relax, panel_size, sp_ienv(3), sp_ienv(4));
	fflush(stdout); */
	
	/* Compute the LU factorization of A*Pc. */
	t0 = SuperLU_timer_();
	dgstrf(refact, &AC, diag_pivot_thresh, drop_tol, relax, panel_size,
	       etree, work, lwork, perm_r, perm_c, L, U, info);
	utime[FACT] = SuperLU_timer_() - t0;
	
	if ( lwork == -1 ) {
	    mem_usage->total_needed = *info - A->ncol;
	    return;
	}
    }

    if ( *info > 0 ) {
	if ( *info <= A->ncol ) {
	    /* Compute the reciprocal pivot growth factor of the leading
	       rank-deficient *info columns of A. */
	    *recip_pivot_growth = dPivotGrowth(*info, AA, perm_c, L, U);
	}
	return;
    }

    /* Compute the reciprocal pivot growth factor *recip_pivot_growth. */
    *recip_pivot_growth = dPivotGrowth(A->ncol, AA, perm_c, L, U);

    /* Estimate the reciprocal of the condition number of A. */
    t0 = SuperLU_timer_();
    if ( notran ) {
	*(unsigned char *)norm = '1';
    } else {
	*(unsigned char *)norm = 'I';
    }
    anorm = dlangs(norm, AA);
    dgscon(norm, L, U, anorm, rcond, info);
    utime[RCOND] = SuperLU_timer_() - t0;
    
    /* Compute the solution matrix X. */
    for (j = 0; j < nrhs; j++)    /* Save a copy of the right hand sides */
	for (i = 0; i < B->nrow; i++)
	    Xmat[i + j*ldx] = Bmat[i + j*ldb];
    
    t0 = SuperLU_timer_();
    dgstrs (trant, L, U, perm_r, perm_c, X, info);
    utime[SOLVE] = SuperLU_timer_() - t0;
    
    /* Use iterative refinement to improve the computed solution and compute
       error bounds and backward error estimates for it. */
    t0 = SuperLU_timer_();
    dgsrfs(trant, AA, L, U, perm_r, perm_c, equed, R, C, B,
	      X, ferr, berr, info);
    utime[REFINE] = SuperLU_timer_() - t0;

    /* Transform the solution matrix X to a solution of the original system. */
    if ( notran ) {
	if ( colequ ) {
	    for (j = 0; j < nrhs; ++j)
		for (i = 0; i < A->nrow; ++i) {
                  Xmat[i + j*ldx] *= C[i];
	        }
	}
    } else if ( rowequ ) {
	for (j = 0; j < nrhs; ++j)
	    for (i = 0; i < A->nrow; ++i) {
	      Xmat[i + j*ldx] *= R[i];
            }
    }

    /* Set INFO = A->ncol+1 if the matrix is singular to working precision. */
    if ( *rcond < dlamch_("E") ) *info = A->ncol + 1;

    dQuerySpace(L, U, panel_size, mem_usage);

    if ( nofact || equil ) Destroy_CompCol_Permuted(&AC);
    if ( A->Stype == NR ) {
	Destroy_SuperMatrix_Store(AA);
	SUPERLU_FREE(AA);
    }

/*     PrintStat( &SuperLUStat ); */
    StatFree();
}
