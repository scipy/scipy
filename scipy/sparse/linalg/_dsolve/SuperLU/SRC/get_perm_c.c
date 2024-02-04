/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*! @file get_perm_c.c
 * \brief Matrix permutation operations
 *
 * <pre>
 * -- SuperLU routine (version 3.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 1, 2008
 * March 25, 2023  add METIS option
 * </pre>
 */
#include "slu_ddefs.h"
#include "colamd.h"


extern int genmmd_(int *neqns, int_t *xadj, int_t *adjncy, 
		   int *invp, int *perm, int_t *delta, int_t *dhead, 
		   int_t *qsize, int_t *llist, int_t *marker, int_t *maxint, 
		   int_t *nofsub);

void
get_colamd(
	   const int m,  /* number of rows in matrix A. */
	   const int n,  /* number of columns in matrix A. */
	   const int_t nnz,/* number of nonzeros in matrix A. */
	   int_t *colptr,  /* column pointer of size n+1 for matrix A. */
	   int_t *rowind,  /* row indices of size nz for matrix A. */
	   int *perm_c   /* out - the column permutation vector. */
	   )
{
    size_t Alen;
    int_t *A, i, *p;
    int info;
    double knobs[COLAMD_KNOBS];
    int_t stats[COLAMD_STATS];

    Alen = COLAMD_recommended(nnz, m, n);

    COLAMD_set_defaults(knobs);

    if ( !(A = intMalloc(Alen)) ) ABORT("Malloc fails for A[]");
    if ( !(p = intMalloc(n+1)) )  ABORT("Malloc fails for p[]");
    for (i = 0; i <= n; ++i) p[i] = colptr[i];
    for (i = 0; i < nnz; ++i) A[i] = rowind[i];
    
    info = COLAMD_MAIN(m, n, Alen, A, p, knobs, stats);

    //printf("after COLAMD_MAIN info %d\n", info);
    if ( info == FALSE ) ABORT("COLAMD failed");

    for (i = 0; i < n; ++i) perm_c[p[i]] = i;

    SUPERLU_FREE(A);
    SUPERLU_FREE(p);
} /* end get_colamd */

void
get_metis(
	  int n,           /* dimension of matrix B */
	  int_t bnz,       /* number of nonzeros in matrix A. */
	  int_t *b_colptr, /* column pointer of size n+1 for matrix B. */
	  int_t *b_rowind, /* row indices of size bnz for matrix B. */
	  int *perm_c      /* out - the column permutation vector. */
	  )
{
#ifdef HAVE_METIS
    /*#define METISOPTIONS 8*/
#define METISOPTIONS 40
    int_t metis_options[METISOPTIONS];
    int numflag = 0; /* C-Style ordering */
    int_t i, nm;
    int_t *perm, *iperm;

    extern int METIS_NodeND(int_t*, int_t*, int_t*, int_t*, int_t*,
			    int_t*, int_t*);

    metis_options[0] = 0; /* Use Defaults for now */

    perm = intMalloc(2*n);
    if (!perm) ABORT("intMalloc fails for perm.");
    iperm = perm + n;
    nm = n;

    /* Call metis */
#undef USEEND
#ifdef USEEND
    METIS_EdgeND(&nm, b_colptr, b_rowind, &numflag, metis_options,
		 perm, iperm);
#else

    /* Earlier version 3.x.x */
    /* METIS_NodeND(&nm, b_colptr, b_rowind, &numflag, metis_options,
       perm, iperm);*/

    /* Latest version 4.x.x */
    METIS_NodeND(&nm, b_colptr, b_rowind, NULL, NULL, perm, iperm);

    /*check_perm_dist("metis perm",  n, perm);*/
#endif

    /* Copy the permutation vector into SuperLU data structure. */
    for (i = 0; i < n; ++i) perm_c[i] = iperm[i];

    SUPERLU_FREE(b_colptr);
    SUPERLU_FREE(b_rowind);
    SUPERLU_FREE(perm);
#endif /* HAVE_METIS */
}

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * Form the structure of A'*A. A is an m-by-n matrix in column oriented
 * format represented by (colptr, rowind). The output A'*A is in column
 * oriented format (symmetrically, also row oriented), represented by
 * (ata_colptr, ata_rowind).
 *
 * This routine is modified from GETATA routine by Tim Davis.
 * The complexity of this algorithm is: SUM_{i=1,m} r(i)^2,
 * i.e., the sum of the square of the row counts.
 *
 * Questions
 * =========
 *     o  Do I need to withhold the *dense* rows?
 *     o  How do I know the number of nonzeros in A'*A?
 * </pre>
 */
void
getata(
       const int m,      /* number of rows in matrix A. */
       const int n,      /* number of columns in matrix A. */
       const int_t nz,     /* number of nonzeros in matrix A */
       int_t *colptr,      /* column pointer of size n+1 for matrix A. */
       int_t *rowind,      /* row indices of size nz for matrix A. */
       int_t *atanz,       /* out - on exit, returns the actual number of
                            nonzeros in matrix A'*A. */
       int_t **ata_colptr, /* out - size n+1 */
       int_t **ata_rowind  /* out - size *atanz */
       )
{
    register int_t i, j, k, col, num_nz, ti, trow;
    int_t *marker, *b_colptr, *b_rowind;
    int_t *t_colptr, *t_rowind; /* a column oriented form of T = A' */

    if ( !(marker = (int_t*) SUPERLU_MALLOC((SUPERLU_MAX(m,n)+1)*sizeof(int_t))) )
	ABORT("SUPERLU_MALLOC fails for marker[]");
    if ( !(t_colptr = (int_t*) SUPERLU_MALLOC((m+1) * sizeof(int_t))) )
	ABORT("SUPERLU_MALLOC t_colptr[]");
    if ( !(t_rowind = (int_t*) SUPERLU_MALLOC(nz * sizeof(int_t))) )
	ABORT("SUPERLU_MALLOC fails for t_rowind[]");

    
    /* Get counts of each column of T, and set up column pointers */
    for (i = 0; i < m; ++i) marker[i] = 0;
    for (j = 0; j < n; ++j) {
	for (i = colptr[j]; i < colptr[j+1]; ++i)
	    ++marker[rowind[i]];
    }
    t_colptr[0] = 0;
    for (i = 0; i < m; ++i) {
	t_colptr[i+1] = t_colptr[i] + marker[i];
	marker[i] = t_colptr[i];
    }

    /* Transpose the matrix from A to T */
    for (j = 0; j < n; ++j)
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    col = rowind[i];
	    t_rowind[marker[col]] = j;
	    ++marker[col];
	}

    
    /* ----------------------------------------------------------------
       compute B = T * A, where column j of B is:

       Struct (B_*j) =    UNION   ( Struct (T_*k) )
                        A_kj != 0

       do not include the diagonal entry
   
       ( Partition A as: A = (A_*1, ..., A_*n)
         Then B = T * A = (T * A_*1, ..., T * A_*n), where
         T * A_*j = (T_*1, ..., T_*m) * A_*j.  )
       ---------------------------------------------------------------- */

    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;

    /* First pass determines number of nonzeros in B */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    /* A_kj is nonzero, add pattern of column T_*k to B_*j */
	    k = rowind[i];
	    for (ti = t_colptr[k]; ti < t_colptr[k+1]; ++ti) {
		trow = t_rowind[ti];
		if ( marker[trow] != j ) {
		    marker[trow] = j;
		    num_nz++;
		}
	    }
	}
    }
    *atanz = num_nz;
    
    /* Allocate storage for A'*A */
    if ( !(*ata_colptr = (int_t*) SUPERLU_MALLOC( (n+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for ata_colptr[]");
    if ( *atanz ) {
	if ( !(*ata_rowind = (int_t*) SUPERLU_MALLOC( *atanz * sizeof(int_t)) ) )
	    ABORT("SUPERLU_MALLOC fails for ata_rowind[]");
    }
    b_colptr = *ata_colptr; /* aliasing */
    b_rowind = *ata_rowind;
    
    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;
    
    /* Compute each column of B, one at a time */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	b_colptr[j] = num_nz;
	
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    /* A_kj is nonzero, add pattern of column T_*k to B_*j */
	    k = rowind[i];
	    for (ti = t_colptr[k]; ti < t_colptr[k+1]; ++ti) {
		trow = t_rowind[ti];
		if ( marker[trow] != j ) {
		    marker[trow] = j;
		    b_rowind[num_nz++] = trow;
		}
	    }
	}
    }
    b_colptr[n] = num_nz;
       
    SUPERLU_FREE(marker);
    SUPERLU_FREE(t_colptr);
    SUPERLU_FREE(t_rowind);
} /* end getata */


/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * Form the structure of A'+A. A is an n-by-n matrix in column oriented
 * format represented by (colptr, rowind). The output A'+A is in column
 * oriented format (symmetrically, also row oriented), represented by
 * (b_colptr, b_rowind).
 * </pre>
 */
void
at_plus_a(
	  const int n,      /* number of columns in matrix A. */
	  const int_t nz,   /* number of nonzeros in matrix A */
	  int_t *colptr,      /* column pointer of size n+1 for matrix A. */
	  int_t *rowind,      /* row indices of size nz for matrix A. */
	  int_t *bnz,         /* out - on exit, returns the actual number of
                               nonzeros in matrix A'*A. */
	  int_t **b_colptr,   /* out - size n+1 */
	  int_t **b_rowind    /* out - size *bnz */
	  )
{
    register int_t i, j, k, col, num_nz;
    int_t *t_colptr, *t_rowind; /* a column oriented form of T = A' */
    int_t *marker;

    if ( !(marker = (int_t*) SUPERLU_MALLOC( n * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for marker[]");
    if ( !(t_colptr = (int_t*) SUPERLU_MALLOC( (n+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for t_colptr[]");
    if ( !(t_rowind = (int_t*) SUPERLU_MALLOC( nz * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails t_rowind[]");

    
    /* Get counts of each column of T, and set up column pointers */
    for (i = 0; i < n; ++i) marker[i] = 0;
    for (j = 0; j < n; ++j) {
	for (i = colptr[j]; i < colptr[j+1]; ++i)
	    ++marker[rowind[i]];
    }
    t_colptr[0] = 0;
    for (i = 0; i < n; ++i) {
	t_colptr[i+1] = t_colptr[i] + marker[i];
	marker[i] = t_colptr[i];
    }

    /* Transpose the matrix from A to T */
    for (j = 0; j < n; ++j)
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    col = rowind[i];
	    t_rowind[marker[col]] = j;
	    ++marker[col];
	}


    /* ----------------------------------------------------------------
       compute B = A + T, where column j of B is:

       Struct (B_*j) = Struct (A_*k) UNION Struct (T_*k)

       do not include the diagonal entry
       ---------------------------------------------------------------- */

    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;

    /* First pass determines number of nonzeros in B */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	/* Add pattern of column A_*k to B_*j */
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    k = rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		++num_nz;
	    }
	}

	/* Add pattern of column T_*k to B_*j */
	for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
	    k = t_rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		++num_nz;
	    }
	}
    }
    *bnz = num_nz;
    
    /* Allocate storage for A+A' */
    if ( !(*b_colptr = (int_t*) SUPERLU_MALLOC( (n+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for b_colptr[]");
    if ( *bnz) {
      if ( !(*b_rowind = (int_t*) SUPERLU_MALLOC( *bnz * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for b_rowind[]");
    }
    
    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;
    
    /* Compute each column of B, one at a time */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	(*b_colptr)[j] = num_nz;
	
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	/* Add pattern of column A_*k to B_*j */
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    k = rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		(*b_rowind)[num_nz++] = k;
	    }
	}

	/* Add pattern of column T_*k to B_*j */
	for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
	    k = t_rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		(*b_rowind)[num_nz++] = k;
	    }
	}
    }
    (*b_colptr)[n] = num_nz;
       
    SUPERLU_FREE(marker);
    SUPERLU_FREE(t_colptr);
    SUPERLU_FREE(t_rowind);
} /* end at_plus_a */

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * GET_PERM_C obtains a permutation matrix Pc, by applying the multiple
 * minimum degree ordering code by Joseph Liu to matrix A'*A or A+A'.
 * or using approximate minimum degree column ordering by Davis et. al.
 * The LU factorization of A*Pc tends to have less fill than the LU 
 * factorization of A.
 *
 * Arguments
 * =========
 *
 * ispec   (input) int
 *         Specifies the type of column ordering to reduce fill:
 *         = 1: minimum degree on the structure of A^T * A
 *         = 2: minimum degree on the structure of A^T + A
 *         = 3: approximate minimum degree for unsymmetric matrices
 *         If ispec == 0, the natural ordering (i.e., Pc = I) is returned.
 * 
 * A       (input) SuperMatrix*
 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *         of the linear equations is A->nrow. Currently, the type of A 
 *         can be: Stype = NC; Dtype = _D; Mtype = GE. In the future,
 *         more general A can be handled.
 *
 * perm_c  (output) int*
 *	   Column permutation vector of size A->ncol, which defines the 
 *         permutation matrix Pc; perm_c[i] = j means column i of A is 
 *         in position j in A*Pc.
 * </pre>
 */
void
get_perm_c(int ispec, SuperMatrix *A, int *perm_c)
{
    NCformat *Astore = A->Store;
    int m, n;
    int_t i, bnz = 0, *b_colptr, *b_rowind;
    int_t delta, maxint, nofsub;
    int *invp;
    int_t *dhead, *qsize, *llist, *marker;
    double t;
    
    m = A->nrow;
    n = A->ncol;

    t = SuperLU_timer_();
    switch ( ispec ) {
    case (NATURAL): /* Natural ordering */
	for (i = 0; i < n; ++i) perm_c[i] = i;
#if ( PRNTlevel>=1 )
	printf("Use natural column ordering.\n");
#endif
	return;
    case (MMD_ATA): /* Minimum degree ordering on A'*A */
	getata(m, n, Astore->nnz, Astore->colptr, Astore->rowind,
		     &bnz, &b_colptr, &b_rowind);
#if ( PRNTlevel>=1 )
	printf("Use minimum degree ordering on A'*A.\n");
#endif
	t = SuperLU_timer_() - t;
	/*printf("Form A'*A time = %8.3f\n", t);*/
	break;
    case (MMD_AT_PLUS_A): /* Minimum degree ordering on A'+A */
	if ( m != n ) ABORT("Matrix is not square");
	at_plus_a(n, Astore->nnz, Astore->colptr, Astore->rowind,
		  &bnz, &b_colptr, &b_rowind);
#if ( PRNTlevel>=1 )
	printf("Use minimum degree ordering on A'+A.\n");
#endif
	t = SuperLU_timer_() - t;
	/*printf("Form A'+A time = %8.3f\n", t);*/
	break;
    case (COLAMD): /* Approximate minimum degree column ordering. */
	get_colamd(m, n, Astore->nnz, Astore->colptr, Astore->rowind, perm_c);
#if ( PRNTlevel>=1 )
	printf(".. Use approximate minimum degree column ordering.\n");
#endif
	return;
#ifdef HAVE_METIS
        case METIS_ATA: /* METIS ordering on A'*A */
	    getata(m, n, Astore->nnz, Astore->colptr, Astore->rowind,
		     &bnz, &b_colptr, &b_rowind);

	    if ( bnz ) { /* non-empty adjacency structure */
		  get_metis(n, bnz, b_colptr, b_rowind, perm_c);
	    } else { /* e.g., diagonal matrix */
		for (i = 0; i < n; ++i) perm_c[i] = i;
		SUPERLU_FREE(b_colptr);
		/* b_rowind is not allocated in this case */
	    }

#if ( PRNTlevel>=1 )
	    printf(".. Use METIS ordering on A'*A\n");
#endif
	    return;
#endif
	
    default:
	ABORT("Invalid ISPEC");
    }

    if ( bnz != 0 ) {
	t = SuperLU_timer_();

	/* Initialize and allocate storage for GENMMD. */
	delta = 0; /* DELTA is a parameter to allow the choice of nodes
		      whose degree <= min-degree + DELTA. */
	maxint = 2147483647; /* 2**31 - 1 */
	invp = (int *) SUPERLU_MALLOC((n+delta)*sizeof(int));
	if ( !invp ) ABORT("SUPERLU_MALLOC fails for invp.");
	dhead = intMalloc(n+delta);
	if ( !dhead ) ABORT("SUPERLU_MALLOC fails for dhead.");
	qsize = intMalloc(n+delta);
	if ( !qsize ) ABORT("SUPERLU_MALLOC fails for qsize.");
	llist = intMalloc(n);
	if ( !llist ) ABORT("SUPERLU_MALLOC fails for llist.");
	marker = intMalloc(n);
	if ( !marker ) ABORT("SUPERLU_MALLOC fails for marker.");

	/* Transform adjacency list into 1-based indexing required by GENMMD.*/
	for (i = 0; i <= n; ++i) ++b_colptr[i];
	for (i = 0; i < bnz; ++i) ++b_rowind[i];
	
	genmmd_(&n, b_colptr, b_rowind, perm_c, invp, &delta, dhead, 
		qsize, llist, marker, &maxint, &nofsub);

	/* Transform perm_c into 0-based indexing. */
	for (i = 0; i < n; ++i) --perm_c[i];

	SUPERLU_FREE(invp);
	SUPERLU_FREE(dhead);
	SUPERLU_FREE(qsize);
	SUPERLU_FREE(llist);
	SUPERLU_FREE(marker);
	SUPERLU_FREE(b_rowind);

	t = SuperLU_timer_() - t;
	/*  printf("call GENMMD time = %8.3f\n", t);*/

    } else { /* Empty adjacency structure */
	for (i = 0; i < n; ++i) perm_c[i] = i;
    }

    SUPERLU_FREE(b_colptr);
}
