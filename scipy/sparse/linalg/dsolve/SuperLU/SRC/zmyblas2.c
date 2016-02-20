
/*! @file zmyblas2.c
 * \brief Level 2 Blas operations
 * 
 * <pre>
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 * </pre>
 * Purpose:
 *     Level 2 BLAS operations: solves and matvec, written in C.
 * Note:
 *     This is only used when the system lacks an efficient BLAS library.
 * </pre>
 */
/*
 * File name:		zmyblas2.c
 */
#include "slu_dcomplex.h"

/*! \brief Solves a dense UNIT lower triangular system
 * 
 * The unit lower 
 * triangular matrix is stored in a 2D array M(1:nrow,1:ncol). 
 * The solution will be returned in the rhs vector.
 */
void zlsolve ( int ldm, int ncol, doublecomplex *M, doublecomplex *rhs )
{
    int k;
    doublecomplex x0, x1, x2, x3, temp;
    doublecomplex *M0;
    doublecomplex *Mki0, *Mki1, *Mki2, *Mki3;
    register int firstcol = 0;

    M0 = &M[0];


    while ( firstcol < ncol - 3 ) { /* Do 4 columns */
      	Mki0 = M0 + 1;
      	Mki1 = Mki0 + ldm + 1;
      	Mki2 = Mki1 + ldm + 1;
      	Mki3 = Mki2 + ldm + 1;

      	x0 = rhs[firstcol];
      	zz_mult(&temp, &x0, Mki0); Mki0++;
      	z_sub(&x1, &rhs[firstcol+1], &temp);
      	zz_mult(&temp, &x0, Mki0); Mki0++;
	z_sub(&x2, &rhs[firstcol+2], &temp);
	zz_mult(&temp, &x1, Mki1); Mki1++;
	z_sub(&x2, &x2, &temp);
      	zz_mult(&temp, &x0, Mki0); Mki0++;
	z_sub(&x3, &rhs[firstcol+3], &temp);
	zz_mult(&temp, &x1, Mki1); Mki1++;
	z_sub(&x3, &x3, &temp);
	zz_mult(&temp, &x2, Mki2); Mki2++;
	z_sub(&x3, &x3, &temp);

 	rhs[++firstcol] = x1;
      	rhs[++firstcol] = x2;
      	rhs[++firstcol] = x3;
      	++firstcol;
    
      	for (k = firstcol; k < ncol; k++) {
	    zz_mult(&temp, &x0, Mki0); Mki0++;
	    z_sub(&rhs[k], &rhs[k], &temp);
	    zz_mult(&temp, &x1, Mki1); Mki1++;
	    z_sub(&rhs[k], &rhs[k], &temp);
	    zz_mult(&temp, &x2, Mki2); Mki2++;
	    z_sub(&rhs[k], &rhs[k], &temp);
	    zz_mult(&temp, &x3, Mki3); Mki3++;
	    z_sub(&rhs[k], &rhs[k], &temp);
	}

        M0 += 4 * ldm + 4;
    }

    if ( firstcol < ncol - 1 ) { /* Do 2 columns */
        Mki0 = M0 + 1;
        Mki1 = Mki0 + ldm + 1;

        x0 = rhs[firstcol];
	zz_mult(&temp, &x0, Mki0); Mki0++;
	z_sub(&x1, &rhs[firstcol+1], &temp);

      	rhs[++firstcol] = x1;
      	++firstcol;
    
      	for (k = firstcol; k < ncol; k++) {
	    zz_mult(&temp, &x0, Mki0); Mki0++;
	    z_sub(&rhs[k], &rhs[k], &temp);
	    zz_mult(&temp, &x1, Mki1); Mki1++;
	    z_sub(&rhs[k], &rhs[k], &temp);
	} 
    }
    
}

/*! \brief Solves a dense upper triangular system. 
 *
 * The upper triangular matrix is
 * stored in a 2-dim array M(1:ldm,1:ncol). The solution will be returned
 * in the rhs vector.
 */
void
zusolve ( ldm, ncol, M, rhs )
int ldm;	/* in */
int ncol;	/* in */
doublecomplex *M;	/* in */
doublecomplex *rhs;	/* modified */
{
    doublecomplex xj, temp;
    int jcol, j, irow;

    jcol = ncol - 1;

    for (j = 0; j < ncol; j++) {

	z_div(&xj, &rhs[jcol], &M[jcol + jcol*ldm]); /* M(jcol, jcol) */
	rhs[jcol] = xj;
	
	for (irow = 0; irow < jcol; irow++) {
	    zz_mult(&temp, &xj, &M[irow+jcol*ldm]); /* M(irow, jcol) */
	    z_sub(&rhs[irow], &rhs[irow], &temp);
	}

	jcol--;

    }
}


/*! \brief Performs a dense matrix-vector multiply: Mxvec = Mxvec + M * vec.
 *
 * The input matrix is M(1:nrow,1:ncol); The product is returned in Mxvec[].
 */
void zmatvec ( ldm, nrow, ncol, M, vec, Mxvec )
int ldm;	/* in -- leading dimension of M */
int nrow;	/* in */ 
int ncol;	/* in */
doublecomplex *M;	/* in */
doublecomplex *vec;	/* in */
doublecomplex *Mxvec;	/* in/out */
{
    doublecomplex vi0, vi1, vi2, vi3;
    doublecomplex *M0, temp;
    doublecomplex *Mki0, *Mki1, *Mki2, *Mki3;
    register int firstcol = 0;
    int k;

    M0 = &M[0];

    while ( firstcol < ncol - 3 ) {	/* Do 4 columns */
	Mki0 = M0;
	Mki1 = Mki0 + ldm;
	Mki2 = Mki1 + ldm;
	Mki3 = Mki2 + ldm;

	vi0 = vec[firstcol++];
	vi1 = vec[firstcol++];
	vi2 = vec[firstcol++];
	vi3 = vec[firstcol++];	
	for (k = 0; k < nrow; k++) {
	    zz_mult(&temp, &vi0, Mki0); Mki0++;
	    z_add(&Mxvec[k], &Mxvec[k], &temp);
	    zz_mult(&temp, &vi1, Mki1); Mki1++;
	    z_add(&Mxvec[k], &Mxvec[k], &temp);
	    zz_mult(&temp, &vi2, Mki2); Mki2++;
	    z_add(&Mxvec[k], &Mxvec[k], &temp);
	    zz_mult(&temp, &vi3, Mki3); Mki3++;
	    z_add(&Mxvec[k], &Mxvec[k], &temp);
	}

	M0 += 4 * ldm;
    }

    while ( firstcol < ncol ) {		/* Do 1 column */
 	Mki0 = M0;
	vi0 = vec[firstcol++];
	for (k = 0; k < nrow; k++) {
	    zz_mult(&temp, &vi0, Mki0); Mki0++;
	    z_add(&Mxvec[k], &Mxvec[k], &temp);
	}
	M0 += ldm;
    }
	
}

