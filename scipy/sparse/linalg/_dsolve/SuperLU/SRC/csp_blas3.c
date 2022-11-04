/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file csp_blas3.c
 * \brief Sparse BLAS3, using some dense BLAS3 operations
 *
 * <pre>
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 * </pre>
 */
/*
 * File name:		sp_blas3.c
 * Purpose:		Sparse BLAS3, using some dense BLAS3 operations.
 */

#include "slu_cdefs.h"

/*! \brief
 *
 * <pre>
 * Purpose   
 *   =======   
 * 
 *   sp_c performs one of the matrix-matrix operations   
 * 
 *      C := alpha*op( A )*op( B ) + beta*C,   
 * 
 *   where  op( X ) is one of 
 * 
 *      op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
 * 
 *   alpha and beta are scalars, and A, B and C are matrices, with op( A ) 
 *   an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. 
 *   
 * 
 *   Parameters   
 *   ==========   
 * 
 *   TRANSA - (input) char*
 *            On entry, TRANSA specifies the form of op( A ) to be used in 
 *            the matrix multiplication as follows:   
 *               TRANSA = 'N' or 'n',  op( A ) = A.   
 *               TRANSA = 'T' or 't',  op( A ) = A'.   
 *               TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).   
 *            Unchanged on exit.   
 * 
 *   TRANSB - (input) char*
 *            On entry, TRANSB specifies the form of op( B ) to be used in 
 *            the matrix multiplication as follows:   
 *               TRANSB = 'N' or 'n',  op( B ) = B.   
 *               TRANSB = 'T' or 't',  op( B ) = B'.   
 *               TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).   
 *            Unchanged on exit.   
 * 
 *   M      - (input) int   
 *            On entry,  M  specifies  the number of rows of the matrix 
 *	     op( A ) and of the matrix C.  M must be at least zero. 
 *	     Unchanged on exit.   
 * 
 *   N      - (input) int
 *            On entry,  N specifies the number of columns of the matrix 
 *	     op( B ) and the number of columns of the matrix C. N must be 
 *	     at least zero.
 *	     Unchanged on exit.   
 * 
 *   K      - (input) int
 *            On entry, K specifies the number of columns of the matrix 
 *	     op( A ) and the number of rows of the matrix op( B ). K must 
 *	     be at least  zero.   
 *           Unchanged on exit.
 *      
 *   ALPHA  - (input) complex
 *            On entry, ALPHA specifies the scalar alpha.   
 * 
 *   A      - (input) SuperMatrix*
 *            Matrix A with a sparse format, of dimension (A->nrow, A->ncol).
 *            Currently, the type of A can be:
 *                Stype = NC or NCP; Dtype = SLU_C; Mtype = GE. 
 *            In the future, more general A can be handled.
 * 
 *   B      - COMPLEX PRECISION array of DIMENSION ( LDB, kb ), where kb is 
 *            n when TRANSB = 'N' or 'n',  and is  k otherwise.   
 *            Before entry with  TRANSB = 'N' or 'n',  the leading k by n 
 *            part of the array B must contain the matrix B, otherwise 
 *            the leading n by k part of the array B must contain the 
 *            matrix B.   
 *            Unchanged on exit.   
 * 
 *   LDB    - (input) int
 *            On entry, LDB specifies the first dimension of B as declared 
 *            in the calling (sub) program. LDB must be at least max( 1, n ).  
 *            Unchanged on exit.   
 * 
 *   BETA   - (input) complex
 *            On entry, BETA specifies the scalar beta. When BETA is   
 *            supplied as zero then C need not be set on input.   
 *  
 *   C      - COMPLEX PRECISION array of DIMENSION ( LDC, n ).   
 *            Before entry, the leading m by n part of the array C must 
 *            contain the matrix C,  except when beta is zero, in which 
 *            case C need not be set on entry.   
 *            On exit, the array C is overwritten by the m by n matrix 
 *	     ( alpha*op( A )*B + beta*C ).   
 *  
 *   LDC    - (input) int
 *            On entry, LDC specifies the first dimension of C as declared 
 *            in the calling (sub)program. LDC must be at least max(1,m).   
 *            Unchanged on exit.   
 *  
 *   ==== Sparse Level 3 Blas routine.   
 * </pre>
 */

int
sp_cgemm(char *transa, char *transb, int m, int n, int k, 
         complex alpha, SuperMatrix *A, complex *b, int ldb, 
         complex beta, complex *c, int ldc)
{
    int    incx = 1, incy = 1;
    int    j;

    for (j = 0; j < n; ++j) {
	sp_cgemv(transa, alpha, A, &b[ldb*j], incx, beta, &c[ldc*j], incy);
    }
    return 0;    
}
