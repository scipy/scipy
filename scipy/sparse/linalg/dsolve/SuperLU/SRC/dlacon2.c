/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file dlacon2.c
 * \brief Estimates the 1-norm
 *
 * <pre>
 * -- SuperLU routine (version 5.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * July 25, 2015
 * </pre>
 */
#include <math.h>
#include "slu_Cnames.h"

/*! \brief
 *
 * <pre>
 *   Purpose   
 *   =======   
 *
 *   DLACON2 estimates the 1-norm of a square matrix A.   
 *   Reverse communication is used for evaluating matrix-vector products. 
 * 
 *   This is a thread safe version of DLACON, which uses the array ISAVE
 *   in place of a STATIC variables, as follows:
 *
 *     DLACON     DLACON2
 *      jump     isave[0]
 *      j        isave[1]
 *      iter     isave[2]
 *
 *
 *   Arguments   
 *   =========   
 *
 *   N      (input) INT
 *          The order of the matrix.  N >= 1.   
 *
 *   V      (workspace) DOUBLE PRECISION array, dimension (N)   
 *          On the final return, V = A*W,  where  EST = norm(V)/norm(W)   
 *          (W is not returned).   
 *
 *   X      (input/output) DOUBLE PRECISION array, dimension (N)   
 *          On an intermediate return, X should be overwritten by   
 *                A * X,   if KASE=1,   
 *                A' * X,  if KASE=2,
 *         and DLACON must be re-called with all the other parameters   
 *          unchanged.   
 *
 *   ISGN   (workspace) INT array, dimension (N)
 *
 *   EST    (output) DOUBLE PRECISION   
 *          An estimate (a lower bound) for norm(A).   
 *
 *   KASE   (input/output) INT
 *          On the initial call to DLACON, KASE should be 0.   
 *          On an intermediate return, KASE will be 1 or 2, indicating   
 *          whether X should be overwritten by A * X  or A' * X.   
 *          On the final return from DLACON, KASE will again be 0.   
 *
 *   isave  (input/output) int [3]
 *          ISAVE is INTEGER array, dimension (3)
 *          ISAVE is used to save variables between calls to DLACON2
 *
 *   Further Details   
 *   ===============   
 *
 *   Contributed by Nick Higham, University of Manchester.   
 *   Originally named CONEST, dated March 16, 1988.   
 *
 *   Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of 
 *   a real or complex matrix, with applications to condition estimation", 
 *   ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.   
 *   ===================================================================== 
 * </pre>
 */

int
dlacon2_(int *n, double *v, double *x, int *isgn, double *est, int *kase, int isave[3])
{
    /* Table of constant values */
    int c__1 = 1;
    double      zero = 0.0;
    double      one = 1.0;
    
    /* Local variables */
    int jlast;
    double altsgn, estold;
    int i;
    double temp;
#ifdef _CRAY
    extern int ISAMAX(int *, double *, int *);
    extern double SASUM(int *, double *, int *);
    extern int SCOPY(int *, double *, int *, double *, int *);
#else
    extern int idamax_(int *, double *, int *);
    extern double dasum_(int *, double *, int *);
    extern int dcopy_(int *, double *, int *, double *, int *);
#endif
#define d_sign(a, b) (b >= 0 ? fabs(a) : -fabs(a))    /* Copy sign */
#define i_dnnt(a) \
	( a>=0 ? floor(a+.5) : -floor(.5-a) ) /* Round to nearest integer */

    if ( *kase == 0 ) {
	for (i = 0; i < *n; ++i) {
	    x[i] = 1. / (double) (*n);
	}
	*kase = 1;
	isave[0] = 1;	/* jump = 1; */
	return 0;
    }

    switch (isave[0]) {
	case 1:  goto L20;
	case 2:  goto L40;
	case 3:  goto L70;
	case 4:  goto L110;
	case 5:  goto L140;
    }

    /*     ................ ENTRY   (isave[0] = 1)   
	   FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */
  L20:
    if (*n == 1) {
	v[0] = x[0];
	*est = fabs(v[0]);
	/*        ... QUIT */
	goto L150;
    }
#ifdef _CRAY
    *est = SASUM(n, x, &c__1);
#else
    *est = dasum_(n, x, &c__1);
#endif

    for (i = 0; i < *n; ++i) {
	x[i] = d_sign(one, x[i]);
	isgn[i] = i_dnnt(x[i]);
    }
    *kase = 2;
    isave[0] = 2;  /* jump = 2; */
    return 0;

    /*     ................ ENTRY   (isave[0] = 2)
	   FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */
L40:
#ifdef _CRAY
    isave[1] = ISAMAX(n, &x[0], &c__1);   /* j */
#else
    isave[1] = idamax_(n, &x[0], &c__1);  /* j */
#endif
    --isave[1];  /* --j; */
    isave[2] = 2;  /* iter = 2; */

    /*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */
L50:
    for (i = 0; i < *n; ++i) x[i] = zero;
    x[isave[1]] = one;
    *kase = 1;
    isave[0] = 3;  /* jump = 3; */
    return 0;

    /*     ................ ENTRY   (isave[0] = 3)   
	   X HAS BEEN OVERWRITTEN BY A*X. */
L70:
#ifdef _CRAY
    SCOPY(n, x, &c__1, v, &c__1);
#else
    dcopy_(n, x, &c__1, v, &c__1);
#endif
    estold = *est;
#ifdef _CRAY
    *est = SASUM(n, v, &c__1);
#else
    *est = dasum_(n, v, &c__1);
#endif

    for (i = 0; i < *n; ++i)
	if (i_dnnt(d_sign(one, x[i])) != isgn[i])
	    goto L90;

    /*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
    goto L120;

L90:
    /*     TEST FOR CYCLING. */
    if (*est <= estold) goto L120;

    for (i = 0; i < *n; ++i) {
	x[i] = d_sign(one, x[i]);
	isgn[i] = i_dnnt(x[i]);
    }
    *kase = 2;
    isave[0] = 4;  /* jump = 4; */
    return 0;

    /*     ................ ENTRY   (isave[0] = 4)   
	   X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */
L110:
    jlast = isave[1];  /* j; */
#ifdef _CRAY
    isave[1] = ISAMAX(n, &x[0], &c__1);  /* j */
#else
    isave[1] = idamax_(n, &x[0], &c__1);  /* j */
#endif
    isave[1] = isave[1] - 1;  /* --j; */
    if (x[jlast] != fabs(x[isave[1]]) && isave[2] < 5) {
	isave[2] = isave[2] + 1;  /* ++iter; */
	goto L50;
    }

    /*     ITERATION COMPLETE.  FINAL STAGE. */
L120:
    altsgn = 1.;
    for (i = 1; i <= *n; ++i) {
	x[i-1] = altsgn * ((double)(i - 1) / (double)(*n - 1) + 1.);
	altsgn = -altsgn;
    }
    *kase = 1;
    isave[0] = 5;  /* jump = 5; */
    return 0;
    
    /*     ................ ENTRY   (isave[0] = 5)   
	   X HAS BEEN OVERWRITTEN BY A*X. */
L140:
#ifdef _CRAY
    temp = SASUM(n, x, &c__1) / (double)(*n * 3) * 2.;
#else
    temp = dasum_(n, x, &c__1) / (double)(*n * 3) * 2.;
#endif
    if (temp > *est) {
#ifdef _CRAY
	SCOPY(n, &x[0], &c__1, &v[0], &c__1);
#else
	dcopy_(n, &x[0], &c__1, &v[0], &c__1);
#endif
	*est = temp;
    }

L150:
    *kase = 0;
    return 0;

} /* dlacon_ */
