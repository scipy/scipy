/* mc64ad.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "slu_ddefs.h"

#define abs(a) ((a) >= 0) ? (a) : -(a)
#define min(a,b) ((a) < (b)) ? (a) : (b)

#if 0
/* Table of constant values */
static int_t c__1 = 1;
static int_t c__2 = 2;
#endif

/* CCCC COPYRIGHT (c) 1999  Council for the Central Laboratory of the */
/* CCCC Research Councils.    All rights reserved. */
/* CCCC PACKAGE MC64A/AD */
/* CCCC AUTHORS Iain Duff (i.duff@rl.ac.uk) and Jacko Koster (jak@ii.uib.no) */
/* CCCC LAST UPDATE 20/09/99 */
/* CCCC */
/* *** Conditions on external use *** */

/* The user shall acknowledge the contribution of this */
/* package in any publication of material dependent upon the use of */
/* the package. The user shall use reasonable endeavours to notify */
/* the authors of the package of this publication. */

/* The user can modify this code but, at no time */
/* shall the right or title to all or any part of this package pass */
/* to the user. The user shall make available free of charge */
/* to the authors for any purpose all information relating to any */
/* alteration or addition made to this package for the purposes of */
/* extending the capabilities or enhancing the performance of this */
/* package. */

/* The user shall not pass this code directly to a third party without the */
/* express prior consent of the authors.  Users wanting to licence their */
/* own copy of these routines should send email to hsl@stfc.ac.uk */

/* None of the comments from the Copyright notice up to and including this */
/* one shall be removed or altered in any way. */
/* ********************************************************************** */
/* Subroutine */ int_t mc64id_(int_t *icntl)
{
    int_t i__;


/* *** Copyright (c) 1999  Council for the Central Laboratory of the */
/*     Research Councils                                             *** */
/* *** Although every effort has been made to ensure robustness and  *** */
/* *** reliability of the subroutines in this MC64 suite, we         *** */
/* *** disclaim any liability arising through the use or misuse of   *** */
/* *** any of the subroutines.                                       *** */
/* *** Any problems?   Contact ... */
/*     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   *** */

/*  Purpose */
/*  ======= */

/*  The components of the array ICNTL control the action of MC64A/AD. */
/*  Default values for these are set in this subroutine. */

/*  Parameters */
/*  ========== */


/*  Local variables */

/*    ICNTL(1) has default value 6. */
/*     It is the output stream for error messages. If it */
/*     is negative, these messages will be suppressed. */

/*    ICNTL(2) has default value 6. */
/*     It is the output stream for warning messages. */
/*     If it is negative, these messages are suppressed. */

/*    ICNTL(3) has default value -1. */
/*     It is the output stream for monitoring printing. */
/*     If it is negative, these messages are suppressed. */

/*    ICNTL(4) has default value 0. */
/*     If left at the defaut value, the incoming data is checked for */
/*     out-of-range indices and duplicates.  Setting ICNTL(4) to any */
/*     other will avoid the checks but is likely to cause problems */
/*     later if out-of-range indices or duplicates are present. */
/*     The user should only set ICNTL(4) non-zero, if the data is */
/*     known to avoid these problems. */

/*    ICNTL(5) to ICNTL(10) are not used by MC64A/AD but are set to */
/*     zero in this routine. */
/* Initialization of the ICNTL array. */
    /* Parameter adjustments */
    --icntl;

    /* Function Body */
    icntl[1] = 6;
    icntl[2] = 6;
    icntl[3] = -1;
    for (i__ = 4; i__ <= 10; ++i__) {
	icntl[i__] = 0;
/* L10: */
    }
    return 0;
} /* mc64id_ */

/* ********************************************************************** */
/* Subroutine */ int_t mc64ad_(int_t *job, int_t *n, int_t *ne, int_t *
	ip, int_t *irn, double *a, int_t *num, int *cperm, 
	int_t *liw, int_t *iw, int_t *ldw, double *dw, int_t *
	icntl, int_t *info)
{
    /* System generated locals */
    int_t i__1, i__2;
    double d__1, d__2;

    /* Builtin functions */
    double log(double);

    /* Local variables */
    int_t i__, j, k;
    double fact, rinf;

    extern /* Subroutine */ int_t mc21ad_(int_t *, int_t *, int_t *, 
	    int_t *, int_t *, int_t *, int_t *, int_t *),
	mc64bd_(int_t *n, int_t *ne, int_t *ip, int_t *irn, double *a,
		int *iperm, int_t *num, int_t *jperm, 
		int_t *pr, int_t *q, int_t *l, double *d__),
	mc64rd_(int_t *n, int_t *ne, int_t *ip, int_t *irn, double *a),
	mc64sd_(int_t *n, int_t *ne, int_t *ip, int_t *irn, double *a,
		int *iperm, int_t *numx, int_t *w, int_t *len, int_t *lenl,
		int_t *lenh, int_t *fc, int_t *iw, int_t *iw4),
	mc64wd_(int_t *n, int_t *ne, int_t *ip, int_t *irn, double *a,
		int *iperm, int_t *num, int_t *jperm, int_t *out,
		int_t *pr, int_t *q, int_t *l, double *u, double *d__);

/* *** Copyright (c) 1999  Council for the Central Laboratory of the */
/*     Research Councils                                             *** */
/* *** Although every effort has been made to ensure robustness and  *** */
/* *** reliability of the subroutines in this MC64 suite, we         *** */
/* *** disclaim any liability arising through the use or misuse of   *** */
/* *** any of the subroutines.                                       *** */
/* *** Any problems?   Contact ... */
/*     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   *** */

/*  Purpose */
/*  ======= */

/* This subroutine attempts to find a column permutation for an NxN */
/* sparse matrix A = {a_ij} that makes the permuted matrix have N */
/* entries on its diagonal. */
/* If the matrix is structurally nonsingular, the subroutine optionally */
/* returns a column permutation that maximizes the smallest element */
/* on the diagonal, maximizes the sum of the diagonal entries, or */
/* maximizes the product of the diagonal entries of the permuted matrix. */
/* For the latter option, the subroutine also finds scaling factors */
/* that may be used to scale the matrix so that the nonzero diagonal */
/* entries of the permuted matrix are one in absolute value and all the */
/* off-diagonal entries are less than or equal to one in absolute value. */
/* The natural logarithms of the scaling factors u(i), i=1..N, for the */
/* rows and v(j), j=1..N, for the columns are returned so that the */
/* scaled matrix B = {b_ij} has entries b_ij = a_ij * EXP(u_i + v_j). */

/*  Parameters */
/*  ========== */


/* JOB is an INT_T variable which must be set by the user to */
/* control the action. It is not altered by the subroutine. */
/* Possible values for JOB are: */
/*   1 Compute a column permutation of the matrix so that the */
/*     permuted matrix has as many entries on its diagonal as possible. */
/*     The values on the diagonal are of arbitrary size. HSL subroutine */
/*     MC21A/AD is used for this. See [1]. */
/*   2 Compute a column permutation of the matrix so that the smallest */
/*     value on the diagonal of the permuted matrix is maximized. */
/*     See [3]. */
/*   3 Compute a column permutation of the matrix so that the smallest */
/*     value on the diagonal of the permuted matrix is maximized. */
/*     The algorithm differs from the one used for JOB = 2 and may */
/*     have quite a different performance. See [2]. */
/*   4 Compute a column permutation of the matrix so that the sum */
/*     of the diagonal entries of the permuted matrix is maximized. */
/*     See [3]. */
/*   5 Compute a column permutation of the matrix so that the product */
/*     of the diagonal entries of the permuted matrix is maximized */
/*     and vectors to scale the matrix so that the nonzero diagonal */
/*     entries of the permuted matrix are one in absolute value and */
/*     all the off-diagonal entries are less than or equal to one in */
/*     absolute value. See [3]. */
/*  Restriction: 1 <= JOB <= 5. */

/* N is an INT_T variable which must be set by the user to the */
/*   order of the matrix A. It is not altered by the subroutine. */
/*   Restriction: N >= 1. */

/* NE is an INT_T variable which must be set by the user to the */
/*   number of entries in the matrix. It is not altered by the */
/*   subroutine. */
/*   Restriction: NE >= 1. */

/* IP is an INT_T array of length N+1. */
/*   IP(J), J=1..N, must be set by the user to the position in array IRN */
/*   of the first row index of an entry in column J. IP(N+1) must be set */
/*   to NE+1. It is not altered by the subroutine. */

/* IRN is an INT_T array of length NE. */
/*   IRN(K), K=1..NE, must be set by the user to hold the row indices of */
/*   the entries of the matrix. Those belonging to column J must be */
/*   stored contiguously in the positions IP(J)..IP(J+1)-1. The ordering */
/*   of the row indices within each column is unimportant. Repeated */
/*   entries are not allowed. The array IRN is not altered by the */
/*   subroutine. */

/* A is a REAL (DOUBLE PRECISION in the D-version) array of length NE. */
/*   The user must set A(K), K=1..NE, to the numerical value of the */
/*   entry that corresponds to IRN(K). */
/*   It is not used by the subroutine when JOB = 1. */
/*   It is not altered by the subroutine. */

/* NUM is an INT_T variable that need not be set by the user. */
/*   On successful exit, NUM will be the number of entries on the */
/*   diagonal of the permuted matrix. */
/*   If NUM < N, the matrix is structurally singular. */

/* CPERM is an INT_T array of length N that need not be set by the */
/*   user. On successful exit, CPERM contains the column permutation. */
/*   Column CPERM(J) of the original matrix is column J in the permuted */
/*   matrix, J=1..N. */

/* LIW is an INT_T variable that must be set by the user to */
/*   the dimension of array IW. It is not altered by the subroutine. */
/*   Restriction: */
/*     JOB = 1 :  LIW >= 5N */
/*     JOB = 2 :  LIW >= 4N */
/*     JOB = 3 :  LIW >= 10N + NE */
/*     JOB = 4 :  LIW >= 5N */
/*     JOB = 5 :  LIW >= 5N */

/* IW is an INT_T array of length LIW that is used for workspace. */

/* LDW is an INT_T variable that must be set by the user to the */
/*   dimension of array DW. It is not altered by the subroutine. */
/*   Restriction: */
/*     JOB = 1 :  LDW is not used */
/*     JOB = 2 :  LDW >= N */
/*     JOB = 3 :  LDW >= NE */
/*     JOB = 4 :  LDW >= 2N + NE */
/*     JOB = 5 :  LDW >= 3N + NE */

/* DW is a REAL (DOUBLE PRECISION in the D-version) array of length LDW */
/*   that is used for workspace. If JOB = 5, on return, */
/*   DW(i) contains u_i, i=1..N, and DW(N+j) contains v_j, j=1..N. */

/* ICNTL is an INT_T array of length 10. Its components control the */
/*   output of MC64A/AD and must be set by the user before calling */
/*   MC64A/AD. They are not altered by the subroutine. */

/*   ICNTL(1) must be set to specify the output stream for */
/*   error messages. If ICNTL(1) < 0, messages are suppressed. */
/*   The default value set by MC46I/ID is 6. */

/*   ICNTL(2) must be set by the user to specify the output stream for */
/*   warning messages. If ICNTL(2) < 0, messages are suppressed. */
/*   The default value set by MC46I/ID is 6. */

/*   ICNTL(3) must be set by the user to specify the output stream for */
/*   diagnostic messages. If ICNTL(3) < 0, messages are suppressed. */
/*   The default value set by MC46I/ID is -1. */

/*   ICNTL(4) must be set by the user to a value other than 0 to avoid */
/*   checking of the input data. */
/*   The default value set by MC46I/ID is 0. */

/* INFO is an INT_T array of length 10 which need not be set by the */
/*   user. INFO(1) is set non-negative to indicate success. A negative */
/*   value is returned if an error occurred, a positive value if a */
/*   warning occurred. INFO(2) holds further information on the error. */
/*   On exit from the subroutine, INFO(1) will take one of the */
/*   following values: */
/*    0 : successful entry (for structurally nonsingular matrix). */
/*   +1 : successful entry (for structurally singular matrix). */
/*   +2 : the returned scaling factors are large and may cause */
/*        overflow when used to scale the matrix. */
/*        (For JOB = 5 entry only.) */
/*   -1 : JOB < 1 or JOB > 5.  Value of JOB held in INFO(2). */
/*   -2 : N < 1.  Value of N held in INFO(2). */
/*   -3 : NE < 1. Value of NE held in INFO(2). */
/*   -4 : the defined length LIW violates the restriction on LIW. */
/*        Value of LIW required given by INFO(2). */
/*   -5 : the defined length LDW violates the restriction on LDW. */
/*        Value of LDW required given by INFO(2). */
/*   -6 : entries are found whose row indices are out of range. INFO(2) */
/*        contains the index of a column in which such an entry is found. */
/*   -7 : repeated entries are found. INFO(2) contains the index of a */
/*        column in which such entries are found. */
/*  INFO(3) to INFO(10) are not currently used and are set to zero by */
/*        the routine. */

/* References: */
/*  [1]  I. S. Duff, (1981), */
/*       "Algorithm 575. Permutations for a zero-free diagonal", */
/*       ACM Trans. Math. Software 7(3), 387-390. */
/*  [2]  I. S. Duff and J. Koster, (1998), */
/*       "The design and use of algorithms for permuting large */
/*       entries to the diagonal of sparse matrices", */
/*       SIAM J. Matrix Anal. Appl., vol. 20, no. 4, pp. 889-901. */
/*  [3]  I. S. Duff and J. Koster, (1999), */
/*       "On algorithms for permuting large entries to the diagonal */
/*       of sparse matrices", */
/*       Technical Report RAL-TR-1999-030, RAL, Oxfordshire, England. */
/* Local variables and parameters */
/* External routines and functions */
/*     EXTERNAL FD05AD */
/*     DOUBLE PRECISION FD05AD */
/* Intrinsic functions */
/* Set RINF to largest positive real number (infinity) */
/* XSL    RINF = FD05AD(5) */
    /* Parameter adjustments */
    --cperm;
    --ip;
    --a;
    --irn;
    --iw;
    --dw;
    --icntl;
    --info;

    /* Function Body */
    rinf = dmach("Overflow");
/* Check value of JOB */
    if (*job < 1 || *job > 5) {
	info[1] = -1;
	info[2] = *job;
	if (icntl[1] >= 0) {
	    printf(" ****** Error in MC64A/AD. INFO(1) = %2d"
		   " because JOB = %d\n", (int) info[1], (int) *job);
	}
	goto L99;
    }
/* Check value of N */
    if (*n < 1) {
	info[1] = -2;
	info[2] = *n;
	if (icntl[1] >= 0) {
	    printf(" ****** Error in MC64A/AD. INFO(1) = %2d"
		   " because N = %d\n", (int) info[1], (int) *job);
	}
	goto L99;
    }
/* Check value of NE */
    if (*ne < 1) {
	info[1] = -3;
	info[2] = *ne;
	if (icntl[1] >= 0) {
	    printf(" ****** Error in MC64A/AD. INFO(1) = %2d"
		   " because NE = %d\n", (int) info[1], (int) *job);
	}
	goto L99;
    }
/* Check LIW */
    if (*job == 1) {
	k = *n * 5;
    }
    if (*job == 2) {
	k = *n << 2;
    }
    if (*job == 3) {
	k = *n * 10 + *ne;
    }
    if (*job == 4) {
	k = *n * 5;
    }
    if (*job == 5) {
	k = *n * 5;
    }
    if (*liw < k) {
	info[1] = -4;
	info[2] = k;
	if (icntl[1] >= 0) {
	    printf(" ****** Error in MC64A/AD. INFO(1) = %2d"
		   " LIW too small, must be at least %8d\n", (int) info[1], (int) k);
	}
	goto L99;
    }
/* Check LDW */
/* If JOB = 1, do not check */
    if (*job > 1) {
	if (*job == 2) {
	    k = *n;
	}
	if (*job == 3) {
	    k = *ne;
	}
	if (*job == 4) {
	    k = (*n << 1) + *ne;
	}
	if (*job == 5) {
	    k = *n * 3 + *ne;
	}
	if (*ldw < k) {
	    info[1] = -5;
	    info[2] = k;
	    if (icntl[1] >= 0) {
		printf(" ****** Error in MC64A/AD. INFO(1) = %2d"
		       " LDW too small, must be at least %8d\n", (int) info[1], (int) k);
	    }
	    goto L99;
	}
    }
    if (icntl[4] == 0) {
/* Check row indices. Use IW(1:N) as workspace */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iw[i__] = 0;
/* L3: */
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = ip[j + 1] - 1;
	    for (k = ip[j]; k <= i__2; ++k) {
		i__ = irn[k];
/* Check for row indices that are out of range */
		if (i__ < 1 || i__ > *n) {
		    info[1] = -6;
		    info[2] = j;
		    if (icntl[1] >= 0) {
			printf(" ****** Error in MC64A/AD. INFO(1) = %2d Column %8d"
			       " contains an entry with invalid row index %8d\n",
			       (int) info[1], (int) j, (int) i__);
		    }
		    goto L99;
		}
/* Check for repeated row indices within a column */
		if (iw[i__] == j) {
		    info[1] = -7;
		    info[2] = j;
		    if (icntl[1] >= 0) {
			printf(" ****** Error in MC64A/AD. INFO(1) = %2d"
			       "        Column %8d"
			       " contains two or more entries with row index %8d\n",
			       (int) info[1], (int) j, (int) i__);
		    }
		    goto L99;
		} else {
		    iw[i__] = j;
		}
/* L4: */
	    }
/* L6: */
	    }
    }
/* Print diagnostics on input */
    if (icntl[3] >= 0) {
	printf("  ****** Input parameters for MC64A/AD: JOB = %8d,"
	       " N = %d, NE = %8d\n", (int) *job, (int) *n, (int) *ne);
	printf(" IP(1:N+1)   = ");
	for (j=1; j<=(*n+1); ++j) {
	    printf("%8d", (int) ip[j]);
	    if (j%8 == 0) printf("\n");
	}
	printf("\n IRN(1:NE) = ");
	for (j=1; j<=(*ne); ++j) {
	    printf("%8d", (int) irn[j]);
	    if (j%8 == 0) printf("\n");
	}
	printf("\n");

	if (*job > 1) {
	    printf(" A(1:NE)     = ");
	    for (j=1; j<=(*ne); ++j) {
		printf("%f14.4", a[j]);
		if (j%4 == 0) printf("\n");
	    }
	    printf("\n");
	}
    }
/* Set components of INFO to zero */
    for (i__ = 1; i__ <= 10; ++i__) {
	info[i__] = 0;
/* L8: */
    }
/* Compute maximum matching with MC21A/AD */
    if (*job == 1) {
/* Put length of column J in IW(J) */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    iw[j] = ip[j + 1] - ip[j];
/* L10: */
	}
/* IW(N+1:5N) is workspace */
#if 0
	mc21ad_(n, &irn[1], ne, &ip[1], &iw[1], &cperm[1], num, &iw[*n+1]);
#else
	printf(" ****** Warning from MC64A/AD. Need to link mc21ad.\n");
#endif
	goto L90;
    }
/* Compute bottleneck matching */
    if (*job == 2) {
/* IW(1:5N), DW(1:N) are workspaces */
	mc64bd_(n, ne, &ip[1], &irn[1], &a[1], &cperm[1], num, &iw[1], &iw[*n 
		+ 1], &iw[(*n << 1) + 1], &iw[*n * 3 + 1], &dw[1]);
	goto L90;
    }
/* Compute bottleneck matching */
    if (*job == 3) {
/* Copy IRN(K) into IW(K), ABS(A(K)) into DW(K), K=1..NE */
	i__1 = *ne;
	for (k = 1; k <= i__1; ++k) {
	    iw[k] = irn[k];
	    dw[k] = (d__1 = a[k], abs(d__1));
/* L20: */
	}
/* Sort entries in each column by decreasing value. */
	mc64rd_(n, ne, &ip[1], &iw[1], &dw[1]);
/* IW(NE+1:NE+10N) is workspace */
	mc64sd_(n, ne, &ip[1], &iw[1], &dw[1], &cperm[1], num, &iw[*ne + 1], &
		iw[*ne + *n + 1], &iw[*ne + (*n << 1) + 1], &iw[*ne + *n * 3 
		+ 1], &iw[*ne + (*n << 2) + 1], &iw[*ne + *n * 5 + 1], &iw[*
		ne + *n * 6 + 1]);
	goto L90;
    }
    if (*job == 4) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    fact = 0.;
	    i__2 = ip[j + 1] - 1;
	    for (k = ip[j]; k <= i__2; ++k) {
		if ((d__1 = a[k], abs(d__1)) > fact) {
		    fact = (d__2 = a[k], abs(d__2));
		}
/* L30: */
	    }
	    i__2 = ip[j + 1] - 1;
	    for (k = ip[j]; k <= i__2; ++k) {
		dw[(*n << 1) + k] = fact - (d__1 = a[k], abs(d__1));
/* L40: */
	    }
/* L50: */
	}
/* B = DW(2N+1:2N+NE); IW(1:5N) and DW(1:2N) are workspaces */
	mc64wd_(n, ne, &ip[1], &irn[1], &dw[(*n << 1) + 1], &cperm[1], num, &
		iw[1], &iw[*n + 1], &iw[(*n << 1) + 1], &iw[*n * 3 + 1], &iw[(
		*n << 2) + 1], &dw[1], &dw[*n + 1]);
	goto L90;
    }
    if (*job == 5) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    fact = 0.;
	    i__2 = ip[j + 1] - 1;
	    for (k = ip[j]; k <= i__2; ++k) {
		dw[*n * 3 + k] = (d__1 = a[k], abs(d__1));
		if (dw[*n * 3 + k] > fact) {
		    fact = dw[*n * 3 + k];
		}
/* L60: */
	    }
	    dw[(*n << 1) + j] = fact;
	    if (fact != 0.) {
		fact = log(fact);
	    } else {
		fact = rinf / *n;
	    }
	    i__2 = ip[j + 1] - 1;
	    for (k = ip[j]; k <= i__2; ++k) {
		if (dw[*n * 3 + k] != 0.) {
		    dw[*n * 3 + k] = fact - log(dw[*n * 3 + k]);
		} else {
		    dw[*n * 3 + k] = rinf / *n;
		}
/* L70: */
	    }
/* L75: */
	}
/* B = DW(3N+1:3N+NE); IW(1:5N) and DW(1:2N) are workspaces */
	mc64wd_(n, ne, &ip[1], &irn[1], &dw[*n * 3 + 1], &cperm[1], num, &iw[
		1], &iw[*n + 1], &iw[(*n << 1) + 1], &iw[*n * 3 + 1], &iw[(*n 
		<< 2) + 1], &dw[1], &dw[*n + 1]);
	if (*num == *n) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (dw[(*n << 1) + j] != 0.) {
		    dw[*n + j] -= log(dw[(*n << 1) + j]);
		} else {
		    dw[*n + j] = 0.;
		}
/* L80: */
	    }
	}
/* Check size of scaling factors */
	fact = log(rinf) * .5f;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (dw[j] < fact && dw[*n + j] < fact) {
		goto L86;
	    }
	    info[1] = 2;
	    goto L90;
L86:
	    ;
	}
/*       GO TO 90 */
    }
L90:
    if (info[1] == 0 && *num < *n) {
/* Matrix is structurally singular, return with warning */
	info[1] = 1;
	if (icntl[2] >= 0) {
	    printf(" ****** Warning from MC64A/AD. INFO(1) = %2d"
		   " The matrix is structurally singular.\n", (int)info[1]);
	}
    }
    if (info[1] == 2) {
/* Scaling factors are large, return with warning */
	if (icntl[2] >= 0) {
	    printf(" ****** Warning from MC64A/AD. INFO(1) = %2d\n"
		   "        Some scaling factors may be too large.\n", (int) info[1]);
	}
    }
/* Print diagnostics on output */
    if (icntl[3] >= 0) {
	printf(" ****** Output parameters for MC64A/AD: INFO(1:2)  = %8d%8d\n",
	       (int) info[1], (int) info[2]);
	printf(" NUM        = %8d", (int) *num);
	printf(" CPERM(1:N) = ");
	for (j=1; j<=*n; ++j) {
	    printf("%8d", (int) cperm[j]);
	    if (j%8 == 0) printf("\n");
	}
	if (*job == 5) {
	    printf("\n DW(1:N)    = ");
	    for (j=1; j<=*n; ++j) {
		printf("%11.3f", dw[j]);
		if (j%5 == 0) printf("\n");
	    }
	    printf("\n DW(N+1:2N) = ");
	    for (j=1; j<=*n; ++j) {
		printf("%11.3f", dw[*n+j]);
		if (j%5 == 0) printf("\n");
	    }
	    printf("\n");
	}
    }
/* Return from subroutine. */
L99:
    return 0;
} /* mc64ad_ */

/* ********************************************************************** */
/* Subroutine */ int_t mc64bd_(int_t *n, int_t *ne, int_t *ip, int_t *
	irn, double *a, int *iperm, int_t *num, int_t *jperm, 
	int_t *pr, int_t *q, int_t *l, double *d__)
{
    /* System generated locals */
    int_t i__1, i__2, i__3;
    double d__1, d__2, d__3;
    int_t c__1 = 1;

    /* Local variables */
    int_t i__, j, k;
    double a0;
    int_t i0, q0;
    double ai, di;
    int_t ii, jj, kk;
    double bv;
    int_t up;
    double dq0;
    int_t kk1, kk2;
    double csp;
    int_t isp, jsp, low;
    double dnew;
    int_t jord, qlen, idum, jdum;
    double rinf;
    extern /* Subroutine */ int_t
	mc64dd_(int_t *, int_t *, int_t *, double *, int_t *, int_t *),
	mc64ed_(int_t *, int_t *, int_t *, double *, int_t *, int_t *),
	mc64fd_(int_t *, int_t *, int_t *, int_t *, double *, int_t *, int_t *);


/* *** Copyright (c) 1999  Council for the Central Laboratory of the */
/*     Research Councils                                             *** */
/* *** Although every effort has been made to ensure robustness and  *** */
/* *** reliability of the subroutines in this MC64 suite, we         *** */
/* *** disclaim any liability arising through the use or misuse of   *** */
/* *** any of the subroutines.                                       *** */
/* *** Any problems?   Contact ... */
/*     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   *** */

/* N, NE, IP, IRN are described in MC64A/AD. */
/* A is a REAL (DOUBLE PRECISION in the D-version) array of length */
/*   NE. A(K), K=1..NE, must be set to the value of the entry */
/*   that corresponds to IRN(K). It is not altered. */
/* IPERM is an INT_T array of length N. On exit, it contains the */
/*    matching: IPERM(I) = 0 or row I is matched to column IPERM(I). */
/* NUM is INT_T variable. On exit, it contains the cardinality of the */
/*    matching stored in IPERM. */
/* IW is an INT_T work array of length 4N. */
/* DW is a REAL (DOUBLE PRECISION in D-version) work array of length N. */
/* Local variables */
/* Local parameters */
/* Intrinsic functions */
/* External subroutines and/or functions */
/*      EXTERNAL FD05AD,MC64DD,MC64ED,MC64FD, DMACH */
/*      DOUBLE PRECISION FD05AD, DMACH */
/* Set RINF to largest positive real number */
/* XSL  RINF = FD05AD(5) */
    /* Parameter adjustments */
    --d__;
    --l;
    --q;
    --pr;
    --jperm;
    --iperm;
    --ip;
    --a;
    --irn;

    /* Function Body */
    rinf = dmach("Overflow");
/* Initialization */
    *num = 0;
    bv = rinf;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	iperm[k] = 0;
	jperm[k] = 0;
	pr[k] = ip[k];
	d__[k] = 0.;
/* L10: */
    }
/* Scan columns of matrix; */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	a0 = -1.;
	i__2 = ip[j + 1] - 1;
	for (k = ip[j]; k <= i__2; ++k) {
	    i__ = irn[k];
	    ai = (d__1 = a[k], abs(d__1));
	    if (ai > d__[i__]) {
		d__[i__] = ai;
	    }
	    if (jperm[j] != 0) {
		goto L30;
	    }
	    if (ai >= bv) {
		a0 = bv;
		if (iperm[i__] != 0) {
		    goto L30;
		}
		jperm[j] = i__;
		iperm[i__] = j;
		++(*num);
	    } else {
		if (ai <= a0) {
		    goto L30;
		}
		a0 = ai;
		i0 = i__;
	    }
L30:
	    ;
	}
	if (a0 != -1. && a0 < bv) {
	    bv = a0;
	    if (iperm[i0] != 0) {
		goto L20;
	    }
	    iperm[i0] = j;
	    jperm[j] = i0;
	    ++(*num);
	}
L20:
	;
    }
/* Update BV with smallest of all the largest maximum absolute values */
/* of the rows. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	d__1 = bv, d__2 = d__[i__];
	bv = min(d__1,d__2);
/* L25: */
    }
    if (*num == *n) {
	goto L1000;
    }
/* Rescan unassigned columns; improve initial assignment */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (jperm[j] != 0) {
	    goto L95;
	}
	i__2 = ip[j + 1] - 1;
	for (k = ip[j]; k <= i__2; ++k) {
	    i__ = irn[k];
	    ai = (d__1 = a[k], abs(d__1));
	    if (ai < bv) {
		goto L50;
	    }
	    if (iperm[i__] == 0) {
		goto L90;
	    }
	    jj = iperm[i__];
	    kk1 = pr[jj];
	    kk2 = ip[jj + 1] - 1;
	    if (kk1 > kk2) {
		goto L50;
	    }
	    i__3 = kk2;
	    for (kk = kk1; kk <= i__3; ++kk) {
		ii = irn[kk];
		if (iperm[ii] != 0) {
		    goto L70;
		}
		if ((d__1 = a[kk], abs(d__1)) >= bv) {
		    goto L80;
		}
L70:
		;
	    }
	    pr[jj] = kk2 + 1;
L50:
	    ;
	}
	goto L95;
L80:
	jperm[jj] = ii;
	iperm[ii] = jj;
	pr[jj] = kk + 1;
L90:
	++(*num);
	jperm[j] = i__;
	iperm[i__] = j;
	pr[j] = k + 1;
L95:
	;
    }
    if (*num == *n) {
	goto L1000;
    }
/* Prepare for main loop */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = -1.;
	l[i__] = 0;
/* L99: */
    }
/* Main loop ... each pass round this loop is similar to Dijkstra's */
/* algorithm for solving the single source shortest path problem */
    i__1 = *n;
    for (jord = 1; jord <= i__1; ++jord) {
	if (jperm[jord] != 0) {
	    goto L100;
	}
	qlen = 0;
	low = *n + 1;
	up = *n + 1;
/* CSP is cost of shortest path to any unassigned row */
/* ISP is matrix position of unassigned row element in shortest path */
/* JSP is column index of unassigned row element in shortest path */
	csp = -1.;
/* Build shortest path tree starting from unassigned column JORD */
	j = jord;
	pr[j] = -1;
/* Scan column J */
	i__2 = ip[j + 1] - 1;
	for (k = ip[j]; k <= i__2; ++k) {
	    i__ = irn[k];
	    dnew = (d__1 = a[k], abs(d__1));
	    if (csp >= dnew) {
		goto L115;
	    }
	    if (iperm[i__] == 0) {
/* Row I is unassigned; update shortest path info */
		csp = dnew;
		isp = i__;
		jsp = j;
		if (csp >= bv) {
		    goto L160;
		}
	    } else {
		d__[i__] = dnew;
		if (dnew >= bv) {
/* Add row I to Q2 */
		    --low;
		    q[low] = i__;
		} else {
/* Add row I to Q, and push it */
		    ++qlen;
		    l[i__] = qlen;
		    mc64dd_(&i__, n, &q[1], &d__[1], &l[1], &c__1);
		}
		jj = iperm[i__];
		pr[jj] = j;
	    }
L115:
	    ;
	}
	i__2 = *num;
	for (jdum = 1; jdum <= i__2; ++jdum) {
/* If Q2 is empty, extract new rows from Q */
	    if (low == up) {
		if (qlen == 0) {
		    goto L160;
		}
		i__ = q[1];
		if (csp >= d__[i__]) {
		    goto L160;
		}
		bv = d__[i__];
		i__3 = *n;
		for (idum = 1; idum <= i__3; ++idum) {
		    mc64ed_(&qlen, n, &q[1], &d__[1], &l[1], &c__1);
		    l[i__] = 0;
		    --low;
		    q[low] = i__;
		    if (qlen == 0) {
			goto L153;
		    }
		    i__ = q[1];
		    if (d__[i__] != bv) {
			goto L153;
		    }
/* L152: */
		}
/* End of dummy loop; this point is never reached */
	    }
/* Move row Q0 */
L153:
	    --up;
	    q0 = q[up];
	    dq0 = d__[q0];
	    l[q0] = up;
/* Scan column that matches with row Q0 */
	    j = iperm[q0];
	    i__3 = ip[j + 1] - 1;
	    for (k = ip[j]; k <= i__3; ++k) {
		i__ = irn[k];
/* Update D(I) */
		if (l[i__] >= up) {
		    goto L155;
		}
/* Computing MIN */
		d__2 = dq0, d__3 = (d__1 = a[k], abs(d__1));
		dnew = min(d__2,d__3);
		if (csp >= dnew) {
		    goto L155;
		}
		if (iperm[i__] == 0) {
/* Row I is unassigned; update shortest path info */
		    csp = dnew;
		    isp = i__;
		    jsp = j;
		    if (csp >= bv) {
			goto L160;
		    }
		} else {
		    di = d__[i__];
		    if (di >= bv || di >= dnew) {
			goto L155;
		    }
		    d__[i__] = dnew;
		    if (dnew >= bv) {
/* Delete row I from Q (if necessary); add row I to Q2 */
			if (di != -1.) {
			    mc64fd_(&l[i__], &qlen, n, &q[1], &d__[1], &l[1], 
				    &c__1);
			}
			l[i__] = 0;
			--low;
			q[low] = i__;
		    } else {
/* Add row I to Q (if necessary); push row I up Q */
			if (di == -1.) {
			    ++qlen;
			    l[i__] = qlen;
			}
			mc64dd_(&i__, n, &q[1], &d__[1], &l[1], &c__1);
		    }
/* Update tree */
		    jj = iperm[i__];
		    pr[jj] = j;
		}
L155:
		;
	    }
/* L150: */
	}
/* If CSP = MINONE, no augmenting path is found */
L160:
	if (csp == -1.) {
	    goto L190;
	}
/* Update bottleneck value */
	bv = min(bv,csp);
/* Find augmenting path by tracing backward in PR; update IPERM,JPERM */
	++(*num);
	i__ = isp;
	j = jsp;
	i__2 = *num + 1;
	for (jdum = 1; jdum <= i__2; ++jdum) {
	    i0 = jperm[j];
	    jperm[j] = i__;
	    iperm[i__] = j;
	    j = pr[j];
	    if (j == -1) {
		goto L190;
	    }
	    i__ = i0;
/* L170: */
	}
/* End of dummy loop; this point is never reached */
L190:
	i__2 = *n;
	for (kk = up; kk <= i__2; ++kk) {
	    i__ = q[kk];
	    d__[i__] = -1.;
	    l[i__] = 0;
/* L191: */
	}
	i__2 = up - 1;
	for (kk = low; kk <= i__2; ++kk) {
	    i__ = q[kk];
	    d__[i__] = -1.;
/* L192: */
	}
	i__2 = qlen;
	for (kk = 1; kk <= i__2; ++kk) {
	    i__ = q[kk];
	    d__[i__] = -1.;
	    l[i__] = 0;
/* L193: */
	}
L100:
	;
    }
/* End of main loop */
/* BV is bottleneck value of final matching */
    if (*num == *n) {
	goto L1000;
    }
/* Matrix is structurally singular, complete IPERM. */
/* JPERM, PR are work arrays */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jperm[j] = 0;
/* L300: */
    }
    k = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iperm[i__] == 0) {
	    ++k;
	    pr[k] = i__;
	} else {
	    j = iperm[i__];
	    jperm[j] = i__;
	}
/* L310: */
    }
    k = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (jperm[i__] != 0) {
	    goto L320;
	}
	++k;
	jdum = pr[k];
	iperm[jdum] = i__;
L320:
	;
    }
L1000:
    return 0;
} /* mc64bd_ */

/* ********************************************************************** */
/* Subroutine */ int_t mc64dd_(int_t *i__, int_t *n, int_t *q, double 
	*d__, int_t *l, int_t *iway)
{
    /* System generated locals */
    int_t i__1;

    /* Local variables */
    double di;
    int_t qk, pos, idum, posk;


/* *** Copyright (c) 1999  Council for the Central Laboratory of the */
/*     Research Councils                                             *** */
/* *** Although every effort has been made to ensure robustness and  *** */
/* *** reliability of the subroutines in this MC64 suite, we         *** */
/* *** disclaim any liability arising through the use or misuse of   *** */
/* *** any of the subroutines.                                       *** */
/* *** Any problems?   Contact ... */
/*     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   *** */

/* Variables N,Q,D,L are described in MC64B/BD */
/* IF IWAY is equal to 1, then */
/* node I is pushed from its current position upwards */
/* IF IWAY is not equal to 1, then */
/* node I is pushed from its current position downwards */
/* Local variables and parameters */
    /* Parameter adjustments */
    --l;
    --d__;
    --q;

    /* Function Body */
    di = d__[*i__];
    pos = l[*i__];
/* POS is index of current position of I in the tree */
    if (*iway == 1) {
	i__1 = *n;
	for (idum = 1; idum <= i__1; ++idum) {
	    if (pos <= 1) {
		goto L20;
	    }
	    posk = pos / 2;
	    qk = q[posk];
	    if (di <= d__[qk]) {
		goto L20;
	    }
	    q[pos] = qk;
	    l[qk] = pos;
	    pos = posk;
/* L10: */
	}
/* End of dummy loop; this point is never reached */
    } else {
	i__1 = *n;
	for (idum = 1; idum <= i__1; ++idum) {
	    if (pos <= 1) {
		goto L20;
	    }
	    posk = pos / 2;
	    qk = q[posk];
	    if (di >= d__[qk]) {
		goto L20;
	    }
	    q[pos] = qk;
	    l[qk] = pos;
	    pos = posk;
/* L15: */
	}
/* End of dummy loop; this point is never reached */
    }
/* End of dummy if; this point is never reached */
L20:
    q[pos] = *i__;
    l[*i__] = pos;
    return 0;
} /* mc64dd_ */

/* ********************************************************************** */
/* Subroutine */ int_t mc64ed_(int_t *qlen, int_t *n, int_t *q, 
	double *d__, int_t *l, int_t *iway)
{
    /* System generated locals */
    int_t i__1;

    /* Local variables */
    int_t i__;
    double di, dk, dr;
    int_t pos, idum, posk;


/* *** Copyright (c) 1999  Council for the Central Laboratory of the */
/*     Research Councils                                             *** */
/* *** Although every effort has been made to ensure robustness and  *** */
/* *** reliability of the subroutines in this MC64 suite, we         *** */
/* *** disclaim any liability arising through the use or misuse of   *** */
/* *** any of the subroutines.                                       *** */
/* *** Any problems?   Contact ... */
/*     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   *** */

/* Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or */
/*     MC64W/WD (IWAY = 2) */
/* The root node is deleted from the binary heap. */
/* Local variables and parameters */
/* Move last element to begin of Q */
    /* Parameter adjustments */
    --l;
    --d__;
    --q;

    /* Function Body */
    i__ = q[*qlen];
    di = d__[i__];
    --(*qlen);
    pos = 1;
    if (*iway == 1) {
	i__1 = *n;
	for (idum = 1; idum <= i__1; ++idum) {
	    posk = pos << 1;
	    if (posk > *qlen) {
		goto L20;
	    }
	    dk = d__[q[posk]];
	    if (posk < *qlen) {
		dr = d__[q[posk + 1]];
		if (dk < dr) {
		    ++posk;
		    dk = dr;
		}
	    }
	    if (di >= dk) {
		goto L20;
	    }
/* Exchange old last element with larger priority child */
	    q[pos] = q[posk];
	    l[q[pos]] = pos;
	    pos = posk;
/* L10: */
	}
/* End of dummy loop; this point is never reached */
    } else {
	i__1 = *n;
	for (idum = 1; idum <= i__1; ++idum) {
	    posk = pos << 1;
	    if (posk > *qlen) {
		goto L20;
	    }
	    dk = d__[q[posk]];
	    if (posk < *qlen) {
		dr = d__[q[posk + 1]];
		if (dk > dr) {
		    ++posk;
		    dk = dr;
		}
	    }
	    if (di <= dk) {
		goto L20;
	    }
/* Exchange old last element with smaller child */
	    q[pos] = q[posk];
	    l[q[pos]] = pos;
	    pos = posk;
/* L15: */
	}
/* End of dummy loop; this point is never reached */
    }
/* End of dummy if; this point is never reached */
L20:
    q[pos] = i__;
    l[i__] = pos;
    return 0;
} /* mc64ed_ */

/* ********************************************************************** */
/* Subroutine */ int_t mc64fd_(int_t *pos0, int_t *qlen, int_t *n, 
	int_t *q, double *d__, int_t *l, int_t *iway)
{
    /* System generated locals */
    int_t i__1;

    /* Local variables */
    int_t i__;
    double di, dk, dr;
    int_t qk, pos, idum, posk;


/* *** Copyright (c) 1999  Council for the Central Laboratory of the */
/*     Research Councils                                             *** */
/* *** Although every effort has been made to ensure robustness and  *** */
/* *** reliability of the subroutines in this MC64 suite, we         *** */
/* *** disclaim any liability arising through the use or misuse of   *** */
/* *** any of the subroutines.                                       *** */
/* *** Any problems?   Contact ... */
/*     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   *** */

/* Variables QLEN,N,Q,D,L are described in MC64B/BD (IWAY = 1) or */
/*     MC64WD (IWAY = 2). */
/* Move last element in the heap */
/* Quick return, if possible */
    /* Parameter adjustments */
    --l;
    --d__;
    --q;

    /* Function Body */
    if (*qlen == *pos0) {
	--(*qlen);
	return 0;
    }
/* Move last element from queue Q to position POS0 */
/* POS is current position of node I in the tree */
    i__ = q[*qlen];
    di = d__[i__];
    --(*qlen);
    pos = *pos0;
    if (*iway == 1) {
	i__1 = *n;
	for (idum = 1; idum <= i__1; ++idum) {
	    if (pos <= 1) {
		goto L20;
	    }
	    posk = pos / 2;
	    qk = q[posk];
	    if (di <= d__[qk]) {
		goto L20;
	    }
	    q[pos] = qk;
	    l[qk] = pos;
	    pos = posk;
/* L10: */
	}
/* End of dummy loop; this point is never reached */
L20:
	q[pos] = i__;
	l[i__] = pos;
	i__1 = *n;
	for (idum = 1; idum <= i__1; ++idum) {
	    posk = pos << 1;
	    if (posk > *qlen) {
		goto L40;
	    }
	    dk = d__[q[posk]];
	    if (posk < *qlen) {
		dr = d__[q[posk + 1]];
		if (dk < dr) {
		    ++posk;
		    dk = dr;
		}
	    }
	    if (di >= dk) {
		goto L40;
	    }
	    qk = q[posk];
	    q[pos] = qk;
	    l[qk] = pos;
	    pos = posk;
/* L30: */
	}
/* End of dummy loop; this point is never reached */
    } else {
	i__1 = *n;
	for (idum = 1; idum <= i__1; ++idum) {
	    if (pos <= 1) {
		goto L34;
	    }
	    posk = pos / 2;
	    qk = q[posk];
	    if (di >= d__[qk]) {
		goto L34;
	    }
	    q[pos] = qk;
	    l[qk] = pos;
	    pos = posk;
/* L32: */
	}
/* End of dummy loop; this point is never reached */
L34:
	q[pos] = i__;
	l[i__] = pos;
	i__1 = *n;
	for (idum = 1; idum <= i__1; ++idum) {
	    posk = pos << 1;
	    if (posk > *qlen) {
		goto L40;
	    }
	    dk = d__[q[posk]];
	    if (posk < *qlen) {
		dr = d__[q[posk + 1]];
		if (dk > dr) {
		    ++posk;
		    dk = dr;
		}
	    }
	    if (di <= dk) {
		goto L40;
	    }
	    qk = q[posk];
	    q[pos] = qk;
	    l[qk] = pos;
	    pos = posk;
/* L36: */
	}
/* End of dummy loop; this point is never reached */
    }
/* End of dummy if; this point is never reached */
L40:
    q[pos] = i__;
    l[i__] = pos;
    return 0;
} /* mc64fd_ */

/* ********************************************************************** */
/* Subroutine */ int_t mc64rd_(int_t *n, int_t *ne, int_t *ip, int_t *
	irn, double *a)
{
    /* System generated locals */
    int_t i__1, i__2, i__3;

    /* Local variables */
    int_t j, k, r__, s;
    double ha;
    int_t hi, td, mid, len, ipj;
    double key;
    int_t last, todo[50], first;


/* *** Copyright (c) 1999  Council for the Central Laboratory of the */
/*     Research Councils                                             *** */
/* *** Although every effort has been made to ensure robustness and  *** */
/* *** reliability of the subroutines in this MC64 suite, we         *** */
/* *** disclaim any liability arising through the use or misuse of   *** */
/* *** any of the subroutines.                                       *** */
/* *** Any problems?   Contact ... */
/*     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   *** */

/* This subroutine sorts the entries in each column of the */
/* sparse matrix (defined by N,NE,IP,IRN,A) by decreasing */
/* numerical value. */
/* Local constants */
/* Local variables */
/* Local arrays */
    /* Parameter adjustments */
    --ip;
    --a;
    --irn;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	len = ip[j + 1] - ip[j];
	if (len <= 1) {
	    goto L100;
	}
	ipj = ip[j];
/* Sort array roughly with partial quicksort */
	if (len < 15) {
	    goto L400;
	}
	todo[0] = ipj;
	todo[1] = ipj + len;
	td = 2;
L500:
	first = todo[td - 2];
	last = todo[td - 1];
/* KEY is the smallest of two values present in interval [FIRST,LAST) */
	key = a[(first + last) / 2];
	i__2 = last - 1;
	for (k = first; k <= i__2; ++k) {
	    ha = a[k];
	    if (ha == key) {
		goto L475;
	    }
	    if (ha > key) {
		goto L470;
	    }
	    key = ha;
	    goto L470;
L475:
	    ;
	}
/* Only one value found in interval, so it is already sorted */
	td += -2;
	goto L425;
/* Reorder interval [FIRST,LAST) such that entries before MID are gt KEY */
L470:
	mid = first;
	i__2 = last - 1;
	for (k = first; k <= i__2; ++k) {
	    if (a[k] <= key) {
		goto L450;
	    }
	    ha = a[mid];
	    a[mid] = a[k];
	    a[k] = ha;
	    hi = irn[mid];
	    irn[mid] = irn[k];
	    irn[k] = hi;
	    ++mid;
L450:
	    ;
	}
/* Both subintervals [FIRST,MID), [MID,LAST) are nonempty */
/* Stack the longest of the two subintervals first */
	if (mid - first >= last - mid) {
	    todo[td + 1] = last;
	    todo[td] = mid;
	    todo[td - 1] = mid;
/*          TODO(TD-1) = FIRST */
	} else {
	    todo[td + 1] = mid;
	    todo[td] = first;
	    todo[td - 1] = last;
	    todo[td - 2] = mid;
	}
	td += 2;
L425:
	if (td == 0) {
	    goto L400;
	}
/* There is still work to be done */
	if (todo[td - 1] - todo[td - 2] >= 15) {
	    goto L500;
	}
/* Next interval is already short enough for straightforward insertion */
	td += -2;
	goto L425;
/* Complete sorting with straightforward insertion */
L400:
	i__2 = ipj + len - 1;
	for (r__ = ipj + 1; r__ <= i__2; ++r__) {
	    if (a[r__ - 1] < a[r__]) {
		ha = a[r__];
		hi = irn[r__];
		a[r__] = a[r__ - 1];
		irn[r__] = irn[r__ - 1];
		i__3 = ipj + 1;
		for (s = r__ - 1; s >= i__3; --s) {
		    if (a[s - 1] < ha) {
			a[s] = a[s - 1];
			irn[s] = irn[s - 1];
		    } else {
			a[s] = ha;
			irn[s] = hi;
			goto L200;
		    }
/* L300: */
		}
		a[ipj] = ha;
		irn[ipj] = hi;
	    }
L200:
	    ;
	}
L100:
	;
    }
    return 0;
} /* mc64rd_ */

/* ********************************************************************** */
/* Subroutine */ int_t mc64sd_(int_t *n, int_t *ne, int_t *ip, int_t *
	irn, double *a, int *iperm, int_t *numx, int_t *w, 
	int_t *len, int_t *lenl, int_t *lenh, int_t *fc, int_t *iw, 
	int_t *iw4)
{
    /* System generated locals */
    int_t i__1, i__2, i__3, i__4;

    /* Local variables */
    int_t i__, j, k, l, ii, mod, cnt, num;
    double bval, bmin, bmax, rinf;
    int_t nval, wlen, idum1, idum2, idum3;
    extern /* Subroutine */ int_t
	mc64qd_(int_t *, int_t *, int_t *, int_t *, int_t *, double *,
		int_t *, double *), 
	mc64ud_(int_t *, int_t *, int_t *, int_t *, int_t *, 
		int_t *, int_t *, int_t *, int_t *, int_t *, int_t *, 
		int_t *, int_t *, int_t *, int_t *);

/* *** Copyright (c) 1999  Council for the Central Laboratory of the */
/*     Research Councils                                             *** */
/* *** Although every effort has been made to ensure robustness and  *** */
/* *** reliability of the subroutines in this MC64 suite, we         *** */
/* *** disclaim any liability arising through the use or misuse of   *** */
/* *** any of the subroutines.                                       *** */
/* *** Any problems?   Contact ... */
/*     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   *** */

/* N, NE, IP, IRN, are described in MC64A/AD. */
/* A is a REAL (DOUBLE PRECISION in the D-version) array of length NE. */
/*   A(K), K=1..NE, must be set to the value of the entry that */
/*   corresponds to IRN(k). The entries in each column must be */
/*   non-negative and ordered by decreasing value. */
/* IPERM is an INT_T array of length N. On exit, it contains the */
/*   bottleneck matching: IPERM(I) - 0 or row I is matched to column */
/*   IPERM(I). */
/* NUMX is an INT_T variable. On exit, it contains the cardinality */
/*   of the matching stored in IPERM. */
/* IW is an INT_T work array of length 10N. */
/* FC is an int_t array of length N that contains the list of */
/*   unmatched columns. */
/* LEN(J), LENL(J), LENH(J) are int_t arrays of length N that point */
/*   to entries in matrix column J. */
/*   In the matrix defined by the column parts IP(J)+LENL(J) we know */
/*   a matching does not exist; in the matrix defined by the column */
/*   parts IP(J)+LENH(J) we know one exists. */
/*   LEN(J) lies between LENL(J) and LENH(J) and determines the matrix */
/*   that is tested for a maximum matching. */
/* W is an int_t array of length N and contains the indices of the */
/*   columns for which LENL ne LENH. */
/* WLEN is number of indices stored in array W. */
/* IW is int_t work array of length N. */
/* IW4 is int_t work array of length 4N used by MC64U/UD. */
/*      EXTERNAL FD05AD,MC64QD,MC64UD */
/*      DOUBLE PRECISION FD05AD */
/* BMIN and BMAX are such that a maximum matching exists for the input */
/*   matrix in which all entries smaller than BMIN are dropped. */
/*   For BMAX, a maximum matching does not exist. */
/* BVAL is a value between BMIN and BMAX. */
/* CNT is the number of calls made to MC64U/UD so far. */
/* NUM is the cardinality of last matching found. */
/* Set RINF to largest positive real number */
/* XSL      RINF = FD05AD(5) */
    /* Parameter adjustments */
    --iw4;
    --iw;
    --fc;
    --lenh;
    --lenl;
    --len;
    --w;
    --iperm;
    --ip;
    --a;
    --irn;

    /* Function Body */
    rinf = dmach("Overflow");
/* Compute a first maximum matching from scratch on whole matrix. */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	fc[j] = j;
	iw[j] = 0;
	len[j] = ip[j + 1] - ip[j];
/* L20: */
    }
/* The first call to MC64U/UD */
    cnt = 1;
    mod = 1;
    *numx = 0;
    mc64ud_(&cnt, &mod, n, &irn[1], ne, &ip[1], &len[1], &fc[1], &iw[1], numx,
	     n, &iw4[1], &iw4[*n + 1], &iw4[(*n << 1) + 1], &iw4[*n * 3 + 1]);
/* IW contains a maximum matching of length NUMX. */
    num = *numx;
    if (num != *n) {
/* Matrix is structurally singular */
	bmax = rinf;
    } else {
/* Matrix is structurally nonsingular, NUM=NUMX=N; */
/* Set BMAX just above the smallest of all the maximum absolute */
/* values of the columns */
	bmax = rinf;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    bval = 0.f;
	    i__2 = ip[j + 1] - 1;
	    for (k = ip[j]; k <= i__2; ++k) {
		if (a[k] > bval) {
		    bval = a[k];
		}
/* L25: */
	    }
	    if (bval < bmax) {
		bmax = bval;
	    }
/* L30: */
	}
	bmax *= 1.001f;
    }
/* Initialize BVAL,BMIN */
    bval = 0.f;
    bmin = 0.f;
/* Initialize LENL,LEN,LENH,W,WLEN according to BMAX. */
/* Set LEN(J), LENH(J) just after last entry in column J. */
/* Set LENL(J) just after last entry in column J with value ge BMAX. */
    wlen = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = ip[j + 1] - ip[j];
	lenh[j] = l;
	len[j] = l;
	i__2 = ip[j + 1] - 1;
	for (k = ip[j]; k <= i__2; ++k) {
	    if (a[k] < bmax) {
		goto L46;
	    }
/* L45: */
	}
/* Column J is empty or all entries are ge BMAX */
	k = ip[j + 1];
L46:
	lenl[j] = k - ip[j];
/* Add J to W if LENL(J) ne LENH(J) */
	if (lenl[j] == l) {
	    goto L48;
	}
	++wlen;
	w[wlen] = j;
L48:
	;
    }
/* Main loop */
    i__1 = *ne;
    for (idum1 = 1; idum1 <= i__1; ++idum1) {
	if (num == *numx) {
/* We have a maximum matching in IW; store IW in IPERM */
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		iperm[i__] = iw[i__];
/* L50: */
	    }
/* Keep going round this loop until matching IW is no longer maximum. */
	    i__2 = *ne;
	    for (idum2 = 1; idum2 <= i__2; ++idum2) {
		bmin = bval;
		if (bmax == bmin) {
		    goto L99;
		}
/* Find splitting value BVAL */
		mc64qd_(&ip[1], &lenl[1], &len[1], &w[1], &wlen, &a[1], &nval,
			 &bval);
		if (nval <= 1) {
		    goto L99;
		}
/* Set LEN such that all matrix entries with value lt BVAL are */
/* discarded. Store old LEN in LENH. Do this for all columns W(K). */
/* Each step, either K is incremented or WLEN is decremented. */
		k = 1;
		i__3 = *n;
		for (idum3 = 1; idum3 <= i__3; ++idum3) {
		    if (k > wlen) {
			goto L71;
		    }
		    j = w[k];
		    i__4 = ip[j] + lenl[j];
		    for (ii = ip[j] + len[j] - 1; ii >= i__4; --ii) {
			if (a[ii] >= bval) {
			    goto L60;
			}
			i__ = irn[ii];
			if (iw[i__] != j) {
			    goto L55;
			}
/* Remove entry from matching */
			iw[i__] = 0;
			--num;
			fc[*n - num] = j;
L55:
			;
		    }
L60:
		    lenh[j] = len[j];
/* IP(J)+LEN(J)-1 is last entry in column ge BVAL */
		    len[j] = ii - ip[j] + 1;
/* If LENH(J) = LENL(J), remove J from W */
		    if (lenl[j] == lenh[j]) {
			w[k] = w[wlen];
			--wlen;
		    } else {
			++k;
		    }
/* L70: */
		}
L71:
		if (num < *numx) {
		    goto L81;
		}
/* L80: */
	    }
/* End of dummy loop; this point is never reached */
/* Set mode for next call to MC64U/UD */
L81:
	    mod = 1;
	} else {
/* We do not have a maximum matching in IW. */
	    bmax = bval;
/* BMIN is the bottleneck value of a maximum matching; */
/* for BMAX the matching is not maximum, so BMAX>BMIN */
/*          IF (BMAX .EQ. BMIN) GO TO 99 */
/* Find splitting value BVAL */
	    mc64qd_(&ip[1], &len[1], &lenh[1], &w[1], &wlen, &a[1], &nval, &
		    bval);
	    if (nval == 0 || bval == bmin) {
		goto L99;
	    }
/* Set LEN such that all matrix entries with value ge BVAL are */
/* inside matrix. Store old LEN in LENL. Do this for all columns W(K). */
/* Each step, either K is incremented or WLEN is decremented. */
	    k = 1;
	    i__2 = *n;
	    for (idum3 = 1; idum3 <= i__2; ++idum3) {
		if (k > wlen) {
		    goto L88;
		}
		j = w[k];
		i__3 = ip[j] + lenh[j] - 1;
		for (ii = ip[j] + len[j]; ii <= i__3; ++ii) {
		    if (a[ii] < bval) {
			goto L86;
		    }
/* L85: */
		}
L86:
		lenl[j] = len[j];
		len[j] = ii - ip[j];
		if (lenl[j] == lenh[j]) {
		    w[k] = w[wlen];
		    --wlen;
		} else {
		    ++k;
		}
/* L87: */
	    }
/* End of dummy loop; this point is never reached */
/* Set mode for next call to MC64U/UD */
L88:
	    mod = 0;
	}
	++cnt;
	mc64ud_(&cnt, &mod, n, &irn[1], ne, &ip[1], &len[1], &fc[1], &iw[1], &
		num, numx, &iw4[1], &iw4[*n + 1], &iw4[(*n << 1) + 1], &iw4[*
		n * 3 + 1]);
/* IW contains maximum matching of length NUM */
/* L90: */
    }
/* End of dummy loop; this point is never reached */
/* BMIN is bottleneck value of final matching */
L99:
    if (*numx == *n) {
	goto L1000;
    }
/* The matrix is structurally singular, complete IPERM */
/* W, IW are work arrays */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	w[j] = 0;
/* L300: */
    }
    k = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iperm[i__] == 0) {
	    ++k;
	    iw[k] = i__;
	} else {
	    j = iperm[i__];
	    w[j] = i__;
	}
/* L310: */
    }
    k = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (w[j] != 0) {
	    goto L320;
	}
	++k;
	idum1 = iw[k];
	iperm[idum1] = j;
L320:
	;
    }
L1000:
    return 0;
} /* mc64sd_ */

/* ********************************************************************** */
/* Subroutine */ int_t mc64qd_(int_t *ip, int_t *lenl, int_t *lenh, 
	int_t *w, int_t *wlen, double *a, int_t *nval, double *
	val)
{
    /* System generated locals */
    int_t i__1, i__2, i__3;

    /* Local variables */
    int_t j, k, s;
    double ha;
    int_t ii, pos;
    double split[10];


/* *** Copyright (c) 1999  Council for the Central Laboratory of the */
/*     Research Councils                                             *** */
/* *** Although every effort has been made to ensure robustness and  *** */
/* *** reliability of the subroutines in this MC64 suite, we         *** */
/* *** disclaim any liability arising through the use or misuse of   *** */
/* *** any of the subroutines.                                       *** */
/* *** Any problems?   Contact ... */
/*     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   *** */

/* This routine searches for at most XX different numerical values */
/* in the columns W(1:WLEN). XX>=2. */
/* Each column J is scanned between IP(J)+LENL(J) and IP(J)+LENH(J)-1 */
/* until XX values are found or all columns have been considered. */
/* On output, NVAL is the number of different values that is found */
/* and SPLIT(1:NVAL) contains the values in decreasing order. */
/* If NVAL > 0, the routine returns VAL = SPLIT((NVAL+1)/2). */

/* Scan columns in W(1:WLEN). For each encountered value, if value not */
/* already present in SPLIT(1:NVAL), insert value such that SPLIT */
/* remains sorted by decreasing value. */
/* The sorting is done by straightforward insertion; therefore the use */
/* of this routine should be avoided for large XX (XX < 20). */
    /* Parameter adjustments */
    --a;
    --w;
    --lenh;
    --lenl;
    --ip;

    /* Function Body */
    *nval = 0;
    i__1 = *wlen;
    for (k = 1; k <= i__1; ++k) {
	j = w[k];
	i__2 = ip[j] + lenh[j] - 1;
	for (ii = ip[j] + lenl[j]; ii <= i__2; ++ii) {
	    ha = a[ii];
	    if (*nval == 0) {
		split[0] = ha;
		*nval = 1;
	    } else {
/* Check presence of HA in SPLIT */
		for (s = *nval; s >= 1; --s) {
		    if (split[s - 1] == ha) {
			goto L15;
		    }
		    if (split[s - 1] > ha) {
			pos = s + 1;
			goto L21;
		    }
/* L20: */
		}
		pos = 1;
/* The insertion */
L21:
		i__3 = pos;
		for (s = *nval; s >= i__3; --s) {
		    split[s] = split[s - 1];
/* L22: */
		}
		split[pos - 1] = ha;
		++(*nval);
	    }
/* Exit loop if XX values are found */
	    if (*nval == 10) {
		goto L11;
	    }
L15:
	    ;
	}
/* L10: */
    }
/* Determine VAL */
L11:
    if (*nval > 0) {
	*val = split[(*nval + 1) / 2 - 1];
    }
    return 0;
} /* mc64qd_ */

/* ********************************************************************** */
/* Subroutine */ int_t mc64ud_(int_t *id, int_t *mod, int_t *n, int_t *
	irn, int_t *lirn, int_t *ip, int_t *lenc, int_t *fc, int_t *
	iperm, int_t *num, int_t *numx, int_t *pr, int_t *arp, 
	int_t *cv, int_t *out)
{
    /* System generated locals */
    int_t i__1, i__2, i__3, i__4;

    /* Local variables */
    int_t i__, j, k, j1, ii, kk, id0, id1, in1, in2, nfc, num0, num1, num2, 
	    jord, last;


/* *** Copyright (c) 1999  Council for the Central Laboratory of the */
/*     Research Councils                                             *** */
/* *** Although every effort has been made to ensure robustness and  *** */
/* *** reliability of the subroutines in this MC64 suite, we         *** */
/* *** disclaim any liability arising through the use or misuse of   *** */
/* *** any of the subroutines.                                       *** */
/* *** Any problems?   Contact ... */
/*     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   *** */

/* PR(J) is the previous column to J in the depth first search. */
/*   Array PR is used as workspace in the sorting algorithm. */
/* Elements (I,IPERM(I)) I=1,..,N are entries at the end of the */
/*   algorithm unless N assignments have not been made in which case */
/*   N-NUM pairs (I,IPERM(I)) will not be entries in the matrix. */
/* CV(I) is the most recent loop number (ID+JORD) at which row I */
/*   was visited. */
/* ARP(J) is the number of entries in column J which have been scanned */
/*   when looking for a cheap assignment. */
/* OUT(J) is one less than the number of entries in column J which have */
/*   not been scanned during one pass through the main loop. */
/* NUMX is maximum possible size of matching. */
    /* Parameter adjustments */
    --out;
    --cv;
    --arp;
    --pr;
    --iperm;
    --fc;
    --lenc;
    --ip;
    --irn;

    /* Function Body */
    if (*id == 1) {
/* The first call to MC64U/UD. */
/* Initialize CV and ARP; parameters MOD, NUMX are not accessed */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    cv[i__] = 0;
	    arp[i__] = 0;
/* L5: */
	}
	num1 = *n;
	num2 = *n;
    } else {
/* Not the first call to MC64U/UD. */
/* Re-initialize ARP if entries were deleted since last call to MC64U/UD */
	if (*mod == 1) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		arp[i__] = 0;
/* L8: */
	    }
	}
	num1 = *numx;
	num2 = *n - *numx;
    }
    num0 = *num;
/* NUM0 is size of input matching */
/* NUM1 is maximum possible size of matching */
/* NUM2 is maximum allowed number of unassigned rows/columns */
/* NUM is size of current matching */
/* Quick return if possible */
/*      IF (NUM.EQ.N) GO TO 199 */
/* NFC is number of rows/columns that could not be assigned */
    nfc = 0;
/* Integers ID0+1 to ID0+N are unique numbers for call ID to MC64U/UD, */
/* so 1st call uses 1..N, 2nd call uses N+1..2N, etc */
    id0 = (*id - 1) * *n;
/* Main loop. Each pass round this loop either results in a new */
/* assignment or gives a column with no assignment */
    i__1 = *n;
    for (jord = num0 + 1; jord <= i__1; ++jord) {
/* Each pass uses unique number ID1 */
	id1 = id0 + jord;
/* J is unmatched column */
	j = fc[jord - num0];
	pr[j] = -1;
	i__2 = jord;
	for (k = 1; k <= i__2; ++k) {
/* Look for a cheap assignment */
	    if (arp[j] >= lenc[j]) {
		goto L30;
	    }
	    in1 = ip[j] + arp[j];
	    in2 = ip[j] + lenc[j] - 1;
	    i__3 = in2;
	    for (ii = in1; ii <= i__3; ++ii) {
		i__ = irn[ii];
		if (iperm[i__] == 0) {
		    goto L80;
		}
/* L20: */
	    }
/* No cheap assignment in row */
	    arp[j] = lenc[j];
/* Begin looking for assignment chain starting with row J */
L30:
	    out[j] = lenc[j] - 1;
/* Inner loop.  Extends chain by one or backtracks */
	    i__3 = jord;
	    for (kk = 1; kk <= i__3; ++kk) {
		in1 = out[j];
		if (in1 < 0) {
		    goto L50;
		}
		in2 = ip[j] + lenc[j] - 1;
		in1 = in2 - in1;
/* Forward scan */
		i__4 = in2;
		for (ii = in1; ii <= i__4; ++ii) {
		    i__ = irn[ii];
		    if (cv[i__] == id1) {
			goto L40;
		    }
/* Column J has not yet been accessed during this pass */
		    j1 = j;
		    j = iperm[i__];
		    cv[i__] = id1;
		    pr[j] = j1;
		    out[j1] = in2 - ii - 1;
		    goto L70;
L40:
		    ;
		}
/* Backtracking step. */
L50:
		j1 = pr[j];
		if (j1 == -1) {
/* No augmenting path exists for column J. */
		    ++nfc;
		    fc[nfc] = j;
		    if (nfc > num2) {
/* A matching of maximum size NUM1 is not possible */
			last = jord;
			goto L101;
		    }
		    goto L100;
		}
		j = j1;
/* L60: */
	    }
/* End of dummy loop; this point is never reached */
L70:
	    ;
	}
/* End of dummy loop; this point is never reached */
/* New assignment is made. */
L80:
	iperm[i__] = j;
	arp[j] = ii - ip[j] + 1;
	++(*num);
	i__2 = jord;
	for (k = 1; k <= i__2; ++k) {
	    j = pr[j];
	    if (j == -1) {
		goto L95;
	    }
	    ii = ip[j] + lenc[j] - out[j] - 2;
	    i__ = irn[ii];
	    iperm[i__] = j;
/* L90: */
	}
/* End of dummy loop; this point is never reached */
L95:
	if (*num == num1) {
/* A matching of maximum size NUM1 is found */
	    last = jord;
	    goto L101;
	}

L100:
	;
    }
/* All unassigned columns have been considered */
    last = *n;
/* Now, a transversal is computed or is not possible. */
/* Complete FC before returning. */
L101:
    i__1 = *n;
    for (jord = last + 1; jord <= i__1; ++jord) {
	++nfc;
	fc[nfc] = fc[jord - num0];
/* L110: */
    }
/*  199 RETURN */
    return 0;
} /* mc64ud_ */

/* ********************************************************************** */
/* Subroutine */ int_t mc64wd_(int_t *n, int_t *ne, int_t *ip, int_t *
	irn, double *a, int *iperm, int_t *num, int_t *jperm, 
	int_t *out, int_t *pr, int_t *q, int_t *l, double *u, 
	double *d__)
{
    /* System generated locals */
    int_t i__1, i__2, i__3, c__2 = 2;

    /* Local variables */
    int_t i__, j, k, i0, k0, k1, k2, q0;
    double di;
    int_t ii, jj, kk;
    double vj;
    int_t up;
    double dq0;
    int_t kk1, kk2;
    double csp;
    int_t isp, jsp, low;
    double dmin__, dnew;
    int_t jord, qlen, jdum;
    double rinf;
    extern /* Subroutine */ int_t
	mc64dd_(int_t *, int_t *, int_t *, double *, int_t *, int_t *),
	mc64ed_(int_t *, int_t *, int_t *, double *, int_t *, int_t *),
	mc64fd_(int_t *, int_t *, int_t *, int_t *, double *, int_t *, int_t *);


/* *** Copyright (c) 1999  Council for the Central Laboratory of the */
/*     Research Councils                                             *** */
/* *** Although every effort has been made to ensure robustness and  *** */
/* *** reliability of the subroutines in this MC64 suite, we         *** */
/* *** disclaim any liability arising through the use or misuse of   *** */
/* *** any of the subroutines.                                       *** */
/* *** Any problems?   Contact ... */
/*     Iain Duff (I.Duff@rl.ac.uk) or Jacko Koster (jak@ii.uib.no)   *** */

/* N, NE, IP, IRN are described in MC64A/AD. */
/* A is a REAL (DOUBLE PRECISION in the D-version) array of length NE. */
/*   A(K), K=1..NE, must be set to the value of the entry that */
/*   corresponds to IRN(K). It is not altered. */
/*   All values A(K) must be non-negative. */
/* IPERM is an INT_T array of length N. On exit, it contains the */
/*   weighted matching: IPERM(I) = 0 or row I is matched to column */
/*   IPERM(I). */
/* NUM is an INT_T variable. On exit, it contains the cardinality of */
/*   the matching stored in IPERM. */
/* IW is an INT_T work array of length 5N. */
/* DW is a REAL (DOUBLE PRECISION in the D-version) array of length 2N. */
/*   On exit, U = D(1:N) contains the dual row variable and */
/*   V = D(N+1:2N) contains the dual column variable. If the matrix */
/*   is structurally nonsingular (NUM = N), the following holds: */
/*      U(I)+V(J) <= A(I,J)  if IPERM(I) |= J */
/*      U(I)+V(J)  = A(I,J)  if IPERM(I)  = J */
/*      U(I) = 0  if IPERM(I) = 0 */
/*      V(J) = 0  if there is no I for which IPERM(I) = J */
/* Local variables */
/* Local parameters */
/* External subroutines and/or functions */
/*      EXTERNAL FD05AD,MC64DD,MC64ED,MC64FD */
/*      DOUBLE PRECISION FD05AD */
/* Set RINF to largest positive real number */
/* XSL      RINF = FD05AD(5) */
    /* Parameter adjustments */
    --d__;
    --u;
    --l;
    --q;
    --pr;
    --out;
    --jperm;
    --iperm;
    --ip;
    --a;
    --irn;

    /* Function Body */
    rinf = dmach("Overflow");
/* Initialization */
    *num = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	u[k] = rinf;
	d__[k] = 0.;
	iperm[k] = 0;
	jperm[k] = 0;
	pr[k] = ip[k];
	l[k] = 0;
/* L10: */
    }
/* Initialize U(I) */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = ip[j + 1] - 1;
	for (k = ip[j]; k <= i__2; ++k) {
	    i__ = irn[k];
	    if (a[k] > u[i__]) {
		goto L20;
	    }
	    u[i__] = a[k];
	    iperm[i__] = j;
	    l[i__] = k;
L20:
	    ;
	}
/* L30: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = iperm[i__];
	if (j == 0) {
	    goto L40;
	}
/* Row I is not empty */
	iperm[i__] = 0;
	if (jperm[j] != 0) {
	    goto L40;
	}
/* Assignment of column J to row I */
	++(*num);
	iperm[i__] = j;
	jperm[j] = l[i__];
L40:
	;
    }
    if (*num == *n) {
	goto L1000;
    }
/* Scan unassigned columns; improve assignment */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* JPERM(J) ne 0 iff column J is already assigned */
	if (jperm[j] != 0) {
	    goto L95;
	}
	k1 = ip[j];
	k2 = ip[j + 1] - 1;
/* Continue only if column J is not empty */
	if (k1 > k2) {
	    goto L95;
	}
	vj = rinf;
	i__2 = k2;
	for (k = k1; k <= i__2; ++k) {
	    i__ = irn[k];
	    di = a[k] - u[i__];
	    if (di > vj) {
		goto L50;
	    }
	    if (di < vj || di == rinf) {
		goto L55;
	    }
	    if (iperm[i__] != 0 || iperm[i0] == 0) {
		goto L50;
	    }
L55:
	    vj = di;
	    i0 = i__;
	    k0 = k;
L50:
	    ;
	}
	d__[j] = vj;
	k = k0;
	i__ = i0;
	if (iperm[i__] == 0) {
	    goto L90;
	}
	i__2 = k2;
	for (k = k0; k <= i__2; ++k) {
	    i__ = irn[k];
	    if (a[k] - u[i__] > vj) {
		goto L60;
	    }
	    jj = iperm[i__];
/* Scan remaining part of assigned column JJ */
	    kk1 = pr[jj];
	    kk2 = ip[jj + 1] - 1;
	    if (kk1 > kk2) {
		goto L60;
	    }
	    i__3 = kk2;
	    for (kk = kk1; kk <= i__3; ++kk) {
		ii = irn[kk];
		if (iperm[ii] > 0) {
		    goto L70;
		}
		if (a[kk] - u[ii] <= d__[jj]) {
		    goto L80;
		}
L70:
		;
	    }
	    pr[jj] = kk2 + 1;
L60:
	    ;
	}
	goto L95;
L80:
	jperm[jj] = kk;
	iperm[ii] = jj;
	pr[jj] = kk + 1;
L90:
	++(*num);
	jperm[j] = k;
	iperm[i__] = j;
	pr[j] = k + 1;
L95:
	;
    }
    if (*num == *n) {
	goto L1000;
    }
/* Prepare for main loop */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = rinf;
	l[i__] = 0;
/* L99: */
    }
/* Main loop ... each pass round this loop is similar to Dijkstra's */
/* algorithm for solving the single source shortest path problem */
    i__1 = *n;
    for (jord = 1; jord <= i__1; ++jord) {
	if (jperm[jord] != 0) {
	    goto L100;
	}
/* JORD is next unmatched column */
/* DMIN is the length of shortest path in the tree */
	dmin__ = rinf;
	qlen = 0;
	low = *n + 1;
	up = *n + 1;
/* CSP is the cost of the shortest augmenting path to unassigned row */
/* IRN(ISP). The corresponding column index is JSP. */
	csp = rinf;
/* Build shortest path tree starting from unassigned column (root) JORD */
	j = jord;
	pr[j] = -1;
/* Scan column J */
	i__2 = ip[j + 1] - 1;
	for (k = ip[j]; k <= i__2; ++k) {
	    i__ = irn[k];
	    dnew = a[k] - u[i__];
	    if (dnew >= csp) {
		goto L115;
	    }
	    if (iperm[i__] == 0) {
		csp = dnew;
		isp = k;
		jsp = j;
	    } else {
		if (dnew < dmin__) {
		    dmin__ = dnew;
		}
		d__[i__] = dnew;
		++qlen;
		q[qlen] = k;
	    }
L115:
	    ;
	}
/* Initialize heap Q and Q2 with rows held in Q(1:QLEN) */
	q0 = qlen;
	qlen = 0;
	i__2 = q0;
	for (kk = 1; kk <= i__2; ++kk) {
	    k = q[kk];
	    i__ = irn[k];
	    if (csp <= d__[i__]) {
		d__[i__] = rinf;
		goto L120;
	    }
	    if (d__[i__] <= dmin__) {
		--low;
		q[low] = i__;
		l[i__] = low;
	    } else {
		++qlen;
		l[i__] = qlen;
		mc64dd_(&i__, n, &q[1], &d__[1], &l[1], &c__2);
	    }
/* Update tree */
	    jj = iperm[i__];
	    out[jj] = k;
	    pr[jj] = j;
L120:
	    ;
	}
	i__2 = *num;
	for (jdum = 1; jdum <= i__2; ++jdum) {
/* If Q2 is empty, extract rows from Q */
	    if (low == up) {
		if (qlen == 0) {
		    goto L160;
		}
		i__ = q[1];
		if (d__[i__] >= csp) {
		    goto L160;
		}
		dmin__ = d__[i__];
L152:
		mc64ed_(&qlen, n, &q[1], &d__[1], &l[1], &c__2);
		--low;
		q[low] = i__;
		l[i__] = low;
		if (qlen == 0) {
		    goto L153;
		}
		i__ = q[1];
		if (d__[i__] > dmin__) {
		    goto L153;
		}
		goto L152;
	    }
/* Q0 is row whose distance D(Q0) to the root is smallest */
L153:
	    q0 = q[up - 1];
	    dq0 = d__[q0];
/* Exit loop if path to Q0 is longer than the shortest augmenting path */
	    if (dq0 >= csp) {
		goto L160;
	    }
	    --up;
/* Scan column that matches with row Q0 */
	    j = iperm[q0];
	    vj = dq0 - a[jperm[j]] + u[q0];
	    i__3 = ip[j + 1] - 1;
	    for (k = ip[j]; k <= i__3; ++k) {
		i__ = irn[k];
		if (l[i__] >= up) {
		    goto L155;
		}
/* DNEW is new cost */
		dnew = vj + a[k] - u[i__];
/* Do not update D(I) if DNEW ge cost of shortest path */
		if (dnew >= csp) {
		    goto L155;
		}
		if (iperm[i__] == 0) {
/* Row I is unmatched; update shortest path info */
		    csp = dnew;
		    isp = k;
		    jsp = j;
		} else {
/* Row I is matched; do not update D(I) if DNEW is larger */
		    di = d__[i__];
		    if (di <= dnew) {
			goto L155;
		    }
		    if (l[i__] >= low) {
			goto L155;
		    }
		    d__[i__] = dnew;
		    if (dnew <= dmin__) {
			if (l[i__] != 0) {
			    mc64fd_(&l[i__], &qlen, n, &q[1], &d__[1], &l[1], 
				    &c__2);
			}
			--low;
			q[low] = i__;
			l[i__] = low;
		    } else {
			if (l[i__] == 0) {
			    ++qlen;
			    l[i__] = qlen;
			}
			mc64dd_(&i__, n, &q[1], &d__[1], &l[1], &c__2);
		    }
/* Update tree */
		    jj = iperm[i__];
		    out[jj] = k;
		    pr[jj] = j;
		}
L155:
		;
	    }
/* L150: */
	}
/* If CSP = RINF, no augmenting path is found */
L160:
	if (csp == rinf) {
	    goto L190;
	}
/* Find augmenting path by tracing backward in PR; update IPERM,JPERM */
	++(*num);
	i__ = irn[isp];
	iperm[i__] = jsp;
	jperm[jsp] = isp;
	j = jsp;
	i__2 = *num;
	for (jdum = 1; jdum <= i__2; ++jdum) {
	    jj = pr[j];
	    if (jj == -1) {
		goto L180;
	    }
	    k = out[j];
	    i__ = irn[k];
	    iperm[i__] = jj;
	    jperm[jj] = k;
	    j = jj;
/* L170: */
	}
/* End of dummy loop; this point is never reached */
/* Update U for rows in Q(UP:N) */
L180:
	i__2 = *n;
	for (kk = up; kk <= i__2; ++kk) {
	    i__ = q[kk];
	    u[i__] = u[i__] + d__[i__] - csp;
/* L185: */
	}
L190:
	i__2 = *n;
	for (kk = low; kk <= i__2; ++kk) {
	    i__ = q[kk];
	    d__[i__] = rinf;
	    l[i__] = 0;
/* L191: */
	}
	i__2 = qlen;
	for (kk = 1; kk <= i__2; ++kk) {
	    i__ = q[kk];
	    d__[i__] = rinf;
	    l[i__] = 0;
/* L193: */
	}
L100:
	;
    }
/* End of main loop */
/* Set dual column variable in D(1:N) */
L1000:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = jperm[j];
	if (k != 0) {
	    d__[j] = a[k] - u[irn[k]];
	} else {
	    d__[j] = 0.;
	}
	if (iperm[j] == 0) {
	    u[j] = 0.;
	}
/* L200: */
    }
    if (*num == *n) {
	goto L1100;
    }
/* The matrix is structurally singular, complete IPERM. */
/* JPERM, OUT are work arrays */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jperm[j] = 0;
/* L300: */
    }
    k = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iperm[i__] == 0) {
	    ++k;
	    out[k] = i__;
	} else {
	    j = iperm[i__];
	    jperm[j] = i__;
	}
/* L310: */
    }
    k = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (jperm[j] != 0) {
	    goto L320;
	}
	++k;
	jdum = out[k];
	iperm[jdum] = j;
L320:
	;
    }
L1100:
    return 0;
} /* mc64wd_ */


