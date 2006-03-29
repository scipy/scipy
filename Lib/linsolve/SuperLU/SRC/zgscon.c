
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
/*
 * File name:	zgscon.c
 * History:     Modified from lapack routines ZGECON.
 */
#include <math.h>
#include "zsp_defs.h"

void
zgscon(char *norm, SuperMatrix *L, SuperMatrix *U,
       double anorm, double *rcond, SuperLUStat_t *stat, int *info)
{
/*
    Purpose   
    =======   

    ZGSCON estimates the reciprocal of the condition number of a general 
    real matrix A, in either the 1-norm or the infinity-norm, using   
    the LU factorization computed by ZGETRF.   

    An estimate is obtained for norm(inv(A)), and the reciprocal of the   
    condition number is computed as   
       RCOND = 1 / ( norm(A) * norm(inv(A)) ).   

    See supermatrix.h for the definition of 'SuperMatrix' structure.
 
    Arguments   
    =========   

    NORM    (input) char*
            Specifies whether the 1-norm condition number or the   
            infinity-norm condition number is required:   
            = '1' or 'O':  1-norm;   
            = 'I':         Infinity-norm.
	    
    L       (input) SuperMatrix*
            The factor L from the factorization Pr*A*Pc=L*U as computed by
            zgstrf(). Use compressed row subscripts storage for supernodes,
            i.e., L has types: Stype = SLU_SC, Dtype = SLU_Z, Mtype = SLU_TRLU.
 
    U       (input) SuperMatrix*
            The factor U from the factorization Pr*A*Pc=L*U as computed by
            zgstrf(). Use column-wise storage scheme, i.e., U has types:
            Stype = SLU_NC, Dtype = SLU_Z, Mtype = TRU.
	    
    ANORM   (input) double
            If NORM = '1' or 'O', the 1-norm of the original matrix A.   
            If NORM = 'I', the infinity-norm of the original matrix A.
	    
    RCOND   (output) double*
            The reciprocal of the condition number of the matrix A,   
            computed as RCOND = 1/(norm(A) * norm(inv(A))).
	    
    INFO    (output) int*
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
*/

    /* Local variables */
    int    kase, kase1, onenrm, i;
    double ainvnm;
    doublecomplex *work;
    extern int zrscl_(int *, doublecomplex *, doublecomplex *, int *);

    extern int zlacon_(int *, doublecomplex *, doublecomplex *, double *, int *);

    
    /* Test the input parameters. */
    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O");
    if (! onenrm && ! lsame_(norm, "I")) *info = -1;
    else if (L->nrow < 0 || L->nrow != L->ncol ||
             L->Stype != SLU_SC || L->Dtype != SLU_Z || L->Mtype != SLU_TRLU)
	 *info = -2;
    else if (U->nrow < 0 || U->nrow != U->ncol ||
             U->Stype != SLU_NC || U->Dtype != SLU_Z || U->Mtype != SLU_TRU) 
	*info = -3;
    if (*info != 0) {
	i = -(*info);
	xerbla_("zgscon", &i);
	return;
    }

    /* Quick return if possible */
    *rcond = 0.;
    if ( L->nrow == 0 || U->nrow == 0) {
	*rcond = 1.;
	return;
    }

    work = doublecomplexCalloc( 3*L->nrow );


    if ( !work )
	ABORT("Malloc fails for work arrays in zgscon.");
    
    /* Estimate the norm of inv(A). */
    ainvnm = 0.;
    if ( onenrm ) kase1 = 1;
    else kase1 = 2;
    kase = 0;

    do {
	zlacon_(&L->nrow, &work[L->nrow], &work[0], &ainvnm, &kase);

	if (kase == 0) break;

	if (kase == kase1) {
	    /* Multiply by inv(L). */
	    sp_ztrsv("L", "No trans", "Unit", L, U, &work[0], stat, info);

	    /* Multiply by inv(U). */
	    sp_ztrsv("U", "No trans", "Non-unit", L, U, &work[0], stat, info);
	    
	} else {

	    /* Multiply by inv(U'). */
	    sp_ztrsv("U", "Transpose", "Non-unit", L, U, &work[0], stat, info);

	    /* Multiply by inv(L'). */
	    sp_ztrsv("L", "Transpose", "Unit", L, U, &work[0], stat, info);
	    
	}

    } while ( kase != 0 );

    /* Compute the estimate of the reciprocal condition number. */
    if (ainvnm != 0.) *rcond = (1. / ainvnm) / anorm;

    SUPERLU_FREE (work);
    return;

} /* zgscon */

