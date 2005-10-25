

/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
/*
 * File name:	dgscon.c
 * History:     Modified from lapack routines DGECON.
 */
#include <math.h>
#include "util.h"
#include "dsp_defs.h"

void
dgscon(char *norm, SuperMatrix *L, SuperMatrix *U,
       double anorm, double *rcond, int *info)
{
/*
    Purpose   
    =======   

    DGSCON estimates the reciprocal of the condition number of a general 
    real matrix A, in either the 1-norm or the infinity-norm, using   
    the LU factorization computed by DGETRF.   

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
            dgstrf(). Use compressed row subscripts storage for supernodes,
            i.e., L has types: Stype = SC, Dtype = D_, Mtype = TRLU.
 
    U       (input) SuperMatrix*
            The factor U from the factorization Pr*A*Pc=L*U as computed by
            dgstrf(). Use column-wise storage scheme, i.e., U has types:
            Stype = NC, Dtype = D_, Mtype = TRU.
	    
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
    double *work;
    int    *iwork;
    extern int drscl_(int *, double *, double *, int *);

    extern int dlacon_(int *, double *, double *, int *, double *, int *);

    
    /* Test the input parameters. */
    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O");
    if (! onenrm && ! lsame_(norm, "I")) *info = -1;
    else if (L->nrow < 0 || L->nrow != L->ncol ||
             L->Stype != SC || L->Dtype != D_ || L->Mtype != TRLU)
	 *info = -2;
    else if (U->nrow < 0 || U->nrow != U->ncol ||
             U->Stype != NC || U->Dtype != D_ || U->Mtype != TRU) 
	*info = -3;
    if (*info != 0) {
	i = -(*info);
	xerbla_("dgscon", &i);
	return;
    }

    /* Quick return if possible */
    *rcond = 0.;
    if ( L->nrow == 0 || U->nrow == 0) {
	*rcond = 1.;
	return;
    }

    work = doubleCalloc( 3*L->nrow );
    iwork = intMalloc( L->nrow );


    if ( !work || !iwork )
	ABORT("Malloc fails for work arrays in dgscon.");
    
    /* Estimate the norm of inv(A). */
    ainvnm = 0.;
    if ( onenrm ) kase1 = 1;
    else kase1 = 2;
    kase = 0;

    do {
	dlacon_(&L->nrow, &work[L->nrow], &work[0], &iwork[0], &ainvnm, &kase);

	if (kase == 0) break;

	if (kase == kase1) {
	    /* Multiply by inv(L). */
	    sp_dtrsv("Lower", "No transpose", "Unit", L, U, &work[0], info);

	    /* Multiply by inv(U). */
	    sp_dtrsv("Upper", "No transpose", "Non-unit", L, U, &work[0],info);
	    
	} else {

	    /* Multiply by inv(U'). */
	    sp_dtrsv("Upper", "Transpose", "Non-unit", L, U, &work[0], info);

	    /* Multiply by inv(L'). */
	    sp_dtrsv("Lower", "Transpose", "Unit", L, U, &work[0], info);
	    
	}

    } while ( kase != 0 );

    /* Compute the estimate of the reciprocal condition number. */
    if (ainvnm != 0.) *rcond = (1. / ainvnm) / anorm;

    SUPERLU_FREE (work);
    SUPERLU_FREE (iwork);
    return;

} /* dgscon */

