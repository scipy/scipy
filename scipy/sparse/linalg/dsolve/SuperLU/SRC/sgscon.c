
/*! @file sgscon.c
 * \brief Estimates reciprocal of the condition number of a general matrix
 * 
 * <pre>
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 * Modified from lapack routines SGECON.
 * </pre> 
 */

/*
 * File name:	sgscon.c
 * History:     Modified from lapack routines SGECON.
 */
#include <math.h>
#include "slu_sdefs.h"

/*! \brief
 *
 * <pre>
 *   Purpose   
 *   =======   
 *
 *   SGSCON estimates the reciprocal of the condition number of a general 
 *   real matrix A, in either the 1-norm or the infinity-norm, using   
 *   the LU factorization computed by SGETRF.   *
 *
 *   An estimate is obtained for norm(inv(A)), and the reciprocal of the   
 *   condition number is computed as   
 *      RCOND = 1 / ( norm(A) * norm(inv(A)) ).   
 *
 *   See supermatrix.h for the definition of 'SuperMatrix' structure.
 * 
 *   Arguments   
 *   =========   
 *
 *    NORM    (input) char*
 *            Specifies whether the 1-norm condition number or the   
 *            infinity-norm condition number is required:   
 *            = '1' or 'O':  1-norm;   
 *            = 'I':         Infinity-norm.
 *	    
 *    L       (input) SuperMatrix*
 *            The factor L from the factorization Pr*A*Pc=L*U as computed by
 *            sgstrf(). Use compressed row subscripts storage for supernodes,
 *            i.e., L has types: Stype = SLU_SC, Dtype = SLU_S, Mtype = SLU_TRLU.
 * 
 *    U       (input) SuperMatrix*
 *            The factor U from the factorization Pr*A*Pc=L*U as computed by
 *            sgstrf(). Use column-wise storage scheme, i.e., U has types:
 *            Stype = SLU_NC, Dtype = SLU_S, Mtype = SLU_TRU.
 *	    
 *    ANORM   (input) float
 *            If NORM = '1' or 'O', the 1-norm of the original matrix A.   
 *            If NORM = 'I', the infinity-norm of the original matrix A.
 *	    
 *    RCOND   (output) float*
 *           The reciprocal of the condition number of the matrix A,   
 *           computed as RCOND = 1/(norm(A) * norm(inv(A))).
 *	    
 *    INFO    (output) int*
 *           = 0:  successful exit   
 *           < 0:  if INFO = -i, the i-th argument had an illegal value   
 *
 *    ===================================================================== 
 * </pre>
 */

void
sgscon(char *norm, SuperMatrix *L, SuperMatrix *U,
       float anorm, float *rcond, SuperLUStat_t *stat, int *info)
{


    /* Local variables */
    int    kase, kase1, onenrm, i;
    float ainvnm;
    float *work;
    int    *iwork;
    extern int srscl_(int *, float *, float *, int *);

    extern int slacon_(int *, float *, float *, int *, float *, int *);

    
    /* Test the input parameters. */
    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O");
    if (! onenrm && ! lsame_(norm, "I")) *info = -1;
    else if (L->nrow < 0 || L->nrow != L->ncol ||
             L->Stype != SLU_SC || L->Dtype != SLU_S || L->Mtype != SLU_TRLU)
	 *info = -2;
    else if (U->nrow < 0 || U->nrow != U->ncol ||
             U->Stype != SLU_NC || U->Dtype != SLU_S || U->Mtype != SLU_TRU) 
	*info = -3;
    if (*info != 0) {
	i = -(*info);
	xerbla_("sgscon", &i);
	return;
    }

    /* Quick return if possible */
    *rcond = 0.;
    if ( L->nrow == 0 || U->nrow == 0) {
	*rcond = 1.;
	return;
    }

    work = floatCalloc( 3*L->nrow );
    iwork = intMalloc( L->nrow );


    if ( !work || !iwork )
	ABORT("Malloc fails for work arrays in sgscon.");
    
    /* Estimate the norm of inv(A). */
    ainvnm = 0.;
    if ( onenrm ) kase1 = 1;
    else kase1 = 2;
    kase = 0;

    do {
	slacon_(&L->nrow, &work[L->nrow], &work[0], &iwork[0], &ainvnm, &kase);

	if (kase == 0) break;

	if (kase == kase1) {
	    /* Multiply by inv(L). */
	    sp_strsv("L", "No trans", "Unit", L, U, &work[0], stat, info);

	    /* Multiply by inv(U). */
	    sp_strsv("U", "No trans", "Non-unit", L, U, &work[0], stat, info);
	    
	} else {

	    /* Multiply by inv(U'). */
	    sp_strsv("U", "Transpose", "Non-unit", L, U, &work[0], stat, info);

	    /* Multiply by inv(L'). */
	    sp_strsv("L", "Transpose", "Unit", L, U, &work[0], stat, info);
	    
	}

    } while ( kase != 0 );

    /* Compute the estimate of the reciprocal condition number. */
    if (ainvnm != 0.) *rcond = (1. / ainvnm) / anorm;

    SUPERLU_FREE (work);
    SUPERLU_FREE (iwork);
    return;

} /* sgscon */

