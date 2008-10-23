

/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#include <math.h>
#include "csp_defs.h"
#include "util.h"

float
cPivotGrowth(int ncols, SuperMatrix *A, int *perm_c, 
             SuperMatrix *L, SuperMatrix *U)
{
/*
 * Purpose
 * =======
 *
 * Compute the reciprocal pivot growth factor of the leading ncols columns
 * of the matrix, using the formula:
 *     min_j ( max_i(abs(A_ij)) / max_i(abs(U_ij)) )
 *
 * Arguments
 * =========
 *
 * ncols    (input) int
 *          The number of columns of matrices A, L and U.
 *
 * A        (input) SuperMatrix*
 *	    Original matrix A, permuted by columns, of dimension
 *          (A->nrow, A->ncol). The type of A can be:
 *          Stype = NC; Dtype = SLU_C; Mtype = GE.
 *
 * L        (output) SuperMatrix*
 *          The factor L from the factorization Pr*A=L*U; use compressed row 
 *          subscripts storage for supernodes, i.e., L has type: 
 *          Stype = SC; Dtype = SLU_C; Mtype = TRLU.
 *
 * U        (output) SuperMatrix*
 *	    The factor U from the factorization Pr*A*Pc=L*U. Use column-wise
 *          storage scheme, i.e., U has types: Stype = NC;
 *          Dtype = SLU_C; Mtype = TRU.
 *
 */
    NCformat *Astore;
    SCformat *Lstore;
    NCformat *Ustore;
    complex  *Aval, *Lval, *Uval;
    int      fsupc, nsupr, luptr, nz_in_U;
    int      i, j, k, oldcol;
    int      *inv_perm_c;
    float   rpg, maxaj, maxuj;
    extern   double slamch_(char *);
    float   smlnum;
    complex   *luval;
    complex   temp_comp;
   
    /* Get machine constants. */
    smlnum = slamch_("S");
    rpg = 1. / smlnum;

    Astore = A->Store;
    Lstore = L->Store;
    Ustore = U->Store;
    Aval = Astore->nzval;
    Lval = Lstore->nzval;
    Uval = Ustore->nzval;
    
    inv_perm_c = (int *) SUPERLU_MALLOC(A->ncol*sizeof(int));
    for (j = 0; j < A->ncol; ++j) inv_perm_c[perm_c[j]] = j;

    for (k = 0; k <= Lstore->nsuper; ++k) {
	fsupc = L_FST_SUPC(k);
	nsupr = L_SUB_START(fsupc+1) - L_SUB_START(fsupc);
	luptr = L_NZ_START(fsupc);
	luval = &Lval[luptr];
	nz_in_U = 1;
	
	for (j = fsupc; j < L_FST_SUPC(k+1) && j < ncols; ++j) {
	    maxaj = 0.;
            oldcol = inv_perm_c[j];
	    for (i = Astore->colptr[oldcol]; i < Astore->colptr[oldcol+1]; ++i)
		maxaj = SUPERLU_MAX( maxaj, slu_c_abs1( &Aval[i]) );
	
	    maxuj = 0.;
	    for (i = Ustore->colptr[j]; i < Ustore->colptr[j+1]; i++)
		maxuj = SUPERLU_MAX( maxuj, slu_c_abs1( &Uval[i]) );
	    
	    /* Supernode */
	    for (i = 0; i < nz_in_U; ++i)
		maxuj = SUPERLU_MAX( maxuj, slu_c_abs1( &luval[i]) );

	    ++nz_in_U;
	    luval += nsupr;

	    if ( maxuj == 0. )
		rpg = SUPERLU_MIN( rpg, 1.);
	    else
		rpg = SUPERLU_MIN( rpg, maxaj / maxuj );
	}
	
	if ( j >= ncols ) break;
    }

    SUPERLU_FREE(inv_perm_c);
    return (rpg);
}
