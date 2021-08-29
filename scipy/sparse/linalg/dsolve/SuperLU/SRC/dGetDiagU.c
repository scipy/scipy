/*! @file dGetDiagU.c
 * \brief Extracts main diagonal of matrix
 *
 * <pre> 
 * -- Auxiliary routine in SuperLU (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * Xiaoye S. Li
 * September 11, 2003
 *
 *  Purpose
 * =======
 *
 * GetDiagU extracts the main diagonal of matrix U of the LU factorization.
 *  
 * Arguments
 * =========
 *
 * L      (input) SuperMatrix*
 *        The factor L from the factorization Pr*A*Pc=L*U as computed by
 *        dgstrf(). Use compressed row subscripts storage for supernodes,
 *        i.e., L has types: Stype = SLU_SC, Dtype = SLU_D, Mtype = SLU_TRLU.
 *
 * diagU  (output) double*, dimension (n)
 *        The main diagonal of matrix U.
 *
 * Note
 * ====
 * The diagonal blocks of the L and U matrices are stored in the L
 * data structures.
 * </pre> 
*/
#include <slu_ddefs.h>

void dGetDiagU(SuperMatrix *L, double *diagU)
{
    int_t i, k, nsupers;
    int_t fsupc, nsupr, nsupc, luptr;
    double *dblock, *Lval;
    SCformat *Lstore;

    Lstore = L->Store;
    Lval = Lstore->nzval;
    nsupers = Lstore->nsuper + 1;

    for (k = 0; k < nsupers; ++k) {
      fsupc = L_FST_SUPC(k);
      nsupc = L_FST_SUPC(k+1) - fsupc;
      nsupr = L_SUB_START(fsupc+1) - L_SUB_START(fsupc);
      luptr = L_NZ_START(fsupc);

      dblock = &diagU[fsupc];
      for (i = 0; i < nsupc; ++i) {
	dblock[i] = Lval[luptr];
	luptr += nsupr + 1;
      }
    }
}

