#ifndef __CSC_H__
#define __CSC_H__

#include "csr.h"

/*
 * Compute Y += A*X for CSC matrix A and dense vectors X,Y
 *
 *
 * Input Arguments:
 *   I  n_row         - number of rows in A
 *   I  n_col         - number of columns in A
 *   I  Ap[n_row+1]   - column pointer
 *   I  Ai[nnz(A)]    - row indices
 *   T  Ax[n_col]     - nonzeros
 *   T  Xx[n_col]     - input vector
 *
 * Output Arguments:
 *   T  Yx[n_row]     - output vector
 *
 * Note:
 *   Output array Yx must be preallocated
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + n_col)
 *
 */
template <class I, class T>
void csc_matvec(const I n_row,
                const I n_col,
                const I Ap[],
                const I Ai[],
                const T Ax[],
                const T Xx[],
                      T Yx[])
{
    for(I j = 0; j < n_col; j++){
        I col_start = Ap[j];
        I col_end   = Ap[j+1];

        for(I ii = col_start; ii < col_end; ii++){
            I i    = Ai[ii];
            Yx[i] += Ax[ii] * Xx[j];
        }
    }
}


/*
 * Compute Y += A*X for CSC matrix A and dense block vectors X,Y
 *
 *
 * Input Arguments:
 *   I  n_row            - number of rows in A
 *   I  n_col            - number of columns in A
 *   I  n_vecs           - number of column vectors in X and Y
 *   I  Ap[n_row+1]      - row pointer
 *   I  Aj[nnz(A)]       - column indices
 *   T  Ax[nnz(A)]       - nonzeros
 *   T  Xx[n_col,n_vecs] - input vector
 *
 * Output Arguments:
 *   T  Yx[n_row,n_vecs] - output vector
 *
 * Note:
 *   Output array Yx must be preallocated
 *
 */
template <class I, class T>
void csc_matvecs(const I n_row,
                 const I n_col,
                 const I n_vecs,
                 const I Ap[],
                 const I Ai[],
                 const T Ax[],
                 const T Xx[],
                       T Yx[])
{
    for(I j = 0; j < n_col; j++){
        for(I ii = Ap[j]; ii < Ap[j+1]; ii++){
            const I i = Ai[ii];
            axpy(n_vecs, Ax[ii], Xx + (npy_intp)n_vecs * j, Yx + (npy_intp)n_vecs * i);
        }
    }
}

#endif
