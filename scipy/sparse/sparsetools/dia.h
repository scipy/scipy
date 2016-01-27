#ifndef __DIA_H__
#define __DIA_H__

#include <algorithm>
#include <cstdlib>
#include "dense.h"


/*
 * Compute Y = alpha*A*X + beta*Y for DIA matrix A and dense vectors X,Y
 * and scalars alpha and beta
 *
 *
 * Input Arguments:
 *   I  n_row            - number of rows in A
 *   I  n_col            - number of columns in A
 *   I  n_diags          - number of diagonals
 *   I  L                - length of each diagonal
 *   I  offsets[n_diags] - diagonal offsets
 *   T  diags[n_diags,L] - nonzeros
 *   T  alpha            - scalar factor multiplying A*X
 *   T  Xx[n_col]        - input vector
 *   T  beta             - scalar factor multiplying Y
 *
 * Output Arguments:
 *   T  Yx[n_row]        - output vector
 *
 * Note:
 *   Output array Yx must be preallocated
 *   Negative offsets correspond to lower diagonals
 *   Positive offsets correspond to upper diagonals
 *
 */
template <class I, class T>
void dia_matvec(const I n_row,
                const I n_col,
                const I n_diags,
                const I L,
                const I offsets[],
                const T diags[],
                const T alpha,
                const T Xx[],
                const T beta,
                      T Yx[])
{
    // Find the diagonal closest to 0 offest and do it first.
    I longest_diag = std::max(n_row, n_col) + 1;
    I longest_diag_index = 0;
    for(I i = 0; i < n_diags; i++){
        if(std::abs(offsets[i]) < longest_diag){
            longest_diag = offsets[i];
            longest_diag_index = i;
            if(longest_diag == 0)
                break;
        }
    }
    const I i_start = std::max<I>(0,-longest_diag);
    const I j_start = std::max<I>(0, longest_diag);
    const I j_end   = std::min<I>(std::min<I>(n_row + longest_diag, n_col), L);
    const I N = j_end - j_start;
    const T * diag = diags + (npy_intp)longest_diag_index*L + j_start;
    const T * x = Xx + j_start;
          T * y = Yx + i_start;

    // needs explicit handling of the beta == 0 case since Y can contains NaN.
    if(beta == 0){
        for(I n = 0; n < i_start; n++){
            Yx[n] = 0;
        }
        for(I n = 0; n < N; n++){
            y[n] = alpha * diag[n] * x[n];
        }
        for(I n = i_start + N; n < n_row; n++){
            Yx[n] = 0;
        }
    }
    else if(beta == 1){
        for(I n = 0; n < N; n++){
            y[n] += alpha * diag[n] * x[n];
        }
    }
    else {
        for(I n = 0; n < i_start; n++){
            Yx[n] *= beta;
        }
        for(I n = 0; n < N; n++){
            y[n] = alpha * diag[n] * x[n] + beta*y[n];
        }
        for(I n = i_start + N; n < n_row; n++){
            Yx[n] *= beta;
        }
    }

    for(I i = 0; i < n_diags; i++){
        if(i == longest_diag_index)
            continue;

        const I k = offsets[i];  //diagonal offset

        const I i_start = std::max<I>(0,-k);
        const I j_start = std::max<I>(0, k);
        const I j_end   = std::min<I>(std::min<I>(n_row + k, n_col),L);

        const I N = j_end - j_start;  //number of elements to process

        const T * diag = diags + (npy_intp)i*L + j_start;
        const T * x = Xx + j_start;
              T * y = Yx + i_start;

        for(I n = 0; n < N; n++){
            y[n] += alpha * diag[n] * x[n];
        }
    }
}


#endif
