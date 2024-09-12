#ifndef __COO_H__
#define __COO_H__

#include <algorithm>

/*
 * Compute B = A for COO matrix A, CSR matrix B
 *
 *
 * Input Arguments:
 *   I  n_row      - number of rows in A
 *   I  n_col      - number of columns in A
 *   I  nnz        - number of nonzeros in A
 *   I  Ai[nnz(A)] - row indices
 *   I  Aj[nnz(A)] - column indices
 *   T  Ax[nnz(A)] - nonzeros
 * Output Arguments:
 *   I Bp  - row pointer
 *   I Bj  - column indices
 *   T Bx  - nonzeros
 *
 * Note:
 *   Output arrays Bp, Bj, and Bx must be preallocated
 *
 * Note: 
 *   Input:  row and column indices *are not* assumed to be ordered
 *           
 *   Note: duplicate entries are carried over to the CSR representation
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col))
 * 
 */
template <class I, class T>
void coo_tocsr(const I n_row,
               const I n_col,
               const I nnz,
               const I Ai[],
               const I Aj[],
               const T Ax[],
                     I Bp[],
                     I Bj[],
                     T Bx[])
{
    //compute number of non-zero entries per row of A 
    std::fill(Bp, Bp + n_row, 0);

    for (I n = 0; n < nnz; n++){            
        Bp[Ai[n]]++;
    }

    //cumsum the nnz per row to get Bp[]
    for(I i = 0, cumsum = 0; i < n_row; i++){     
        I temp = Bp[i];
        Bp[i] = cumsum;
        cumsum += temp;
    }
    Bp[n_row] = nnz; 

    //write Aj,Ax into Bj,Bx
    for(I n = 0; n < nnz; n++){
        I row  = Ai[n];
        I dest = Bp[row];

        Bj[dest] = Aj[n];
        Bx[dest] = Ax[n];

        Bp[row]++;
    }

    for(I i = 0, last = 0; i <= n_row; i++){
        I temp = Bp[i];
        Bp[i]  = last;
        last   = temp;
    }

    //now Bp,Bj,Bx form a CSR representation (with possible duplicates)
}

/*
 * Compute B += A for COO matrix A, dense matrix B
 *
 * Input Arguments:
 *   I  n_row           - number of rows in A
 *   I  n_col           - number of columns in A
 *   npy_int64  nnz     - number of nonzeros in A
 *   I  Ai[nnz(A)]      - row indices
 *   I  Aj[nnz(A)]      - column indices
 *   T  Ax[nnz(A)]      - nonzeros 
 *   T  Bx[n_row*n_col] - dense matrix
 *
 */
template <class I, class T>
void coo_todense(const I n_row,
                 const I n_col,
                 const npy_int64 nnz,
                 const I Ai[],
                 const I Aj[],
                 const T Ax[],
                       T Bx[],
                 const int fortran)
{
    if (!fortran) {
        for(npy_int64 n = 0; n < nnz; n++){
            Bx[ (npy_intp)n_col * Ai[n] + Aj[n] ] += Ax[n];
        }
    }
    else {
        for(npy_int64 n = 0; n < nnz; n++){
            Bx[ (npy_intp)n_row * Aj[n] + Ai[n] ] += Ax[n];
        }
    }
}


/*
 * Input Arguments:
 *   I  strides[num_dims]   - strides array
 *   npy_int64  nnz         - number of nonzeros in A
 *   npy_int64  num_dims    - number of dimensions of A
 *   I  Aijk[nnz(A)]        - coords for nonzeros in A
 *   T  Ax[nnz(A)]          - nonzeros in A
 *   T  Bx[]                - dense array
 *
 */
template <class I, class T>
void coo_todense_nd(const I strides[],
                    const npy_int64 nnz,
                    const npy_int64 num_dims,
                    const I Aijk[],
                    const T Ax[],
                          T Bx[],
                    const int fortran)
{

    if (!fortran) {
        for(npy_int64 n = 0; n < nnz; n++) {
            npy_intp index = 0;
            for(npy_int64 d = num_dims - 1; d >= 0; d--) {
                index += Aijk[d * nnz + n] * strides[d];
            }
            Bx[index] += Ax[n];
        }
    }
    else {
        for(npy_int64 n = 0; n < nnz; n++) {
            npy_intp index = 0;
            for(npy_int64 d = 0; d < num_dims; d++) {
                index += Aijk[d * nnz + n] * strides[d];
            }
            Bx[index] += Ax[n];
        }
    }
}


/*
 * Compute Y += A*X for COO matrix A and dense vectors X,Y
 *
 *
 * Input Arguments:
 *   npy_int64  nnz   - number of nonzeros in A
 *   I  Ai[nnz]       - row indices
 *   I  Aj[nnz]       - column indices
 *   T  Ax[nnz]       - nonzero values
 *   T  Xx[n_col]     - input vector
 *
 * Output Arguments:
 *   T  Yx[n_row]     - output vector
 *
 * Notes:
 *   Output array Yx must be preallocated
 *
 *   Complexity: Linear.  Specifically O(nnz(A))
 * 
 */
template <class I, class T>
void coo_matvec(const npy_int64 nnz,
                const I Ai[],
                const I Aj[],
                const T Ax[],
                const T Xx[],
                      T Yx[])
{
    for(npy_int64 n = 0; n < nnz; n++){
        Yx[Ai[n]] += Ax[n] * Xx[Aj[n]];
    }
}


/*
 * Input Arguments:
 *   npy_int64  nnz         - number of nonzeros in A
 *   npy_int64 num_dims     - number of dimensions
 *   I  strides[num_dims-1] - strides array
 *   I  Aijk[nnz * num_dims]- coords for nonzeros in A
 *   T  Ax[nnz]             - nonzero values in A
 *   T  Xx[n_col]           - input vector
 *
 * Output Arguments:
 *   T  Yx[n_row]     - output vector
 * 
 */
template <class I, class T>
void coo_matvec_nd(const npy_int64 nnz,
                const npy_int64 num_dims,
                const I strides[],
                const I Aijk[],
                const T Ax[],
                const T Xx[],
                      T Yx[])
{
    for(npy_int64 n = 0; n < nnz; n++){
        npy_intp index = 0;
        for(npy_int64 d = num_dims - 2; d >= 0; d--) {
            index += Aijk[d * nnz + n] * strides[d];
        }
        Yx[index] += Ax[n] * Xx[Aijk[nnz * (num_dims - 1) + n]];
    }
}



/*
 *
 * Input Arguments:
 *   npy_int64  nnz           - number of nonzeros in A
 *   npy_int64  n_col_B       - number of columns in B
 *   I  Ai[nnz]               - row indices
 *   I  Aj[nnz]               - column indices
 *   T  Ax[nnz]               - nonzero values
 *   T  Bx[n_row_B * n_col_B] - input matrix flattened
 *
 * Output Arguments:
 *   T  Yx[n_row_A * n_col_B] - output matrix flattened
 *
 * Notes:
 *   Output array Yx must be preallocated
 * 
 */
template <class I, class T>
void coo_matmat_dense(const npy_int64 nnz,
                const npy_int64 n_col_B,
                const I Ai[],
                const I Aj[],
                const T Ax[],
                const T Bx[],
                      T Yx[])
{
    for (npy_int64 n = 0; n < nnz; n++) {
        const T x = Ax[n];
        if (x != 0) {
            const npy_int64 dst_offset = Ai[n] * n_col_B;
            const npy_int64 src_offset = Aj[n] * n_col_B;
            for (npy_int64 i = 0; i < n_col_B; i++) {
                Yx[dst_offset + i] += x * Bx[src_offset + i];
            }
        }
    }
}

#endif
