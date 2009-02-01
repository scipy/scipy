#ifndef __BSR_H__
#define __BSR_H__

#include <vector>
#include <algorithm>
#include <functional>

#include "csr.h"
#include "dense.h"


template <class I, class T>
void bsr_diagonal(const I n_brow,
                  const I n_bcol, 
                  const I R,
                  const I C,
	              const I Ap[], 
	              const I Aj[], 
	              const T Ax[],
	                    T Yx[])
{
    const I D  = std::min(R*n_brow, C*n_bcol);
    const I RC = R*C;

    for(I i = 0; i < D; i++){
        Yx[i] = 0;
    }

    if ( R == C ){
        //main diagonal with square blocks
        const I end = std::min(n_brow,n_bcol);
        for(I i = 0; i < end; i++){
            for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
                if (i == Aj[jj]){
                    I row = R*i;
                    const T * val = Ax + RC*jj;
                    for(I bi = 0; bi < R; bi++){
                        Yx[row + bi] = *val;
                        val += C + 1;
                    }
                }
            }
        }
    } 
    else 
    {
        //This could be made faster
        const I end = (D/R) + (D % R == 0 ? 0 : 1);
        for(I i = 0; i < end; i++){
            for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
                const I base_row = R*i;
                const I base_col = C*Aj[jj];
                const T * base_val = Ax + RC*jj;

                for(I bi = 0; bi < R; bi++){
                    const I row = base_row + bi;
                    if (row >= D) break;

                    for(I bj = 0; bj < C; bj++){
                        const I col = base_col + bj;
                        if (row == col){
                            Yx[row] = base_val[bi*C + bj];
                        }
                    }
                }
            }
        }
    }
}



/*
 * Scale the rows of a BSR matrix *in place*
 *
 *   A[i,:] *= X[i]
 *
 */
template <class I, class T>
void bsr_scale_rows(const I n_brow,
                    const I n_bcol, 
                    const I R,
                    const I C,
	                const I Ap[], 
	                const I Aj[], 
	                      T Ax[],
	                const T Xx[])
{
    const I RC = R*C;

    for(I i = 0; i < n_brow; i++){
        const T * row_scales = Xx + R*i;

        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
           T * block = Ax + RC*jj;

            for(I bi = 0; bi < R; bi++){
                scal(C, row_scales[bi], block + C*bi);
            }
        }
    }
}

/*
 * Scale the columns of a BSR matrix *in place*
 *
 *   A[:,i] *= X[i]
 *
 */
template <class I, class T>
void bsr_scale_columns(const I n_brow,
                       const I n_bcol, 
                       const I R,
                       const I C,
	                   const I Ap[], 
	                   const I Aj[], 
	                         T Ax[],
	                   const T Xx[])
{
    const I bnnz = Ap[n_brow];
    const I RC  = R*C;
    for(I i = 0; i < bnnz; i++){
        const T * scales = Xx + C*Aj[i] ;
        T * block = Ax + RC*i;
        
        for(I bi = 0; bi < R; bi++){
            for(I bj = 0; bj < C; bj++){
                block[C*bi + bj] *= scales[bj];
            }
        }

    }
}



/*
 * Sort the column block indices of a BSR matrix inplace
 *
 * Input Arguments:
 *   I  n_brow        - number of row blocks in A
 *   I  n_bcol        - number of column blocks in A
 *   I  R             - rows per block
 *   I  C             - columns per block
 *   I  Ap[n_brow+1]  - row pointer
 *   I  Aj[nblk(A)]   - column indices
 *   T  Ax[nnz(A)]    - nonzeros
 *
 */
template <class I, class T>
void bsr_sort_indices(const I n_brow,
	                  const I n_bcol, 
                      const I R,
                      const I C,
	                        I Ap[], 
	                        I Aj[], 
	                        T Ax[])
{  
    if( R == 1 && C == 1 ){
        csr_sort_indices(n_brow, Ap, Aj, Ax);
        return;
    }
    
    
    const I nblks = Ap[n_brow];
    const I RC    = R*C;
    const I nnz   = RC*nblks;

    //compute permutation of blocks using CSR
    std::vector<I> perm(nblks);

    for(I i = 0; i < nblks; i++)
        perm[i] = i;

    csr_sort_indices(n_brow, Ap, Aj, &perm[0]);

    std::vector<T> Ax_copy(nnz);
    std::copy(Ax, Ax + nnz, Ax_copy.begin());

    for(I i = 0; i < nblks; i++){
        const T * input = &Ax_copy[RC * perm[i]];
              T * output = Ax + RC*i;
        std::copy(input, input + RC, output);
    }
}


/*
 * Compute transpose(A) BSR matrix A
 *
 * Input Arguments:
 *   I  n_brow        - number of row blocks in A
 *   I  n_bcol        - number of column blocks in A
 *   I  R             - rows per block
 *   I  C             - columns per block
 *   I  Ap[n_brow+1]  - row pointer
 *   I  Aj[nblk(A)]   - column indices
 *   T  Ax[nnz(A)]    - nonzeros
 *
 * Output Arguments:
 *   I  Bp[n_col+1]   - row pointer
 *   I  Bj[nblk(A)]   - column indices
 *   T  Bx[nnz(A)]    - nonzeros
 *
 * Note:
 *   Output arrays Bp, Bj, Bx must be preallocated
 *
 * Note: 
 *   Input:  column indices *are not* assumed to be in sorted order
 *   Output: row indices *will be* in sorted order
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col))
 * 
 */
template <class I, class T>
void bsr_transpose(const I n_brow,
	               const I n_bcol, 
                   const I R,
                   const I C,
	               const I Ap[], 
	               const I Aj[], 
	               const T Ax[],
	                     I Bp[],
	                     I Bj[],
	                     T Bx[])
{  
    const I nblks = Ap[n_brow];
    const I RC    = R*C;

    //compute permutation of blocks using tranpose(CSR)
    std::vector<I> perm_in (nblks);
    std::vector<I> perm_out(nblks);

    for(I i = 0; i < nblks; i++)
        perm_in[i] = i;

    csr_tocsc(n_brow, n_bcol, Ap, Aj, &perm_in[0], Bp, Bj, &perm_out[0]);

    for(I i = 0; i < nblks; i++){
        const T * Ax_blk = Ax + RC * perm_out[i];
              T * Bx_blk = Bx + RC * i;
        for(I r = 0; r < R; r++){
            for(I c = 0; c < C; c++){
                Bx_blk[c * R + r] = Ax_blk[r * C + c];
            }
        }
    }
}



template <class I, class T>
void bsr_matmat_pass2(const I n_brow,  const I n_bcol, 
                      const I R,       const I C,       const I N,
      	              const I Ap[],    const I Aj[],    const T Ax[],
      	              const I Bp[],    const I Bj[],    const T Bx[],
      	                    I Cp[],          I Cj[],          T Cx[])
{
    assert(R > 0 && C > 0 && N > 0);

    if( R == 1 && N == 1 && C == 1 ){
        // Use CSR for 1x1 blocksize
        csr_matmat_pass2(n_brow, n_bcol, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx);
        return;
    }

    const I RC = R*C;
    const I RN = R*N;
    const I NC = N*C;

    std::fill( Cx, Cx + RC * Cp[n_brow], 0 ); //clear output array
 
    std::vector<I>  next(n_bcol,-1);
    std::vector<T*> mats(n_bcol);
    
    I nnz = 0;
    Cp[0] = 0;

    for(I i = 0; i < n_brow; i++){
        I head   = -2;
        I length =  0;

        I jj_start = Ap[i];
        I jj_end   = Ap[i+1];
        for(I jj = jj_start; jj < jj_end; jj++){
            I j = Aj[jj];

            I kk_start = Bp[j];
            I kk_end   = Bp[j+1];
            for(I kk = kk_start; kk < kk_end; kk++){
                I k = Bj[kk];

                if(next[k] == -1){
                    next[k] = head;                        
                    head = k;
                    Cj[nnz] = k;
                    mats[k] = Cx + RC*nnz;
                    nnz++;
                    length++;
                }

                const T * A = Ax + jj*RN;
                const T * B = Bx + kk*NC;

                gemm(R, C, N, A, B, mats[k]);
            }
        }         

        for(I jj = 0; jj < length; jj++){
            I temp = head;                
            head = next[head];
            next[temp] = -1; //clear arrays
        }

    }
}




template <class I, class T>
bool is_nonzero_block(const T block[], const I blocksize){
    for(I i = 0; i < blocksize; i++){
        if(block[i] != 0){
            return true;
        }
    }
    return false;
}



/*
 * Compute C = A (binary_op) B for BSR matrices that are not
 * necessarily canonical BSR format.  Specifically, this method
 * works even when the input matrices have duplicate and/or
 * unsorted column indices within a given row.
 *
 * Refer to bsr_binop_bsr() for additional information
 *   
 * Note:
 *   Output arrays Cp, Cj, and Cx must be preallocated
 *   If nnz(C) is not known a priori, a conservative bound is:
 *          nnz(C) <= nnz(A) + nnz(B)
 *
 * Note: 
 *   Input:  A and B column indices are not assumed to be in sorted order 
 *   Output: C column indices are not generally in sorted order
 *           C will not contain any duplicate entries or explicit zeros.
 *
 */
template <class I, class T, class bin_op>
void bsr_binop_bsr_general(const I n_brow, const I n_bcol,
                           const I R,      const I C,
                           const I Ap[],  const I Aj[],  const T Ax[],
                           const I Bp[],  const I Bj[],  const T Bx[],
                                 I Cp[],        I Cj[],        T Cx[],
                           const bin_op& op)
{
    //Method that works for duplicate and/or unsorted indices
    const I RC = R*C;

    Cp[0] = 0;
    I nnz = 0;

    std::vector<I>  next(n_bcol,     -1);
    std::vector<T> A_row(n_bcol * RC, 0);   // this approach can be problematic for large R
    std::vector<T> B_row(n_bcol * RC, 0);

    for(I i = 0; i < n_brow; i++){
        I head   = -2;
        I length =  0;

        //add a row of A to A_row
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            I j = Aj[jj];

            for(I n = 0; n < RC; n++)
                A_row[RC*j + n] += Ax[RC*jj + n];

            if(next[j] == -1){
                next[j] = head;                       
                head = j;
                length++;
            }
        }

        //add a row of B to B_row
        for(I jj = Bp[i]; jj < Bp[i+1]; jj++){
            I j = Bj[jj];

            for(I n = 0; n < RC; n++)
                B_row[RC*j + n] += Bx[RC*jj + n];

            if(next[j] == -1){
                next[j] = head;                       
                head = j;
                length++;
            }
        }


        for(I jj = 0; jj < length; jj++){
            // compute op(block_A, block_B)
            for(I n = 0; n < RC; n++)
                Cx[RC * nnz + n] = op(A_row[RC*head + n], B_row[RC*head + n]);

            // advance counter if block is nonzero
            if( is_nonzero_block(Cx + (RC * nnz), RC) )
                Cj[nnz++] = head;

            // clear block_A and block_B values
            for(I n = 0; n < RC; n++){
                A_row[RC*head + n] = 0;
                B_row[RC*head + n] = 0;
            }

            I temp = head;               
            head = next[head];
            next[temp] = -1;
        }
        
        Cp[i + 1] = nnz;
    }
}


/*
 * Compute C = A (binary_op) B for BSR matrices that are in the 
 * canonical BSR format.  Specifically, this method requires that
 * the rows of the input matrices are free of duplicate column indices
 * and that the column indices are in sorted order.
 *
 * Refer to bsr_binop_bsr() for additional information
 *
 * Note: 
 *   Input:  A and B column indices are assumed to be in sorted order 
 *   Output: C column indices will be in sorted order
 *           Cx will not contain any zero entries
 *
 */
template <class I, class T, class bin_op>
void bsr_binop_bsr_canonical(const I n_brow, const I n_bcol, 
                             const I R,      const I C, 
                             const I Ap[],  const I Aj[],  const T Ax[],
                             const I Bp[],  const I Bj[],  const T Bx[],
                                   I Cp[],        I Cj[],        T Cx[],
                             const bin_op& op)
{
    const I RC = R*C;
    T * result = Cx;

    Cp[0] = 0;
    I nnz = 0;

    for(I i = 0; i < n_brow; i++){
        I A_pos = Ap[i];
        I B_pos = Bp[i];
        I A_end = Ap[i+1];
        I B_end = Bp[i+1];

        //while not finished with either row
        while(A_pos < A_end && B_pos < B_end){
            I A_j = Aj[A_pos];
            I B_j = Bj[B_pos];

            if(A_j == B_j){
                for(I n = 0; n < RC; n++){
                    result[n] = op(Ax[RC*A_pos + n], Bx[RC*B_pos + n]);
                }

                if( is_nonzero_block(result,RC) ){
                    Cj[nnz] = A_j;
                    result += RC;
                    nnz++;
                }

                A_pos++; 
                B_pos++;
            } else if (A_j < B_j) {
                for(I n = 0; n < RC; n++){
                    result[n] = op(Ax[RC*A_pos + n], 0);
                }

                if(is_nonzero_block(result,RC)){
                    Cj[nnz] = A_j;
                    result += RC;
                    nnz++;
                }

                A_pos++; 
            } else {
                //B_j < A_j
                for(I n = 0; n < RC; n++){
                    result[n] = op(0, Bx[RC*B_pos + n]);
                }
                if(is_nonzero_block(result,RC)){
                    Cj[nnz] = B_j;
                    result += RC;
                    nnz++;
                }

                B_pos++;
            }
        }

        //tail
        while(A_pos < A_end){
            for(I n = 0; n < RC; n++){
                result[n] = op(Ax[RC*A_pos + n], 0);
            }

            if(is_nonzero_block(result, RC)){
                Cj[nnz] = Aj[A_pos];
                result += RC;
                nnz++;
            }

            A_pos++; 
        }
        while(B_pos < B_end){
            for(I n = 0; n < RC; n++){
                result[n] = op(0,Bx[RC*B_pos + n]);
            }

            if(is_nonzero_block(result, RC)){
                Cj[nnz] = Bj[B_pos];
                result += RC;
                nnz++;
            }

            B_pos++;
        }

        Cp[i+1] = nnz;
    }
}


/*
 * Compute C = A (binary_op) B for CSR matrices A,B where the column 
 * indices with the rows of A and B are known to be sorted.
 *
 *   binary_op(x,y) - binary operator to apply elementwise
 *
 * Input Arguments:
 *   I    n_row       - number of rows in A (and B)
 *   I    n_col       - number of columns in A (and B)
 *   I    Ap[n_row+1] - row pointer
 *   I    Aj[nnz(A)]  - column indices
 *   T    Ax[nnz(A)]  - nonzeros
 *   I    Bp[n_row+1] - row pointer
 *   I    Bj[nnz(B)]  - column indices
 *   T    Bx[nnz(B)]  - nonzeros
 * Output Arguments:
 *   I    Cp[n_row+1] - row pointer
 *   I    Cj[nnz(C)]  - column indices
 *   T    Cx[nnz(C)]  - nonzeros
 *   
 * Note:
 *   Output arrays Cp, Cj, and Cx must be preallocated
 *   If nnz(C) is not known a priori, a conservative bound is:
 *          nnz(C) <= nnz(A) + nnz(B)
 *
 * Note: 
 *   Input:  A and B column indices are not assumed to be in sorted order.
 *   Output: C column indices will be in sorted if both A and B have sorted indices.
 *           Cx will not contain any zero entries
 *
 */
template <class I, class T, class bin_op>
void bsr_binop_bsr(const I n_brow, const I n_bcol, 
                   const I R,     const I C, 
                   const I Ap[],  const I Aj[],  const T Ax[],
                   const I Bp[],  const I Bj[],  const T Bx[],
                         I Cp[],        I Cj[],        T Cx[],
                   const bin_op& op)
{
    assert( R > 0 && C > 0);
    
    if( R == 1 && C == 1 ){
        //use CSR for 1x1 blocksize
        csr_binop_csr(n_brow, n_bcol, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx, op);
    }
    else if ( csr_has_canonical_format(n_brow, Ap, Aj) && csr_has_canonical_format(n_brow, Bp, Bj) ){
        // prefer faster implementation
        bsr_binop_bsr_canonical(n_brow, n_bcol, R, C, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx, op);
    }
    else {
        // slower fallback method
        bsr_binop_bsr_general(n_brow, n_bcol, R, C, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx, op);
    }
}

/* element-wise binary operations */
template <class I, class T>
void bsr_elmul_bsr(const I n_row, const I n_col, const I R, const I C, 
                   const I Ap[], const I Aj[], const T Ax[],
                   const I Bp[], const I Bj[], const T Bx[],
                         I Cp[],       I Cj[],       T Cx[])
{
    bsr_binop_bsr(n_row,n_col,R,C,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::multiplies<T>());
}

template <class I, class T>
void bsr_eldiv_bsr(const I n_row, const I n_col, const I R, const I C,
                   const I Ap[], const I Aj[], const T Ax[],
                   const I Bp[], const I Bj[], const T Bx[],
                         I Cp[],       I Cj[],       T Cx[])
{
    bsr_binop_bsr(n_row,n_col,R,C,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::divides<T>());
}


template <class I, class T>
void bsr_plus_bsr(const I n_row, const I n_col, const I R, const I C, 
                  const I Ap[], const I Aj[], const T Ax[],
                  const I Bp[], const I Bj[], const T Bx[],
                        I Cp[],       I Cj[],       T Cx[])
{
    bsr_binop_bsr(n_row,n_col,R,C,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::plus<T>());
}

template <class I, class T>
void bsr_minus_bsr(const I n_row, const I n_col, const I R, const I C, 
                   const I Ap[], const I Aj[], const T Ax[],
                   const I Bp[], const I Bj[], const T Bx[],
                         I Cp[],       I Cj[],       T Cx[])
{
    bsr_binop_bsr(n_row,n_col,R,C,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::minus<T>());
}





//template <class I, class T>
//void bsr_tocsr(const I n_brow,
//	           const I n_bcol, 
//	           const I R, 
//	           const I C, 
//	           const I Ap[], 
//	           const I Aj[], 
//	           const T Ax[],
//	                 I Bp[],
//                     I Bj[]
//	                 T Bx[])
//{
//    const I RC = R*C;
//
//    for(I brow = 0; brow < n_brow; brow++){
//        I row_size = C * (Ap[brow + 1] - Ap[brow]);
//        for(I r = 0; r < R; r++){
//            Bp[R*brow + r] = RC * Ap[brow] + r * row_size
//        }
//    }
//}


template <class I, class T>
void bsr_matvec(const I n_brow,
	            const I n_bcol, 
	            const I R, 
	            const I C, 
	            const I Ap[], 
	            const I Aj[], 
	            const T Ax[],
	            const T Xx[],
	                  T Yx[])
{
    assert(R > 0 && C > 0);

    if( R == 1 && C == 1 ){
        //use CSR for 1x1 blocksize 
        csr_matvec(n_brow, n_bcol, Ap, Aj, Ax, Xx, Yx);
        return;
    }

    const I RC = R*C;
    for(I i = 0; i < n_brow; i++){
        T * y = Yx + R * i;
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            const I j = Aj[jj];
            const T * A = Ax + RC * jj;
            const T * x = Xx + C * j;
            gemv(R, C, A, x, y); // y += A*x
        }
    }
}


/*
 * Compute Y += A*X for BSR matrix A and dense block vectors X,Y
 *
 *
 * Input Arguments:
 *   I  n_brow              - number of row blocks in A
 *   I  n_bcol              - number of column blocks in A
 *   I  n_vecs              - number of column vectors in X and Y
 *   I  R                   - rows per block
 *   I  C                   - columns per block
 *   I  Ap[n_brow+1]        - row pointer
 *   I  Aj[nblks(A)]        - column indices
 *   T  Ax[nnz(A)]          - nonzeros
 *   T  Xx[C*n_bcol,n_vecs] - input vector
 *
 * Output Arguments:
 *   T  Yx[R*n_brow,n_vecs] - output vector
 *
 */
template <class I, class T>
void bsr_matvecs(const I n_brow,
	             const I n_bcol, 
                 const I n_vecs,
	             const I R, 
	             const I C, 
	             const I Ap[], 
	             const I Aj[], 
	             const T Ax[],
	             const T Xx[],
	                   T Yx[])
{
    assert(R > 0 && C > 0);

    if( R == 1 && C == 1 ){
        //use CSR for 1x1 blocksize 
        csr_matvecs(n_brow, n_bcol, n_vecs, Ap, Aj, Ax, Xx, Yx);
        return;
    }

    const I A_bs = R*C;      //Ax blocksize
    const I Y_bs = n_vecs*R; //Yx blocksize
    const I X_bs = C*n_vecs; //Xx blocksize

    for(I i = 0; i < n_brow; i++){
        T * y = Yx + Y_bs * i;
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            const I j = Aj[jj];
            const T * A = Ax + A_bs * jj;
            const T * x = Xx + X_bs * j;
            gemm(R, n_vecs, C, A, x, y); // y += A*x
        }
    }
}


#endif
