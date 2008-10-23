#ifndef __CSR_H__
#define __CSR_H__

#include <set>
#include <vector>
#include <algorithm>
#include <functional>

#include "dense.h"

/*
 * Extract main diagonal of CSR matrix A
 *
 * Input Arguments:
 *   I  n_row         - number of rows in A
 *   I  n_col         - number of columns in A
 *   I  Ap[n_row+1]   - row pointer
 *   I  Aj[nnz(A)]    - column indices
 *   T  Ax[n_col]     - nonzeros
 *
 * Output Arguments:
 *   T  Yx[min(n_row,n_col)] - diagonal entries
 *
 * Note:
 *   Output array Yx must be preallocated
 *
 *   Duplicate entries will be summed.
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + min(n_row,n_col))
 * 
 */
template <class I, class T>
void csr_diagonal(const I n_row,
                  const I n_col, 
	              const I Ap[], 
	              const I Aj[], 
	              const T Ax[],
	                    T Yx[])
{
    const I N = std::min(n_row, n_col);

    for(I i = 0; i < N; i++){
        const I row_start = Ap[i];
        const I row_end   = Ap[i+1];

        T diag = 0;
        for(I jj = row_start; jj < row_end; jj++){
            if (Aj[jj] == i)
                diag += Ax[jj];
        }

        Yx[i] = diag;
    }
}


/*
 * Expand a compressed row pointer into a row array
 *
 * Input Arguments:
 *   I  n_row         - number of rows in A
 *   I  Ap[n_row+1]   - row pointer
 *
 * Output Arguments:
 *   Bi  - row indices
 *
 * Note:
 *   Output array Bi must be preallocated
 *
 * Note: 
 *   Complexity: Linear
 * 
 */
template <class I>
void expandptr(const I n_row,
               const I Ap[], 
                     I Bi[])
{
    for(I i = 0; i < n_row; i++){
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            Bi[jj] = i;
        }
    }
}

/*
 * Scale the rows of a CSR matrix *in place*
 *
 *   A[i,:] *= X[i]
 *
 */
template <class I, class T>
void csr_scale_rows(const I n_row,
                    const I n_col, 
	                const I Ap[], 
	                const I Aj[], 
	                      T Ax[],
	                const T Xx[])
{
    for(I i = 0; i < n_row; i++){
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            Ax[jj] *= Xx[i];
        }
    }
}

/*
 * Scale the columns of a CSR matrix *in place*
 *
 *   A[:,i] *= X[i]
 *
 */
template <class I, class T>
void csr_scale_columns(const I n_row,
                       const I n_col, 
	                   const I Ap[], 
	                   const I Aj[], 
	                         T Ax[],
	                   const T Xx[])
{
    const I nnz = Ap[n_row];
    for(I i = 0; i < nnz; i++){
        Ax[i] *= Xx[Aj[i]];
    }
}


/*
 * Compute the number of occupied RxC blocks in a matrix
 *
 * Input Arguments:
 *   I  n_row         - number of rows in A
 *   I  R             - row blocksize
 *   I  C             - column blocksize
 *   I  Ap[n_row+1]   - row pointer
 *   I  Aj[nnz(A)]    - column indices
 *
 * Output Arguments:
 *   I  num_blocks    - number of blocks
 *
 * Note: 
 *   Complexity: Linear
 * 
 */
template <class I>
I csr_count_blocks(const I n_row,
                   const I n_col,
                   const I R,
                   const I C,
                   const I Ap[], 
                   const I Aj[])
{
    std::vector<I> mask(n_col/C + 1,-1);
    I n_blks = 0;
    for(I i = 0; i < n_row; i++){
        I bi = i/R;
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            I bj = Aj[jj]/C;
            if(mask[bj] != bi){
                mask[bj] = bi;
                n_blks++;
            }
        }
    }
    return n_blks;
}


/*
 * Convert a CSR matrix to BSR format
 *
 * Input Arguments:
 *   I  n_row           - number of rows in A
 *   I  n_col           - number of columns in A
 *   I  R               - row blocksize
 *   I  C               - column blocksize
 *   I  Ap[n_row+1]     - row pointer
 *   I  Aj[nnz(A)]      - column indices
 *   T  Ax[nnz(A)]      - nonzero values
 *
 * Output Arguments:
 *   I  Bp[n_row/R + 1] - block row pointer
 *   I  Bj[nnz(B)]      - column indices
 *   T  Bx[nnz(B)]      - nonzero blocks
 *
 * Note:
 *   Complexity: Linear
 *   Output arrays must be preallocated (with Bx initialized to zero)
 *
 * 
 */

template <class I, class T>
void csr_tobsr(const I n_row,
	           const I n_col, 
	           const I R, 
	           const I C, 
	           const I Ap[], 
	           const I Aj[], 
	           const T Ax[],
	                 I Bp[],
                     I Bj[],
	                 T Bx[])
{
    std::vector<T*> blocks(n_col/C + 1, (T*)0 );

    assert( n_row % R == 0 );
    assert( n_col % C == 0 );

    I n_brow = n_row / R;
    //I n_bcol = n_col / C;

    I RC = R*C;
    I n_blks = 0;

    Bp[0] = 0;

    for(I bi = 0; bi < n_brow; bi++){
        for(I r = 0; r < R; r++){
            I i = R*bi + r;  //row index
            for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
                I j = Aj[jj]; //column index

                I bj = j / C;
                I c  = j % C;
                
                if( blocks[bj] == 0 ){
                    blocks[bj] = Bx + RC*n_blks;
                    Bj[n_blks] = bj;
                    n_blks++;
                }

                *(blocks[bj] + C*r + c) += Ax[jj];
            }
        }

        for(I jj = Ap[R*bi]; jj < Ap[R*(bi+1)]; jj++){
            blocks[Aj[jj] / C] = 0;
        }

        Bp[bi+1] = n_blks;
    }
}



/*
 * Sort CSR column indices inplace
 *
 * Input Arguments:
 *   I  n_row           - number of rows in A
 *   I  Ap[n_row+1]     - row pointer
 *   I  Aj[nnz(A)]      - column indices
 *   T  Ax[nnz(A)]      - nonzeros 
 *
 */
template <class I>
bool csr_has_sorted_indices(const I n_row, 
                            const I Ap[],
                            const I Aj[])
{
  for(I i = 0; i < n_row; i++){
      for(I jj = Ap[i]; jj < Ap[i+1] - 1; jj++){
          if(Aj[jj] > Aj[jj+1]){
              return false;
          }
      }
  }
  return true;
}
template< class T1, class T2 >
bool kv_pair_less(const std::pair<T1,T2>& x, const std::pair<T1,T2>& y){
    return x.first < y.first;
}

template<class I, class T>
void csr_sort_indices(const I n_row,
                      const I Ap[], 
                            I Aj[], 
                            T Ax[])
{
    std::vector< std::pair<I,T> > temp;

    for(I i = 0; i < n_row; i++){
        I row_start = Ap[i];
        I row_end   = Ap[i+1];

        temp.clear();

        for(I jj = row_start; jj < row_end; jj++){
            temp.push_back(std::make_pair(Aj[jj],Ax[jj]));
        }

        std::sort(temp.begin(),temp.end(),kv_pair_less<I,T>);

        for(I jj = row_start, n = 0; jj < row_end; jj++, n++){
            Aj[jj] = temp[n].first;
            Ax[jj] = temp[n].second;
        }
    }    
}




/*
 * Compute B = A for CSR matrix A, CSC matrix B
 *
 * Also, with the appropriate arguments can also be used to:
 *   - compute B = A^t for CSR matrix A, CSR matrix B
 *   - compute B = A^t for CSC matrix A, CSC matrix B
 *   - convert CSC->CSR
 *
 * Input Arguments:
 *   I  n_row         - number of rows in A
 *   I  n_col         - number of columns in A
 *   I  Ap[n_row+1]   - row pointer
 *   I  Aj[nnz(A)]    - column indices
 *   T  Ax[nnz(A)]    - nonzeros
 *
 * Output Arguments:
 *   I  Bp[n_col+1] - column pointer
 *   I  Bj[nnz(A)]  - row indices
 *   T  Bx[nnz(A)]  - nonzeros
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
void csr_tocsc(const I n_row,
	           const I n_col, 
	           const I Ap[], 
	           const I Aj[], 
	           const T Ax[],
	                 I Bp[],
	                 I Bi[],
	                 T Bx[])
{  
    const I nnz = Ap[n_row];

    //compute number of non-zero entries per column of A 
    std::fill(Bp, Bp + n_col, 0);

    for (I n = 0; n < nnz; n++){            
        Bp[Aj[n]]++;
    }

    //cumsum the nnz per column to get Bp[]
    for(I col = 0, cumsum = 0; col < n_col; col++){     
        I temp  = Bp[col];
        Bp[col] = cumsum;
        cumsum += temp;
    }
    Bp[n_col] = nnz; 

    for(I row = 0; row < n_row; row++){
        for(I jj = Ap[row]; jj < Ap[row+1]; jj++){
            I col  = Aj[jj];
            I dest = Bp[col];

            Bi[dest] = row;
            Bx[dest] = Ax[jj];

            Bp[col]++;
        }
    }  

    for(I col = 0, last = 0; col <= n_col; col++){
        I temp  = Bp[col];
        Bp[col] = last;
        last    = temp;
    }
}   



/*
 * Compute B = A for CSR matrix A, ELL matrix B
 *
 * Input Arguments:
 *   I  n_row         - number of rows in A
 *   I  n_col         - number of columns in A
 *   I  Ap[n_row+1]   - row pointer
 *   I  Aj[nnz(A)]    - column indices
 *   T  Ax[nnz(A)]    - nonzeros
 *   I  row_length    - maximum nnz in a row of A
 *
 * Output Arguments:
 *   I  Bj[n_row * row_length]  - column indices
 *   T  Bx[n_row * row_length]  - nonzeros
 *
 * Note:
 *   Output arrays Bj, Bx must be preallocated
 *   Duplicate entries in A are not merged.
 *   Explicit zeros in A are carried over to B.
 *   Rows with fewer than row_length columns are padded with zeros.
 *
 */
template <class I, class T>
void csr_toell(const I n_row,
	           const I n_col, 
	           const I Ap[], 
	           const I Aj[], 
	           const T Ax[],
               const I row_length,
	                 I Bj[],
	                 T Bx[])
{
    const I ell_nnz = row_length * n_row;
    std::fill(Bj, Bj + ell_nnz, 0);
    std::fill(Bx, Bx + ell_nnz, 0);

    for(I i = 0; i < n_row; i++){
        I * Bj_row = Bj + row_length * i;
        T * Bx_row = Bx + row_length * i;
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            *Bj_row = Aj[jj];            
            *Bx_row = Ax[jj];            
            Bj_row++;
            Bx_row++;
        }
    }
}


/*
 * Compute C = A*B for CSR matrices A,B
 *
 *
 * Input Arguments:
 *   I  n_row       - number of rows in A
 *   I  n_col       - number of columns in B (hence C is n_row by n_col)
 *   I  Ap[n_row+1] - row pointer
 *   I  Aj[nnz(A)]  - column indices
 *   T  Ax[nnz(A)]  - nonzeros
 *   I  Bp[?]       - row pointer
 *   I  Bj[nnz(B)]  - column indices
 *   T  Bx[nnz(B)]  - nonzeros
 * Output Arguments:
 *   I  Cp[n_row+1] - row pointer
 *   I  Cj[nnz(C)]  - column indices
 *   T  Cx[nnz(C)]  - nonzeros
 *   
 * Note:
 *   Output arrays Cp, Cj, and Cx must be preallocated
 *   The value of nnz(C) will be stored in Ap[n_row] after the first pass.
 *
 * Note: 
 *   Input:  A and B column indices *are not* assumed to be in sorted order 
 *   Output: C column indices *are not* assumed to be in sorted order
 *           Cx will not contain any zero entries
 *
 *   Complexity: O(n_row*K^2 + max(n_row,n_col)) 
 *                 where K is the maximum nnz in a row of A
 *                 and column of B.
 *
 *
 *  This is an implementation of the SMMP algorithm:
 *
 *    "Sparse Matrix Multiplication Package (SMMP)"
 *      Randolph E. Bank and Craig C. Douglas
 *
 *    http://citeseer.ist.psu.edu/445062.html
 *    http://www.mgnet.org/~douglas/ccd-codes.html
 *
 */


/*
 * Pass 1 computes CSR row pointer for the matrix product C = A * B
 *
 */
template <class I>
void csr_matmat_pass1(const I n_row,
                      const I n_col, 
                      const I Ap[], 
                      const I Aj[], 
                      const I Bp[],
                      const I Bj[],
                            I Cp[])
{
//    // method that uses O(1) temp storage
//    const I hash_size = 1 << 5;
//    I vals[hash_size];
//    I mask[hash_size];
//
//    std::set<I> spill;    
//    
//    for(I i = 0; i < hash_size; i++){
//        vals[i] = -1;
//        mask[i] = -1;
//    }
//
//    Cp[0] = 0;
//
//    I slow_inserts = 0;
//    I total_inserts = 0;
//    I nnz = 0;
//    for(I i = 0; i < n_row; i++){
//        spill.clear();
//        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
//            I j = Aj[jj];
//            for(I kk = Bp[j]; kk < Bp[j+1]; kk++){
//                I k = Bj[kk];
//                // I hash = k & (hash_size - 1);
//                I hash = ((I)2654435761 * k) & (hash_size -1 );
//                total_inserts++;
//                if(mask[hash] != i){
//                    mask[hash] = i;                        
//                    vals[hash] = k;
//                    nnz++;
//                } else {
//                    if (vals[hash] != k){
//                        slow_inserts++;
//                        spill.insert(k);
//                    }
//                }
//            }
//        }       
//        nnz += spill.size();
//        Cp[i+1] = nnz;
//    }
//
//    std::cout << "slow fraction " << ((float) slow_inserts)/ ((float) total_inserts) << std::endl;

    // method that uses O(n) temp storage
    std::vector<I> mask(n_col,-1);
    Cp[0] = 0;

    I nnz = 0;
    for(I i = 0; i < n_row; i++){
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            I j = Aj[jj];
            for(I kk = Bp[j]; kk < Bp[j+1]; kk++){
                I k = Bj[kk];
                if(mask[k] != i){
                    mask[k] = i;                        
                    nnz++;
                }
            }
        }         
        Cp[i+1] = nnz;
    }
}

/*
 * Pass 2 computes CSR entries for matrix C = A*B using the 
 * row pointer Cp[] computed in Pass 1.
 *
 */
template <class I, class T>
void csr_matmat_pass2(const I n_row,
      	              const I n_col, 
      	              const I Ap[], 
      	              const I Aj[], 
      	              const T Ax[],
      	              const I Bp[],
      	              const I Bj[],
      	              const T Bx[],
      	                    I Cp[],
      	                    I Cj[],
      	                    T Cx[])
{
    std::vector<I> next(n_col,-1);
    std::vector<T> sums(n_col, 0);

    I nnz = 0;

    Cp[0] = 0;

    for(I i = 0; i < n_row; i++){
        I head   = -2;
        I length =  0;

        I jj_start = Ap[i];
        I jj_end   = Ap[i+1];
        for(I jj = jj_start; jj < jj_end; jj++){
            I j = Aj[jj];
            T v = Ax[jj];

            I kk_start = Bp[j];
            I kk_end   = Bp[j+1];
            for(I kk = kk_start; kk < kk_end; kk++){
                I k = Bj[kk];

                sums[k] += v*Bx[kk];

                if(next[k] == -1){
                    next[k] = head;                        
                    head = k;
                    length++;
                }
            }
        }         

        for(I jj = 0; jj < length; jj++){

            if(sums[head] != 0){
                Cj[nnz] = head;
                Cx[nnz] = sums[head];
                nnz++;
            }

            I temp = head;                
            head = next[head];

            next[temp] = -1; //clear arrays
            sums[temp] =  0;                              
        }

        Cp[i+1] = nnz;
    }
}





/*
 * Compute C = A (bin_op) B for CSR matrices A,B
 *
 *   bin_op(x,y) - binary operator to apply elementwise
 *
 *   
 * Input Arguments:
 *   I    n_row       - number of rows in A (and B)
 *   I    n_col       - number of columns in A (and B)
 *   I    Ap[n_row+1] - row pointer
 *   I    Aj[nnz(A)]  - column indices
 *   T    Ax[nnz(A)]  - nonzeros
 *   I    Bp[?]       - row pointer
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
 *   Input:  A and B column indices are assumed to be in sorted order 
 *   Output: C column indices are assumed to be in sorted order
 *           Cx will not contain any zero entries
 *
 */
template <class I, class T, class bin_op>
void csr_binop_csr(const I n_row,
                   const I n_col, 
                   const I Ap[], 
                   const I Aj[], 
                   const T Ax[],
                   const I Bp[],
                   const I Bj[],
                   const T Bx[],
                         I Cp[],
                         I Cj[],
                         T Cx[],
                   const bin_op& op)
{
    //Method that works for sorted indices
    // assert( csr_has_sorted_indices(n_row,Ap,Aj) );
    // assert( csr_has_sorted_indices(n_row,Bp,Bj) );

    Cp[0] = 0;
    I nnz = 0;

    for(I i = 0; i < n_row; i++){
        I A_pos = Ap[i];
        I B_pos = Bp[i];
        I A_end = Ap[i+1];
        I B_end = Bp[i+1];

        //while not finished with either row
        while(A_pos < A_end && B_pos < B_end){
            I A_j = Aj[A_pos];
            I B_j = Bj[B_pos];

            if(A_j == B_j){
                T result = op(Ax[A_pos],Bx[B_pos]);
                if(result != 0){
                    Cj[nnz] = A_j;
                    Cx[nnz] = result;
                    nnz++;
                }
                A_pos++;
                B_pos++;
            } else if (A_j < B_j) {
                T result = op(Ax[A_pos],0);
                if (result != 0){
                    Cj[nnz] = A_j;
                    Cx[nnz] = result;
                    nnz++;
                }
                A_pos++; 
            } else {
                //B_j < A_j
                T result = op(0,Bx[B_pos]);
                if (result != 0){
                    Cj[nnz] = B_j;
                    Cx[nnz] = result;
                    nnz++;
                }
                B_pos++;
            }
        }

        //tail
        while(A_pos < A_end){
            T result = op(Ax[A_pos],0);
            if (result != 0){
                Cj[nnz] = Aj[A_pos];
                Cx[nnz] = result;
                nnz++;
            }
            A_pos++;
        }
        while(B_pos < B_end){
            T result = op(0,Bx[B_pos]);
            if (result != 0){
                Cj[nnz] = Bj[B_pos];
                Cx[nnz] = result;
                nnz++;
            }
            B_pos++;
        }
        Cp[i+1] = nnz;
    }
}

/* element-wise binary operations*/
template <class I, class T>
void csr_elmul_csr(const I n_row, const I n_col, 
                   const I Ap[], const I Aj[], const T Ax[],
                   const I Bp[], const I Bj[], const T Bx[],
                         I Cp[],       I Cj[],       T Cx[])
{
    csr_binop_csr(n_row,n_col,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::multiplies<T>());
}

template <class I, class T>
void csr_eldiv_csr(const I n_row, const I n_col, 
                   const I Ap[], const I Aj[], const T Ax[],
                   const I Bp[], const I Bj[], const T Bx[],
                         I Cp[],       I Cj[],       T Cx[])
{
    csr_binop_csr(n_row,n_col,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::divides<T>());
}


template <class I, class T>
void csr_plus_csr(const I n_row, const I n_col, 
                  const I Ap[], const I Aj[], const T Ax[],
                  const I Bp[], const I Bj[], const T Bx[],
                        I Cp[],       I Cj[],       T Cx[])
{
    csr_binop_csr(n_row,n_col,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::plus<T>());
}

template <class I, class T>
void csr_minus_csr(const I n_row, const I n_col, 
                   const I Ap[], const I Aj[], const T Ax[],
                   const I Bp[], const I Bj[], const T Bx[],
                         I Cp[],       I Cj[],       T Cx[])
{
    csr_binop_csr(n_row,n_col,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::minus<T>());
}



/*
 * Sum together duplicate column entries in each row of CSR matrix A
 *
 *   
 * Input Arguments:
 *   I    n_row       - number of rows in A (and B)
 *   I    n_col       - number of columns in A (and B)
 *   I    Ap[n_row+1] - row pointer
 *   I    Aj[nnz(A)]  - column indices
 *   T    Ax[nnz(A)]  - nonzeros
 *   
 * Note:
 *   The column indicies within each row must be in sorted order.
 *   Explicit zeros are retained.
 *   Ap, Aj, and Ax will be modified *inplace*
 *
 */
template <class I, class T>
void csr_sum_duplicates(const I n_row,
                        const I n_col, 
                              I Ap[], 
                              I Aj[], 
                              T Ax[])
{
    I nnz = 0;
    I row_end = 0;
    for(I i = 0; i < n_row; i++){
        I jj = row_end;
        row_end = Ap[i+1];
        while( jj < row_end ){
            I j = Aj[jj];
            T x = Ax[jj];
            jj++;
            while( jj < row_end && Aj[jj] == j ){
                x += Ax[jj];
                jj++;
            }
            Aj[nnz] = j;
            Ax[nnz] = x;
            nnz++;
        }
        Ap[i+1] = nnz;
    }
}

/*
 * Eliminate zero entries from CSR matrix A
 *
 *   
 * Input Arguments:
 *   I    n_row       - number of rows in A (and B)
 *   I    n_col       - number of columns in A (and B)
 *   I    Ap[n_row+1] - row pointer
 *   I    Aj[nnz(A)]  - column indices
 *   T    Ax[nnz(A)]  - nonzeros
 *   
 * Note:
 *   Ap, Aj, and Ax will be modified *inplace*
 *
 */
template <class I, class T>
void csr_eliminate_zeros(const I n_row,
                         const I n_col, 
                               I Ap[], 
                               I Aj[], 
                               T Ax[])
{
    I nnz = 0;
    I row_end = 0;
    for(I i = 0; i < n_row; i++){
        I jj = row_end;
        row_end = Ap[i+1];
        while( jj < row_end ){
            I j = Aj[jj];
            T x = Ax[jj];
            if(x != 0){
                Aj[nnz] = j;
                Ax[nnz] = x;
                nnz++;
            }
            jj++;
        }
        Ap[i+1] = nnz;
    }
}



/*
 * Compute Y += A*X for CSR matrix A and dense vectors X,Y
 *
 *
 * Input Arguments:
 *   I  n_row         - number of rows in A
 *   I  n_col         - number of columns in A
 *   I  Ap[n_row+1]   - row pointer
 *   I  Aj[nnz(A)]    - column indices
 *   T  Ax[nnz(A)]    - nonzeros
 *   T  Xx[n_col]     - input vector
 *
 * Output Arguments:
 *   T  Yx[n_row]     - output vector
 *
 * Note:
 *   Output array Yx must be preallocated
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + n_row)
 * 
 */
template <class I, class T>
void csr_matvec(const I n_row,
	            const I n_col, 
	            const I Ap[], 
	            const I Aj[], 
	            const T Ax[],
	            const T Xx[],
	                  T Yx[])
{
    for(I i = 0; i < n_row; i++){
        T sum = Yx[i];
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            sum += Ax[jj] * Xx[Aj[jj]];
        }
        Yx[i] = sum;
    }
}


/*
 * Compute Y += A*X for CSR matrix A and dense block vectors X,Y
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
 */
template <class I, class T>
void csr_matvecs(const I n_row,
	             const I n_col, 
                 const I n_vecs,
	             const I Ap[], 
	             const I Aj[], 
	             const T Ax[],
	             const T Xx[],
	                   T Yx[])
{
    for(I i = 0; i < n_row; i++){
        T * y = Yx + n_vecs * i;
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            const I j = Aj[jj];
            const T a = Ax[jj];
            const T * x = Xx + n_vecs * j;
            axpy(n_vecs, a, x, y);
        }
    }
}




template<class I, class T>
void get_csr_submatrix(const I n_row,
		               const I n_col,
		               const I Ap[], 
		               const I Aj[], 
		               const T Ax[],
		               const I ir0,
		               const I ir1,
		               const I ic0,
		               const I ic1,
		               std::vector<I>* Bp,
		               std::vector<I>* Bj,
		               std::vector<T>* Bx)
{
    I new_n_row = ir1 - ir0;
    //I new_n_col = ic1 - ic0;  //currently unused
    I new_nnz = 0;
    I kk = 0;

    // Count nonzeros total/per row.
    for(I i = 0; i < new_n_row; i++){
        I row_start = Ap[ir0+i];
        I row_end   = Ap[ir0+i+1];

        for(I jj = row_start; jj < row_end; jj++){
            if ((Aj[jj] >= ic0) && (Aj[jj] < ic1)) {
                new_nnz++;
            }
        }
    }

    // Allocate.
    Bp->resize(new_n_row+1);
    Bj->resize(new_nnz);
    Bx->resize(new_nnz);

    // Assign.
    (*Bp)[0] = 0;
    for(I i = 0; i < new_n_row; i++){
        I row_start = Ap[ir0+i];
        I row_end   = Ap[ir0+i+1];

        for(I jj = row_start; jj < row_end; jj++){
            if ((Aj[jj] >= ic0) && (Aj[jj] < ic1)) {
                (*Bj)[kk] = Aj[jj] - ic0;
                (*Bx)[kk] = Ax[jj];
                kk++;
            }
        }
        (*Bp)[i+1] = kk;
    }
}


/*
 * Count the number of occupied diagonals in CSR matrix A
 *
 * Input Arguments:
 *   I  nnz             - number of nonzeros in A
 *   I  Ai[nnz(A)]      - row indices
 *   I  Aj[nnz(A)]      - column indices
 *
 */
template <class I>
I csr_count_diagonals(const I n_row,
                      const I Ap[],
                      const I Aj[])
{
    std::set<I> diagonals;
    
    for(I i = 0; i < n_row; i++){
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            diagonals.insert(Aj[jj] - i);
        }
    }
    return diagonals.size();
}


#endif
