#ifndef SPARSETOOLS_H
#define SPARSETOOLS_H


/*
 * sparsetools.h 
 *   A collection of CSR/CSC/COO/BSR matrix conversion and arithmetic functions.
 *  
 * Original code by Nathan Bell
 *
 */




#include <vector>
#include <algorithm>
#include <functional>
#include "fixed_size.h"

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
        I row_start = Ap[i];
        I row_end   = Ap[i+1];

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
 *   Complexity: Linear.
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
 *   Complexity: Linear.
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
 * Sort CSR column indices inplace
 *
 * Input Arguments:
 *   I  n_row           - number of rows in A
 *   I  n_col           - number of columns in A
 *   I  Ap[n_row+1]     - row pointer
 *   I  Aj[nnz(A)]      - column indices
 *   T  Ax[nnz(A)]      - nonzeros 
 *
 */
template <class I>
bool has_sorted_indices(const I n_row, 
                        const I n_col,
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
                      const I n_col,
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
 *   vec<I> Cp - row pointer
 *   vec<I> Cj - column indices
 *   vec<T> Cx - nonzeros
 *   
 * Note:
 *   Output arrays Cp, Cj, and Cx must be preallocated
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
 * Pass 2 computes CSR entries for C using the row pointer computed in Pass 1
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








template <class I, class T, class bin_op>
void bsr_binop_bsr(const I n_brow, const I n_bcol, 
                   const I R,      const I C, 
                   const I Ap[],   const I Aj[],    const T Ax[],
                   const I Bp[],   const I Bj[],    const T Bx[],
                   std::vector<I>* Cp,
                   std::vector<I>* Cj,
                   std::vector<T>* Cx,
                   const bin_op& op)
{
    Cp->resize(n_brow + 1, 0);

    const I RC = R*C;

    std::vector<I>  next(n_bcol,-1);
    std::vector<T> A_row(n_bcol*RC, 0);
    std::vector<T> B_row(n_bcol*RC, 0);

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
            bool nonzero_block = false;
            for(I n = 0; n < RC; n++){
                T result = op(A_row[RC*head + n],B_row[RC*head + n]);
                A_row[RC*head + n] = result;
                if(result != 0)
                    nonzero_block = true;
            }


            if(nonzero_block){
                Cj->push_back(head);
                for(I n = 0; n < RC; n++){
                    Cx->push_back(A_row[RC*head + n]);
                }
            }

            for(I n = 0; n < RC; n++){
                A_row[RC*head + n] = 0;
                B_row[RC*head + n] = 0;
            }

            I temp = head;                
            head = next[head];
            next[temp] = -1;
        }

        (*Cp)[i+1] = Cj->size();
    }
}

/* element-wise binary operations*/
template <class I, class T>
void bsr_elmul_bsr(const I n_row, const I n_col, const I R, const I C, 
                   const I Ap [], const I Aj [], const T Ax [],
                   const I Bp [], const I Bj [], const T Bx [],
                   std::vector<I>* Cp, std::vector<I>* Cj, std::vector<T>* Cx)
{
    bsr_binop_bsr(n_row,n_col,R,C,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::multiplies<T>());
}

template <class I, class T>
void bsr_eldiv_bsr(const I n_row, const I n_col, const I R, const I C,
                   const I Ap [], const I Aj [], const T Ax [],
                   const I Bp [], const I Bj [], const T Bx [],
                   std::vector<I>* Cp, std::vector<I>* Cj, std::vector<T>* Cx)
{
    bsr_binop_bsr(n_row,n_col,R,C,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::divides<T>());
}


template <class I, class T>
void bsr_plus_bsr(const I n_row, const I n_col, const I R, const I C, 
                  const I Ap [], const I Aj [], const T Ax [],
                  const I Bp [], const I Bj [], const T Bx [],
                  std::vector<I>* Cp, std::vector<I>* Cj, std::vector<T>* Cx)
{
    bsr_binop_bsr(n_row,n_col,R,C,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::plus<T>());
}

template <class I, class T>
void bsr_minus_bsr(const I n_row, const I n_col, const I R, const I C, 
                   const I Ap [], const I Aj [], const T Ax [],
                   const I Bp [], const I Bj [], const T Bx [],
                   std::vector<I>* Cp, std::vector<I>* Cj, std::vector<T>* Cx)
{
    bsr_binop_bsr(n_row,n_col,R,C,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::minus<T>());
}




/*
 * Compute C = A (bin_op) B for CSR matrices A,B
 *
 *   (bin_op) - binary operator to apply elementwise
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
 *   vec<I> Cp  - row pointer
 *   vec<I> Cj  - column indices
 *   vec<T> Cx  - nonzeros
 *   
 * Note:
 *   Output arrays Cp, Cj, and Cx will be allocated within in the method
 *
 * Note: 
 *   Input:  A and B column indices *are not* assumed to be in sorted order 
 *   Output: C column indices *are not* assumed to be in sorted order
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
                   std::vector<I>* Cp,
                   std::vector<I>* Cj,
                   std::vector<T>* Cx,
                   const bin_op& op)
{
//   //Method that works for unsorted indices
//    Cp->resize(n_row + 1, 0);
//    (*Cp)[0] = 0;
//
//    for(I i = 0; i < n_row; i++){
//        I A_pos = Ap[i];
//        I B_pos = Bp[i];
//        I A_end = Ap[i+1];
//        I B_end = Bp[i+1];
//
//        I A_j = Aj[A_pos];
//        I B_j = Bj[B_pos];
//            
//        //while not finished with either row
//        while(A_pos < A_end && B_pos < B_end){
//            if(A_j == B_j){
//                T result = op(Ax[A_pos],Bx[B_pos]);
//                if(result != 0){
//                    Cj->push_back(A_j);
//                    Cx->push_back(result);
//                }
//                A_j = Aj[++A_pos]; 
//                B_j = Bj[++B_pos];
//            } else if (A_j < B_j) {
//                T result = op(Ax[A_pos],0);
//                if (result != 0){
//                    Cj->push_back(A_j);
//                    Cx->push_back(result);
//                }
//                A_j = Aj[++A_pos]; 
//            } else {
//                //B_j < A_j
//                T result = op(0,Bx[B_pos]);
//                if (result != 0){
//                    Cj->push_back(B_j);
//                    Cx->push_back(result);
//                }
//                B_j = Bj[++B_pos];
//            }
//        }
//
//        //tail
//        while(A_pos < A_end){
//            T result = op(Ax[A_pos],0);
//            if (result != 0){
//                Cj->push_back(A_j);
//                Cx->push_back(result);
//            }
//            A_j = Aj[++A_pos]; 
//        }
//        while(B_pos < B_end){
//            T result = op(0,Bx[B_pos]);
//            if (result != 0){
//                Cj->push_back(B_j);
//                Cx->push_back(result);
//            }
//            B_j = Bj[++B_pos];
//        }
//        (*Cp)[i+1] = Cx->size();
//    }


   //Method that works for unsorted indices

    Cp->resize(n_row + 1, 0);

    std::vector<I>  next(n_col,-1);
    std::vector<T> A_row(n_col, 0);
    std::vector<T> B_row(n_col, 0);

    for(I i = 0; i < n_row; i++){
        I head   = -2;
        I length =  0;

        //add a row of A to A_row
        I i_start = Ap[i];
        I i_end   = Ap[i+1];
        for(I jj = i_start; jj < i_end; jj++){
            I j = Aj[jj];

            A_row[j] += Ax[jj];

            if(next[j] == -1){
                next[j] = head;                        
                head = j;
                length++;
            }
        }

        //add a row of B to B_row
        i_start = Bp[i];
        i_end   = Bp[i+1];
        for(I jj = i_start; jj < i_end; jj++){
            I j = Bj[jj];

            B_row[j] += Bx[jj];

            if(next[j] == -1){
                next[j] = head;                        
                head = j;
                length++;
            }
        }


        for(I jj = 0; jj < length; jj++){
            T result = op(A_row[head],B_row[head]);

            if(result != 0){
                Cj->push_back(head);
                Cx->push_back(result);
            }

            I temp = head;                
            head = next[head];

            next[temp] = -1;
            A_row[temp] =  0;                              
            B_row[temp] =  0;
        }

        (*Cp)[i+1] = Cj->size();
    }
}

/* element-wise binary operations*/
template <class I, class T>
void csr_elmul_csr(const I n_row, const I n_col, 
                   const I Ap [], const I Aj [], const T Ax [],
                   const I Bp [], const I Bj [], const T Bx [],
                   std::vector<I>* Cp, std::vector<I>* Cj, std::vector<T>* Cx)
{
    csr_binop_csr(n_row,n_col,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::multiplies<T>());
}

template <class I, class T>
void csr_eldiv_csr(const I n_row, const I n_col, 
                   const I Ap [], const I Aj [], const T Ax [],
                   const I Bp [], const I Bj [], const T Bx [],
                   std::vector<I>* Cp, std::vector<I>* Cj, std::vector<T>* Cx)
{
    csr_binop_csr(n_row,n_col,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::divides<T>());
}


template <class I, class T>
void csr_plus_csr(const I n_row, const I n_col, 
                 const I Ap [], const I Aj [], const T Ax [],
                 const I Bp [], const I Bj [], const T Bx [],
                 std::vector<I>* Cp, std::vector<I>* Cj, std::vector<T>* Cx)
{
    csr_binop_csr(n_row,n_col,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::plus<T>());
}

template <class I, class T>
void csr_minus_csr(const I n_row, const I n_col, 
                   const I Ap [], const I Aj [], const T Ax [],
                   const I Bp [], const I Bj [], const T Bx [],
                   std::vector<I>* Cp, std::vector<I>* Cj, std::vector<T>* Cx)
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
 *   Ap,Aj, and Ax will be modified *inplace*
 *
 */
template <class I, class T>
void csr_sum_duplicates(const I n_row,
                        const I n_col, 
                              I Ap[], 
                              I Aj[], 
                              T Ax[])
{
  std::vector<I>  next(n_col,-1);
  std::vector<T>  sums(n_col, 0);

  I nnz = 0;

  I row_start = 0;
  I row_end   = 0;
  
  for(I i = 0; i < n_row; i++){
    I head = -2;
    
    row_start = row_end; //Ap[i] may have been changed
    row_end   = Ap[i+1]; //Ap[i+1] is safe
    
    for(I jj = row_start; jj < row_end; jj++){
      I j = Aj[jj];

      sums[j] += Ax[jj];
      
      if(next[j] == -1){
	    next[j] = head;                        
	    head    = j;
      }
    }

    while(head != -2){
        I curr = head; //current column
        head   = next[curr];
        
        if(sums[curr] != 0){
            Aj[nnz] = curr;
            Ax[nnz] = sums[curr];
            nnz++;
        }
        
        next[curr] = -1;
        sums[curr] =  0;
    }
    Ap[i+1] = nnz;
  }
}



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
 *   Output: CSR column indices *will be* in sorted order
 *
 *   Note: duplicate entries are carried over to the CSR represention
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
 * Compute Y = A*X for CSR matrix A and dense vectors X,Y
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
        T sum = 0;
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            sum += Ax[jj] * Xx[Aj[jj]];
        }
        Yx[i] = sum;
    }
}



/*
 * Compute Y = A*X for CSC matrix A and dense vectors X,Y
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
    for(I i = 0; i < n_row; i++){
        Yx[i] = 0;
    }

    for(I j = 0; j < n_col; j++){
        I col_start = Ap[j];
        I col_end   = Ap[j+1];

        for(I ii = col_start; ii < col_end; ii++){
            I row  = Ai[ii];
            Yx[row] += Ax[ii] * Xx[j];
        }
    }
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

template <class I, class T, int R, int C>
void bsr_matvec_fixed(const I n_brow,
	                  const I n_bcol, 
	                  const I Ap[], 
	                  const I Aj[], 
	                  const T Ax[],
	                  const T Xx[],
	                        T Yx[])
{
    for(I i = 0; i < n_brow; i++) {
        T r0 = 0;
        T r1 = 0;
        T r2 = 0;
        T r3 = 0;
        T r4 = 0;
        T r5 = 0;
        T r6 = 0;
        T r7 = 0;

        for(I jj = Ap[i]; jj < Ap[i+1]; jj++) {
            I j = Aj[jj];
            const T * base = Ax + jj*(R*C);
            if (R > 0) r0 += dot<C>(base + 0*C, Xx + j*C);
            if (R > 1) r1 += dot<C>(base + 1*C, Xx + j*C);
            if (R > 2) r2 += dot<C>(base + 2*C, Xx + j*C);
            if (R > 3) r3 += dot<C>(base + 3*C, Xx + j*C);
            if (R > 4) r4 += dot<C>(base + 4*C, Xx + j*C);
            if (R > 5) r5 += dot<C>(base + 5*C, Xx + j*C);
            if (R > 6) r6 += dot<C>(base + 6*C, Xx + j*C);
            if (R > 7) r7 += dot<C>(base + 7*C, Xx + j*C);
        }

        if (R > 0) Yx[R*i+0] = r0; 
        if (R > 1) Yx[R*i+1] = r1;
        if (R > 2) Yx[R*i+2] = r2;
        if (R > 3) Yx[R*i+3] = r3;
        if (R > 4) Yx[R*i+4] = r4;
        if (R > 5) Yx[R*i+5] = r5;
        if (R > 6) Yx[R*i+6] = r6;
        if (R > 7) Yx[R*i+7] = r7;
    }
}
#define F(X,Y) bsr_matvec_fixed<I,T,X,Y>

/*
 * Generate the table below with:
 *   out = ''
 *   N = 8
 *   for i in range(N):
 *       out += '{'
 *       for j in range(N-1):
 *           out += ' F(%d,%d),' % (i+1,j+1)
 *       out += ' F(%d,%d) },\n' % (i+1,j+2)
 *   out = out[:-2]
 *
 */

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

    void (*dispatch[8][8])(I,I,const I*,const I*,const T*,const T*,T*) = \
        {
          { F(1,1), F(1,2), F(1,3), F(1,4), F(1,5), F(1,6), F(1,7), F(1,8) },
          { F(2,1), F(2,2), F(2,3), F(2,4), F(2,5), F(2,6), F(2,7), F(2,8) },
          { F(3,1), F(3,2), F(3,3), F(3,4), F(3,5), F(3,6), F(3,7), F(3,8) },
          { F(4,1), F(4,2), F(4,3), F(4,4), F(4,5), F(4,6), F(4,7), F(4,8) },
          { F(5,1), F(5,2), F(5,3), F(5,4), F(5,5), F(5,6), F(5,7), F(5,8) },
          { F(6,1), F(6,2), F(6,3), F(6,4), F(6,5), F(6,6), F(6,7), F(6,8) },
          { F(7,1), F(7,2), F(7,3), F(7,4), F(7,5), F(7,6), F(7,7), F(7,8) },
          { F(8,1), F(8,2), F(8,3), F(8,4), F(8,5), F(8,6), F(8,7), F(8,8) }
        };

    if (R <= 8 && C <= 8){
        dispatch[R-1][C-1](n_brow,n_bcol,Ap,Aj,Ax,Xx,Yx);
        return;
    }


    //otherwise use general method
    for(I i = 0; i < n_brow; i++){
        Yx[i] = 0;
    }

    for(I i = 0; i < n_brow; i++){
        const T * A = Ax + R * C * Ap[i];
              T * y = Yx + R * i;
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            const T * x = Xx + C*Aj[jj];

            //TODO replace this with a proper matvec
            for( I r = 0; r < R; r++ ){
                T sum = 0;
                for( I c = 0; c < C; c++ ){
                    sum += (*A) * x[c];
                    A++;
                }
                y[r] += sum;
            }

        }
    }
}
#undef F










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
 * Derived methods
 */
template <class I, class T>
void csc_diagonal(const I n_row,
                  const I n_col, 
	              const I Ap[], 
	              const I Aj[], 
	              const T Ax[],
	                    T Yx[])
{ csr_diagonal(n_col, n_row, Ap, Aj, Ax, Yx); }


template <class I, class T>
void csc_tocsr(const I n_row,
               const I n_col, 
               const I Ap[], 
               const I Ai[], 
               const T Ax[],
                     I Bp[],
                     I Bj[],
                     T Bx[])
{ csr_tocsc<I,T>(n_col, n_row, Ap, Ai, Ax, Bp, Bj, Bx); }

    
template <class I>
void csc_matmat_pass1(const I n_row,
                      const I n_col, 
                      const I Ap[], 
                      const I Ai[], 
                      const I Bp[],
                      const I Bi[],
                            I Cp[])
{ csr_matmat_pass1(n_col, n_row, Bp, Bi, Ap, Ai, Cp); }
    
template <class I, class T>
void csc_matmat_pass2(const I n_row,
      	              const I n_col, 
      	              const I Ap[], 
      	              const I Ai[], 
      	              const T Ax[],
      	              const I Bp[],
      	              const I Bi[],
      	              const T Bx[],
      	                    I Cp[],
      	                    I Ci[],
      	                    T Cx[])
{ csr_matmat_pass2(n_col, n_row, Bp, Bi, Bx, Ap, Ai, Ax, Cp, Ci, Cx); }


template<class I, class T>
void coo_tocsc(const I n_row,
      	       const I n_col,
      	       const I nnz,
      	       const I Ai[],
      	       const I Aj[],
      	       const T Ax[],
      	             I Bp[],
      	             I Bi[],
      	             T Bx[])
{ coo_tocsr<I,T>(n_col, n_row, nnz, Aj, Ai, Ax, Bp, Bi, Bx); }



template <class I, class T>
void csc_elmul_csc(const I n_row, const I n_col, 
                   const I Ap [], const I Ai [], const T Ax [],
                   const I Bp [], const I Bi [], const T Bx [],
                   std::vector<I>* Cp, std::vector<I>* Ci, std::vector<T>* Cx)
{
    csr_elmul_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}

template <class I, class T>
void csc_eldiv_csc(const I n_row, const I n_col, 
                   const I Ap [], const I Ai [], const T Ax [],
                   const I Bp [], const I Bi [], const T Bx [],
                   std::vector<I>* Cp, std::vector<I>* Ci, std::vector<T>* Cx)
{
    csr_eldiv_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}


template <class I, class T>
void csc_plus_csc(const I n_row, const I n_col, 
                  const I Ap [], const I Ai [], const T Ax [],
                  const I Bp [], const I Bi [], const T Bx [],
                  std::vector<I>* Cp, std::vector<I>* Ci, std::vector<T>* Cx)
{
    csr_plus_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}

template <class I, class T>
void csc_minus_csc(const I n_row, const I n_col, 
                   const I Ap [], const I Ai [], const T Ax [],
                   const I Bp [], const I Bi [], const T Bx [],
                   std::vector<I>* Cp, std::vector<I>* Ci, std::vector<T>* Cx)
{
    csr_minus_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}



template <class I, class T>
void csc_sum_duplicates(const I n_row,
                        const I n_col, 
                              I Ap[], 
                              I Ai[], 
                              T Ax[])
{ csr_sum_duplicates(n_col, n_row, Ap, Ai, Ax); }


template<class I, class T>
void csc_sort_indices(const I n_row,
                      const I n_col,
                      const I Ap[], 
                      I       Ai[], 
                      T       Ax[])
{ csr_sort_indices(n_col, n_row, Ap, Ai, Ax); }






/* 
 * These are sparsetools functions that are not currently used
 * 
 */

///*
// * Compute C = A*B for CSR matrices A,B
// *
// *
// * Input Arguments:
// *   I  n_row       - number of rows in A
// *   I  n_col       - number of columns in B (hence C is n_row by n_col)
// *   I  Ap[n_row+1] - row pointer
// *   I  Aj[nnz(A)]  - column indices
// *   T  Ax[nnz(A)]  - nonzeros
// *   I  Bp[?]       - row pointer
// *   I  Bj[nnz(B)]  - column indices
// *   T  Bx[nnz(B)]  - nonzeros
// * Output Arguments:
// *   vec<I> Cp - row pointer
// *   vec<I> Cj - column indices
// *   vec<T> Cx - nonzeros
// *   
// * Note:
// *   Output arrays Cp, Cj, and Cx will be allocated within in the method
// *
// * Note: 
// *   Input:  A and B column indices *are not* assumed to be in sorted order 
// *   Output: C column indices *are not* assumed to be in sorted order
// *           Cx will not contain any zero entries
// *
// *   Complexity: O(n_row*K^2 + max(n_row,n_col)) 
// *                 where K is the maximum nnz in a row of A
// *                 and column of B.
// *
// *
// *  This implementation closely follows the SMMP algorithm:
// *
// *    "Sparse Matrix Multiplication Package (SMMP)"
// *      Randolph E. Bank and Craig C. Douglas
// *
// *    http://citeseer.ist.psu.edu/445062.html
// *    http://www.mgnet.org/~douglas/ccd-codes.html
// *
// */
//template <class I, class T>
//void csrmucsr(const I n_row,
//      	      const I n_col, 
//      	      const I Ap[], 
//      	      const I Aj[], 
//      	      const T Ax[],
//      	      const I Bp[],
//      	      const I Bj[],
//      	      const T Bx[],
//      	      std::vector<I>* Cp,
//      	      std::vector<I>* Cj,
//      	      std::vector<T>* Cx)
//{
//    Cp->resize(n_row+1,0);
//
//    std::vector<I> next(n_col,-1);
//    std::vector<T> sums(n_col, 0);
//
//    for(I i = 0; i < n_row; i++){
//        I head = -2;
//        I length =  0;
//
//        I jj_start = Ap[i];
//        I jj_end   = Ap[i+1];
//        for(I jj = jj_start; jj < jj_end; jj++){
//            I j = Aj[jj];
//
//            I kk_start = Bp[j];
//            I kk_end   = Bp[j+1];
//            for(I kk = kk_start; kk < kk_end; kk++){
//                I k = Bj[kk];
//
//                sums[k] += Ax[jj]*Bx[kk];
//
//                if(next[k] == -1){
//                    next[k] = head;                        
//                    head = k;
//                    length++;
//                }
//            }
//        }         
//
//        for(I jj = 0; jj < length; jj++){
//            if(sums[head] != 0){
//                Cj->push_back(head);
//                Cx->push_back(sums[head]);
//            }
//
//            I temp = head;                
//            head = next[head];
//
//            next[temp] = -1; //clear arrays
//            sums[temp]  =  0;                              
//        }
//
//        (*Cp)[i+1] = Cx->size();
//    }
//}
//
//
//
//
//
//
//
//
//
//
//
//
///*
// * Compute A = M for CSR matrix A, dense matrix M
// *
// * Input Arguments:
// *   I  n_row           - number of rows in A
// *   I  n_col           - number of columns in A
// *   T  Mx[n_row*n_col] - dense matrix
// *   I  Ap[n_row+1]     - row pointer
// *   I  Aj[nnz(A)]      - column indices
// *   T  Ax[nnz(A)]      - nonzeros 
// *
// * Note:
// *    Output arrays Ap, Aj, and Ax will be allocated within the method
// *
// */
//template <class I, class T>
//void dense_tocsr(const I n_row,
//                 const I n_col,
//                 const T Mx[],
//                 std::vector<I>* Ap,
//                 std::vector<I>* Aj,
//                 std::vector<T>* Ax)
//{
//  const T* x_ptr = Mx;
//
//  Ap->push_back(0);
//  for(I i = 0; i < n_row; i++){
//    for(I j = 0; j < n_col; j++){
//      if(*x_ptr != 0){
//	    Aj->push_back(j);
//	    Ax->push_back(*x_ptr);
//      }
//      x_ptr++;
//    }
//    Ap->push_back(Aj->size());
//  }
//}
//
//
///*
// * Compute M = A for CSR matrix A, dense matrix M
// *
// * Input Arguments:
// *   I  n_row           - number of rows in A
// *   I  n_col           - number of columns in A
// *   I  Ap[n_row+1]     - row pointer
// *   I  Aj[nnz(A)]      - column indices
// *   T  Ax[nnz(A)]      - nonzeros 
// *   T  Mx[n_row*n_col] - dense matrix
// *
// * Note:
// *   Output array Mx is assumed to be allocated and
// *   initialized to 0 by the caller.
// *
// */
//template <class I, class T>
//void csr_todense(const I  n_row,
//                 const I  n_col,
//                 const I  Ap[],
//                 const I  Aj[],
//                 const T  Ax[],
//                       T  Mx[])
//{
//    I row_base = 0;
//    for(I i = 0; i < n_row; i++){
//        I row_start = Ap[i];
//        I row_end   = Ap[i+1];
//        for(I jj = row_start; jj < row_end; jj++){
//            I j = Aj[jj];
//            Mx[row_base + j] = Ax[jj];
//        }	
//        row_base += n_col;
//    }
//}
///*
// * Compute B = A for CSR matrix A, COO matrix B
// *
// * Also, with the appropriate arguments can also be used to:
// *   - convert CSC->COO
// *
// * Input Arguments:
// *   I  n_row         - number of rows in A
// *   I  n_col         - number of columns in A
// *   I  Ap[n_row+1]   - row pointer
// *   I  Aj[nnz(A)]    - column indices
// *   T  Ax[nnz(A)]    - nonzeros
// *
// * Output Arguments:
// *   vec<I> Bi  - row indices
// *   vec<I> Bj  - column indices
// *   vec<T> Bx  - nonzeros
// *
// * Note:
// *   Output arrays Bi, Bj, Bx will be allocated within in the method
// *
// * Note: 
// *   Complexity: Linear.
// * 
// */
//template <class I, class T>
//void csr_tocoo(const I n_row,
//	           const I n_col, 
//               const I Ap[], 
//               const I Aj[], 
//               const T Ax[],
//               std::vector<I>* Bi,
//               std::vector<I>* Bj,
//               std::vector<T>* Bx)
//{
//  I nnz = Ap[n_row];
//  Bi->reserve(nnz);
//  Bi->reserve(nnz);
//  Bx->reserve(nnz);
//  for(I i = 0; i < n_row; i++){
//    I row_start = Ap[i];
//    I row_end   = Ap[i+1];
//    for(I jj = row_start; jj < row_end; jj++){
//      Bi->push_back(i);
//      Bj->push_back(Aj[jj]);
//      Bx->push_back(Ax[jj]);
//    }
//  }
//}
//
//
///*
// * Construct CSC matrix A from diagonals
// *
// * Input Arguments:
// *   I  n_row                            - number of rows in A
// *   I  n_col                            - number of columns in A
// *   I  n_diags                          - number of diagonals
// *   I  diags_indx[n_diags]              - where to place each diagonal 
// *   T  diags[n_diags][min(n_row,n_col)] - diagonals
// *
// * Output Arguments:
// *   vec<I> Ap  - row pointer
// *   vec<I> Aj  - column indices
// *   vec<T> Ax  - nonzeros
// *
// * Note:
// *   Output arrays Ap, Aj, Ax will be allocated within in the method
// *
// * Note: 
// *   Output: row indices are not in sorted order
// *
// *   Complexity: Linear
// * 
// */
//template <class I, class T>
//void spdiags(const I n_row,
//             const I n_col,
//             const I n_diag,
//             const I offsets[],
//             const T diags[],
//             std::vector<I> * Ap,
//             std::vector<I> * Ai,
//             std::vector<T> * Ax)
//{
//    const I diags_length = std::min(n_row,n_col);
//    Ap->push_back(0);
//
//    for(I i = 0; i < n_col; i++){
//        for(I j = 0; j < n_diag; j++){
//            if(offsets[j] <= 0){              //sub-diagonal
//                I row = i - offsets[j];
//                if (row >= n_row){ continue; }
//
//                Ai->push_back(row);
//                Ax->push_back(diags[j*diags_length + i]);
//            } else {                          //super-diagonal
//                I row = i - offsets[j];
//                if (row < 0 || row >= n_row){ continue; }
//                Ai->push_back(row);
//                Ax->push_back(diags[j*diags_length + row]);
//            }
//        }
//        Ap->push_back(Ai->size());
//    }
//}
//

#endif
