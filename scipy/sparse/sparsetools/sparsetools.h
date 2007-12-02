#ifndef SPARSETOOLS_H
#define SPARSETOOLS_H


/*
 * sparsetools.h 
 *   A collection of CSR/CSC/COO matrix conversion and arithmetic functions.
 *  
 * Authors:  
 *    Nathan Bell
 *
 * Revisions:
 *    07/14/2007 - added sum_csr_duplicates
 *    07/12/2007 - added templated function for binary arithmetic ops
 *    01/09/2007 - index type is now templated
 *    01/06/2007 - initial inclusion into SciPy
 *
 */




#include <vector>
#include <algorithm>
#include <functional>

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
 *   vec<T> Yx - diagonal entries
 *
 * Note:
 *   Output array Yx will be allocated within in the method
 *   Duplicate entries will be summed.
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + min(n_row,n_col))
 * 
 */
template <class I, class T>
void extract_csr_diagonal(const I n_row,
                          const I n_col, 
	                      const I Ap[], 
	                      const I Aj[], 
	                      const T Ax[],
	                      std::vector<T>*  Yx)
{
  const I N = std::min(n_row, n_col);
  
  Yx->resize(N);

  for(I i = 0; i < N; i++){
    I row_start = Ap[i];
    I row_end   = Ap[i+1];
    
    T diag = 0;
    for(I jj = row_start; jj < row_end; jj++){
        if (Aj[jj] == i)
            diag += Ax[jj];
    }
    
    (*Yx)[i] = diag;
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
 *   vec<I>  Bp  - row pointer
 *   vec<I>  Bj  - column indices
 *   vec<T>  Bx  - nonzeros
 *
 * Note:
 *   Output arrays Bp,Bj,Bx will be allocated within in the method
 *
 * Note: 
 *   Input:  column indices *are not* assumed to be in sorted order
 *   Output: row indices *will be* in sorted order
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col))
 * 
 */
template <class I, class T>
void csrtocsc(const I n_row,
	          const I n_col, 
	          const I Ap[], 
	          const I Aj[], 
	          const T Ax[],
	          std::vector<I>* Bp,
	          std::vector<I>* Bi,
	          std::vector<T>* Bx)
{  
  I NNZ = Ap[n_row];
  
  Bp->resize(n_col+1);
  Bi->resize(NNZ);
  Bx->resize(NNZ);
 
  std::vector<I> nnz_per_col(n_col,0); //temp array
 
  //compute number of non-zero entries per column of A 
  for (I i = 0; i < NNZ; i++){            
    nnz_per_col[Aj[i]]++;
  }
        
  //cumsum the nnz_per_col to get Bp[]
  for(I i = 0, cumsum = 0; i < n_col; i++){     
    (*Bp)[i]   = cumsum; 
    cumsum += nnz_per_col[i];
    nnz_per_col[i] = 0;              //reset count
  }
  (*Bp)[n_col] = NNZ;
  
  for(I i = 0; i < n_row; i++){
    I row_start = Ap[i];
    I row_end   = Ap[i+1];
    for(I j = row_start; j < row_end; j++){
      I col = Aj[j];
      I k   = (*Bp)[col] + nnz_per_col[col];

      (*Bi)[k] = i;
      (*Bx)[k] = Ax[j];

      nnz_per_col[col]++;
    }
  }  
}   




/*
 * Compute B = A for CSR matrix A, COO matrix B
 *
 * Also, with the appropriate arguments can also be used to:
 *   - convert CSC->COO
 *
 * Input Arguments:
 *   I  n_row         - number of rows in A
 *   I  n_col         - number of columns in A
 *   I  Ap[n_row+1]   - row pointer
 *   I  Aj[nnz(A)]    - column indices
 *   T  Ax[nnz(A)]    - nonzeros
 *
 * Output Arguments:
 *   vec<I> Bi  - row indices
 *   vec<I> Bj  - column indices
 *   vec<T> Bx  - nonzeros
 *
 * Note:
 *   Output arrays Bi,Bj,Bx will be allocated within in the method
 *
 * Note: 
 *   Complexity: Linear.
 * 
 */
template <class I, class T>
void csrtocoo(const I n_row,
	          const I n_col, 
              const I Ap[], 
              const I Aj[], 
              const T Ax[],
              std::vector<I>* Bi,
              std::vector<I>* Bj,
              std::vector<T>* Bx)
{
  I NNZ = Ap[n_row];
  Bi->reserve(NNZ);
  Bi->reserve(NNZ);
  Bx->reserve(NNZ);
  for(I i = 0; i < n_row; i++){
    I row_start = Ap[i];
    I row_end   = Ap[i+1];
    for(I jj = row_start; jj < row_end; jj++){
      Bi->push_back(i);
      Bj->push_back(Aj[jj]);
      Bx->push_back(Ax[jj]);
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
 *   vec<I> Cp - row pointer
 *   vec<I> Cj - column indices
 *   vec<T> Cx - nonzeros
 *   
 * Note:
 *   Output arrays Cp,Cj, and Cx will be allocated within in the method
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
 *  This implementation closely follows the SMMP algorithm:
 *
 *    "Sparse Matrix Multiplication Package (SMMP)"
 *      Randolph E. Bank and Craig C. Douglas
 *
 *    http://citeseer.ist.psu.edu/445062.html
 *    http://www.mgnet.org/~douglas/ccd-codes.html
 *
 */
template <class I, class T>
void csrmucsr(const I n_row,
      	      const I n_col, 
      	      const I Ap[], 
      	      const I Aj[], 
      	      const T Ax[],
      	      const I Bp[],
      	      const I Bj[],
      	      const T Bx[],
      	      std::vector<I>* Cp,
      	      std::vector<I>* Cj,
      	      std::vector<T>* Cx)
{
  Cp->resize(n_row+1,0);
  
  std::vector<I> index(n_col,-1);
  std::vector<T> sums(n_col,0);

  for(I i = 0; i < n_row; i++){
    I istart = -2;
    I length =  0;
    
    for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
      I j = Aj[jj];
      for(I kk = Bp[j]; kk < Bp[j+1]; kk++){
    	I k = Bj[kk];
        
	    sums[k] += Ax[jj]*Bx[kk];
        
	    if(index[k] == -1){
	        index[k] = istart;                        
	        istart = k;
	        length++;
	    }
      }
    }         

    for(I jj = 0; jj < length; jj++){
      if(sums[istart] != 0){
	    Cj->push_back(istart);
	    Cx->push_back(sums[istart]);
      }
	
      I temp = istart;                
      istart = index[istart];
      
      index[temp] = -1; //clear arrays
      sums[temp]  =  0;                              
    }
    
    (*Cp)[i+1] = Cx->size();
  }
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
 *   Output arrays Cp,Cj, and Cx will be allocated within in the method
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
  Cp->resize(n_row+1,0);
  
  std::vector<I>   index(n_col,-1);
  std::vector<T> A_row(n_col,0);
  std::vector<T> B_row(n_col,0);

  for(I i = 0; i < n_row; i++){
    I istart = -2;
    I length =  0;
    
    //add a row of A to A_row
    for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
      I j = Aj[jj];

      A_row[j] += Ax[jj];
      
      if(index[j] == -1){
	    index[j] = istart;                        
	    istart = j;
	    length++;
      }
    }
    
    //add a row of B to B_row
    for(I jj = Bp[i]; jj < Bp[i+1]; jj++){
      I j = Bj[jj];

      B_row[j] += Bx[jj];

      if(index[j] == -1){
          index[j] = istart;                        
	      istart = j;
	      length++;
      }
    }


    for(I jj = 0; jj < length; jj++){
      T result = op(A_row[istart],B_row[istart]);
      
      if(result != 0){
	    Cj->push_back(istart);
	    Cx->push_back(result);
      }
      
      I temp = istart;                
      istart = index[istart];
      
      index[temp] = -1;
      A_row[temp] =  0;                              
      B_row[temp] =  0;
    }
    
    (*Cp)[i+1] = Cx->size();
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
void sum_csr_duplicates(const I n_row,
                        const I n_col, 
                              I Ap[], 
                              I Aj[], 
                              T Ax[])
{
  std::vector<I>  next(n_col,-1);
  std::vector<T>  sums(n_col, 0);

  I NNZ = 0;

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
            Aj[NNZ] = curr;
            Ax[NNZ] = sums[curr];
            NNZ++;
        }
        
        next[curr] = -1;
        sums[curr] =  0;
    }
    Ap[i+1] = NNZ;
  }
}



/*
 * Compute B = A for COO matrix A, CSR matrix B
 *
 *
 * Input Arguments:
 *   I  n_row         - number of rows in A
 *   I  n_col         - number of columns in A
 *   I  Ai[nnz(A)]    - row indices
 *   I  Aj[nnz(A)]    - column indices
 *   T  Ax[nnz(A)]    - nonzeros
 * Output Arguments:
 *   vec<I> Bp        - row pointer
 *   vec<I> Bj        - column indices
 *   vec<T> Bx        - nonzeros
 *
 * Note:
 *   Output arrays Bp,Bj,Bx will be allocated within in the method
 *
 * Note: 
 *   Input:  row and column indices *are not* assumed to be ordered
 *           duplicate (i,j) entries will be summed together
 *
 *   Output: CSR column indices *will be* in sorted order
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col))
 * 
 */
template <class I, class T>
void cootocsr(const I n_row,
              const I n_col,
              const I NNZ,
              const I Ai[],
              const I Aj[],
              const T Ax[],
              std::vector<I>* Bp,
              std::vector<I>* Bj,
              std::vector<T>* Bx)
{
  Bp->resize(n_row+1,0);
  Bj->resize(NNZ);
  Bx->resize(NNZ);
  
  std::vector<I> nnz_per_row(n_row,0); //temp array

  //compute nnz per row, then compute Bp
  for(I i = 0; i < NNZ; i++){
    nnz_per_row[Ai[i]]++;
  }
  for(I i = 0, cumsum = 0; i < n_row; i++){
    (*Bp)[i]          = cumsum;
    cumsum        += nnz_per_row[i];
    nnz_per_row[i] = 0; //reset count
  }
  (*Bp)[n_row] = NNZ;


  //write Aj,Ax Io tempBj,tempBx
  for(I i = 0; i < NNZ; i++){
    I row = Ai[i];
    I n   = (*Bp)[row] + nnz_per_row[row];

    (*Bj)[n] = Aj[i];
    (*Bx)[n] = Ax[i];

    nnz_per_row[row]++;
  }
  //now tempBp,tempBj,tempBx form a CSR representation (with duplicates)

  sum_csr_duplicates(n_row,n_col,&(*Bp)[0],&(*Bj)[0],&(*Bx)[0]);

  //trim unused space at the end
  Bj->resize(Bp->back());
  Bx->resize(Bp->back());
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
 *   T  Ax[n_col]     - nonzeros
 *   T  Xx[n_col]     - nonzeros
 *
 * Output Arguments:
 *   vec<T> Yx - nonzeros 
 *
 * Note:
 *   Output array Xx will be allocated within in the method
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col))
 * 
 */
template <class I, class T>
void csrmux(const I n_row,
	        const I n_col, 
	        const I Ap[], 
	        const I Aj[], 
	        const T Ax[],
	        const T Xx[],
	        std::vector<T>*  Yx)
{
  Yx->resize(n_row);

  for(I i = 0; i < n_row; i++){
    I row_start = Ap[i];
    I row_end   = Ap[i+1];
    
    T sum = 0;
    for(I jj = row_start; jj < row_end; jj++){
      sum += Ax[jj] * Xx[Aj[jj]];
    }
    (*Yx)[i] = sum;
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
 *   T  Xx[n_col]     - nonzeros
 *
 * Output Arguments:
 *   vec<T> Yx - nonzeros 
 *
 * Note:
 *   Output arrays Xx will be allocated within in the method
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col))
 * 
 */
template <class I, class T>
void cscmux(const I n_row,
	        const I n_col, 
	        const I Ap[], 
	        const I Ai[], 
	        const T Ax[],
	        const T Xx[],
	        std::vector<T>*  Yx)
{
  Yx->resize(n_row,0);
  
  for(I j = 0; j < n_col; j++){
    I col_start = Ap[j];
    I col_end   = Ap[j+1];
    
    for(I ii = col_start; ii < col_end; ii++){
      I row  = Ai[ii];
      (*Yx)[row] += Ax[ii] * Xx[j];
    }
  }
}




/*
 * Construct CSC matrix A from diagonals
 *
 * Input Arguments:
 *   I  n_row                            - number of rows in A
 *   I  n_col                            - number of columns in A
 *   I  n_diags                          - number of diagonals
 *   I  diags_indx[n_diags]              - where to place each diagonal 
 *   T  diags[n_diags][min(n_row,n_col)] - diagonals
 *
 * Output Arguments:
 *   vec<I> Ap  - row pointer
 *   vec<I> Aj  - column indices
 *   vec<T> Ax  - nonzeros
 *
 * Note:
 *   Output arrays Ap,Aj,Ax will be allocated within in the method
 *
 * Note: 
 *   Output: row indices are not in sorted order
 *
 *   Complexity: Linear
 * 
 */
template <class I, class T>
void spdiags(const I n_row,
             const I n_col,
             const I n_diag,
             const I offsets[],
             const T diags[],
             std::vector<I> * Ap,
             std::vector<I> * Ai,
             std::vector<T> * Ax)
{
  const I diags_length = std::min(n_row,n_col);
  Ap->push_back(0);

  for(I i = 0; i < n_col; i++){
    for(I j = 0; j < n_diag; j++){
      if(offsets[j] <= 0){              //sub-diagonal
    	I row = i - offsets[j];
	    if (row >= n_row){ continue; }
        
	    Ai->push_back(row);
        Ax->push_back(diags[j*diags_length + i]);
      } else {                          //super-diagonal
	    I row = i - offsets[j];
	    if (row < 0 || row >= n_row){ continue; }
        Ai->push_back(row);
        Ax->push_back(diags[j*diags_length + row]);
      }
    }
    Ap->push_back(Ai->size());
  }
}



/*
 * Compute M = A for CSR matrix A, dense matrix M
 *
 * Input Arguments:
 *   I  n_row           - number of rows in A
 *   I  n_col           - number of columns in A
 *   I  Ap[n_row+1]     - row pointer
 *   I  Aj[nnz(A)]      - column indices
 *   T    Ax[nnz(A)]      - nonzeros 
 *   T    Mx[n_row*n_col] - dense matrix
 *
 * Note:
 *   Output array Mx is assumed to be allocated and
 *   initialized to 0 by the caller.
 *
 */
template <class I, class T>
void csrtodense(const I  n_row,
                const I  n_col,
                const I  Ap[],
                const I  Aj[],
                const T  Ax[],
                      T  Mx[])
{
  I row_base = 0;
  for(I i = 0; i < n_row; i++){
    I row_start = Ap[i];
    I row_end   = Ap[i+1];
    for(I jj = row_start; jj < row_end; jj++){
      I j = Aj[jj];
      Mx[row_base + j] = Ax[jj];
    }	
    row_base += n_col;
  }
}



/*
 * Compute A = M for CSR matrix A, dense matrix M
 *
 * Input Arguments:
 *   I  n_row           - number of rows in A
 *   I  n_col           - number of columns in A
 *   T  Mx[n_row*n_col] - dense matrix
 *   I  Ap[n_row+1]     - row pointer
 *   I  Aj[nnz(A)]      - column indices
 *   T  Ax[nnz(A)]      - nonzeros 
 *
 * Note:
 *    Output arrays Ap,Aj,Ax will be allocated within the method
 *
 */
template <class I, class T>
void densetocsr(const I n_row,
                const I n_col,
                const T Mx[],
                std::vector<I>* Ap,
                std::vector<I>* Aj,
                std::vector<T>* Ax)
{
  const T* x_ptr = Mx;

  Ap->push_back(0);
  for(I i = 0; i < n_row; i++){
    for(I j = 0; j < n_col; j++){
      if(*x_ptr != 0){
	    Aj->push_back(j);
	    Ax->push_back(*x_ptr);
      }
      x_ptr++;
    }
    Ap->push_back(Aj->size());
  }
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
template< class T1, class T2 >
bool kv_pair_less(const std::pair<T1,T2>& x, const std::pair<T1,T2>& y){
    return x.first < y.first;
}

template<class I, class T>
void sort_csr_indices(const I n_row,
                      const I n_col,
                      const I Ap[], 
                      I       Aj[], 
                      T       Ax[])
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
void extract_csc_diagonal(const I n_row,
                          const I n_col, 
	                      const I Ap[], 
	                      const I Aj[], 
	                      const T Ax[],
	                      std::vector<T>*  Yx){
    extract_csr_diagonal(n_col, n_row, Ap, Aj, Ax, Yx);
}


template <class I, class T>
void csctocsr(const I n_row,
              const I n_col, 
              const I Ap[], 
              const I Ai[], 
              const T Ax[],
              std::vector<I>* Bp,
              std::vector<I>* Bj,
              std::vector<T>* Bx)
{ csrtocsc<I,T>(n_col,n_row,Ap,Ai,Ax,Bp,Bj,Bx); }

template <class I, class T>
void csctocoo(const I n_row,
              const I n_col, 
              const I Ap[], 
              const I Ai[], 
              const T Ax[],
              std::vector<I>* Bi,
              std::vector<I>* Bj,
              std::vector<T>* Bx)
{ csrtocoo<I,T>(n_col,n_row,Ap,Ai,Ax,Bj,Bi,Bx); }

template <class I, class T>
void cscmucsc(const I n_row,
              const I n_col, 
              const I Ap[], 
              const I Ai[], 
              const T Ax[],
              const I Bp[],
              const I Bi[],
              const T Bx[],
              std::vector<I>* Cp,
              std::vector<I>* Ci,
              std::vector<T>* Cx)
{ csrmucsr<I,T>(n_col,n_row,Bp,Bi,Bx,Ap,Ai,Ax,Cp,Ci,Cx); }

template<class I, class T>
void cootocsc(const I n_row,
      	      const I n_col,
      	      const I NNZ,
      	      const I Ai[],
      	      const I Aj[],
      	      const T Ax[],
      	      std::vector<I>* Bp,
      	      std::vector<I>* Bi,
      	      std::vector<T>* Bx)
{ cootocsr<I,T>(n_col,n_row,NNZ,Aj,Ai,Ax,Bp,Bi,Bx); }



template <class I, class T>
void csc_elmul_csc(const I n_row, const I n_col, 
                   const I Ap [], const I Ai [], const T Ax [],
                   const I Bp [], const I Bi [], const T Bx [],
                   std::vector<I>* Cp, std::vector<I>* Ci, std::vector<T>* Cx)
{
    csr_elmul_csr(n_col,n_row,Ap,Ai,Ax,Bp,Bi,Bx,Cp,Ci,Cx);
}

template <class I, class T>
void csc_eldiv_csc(const I n_row, const I n_col, 
                   const I Ap [], const I Ai [], const T Ax [],
                   const I Bp [], const I Bi [], const T Bx [],
                   std::vector<I>* Cp, std::vector<I>* Ci, std::vector<T>* Cx)
{
    csr_eldiv_csr(n_col,n_row,Ap,Ai,Ax,Bp,Bi,Bx,Cp,Ci,Cx);
}


template <class I, class T>
void csc_plus_csc(const I n_row, const I n_col, 
                  const I Ap [], const I Ai [], const T Ax [],
                  const I Bp [], const I Bi [], const T Bx [],
                  std::vector<I>* Cp, std::vector<I>* Ci, std::vector<T>* Cx)
{
    csr_plus_csr(n_col,n_row,Ap,Ai,Ax,Bp,Bi,Bx,Cp,Ci,Cx);
}

template <class I, class T>
void csc_minus_csc(const I n_row, const I n_col, 
                   const I Ap [], const I Ai [], const T Ax [],
                   const I Bp [], const I Bi [], const T Bx [],
                   std::vector<I>* Cp, std::vector<I>* Ci, std::vector<T>* Cx)
{
    csr_minus_csr(n_col,n_row,Ap,Ai,Ax,Bp,Bi,Bx,Cp,Ci,Cx);
}



template <class I, class T>
void sum_csc_duplicates(const I n_row,
                        const I n_col, 
                              I Ap[], 
                              I Ai[], 
                              T Ax[])
{ sum_csr_duplicates(n_col,n_row,Ap,Ai,Ax); }


template<class I, class T>
void sort_csc_indices(const I n_row,
                      const I n_col,
                      const I Ap[], 
                      I       Ai[], 
                      T       Ax[])
{ sort_csr_indices(n_col,n_row,Ap,Ai,Ax); }

#endif
