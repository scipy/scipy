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
 *    01/09/2007 - index type is now templated
 *    01/06/2007 - initial inclusion into SciPy
 *
 */




#include <vector>





/*
 * Return zero of the appropriate type
 *
 *  this is a workaround for NumPy complex types 
 *  where T x = 0; doesn't make sense.
 *
 */
template <class T> 
T ZERO(){
  T temp = {0};
  return temp;
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
  
  *Bp = std::vector<I>(n_col+1);
  *Bi = std::vector<I>(NNZ);
  *Bx = std::vector<T>(NNZ);
 
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
	      const I Ap [], 
	      const I Aj[], 
	      const T Ax[],
	      std::vector<I>* Bi,
	      std::vector<I>* Bj,
	      std::vector<T>* Bx)
{
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
  *Cp = std::vector<I>(n_row+1,0);
  
  const T zero = ZERO<T>();

  std::vector<I> index(n_col,-1);
  std::vector<T> sums(n_col,zero);

  for(I i = 0; i < n_row; i++){
    I istart = -1;
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
      if(sums[istart] != zero){
	Cj->push_back(istart);
	Cx->push_back(sums[istart]);
      }
	
      I temp = istart;                
      istart = index[istart];
      
      index[temp] = -1; //clear arrays
      sums[temp]  = zero;                              
    }
    
    (*Cp)[i+1] = Cx->size();
  }
}




/*
 * Compute C = A+B for CSR matrices A,B
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
template <class I, class T>
void csrplcsr(const I n_row,
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

  *Cp = std::vector<I>(n_row+1,0);
  
  const T zero = ZERO<T>();

  std::vector<I> index(n_col,-1);
  std::vector<T>  sums(n_col,zero);

  for(I i = 0; i < n_row; i++){
    I istart = -1;
    I length =  0;
    
    //add a row of A to sums
    for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
      I j = Aj[jj];
      sums[j] += Ax[jj];
              
      if(index[j] == -1){
	index[j] = istart;                        
	istart = j;
	length++;
      }
    }
    
    //add a row of B to sums
    for(I jj = Bp[i]; jj < Bp[i+1]; jj++){
      I j = Bj[jj];
      sums[j] += Bx[jj];

      if(index[j] == -1){
	index[j] = istart;                        
	istart = j;
	length++;
      }
    }


    for(I jj = 0; jj < length; jj++){
      if(sums[istart] != zero){
	Cj->push_back(istart);
	Cx->push_back(sums[istart]);
      }
      
      I temp = istart;                
      istart = index[istart];
      
      index[temp] = -1;
      sums[temp]  = zero;                              
    }
    
    (*Cp)[i+1] = Cx->size();
  }
}

/*
 * Compute C = A (elmul) B for CSR matrices A,B
 *
 *   (elmul) - elementwise multiplication
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
template <class I, class T>
void csrelmulcsr(const I n_row,
		 const I n_col, 
		 const I Ap [], 
		 const I Aj[], 
		 const T Ax[],
		 const I Bp[],
		 const I Bj[],
		 const T Bx[],
		 std::vector<I>* Cp,
		 std::vector<I>* Cj,
		 std::vector<T>* Cx)
{
  *Cp = std::vector<I>(n_row+1,0);
  
  const T zero = ZERO<T>();

  std::vector<I>   index(n_col,-1);
  std::vector<T> A_row(n_col,zero);
  std::vector<T> B_row(n_col,zero);

  for(I i = 0; i < n_row; i++){
    I istart = -1;
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
      T prod = A_row[istart] * B_row[istart];
      
      if(prod != zero){
	Cj->push_back(istart);
	Cx->push_back(prod);
      }
      
      I temp = istart;                
      istart = index[istart];
      
      index[temp] = -1;
      A_row[temp] = zero;                              
      B_row[temp] = zero;
    }
    
    (*Cp)[i+1] = Cx->size();
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
  std::vector<I> tempBp(n_row+1,0);
  std::vector<I> tempBj(NNZ);
  std::vector<T> tempBx(NNZ);

  std::vector<I> nnz_per_row(n_row,0); //temp array

  //compute nnz per row, then compute Bp
  for(I i = 0; i < NNZ; i++){
    nnz_per_row[Ai[i]]++;
  }
  for(I i = 0, cumsum = 0; i < n_row; i++){
    tempBp[i]      = cumsum;
    cumsum        += nnz_per_row[i];
    nnz_per_row[i] = 0; //reset count
  }
  tempBp[n_row] = NNZ;


  //write Aj,Ax Io tempBj,tempBx
  for(I i = 0; i < NNZ; i++){
    I row = Ai[i];
    I n   = tempBp[row] + nnz_per_row[row];

    tempBj[n] = Aj[i];
    tempBx[n] = Ax[i];

    nnz_per_row[row]++;
  }
  //now tempBp,tempBj,tempBx form a CSR representation (with duplicates)


  //use (tempB + 0) to sum duplicates
  std::vector<I> Xp(n_row+1,0); //row pointer for an empty matrix

  csrplcsr<I,T>(n_row,n_col,
		&tempBp[0],&tempBj[0],&tempBx[0],
		&Xp[0],NULL,NULL,
		Bp,Bj,Bx);    	   
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
 *   vec<T> Yx - nonzeros (real part)
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
	    const I Ap [], 
	    const I Aj[], 
	    const T Ax[],
	    const T Xx[],
	    std::vector<T>*  Yx)
{
  const T zero = ZERO<T>();

  *Yx = std::vector<T>(n_row,zero);

  for(I i = 0; i < n_row; i++){
    I row_start = Ap[i];
    I row_end   = Ap[i+1];
    
    T& Yx_i = (*Yx)[i];
    for(I jj = row_start; jj < row_end; jj++){
      Yx_i += Ax[jj] * Xx[Aj[jj]];
    }
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
 *   T    Ax[n_col]     - nonzeros (real part)
 *   T    Xx[n_col]     - nonzeros (real part)
 *   bool do_complex    - switch scalar/complex modes
 *
 * Output Arguments:
 *   vec<T> Yx - nonzeros (real part)
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
  const T zero = ZERO<T>();

  *Yx = std::vector<T>(n_row,zero);
  
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
  const T  zero  = ZERO<T>();
  const T* x_ptr = Mx;

  Ap->push_back(0);
  for(I i = 0; i < n_row; i++){
    for(I j = 0; j < n_col; j++){
      if(*x_ptr != zero){
	Aj->push_back(j);
	Ax->push_back(*x_ptr);
      }
      x_ptr++;
    }
    Ap->push_back(Aj->size());
  }
}



/*
 * Derived methods
 */
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

template <class I, class T>
void cscplcsc(const I n_row,
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
{ csrplcsr<I,T>(n_col,n_row,Ap,Ai,Ax,Bp,Bi,Bx,Cp,Ci,Cx); }

template <class I, class T>
void cscelmulcsc(const I n_row,
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
{ csrelmulcsr<I,T>(n_col,n_row,Ap,Ai,Ax,Bp,Bi,Bx,Cp,Ci,Cx); }

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

/* Taken from numpy. */
#define PYA_QS_STACK 100
#define SMALL_QUICKSORT 15
#define STDC_LT(a,b) ((a) < (b))
#define STDC_LE(a,b) ((a) <= (b))
#define STDC_EQ(a,b) ((a) == (b))
#define SWAP(a,b) {SWAP_temp = (b); (b)=(a); (a) = SWAP_temp;}
template<class I, class Ip>
void int_aquicksort(I *v, Ip* tosort, Ip num, void *unused)
{
  I vp;
  Ip *pl, *pr, SWAP_temp;
  Ip *stack[PYA_QS_STACK], **sptr=stack, *pm, *pi, *pj, *pt, vi;

  pl = tosort;
  pr = tosort + num - 1;

  for(;;) {
    while ((pr - pl) > SMALL_QUICKSORT) {
      /* quicksort partition */
      pm = pl + ((pr - pl) >> 1);
      if (STDC_LT(v[*pm],v[*pl])) SWAP(*pm,*pl);
      if (STDC_LT(v[*pr],v[*pm])) SWAP(*pr,*pm);
      if (STDC_LT(v[*pm],v[*pl])) SWAP(*pm,*pl);
      vp = v[*pm];
      pi = pl;
      pj = pr - 1;
      SWAP(*pm,*pj);
      for(;;) {
	do ++pi; while (STDC_LT(v[*pi],vp));
	do --pj; while (STDC_LT(vp,v[*pj]));
	if (pi >= pj)  break;
	SWAP(*pi,*pj);
      }
      SWAP(*pi,*(pr-1));
      /* push largest partition on stack */
      if (pi - pl < pr - pi) {
	*sptr++ = pi + 1;
	*sptr++ = pr;
	pr = pi - 1;
      }else{
	*sptr++ = pl;
	*sptr++ = pi - 1;
	pl = pi + 1;
      }
    }
    /* insertion sort */
    for(pi = pl + 1; pi <= pr; ++pi) {
      vi = *pi;
      vp = v[vi];
      for(pj = pi, pt = pi - 1; \
	    pj > pl && STDC_LT(vp, v[*pt]);)
	{
	  *pj-- = *pt--;
	}
      *pj = vi;
    }
    if (sptr == stack) break;
    pr = *(--sptr);
    pl = *(--sptr);
  }
}

template<class I, class T>
void ensure_sorted_indices(const I n_row,
			   const I n_col,
			   const I Ap[], 
			   I Aj[], 
			   T Ax[])
{
  const T zero = ZERO<T>();
  I isort[ n_col ];
  std::vector<I> itemp(n_col,0);
  std::vector<T> atemp(n_col,zero);

  for(I i = 0; i < n_row; i++){
    I row_start = Ap[i];
    I row_end   = Ap[i+1];
    I ncol = row_end - row_start;
    I ii;

    for(I jj = 0; jj < ncol; jj++){
      isort[jj] = jj;
      atemp[jj] = Ax[row_start + jj];
      itemp[jj] = Aj[row_start + jj];
    }    
    int_aquicksort( Aj + row_start, isort, ncol, 0 );
    
    /* Permute in-place both Aj and Ax of row i. */
    for(I jj = row_start; jj < row_end; jj++){
      ii = isort[jj-row_start];
      Aj[jj] = itemp[ii];
      Ax[jj] = atemp[ii];
    }
  }
}
			   

#endif
