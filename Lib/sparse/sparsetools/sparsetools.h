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
 *    01/06/2007 - initial inclusion into SciPy
 *
 */




#include <vector>





/*
 * Return zero of the appropriate type
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
 *   int  n_row         - number of rows in A
 *   int  n_col         - number of columns in A
 *   int  Ap[n_row+1]   - row pointer
 *   int  Aj[nnz(A)]    - column indices
 *   T    Ax[nnz(A)]    - nonzeros
 *
 * Output Arguments:
 *   vec<int> Bp  - row pointer
 *   vec<int> Bj  - column indices
 *   vec<T>   Bx  - nonzeros
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
template <class T>
void csrtocsc(const int n_row,
	      const int n_col, 
	      const int Ap[], 
	      const int Aj[], 
	      const T   Ax[],
	      std::vector<int>* Bp,
	      std::vector<int>* Bi,
	      std::vector<T>*   Bx)
{  
  int NNZ = Ap[n_row];
  
  *Bp = std::vector<int>(n_col+1);
  *Bi = std::vector<int>(NNZ);
  *Bx = std::vector<T>(NNZ);
 
  std::vector<int> nnz_per_col(n_col,0); //temp array
 
  //compute number of non-zero entries per column of A 
  for (int i = 0; i < NNZ; i++){            
    nnz_per_col[Aj[i]]++;
  }
        
  //cumsum the nnz_per_col to get Bp[]
  for(int i = 0, cumsum = 0; i < n_col; i++){     
    (*Bp)[i]   = cumsum; 
    cumsum += nnz_per_col[i];
    nnz_per_col[i] = 0;              //reset count
  }
  (*Bp)[n_col] = NNZ;
  
  for(int i = 0; i < n_row; i++){
    int row_start = Ap[i];
    int row_end   = Ap[i+1];
    for(int j = row_start; j < row_end; j++){
      int col = Aj[j];
      int k   = (*Bp)[col] + nnz_per_col[col];

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
 *   int  n_row         - number of rows in A
 *   int  n_col         - number of columns in A
 *   int  Ap[n_row+1]   - row pointer
 *   int  Aj[nnz(A)]    - column indices
 *   T    Ax[nnz(A)]    - nonzeros
 *
 * Output Arguments:
 *   vec<int> Bi  - row indices
 *   vec<int> Bj  - column indices
 *   vec<T>   Bx  - nonzeros
 *
 * Note:
 *   Output arrays Bi,Bj,Bx will be allocated within in the method
 *
 * Note: 
 *   Complexity: Linear.
 * 
 */
template<class T>
void csrtocoo(const int n_row,
	      const int n_col, 
	      const int Ap [], 
	      const int Aj[], 
	      const T   Ax[],
	      std::vector<int>*    Bi,
	      std::vector<int>*    Bj,
	      std::vector<T>* Bx)
{
  for(int i = 0; i < n_row; i++){
    int row_start = Ap[i];
    int row_end   = Ap[i+1];
    for(int jj = row_start; jj < row_end; jj++){
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
 *   int  n_row       - number of rows in A
 *   int  n_col       - number of columns in B (hence C is n_row by n_col)
 *   int  Ap[n_row+1] - row pointer
 *   int  Aj[nnz(A)]  - column indices
 *   T    Ax[nnz(A)]  - nonzeros
 *   int  Bp[?]       - row pointer
 *   int  Bj[nnz(B)]  - column indices
 *   T    Bx[nnz(B)]  - nonzeros
 * Output Arguments:
 *   vec<int> Cp - row pointer
 *   vec<int> Cj - column indices
 *   vec<T>   Cx - nonzeros
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
template<class T>
void csrmucsr(const int n_row,
	      const int n_col, 
	      const int Ap [], 
	      const int Aj[], 
	      const T Ax[],
	      const int Bp[],
	      const int Bj[],
	      const T Bx[],
	      std::vector<int>* Cp,
	      std::vector<int>* Cj,
	      std::vector<T>* Cx)
{
  *Cp = std::vector<int>(n_row+1,0);
  Cj->clear();        
  Cx->clear();
  
  const T zero = ZERO<T>();

  std::vector<int>    index(n_col,-1);
  std::vector<T> sums(n_col,zero);

  for(int i = 0; i < n_row; i++){
    int istart = -1;
    int length =  0;
    
    for(int jj = Ap[i]; jj < Ap[i+1]; jj++){
      int j = Aj[jj];
      for(int kk = Bp[j]; kk < Bp[j+1]; kk++){
	int k = Bj[kk];
        
	sums[k] += Ax[jj]*Bx[kk];
        
	if(index[k] == -1){
	  index[k] = istart;                        
	  istart = k;
	  length++;
	}
      }
    }         

    for(int jj = 0; jj < length; jj++){
      if(sums[istart] != zero){
	Cj->push_back(istart);
	Cx->push_back(sums[istart]);
      }
	
      int temp = istart;                
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
 *   int    n_row       - number of rows in A (and B)
 *   int    n_col       - number of columns in A (and B)
 *   int    Ap[n_row+1] - row pointer
 *   int    Aj[nnz(A)]  - column indices
 *   T      Ax[nnz(A)]  - nonzeros
 *   int    Bp[?]       - row pointer
 *   int    Bj[nnz(B)]  - column indices
 *   T      Bx[nnz(B)]  - nonzeros
 * Output Arguments:
 *   vec<int> Cp  - row pointer
 *   vec<int> Cj  - column indices
 *   vec<T>   Cx  - nonzeros
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
template <class T>
void csrplcsr(const int n_row,
	      const int n_col, 
	      const int Ap[], 
	      const int Aj[], 
	      const T   Ax[],
	      const int Bp[],
	      const int Bj[],
	      const T   Bx[],
	      std::vector<int>* Cp,
	      std::vector<int>* Cj,
	      std::vector<T>  * Cx)
{

  *Cp = std::vector<int>(n_row+1,0);
  Cj->clear();        
  Cx->clear();
  
  const T zero = ZERO<T>();

  std::vector<int> index(n_col,-1);
  std::vector<T>   sums(n_col,zero);

  for(int i = 0; i < n_row; i++){
    int istart = -1;
    int length =  0;
    
    //add a row of A to sums
    for(int jj = Ap[i]; jj < Ap[i+1]; jj++){
      int j = Aj[jj];
      sums[j] += Ax[jj];
              
      if(index[j] == -1){
	index[j] = istart;                        
	istart = j;
	length++;
      }
    }
    
    //add a row of B to sums
    for(int jj = Bp[i]; jj < Bp[i+1]; jj++){
      int j = Bj[jj];
      sums[j] += Bx[jj];

      if(index[j] == -1){
	index[j] = istart;                        
	istart = j;
	length++;
      }
    }


    for(int jj = 0; jj < length; jj++){
      if(sums[istart] != zero){
	Cj->push_back(istart);
	Cx->push_back(sums[istart]);
      }
      
      int temp = istart;                
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
 *   int    n_row       - number of rows in A (and B)
 *   int    n_col       - number of columns in A (and B)
 *   int    Ap[n_row+1] - row pointer
 *   int    Aj[nnz(A)]  - column indices
 *   T      Ax[nnz(A)]  - nonzeros
 *   int    Bp[?]       - row pointer
 *   int    Bj[nnz(B)]  - column indices
 *   T      Bx[nnz(B)]  - nonzeros
 * Output Arguments:
 *   vec<int> Cp  - row pointer
 *   vec<int> Cj  - column indices
 *   vec<T>   Cx  - nonzeros
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
template <class T>
void csrelmulcsr(const int n_row,
		 const int n_col, 
		 const int Ap [], 
		 const int Aj[], 
		 const T   Ax[],
		 const int Bp[],
		 const int Bj[],
		 const T   Bx[],
		 std::vector<int>* Cp,
		 std::vector<int>* Cj,
		 std::vector<T>*   Cx)
{
  *Cp = std::vector<int>(n_row+1,0);
  Cj->clear();        
  Cx->clear();
  
  const T zero = ZERO<T>();

  std::vector<int>   index(n_col,-1);
  std::vector<T> A_row(n_col,zero);
  std::vector<T> B_row(n_col,zero);

  for(int i = 0; i < n_row; i++){
    int istart = -1;
    int length =  0;
    
    //add a row of A to A_row
    for(int jj = Ap[i]; jj < Ap[i+1]; jj++){
      int j = Aj[jj];

      A_row[j] += Ax[jj];
      
      if(index[j] == -1){
	index[j] = istart;                        
	istart = j;
	length++;
      }
    }
    
    //add a row of B to B_row
    for(int jj = Bp[i]; jj < Bp[i+1]; jj++){
      int j = Bj[jj];

      B_row[j] += Bx[jj];

      if(index[j] == -1){
	index[j] = istart;                        
	istart = j;
	length++;
      }
    }


    for(int jj = 0; jj < length; jj++){
      T prod = A_row[istart] * B_row[istart];
      
      if(prod != zero){
	Cj->push_back(istart);
	Cx->push_back(prod);
      }
      
      int temp = istart;                
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
 *   int  n_row         - number of rows in A
 *   int  n_col         - number of columns in A
 *   int  Ai[nnz(A)]    - row indices
 *   int  Aj[nnz(A)]    - column indices
 *   T    Ax[nnz(A)]    - nonzeros
 * Output Arguments:
 *   vec<int> Bp        - row pointer
 *   vec<int> Bj        - column indices
 *   vec<T>   Bx        - nonzeros
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
template<class T>
void cootocsr(const int n_row,
	      const int n_col,
	      const int NNZ,
	      const int Ai[],
	      const int Aj[],
	      const T   Ax[],
	      std::vector<int>* Bp,
	      std::vector<int>* Bj,
	      std::vector<T>* Bx)
{
  std::vector<int> tempBp(n_row+1,0);
  std::vector<int> tempBj(NNZ);
  std::vector<T>   tempBx(NNZ);

  std::vector<int> nnz_per_row(n_row,0); //temp array

  //compute nnz per row, then compute Bp
  for(int i = 0; i < NNZ; i++){
    nnz_per_row[Ai[i]]++;
  }
  for(int i = 0, cumsum = 0; i < n_row; i++){
    tempBp[i]      = cumsum;
    cumsum        += nnz_per_row[i];
    nnz_per_row[i] = 0; //reset count
  }
  tempBp[n_row] = NNZ;


  //write Aj,Ax into tempBj,tempBx
  for(int i = 0; i < NNZ; i++){
    int row = Ai[i];
    int n   = tempBp[row] + nnz_per_row[row];

    tempBj[n] = Aj[i];
    tempBx[n] = Ax[i];

    nnz_per_row[row]++;
  }
  //now tempBp,tempBj,tempBx form a CSR representation (with duplicates)


  //use (tempB + 0) to sum duplicates
  std::vector<int> Xp(n_row+1,0); //row pointer for an empty matrix

  csrplcsr<T>(n_row,n_col,
	      &tempBp[0],&tempBj[0],&tempBx[0],
	      &Xp[0],NULL,NULL,
	      Bp,Bj,Bx);    	   
}
	    




/*
 * Compute Y = A*X for CSR matrix A and dense vectors X,Y
 *
 *
 * Input Arguments:
 *   int  n_row         - number of rows in A
 *   int  n_col         - number of columns in A
 *   int  Ap[n_row+1]   - row pointer
 *   int  Aj[nnz(A)]    - column indices
 *   T    Ax[n_col]     - nonzeros
 *   T    Xx[n_col]     - nonzeros
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
template <class T>
void csrmux(const int n_row,
	    const int n_col, 
	    const int Ap [], 
	    const int Aj[], 
	    const T   Ax[],
	    const T   Xx[],
	    std::vector<T>*  Yx)
{
  const T zero = ZERO<T>();

  *Yx = std::vector<T>(n_row,zero);

  for(int i = 0; i < n_row; i++){
    int row_start = Ap[i];
    int row_end   = Ap[i+1];
    
    T& Yx_i = (*Yx)[i];
    for(int jj = row_start; jj < row_end; jj++){
      Yx_i += Ax[jj] * Xx[Aj[jj]];
    }
  }
}



/*
 * Compute Y = A*X for CSC matrix A and dense vectors X,Y
 *
 *
 * Input Arguments:
 *   int  n_row         - number of rows in A
 *   int  n_col         - number of columns in A
 *   int  Ap[n_row+1]   - column pointer
 *   int  Ai[nnz(A)]    - row indices
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
template <class T>
void cscmux(const int n_row,
	    const int n_col, 
	    const int Ap[], 
	    const int Ai[], 
	    const T   Ax[],
	    const T   Xx[],
	    std::vector<T>*  Yx)
{
  const T zero = ZERO<T>();

  *Yx = std::vector<T>(n_row,zero);
  
  for(int j = 0; j < n_col; j++){
    int col_start = Ap[j];
    int col_end   = Ap[j+1];
    
    for(int ii = col_start; ii < col_end; ii++){
      int row  = Ai[ii];
      (*Yx)[row] += Ax[ii] * Xx[j];
    }
  }
}




/*
 * Construct CSC matrix A from diagonals
 *
 * Input Arguments:
 *   int  n_row                            - number of rows in A
 *   int  n_col                            - number of columns in A
 *   int  n_diags                          - number of diagonals
 *   int  diags_indx[n_diags]              - where to place each diagonal 
 *   T    diags[n_diags][min(n_row,n_col)] - diagonals
 *
 * Output Arguments:
 *   vec<int> Ap  - row pointer
 *   vec<int> Aj  - column indices
 *   vec<T>   Ax  - nonzeros
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
template<class T>
void spdiags(const int n_row,
	     const int n_col,
	     const int n_diag,
	     const int offsets[],
	     const T   diags[],
	     std::vector<int> * Ap,
	     std::vector<int> * Ai,
	     std::vector<T>   * Ax)
{
  const int diags_length = std::min(n_row,n_col);
  Ap->push_back(0);

  for(int i = 0; i < n_col; i++){
    for(int j = 0; j < n_diag; j++){
      if(offsets[j] <= 0){              //sub-diagonal
	int row = i - offsets[j];
	if (row >= n_row){ continue; }
	
	Ai->push_back(row);
	Ax->push_back(diags[j*diags_length + i]);

      } else {                          //super-diagonal
	int row = i - offsets[j];
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
 *   int  n_row           - number of rows in A
 *   int  n_col           - number of columns in A
 *   int  Ap[n_row+1]     - row pointer
 *   int  Aj[nnz(A)]      - column indices
 *   T    Ax[nnz(A)]      - nonzeros 
 *   T    Mx[n_row*n_col] - dense matrix
 *
 * Note:
 *   Output arrays Mx are assumed to be allocated and
 *   initialized to 0 by the caller.
 *
 */
template<class T>
void csrtodense(const int  n_row,
		const int  n_col,
		const int  Ap[],
		const int  Aj[],
		const T    Ax[],
		      T    Mx[])
{
  int row_base = 0;
  for(int i = 0; i < n_row; i++){
    int row_start = Ap[i];
    int row_end   = Ap[i+1];
    for(int jj = row_start; jj < row_end; jj++){
      int j = Aj[jj];

      Mx[row_base + j] = Ax[jj];
    }	
    row_base +=  n_col;
  }
}



/*
 * Compute A = M for CSR matrix A, dense matrix M
 *
 * Input Arguments:
 *   int  n_row           - number of rows in A
 *   int  n_col           - number of columns in A
 *   T    Mx[n_row*n_col] - dense matrix
 *   int  Ap[n_row+1]     - row pointer
 *   int  Aj[nnz(A)]      - column indices
 *   T    Ax[nnz(A)]      - nonzeros 
 *
 * Note:
 *    Output arrays Ap,Aj,Ax will be allocated within the method
 *
 */
template<class T>
void densetocsr(const int  n_row,
		const int  n_col,
		const T    Mx[],
		std::vector<int>* Ap,
		std::vector<int>* Aj,
		std::vector<T>*   Ax)
{
  const T zero = ZERO<T>();
  const T* x_ptr = Mx;

  Ap->push_back(0);
  for(int i = 0; i < n_row; i++){
    for(int j = 0; j < n_col; j++){
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
template <class T>
void csctocsr(const int n_row,
	      const int n_col, 
	      const int Ap[], 
	      const int Ai[], 
	      const T   Ax[],
	      std::vector<int>* Bp,
	      std::vector<int>* Bj,
	      std::vector<T>*   Bx)
{ csrtocsc<T>(n_col,n_row,Ap,Ai,Ax,Bp,Bj,Bx); }

template<class T>
void csctocoo(const int n_row,
	      const int n_col, 
	      const int Ap[], 
	      const int Ai[], 
	      const T   Ax[],
	      std::vector<int>*    Bi,
	      std::vector<int>*    Bj,
	      std::vector<T>*      Bx)
{ csrtocoo<T>(n_col,n_row,Ap,Ai,Ax,Bj,Bi,Bx); }

template<class T>
void cscmucsc(const int n_row,
	      const int n_col, 
	      const int Ap [], 
	      const int Ai[], 
	      const T   Ax[],
	      const int Bp[],
	      const int Bi[],
	      const T   Bx[],
	      std::vector<int>* Cp,
	      std::vector<int>* Ci,
	      std::vector<T>  * Cx)
{ csrmucsr<T>(n_col,n_row,Bp,Bi,Bx,Ap,Ai,Ax,Cp,Ci,Cx); }

template <class T>
void cscplcsc(const int n_row,
	      const int n_col, 
	      const int Ap [], 
	      const int Ai[], 
	      const T   Ax[],
	      const int Bp[],
	      const int Bi[],
	      const T   Bx[],
	      std::vector<int>* Cp,
	      std::vector<int>* Ci,
	      std::vector<T>*   Cx)
{ csrplcsr<T>(n_col,n_row,Ap,Ai,Ax,Bp,Bi,Bx,Cp,Ci,Cx); }

template <class T>
void cscelmulcsc(const int n_row,
		 const int n_col, 
		 const int Ap [], 
		 const int Ai[], 
		 const T   Ax[],
		 const int Bp[],
		 const int Bi[],
		 const T   Bx[],
		 std::vector<int>* Cp,
		 std::vector<int>* Ci,
		 std::vector<T>*   Cx)
{ csrelmulcsr<T>(n_col,n_row,Ap,Ai,Ax,Bp,Bi,Bx,Cp,Ci,Cx); }

template<class T>
void cootocsc(const int n_row,
	      const int n_col,
	      const int NNZ,
	      const int Ai[],
	      const int Aj[],
	      const T   Ax[],
	      std::vector<int>* Bp,
	      std::vector<int>* Bi,
	      std::vector<T>*   Bx)
{ cootocsr<T>(n_col,n_row,NNZ,Aj,Ai,Ax,Bp,Bi,Bx); }



#endif
