/* 
 * These are sparsetools functions that are not currently used
 * 
 */

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
 *   Output arrays Cp, Cj, and Cx will be allocated within in the method
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

    std::vector<I> next(n_col,-1);
    std::vector<T> sums(n_col, 0);

    for(I i = 0; i < n_row; i++){
        I head = -2;
        I length =  0;

        I jj_start = Ap[i];
        I jj_end   = Ap[i+1];
        for(I jj = jj_start; jj < jj_end; jj++){
            I j = Aj[jj];

            I kk_start = Bp[j];
            I kk_end   = Bp[j+1];
            for(I kk = kk_start; kk < kk_end; kk++){
                I k = Bj[kk];

                sums[k] += Ax[jj]*Bx[kk];

                if(next[k] == -1){
                    next[k] = head;                        
                    head = k;
                    length++;
                }
            }
        }         

        for(I jj = 0; jj < length; jj++){
            if(sums[head] != 0){
                Cj->push_back(head);
                Cx->push_back(sums[head]);
            }

            I temp = head;                
            head = next[head];

            next[temp] = -1; //clear arrays
            sums[temp]  =  0;                              
        }

        (*Cp)[i+1] = Cx->size();
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
 *   T  Ax[nnz(A)]      - nonzeros 
 *   T  Mx[n_row*n_col] - dense matrix
 *
 * Note:
 *   Output array Mx is assumed to be allocated and
 *   initialized to 0 by the caller.
 *
 */
template <class I, class T>
void csr_todense(const I  n_row,
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
 *   Output arrays Bi, Bj, Bx will be allocated within in the method
 *
 * Note: 
 *   Complexity: Linear.
 * 
 */
template <class I, class T>
void csr_tocoo(const I n_row,
	           const I n_col, 
               const I Ap[], 
               const I Aj[], 
               const T Ax[],
               std::vector<I>* Bi,
               std::vector<I>* Bj,
               std::vector<T>* Bx)
{
  I nnz = Ap[n_row];
  Bi->reserve(nnz);
  Bi->reserve(nnz);
  Bx->reserve(nnz);
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
 *   Output arrays Ap, Aj, Ax will be allocated within in the method
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

template <class I, class T, int R, int C, class bin_op>
void bsr_binop_bsr_fixed(const I n_brow, const I n_bcol, 
                         const I Ap[],   const I Aj[],    const T Ax[],
                         const I Bp[],   const I Bj[],    const T Bx[],
                               I Cp[],         I Cj[],          T Cx[],
                         const bin_op& op)
{
   //Method that works for unsorted indices
    const I RC = R*C;
    T zeros[RC] = {0};
    Cp[0] = 0;
    I nnz = 0;

    std::cout << "using bsr_ fixed" << std::endl;
    for(I i = 0; i < n_brow; i++){
        I A_pos = Ap[i];
        I B_pos = Bp[i];
        I A_end = Ap[i+1];
        I B_end = Bp[i+1];

        I A_j = Aj[A_pos];
        I B_j = Bj[B_pos];
            
        //while not finished with either row
        while(A_pos < A_end && B_pos < B_end){
            if(A_j == B_j){
                Cj[nnz] = A_j;
                vec_binop_vec<RC> (Ax + RC*A_pos, Bx + RC*B_pos, Cx + RC*nnz, op);
                if( is_nonzero_block(Cx + RC*nnz,RC) ){
                    nnz++;
                }
                A_j = Aj[++A_pos]; 
                B_j = Bj[++B_pos];
            } else if (A_j < B_j) {
                Cj[nnz] = A_j;
                vec_binop_vec<RC> (Ax + RC*A_pos, zeros, Cx + RC*nnz, op);
                if( is_nonzero_block(Cx + RC*nnz,RC) ){
                    nnz++;
                }
                A_j = Aj[++A_pos]; 
            } else {
                //B_j < A_j
                Cj[nnz] = B_j;
                vec_binop_vec<RC> (zeros, Bx + RC*A_pos, Cx + RC*nnz, op);
                if( is_nonzero_block(Cx + RC*nnz,RC) ){
                    nnz++;
                }
                B_j = Bj[++B_pos];
            }
        }

        //tail
        while(A_pos < A_end){
            Cj[nnz] = A_j;
            vec_binop_vec<RC> (Ax + RC*A_pos, zeros, Cx + RC*nnz, op);
            if( is_nonzero_block(Cx + RC*nnz,RC) ){
                nnz++;
            }
            A_j = Aj[++A_pos]; 
        }
        while(B_pos < B_end){
            Cj[nnz] = B_j;
            vec_binop_vec<RC> (zeros, Bx + RC*A_pos, Cx + RC*nnz, op);
            if( is_nonzero_block(Cx + RC*nnz,RC) ){
                nnz++;
            }
            B_j = Bj[++B_pos];
        }

        Cp[i+1] = nnz;
    }
}


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
    // method that uses O(1) temp storage
    const I hash_size = 1 << 5;
    I vals[hash_size];
    I mask[hash_size];

    std::set<I> spill;    
    
    for(I i = 0; i < hash_size; i++){
        vals[i] = -1;
        mask[i] = -1;
    }

    Cp[0] = 0;

    I slow_inserts = 0;
    I total_inserts = 0;
    I nnz = 0;
    for(I i = 0; i < n_row; i++){
        spill.clear();
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            I j = Aj[jj];
            for(I kk = Bp[j]; kk < Bp[j+1]; kk++){
                I k = Bj[kk];
                // I hash = k & (hash_size - 1);
                I hash = ((I)2654435761 * k) & (hash_size -1 );
                total_inserts++;
                if(mask[hash] != i){
                    mask[hash] = i;                        
                    vals[hash] = k;
                    nnz++;
                } else {
                    if (vals[hash] != k){
                        slow_inserts++;
                        spill.insert(k);
                    }
                }
            }
        }       
        nnz += spill.size();
        Cp[i+1] = nnz;
    }

    std::cout << "slow fraction " << ((float) slow_inserts)/ ((float) total_inserts) << std::endl;
}


