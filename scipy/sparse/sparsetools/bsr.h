#ifndef __BSR_H__
#define __BSR_H__

#include <vector>
#include <algorithm>
#include <functional>

#include "csr.h"
#include "fixed_size.h"


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
    const I N  = std::min(R*n_brow, C*n_bcol);
    const I RC = R*C;

    for(I i = 0; i < N; i++){
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
        const I end = (N/R) + (N % R == 0 ? 0 : 1);
        for(I i = 0; i < end; i++){
            for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
                const I base_row = R*i;
                const I base_col = C*Aj[jj];
                const T * base_val = Ax + RC*jj;

                for(I bi = 0; bi < R; bi++){
                    const I row = base_row + bi;
                    if (row >= N) break;

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
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            for(I bi = 0; bi < R; bi++){
                const T s = Xx[R*i + bi];
                T * block_row = Ax + RC*jj + C*bi;
                
                for(I bj = 0; bj < C; bj++){
                    block_row[bj] *= s;
                }

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

    //compute permutation of blocks using tranpose(CSR)
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


template <class I, class T, int R, int C, int N>
void bsr_matmat_pass2_fixed(const I n_brow,  const I n_bcol, 
      	                    const I Ap[],    const I Aj[],    const T Ax[],
      	                    const I Bp[],    const I Bj[],    const T Bx[],
      	                          I Cp[],          I Cj[],          T Cx[])
{
    const I RC = R*C;
    const I RN = R*N;
    const I NC = N*C;
    const I SIZE = RC*Cp[n_brow];

    for(I i = 0; i < SIZE; i++){
        Cx[i] = 0;
    }
 
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
                T * result = mats[k];
                matmat<R,C,N>(A,B,result);
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
void bsr_matmat_pass2(const I n_brow,  const I n_bcol, 
                      const I R,       const I C,       const I N,
      	              const I Ap[],    const I Aj[],    const T Ax[],
      	              const I Bp[],    const I Bj[],    const T Bx[],
      	                    I Cp[],          I Cj[],          T Cx[])
{
    assert(R > 0 && C > 0 && N > 0);

#ifdef SPARSETOOLS_TESTING
#define F(X,Y,Z) bsr_matmat_pass2_fixed<I,T,X,Y,Z>

    void (*dispatch[4][4][4])(I,I,const I*,const I*,const T*,
                                  const I*,const I*,const T*,
                                        I*,      I*,      T*) = \
    {
        { { F(1,1,1), F(1,1,2), F(1,1,3), F(1,1,4) },
          { F(1,2,1), F(1,2,2), F(1,2,3), F(1,2,4) },
          { F(1,3,1), F(1,3,2), F(1,3,3), F(1,3,4) },
          { F(1,4,1), F(1,4,2), F(1,4,3), F(1,4,4) },
        },
        { { F(2,1,1), F(2,1,2), F(2,1,3), F(2,1,4) },
          { F(2,2,1), F(2,2,2), F(2,2,3), F(2,2,4) },
          { F(2,3,1), F(2,3,2), F(2,3,3), F(2,3,4) },
          { F(2,4,1), F(2,4,2), F(2,4,3), F(2,4,4) },
        },
        { { F(3,1,1), F(3,1,2), F(3,1,3), F(3,1,4) },
          { F(3,2,1), F(3,2,2), F(3,2,3), F(3,2,4) },
          { F(3,3,1), F(3,3,2), F(3,3,3), F(3,3,4) },
          { F(3,4,1), F(3,4,2), F(3,4,3), F(3,4,4) },
        },
        { { F(4,1,1), F(4,1,2), F(4,1,3), F(4,1,4) },
          { F(4,2,1), F(4,2,2), F(4,2,3), F(4,2,4) },
          { F(4,3,1), F(4,3,2), F(4,3,3), F(4,3,4) },
          { F(4,4,1), F(4,4,2), F(4,4,3), F(4,4,4) },
        }
    };
    
    if (R <= 4 && C <= 4 && N <= 4){
        dispatch[R-1][N-1][C-1](n_brow,n_bcol,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx);
        return;
    }

#undef F
#endif

    if( R == 1 && N == 1 && C == 1 ){
        csr_matmat_pass2(n_brow, n_bcol, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx);
        return;
    }

    const I RC = R*C;
    const I RN = R*N;
    const I NC = N*C;
    const I SIZE = RC*Cp[n_brow];


    for(I i = 0; i < SIZE; i++){
        Cx[i] = 0;
    }
 
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
                T * result = mats[k];
                for(I r = 0; r < R; r++){
                    for(I c = 0; c < C; c++){
                        for(I n = 0; n < N; n++){
                            result[C*r + c] += A[N*r + n] * B[C*n + c];

                        }
                    }
                }
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



template <class I, class T, class bin_op>
void bsr_binop_bsr(const I n_brow, const I n_bcol, 
                   const I R,      const I C, 
                   const I Ap[],   const I Aj[],    const T Ax[],
                   const I Bp[],   const I Bj[],    const T Bx[],
                         I Cp[],         I Cj[],          T Cx[],
                   const bin_op& op)
{
    assert( R > 0 && C > 0);
    
    if( R == 1 && C == 1 ){
        csr_binop_csr(n_brow, n_bcol, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx, op); //use CSR for 1x1 blocksize 
        return;
    }

    const I RC = R*C;
    T * result = Cx;

    Cp[0] = 0;
    I nnz = 0;

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
                for(I n = 0; n < RC; n++){
                    result[n] = op(Ax[RC*A_pos + n],Bx[RC*B_pos + n]);
                }

                if( is_nonzero_block(result,RC) ){
                    Cj[nnz] = A_j;
                    result += RC;
                    nnz++;
                }

                A_j = Aj[++A_pos]; 
                B_j = Bj[++B_pos];

            } else if (A_j < B_j) {
                for(I n = 0; n < RC; n++){
                    result[n] = op(Ax[RC*A_pos + n],0);
                }

                if(is_nonzero_block(result,RC)){
                    Cj[nnz] = A_j;
                    result += RC;
                    nnz++;
                }

                A_j = Aj[++A_pos]; 

            } else {
                //B_j < A_j
                for(I n = 0; n < RC; n++){
                    result[n] = op(0,Bx[RC*B_pos + n]);
                }
                if(is_nonzero_block(result,RC)){
                    Cj[nnz] = B_j;
                    result += RC;
                    nnz++;
                }

                B_j = Bj[++B_pos];

            }
        }

        //tail
        while(A_pos < A_end){

            for(I n = 0; n < RC; n++){
                result[n] = op(Ax[RC*A_pos + n],0);
            }

            if(is_nonzero_block(result,RC)){
                Cj[nnz] = A_j;
                result += RC;
                nnz++;
            }

            A_j = Aj[++A_pos]; 

        }
        while(B_pos < B_end){
            for(I n = 0; n < RC; n++){
                result[n] = op(0,Bx[RC*B_pos + n]);
            }

            if(is_nonzero_block(result,RC)){
                Cj[nnz] = B_j;
                result += RC;
                nnz++;
            }

            B_j = Bj[++B_pos];

        }

        Cp[i+1] = nnz;
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

template <class I, class T, int R, int C>
void bsr_matvec_fixed(const I n_brow,
	                  const I n_bcol, 
	                  const I Ap[], 
	                  const I Aj[], 
	                  const T Ax[],
	                  const T Xx[],
	                        T Yx[])
{
    for(I i = 0; i < R*n_brow; i++){
        Yx[i] = 0;
    }

    for(I i = 0; i < n_brow; i++) {
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++) {
            I j = Aj[jj];
            matvec<R,C,1,1>(Ax + jj*R*C, Xx + j*C, Yx + i*R);
        }
    }
}

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

    if( R == 1 && C == 1 ){
        csr_matvec(n_brow, n_bcol, Ap, Aj, Ax, Xx, Yx); //use CSR for 1x1 blocksize 
        return;
    }

#ifdef SPARSETOOLS_TESTING
#define F(X,Y) bsr_matvec_fixed<I,T,X,Y>

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

#undef F
#endif

    //otherwise use general method

    for(I i = 0; i < R*n_brow; i++){
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


#endif
