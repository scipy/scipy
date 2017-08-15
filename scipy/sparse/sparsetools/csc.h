#ifndef __CSC_H__
#define __CSC_H__

#include "csr.h"

/*
 * Compute Y += A*X for CSC matrix A and dense vectors X,Y
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
    for(I j = 0; j < n_col; j++){
        I col_start = Ap[j];
        I col_end   = Ap[j+1];

        for(I ii = col_start; ii < col_end; ii++){
            I i    = Ai[ii];
            Yx[i] += Ax[ii] * Xx[j];
        }
    }
}


/*
 * Compute Y += A*X for CSC matrix A and dense block vectors X,Y
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
 * Note:
 *   Output array Yx must be preallocated
 *
 */
template <class I, class T>
void csc_matvecs(const I n_row,
                 const I n_col,
                 const I n_vecs,
                 const I Ap[],
                 const I Ai[],
                 const T Ax[],
                 const T Xx[],
                       T Yx[])
{
    for(I j = 0; j < n_col; j++){
        for(I ii = Ap[j]; ii < Ap[j+1]; ii++){
            const I i = Ai[ii];
            axpy(n_vecs, Ax[ii], Xx + (npy_intp)n_vecs * j, Yx + (npy_intp)n_vecs * i);
        }
    }
}




/*
 * Derived methods
 */
template <class I, class T>
void csc_diagonal(const I k,
                  const I n_row,
                  const I n_col,
                  const I Ap[],
                  const I Aj[],
                  const T Ax[],
                        T Yx[])
{ csr_diagonal(-k, n_col, n_row, Ap, Aj, Ax, Yx); }


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

template <class I, class T, class T2>
void csc_ne_csc(const I n_row, const I n_col,
                const I Ap[], const I Ai[], const T Ax[],
                const I Bp[], const I Bi[], const T Bx[],
                      I Cp[],       I Ci[],      T2 Cx[])
{
    csr_ne_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}

template <class I, class T, class T2>
void csc_lt_csc(const I n_row, const I n_col,
                const I Ap[], const I Ai[], const T Ax[],
                const I Bp[], const I Bi[], const T Bx[],
                      I Cp[],       I Ci[],      T2 Cx[])
{
    csr_lt_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}

template <class I, class T, class T2>
void csc_gt_csc(const I n_row, const I n_col,
                const I Ap[], const I Ai[], const T Ax[],
                const I Bp[], const I Bi[], const T Bx[],
                      I Cp[],       I Ci[],      T2 Cx[])
{
    csr_gt_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}

template <class I, class T, class T2>
void csc_le_csc(const I n_row, const I n_col,
                const I Ap[], const I Ai[], const T Ax[],
                const I Bp[], const I Bi[], const T Bx[],
                      I Cp[],       I Ci[],      T2 Cx[])
{
    csr_le_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}

template <class I, class T, class T2>
void csc_ge_csc(const I n_row, const I n_col,
                const I Ap[], const I Ai[], const T Ax[],
                const I Bp[], const I Bi[], const T Bx[],
                      I Cp[],       I Ci[],      T2 Cx[])
{
    csr_ge_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}

template <class I, class T>
void csc_elmul_csc(const I n_row, const I n_col,
                   const I Ap[], const I Ai[], const T Ax[],
                   const I Bp[], const I Bi[], const T Bx[],
                         I Cp[],       I Ci[],       T Cx[])
{
    csr_elmul_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}

template <class I, class T>
void csc_eldiv_csc(const I n_row, const I n_col,
                   const I Ap[], const I Ai[], const T Ax[],
                   const I Bp[], const I Bi[], const T Bx[],
                         I Cp[],       I Ci[],       T Cx[])
{
    csr_eldiv_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}


template <class I, class T>
void csc_plus_csc(const I n_row, const I n_col,
                  const I Ap[], const I Ai[], const T Ax[],
                  const I Bp[], const I Bi[], const T Bx[],
                        I Cp[],       I Ci[],       T Cx[])
{
    csr_plus_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}

template <class I, class T>
void csc_minus_csc(const I n_row, const I n_col,
                   const I Ap[], const I Ai[], const T Ax[],
                   const I Bp[], const I Bi[], const T Bx[],
                         I Cp[],       I Ci[],       T Cx[])
{
    csr_minus_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}


template <class I, class T>
void csc_maximum_csc(const I n_row, const I n_col,
                     const I Ap[], const I Ai[], const T Ax[],
                     const I Bp[], const I Bi[], const T Bx[],
                           I Cp[],       I Ci[],       T Cx[])
{
    csr_maximum_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}

template <class I, class T>
void csc_minimum_csc(const I n_row, const I n_col,
                     const I Ap[], const I Ai[], const T Ax[],
                     const I Bp[], const I Bi[], const T Bx[],
                           I Cp[],       I Ci[],       T Cx[])
{
    csr_minimum_csr(n_col, n_row, Ap, Ai, Ax, Bp, Bi, Bx, Cp, Ci, Cx);
}


#endif
