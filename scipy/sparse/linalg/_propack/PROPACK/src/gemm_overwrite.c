#include "gemm_overwrite.h"


void sgemm_ovwr(const int transa, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, float alpha,
                float* restrict A, CBLAS_INT lda, float beta,
                float* restrict B, CBLAS_INT ldb,
                float* restrict work, CBLAS_INT blocksize)
{
    CBLAS_INT i, j, l;
    char transchara = (transa ? 'T' : 'N');
    float zero = 0.0f;

    // Early return for degenerate cases
    if ((m <= 0) || (n <= 0) || (k <= 0)) { return; }

    // Pre-calculate block structure
    CBLAS_INT num_full_blocks = n / blocksize;
    CBLAS_INT remainder_cols = n % blocksize;

    // Process full blocks
    for (i = 0; i < num_full_blocks; i++)
    {
        CBLAS_INT block_start = i * blocksize;

        // Compute work = alpha * op(A) * B(:, block_start:block_start+blocksize-1)
        BLAS_FUNC(sgemm)(&transchara, "N", &m, &blocksize, &k, &alpha, A, &lda, &B[block_start * ldb], &ldb, &zero, work, &m);

        // Copy result back to B
        float* B_col = &B[block_start * ldb];
        float* work_ptr = work;

        for (j = 0; j < blocksize; j++)
        {
            float* B_ptr = B_col;
            for (l = 0; l < m; l++)
            {
                // If beta != 0, we need to scale the existing B values
                if (beta == 0.0f)
                {
                    *B_ptr++ = *work_ptr++;
                } else {
                    *B_ptr = *work_ptr++ + beta * (*B_ptr);
                    B_ptr++;
                }
            }
            B_col += ldb;  // Move to next column
        }
    }

    // Handle remainders
    if (remainder_cols > 0) {
        CBLAS_INT block_start = num_full_blocks * blocksize;

        BLAS_FUNC(sgemm)(&transchara, "N", &m, &remainder_cols, &k, &alpha, A, &lda, &B[block_start * ldb], &ldb, &zero, work, &m);

        // Copy remainder results back
        float* B_col = &B[block_start * ldb];
        float* work_ptr = work;

        for (j = 0; j < remainder_cols; j++)
        {
            float* B_ptr = B_col;
            for (l = 0; l < m; l++)
            {
                // If beta != 0, we need to scale the existing B values
                if (beta == 0.0f)
                {
                    *B_ptr++ = *work_ptr++;
                } else {
                    *B_ptr = *work_ptr++ + beta * (*B_ptr);
                    B_ptr++;
                }
            }
            B_col += ldb;  // Move to next column
        }
    }
}


void sgemm_ovwr_left(const CBLAS_INT transb, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, float alpha,
                     float* restrict A, CBLAS_INT lda, float* restrict B, CBLAS_INT ldb,
                     float* restrict work, CBLAS_INT blocksize)
{
    CBLAS_INT i, j, l;
    char transcharb = (transb ? 'T' : 'N');
    float zero = 0.0f;

    if (m <= 0 || n <= 0 || k <= 0) { return; }

    // Pre-calculate block structure
    CBLAS_INT num_full_blocks = m / blocksize;
    CBLAS_INT remainder_rows = m % blocksize;

    // Process full blocks of rows
    for (i = 0; i < num_full_blocks; i++)
    {
        CBLAS_INT block_start = i * blocksize;

        // Compute work = alpha * A(block_start:block_start+blocksize-1, :) * op(B)
        BLAS_FUNC(sgemm)("N", &transcharb, &blocksize, &n, &k, &alpha, &A[block_start], &lda, B, &ldb, &zero, work, &blocksize);

        // Copy result back to A
        float* work_ptr = work;

        for (j = 0; j < n; j++)
        {
            float* A_ptr = &A[block_start + j * lda];
            for (l = 0; l < blocksize; l++)
            {
                *A_ptr++ = *work_ptr++;
            }
        }
    }

    // Handle remainder rows
    if (remainder_rows > 0) {
        CBLAS_INT block_start = num_full_blocks * blocksize;

        BLAS_FUNC(sgemm)("N", &transcharb, &remainder_rows, &n, &k, &alpha, &A[block_start], &lda, B, &ldb, &zero, work, &remainder_rows);

        // Copy remainder results back
        float* work_ptr = work;

        for (j = 0; j < n; j++)
        {
            float* A_ptr = &A[block_start + j * lda];
            for (l = 0; l < remainder_rows; l++)
            {
                *A_ptr++ = *work_ptr++;
            }
        }
    }
}


void dgemm_ovwr(
    const int transa, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, double alpha, double* restrict A, CBLAS_INT lda, double beta,
    double* restrict B, CBLAS_INT ldb, double* restrict work, CBLAS_INT blocksize)
{
    CBLAS_INT i, j, l;
    char transchara = (transa ? 'T' : 'N');
    double zero = 0.0;

    // Early return for degenerate cases
    if ((m <= 0) || (n <= 0) || (k <= 0)) { return; }

    // Pre-calculate block structure
    CBLAS_INT num_full_blocks = n / blocksize;
    CBLAS_INT remainder_cols = n % blocksize;

    // Process full blocks - beta = 0 case
    for (i = 0; i < num_full_blocks; i++)
    {
        CBLAS_INT block_start = i * blocksize;

        // Compute work = alpha * op(A) * B(:, block_start:block_start+blocksize-1)
        BLAS_FUNC(dgemm)(&transchara, "N", &m, &blocksize, &k, &alpha, A, &lda, &B[block_start * ldb], &ldb, &zero, work, &m);

        // Copy result back to B
        double* B_col = &B[block_start * ldb];
        double* work_ptr = work;

        for (j = 0; j < blocksize; j++)
        {
            double* B_ptr = B_col;
            for (l = 0; l < m; l++)
            {
                // If beta != 0, we need to scale the existing B values
                if (beta == 0.0)
                {
                    *B_ptr++ = *work_ptr++;
                } else {
                    *B_ptr = *work_ptr++ + beta * (*B_ptr);
                    B_ptr++;
                }
            }
            B_col += ldb;  // Move to next column
        }
    }

    // Handle remainder - beta = 0 case
    if (remainder_cols > 0) {
        CBLAS_INT block_start = num_full_blocks * blocksize;

        BLAS_FUNC(dgemm)(&transchara, "N", &m, &remainder_cols, &k, &alpha, A, &lda, &B[block_start * ldb], &ldb, &zero, work, &m);

        // Copy remainder results back
        double* B_col = &B[block_start * ldb];
        double* work_ptr = work;

        for (j = 0; j < remainder_cols; j++)
        {
            double* B_ptr = B_col;
            for (l = 0; l < m; l++)
            {
                // If beta != 0, we need to scale the existing B values
                if (beta == 0.0)
                {
                    *B_ptr++ = *work_ptr++;
                } else {
                    *B_ptr = *work_ptr++ + beta * (*B_ptr);
                    B_ptr++;
                }
            }
            B_col += ldb;  // Move to next column
        }
    }
}


void dgemm_ovwr_left(const CBLAS_INT transb, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, double alpha,
                     double* restrict A, CBLAS_INT lda, double* restrict B, CBLAS_INT ldb,
                     double* restrict work, CBLAS_INT blocksize)
{
    CBLAS_INT i, j, l;
    char transcharb = (transb ? 'T' : 'N');
    double zero = 0.0;

    // Early return for degenerate cases
    if (m <= 0 || n <= 0 || k <= 0) return;

    // Pre-calculate block structure
    CBLAS_INT num_full_blocks = m / blocksize;
    CBLAS_INT remainder_rows = m % blocksize;

    // Process full blocks of rows
    for (i = 0; i < num_full_blocks; i++) {
        CBLAS_INT block_start = i * blocksize;

        // Compute work = alpha * A(block_start:block_start+blocksize-1, :) * op(B)
        BLAS_FUNC(dgemm)("N", &transcharb, &blocksize, &n, &k, &alpha, &A[block_start], &lda, B, &ldb, &zero, work, &blocksize);

        // Copy result back to A
        double* work_ptr = work;
        for (j = 0; j < n; j++)
        {
            double* A_ptr = &A[block_start + j * lda];
            for (l = 0; l < blocksize; l++)
            {
                *A_ptr++ = *work_ptr++;
            }
        }
    }

    // Handle remainder rows
    if (remainder_rows > 0) {
        CBLAS_INT block_start = num_full_blocks * blocksize;

        BLAS_FUNC(dgemm)("N", &transcharb, &remainder_rows, &n, &k, &alpha, &A[block_start], &lda, B, &ldb, &zero, work, &remainder_rows);

        // Copy remainder results back
        double* work_ptr = work;
        for (j = 0; j < n; j++)
        {
            double* A_ptr = &A[block_start + j * lda];
            for (l = 0; l < remainder_rows; l++)
            {
                *A_ptr++ = *work_ptr++;
            }
        }
    }
}


/**
 * @brief Naive (complex float) x (float) GEMM kernel.
 *
 * @param transb Whether to transpose B (ignored, always A * B^T)
 * @param m Number of rows of A and C
 * @param n Number of columns of B and C
 * @param k Number of columns of A and rows of B
 * @param A complex float matrix A of size (m, k)
 * @param lda Leading dimension of A
 * @param B float matrix B of size (n, k)
 * @param ldb Leading dimension of B
 * @param C complex float output matrix C of size (m, n)
 * @param ldc Leading dimension of C
 *
 * @note Fortran code is a bit confusing over this kernel since they also provide a
 *       blocked version however left a note that the performance is not promising.
 *
 */
static void csgemm_kernel(
    const CBLAS_INT transb, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, const PROPACK_CPLXF_TYPE* restrict A, CBLAS_INT lda,
    const float* restrict B, CBLAS_INT ldb, PROPACK_CPLXF_TYPE* restrict C, CBLAS_INT ldc)
{
    (void)transb;  // Unused parameter
    // Initialize C to zero
    for (CBLAS_INT j = 0; j < n; j++)
    {
        for (CBLAS_INT i = 0; i < m; i++)
        {
            C[i + j * ldc] = PROPACK_cplxf(0.0f, 0.0f );
        }
    }
    // Always compute C = A * B^T
    for (CBLAS_INT l = 0; l < k; l++)
    {
        for (CBLAS_INT j = 0; j < n; j++)
        {
            float b_val = B[j + l * ldb];
            for (CBLAS_INT i = 0; i < m; i++)
            {
#ifdef _MSC_VER
                PROPACK_CPLXF_TYPE tmp1 = _FCmulcr(A[i + l * lda], b_val);
                PROPACK_CPLXF_TYPE tmp2 = C[i + j * ldc];
                C[i + j * ldc] = PROPACK_cplxf(crealf(tmp2) + crealf(tmp1), cimagf(tmp2) + cimagf(tmp1));
#else
                C[i + j * ldc] += A[i + l * lda] * b_val;
#endif
            }
        }
    }
}


// GEMM operation for mixed dtype, float complex x float
void csgemm_ovwr_left(
    const CBLAS_INT transb, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, PROPACK_CPLXF_TYPE* restrict A, CBLAS_INT lda,
    const float* restrict B, CBLAS_INT ldb, PROPACK_CPLXF_TYPE* restrict work, CBLAS_INT blocksize)
{
    if ((m <= 0) || (n <= 0) || (k <= 0)) { return; }

    CBLAS_INT num_full_blocks = m / blocksize;
    CBLAS_INT remainder_rows = m % blocksize;

    // Process full blocks
    for (CBLAS_INT i = 0; i < num_full_blocks; i++)
    {
        CBLAS_INT block_start = i * blocksize;
        // Compute: work = A[block_start:block_start+blocksize-1, :] * op(B)
        csgemm_kernel(transb, blocksize, n, k, &A[block_start], lda, B, ldb, work, blocksize);

        // Copy result back to A
        for (CBLAS_INT j = 0; j < n; j++)
        {
            for (CBLAS_INT l = 0; l < blocksize; l++)
            {
                A[block_start + l + j * lda] = work[l + j * blocksize];
            }
        }
    }

    // Handle remainder rows
    if (remainder_rows > 0)
    {
        CBLAS_INT block_start = num_full_blocks * blocksize;
        csgemm_kernel(transb, remainder_rows, n, k, &A[block_start], lda, B, ldb, work, remainder_rows);
        for (CBLAS_INT j = 0; j < n; j++)
        {
            for (CBLAS_INT l = 0; l < remainder_rows; l++)
            {
                A[block_start + l + j * lda] = work[l + j * remainder_rows];
            }
        }
    }
}


/**
 * @brief Naive (complex double) x (double) GEMM kernel.
 *
 * @param transb Whether to transpose B (ignored, always A * B^T)
 * @param m Number of rows of A and C
 * @param n Number of columns of B and C
 * @param k Number of columns of A and rows of B
 * @param A Complex double matrix A of size (m, k)
 * @param lda Leading dimension of A
 * @param B Double matrix B of size (n, k)
 * @param ldb Leading dimension of B
 * @param C Complex double output matrix C of size (m, n)
 * @param ldc Leading dimension of C
 *
 * @note Fortran code is a bit confusing over this kernel since they also provide a
 *       blocked version however left a note that the performance is not promising.
 *
 */
static void zdgemm_kernel(
    const CBLAS_INT transb, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, const PROPACK_CPLX_TYPE* restrict A, CBLAS_INT lda,
    const double* restrict B, CBLAS_INT ldb, PROPACK_CPLX_TYPE* restrict C, CBLAS_INT ldc)
{
    (void)transb;  // Unused parameter
    // Initialize C to zero
    for (CBLAS_INT j = 0; j < n; j++)
    {
        for (CBLAS_INT i = 0; i < m; i++)
        {
            C[i + j * ldc] = PROPACK_cplx(0.0, 0.0);
        }
    }
    // Always compute C = A * B^T
    for (CBLAS_INT l = 0; l < k; l++)
    {
        for (CBLAS_INT j = 0; j < n; j++)
        {
            double b_val = B[j + l * ldb];
            for (CBLAS_INT i = 0; i < m; i++)
            {
#ifdef _MSC_VER
                PROPACK_CPLX_TYPE tmp1 = _Cmulcr(A[i + l * lda], b_val);
                PROPACK_CPLX_TYPE tmp2 = C[i + j * ldc];
                C[i + j * ldc] = PROPACK_cplx(creal(tmp2) + creal(tmp1), cimag(tmp2) + cimag(tmp1));
#else
                C[i + j * ldc] += A[i + l * lda] * b_val;
#endif
            }
        }
    }
}


// GEMM operation for mixed dtype, complex double x double
void zdgemm_ovwr_left(
    const CBLAS_INT transb, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, PROPACK_CPLX_TYPE* restrict A, CBLAS_INT lda,
    const double* restrict B, CBLAS_INT ldb, PROPACK_CPLX_TYPE* restrict work, CBLAS_INT blocksize)
{
    if ((m <= 0) || (n <= 0) || (k <= 0)) { return; }

    CBLAS_INT num_full_blocks = m / blocksize;
    CBLAS_INT remainder_rows = m % blocksize;

    // Process full blocks
    for (CBLAS_INT i = 0; i < num_full_blocks; i++)
    {
        CBLAS_INT block_start = i * blocksize;
        // Compute: work = A[block_start:block_start+blocksize-1, :] * op(B)
        zdgemm_kernel(transb, blocksize, n, k, &A[block_start], lda, B, ldb, work, blocksize);

        // Copy result back to A
        for (CBLAS_INT j = 0; j < n; j++)
        {
            for (CBLAS_INT l = 0; l < blocksize; l++)
            {
                A[block_start + l + j * lda] = work[l + j * blocksize];
            }
        }
    }

    // Handle remainder rows
    if (remainder_rows > 0)
    {
        CBLAS_INT block_start = num_full_blocks * blocksize;
        zdgemm_kernel(transb, remainder_rows, n, k, &A[block_start], lda, B, ldb, work, remainder_rows);
        for (CBLAS_INT j = 0; j < n; j++)
        {
            for (CBLAS_INT l = 0; l < remainder_rows; l++)
            {
                A[block_start + l + j * lda] = work[l + j * remainder_rows];
            }
        }
    }
}
