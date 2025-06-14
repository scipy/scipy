#include "gemm_overwrite.h"


void sgemm_ovwr(const int transa, int m, int n, int k, float alpha,
                const float* restrict A, int lda, float beta,
                float* restrict B, int ldb,
                float* restrict work, int blocksize)
{
    int i, j, l;
    const char transb = 'N';
    const char transchara = (transa ? 'T' : 'N');
    const float zero = 0.0f;

    // Early return for degenerate cases
    if ((m <= 0) || (n <= 0) || (k <= 0)) { return; }

    // Pre-calculate block structure
    int num_full_blocks = n / blocksize;
    int remainder_cols = n % blocksize;

    // Process full blocks
    for (i = 0; i < num_full_blocks; i++)
    {
        int block_start = i * blocksize;

        // Compute work = alpha * op(A) * B(:, block_start:block_start+blocksize-1)
        sgemm_(&transchara, &transb, &m, &blocksize, &k, &alpha, A, &lda, &B[block_start * ldb], &ldb, &zero, work, &m);

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
        int block_start = num_full_blocks * blocksize;

        sgemm_(&transchara, &transb, &m, &remainder_cols, &k, &alpha, A, &lda, &B[block_start * ldb], &ldb, &zero, work, &m);

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


void sgemm_ovwr_left(const int transb, int m, int n, int k, float alpha,
                     float* restrict A, int lda, const float* restrict B, int ldb,
                     float* restrict work, int blocksize)
{
    int i, j, l;
    const char transa = 'N';
    const char transcharb = (transb ? 'T' : 'N');
    const float zero = 0.0f;

    if (m <= 0 || n <= 0 || k <= 0) { return; }

    // Pre-calculate block structure
    int num_full_blocks = m / blocksize;
    int remainder_rows = m % blocksize;

    // Process full blocks of rows
    for (i = 0; i < num_full_blocks; i++) {
        int block_start = i * blocksize;

        // Compute work = alpha * A(block_start:block_start+blocksize-1, :) * op(B)
        sgemm_(&transa, &transcharb, &blocksize, &n, &k, &alpha, &A[block_start], &lda, B, &ldb, &zero, work, &blocksize);

        // Copy result back to A
        float* work_ptr = work;

        for (j = 0; j < n; j++) {
            float* A_ptr = &A[block_start + j * lda];
            for (l = 0; l < blocksize; l++) {
                *A_ptr++ = *work_ptr++;
            }
        }
    }

    // Handle remainder rows
    if (remainder_rows > 0) {
        int block_start = num_full_blocks * blocksize;

        sgemm_(&transa, &transcharb, &remainder_rows, &n, &k, &alpha, &A[block_start], &lda, B, &ldb, &zero, work, &remainder_rows);

        // Copy remainder results back
        float* work_ptr = work;

        for (j = 0; j < n; j++) {
            float* A_ptr = &A[block_start + j * lda];
            for (l = 0; l < remainder_rows; l++) {
                *A_ptr++ = *work_ptr++;
            }
        }
    }
}


void dgemm_ovwr(
    const int transa, int m, int n, int k, double alpha, const double* restrict A, int lda, double beta,
    double* restrict B, int ldb, double* restrict work, int blocksize)
{
    int i, j, l;
    const char transb = 'N';
    const char transchara = (transa ? 'T' : 'N');
    const double zero = 0.0;

    // Early return for degenerate cases
    if ((m <= 0) || (n <= 0) || (k <= 0)) { return; }

    // Pre-calculate block structure
    int num_full_blocks = n / blocksize;
    int remainder_cols = n % blocksize;

    // Process full blocks - beta = 0 case
    for (i = 0; i < num_full_blocks; i++)
    {
        int block_start = i * blocksize;

        // Compute work = alpha * op(A) * B(:, block_start:block_start+blocksize-1)
        dgemm_(&transchara, &transb, &m, &blocksize, &k, &alpha, A, &lda, &B[block_start * ldb], &ldb, &zero, work, &m);

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
        int block_start = num_full_blocks * blocksize;

        dgemm_(&transchara, &transb, &m, &remainder_cols, &k, &alpha, A, &lda, &B[block_start * ldb], &ldb, &zero, work, &m);

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


void dgemm_ovwr_left(const int transb, int m, int n, int k, double alpha,
                     double* restrict A, int lda, const double* restrict B, int ldb,
                     double* restrict work, int blocksize)
{
    int i, j, l;
    const char transa = 'N';
    const char transcharb = (transb ? 'T' : 'N');
    const double zero = 0.0;

    // Early return for degenerate cases
    if (m <= 0 || n <= 0 || k <= 0) return;

    // Pre-calculate block structure
    int num_full_blocks = m / blocksize;
    int remainder_rows = m % blocksize;

    // Process full blocks of rows
    for (i = 0; i < num_full_blocks; i++) {
        int block_start = i * blocksize;

        // Compute work = alpha * A(block_start:block_start+blocksize-1, :) * op(B)
        dgemm_(&transa, &transcharb, &blocksize, &n, &k, &alpha, &A[block_start], &lda, B, &ldb, &zero, work, &blocksize);

        // Copy result back to A - optimized pointer arithmetic
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
        int block_start = num_full_blocks * blocksize;

        dgemm_(&transa, &transcharb, &remainder_rows, &n, &k, &alpha, &A[block_start], &lda, B, &ldb, &zero, work, &remainder_rows);

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


// GEMM operation for mixed dtype, complex float x float
void csgemm_ovwr_left(const int transb, int m, int n, int k,
                      complex float* restrict A, int lda,
                      const float* restrict B, int ldb,
                      complex float* restrict work, int blocksize)
{
    if ((m <= 0) || (n <= 0) || (k <= 0)) { return; }

    int num_full_blocks = m / blocksize;
    int remainder_rows = m % blocksize;

    // Process full blocks
    for (int i = 0; i < num_full_blocks; i++)
    {
        int block_start = i * blocksize;
        // Compute: work = A[block_start:block_start+blocksize-1, :] * op(B)
        csgemm_kernel(transb, blocksize, n, k, &A[block_start], lda, B, ldb, work, blocksize);

        // Copy result back to A
        for (int j = 0; j < n; j++)
        {
            for (int l = 0; l < blocksize; l++)
            {
                A[block_start + l + j * lda] = work[l + j * blocksize];
            }
        }
    }

    // Handle remainder rows
    if (remainder_rows > 0)
    {
        int block_start = num_full_blocks * blocksize;
        csgemm_kernel(transb, remainder_rows, n, k, &A[block_start], lda, B, ldb, work, remainder_rows);
        for (int j = 0; j < n; j++)
        {
            for (int l = 0; l < remainder_rows; l++)
            {
                A[block_start + l + j * lda] = work[l + j * remainder_rows];
            }
        }
    }
}


static void csgemm_kernel(
    const int transb, int m, int n, int k, const complex float* restrict A, int lda,
    const float* restrict B, int ldb, complex float* restrict C, int ldc)
{
    // Initialize C to zero
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < m; i++)
        {
            C[i + j * ldc] = PROPACK_cplxf(0.0f, 0.0f );
        }
    }

    if (transb)
    {
        // C = A * B^T
        for (int l = 0; l < k; l++)
        {
            for (int j = 0; j < n; j++)
            {
                float b_val = B[j + l * ldb];  // B^T[j,l] = B[l,j]
                for (int i = 0; i < m; i++)
                {
#ifdef _MSC_VER
                    PROPACK_CPLXF_TYPE tmp1 = _FCMulcr(A[i + l * lda], b_val);
                    PROPACK_CPLXF_TYPE tmp2 = C[i + j * ldc];
                    C[i + j * ldc] = PROPACK_cplxf(creal(tmp2) + creal(tmp1),
                                                   cimag(tmp2) + cimag(tmp1));
#else
                    C[i + j * ldc] += A[i + l * lda] * b_val;
#endif
                }
            }
        }
    } else {
        // C = A * B (no transpose)
        for (int l = 0; l < k; l++) {
            for (int j = 0; j < n; j++) {
                float b_val = B[l + j * ldb];  // B[l,j]
                for (int i = 0; i < m; i++) {
#ifdef _MSC_VER
                    PROPACK_CPLXF_TYPE tmp1 = _FCMulcr(A[i + l * lda], b_val);
                    PROPACK_CPLXF_TYPE tmp2 = C[i + j * ldc];
                    C[i + j * ldc] = PROPACK_cplxf(creal(tmp2) + creal(tmp1),
                                                   cimag(tmp2) + cimag(tmp1));
#else
                    C[i + j * ldc] += A[i + l * lda] * b_val;
#endif
                }
            }
        }
    }
}


// GEMM operation for mixed dtype, complex double x double
void zdgemm_ovwr_left(const int transb, int m, int n, int k,
                      complex double* restrict A, int lda,
                      const double* restrict B, int ldb,
                      complex double* restrict work, int blocksize)
{
    if ((m <= 0) || (n <= 0) || (k <= 0)) { return; }

    int num_full_blocks = m / blocksize;
    int remainder_rows = m % blocksize;

    // Process full blocks
    for (int i = 0; i < num_full_blocks; i++)
    {
        int block_start = i * blocksize;
        // Compute: work = A[block_start:block_start+blocksize-1, :] * op(B)
        zdgemm_kernel(transb, blocksize, n, k, &A[block_start], lda, B, ldb, work, blocksize);

        // Copy result back to A
        for (int j = 0; j < n; j++)
        {
            for (int l = 0; l < blocksize; l++)
            {
                A[block_start + l + j * lda] = work[l + j * blocksize];
            }
        }
    }

    // Handle remainder rows
    if (remainder_rows > 0)
    {
        int block_start = num_full_blocks * blocksize;
        zdgemm_kernel(transb, remainder_rows, n, k, &A[block_start], lda, B, ldb, work, remainder_rows);
        for (int j = 0; j < n; j++)
        {
            for (int l = 0; l < remainder_rows; l++)
            {
                A[block_start + l + j * lda] = work[l + j * remainder_rows];
            }
        }
    }
}


static void zdgemm_kernel(
    const int transb, int m, int n, int k, const complex double* restrict A, int lda,
    const double* restrict B, int ldb, complex double* restrict C, int ldc)
{
    // Initialize C to zero
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < m; i++)
        {
            C[i + j * ldc] = PROPACK_cplx(0.0, 0.0);
        }
    }

    if (transb)
    {
        // C = A * B^T
        for (int l = 0; l < k; l++)
        {
            for (int j = 0; j < n; j++)
            {
                double b_val = B[j + l * ldb];  // B^T[j,l] = B[l,j]
                for (int i = 0; i < m; i++)
                {
#ifdef _MSC_VER
                    PROPACK_CPLX_TYPE tmp1 = _Cmulcr(A[i + l * lda], b_val);
                    PROPACK_CPLX_TYPE tmp2 = C[i + j * ldc];
                    C[i + j * ldc] = PROPACK_cplx(creal(tmp2) + creal(tmp1),
                                                  cimag(tmp2) + cimag(tmp1));
#else
                    C[i + j * ldc] += A[i + l * lda] * b_val;
#endif
                }
            }
        }
    } else {
        // C = A * B (no transpose)
        for (int l = 0; l < k; l++) {
            for (int j = 0; j < n; j++) {
                double b_val = B[l + j * ldb];  // B[l,j]
                for (int i = 0; i < m; i++) {
#ifdef _MSC_VER
                    PROPACK_CPLX_TYPE tmp1 = _Cmulcr(A[i + l * lda], b_val);
                    PROPACK_CPLX_TYPE tmp2 = C[i + j * ldc];
                    C[i + j * ldc] = PROPACK_cplx(creal(tmp2) + creal(tmp1),
                                                  cimag(tmp2) + cimag(tmp1));
#else
                    C[i + j * ldc] += A[i + l * lda] * b_val;
#endif
                }
            }
        }
    }
}
