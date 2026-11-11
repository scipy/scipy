#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "DenseFact.h"
#include "DgemmParallel.h"
#include "ReturnValues.h"
#include "Swaps.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"

namespace hipo {

// Factorisation with "hybrid formats".

Int denseFactFH(char format, Int n, Int k, double* A, double* B,
                const Int* pivot_sign, double thresh, double* totalreg,
                Int* swaps, double* pivot_2x2, DataCollector& data,
                const FHoptions& options) {
  // ===========================================================================
  // Partial blocked factorisation
  // Matrix A is in format FH
  // Matrix B is in format FP, if format == 'P'
  //                       FH, if format == 'H'
  // BLAS calls: dcopy, dscal, daxpy, dgemm, dtrsm
  // ===========================================================================

  HIPO_CLOCK_CREATE;

  // check input
  if (n < 0 || k < 0 || !A || (k < n && !B)) return kRetInvalidInput;

  // quick return
  if (n == 0) return kRetOk;

  const Int nb = options.nb;

  // number of blocks of columns
  const Int n_blocks = (k - 1) / nb + 1;

  // start of diagonal blocks
  std::vector<Int64> diag_start(n_blocks);
  getDiagStart(n, k, nb, n_blocks, diag_start);

  // size of blocks
  const Int diag_size = nb * nb;
  const Int full_size = nb * nb;

  // buffer for copy of block column
  std::vector<double> T(n * nb);

  // number of rows/columns in the Schur complement
  const Int ns = n - k;

  // number of blocks in Schur complement
  const Int s_blocks = (ns - 1) / nb + 1;

  // buffer for full-format of block of columns of Schur complement
  std::vector<double> schur_buf;
  if (format == 'P') schur_buf.resize(ns * nb);

  // ===========================================================================
  // LOOP OVER BLOCKS
  // ===========================================================================
  for (Int j = 0; j < n_blocks; ++j) {
    // j is the index of the block column

    HIPO_CLOCK_START(2);

    // jb is the number of columns
    const Int jb = std::min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    const Int this_diag_size = jb * jb;
    const Int this_full_size = nb * jb;

    // diagonal block j
    double* D = &A[diag_start[j]];

    // number of rows left below block j
    const Int M = n - nb * j - jb;

    // block of columns below diagonal block j
    const Int64 R_pos = diag_start[j] + this_diag_size;
    double* R = &A[R_pos];

    // ===========================================================================
    // FACTORISE DIAGONAL BLOCK
    // ===========================================================================
    double max_in_R = -1.0;
    if (jb == 1) {
      for (Int64 i = 0; i < M * jb; ++i)
        max_in_R = std::max(max_in_R, std::abs(R[i]));
    }

    double* regul_current = &totalreg[j * nb];
    std::vector<Int> pivot_sign_current(&pivot_sign[j * nb],
                                        &pivot_sign[j * nb] + jb);
    Int* swaps_current = &swaps[j * nb];
    double* pivot_2x2_current = &pivot_2x2[j * nb];
    Int info = denseFactK('U', jb, D, jb, pivot_sign_current.data(), thresh,
                          regul_current, swaps_current, pivot_2x2_current, data,
                          options);
    if (info != 0) return info;

    if (options.pivoting) {
      // swap columns in R
      applySwaps(swaps_current, M, jb, R, data);

      // unswap regularisation, to keep it with original ordering
      permuteWithSwaps(regul_current, swaps_current, jb, true);
    }

    if (M > 0) {
      // ===========================================================================
      // SOLVE COLUMNS
      // ===========================================================================
      // solve block R with D
      callAndTime_dtrsm('L', 'U', 'T', 'U', jb, M, 1.0, D, jb, R, jb, data);

      // make copy of partially solved columns
      callAndTime_dcopy(jb * M, R, 1, T.data(), 1, data);

      // solve block R with pivots
      Int step = 1;
      for (Int col = 0; col < jb; col += step) {
        if (pivot_2x2_current[col] == 0.0) {
          // 1x1 pivots
          step = 1;
          callAndTime_dscal(M, D[col + jb * col], &R[col], jb, data);
        } else {
          // 2x2 pivots
          step = 2;

          // columns affected
          double* c1 = &R[col];
          double* c2 = &R[col + 1];

          // inverse of 2x2 pivot
          double i_d1 = D[col + jb * col];
          double i_d2 = D[col + 1 + jb * (col + 1)];
          double i_off = pivot_2x2_current[col];

          // copy of original col1
          std::vector<double> c1_temp(M);
          callAndTime_dcopy(M, c1, jb, c1_temp.data(), 1, data);

          // solve col and col+1
          callAndTime_dscal(M, i_d1, c1, jb, data);
          callAndTime_daxpy(M, i_off, c2, jb, c1, jb, data);
          callAndTime_dscal(M, i_d2, c2, jb, data);
          callAndTime_daxpy(M, i_off, c1_temp.data(), 1, c2, jb, data);
        }
      }

      // check entries of L
      /*double max_in_R = -1.0;
      for (Int64 i = 0; i < M * jb; ++i) {
        max_in_R = std::max(max_in_R, std::abs(R[i]));
      }
      if (max_in_R > 1e8) printf("%.1e, %5d %5d\n", max_in_R, jb, M);*/

      // ===========================================================================
      // UPDATE FRONTAL
      // ===========================================================================
      Int64 offset{};

      // go through remaining blocks of columns
      for (Int jj = j + 1; jj < n_blocks; ++jj) {
        // number of columns in block jj
        const Int col_jj = std::min(nb, k - nb * jj);

        // number of rows in block jj
        const Int row_jj = n - nb * jj;

        const double* P = &T[offset];
        double* Q = &A[diag_start[jj]];
        const double* Rjj = &R[offset];

        // perform gemm (potentially) in parallel
        if (options.parallel_node)
          dgemmParallel(P, Rjj, Q, col_jj, jb, row_jj, nb, 1.0, data);
        else
          callAndTime_dgemm('T', 'N', col_jj, row_jj, jb, -1.0, P, jb, Rjj, jb,
                            1.0, Q, col_jj, data);

        offset += jb * col_jj;
      }
      HIPO_CLOCK_STOP(2, data, kTimeDenseFact_main);

      // ===========================================================================
      // UPDATE SCHUR COMPLEMENT
      // ===========================================================================
      HIPO_CLOCK_START(2);
      if (k < n) {
        Int64 B_offset{};

        // go through blocks of columns of the Schur complement
        for (Int sb = 0; sb < s_blocks; ++sb) {
          // number of rows of the block
          const Int nrow = ns - nb * sb;

          // number of columns of the block
          const Int ncol = std::min(nb, nrow);

          const double* P = &T[offset];
          double* Q = format == 'P' ? schur_buf.data() : &B[B_offset];
          const double* Rjj = &R[offset];

          // beta is 0 to avoid initialising schur_buf if format=='P'
          double beta = format == 'P' ? 0.0 : 1.0;

          // perform gemm (potentially) in parallel
          if (options.parallel_node)
            dgemmParallel(P, Rjj, Q, ncol, jb, nrow, nb, beta, data);
          else
            callAndTime_dgemm('T', 'N', ncol, nrow, jb, -1.0, P, jb, Rjj, jb,
                              beta, Q, ncol, data);

          if (format == 'P') {
            // schur_buf contains Schur complement in hybrid format (with full
            // diagonal blocks). Store it by columns in B (with full diagonal
            // blocks).
            for (Int buf_row = 0; buf_row < nrow; ++buf_row) {
              const Int N = ncol;
              callAndTime_daxpy(N, 1.0, &schur_buf[buf_row * ncol], 1,
                                &B[B_offset + buf_row], nrow, data);
            }
          }

          B_offset += nrow * ncol;
          offset += jb * ncol;
        }
      }
      HIPO_CLOCK_STOP(2, data, kTimeDenseFact_schur);
    }
  }

  return kRetOk;
}

Int denseFactFP2FH(double* A, Int nrow, Int ncol, Int nb, DataCollector& data) {
  // ===========================================================================
  // Packed to Hybrid conversion
  // Matrix A on  input is in format FP
  // Matrix A on output is in format FH
  // BLAS calls: dcopy
  // ===========================================================================

  HIPO_CLOCK_CREATE;

  std::vector<double> buf(nrow * nb);

  Int64 startAtoBuf = 0;
  Int64 startBuftoA = 0;

  for (Int k = 0; k <= (ncol - 1) / nb; ++k) {
    // Number of columns in the block. Can be smaller than nb for last block.
    const Int block_size = std::min(nb, ncol - k * nb);

    // Number of rows in the block
    const Int row_size = nrow - k * nb;

    // Copy block into buf
    callAndTime_dcopy(row_size * block_size, &A[startAtoBuf], 1, buf.data(), 1,
                      data);
    startAtoBuf += row_size * block_size;

    // Copy columns back into A, row by row.
    // One call of dcopy_ for each row of the block of columns.
    for (Int i = 0; i < row_size; ++i) {
      const Int N = block_size;
      callAndTime_dcopy(N, &buf[i], row_size, &A[startBuftoA], 1, data);
      startBuftoA += N;
    }
  }

  HIPO_CLOCK_STOP(2, data, kTimeDenseFact_convert);

  return kRetOk;
}

}  // namespace hipo
