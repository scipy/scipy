#include "Swaps.h"

#include "CallAndTimeBlas.h"
#include "DataCollector.h"

namespace hipo {

void permuteWithSwaps(double* x, const Int* swaps, Int n, bool reverse) {
  // Apply swaps to vector x of length n

  if (!reverse) {
    // apply the swaps in forward order
    for (Int i = 0; i < n; ++i) {
      if (swaps[i] != i) std::swap(x[i], x[swaps[i]]);
    }
  } else {
    // apply the swaps in backward order
    for (Int i = n - 1; i >= 0; --i) {
      if (swaps[i] != i) std::swap(x[i], x[swaps[i]]);
    }
  }
}

void swapCols(char uplo, Int n, double* A, Int lda, Int i, Int j, Int* swaps,
              Int* sign, DataCollector& data) {
  // Exchange rows/cols i and j of symmetric matrix A

  // make sure that i < j
  if (i == j) return;
  if (i > j) std::swap(i, j);

  // swap diagonal elements
  std::swap(A[i + i * lda], A[j + j * lda]);

  // swap rest of rows/cols
  if (uplo == 'L') {
    callAndTime_dswap(i, &A[i], lda, &A[j], lda, data);
    callAndTime_dswap(n - j - 1, &A[j + 1 + i * lda], 1, &A[j + 1 + j * lda], 1,
                      data);
    callAndTime_dswap(j - i - 1, &A[i + 1 + i * lda], 1, &A[j + (i + 1) * lda],
                      lda, data);
  } else {
    callAndTime_dswap(i, &A[i * lda], 1, &A[j * lda], 1, data);
    callAndTime_dswap(n - j - 1, &A[i + (j + 1) * lda], lda,
                      &A[j + (j + 1) * lda], lda, data);
    callAndTime_dswap(j - i - 1, &A[i + (i + 1) * lda], lda,
                      &A[i + 1 + j * lda], 1, data);
  }

  // swap pivot sign
  std::swap(sign[i], sign[j]);

  // keep track of order of swaps
  swaps[i] = j;

  data.countSwap();
}

void applySwaps(const Int* swaps, Int nrow, Int ncol, double* R,
                DataCollector& data) {
  // apply the column swaps to block R
  for (Int i = 0; i < ncol; ++i) {
    if (swaps[i] != i) {
      // swap col i and col swaps[i]
      callAndTime_dswap(nrow, &R[i], ncol, &R[swaps[i]], ncol, data);
    }
  }
}

}  // namespace hipo