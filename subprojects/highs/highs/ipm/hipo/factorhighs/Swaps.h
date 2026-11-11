#ifndef FACTORHIGHS_SWAPS_H
#define FACTORHIGHS_SWAPS_H

#include "DataCollector.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

void permuteWithSwaps(double* x, const Int* swaps, Int n, bool reverse = false);

void swapCols(char uplo, Int n, double* A, Int lda, Int i, Int j, Int* swaps,
              Int* sign, DataCollector& data);

void applySwaps(const Int* swaps, Int nrow, Int ncol, double* R,
                DataCollector& data);

}  // namespace hipo

#endif