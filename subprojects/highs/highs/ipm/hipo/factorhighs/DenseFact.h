#ifndef FACTORHIGHS_DENSE_FACT_H
#define FACTORHIGHS_DENSE_FACT_H

#include "DataCollector.h"
#include "FactorHighsOptions.h"
#include "FactorHighsSettings.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

/*
  Names:
  denseFact:
  - K : factorisation kernel for diagonal blocks
  - F : blocked factorisation in full format
  - FP: blocked factorisation in format FP
  - FH : blocked factorisation in "hybrid formats"

  Formats used:
  - F : Full format
  - P : lower Packed format
  - FP: lower Packed format with Full diagonal blocks
  - H : lower-blocked-Hybrid format
  - FH: lower-blocked-Hybrid format with Full diagonal blocks

  F, P do not use blocks. FP, H, FH use blocks.
  Blocks are always blocks of columns.
  F, P store by columns.
  FP stores by columns within the blocks. H, FH store by rows within the blocks.
  See report for more details.
*/

// dense factorisation kernel
Int denseFactK(char uplo, Int n, double* A, Int lda, Int* pivot_sign,
               double thresh, double* totalreg, Int* swaps, double* pivot_2x2,
               DataCollector& data, const FHoptions& options);

// dense partial factorisation, in "hybrid formats"
Int denseFactFH(char format, Int n, Int k, double* A, double* B,
                const Int* pivot_sign, double thresh, double* totalreg,
                Int* swaps, double* pivot_2x2, DataCollector& data,
                const FHoptions& options);

// function to convert A from lower packed, to lower-blocked-hybrid format
Int denseFactFP2FH(double* A, Int nrow, Int ncol, Int nb, DataCollector& data);

}  // namespace hipo

#endif