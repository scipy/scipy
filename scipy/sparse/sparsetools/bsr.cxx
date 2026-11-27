#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_sparsetools_ARRAY_API

#include "sparsetools.h"
#include "bsr.h"

#define SPTOOLS_BSR_ROUTINES(X_0, X_I, X_IT) \
  X_IT(bsr_diagonal)                         \
  X_IT(bsr_tocsr)                            \
  X_IT(bsr_scale_rows)                       \
  X_IT(bsr_scale_columns)                    \
  X_IT(bsr_sort_indices)                     \
  X_IT(bsr_transpose)                        \
  X_IT(bsr_matmat)                           \
  X_IT(bsr_matvec)                           \
  X_IT(bsr_matvecs)                          \
  X_IT(bsr_elmul_bsr)                        \
  X_IT(bsr_eldiv_bsr)                        \
  X_IT(bsr_plus_bsr)                         \
  X_IT(bsr_minus_bsr)                        \
  X_IT(bsr_maximum_bsr)                      \
  X_IT(bsr_minimum_bsr)                      \
  X_IT(bsr_ne_bsr)                           \
  X_IT(bsr_lt_bsr)                           \
  X_IT(bsr_gt_bsr)                           \
  X_IT(bsr_le_bsr)                           \
  X_IT(bsr_ge_bsr)

SPTOOLS_DEFINE_UNIT(bsr, SPTOOLS_BSR_ROUTINES)
