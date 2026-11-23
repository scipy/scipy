#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_sparsetools_ARRAY_API

#include "sparsetools.h"
#include "csr.h"

#define SPTOOLS_CSR_ROUTINES(X_0, X_I, X_IT) \
  X_I(csr_matmat_maxnnz)                     \
  X_IT(csr_matmat)                           \
  X_IT(csr_diagonal)                         \
  X_IT(csr_tocsc)                            \
  X_IT(csr_tobsr)                            \
  X_IT(csr_todense)                          \
  X_IT(csr_matvec)                           \
  X_IT(csr_matvecs)                          \
  X_IT(csr_elmul_csr)                        \
  X_IT(csr_eldiv_csr)                        \
  X_IT(csr_plus_csr)                         \
  X_IT(csr_minus_csr)                        \
  X_IT(csr_maximum_csr)                      \
  X_IT(csr_minimum_csr)                      \
  X_IT(csr_ne_csr)                           \
  X_IT(csr_lt_csr)                           \
  X_IT(csr_gt_csr)                           \
  X_IT(csr_le_csr)                           \
  X_IT(csr_ge_csr)                           \
  X_IT(csr_scale_rows)                       \
  X_IT(csr_scale_columns)                    \
  X_IT(csr_sort_indices)                     \
  X_IT(csr_eliminate_zeros)                  \
  X_IT(csr_sum_duplicates)                   \
  X_IT(get_csr_submatrix)                    \
  X_IT(csr_row_index)                        \
  X_IT(csr_row_slice)                        \
  X_I(csr_column_index1)                     \
  X_IT(csr_column_index2)                    \
  X_IT(csr_sample_values)                    \
  X_I(csr_count_blocks)                      \
  X_I(csr_sample_offsets)                    \
  X_IT(csr_hstack)                           \
  X_I(expandptr)                             \
  X_0(test_throw_error)                      \
  X_I(csr_has_sorted_indices)                \
  X_I(csr_has_canonical_format)

SPTOOLS_DEFINE_UNIT(csr, SPTOOLS_CSR_ROUTINES)

/* Functions templated on the index type only. */
#define SPTOOLS_CSR_DEFINE_INDEX_TEMPLATE(I) \
  template bool csr_has_canonical_format(const I n_row, const I Ap[], const I Aj[]);

SPTOOLS_FOR_EACH_INDEX_TYPE(SPTOOLS_CSR_DEFINE_INDEX_TEMPLATE)
#undef SPTOOLS_CSR_DEFINE_INDEX_TEMPLATE

/* Functions templated on the index and data types. */
#define SPTOOLS_CSR_DEFINE_TEMPLATE(I, T) \
  template void csr_sort_indices(const I n_row, const I Ap[], I Aj[], T Ax[]); \
  template void csr_tocsc(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], I Bp[], I Bi[], T Bx[]); \
  template void csr_matmat(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[]); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], npy_bool_wrapper Cx[], const std::not_equal_to<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], npy_bool_wrapper Cx[], const std::less<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], npy_bool_wrapper Cx[], const std::less_equal<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const std::multiplies<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const safe_divides<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const std::plus<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const std::minus<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const maximum<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const minimum<T>& op); \
  template void csr_matvec(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const T Xx[], T Yx[]); \
  template void csr_matvecs(const I n_row, const I n_col, const I n_vecs, const I Ap[], const I Aj[], const T Ax[], const T Xx[], T Yx[]);

SPTOOLS_FOR_EACH_INDEX_DATA_TYPE_COMBINATION(SPTOOLS_CSR_DEFINE_TEMPLATE)
#undef SPTOOLS_CSR_DEFINE_TEMPLATE
