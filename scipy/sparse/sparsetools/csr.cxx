#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_sparsetools_ARRAY_API

#include "sparsetools.h"
#include "csr.h"

extern "C" {
#include "csr_impl.h"
}

#define SPTOOLS_CSR_DEFINE_TEMPLATE(I, T) \
  template void csr_diagonal(const I k, const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], T Yx[]); \
  template void csr_scale_rows(const I n_row, const I n_col, const I Ap[], const I Aj[], T Ax[], const T Xx[]); \
  template void csr_scale_columns(const I n_row, const I n_col, const I Ap[], const I Aj[], T Ax[], const T Xx[]); \
  template void csr_tobsr(const I n_row, const I n_col, const I R, const I C, const I Ap[], const I Aj[], const T Ax[], I Bp[], I Bj[], T Bx[]); \
  template void csr_todense(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], T Bx[]); \
  template void csr_sort_indices(const I n_row, const I Ap[], I Aj[], T Ax[]); \
  template void csr_tocsc(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], I Bp[], I Bi[], T Bx[]); \
  template void csr_toell(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I row_length, I Bj[], T Bx[]); \
  template void csr_matmat(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[]); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const std::not_equal_to<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const std::less<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const std::less_equal<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const std::greater_equal<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const std::multiplies<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const safe_divides<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const std::plus<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const std::minus<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const maximum<T>& op); \
  template void csr_binop_csr(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I Bp[], const I Bj[], const T Bx[], I Cp[], I Cj[], T Cx[], const minimum<T>& op); \
  template void csr_sum_duplicates(const I n_row, const I n_col, I Ap[], I Aj[], T Ax[]); \
  template void csr_eliminate_zeros(const I n_row, const I n_col, I Ap[], I Aj[], T Ax[]); \
  template void csr_matvec(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const T Xx[], T Yx[]); \
  template void csr_matvecs(const I n_row, const I n_col, const I n_vecs, const I Ap[], const I Aj[], const T Ax[], const T Xx[], T Yx[]); \
  template void get_csr_submatrix(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I ir0, const I ir1, const I ic0, const I ic1, std::vector<I>* Bp, std::vector<I>* Bj, std::vector<T>* Bx); \
  template void csr_row_index(const I n_row_idx, const I rows[], const I Ap[], const I Aj[], const T Ax[], I Bj[], T Bx[]); \
  template void csr_row_slice(const I start, const I stop, const I step, const I Ap[], const I Aj[], const T Ax[], I Bj[], T Bx[]); \
  template void csr_column_index2(const I col_order[], const I col_offsets[], const I nnz, const I Aj[], const T Ax[], I Bj[], T Bx[]); \
  template void csr_sample_values(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], const I n_samples, const I Bi[], const I Bj[], T Bx[]);

SPTOOLS_FOR_EACH_INDEX_DATA_TYPE_COMBINATION(SPTOOLS_CSR_DEFINE_TEMPLATE)
