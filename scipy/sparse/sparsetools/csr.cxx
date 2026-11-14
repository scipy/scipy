#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_sparsetools_ARRAY_API

#include "sparsetools.h"
#include "csr.h"

extern "C" {
#include "csr_impl.h"
}

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
