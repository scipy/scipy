/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsExtrasExternalDeps.h
 * @brief Defines the set of external features available
 */
#ifndef HIGHS_EXTRAS_EXTERNAL_DEPS_H_
#define HIGHS_EXTRAS_EXTERNAL_DEPS_H_

#include <tuple>

#include "HighsExtrasApiBinding.h"
#include "HighsMetis.h"
#include "amd/amd.h"
#include "blas/mycblas.h"
#include "rcm/rcm.h"

//
// Feature API for import/export dynamic loading
// Each feature needs 3 layers:
//   1. internal functions
//   2. API import/exports (struct of function pointers + binding)
//   3. wrapper class (to access API via static methods)
//
// We use a template wrapper to simplify the API import/export (i.e., creating
// struct and binding to the internal functions), however, we write an explicit
// wrapper for better intellisense etc.
//

//
// 0. Define Features
//
namespace HighsExtras {
struct extras_family {
  static const HighsExtrasFeatureInfo* getInfo();
};

template <>
struct wrapper_storage<extras_family> {
  template <class Methods>
  static feature_api<Methods>& getApi();

  static const HighsExtrasFeatureInfo* getInfo();
};

// feature names need to be available at compile-time for any consumers
template <int Index>
struct extras_feature : feature_base<extras_family, Index> {
  static const char* name() {
    return std::get<Index>(std::make_tuple("amd", "blas", "metis", "rcm"));
  }
};

//
// 1. Methods
//

using amd_methods =
    std::tuple<decltype(&Highs_amd_defaults), decltype(&Highs_amd_order)>;

using blas_methods = std::tuple<
    decltype(&cblas_daxpy), decltype(&cblas_dcopy), decltype(&cblas_dscal),
    decltype(&cblas_dswap), decltype(&cblas_dgemv), decltype(&cblas_dtpsv),
    decltype(&cblas_dtrsv), decltype(&cblas_dger), decltype(&cblas_dgemm),
    decltype(&cblas_dsyrk), decltype(&cblas_dtrsm),
    decltype(&highs_openblas_set_num_threads)>;

using metis_methods = std::tuple<decltype(&Highs_METIS_SetDefaultOptions),
                                 decltype(&Highs_METIS_NodeND)>;

using rcm_methods = std::tuple<decltype(&Highs_genrcm)>;

//
// 2. API import/export structure
//

// define the struct of function pointers for the HighsExtrasApi
struct HighsExtrasApi : feature_api<amd_methods>,
                        feature_api<metis_methods>,
                        feature_api<blas_methods>,
                        feature_api<rcm_methods> {
  template <class Methods>
  feature_api<Methods>& as() {
    return static_cast<feature_api<Methods>&>(*this);
  }

  template <class Methods>
  const feature_api<Methods>& as() const {
    return static_cast<const feature_api<Methods>&>(*this);
  }
};

//
// 3. Wrapper access to methods
//
struct amd : extras_feature<0> {
  using impl = feature_wrapper<extras_family, amd_methods>;

  static void set_defaults(double Control[]) {
    impl::template fn<0>()(Control);
  }

  static int order(amd_int n, const amd_int Ap[], const amd_int Ai[],
                   amd_int P[], double Control[], double Info[]) {
    return impl::template fn<1>()(n, Ap, Ai, P, Control, Info);
  }
};

struct blas : extras_feature<1> {
  using impl = feature_wrapper<extras_family, blas_methods>;

  static void daxpy(blasint n, double da, const double* dx, blasint incx,
                    double* dy, blasint incy) {
    impl::template fn<0>()(n, da, dx, incx, dy, incy);
  }

  static void dcopy(blasint n, const double* dx, blasint incx, double* dy,
                    blasint incy) {
    impl::template fn<1>()(n, dx, incx, dy, incy);
  }

  static void dscal(blasint n, const double da, double* dx, blasint incx) {
    impl::template fn<2>()(n, da, dx, incx);
  }

  static void dswap(blasint n, double* dx, blasint incx, double* dy,
                    blasint incy) {
    impl::template fn<3>()(n, dx, incx, dy, incy);
  }

  static void dgemv(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE transa,
                    blasint M, blasint n, double alpha, const double* A,
                    blasint lda, const double* x, blasint incx, double beta,
                    double* y, blasint incy) {
    impl::template fn<4>()(order, transa, M, n, alpha, A, lda, x, incx, beta, y,
                           incy);
  }

  static void dtpsv(enum CBLAS_ORDER order, enum CBLAS_UPLO uplo,
                    enum CBLAS_TRANSPOSE transa, enum CBLAS_DIAG diag,
                    blasint n, const double* ap, double* x, blasint incx) {
    impl::template fn<5>()(order, uplo, transa, diag, n, ap, x, incx);
  }

  static void dtrsv(enum CBLAS_ORDER order, enum CBLAS_UPLO uplo,
                    enum CBLAS_TRANSPOSE transa, enum CBLAS_DIAG diag,
                    blasint n, const double* A, blasint lda, double* x,
                    blasint incx) {
    impl::template fn<6>()(order, uplo, transa, diag, n, A, lda, x, incx);
  }

  static void dger(enum CBLAS_ORDER order, blasint m, blasint n, double alpha,
                   const double* x, blasint incx, const double* y, blasint incy,
                   double* A, blasint lda) {
    impl::template fn<7>()(order, m, n, alpha, x, incx, y, incy, A, lda);
  }

  static void dgemm(enum CBLAS_ORDER order, enum CBLAS_TRANSPOSE transa,
                    enum CBLAS_TRANSPOSE transb, blasint m, blasint n,
                    blasint k, double alpha, const double* A, blasint lda,
                    const double* B, blasint ldb, double beta, double* C,
                    blasint ldc) {
    impl::template fn<8>()(order, transa, transb, m, n, k, alpha, A, lda, B,
                           ldb, beta, C, ldc);
  }

  static void dsyrk(enum CBLAS_ORDER order, enum CBLAS_UPLO uplo,
                    enum CBLAS_TRANSPOSE trans, blasint n, blasint k,
                    double alpha, const double* A, blasint lda, double beta,
                    double* C, blasint ldc) {
    impl::template fn<9>()(order, uplo, trans, n, k, alpha, A, lda, beta, C,
                           ldc);
  }

  static void dtrsm(enum CBLAS_ORDER order, enum CBLAS_SIDE side,
                    enum CBLAS_UPLO uplo, enum CBLAS_TRANSPOSE transa,
                    enum CBLAS_DIAG diag, blasint m, blasint n, double alpha,
                    const double* A, blasint lda, double* B, blasint ldb) {
    impl::template fn<10>()(order, side, uplo, transa, diag, m, n, alpha, A,
                            lda, B, ldb);
  }

  static int openblas_set_num_threads(int num_threads) {
    return impl::template fn<11>()(num_threads);
  }
};

struct metis : extras_feature<2> {
  using impl = feature_wrapper<extras_family, metis_methods>;

  static int set_default_options(idx_t options[]) {
    return impl::template fn<0>()(options);
  }

  static int nodeND(idx_t* nvtxs, const idx_t* xadj, const idx_t* adjncy,
                    idx_t* vwgt, idx_t* options, idx_t* perm, idx_t* iperm) {
    return impl::template fn<1>()(nvtxs, xadj, adjncy, vwgt, options, perm,
                                  iperm);
  }
};

struct rcm : extras_feature<3> {
  using impl = feature_wrapper<extras_family, rcm_methods>;

  static int genrcm(HighsInt node_num, HighsInt adj_num,
                    const HighsInt adj_row[], const HighsInt adj[],
                    HighsInt perm[]) {
    return impl::template fn<0>()(node_num, adj_num, adj_row, adj, perm);
  }
};

//
// define feature set
//

using hipo = require<amd, blas, metis, rcm>;
using extrasAll = require<hipo>;

}  // namespace HighsExtras

#endif  // HIGHS_EXTRAS_EXTERNAL_DEPS_H_