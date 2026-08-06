/* Compute eigenvectors and eigenvalues of symmetric tridiagonal matrices for
 * use in special function kernels. Wraps LAPACK dstevd.
 */

#pragma once

#include <new>
#include <vector>
#include <xsf/error.h>

#include "blaslapack_declarations.h"

namespace special {

struct eigvalsh_tridiagonal {
  std::vector<double> work;
  std::vector<CBLAS_INT> iwork;
  std::vector<CBLAS_INT> isuppz{2};
  std::vector<double> W;

  sf_error_t operator()(std::vector<double> &D, std::vector<double> &E,
                        decltype(D.size()) idx, double &eigenvalue) {
    if (idx >= D.size()) {
      return SF_ERROR_OTHER;
    }

    auto N = static_cast<CBLAS_INT>(D.size());

    if (N == 0) {
      return SF_ERROR_OK;
    }

    char jobz = 'N';
    char range = 'I';
    CBLAS_INT ldz = 1;
    CBLAS_INT info = 0;

    /* Dummy bounds since we find by index. */
    double vl = 0.0;
    double vu = 0.0;

    auto il = static_cast<CBLAS_INT>(idx) + 1;
    CBLAS_INT iu = il;

    char cmach = 'S';
    double abstol = BLAS_FUNC(dlamch)(&cmach);
    CBLAS_INT m = 0;

    try {
      W.resize(N);
    } catch (const std::bad_alloc &) {
      return SF_ERROR_MEMORY;
    }

    /* Dummy Z since we aren't finding eigenvectors. */
    double Z = 0.0;

    /* Allocate the optimal workspace */
    CBLAS_INT lwork = -1;
    CBLAS_INT liwork = -1;
    double work_query = 0.0;
    CBLAS_INT iwork_query = 0;
    BLAS_FUNC(dstevr)
    (&jobz, &range, &N, D.data(), E.data(), &vl, &vu, &il, &iu, &abstol, &m,
     W.data(), &Z, &ldz, isuppz.data(), &work_query, &lwork, &iwork_query,
     &liwork, &info);

    if (info != 0) {
      return SF_ERROR_OTHER;
    }

    lwork = static_cast<CBLAS_INT>(work_query);
    liwork = iwork_query;
    try {
      work.resize(lwork);
      iwork.resize(liwork);
    } catch (const std::bad_alloc &) {
      return SF_ERROR_MEMORY;
    }

    BLAS_FUNC(dstevr)
    (&jobz, &range, &N, D.data(), E.data(), &vl, &vu, &il, &iu, &abstol, &m,
     W.data(), &Z, &ldz, isuppz.data(), work.data(), &lwork, iwork.data(),
     &liwork, &info);

    if (info < 0) {
      return SF_ERROR_OTHER;
    }

    if (info > 0) {
      return SF_ERROR_NO_RESULT;
    }

    eigenvalue = W[0];
    return SF_ERROR_OK;
  }
};

struct eigh_tridiagonal {
  std::vector<double> work;
  std::vector<CBLAS_INT> iwork;
  std::vector<CBLAS_INT> isuppz{2};
  std::vector<double> W;

  sf_error_t operator()(std::vector<double> &D, std::vector<double> &E,
                        decltype(D.size()) idx, std::vector<double> &Z,
                        double &eigenvalue) {

    if (idx >= D.size()) {
      return SF_ERROR_OTHER;
    }

    auto N = static_cast<CBLAS_INT>(D.size());

    if (N == 0) {
      return SF_ERROR_OK;
    }

    char jobz = 'V';
    char range = 'I';

    double vl = 0.0;
    double vu = 0.0;
    CBLAS_INT il = static_cast<CBLAS_INT>(idx) + 1;
    CBLAS_INT iu = il;

    char cmach = 'S';
    double abstol = BLAS_FUNC(dlamch)(&cmach);
    CBLAS_INT m = 0;

    CBLAS_INT ldz = N;
    CBLAS_INT info = 0;

    try {
      W.resize(N);
      Z.resize(N);
    } catch (const std::bad_alloc &) {
      return SF_ERROR_MEMORY;
    }

    /* Allocate the optimal workspace */
    CBLAS_INT lwork = -1;
    CBLAS_INT liwork = -1;
    double work_query = 0.0;
    CBLAS_INT iwork_query = 0;

    BLAS_FUNC(dstevr)
    (&jobz, &range, &N, D.data(), E.data(), &vl, &vu, &il, &iu, &abstol, &m,
     W.data(), Z.data(), &ldz, isuppz.data(), &work_query, &lwork, &iwork_query,
     &liwork, &info);

    if (info != 0) {
      return SF_ERROR_OTHER;
    }

    lwork = static_cast<CBLAS_INT>(work_query);
    liwork = iwork_query;

    try {
      work.resize(lwork);
      iwork.resize(liwork);
    } catch (const std::bad_alloc &) {
      return SF_ERROR_MEMORY;
    }

    BLAS_FUNC(dstevr)
    (&jobz, &range, &N, D.data(), E.data(), &vl, &vu, &il, &iu, &abstol, &m,
     W.data(), Z.data(), &ldz, isuppz.data(), work.data(), &lwork, iwork.data(),
     &liwork, &info);

    if (info < 0) {
      return SF_ERROR_OTHER;
    }

    if (info > 0) {
      return SF_ERROR_NO_RESULT;
    }

    eigenvalue = W[0];
    return SF_ERROR_OK;
  }
};
} // namespace special
