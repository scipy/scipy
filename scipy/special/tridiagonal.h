/* Compute eigenvectors and eigenvalues of symmetric tridiagonal matrices for use in
 * special function kernels. Wraps LAPACK dstevd.
 */

#pragma once

#include <vector>
#include <xsf/error.h>

#include "blaslapack_declarations.h"

namespace special {

struct eigvalsh_tridiagonal {
    std::vector<double> Z{0.0};
    std::vector<double> work;
    std::vector<CBLAS_INT> iwork{0};
    sf_error_t operator()(std::vector<double> &D, std::vector<double> &E) {
        auto N = static_cast<CBLAS_INT>(D.size());

        if (N == 0) {
            return SF_ERROR_OK;
        }

        char jobz = 'N';
        CBLAS_INT ldz = 1;
        CBLAS_INT info;

        /* Allocate the optimal workspace */
        CBLAS_INT lwork = 1 + 2 * N;
        CBLAS_INT liwork = 1;
        work.resize(lwork);
        BLAS_FUNC(dstevd)
        (&jobz, &N, D.data(), E.data(), Z.data(), &ldz, work.data(), &lwork, iwork.data(), &liwork, &info);

        if (info < 0) {
            return SF_ERROR_OTHER;
        }

        if (info > 0) {
            return SF_ERROR_NO_RESULT;
        }

        return SF_ERROR_OK;
    }
};

struct eigh_tridiagonal {
    std::vector<double> work;
    std::vector<CBLAS_INT> iwork;

    sf_error_t operator()(std::vector<double> &D, std::vector<double> &E, std::vector<double> &Z) {
        auto N = static_cast<CBLAS_INT>(D.size());
        if (N == 0) {
            return SF_ERROR_OK;
        }

        char jobz = 'V';
        CBLAS_INT ldz = N;
        CBLAS_INT info = 0;

        /* Allocate the optimal workspace */
        CBLAS_INT lwork = 1 + 4 * N + N * N;
        CBLAS_INT liwork = 3 + 5 * N;
        work.resize(lwork);
        iwork.resize(liwork);

        BLAS_FUNC(dstevd)
        (&jobz, &N, D.data(), E.data(), Z.data(), &ldz, work.data(), &lwork, iwork.data(), &liwork, &info);

        if (info < 0) {
            return SF_ERROR_OTHER;
        }

        if (info > 0) {
            return SF_ERROR_NO_RESULT;
        }

        return SF_ERROR_OK;
    }
};
} // namespace special
