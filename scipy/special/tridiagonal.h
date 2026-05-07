/* Compute eigenvectors and eigenvalues of symmetric tridiagonal matrices for use in 
 * special function kernels. Wraps LAPACK dstevd.
 */


#pragma once

#include <xsf/error.h>
#include <xsf/numpy.h>

#include <vector>

#include "blaslapack_declarations.h"

template <typename DMat, typename EMat>
sf_error_t eigvalsh_tridiagonal(DMat D, EMat E) {
    auto N = static_cast<CBLAS_INT>(D.extent(0));

    if (N == 0) {
        return SF_ERROR_OK;
    }

    char jobz = 'N';
    std::vector<double> Z(1);
    CBLAS_INT ldz = 1;
    CBLAS_INT info;

    /* Allocate the optimal workspace */
    CBLAS_INT lwork = 1 + 2 * N;
    CBLAS_INT liwork = 1;
    std::vector<double> work(lwork);
    std::vector<CBLAS_INT> iwork(liwork);

    BLAS_FUNC(dstevd)
    (&jobz, &N, D.data_handle(), E.data_handle(), Z.data(), &ldz, work.data(), &lwork, iwork.data(), &liwork, &info);

    if (info < 0) {
        return SF_ERROR_OTHER;
    }

    if (info > 0) {
        return SF_ERROR_NO_RESULT;
    }

    return SF_ERROR_OK;
}

template <typename DMat, typename EMat, typename ZMat>
sf_error_t eigh_tridiagonal(DMat D, EMat E, ZMat Z) {
    auto N = static_cast<CBLAS_INT>(D.extent(0));
    if (N == 0) {
        return SF_ERROR_OK;
    }

    char jobz = 'V';
    CBLAS_INT ldz = N;
    CBLAS_INT info = 0;

    // query for optimal workspace
    double work_query;
    CBLAS_INT liwork_query_res;
    CBLAS_INT lwork_q = -1;
    CBLAS_INT liwork_q = -1;

    BLAS_FUNC(dstevd)
    (&jobz, &N, D.data_handle(), E.data_handle(), Z.data_handle(), &ldz, &work_query, &lwork_q, &liwork_query_res,
     &liwork_q, &info);

    if (info != 0) {
        return SF_ERROR_OTHER;
    }

    CBLAS_INT lwork = static_cast<CBLAS_INT>(work_query);
    CBLAS_INT liwork = liwork_query_res;

    std::vector<double> work(lwork);
    std::vector<CBLAS_INT> iwork(liwork);

    BLAS_FUNC(dstevd)
    (&jobz, &N, D.data_handle(), E.data_handle(), Z.data_handle(), &ldz, work.data(), &lwork, iwork.data(), &liwork,
     &info);

    if (info < 0) {
        return SF_ERROR_OTHER;
    }

    if (info > 0) {
        return SF_ERROR_NO_RESULT;
    }

    return SF_ERROR_OK;
}
