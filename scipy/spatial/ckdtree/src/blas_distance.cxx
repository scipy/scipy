
#include <vector>
#include <algorithm>

#define NPY_NO_DEPRECATED_API NPY_API_VERSION

#include "numpy/arrayobject.h"
#include "npy_cblas.h"

#if defined(__GNUC__)
#define CKDTREE_RESTRICT __restrict__
#elif defined(_MSC_VER)
#define CKDTREE_RESTRICT __restrict
#else
#define CKDTREE_RESTRICT
#endif

// blas
extern "C"
{
    double BLAS_FUNC(ddot)(CBLAS_INT *n, double *x, CBLAS_INT *incx, double *y, CBLAS_INT *incy);
}

#define CKDTREE_BLAS_DIST

#include "ckdtree_decl.h"

double
sqeuclidean_distance_double_blas(const double * CKDTREE_RESTRICT u, 
                                 const double * CKDTREE_RESTRICT v, 
                                 const ckdtree_intp_t n)
{
    CBLAS_INT one = 1;
    ckdtree_intp_t nn = n;
    CBLAS_INT bn = static_cast<CBLAS_INT>(nn); // this is safe because n is the number of dimensions
    double minus_one = -1;
    std::vector<double> tmp(n);
    ckdtree_intp_t i;
    #pragma unroll
    for (i=0; i<n; ++i) {
        tmp[i] = u[i] - v[i];
    }
    return BLAS_FUNC(ddot)(&bn, &tmp[0], &one, &tmp[0], &one);
}


