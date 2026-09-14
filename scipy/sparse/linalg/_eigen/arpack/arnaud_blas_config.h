/*
 * SciPy-specific BLAS/LAPACK configuration for the arnaud library.
 *
 * This header is used when we pass the `ARNAUD_HAS_BLAS_CONFIG` define
 * when building the arnaud sources.
 *
 * The ARNAUD_BLAS() macro uses npy_cblas.h's BLAS_FUNC() for correct symbol
 * mangling on all platforms (MKL, OpenBLAS, Accelerate, etc.).
 *
 * It is conditionally included by blaslapack_declarations.h when
 * ARNAUD_HAS_BLAS_CONFIG is defined.
 */
#include "scipy_blas_defines.h"
#define ARNAUD_BLAS(name) BLAS_FUNC(name)
