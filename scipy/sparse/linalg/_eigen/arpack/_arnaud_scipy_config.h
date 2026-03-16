/*
 * SciPy-specific BLAS/LAPACK configuration for the arnaud library.
 *
 * This header is injected via -include by meson.build so that arnaud's
 * ARNAUD_BLAS() macro uses npy_cblas.h's BLAS_FUNC() for correct symbol
 * mangling on all platforms (MKL, OpenBLAS, Accelerate, etc.).
 *
 * The arnaud library itself never includes this file directly.
 */
#include "npy_cblas.h"
#define ARNAUD_BLAS(name) BLAS_FUNC(name)
