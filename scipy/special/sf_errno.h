/*
 * Define the global error variable for cython_special. Do it in a
 * header file instead of directly in cython_special so that we can
 * use openmp pragmas.
 */

#include "sf_error.h"

sf_error_t sf_errno = SF_ERROR_OK;
#ifdef HAVE_OPENMP
#pragma omp threadprivate(sf_errno)
#endif

