#ifndef SCIPY_SLU_CONFIG_H
#define SCIPY_SLU_CONFIG_H

#include <stdlib.h>

/*
 * Support routines
 */
void superlu_python_module_abort(char *msg);
void *superlu_python_module_malloc(size_t size);
void superlu_python_module_free(void *ptr);

#define USER_ABORT  superlu_python_module_abort
#define USER_MALLOC superlu_python_module_malloc
#define USER_FREE   superlu_python_module_free

#define SCIPY_FIX 1

/*
 * Fortran configuration
 */
#ifndef HAVE_BLAS_ILP64

#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define UpCase 1
#else
#define NoChange 1
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#error Uppercase and trailing slash in Fortran names not supported
#else
#define Add_ 1
#endif
#endif

#else

/* Use custom ILP64 BLAS wrappers */

#define sswap_    sswap_superlu_blas_
#define saxpy_    saxpy_superlu_blas_
#define sasum_    sasum_superlu_blas_
#define isamax_   isamax_superlu_blas_
#define scopy_    scopy_superlu_blas_
#define sscal_    sscal_superlu_blas_
#define sger_     sger_superlu_blas_
#define snrm2_    snrm2_superlu_blas_
#define ssymv_    ssymv_superlu_blas_
#define sdot_     sdot_superlu_blas_
#define saxpy_    saxpy_superlu_blas_
#define ssyr2_    ssyr2_superlu_blas_
#define srot_     srot_superlu_blas_
#define sgemv_    sgemv_superlu_blas_
#define strsv_    strsv_superlu_blas_
#define sgemm_    sgemm_superlu_blas_
#define strsm_    strsm_superlu_blas_

#define dswap_    dswap_superlu_blas_
#define daxpy_    daxpy_superlu_blas_
#define dasum_    dasum_superlu_blas_
#define idamax_   idamax_superlu_blas_
#define dcopy_    dcopy_superlu_blas_
#define dscal_    dscal_superlu_blas_
#define dger_     dger_superlu_blas_
#define dnrm2_    dnrm2_superlu_blas_
#define dsymv_    dsymv_superlu_blas_
#define ddot_     ddot_superlu_blas_
#define dsyr2_    dsyr2_superlu_blas_
#define drot_     drot_superlu_blas_
#define dgemv_    dgemv_superlu_blas_
#define dtrsv_    dtrsv_superlu_blas_
#define dgemm_    dgemm_superlu_blas_
#define dtrsm_    dtrsm_superlu_blas_

#define cswap_    cswap_superlu_blas_
#define caxpy_    caxpy_superlu_blas_
#define scasum_   scasum_superlu_blas_
#define icamax_   icamax_superlu_blas_
#define ccopy_    ccopy_superlu_blas_
#define cscal_    cscal_superlu_blas_
#define scnrm2_   scnrm2_superlu_blas_
#define caxpy_    caxpy_superlu_blas_
#define cgemv_    cgemv_superlu_blas_
#define ctrsv_    ctrsv_superlu_blas_
#define cgemm_    cgemm_superlu_blas_
#define ctrsm_    ctrsm_superlu_blas_
#define cgerc_    cgerc_superlu_blas_
#define chemv_    chemv_superlu_blas_
#define cher2_    cher2_superlu_blas_

#define zswap_    zswap_superlu_blas_
#define zaxpy_    zaxpy_superlu_blas_
#define dzasum_   dzasum_superlu_blas_
#define izamax_   izamax_superlu_blas_
#define zcopy_    zcopy_superlu_blas_
#define zscal_    zscal_superlu_blas_
#define dznrm2_   dznrm2_superlu_blas_
#define zaxpy_    zaxpy_superlu_blas_
#define zgemv_    zgemv_superlu_blas_
#define ztrsv_    ztrsv_superlu_blas_
#define zgemm_    zgemm_superlu_blas_
#define ztrsm_    ztrsm_superlu_blas_
#define zgerc_    zgerc_superlu_blas_
#define zhemv_    zhemv_superlu_blas_
#define zher2_    zher2_superlu_blas_

/* LAPACK */
#define dlacon_   dlacon_superlu_blas_
#define slacon_   slacon_superlu_blas_
#define icmax1_   icmax1_superlu_blas_
#define scsum1_   scsum1_superlu_blas_
#define clacon_   clacon_superlu_blas_
#define dzsum1_   dzsum1_superlu_blas_
#define izmax1_   izmax1_superlu_blas_
#define zlacon_   zlacon_superlu_blas_

#endif

#endif
