#ifndef SCIPY_SLU_BLAS_CONFIG_H
#define SCIPY_SLU_BLAS_CONFIG_H

#include "scipy_slu_config.h"

/*
 * BLAS symbol name mangling for ILP64.
 * BLAS_SYMBOL_SUFFIX and BLAS_FORTRAN_SUFFIX are passed via -D flags from meson.
 * SLU_BLAS_FUNC(name) produces the correctly-mangled BLAS symbol name,
 * e.g., dgemm_64_ for OpenBLAS ILP64, dgemm_ for LP64.
 */
#ifndef BLAS_SYMBOL_PREFIX
#define BLAS_SYMBOL_PREFIX
#endif
#ifndef BLAS_SYMBOL_SUFFIX
#define BLAS_SYMBOL_SUFFIX
#endif
#ifndef BLAS_FORTRAN_SUFFIX
#define BLAS_FORTRAN_SUFFIX _
#endif

/* Accelerate doesn't use an underscore as suffix, so fix that up here */
#ifdef ACCELERATE_NEW_LAPACK
#undef BLAS_FORTRAN_SUFFIX
#define BLAS_FORTRAN_SUFFIX
#endif

#define SLU_CONCAT4(a,b,c,d) a ## b ## c ## d
#define SLU_EXPAND4(a,b,c,d) SLU_CONCAT4(a,b,c,d)

#ifdef OPENBLAS_ILP64_NAMING_SCHEME
/* OpenBLAS: prefix + name + fortran_suffix + symbol_suffix (e.g., dgemm_64_) */
#define SLU_BLAS_FUNC(name) SLU_EXPAND4(BLAS_SYMBOL_PREFIX, name, BLAS_FORTRAN_SUFFIX, BLAS_SYMBOL_SUFFIX)
#else
/* MKL and others: prefix + name + symbol_suffix + fortran_suffix (e.g., dgemm_64_) */
#define SLU_BLAS_FUNC(name) SLU_EXPAND4(BLAS_SYMBOL_PREFIX, name, BLAS_SYMBOL_SUFFIX, BLAS_FORTRAN_SUFFIX)
#endif

/*
 * ILP64 BLAS: override all BLAS symbol names with the ILP64-suffixed
 * variants (e.g., dgemm_ -> dgemm_64_ for OpenBLAS, dgemm_64 for MKL).
 * This must come after the Fortran naming blocks above so it overrides them.
 * SLU_BLAS_FUNC() is defined above.
 *
 * Note: LAPACK-like routines (dlacon2, scsum1, dzsum1, icmax1, izmax1) are
 * implemented internally by SuperLU and do NOT need ILP64 renaming.
 */
/* BLAS symbol renaming: active for both ILP64 (suffix) and symbol-prefixed
 * libraries like scipy-openblas (prefix).  When prefix/suffix are empty,
 * SLU_BLAS_FUNC(name) expands back to name_, which is a harmless no-op. */

/* Single precision BLAS */
#undef sswap_
#define sswap_    SLU_BLAS_FUNC(sswap)
#undef saxpy_
#define saxpy_    SLU_BLAS_FUNC(saxpy)
#undef sasum_
#define sasum_    SLU_BLAS_FUNC(sasum)
#undef isamax_
#define isamax_   SLU_BLAS_FUNC(isamax)
#undef scopy_
#define scopy_    SLU_BLAS_FUNC(scopy)
#undef sscal_
#define sscal_    SLU_BLAS_FUNC(sscal)
#undef sger_
#define sger_     SLU_BLAS_FUNC(sger)
#undef snrm2_
#define snrm2_    SLU_BLAS_FUNC(snrm2)
#undef ssymv_
#define ssymv_    SLU_BLAS_FUNC(ssymv)
#undef sdot_
#define sdot_     SLU_BLAS_FUNC(sdot)
#undef ssyr2_
#define ssyr2_    SLU_BLAS_FUNC(ssyr2)
#undef srot_
#define srot_     SLU_BLAS_FUNC(srot)
#undef sgemv_
#define sgemv_    SLU_BLAS_FUNC(sgemv)
#undef strsv_
#define strsv_    SLU_BLAS_FUNC(strsv)
#undef sgemm_
#define sgemm_    SLU_BLAS_FUNC(sgemm)
#undef strsm_
#define strsm_    SLU_BLAS_FUNC(strsm)

/* Double precision BLAS */
#undef dswap_
#define dswap_    SLU_BLAS_FUNC(dswap)
#undef daxpy_
#define daxpy_    SLU_BLAS_FUNC(daxpy)
#undef dasum_
#define dasum_    SLU_BLAS_FUNC(dasum)
#undef idamax_
#define idamax_   SLU_BLAS_FUNC(idamax)
#undef dcopy_
#define dcopy_    SLU_BLAS_FUNC(dcopy)
#undef dscal_
#define dscal_    SLU_BLAS_FUNC(dscal)
#undef dger_
#define dger_     SLU_BLAS_FUNC(dger)
#undef dnrm2_
#define dnrm2_    SLU_BLAS_FUNC(dnrm2)
#undef dsymv_
#define dsymv_    SLU_BLAS_FUNC(dsymv)
#undef ddot_
#define ddot_     SLU_BLAS_FUNC(ddot)
#undef dsyr2_
#define dsyr2_    SLU_BLAS_FUNC(dsyr2)
#undef drot_
#define drot_     SLU_BLAS_FUNC(drot)
#undef dgemv_
#define dgemv_    SLU_BLAS_FUNC(dgemv)
#undef dtrsv_
#define dtrsv_    SLU_BLAS_FUNC(dtrsv)
#undef dgemm_
#define dgemm_    SLU_BLAS_FUNC(dgemm)
#undef dtrsm_
#define dtrsm_    SLU_BLAS_FUNC(dtrsm)

/* Complex single precision BLAS */
#undef cdotc_
#define cdotc_    SLU_BLAS_FUNC(cdotc)
#undef dcabs1_
#define dcabs1_   SLU_BLAS_FUNC(dcabs1)
#undef cswap_
#define cswap_    SLU_BLAS_FUNC(cswap)
#undef caxpy_
#define caxpy_    SLU_BLAS_FUNC(caxpy)
#undef scasum_
#define scasum_   SLU_BLAS_FUNC(scasum)
#undef icamax_
#define icamax_   SLU_BLAS_FUNC(icamax)
#undef ccopy_
#define ccopy_    SLU_BLAS_FUNC(ccopy)
#undef cscal_
#define cscal_    SLU_BLAS_FUNC(cscal)
#undef scnrm2_
#define scnrm2_   SLU_BLAS_FUNC(scnrm2)
#undef cgemv_
#define cgemv_    SLU_BLAS_FUNC(cgemv)
#undef ctrsv_
#define ctrsv_    SLU_BLAS_FUNC(ctrsv)
#undef cgemm_
#define cgemm_    SLU_BLAS_FUNC(cgemm)
#undef ctrsm_
#define ctrsm_    SLU_BLAS_FUNC(ctrsm)
#undef cgerc_
#define cgerc_    SLU_BLAS_FUNC(cgerc)
#undef chemv_
#define chemv_    SLU_BLAS_FUNC(chemv)
#undef cher2_
#define cher2_    SLU_BLAS_FUNC(cher2)

/* Complex double precision BLAS */
#undef zdotc_
#define zdotc_    SLU_BLAS_FUNC(zdotc)
#undef zswap_
#define zswap_    SLU_BLAS_FUNC(zswap)
#undef zaxpy_
#define zaxpy_    SLU_BLAS_FUNC(zaxpy)
#undef dzasum_
#define dzasum_   SLU_BLAS_FUNC(dzasum)
#undef izamax_
#define izamax_   SLU_BLAS_FUNC(izamax)
#undef zcopy_
#define zcopy_    SLU_BLAS_FUNC(zcopy)
#undef zscal_
#define zscal_    SLU_BLAS_FUNC(zscal)
#undef dznrm2_
#define dznrm2_   SLU_BLAS_FUNC(dznrm2)
#undef zgemv_
#define zgemv_    SLU_BLAS_FUNC(zgemv)
#undef ztrsv_
#define ztrsv_    SLU_BLAS_FUNC(ztrsv)
#undef zgemm_
#define zgemm_    SLU_BLAS_FUNC(zgemm)
#undef ztrsm_
#define ztrsm_    SLU_BLAS_FUNC(ztrsm)
#undef zgerc_
#define zgerc_    SLU_BLAS_FUNC(zgerc)
#undef zhemv_
#define zhemv_    SLU_BLAS_FUNC(zhemv)
#undef zher2_
#define zher2_    SLU_BLAS_FUNC(zher2)

#endif /* SCIPY_SLU_BLAS_CONFIG_H */
