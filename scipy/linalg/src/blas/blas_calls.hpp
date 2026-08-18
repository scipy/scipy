/**
 * @file
 * @brief Typed C++ overloads over the (untyped) Fortran BLAS symbols.
 *
 * A wrapper calls `blas::gemv(...)` and overload resolution picks the right symbol from the
 * scalar type -- a wrong argument list is a compile error, where f2py's `char*`-cast function
 * pointer made it a silent ABI bug.  Scalars and single-character flags are taken by value;
 * each overload passes their addresses to Fortran (everything by reference).  No Python.h,
 * no numpy: this layer speaks only the width aliases below, `CBLAS_INT` and `char`.
 *
 * Everything -- including the `extern "C"` prototypes -- lives in `namespace blas`: C linkage
 * governs only the symbol name (unmangled), not the C++ scope, so nothing here leaks into the
 * global namespace.
 *
 * @note Every routine is enumerated twice here (extern prototype + overload) and once more in
 *       the wrapper module.  This layer is mechanical and is intended to be generated from
 *       `cython_blas_signatures.txt` once the approach is validated; hand-written for now.
 */
#pragma once

#include <complex>
#include "scipy_blas_defines.h"
#include "fortran_defs.h"   /* F_FUNC, for the non-BLAS shim symbols of wrap_*_g77_abi.c */
/**
 * MKL ILP64 exports a few LAPACK auxiliaries as `foo_64` (no trailing Fortran underscore)
 * instead of `foo_64_`; of the routines we bind, only `cspr` is affected.  This header
 * `#define`s `cspr_64_ -> cspr_64` (and its siblings), gated on FIX_MKL_2025_ILP64_MISSING_SYMBOL
 * which meson sets only for MKL, and is a no-op for every other library. It must precede the
 * BLAS_FUNC prototypes below. f2py got this transitively via `_blas64_defines.h`.
 */
#include "_mkl_ilp64_fixes.h"

/**
 * Mangling for the g77-ABI complex-dot shims from wrap_*_g77_abi.c. These are not BLAS symbols,
 * so they carry no ILP64 width in their name, but the shim source still suffixes them in ILP64
 * builds. In LP64 both sides use the plain Fortran mangling, so `F_FUNC` is used.
 */
#ifdef HAVE_BLAS_ILP64
#define WRP_FUNC(f, F) BLAS_FUNC(f)
#else
#define WRP_FUNC(f, F) F_FUNC(f, F)
#endif

namespace blas {

/* numpy-style width aliases (float32, ..., complex128); the s/d/c/z flavor columns below
 * stay short and aligned. */
using f32  = float;
using f64  = double;
using c64  = std::complex<float>;
using c128 = std::complex<double>;

extern "C" {
f32 BLAS_FUNC(sasum)(CBLAS_INT *, f32 *,  CBLAS_INT *);
f64 BLAS_FUNC(dasum)(CBLAS_INT *, f64 *,  CBLAS_INT *);
f32 BLAS_FUNC(scasum)(CBLAS_INT *, c64 *,  CBLAS_INT *);
f64 BLAS_FUNC(dzasum)(CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(saxpy)(CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *);
void BLAS_FUNC(daxpy)(CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *);
void BLAS_FUNC(caxpy)(CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(zaxpy)(CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(scopy)(CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *);
void BLAS_FUNC(dcopy)(CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *);
void BLAS_FUNC(ccopy)(CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(zcopy)(CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

f32 BLAS_FUNC(sdot)(CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *);
f64 BLAS_FUNC(ddot)(CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *);

/* The complex dot products return a complex by value -- these are the source of problematic g77 ABI
 * These subroutine-style shims from scipy/_build_utils/src/wrap_*_g77_abi.c (linked via
 * blas_lapack_wrapper_lib) take the result as their first argument instead.  The uppercase
 * spellings (WCDOTCWRP vs CDOTUWRP) mirror the shim source.  Mangled via WRP_FUNC (see above),
 * which is ILP64-aware unlike the real BLAS symbols' BLAS_FUNC. */
void WRP_FUNC(cdotuwrp, CDOTUWRP)(c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void WRP_FUNC(zdotuwrp, ZDOTUWRP)(c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);
void WRP_FUNC(cdotcwrp, WCDOTCWRP)(c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void WRP_FUNC(zdotcwrp, WZDOTCWRP)(c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(sger)(CBLAS_INT *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *);
void BLAS_FUNC(dger)(CBLAS_INT *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *);
void BLAS_FUNC(cgeru)(CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(zgeru)(CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);
void BLAS_FUNC(cgerc)(CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(zgerc)(CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(sgemm)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *);
void BLAS_FUNC(dgemm)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *);
void BLAS_FUNC(cgemm)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zgemm)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);

/* Level-3 symmetric families are standard four-flavor BLAS (csymm/csyrk/csyr2k are not
 * auxiliaries); the hermitian members follow the reference's REAL scalar rules: herk takes
 * real alpha AND beta, her2k real beta only.  The .pyf declared those complex and let the
 * miscast pointer truncate; the wrappers convert as complex and pass .real(). */
void BLAS_FUNC(ssymm)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *);
void BLAS_FUNC(dsymm)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *);
void BLAS_FUNC(csymm)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zsymm)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);
void BLAS_FUNC(chemm)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zhemm)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);

void BLAS_FUNC(ssyrk)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *);
void BLAS_FUNC(dsyrk)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *);
void BLAS_FUNC(csyrk)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zsyrk)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);
void BLAS_FUNC(cherk)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *,  c64 *,  CBLAS_INT *, f32 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zherk)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *,  c128 *, CBLAS_INT *, f64 *,  c128 *, CBLAS_INT *);

void BLAS_FUNC(ssyr2k)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *);
void BLAS_FUNC(dsyr2k)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *);
void BLAS_FUNC(csyr2k)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zsyr2k)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);
void BLAS_FUNC(cher2k)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, f32 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zher2k)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *,  c128 *, CBLAS_INT *);

void BLAS_FUNC(strmm)(char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *);
void BLAS_FUNC(dtrmm)(char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *);
void BLAS_FUNC(ctrmm)(char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(ztrmm)(char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(sgbmv)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *);
void BLAS_FUNC(dgbmv)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *);
void BLAS_FUNC(cgbmv)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zgbmv)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);

void BLAS_FUNC(sgemv)(char *, CBLAS_INT *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *);
void BLAS_FUNC(dgemv)(char *, CBLAS_INT *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *);
void BLAS_FUNC(cgemv)(char *, CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zgemv)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);

f32 BLAS_FUNC(snrm2)(CBLAS_INT *, f32 *,  CBLAS_INT *);
f64 BLAS_FUNC(dnrm2)(CBLAS_INT *, f64 *,  CBLAS_INT *);
f32 BLAS_FUNC(scnrm2)(CBLAS_INT *, c64 *,  CBLAS_INT *);
f64 BLAS_FUNC(dznrm2)(CBLAS_INT *, c128 *, CBLAS_INT *);

/* iamax returns the (1-based, Fortran) index as the integer type of the build's BLAS ABI */
CBLAS_INT BLAS_FUNC(isamax)(CBLAS_INT *, f32 *,  CBLAS_INT *);
CBLAS_INT BLAS_FUNC(idamax)(CBLAS_INT *, f64 *,  CBLAS_INT *);
CBLAS_INT BLAS_FUNC(icamax)(CBLAS_INT *, c64 *,  CBLAS_INT *);
CBLAS_INT BLAS_FUNC(izamax)(CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(srotg)(f32 *,  f32 *,  f32 *, f32 *);
void BLAS_FUNC(drotg)(f64 *,  f64 *,  f64 *, f64 *);
void BLAS_FUNC(crotg)(c64 *,  c64 *,  f32 *, c64 *);
void BLAS_FUNC(zrotg)(c128 *, c128 *, f64 *, c128 *);

/* The packed families are six-flavored in netlib (like symv): cspmv/zspmv and cspr/zspr
 * are LAPACK auxiliary complex-symmetric routines; chp* are the hermitian members.
 * hpr's alpha is REAL; everything else takes the flavor type.  hpr2/spr2 has no
 * complex-symmetric members. */
void BLAS_FUNC(sspmv)(char *, CBLAS_INT *, f32 *,  f32 *,  f32 *,  CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *);
void BLAS_FUNC(dspmv)(char *, CBLAS_INT *, f64 *,  f64 *,  f64 *,  CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *);
void BLAS_FUNC(cspmv)(char *, CBLAS_INT *, c64 *,  c64 *,  c64 *,  CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zspmv)(char *, CBLAS_INT *, c128 *, c128 *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);
void BLAS_FUNC(chpmv)(char *, CBLAS_INT *, c64 *,  c64 *,  c64 *,  CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zhpmv)(char *, CBLAS_INT *, c128 *, c128 *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);

void BLAS_FUNC(sspr)(char *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *);
void BLAS_FUNC(dspr)(char *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *);
void BLAS_FUNC(cspr)(char *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *);
void BLAS_FUNC(zspr)(char *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *);
void BLAS_FUNC(chpr)(char *, CBLAS_INT *, f32 *,  c64 *,  CBLAS_INT *, c64 *);
void BLAS_FUNC(zhpr)(char *, CBLAS_INT *, f64 *,  c128 *, CBLAS_INT *, c128 *);

void BLAS_FUNC(sspr2)(char *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *);
void BLAS_FUNC(dspr2)(char *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *);
void BLAS_FUNC(chpr2)(char *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *);
void BLAS_FUNC(zhpr2)(char *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *);

void BLAS_FUNC(srot)(CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *, f32 *);
void BLAS_FUNC(drot)(CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *, f64 *);
void BLAS_FUNC(csrot)(CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, f32 *, f32 *);
void BLAS_FUNC(zdrot)(CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *);

void BLAS_FUNC(srotm)(CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *);
void BLAS_FUNC(drotm)(CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *);
void BLAS_FUNC(srotmg)(f32 *, f32 *, f32 *, f32 *, f32 *);
void BLAS_FUNC(drotmg)(f64 *, f64 *, f64 *, f64 *, f64 *);

/* scal has six flavors: the four regular ones plus the real-scalar-on-complex-data pair. */
void BLAS_FUNC(sscal)(CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *);
void BLAS_FUNC(dscal)(CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *);
void BLAS_FUNC(cscal)(CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zscal)(CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);
void BLAS_FUNC(csscal)(CBLAS_INT *, f32 *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(zdscal)(CBLAS_INT *, f64 *, c128 *, CBLAS_INT *);

/* The banded symmetric/hermitian family has NO complex-symmetric auxiliaries in LAPACK
 * (checked against netlib): unlike symv/spmv/spr, sbmv is s/d only and the complex
 * members are the hermitian chbmv/zhbmv.  A csbmv_ symbol exported by some OpenBLAS
 * builds is testing-tree leakage -- do not bind it. */
void BLAS_FUNC(ssbmv)(char *, CBLAS_INT *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *);
void BLAS_FUNC(dsbmv)(char *, CBLAS_INT *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *);
void BLAS_FUNC(chbmv)(char *, CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zhbmv)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);

/* symv is a full four-flavor family: csymv/zsymv are LAPACK auxiliary routines (complex
 * *symmetric*, no conjugation), which scipy's .pyf never exposed.  The hermitian chemv/
 * zhemv are the separate overload set -- they collide with csymv/zsymv on the data type
 * alone, hence the distinct names. */
void BLAS_FUNC(ssymv)(char *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *);
void BLAS_FUNC(dsymv)(char *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *);
void BLAS_FUNC(csymv)(char *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zsymv)(char *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);
void BLAS_FUNC(chemv)(char *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(zhemv)(char *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);

/* syr is six-flavored in netlib (csyr/zsyr are LAPACK auxiliary complex-symmetric); her's
 * alpha is real per the reference, and the wrapper passes the real part of the converted
 * complex (mirroring the .pyf's declared-complex/passed-real arrangement).  syr2 has no
 * complex-symmetric members. */
void BLAS_FUNC(ssyr)(char *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *);
void BLAS_FUNC(dsyr)(char *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *);
void BLAS_FUNC(csyr)(char *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(zsyr)(char *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);
void BLAS_FUNC(cher)(char *, CBLAS_INT *, f32 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(zher)(char *, CBLAS_INT *, f64 *,  c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(ssyr2)(char *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *);
void BLAS_FUNC(dsyr2)(char *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *);
void BLAS_FUNC(cher2)(char *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(zher2)(char *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(sswap)(CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *);
void BLAS_FUNC(dswap)(CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *);
void BLAS_FUNC(cswap)(CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(zswap)(CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(stbmv)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *);
void BLAS_FUNC(dtbmv)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *);
void BLAS_FUNC(ctbmv)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(ztbmv)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(stbsv)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *);
void BLAS_FUNC(dtbsv)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *);
void BLAS_FUNC(ctbsv)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(ztbsv)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(stpmv)(char *, char *, char *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *);
void BLAS_FUNC(dtpmv)(char *, char *, char *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *);
void BLAS_FUNC(ctpmv)(char *, char *, char *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(ztpmv)(char *, char *, char *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);

void BLAS_FUNC(stpsv)(char *, char *, char *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *);
void BLAS_FUNC(dtpsv)(char *, char *, char *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *);
void BLAS_FUNC(ctpsv)(char *, char *, char *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *);
void BLAS_FUNC(ztpsv)(char *, char *, char *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);

void BLAS_FUNC(strmv)(char *, char *, char *, CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *);
void BLAS_FUNC(dtrmv)(char *, char *, char *, CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *);
void BLAS_FUNC(ctrmv)(char *, char *, char *, CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(ztrmv)(char *, char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(strsv)(char *, char *, char *, CBLAS_INT *, f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *);
void BLAS_FUNC(dtrsv)(char *, char *, char *, CBLAS_INT *, f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *);
void BLAS_FUNC(ctrsv)(char *, char *, char *, CBLAS_INT *, c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(ztrsv)(char *, char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

void BLAS_FUNC(strsm)(char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *,  f32 *,  CBLAS_INT *, f32 *,  CBLAS_INT *);
void BLAS_FUNC(dtrsm)(char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *,  f64 *,  CBLAS_INT *, f64 *,  CBLAS_INT *);
void BLAS_FUNC(ctrsm)(char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *,  c64 *,  CBLAS_INT *, c64 *,  CBLAS_INT *);
void BLAS_FUNC(ztrsm)(char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);
}

inline f32 asum(CBLAS_INT n, f32 *x,  CBLAS_INT incx) { return BLAS_FUNC(sasum)(&n, x, &incx); }
inline f64 asum(CBLAS_INT n, f64 *x,  CBLAS_INT incx) { return BLAS_FUNC(dasum)(&n, x, &incx); }
inline f32 asum(CBLAS_INT n, c64 *x,  CBLAS_INT incx) { return BLAS_FUNC(scasum)(&n, x, &incx); }
inline f64 asum(CBLAS_INT n, c128 *x, CBLAS_INT incx) { return BLAS_FUNC(dzasum)(&n, x, &incx); }

inline void axpy(CBLAS_INT n, f32 a,  f32 *x,  CBLAS_INT incx, f32 *y,  CBLAS_INT incy) { BLAS_FUNC(saxpy)(&n, &a, x, &incx, y, &incy); }
inline void axpy(CBLAS_INT n, f64 a,  f64 *x,  CBLAS_INT incx, f64 *y,  CBLAS_INT incy) { BLAS_FUNC(daxpy)(&n, &a, x, &incx, y, &incy); }
inline void axpy(CBLAS_INT n, c64 a,  c64 *x,  CBLAS_INT incx, c64 *y,  CBLAS_INT incy) { BLAS_FUNC(caxpy)(&n, &a, x, &incx, y, &incy); }
inline void axpy(CBLAS_INT n, c128 a, c128 *x, CBLAS_INT incx, c128 *y, CBLAS_INT incy) { BLAS_FUNC(zaxpy)(&n, &a, x, &incx, y, &incy); }

inline void copy(CBLAS_INT n, f32 *x,  CBLAS_INT incx, f32 *y,  CBLAS_INT incy) { BLAS_FUNC(scopy)(&n, x, &incx, y, &incy); }
inline void copy(CBLAS_INT n, f64 *x,  CBLAS_INT incx, f64 *y,  CBLAS_INT incy) { BLAS_FUNC(dcopy)(&n, x, &incx, y, &incy); }
inline void copy(CBLAS_INT n, c64 *x,  CBLAS_INT incx, c64 *y,  CBLAS_INT incy) { BLAS_FUNC(ccopy)(&n, x, &incx, y, &incy); }
inline void copy(CBLAS_INT n, c128 *x, CBLAS_INT incx, c128 *y, CBLAS_INT incy) { BLAS_FUNC(zcopy)(&n, x, &incx, y, &incy); }

inline f32 dot(CBLAS_INT n, f32 *x, CBLAS_INT incx, f32 *y, CBLAS_INT incy) { return BLAS_FUNC(sdot)(&n, x, &incx, y, &incy); }
inline f64 dot(CBLAS_INT n, f64 *x, CBLAS_INT incx, f64 *y, CBLAS_INT incy) { return BLAS_FUNC(ddot)(&n, x, &incx, y, &incy); }

inline void ger(CBLAS_INT m, CBLAS_INT n, f32 alpha,  f32 *x,  CBLAS_INT incx, f32 *y,  CBLAS_INT incy, f32 *a,  CBLAS_INT lda) { BLAS_FUNC(sger)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }
inline void ger(CBLAS_INT m, CBLAS_INT n, f64 alpha,  f64 *x,  CBLAS_INT incx, f64 *y,  CBLAS_INT incy, f64 *a,  CBLAS_INT lda) { BLAS_FUNC(dger)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }
inline void geru(CBLAS_INT m, CBLAS_INT n, c64 alpha,  c64 *x,  CBLAS_INT incx, c64 *y,  CBLAS_INT incy, c64 *a,  CBLAS_INT lda) { BLAS_FUNC(cgeru)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }
inline void geru(CBLAS_INT m, CBLAS_INT n, c128 alpha, c128 *x, CBLAS_INT incx, c128 *y, CBLAS_INT incy, c128 *a, CBLAS_INT lda) { BLAS_FUNC(zgeru)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }
inline void gerc(CBLAS_INT m, CBLAS_INT n, c64 alpha,  c64 *x,  CBLAS_INT incx, c64 *y,  CBLAS_INT incy, c64 *a,  CBLAS_INT lda) { BLAS_FUNC(cgerc)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }
inline void gerc(CBLAS_INT m, CBLAS_INT n, c128 alpha, c128 *x, CBLAS_INT incx, c128 *y, CBLAS_INT incy, c128 *a, CBLAS_INT lda) { BLAS_FUNC(zgerc)(&m, &n, &alpha, x, &incx, y, &incy, a, &lda); }

inline CBLAS_INT iamax(CBLAS_INT n, f32 *x,  CBLAS_INT incx) { return BLAS_FUNC(isamax)(&n, x, &incx); }
inline CBLAS_INT iamax(CBLAS_INT n, f64 *x,  CBLAS_INT incx) { return BLAS_FUNC(idamax)(&n, x, &incx); }
inline CBLAS_INT iamax(CBLAS_INT n, c64 *x,  CBLAS_INT incx) { return BLAS_FUNC(icamax)(&n, x, &incx); }
inline CBLAS_INT iamax(CBLAS_INT n, c128 *x, CBLAS_INT incx) { return BLAS_FUNC(izamax)(&n, x, &incx); }

inline c64  dotu(CBLAS_INT n, c64 *x,  CBLAS_INT incx, c64 *y,  CBLAS_INT incy) { c64 r;  WRP_FUNC(cdotuwrp, CDOTUWRP)(&r, &n, x, &incx, y, &incy);  return r; }
inline c128 dotu(CBLAS_INT n, c128 *x, CBLAS_INT incx, c128 *y, CBLAS_INT incy) { c128 r; WRP_FUNC(zdotuwrp, ZDOTUWRP)(&r, &n, x, &incx, y, &incy);  return r; }
inline c64  dotc(CBLAS_INT n, c64 *x,  CBLAS_INT incx, c64 *y,  CBLAS_INT incy) { c64 r;  WRP_FUNC(cdotcwrp, WCDOTCWRP)(&r, &n, x, &incx, y, &incy); return r; }
inline c128 dotc(CBLAS_INT n, c128 *x, CBLAS_INT incx, c128 *y, CBLAS_INT incy) { c128 r; WRP_FUNC(zdotcwrp, WZDOTCWRP)(&r, &n, x, &incx, y, &incy); return r; }

inline void gemm(char ta, char tb, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, f32 alpha,  f32 *a,  CBLAS_INT lda, f32 *b,  CBLAS_INT ldb, f32 beta,  f32 *c,  CBLAS_INT ldc) { BLAS_FUNC(sgemm)(&ta, &tb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void gemm(char ta, char tb, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, f64 alpha,  f64 *a,  CBLAS_INT lda, f64 *b,  CBLAS_INT ldb, f64 beta,  f64 *c,  CBLAS_INT ldc) { BLAS_FUNC(dgemm)(&ta, &tb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void gemm(char ta, char tb, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 *b,  CBLAS_INT ldb, c64 beta,  c64 *c,  CBLAS_INT ldc) { BLAS_FUNC(cgemm)(&ta, &tb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void gemm(char ta, char tb, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, c128 alpha, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, c128 beta, c128 *c, CBLAS_INT ldc) { BLAS_FUNC(zgemm)(&ta, &tb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }

inline void symm(char side, char uplo, CBLAS_INT m, CBLAS_INT n, f32 alpha,  f32 *a,  CBLAS_INT lda, f32 *b,  CBLAS_INT ldb, f32 beta,  f32 *c,  CBLAS_INT ldc) { BLAS_FUNC(ssymm)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void symm(char side, char uplo, CBLAS_INT m, CBLAS_INT n, f64 alpha,  f64 *a,  CBLAS_INT lda, f64 *b,  CBLAS_INT ldb, f64 beta,  f64 *c,  CBLAS_INT ldc) { BLAS_FUNC(dsymm)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void symm(char side, char uplo, CBLAS_INT m, CBLAS_INT n, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 *b,  CBLAS_INT ldb, c64 beta,  c64 *c,  CBLAS_INT ldc) { BLAS_FUNC(csymm)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void symm(char side, char uplo, CBLAS_INT m, CBLAS_INT n, c128 alpha, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, c128 beta, c128 *c, CBLAS_INT ldc) { BLAS_FUNC(zsymm)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void hemm(char side, char uplo, CBLAS_INT m, CBLAS_INT n, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 *b,  CBLAS_INT ldb, c64 beta,  c64 *c,  CBLAS_INT ldc) { BLAS_FUNC(chemm)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void hemm(char side, char uplo, CBLAS_INT m, CBLAS_INT n, c128 alpha, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, c128 beta, c128 *c, CBLAS_INT ldc) { BLAS_FUNC(zhemm)(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }

inline void syrk(char uplo, char trans, CBLAS_INT n, CBLAS_INT k, f32 alpha,  f32 *a,  CBLAS_INT lda, f32 beta,  f32 *c,  CBLAS_INT ldc) { BLAS_FUNC(ssyrk)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc); }
inline void syrk(char uplo, char trans, CBLAS_INT n, CBLAS_INT k, f64 alpha,  f64 *a,  CBLAS_INT lda, f64 beta,  f64 *c,  CBLAS_INT ldc) { BLAS_FUNC(dsyrk)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc); }
inline void syrk(char uplo, char trans, CBLAS_INT n, CBLAS_INT k, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 beta,  c64 *c,  CBLAS_INT ldc) { BLAS_FUNC(csyrk)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc); }
inline void syrk(char uplo, char trans, CBLAS_INT n, CBLAS_INT k, c128 alpha, c128 *a, CBLAS_INT lda, c128 beta, c128 *c, CBLAS_INT ldc) { BLAS_FUNC(zsyrk)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc); }
inline void herk(char uplo, char trans, CBLAS_INT n, CBLAS_INT k, f32 alpha,  c64 *a,  CBLAS_INT lda, f32 beta,  c64 *c,  CBLAS_INT ldc) { BLAS_FUNC(cherk)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc); }
inline void herk(char uplo, char trans, CBLAS_INT n, CBLAS_INT k, f64 alpha,  c128 *a, CBLAS_INT lda, f64 beta,  c128 *c, CBLAS_INT ldc) { BLAS_FUNC(zherk)(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc); }

inline void syr2k(char uplo, char trans, CBLAS_INT n, CBLAS_INT k, f32 alpha,  f32 *a,  CBLAS_INT lda, f32 *b,  CBLAS_INT ldb, f32 beta,  f32 *c,  CBLAS_INT ldc) { BLAS_FUNC(ssyr2k)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void syr2k(char uplo, char trans, CBLAS_INT n, CBLAS_INT k, f64 alpha,  f64 *a,  CBLAS_INT lda, f64 *b,  CBLAS_INT ldb, f64 beta,  f64 *c,  CBLAS_INT ldc) { BLAS_FUNC(dsyr2k)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void syr2k(char uplo, char trans, CBLAS_INT n, CBLAS_INT k, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 *b,  CBLAS_INT ldb, c64 beta,  c64 *c,  CBLAS_INT ldc) { BLAS_FUNC(csyr2k)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void syr2k(char uplo, char trans, CBLAS_INT n, CBLAS_INT k, c128 alpha, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, c128 beta, c128 *c, CBLAS_INT ldc) { BLAS_FUNC(zsyr2k)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void her2k(char uplo, char trans, CBLAS_INT n, CBLAS_INT k, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 *b,  CBLAS_INT ldb, f32 beta,  c64 *c,  CBLAS_INT ldc) { BLAS_FUNC(cher2k)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }
inline void her2k(char uplo, char trans, CBLAS_INT n, CBLAS_INT k, c128 alpha, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, f64 beta,  c128 *c, CBLAS_INT ldc) { BLAS_FUNC(zher2k)(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); }

inline void trmm(char side, char uplo, char transa, char diag, CBLAS_INT m, CBLAS_INT n, f32 alpha,  f32 *a,  CBLAS_INT lda, f32 *b,  CBLAS_INT ldb) { BLAS_FUNC(strmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }
inline void trmm(char side, char uplo, char transa, char diag, CBLAS_INT m, CBLAS_INT n, f64 alpha,  f64 *a,  CBLAS_INT lda, f64 *b,  CBLAS_INT ldb) { BLAS_FUNC(dtrmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }
inline void trmm(char side, char uplo, char transa, char diag, CBLAS_INT m, CBLAS_INT n, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 *b,  CBLAS_INT ldb) { BLAS_FUNC(ctrmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }
inline void trmm(char side, char uplo, char transa, char diag, CBLAS_INT m, CBLAS_INT n, c128 alpha, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb) { BLAS_FUNC(ztrmm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }

inline void gbmv(char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, f32 alpha,  f32 *a,  CBLAS_INT lda, f32 *x,  CBLAS_INT incx, f32 beta,  f32 *y,  CBLAS_INT incy) { BLAS_FUNC(sgbmv)(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void gbmv(char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, f64 alpha,  f64 *a,  CBLAS_INT lda, f64 *x,  CBLAS_INT incx, f64 beta,  f64 *y,  CBLAS_INT incy) { BLAS_FUNC(dgbmv)(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void gbmv(char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 *x,  CBLAS_INT incx, c64 beta,  c64 *y,  CBLAS_INT incy) { BLAS_FUNC(cgbmv)(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void gbmv(char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, c128 alpha, c128 *a, CBLAS_INT lda, c128 *x, CBLAS_INT incx, c128 beta, c128 *y, CBLAS_INT incy) { BLAS_FUNC(zgbmv)(&trans, &m, &n, &kl, &ku, &alpha, a, &lda, x, &incx, &beta, y, &incy); }

inline void gemv(char trans, CBLAS_INT m, CBLAS_INT n, f32 alpha,  f32 *a,  CBLAS_INT lda, f32 *x,  CBLAS_INT incx, f32 beta,  f32 *y,  CBLAS_INT incy) { BLAS_FUNC(sgemv)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void gemv(char trans, CBLAS_INT m, CBLAS_INT n, f64 alpha,  f64 *a,  CBLAS_INT lda, f64 *x,  CBLAS_INT incx, f64 beta,  f64 *y,  CBLAS_INT incy) { BLAS_FUNC(dgemv)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void gemv(char trans, CBLAS_INT m, CBLAS_INT n, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 *x,  CBLAS_INT incx, c64 beta,  c64 *y,  CBLAS_INT incy) { BLAS_FUNC(cgemv)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void gemv(char trans, CBLAS_INT m, CBLAS_INT n, c128 alpha, c128 *a, CBLAS_INT lda, c128 *x, CBLAS_INT incx, c128 beta, c128 *y, CBLAS_INT incy) { BLAS_FUNC(zgemv)(&trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }

inline f32 nrm2(CBLAS_INT n, f32 *x,  CBLAS_INT incx) { return BLAS_FUNC(snrm2)(&n, x, &incx); }
inline f64 nrm2(CBLAS_INT n, f64 *x,  CBLAS_INT incx) { return BLAS_FUNC(dnrm2)(&n, x, &incx); }
inline f32 nrm2(CBLAS_INT n, c64 *x,  CBLAS_INT incx) { return BLAS_FUNC(scnrm2)(&n, x, &incx); }
inline f64 nrm2(CBLAS_INT n, c128 *x, CBLAS_INT incx) { return BLAS_FUNC(dznrm2)(&n, x, &incx); }

/* rot/rotg take their in-out scalars by reference; Fortran overwrites them in place. */
inline void spmv(char uplo, CBLAS_INT n, f32 alpha,  f32 *ap,  f32 *x,  CBLAS_INT incx, f32 beta,  f32 *y,  CBLAS_INT incy) { BLAS_FUNC(sspmv)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy); }
inline void spmv(char uplo, CBLAS_INT n, f64 alpha,  f64 *ap,  f64 *x,  CBLAS_INT incx, f64 beta,  f64 *y,  CBLAS_INT incy) { BLAS_FUNC(dspmv)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy); }
inline void spmv(char uplo, CBLAS_INT n, c64 alpha,  c64 *ap,  c64 *x,  CBLAS_INT incx, c64 beta,  c64 *y,  CBLAS_INT incy) { BLAS_FUNC(cspmv)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy); }
inline void spmv(char uplo, CBLAS_INT n, c128 alpha, c128 *ap, c128 *x, CBLAS_INT incx, c128 beta, c128 *y, CBLAS_INT incy) { BLAS_FUNC(zspmv)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy); }
inline void hpmv(char uplo, CBLAS_INT n, c64 alpha,  c64 *ap,  c64 *x,  CBLAS_INT incx, c64 beta,  c64 *y,  CBLAS_INT incy) { BLAS_FUNC(chpmv)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy); }
inline void hpmv(char uplo, CBLAS_INT n, c128 alpha, c128 *ap, c128 *x, CBLAS_INT incx, c128 beta, c128 *y, CBLAS_INT incy) { BLAS_FUNC(zhpmv)(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy); }

inline void spr(char uplo, CBLAS_INT n, f32 alpha,  f32 *x,  CBLAS_INT incx, f32 *ap)  { BLAS_FUNC(sspr)(&uplo, &n, &alpha, x, &incx, ap); }
inline void spr(char uplo, CBLAS_INT n, f64 alpha,  f64 *x,  CBLAS_INT incx, f64 *ap)  { BLAS_FUNC(dspr)(&uplo, &n, &alpha, x, &incx, ap); }
inline void spr(char uplo, CBLAS_INT n, c64 alpha,  c64 *x,  CBLAS_INT incx, c64 *ap)  { BLAS_FUNC(cspr)(&uplo, &n, &alpha, x, &incx, ap); }
inline void spr(char uplo, CBLAS_INT n, c128 alpha, c128 *x, CBLAS_INT incx, c128 *ap) { BLAS_FUNC(zspr)(&uplo, &n, &alpha, x, &incx, ap); }
inline void hpr(char uplo, CBLAS_INT n, f32 alpha,  c64 *x,  CBLAS_INT incx, c64 *ap)  { BLAS_FUNC(chpr)(&uplo, &n, &alpha, x, &incx, ap); }
inline void hpr(char uplo, CBLAS_INT n, f64 alpha,  c128 *x, CBLAS_INT incx, c128 *ap) { BLAS_FUNC(zhpr)(&uplo, &n, &alpha, x, &incx, ap); }

inline void spr2(char uplo, CBLAS_INT n, f32 alpha,  f32 *x,  CBLAS_INT incx, f32 *y,  CBLAS_INT incy, f32 *ap)  { BLAS_FUNC(sspr2)(&uplo, &n, &alpha, x, &incx, y, &incy, ap); }
inline void spr2(char uplo, CBLAS_INT n, f64 alpha,  f64 *x,  CBLAS_INT incx, f64 *y,  CBLAS_INT incy, f64 *ap)  { BLAS_FUNC(dspr2)(&uplo, &n, &alpha, x, &incx, y, &incy, ap); }
inline void hpr2(char uplo, CBLAS_INT n, c64 alpha,  c64 *x,  CBLAS_INT incx, c64 *y,  CBLAS_INT incy, c64 *ap)  { BLAS_FUNC(chpr2)(&uplo, &n, &alpha, x, &incx, y, &incy, ap); }
inline void hpr2(char uplo, CBLAS_INT n, c128 alpha, c128 *x, CBLAS_INT incx, c128 *y, CBLAS_INT incy, c128 *ap) { BLAS_FUNC(zhpr2)(&uplo, &n, &alpha, x, &incx, y, &incy, ap); }

inline void rot(CBLAS_INT n, f32 *x,  CBLAS_INT incx, f32 *y,  CBLAS_INT incy, f32 c, f32 s) { BLAS_FUNC(srot)(&n, x, &incx, y, &incy, &c, &s); }
inline void rot(CBLAS_INT n, f64 *x,  CBLAS_INT incx, f64 *y,  CBLAS_INT incy, f64 c, f64 s) { BLAS_FUNC(drot)(&n, x, &incx, y, &incy, &c, &s); }
inline void rot(CBLAS_INT n, c64 *x,  CBLAS_INT incx, c64 *y,  CBLAS_INT incy, f32 c, f32 s) { BLAS_FUNC(csrot)(&n, x, &incx, y, &incy, &c, &s); }
inline void rot(CBLAS_INT n, c128 *x, CBLAS_INT incx, c128 *y, CBLAS_INT incy, f64 c, f64 s) { BLAS_FUNC(zdrot)(&n, x, &incx, y, &incy, &c, &s); }

inline void rotg(f32 &a,  f32 &b,  f32 &c, f32 &s)  { BLAS_FUNC(srotg)(&a, &b, &c, &s); }
inline void rotg(f64 &a,  f64 &b,  f64 &c, f64 &s)  { BLAS_FUNC(drotg)(&a, &b, &c, &s); }
inline void rotg(c64 &a,  c64 &b,  f32 &c, c64 &s)  { BLAS_FUNC(crotg)(&a, &b, &c, &s); }
inline void rotg(c128 &a, c128 &b, f64 &c, c128 &s) { BLAS_FUNC(zrotg)(&a, &b, &c, &s); }

inline void rotm(CBLAS_INT n, f32 *x, CBLAS_INT incx, f32 *y, CBLAS_INT incy, f32 *param) { BLAS_FUNC(srotm)(&n, x, &incx, y, &incy, param); }
inline void rotm(CBLAS_INT n, f64 *x, CBLAS_INT incx, f64 *y, CBLAS_INT incy, f64 *param) { BLAS_FUNC(drotm)(&n, x, &incx, y, &incy, param); }

inline void rotmg(f32 &d1, f32 &d2, f32 &x1, f32 y1, f32 *param) { BLAS_FUNC(srotmg)(&d1, &d2, &x1, &y1, param); }
inline void rotmg(f64 &d1, f64 &d2, f64 &x1, f64 y1, f64 *param) { BLAS_FUNC(drotmg)(&d1, &d2, &x1, &y1, param); }

/* the (scalar type, data type) pair selects among all six scal flavors */
inline void scal(CBLAS_INT n, f32 a,  f32 *x,  CBLAS_INT incx) { BLAS_FUNC(sscal)(&n, &a, x, &incx); }
inline void scal(CBLAS_INT n, f64 a,  f64 *x,  CBLAS_INT incx) { BLAS_FUNC(dscal)(&n, &a, x, &incx); }
inline void scal(CBLAS_INT n, c64 a,  c64 *x,  CBLAS_INT incx) { BLAS_FUNC(cscal)(&n, &a, x, &incx); }
inline void scal(CBLAS_INT n, c128 a, c128 *x, CBLAS_INT incx) { BLAS_FUNC(zscal)(&n, &a, x, &incx); }
inline void scal(CBLAS_INT n, f32 a,  c64 *x,  CBLAS_INT incx) { BLAS_FUNC(csscal)(&n, &a, x, &incx); }
inline void scal(CBLAS_INT n, f64 a,  c128 *x, CBLAS_INT incx) { BLAS_FUNC(zdscal)(&n, &a, x, &incx); }

inline void sbmv(char uplo, CBLAS_INT n, CBLAS_INT k, f32 alpha,  f32 *a,  CBLAS_INT lda, f32 *x,  CBLAS_INT incx, f32 beta,  f32 *y,  CBLAS_INT incy) { BLAS_FUNC(ssbmv)(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void sbmv(char uplo, CBLAS_INT n, CBLAS_INT k, f64 alpha,  f64 *a,  CBLAS_INT lda, f64 *x,  CBLAS_INT incx, f64 beta,  f64 *y,  CBLAS_INT incy) { BLAS_FUNC(dsbmv)(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void hbmv(char uplo, CBLAS_INT n, CBLAS_INT k, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 *x,  CBLAS_INT incx, c64 beta,  c64 *y,  CBLAS_INT incy) { BLAS_FUNC(chbmv)(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void hbmv(char uplo, CBLAS_INT n, CBLAS_INT k, c128 alpha, c128 *a, CBLAS_INT lda, c128 *x, CBLAS_INT incx, c128 beta, c128 *y, CBLAS_INT incy) { BLAS_FUNC(zhbmv)(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy); }

inline void symv(char uplo, CBLAS_INT n, f32 alpha,  f32 *a,  CBLAS_INT lda, f32 *x,  CBLAS_INT incx, f32 beta,  f32 *y,  CBLAS_INT incy) { BLAS_FUNC(ssymv)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void symv(char uplo, CBLAS_INT n, f64 alpha,  f64 *a,  CBLAS_INT lda, f64 *x,  CBLAS_INT incx, f64 beta,  f64 *y,  CBLAS_INT incy) { BLAS_FUNC(dsymv)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void symv(char uplo, CBLAS_INT n, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 *x,  CBLAS_INT incx, c64 beta,  c64 *y,  CBLAS_INT incy) { BLAS_FUNC(csymv)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void symv(char uplo, CBLAS_INT n, c128 alpha, c128 *a, CBLAS_INT lda, c128 *x, CBLAS_INT incx, c128 beta, c128 *y, CBLAS_INT incy) { BLAS_FUNC(zsymv)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void hemv(char uplo, CBLAS_INT n, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 *x,  CBLAS_INT incx, c64 beta,  c64 *y,  CBLAS_INT incy) { BLAS_FUNC(chemv)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }
inline void hemv(char uplo, CBLAS_INT n, c128 alpha, c128 *a, CBLAS_INT lda, c128 *x, CBLAS_INT incx, c128 beta, c128 *y, CBLAS_INT incy) { BLAS_FUNC(zhemv)(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy); }

inline void syr(char uplo, CBLAS_INT n, f32 alpha,  f32 *x,  CBLAS_INT incx, f32 *a,  CBLAS_INT lda) { BLAS_FUNC(ssyr)(&uplo, &n, &alpha, x, &incx, a, &lda); }
inline void syr(char uplo, CBLAS_INT n, f64 alpha,  f64 *x,  CBLAS_INT incx, f64 *a,  CBLAS_INT lda) { BLAS_FUNC(dsyr)(&uplo, &n, &alpha, x, &incx, a, &lda); }
inline void syr(char uplo, CBLAS_INT n, c64 alpha,  c64 *x,  CBLAS_INT incx, c64 *a,  CBLAS_INT lda) { BLAS_FUNC(csyr)(&uplo, &n, &alpha, x, &incx, a, &lda); }
inline void syr(char uplo, CBLAS_INT n, c128 alpha, c128 *x, CBLAS_INT incx, c128 *a, CBLAS_INT lda) { BLAS_FUNC(zsyr)(&uplo, &n, &alpha, x, &incx, a, &lda); }
inline void her(char uplo, CBLAS_INT n, f32 alpha,  c64 *x,  CBLAS_INT incx, c64 *a,  CBLAS_INT lda) { BLAS_FUNC(cher)(&uplo, &n, &alpha, x, &incx, a, &lda); }
inline void her(char uplo, CBLAS_INT n, f64 alpha,  c128 *x, CBLAS_INT incx, c128 *a, CBLAS_INT lda) { BLAS_FUNC(zher)(&uplo, &n, &alpha, x, &incx, a, &lda); }

inline void syr2(char uplo, CBLAS_INT n, f32 alpha,  f32 *x,  CBLAS_INT incx, f32 *y,  CBLAS_INT incy, f32 *a,  CBLAS_INT lda) { BLAS_FUNC(ssyr2)(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda); }
inline void syr2(char uplo, CBLAS_INT n, f64 alpha,  f64 *x,  CBLAS_INT incx, f64 *y,  CBLAS_INT incy, f64 *a,  CBLAS_INT lda) { BLAS_FUNC(dsyr2)(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda); }
inline void her2(char uplo, CBLAS_INT n, c64 alpha,  c64 *x,  CBLAS_INT incx, c64 *y,  CBLAS_INT incy, c64 *a,  CBLAS_INT lda) { BLAS_FUNC(cher2)(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda); }
inline void her2(char uplo, CBLAS_INT n, c128 alpha, c128 *x, CBLAS_INT incx, c128 *y, CBLAS_INT incy, c128 *a, CBLAS_INT lda) { BLAS_FUNC(zher2)(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda); }

inline void swap(CBLAS_INT n, f32 *x,  CBLAS_INT incx, f32 *y,  CBLAS_INT incy) { BLAS_FUNC(sswap)(&n, x, &incx, y, &incy); }
inline void swap(CBLAS_INT n, f64 *x,  CBLAS_INT incx, f64 *y,  CBLAS_INT incy) { BLAS_FUNC(dswap)(&n, x, &incx, y, &incy); }
inline void swap(CBLAS_INT n, c64 *x,  CBLAS_INT incx, c64 *y,  CBLAS_INT incy) { BLAS_FUNC(cswap)(&n, x, &incx, y, &incy); }
inline void swap(CBLAS_INT n, c128 *x, CBLAS_INT incx, c128 *y, CBLAS_INT incy) { BLAS_FUNC(zswap)(&n, x, &incx, y, &incy); }

inline void tbmv(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT k, f32 *a,  CBLAS_INT lda, f32 *x,  CBLAS_INT incx) { BLAS_FUNC(stbmv)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx); }
inline void tbmv(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT k, f64 *a,  CBLAS_INT lda, f64 *x,  CBLAS_INT incx) { BLAS_FUNC(dtbmv)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx); }
inline void tbmv(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT k, c64 *a,  CBLAS_INT lda, c64 *x,  CBLAS_INT incx) { BLAS_FUNC(ctbmv)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx); }
inline void tbmv(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT k, c128 *a, CBLAS_INT lda, c128 *x, CBLAS_INT incx) { BLAS_FUNC(ztbmv)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx); }

inline void tbsv(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT k, f32 *a,  CBLAS_INT lda, f32 *x,  CBLAS_INT incx) { BLAS_FUNC(stbsv)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx); }
inline void tbsv(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT k, f64 *a,  CBLAS_INT lda, f64 *x,  CBLAS_INT incx) { BLAS_FUNC(dtbsv)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx); }
inline void tbsv(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT k, c64 *a,  CBLAS_INT lda, c64 *x,  CBLAS_INT incx) { BLAS_FUNC(ctbsv)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx); }
inline void tbsv(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT k, c128 *a, CBLAS_INT lda, c128 *x, CBLAS_INT incx) { BLAS_FUNC(ztbsv)(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx); }

inline void tpmv(char uplo, char trans, char diag, CBLAS_INT n, f32 *ap,  f32 *x,  CBLAS_INT incx) { BLAS_FUNC(stpmv)(&uplo, &trans, &diag, &n, ap, x, &incx); }
inline void tpmv(char uplo, char trans, char diag, CBLAS_INT n, f64 *ap,  f64 *x,  CBLAS_INT incx) { BLAS_FUNC(dtpmv)(&uplo, &trans, &diag, &n, ap, x, &incx); }
inline void tpmv(char uplo, char trans, char diag, CBLAS_INT n, c64 *ap,  c64 *x,  CBLAS_INT incx) { BLAS_FUNC(ctpmv)(&uplo, &trans, &diag, &n, ap, x, &incx); }
inline void tpmv(char uplo, char trans, char diag, CBLAS_INT n, c128 *ap, c128 *x, CBLAS_INT incx) { BLAS_FUNC(ztpmv)(&uplo, &trans, &diag, &n, ap, x, &incx); }

inline void tpsv(char uplo, char trans, char diag, CBLAS_INT n, f32 *ap,  f32 *x,  CBLAS_INT incx) { BLAS_FUNC(stpsv)(&uplo, &trans, &diag, &n, ap, x, &incx); }
inline void tpsv(char uplo, char trans, char diag, CBLAS_INT n, f64 *ap,  f64 *x,  CBLAS_INT incx) { BLAS_FUNC(dtpsv)(&uplo, &trans, &diag, &n, ap, x, &incx); }
inline void tpsv(char uplo, char trans, char diag, CBLAS_INT n, c64 *ap,  c64 *x,  CBLAS_INT incx) { BLAS_FUNC(ctpsv)(&uplo, &trans, &diag, &n, ap, x, &incx); }
inline void tpsv(char uplo, char trans, char diag, CBLAS_INT n, c128 *ap, c128 *x, CBLAS_INT incx) { BLAS_FUNC(ztpsv)(&uplo, &trans, &diag, &n, ap, x, &incx); }

inline void trmv(char uplo, char trans, char diag, CBLAS_INT n, f32 *a,  CBLAS_INT lda, f32 *x,  CBLAS_INT incx) { BLAS_FUNC(strmv)(&uplo, &trans, &diag, &n, a, &lda, x, &incx); }
inline void trmv(char uplo, char trans, char diag, CBLAS_INT n, f64 *a,  CBLAS_INT lda, f64 *x,  CBLAS_INT incx) { BLAS_FUNC(dtrmv)(&uplo, &trans, &diag, &n, a, &lda, x, &incx); }
inline void trmv(char uplo, char trans, char diag, CBLAS_INT n, c64 *a,  CBLAS_INT lda, c64 *x,  CBLAS_INT incx) { BLAS_FUNC(ctrmv)(&uplo, &trans, &diag, &n, a, &lda, x, &incx); }
inline void trmv(char uplo, char trans, char diag, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *x, CBLAS_INT incx) { BLAS_FUNC(ztrmv)(&uplo, &trans, &diag, &n, a, &lda, x, &incx); }

inline void trsv(char uplo, char trans, char diag, CBLAS_INT n, f32 *a,  CBLAS_INT lda, f32 *x,  CBLAS_INT incx) { BLAS_FUNC(strsv)(&uplo, &trans, &diag, &n, a, &lda, x, &incx); }
inline void trsv(char uplo, char trans, char diag, CBLAS_INT n, f64 *a,  CBLAS_INT lda, f64 *x,  CBLAS_INT incx) { BLAS_FUNC(dtrsv)(&uplo, &trans, &diag, &n, a, &lda, x, &incx); }
inline void trsv(char uplo, char trans, char diag, CBLAS_INT n, c64 *a,  CBLAS_INT lda, c64 *x,  CBLAS_INT incx) { BLAS_FUNC(ctrsv)(&uplo, &trans, &diag, &n, a, &lda, x, &incx); }
inline void trsv(char uplo, char trans, char diag, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *x, CBLAS_INT incx) { BLAS_FUNC(ztrsv)(&uplo, &trans, &diag, &n, a, &lda, x, &incx); }

inline void trsm(char side, char uplo, char transa, char diag, CBLAS_INT m, CBLAS_INT n, f32 alpha,  f32 *a,  CBLAS_INT lda, f32 *b,  CBLAS_INT ldb) { BLAS_FUNC(strsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }
inline void trsm(char side, char uplo, char transa, char diag, CBLAS_INT m, CBLAS_INT n, f64 alpha,  f64 *a,  CBLAS_INT lda, f64 *b,  CBLAS_INT ldb) { BLAS_FUNC(dtrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }
inline void trsm(char side, char uplo, char transa, char diag, CBLAS_INT m, CBLAS_INT n, c64 alpha,  c64 *a,  CBLAS_INT lda, c64 *b,  CBLAS_INT ldb) { BLAS_FUNC(ctrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }
inline void trsm(char side, char uplo, char transa, char diag, CBLAS_INT m, CBLAS_INT n, c128 alpha, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb) { BLAS_FUNC(ztrsm)(&side, &uplo, &transa, &diag, &m, &n, &alpha, a, &lda, b, &ldb); }

}  // namespace blas
