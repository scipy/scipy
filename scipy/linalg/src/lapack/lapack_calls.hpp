/**
 * @file
 * @brief Fortran declarations and typed C++ overloads for the callback-bearing LAPACK routines.
 *
 * The overloads take scalars by value and pass their addresses on, so wrapper bodies never spell
 * `&`. More families are added here alongside their wrapper .cpp files.
 *
 * The declarations use only the shared width aliases, `CBLAS_INT` and `char`, so this header
 * includes neither Python.h nor numpy. A routine's callback is therefore declared elsewhere:
 * the selector trampolines and the per-flavor traits that name their Python keyword arguments
 * are in `lapack_callback.hpp`.
 */
#pragma once

#include "wrapper_types.hpp"       /* f32, f64, c64, c128 types and helper templates */
#include "scipy_blas_defines.h"    /* CBLAS_INT and the BLAS_FUNC symbol mangling    */

/* MKL ILP64 exports a handful of LAPACK auxiliaries as `foo_64` rather than `foo_64_`; this
 * header patches those spellings and is a no-op everywhere else.  It must precede the
 * BLAS_FUNC prototypes below. */
#include "_mkl_ilp64_fixes.h"

namespace lapack {

    /* The shared width aliases, so the s/d/c/z flavor columns below stay short and aligned.
     * `lapack_helpers.hpp` re-exports the same names, which is redundant but legal and lets
     * either header be included on its own. */
    using wrapper::f32;
    using wrapper::f64;
    using wrapper::c64;
    using wrapper::c128;

    extern "C" {
        void BLAS_FUNC(sgees)(char *, char *, CBLAS_INT (*)(f32 *, f32 *), CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgees)(char *, char *, CBLAS_INT (*)(f64 *, f64 *), CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgees)(char *, char *, CBLAS_INT (*)(c64 *), CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgees)(char *, char *, CBLAS_INT (*)(c128 *), CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgges)(char *, char *, char *, CBLAS_INT (*)(f32 *, f32 *, f32 *), CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgges)(char *, char *, char *, CBLAS_INT (*)(f64 *, f64 *, f64 *), CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgges)(char *, char *, char *, CBLAS_INT (*)(c64 *, c64 *), CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, c64 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgges)(char *, char *, char *, CBLAS_INT (*)(c128 *, c128 *), CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgebal)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dgebal)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cgebal)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zgebal)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sgecon)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgecon)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgecon)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, c64 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zgecon)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, c128 *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sgesv)(CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgesv)(CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgesv)(CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgesv)(CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgetrf)(CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgetrf)(CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgetrf)(CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgetrf)(CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgetrs)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgetrs)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgetrs)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgetrs)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgetc2)(CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgetc2)(CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgetc2)(CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgetc2)(CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgesc2)(CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *);
        void BLAS_FUNC(dgesc2)(CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *);
        void BLAS_FUNC(cgesc2)(CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, f32 *);
        void BLAS_FUNC(zgesc2)(CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, f64 *);

        void BLAS_FUNC(sgetri)(CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgetri)(CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgetri)(CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgetri)(CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgesdd)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgesdd)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgesdd)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgesdd)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgesvd)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgesvd)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgesvd)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zgesvd)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sgels)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgels)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgels)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgels)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgeqp3)(CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgeqp3)(CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgeqp3)(CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zgeqp3)(CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sgelsy)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgelsy)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgelsy)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zgelsy)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sgelsd)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgelsd)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgelsd)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgelsd)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgelss)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgelss)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgelss)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zgelss)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sgehrd)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgehrd)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgehrd)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgehrd)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgeqrf)(CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgeqrf)(CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgeqrf)(CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgeqrf)(CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgeqrfp)(CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgeqrfp)(CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgeqrfp)(CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgeqrfp)(CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgerqf)(CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgerqf)(CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgerqf)(CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgerqf)(CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgeequ)(CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dgeequ)(CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cgeequ)(CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zgeequ)(CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sgeequb)(CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dgeequb)(CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cgeequb)(CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zgeequb)(CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sgeev)(char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgeev)(char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgeev)(char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zgeev)(char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sggev)(char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dggev)(char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cggev)(char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zggev)(char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sgesvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, char *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgesvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, char *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgesvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, char *, f32 *, f32 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, f32 *, c64 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zgesvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, char *, f64 *, f64 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, f64 *, c128 *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sgtsv)(CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgtsv)(CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgtsv)(CBLAS_INT *, CBLAS_INT *, c64 *, c64 *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgtsv)(CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgttrf)(CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgttrf)(CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgttrf)(CBLAS_INT *, c64 *, c64 *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgttrf)(CBLAS_INT *, c128 *, c128 *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgttrs)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgttrs)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgttrs)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, c64 *, c64 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgttrs)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        /* The complex flavors take no IWORK: their condition estimator works entirely in the
         * complex WORK buffer.  Same split in `?gtsvx` below, where complex takes RWORK instead. */
        void BLAS_FUNC(sgtcon)(char *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgtcon)(char *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgtcon)(char *, CBLAS_INT *, c64 *, c64 *, c64 *, c64 *, CBLAS_INT *, f32 *, f32 *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zgtcon)(char *, CBLAS_INT *, c128 *, c128 *, c128 *, c128 *, CBLAS_INT *, f64 *, f64 *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(sgtsvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgtsvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgtsvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, c64 *, c64 *, c64 *, c64 *, c64 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, f32 *, c64 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zgtsvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, c128 *, c128 *, c128 *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, f64 *, c128 *, f64 *, CBLAS_INT *);

        /* The symmetric-tridiagonal eigensolvers below are real-only: a Hermitian tridiagonal
         * matrix is unitarily similar to a real symmetric one, so LAPACK ships no c/z flavors. */
        void BLAS_FUNC(sstev)(char *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dstev)(char *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sstebz)(char *, char *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dstebz)(char *, char *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(ssterf)(CBLAS_INT *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dsterf)(CBLAS_INT *, f64 *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sstein)(CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dstein)(CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sstemr)(char *, char *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dstemr)(char *, char *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sstevd)(char *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dstevd)(char *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgbsv)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgbsv)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgbsv)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgbsv)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgbtrf)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgbtrf)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgbtrf)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgbtrf)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgbtrs)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgbtrs)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgbtrs)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgbtrs)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        /* The same real/complex split as `?gtcon`: IWORK for the real flavors, RWORK for the complex. */
        void BLAS_FUNC(sgbcon)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgbcon)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgbcon)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, c64 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zgbcon)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, c128 *, f64 *, CBLAS_INT *);

        /* A value-returning Fortran function, like BLAS's `?nrm2`: the norm comes back as the
         * real counterpart of the flavor, so `clangb` returns f32 and `zlangb` returns f64. */
        f32 BLAS_FUNC(slangb)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *);
        f64 BLAS_FUNC(dlangb)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *);
        f32 BLAS_FUNC(clangb)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *);
        f64 BLAS_FUNC(zlangb)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *);

        /* `tol` and `work` are real beside a complex `a`: a pivot threshold and the diagonal
         * magnitudes the pivoting compares, both of which are magnitudes. */
        void BLAS_FUNC(spstrf)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dpstrf)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cpstrf)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zpstrf)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(spstf2)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dpstf2)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cpstf2)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zpstf2)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sposv)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dposv)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cposv)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zposv)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        /* The last buffer is IWORK for the real flavors and RWORK for the complex ones -- same
         * position, same length `n`, different element type; see the `W` alias in the wrappers. */
        void BLAS_FUNC(sposvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, char *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dposvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, char *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cposvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, char *, f32 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, f32 *, c64 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zposvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, char *, f64 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, f64 *, c128 *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(spocon)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dpocon)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cpocon)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, c64 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zpocon)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, c128 *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(spotrf)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dpotrf)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cpotrf)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zpotrf)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(spotrs)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dpotrs)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cpotrs)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zpotrs)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(spotri)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dpotri)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cpotri)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zpotri)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        /* Throughout this group `d` is real even beside a complex `e`: it is the diagonal of
         * D in the L*D*L**H factorization, which is real for a Hermitian positive definite
         * matrix.  `e` follows the flavor. */
        void BLAS_FUNC(sptsv)(CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dptsv)(CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cptsv)(CBLAS_INT *, CBLAS_INT *, f32 *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zptsv)(CBLAS_INT *, CBLAS_INT *, f64 *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(spttrf)(CBLAS_INT *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dpttrf)(CBLAS_INT *, f64 *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cpttrf)(CBLAS_INT *, f32 *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zpttrf)(CBLAS_INT *, f64 *, c128 *, CBLAS_INT *);

        /* The complex flavors take a leading UPLO the real ones do not have: for a Hermitian
         * matrix it says whether `e` is the sub- or superdiagonal, which for a symmetric one
         * is not a distinction. */
        void BLAS_FUNC(spttrs)(CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dpttrs)(CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cpttrs)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zpttrs)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(spteqr)(char *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dpteqr)(char *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cpteqr)(char *, CBLAS_INT *, f32 *, f32 *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zpteqr)(char *, CBLAS_INT *, f64 *, f64 *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sptsvx)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dptsvx)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cptsvx)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, c64 *, f32 *, c64 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, f32 *, c64 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zptsvx)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, c128 *, f64 *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, f64 *, c128 *, f64 *, CBLAS_INT *);

        /* The symmetric/Hermitian indefinite group.  `?sy*` exists for all four flavors -- a
         * complex *symmetric* matrix (A**T = A) is a real LAPACK subject, distinct from the
         * Hermitian `?he*` (A**H = A), which exists only in c/z.  So `csysv` and `chesv` are
         * two different routines, not two spellings of one. */
        void BLAS_FUNC(ssytrf)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsytrf)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(csytrf)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zsytrf)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(chetrf)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zhetrf)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(ssytf2)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsytf2)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(csytf2)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zsytf2)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(ssytrs)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsytrs)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(csytrs)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zsytrs)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(chetrs)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zhetrs)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        /* `csytri` wants 2n of workspace where the others want n -- see the wrapper. */
        void BLAS_FUNC(ssytri)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dsytri)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(csytri)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zsytri)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *);
        void BLAS_FUNC(chetri)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zhetri)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(ssyconv)(char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dsyconv)(char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(csyconv)(char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zsyconv)(char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *);

        /* `s`, `scond` and `amax` are magnitudes, so real; `work` follows the flavor. */
        void BLAS_FUNC(ssyequb)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dsyequb)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(csyequb)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, f32 *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zsyequb)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, f64 *, c128 *, CBLAS_INT *);
        void BLAS_FUNC(cheequb)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, f32 *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zheequb)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, f64 *, c128 *, CBLAS_INT *);

        /* Real `?sycon` takes IWORK; the complex `?sycon`/`?hecon` do not -- the `?gtcon` split. */
        void BLAS_FUNC(ssycon)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsycon)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(csycon)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zsycon)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, c128 *, CBLAS_INT *);
        void BLAS_FUNC(checon)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zhecon)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(ssysv)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsysv)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(csysv)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zsysv)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(chesv)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zhesv)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        /* `?sysvx`'s last workspace is IWORK for the real flavors and RWORK for the complex --
         * same slot, same length `n`, different element type (the `?gecon` shape).  `?hesvx`
         * is complex only, so its RWORK needs no such alias. */
        void BLAS_FUNC(ssysvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsysvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(csysvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, f32 *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zsysvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, f64 *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(chesvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, f32 *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zhesvx)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, f64 *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        /* The symmetric/Hermitian eigensolvers.  `w` is real in every flavor, and the complex
         * routines add an RWORK the real ones do not have -- sometimes hidden, sometimes with
         * its own `lrwork` in the Python signature (see `heevd`, `heevr`, `hegvd`). */
        void BLAS_FUNC(ssyev)(char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsyev)(char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cheev)(char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zheev)(char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(ssyevd)(char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsyevd)(char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cheevd)(char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zheevd)(char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(ssyevr)(char *, char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsyevr)(char *, char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cheevr)(char *, char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zheevr)(char *, char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(ssyevx)(char *, char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsyevx)(char *, char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cheevx)(char *, char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zheevx)(char *, char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(ssygv)(CBLAS_INT *, char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsygv)(CBLAS_INT *, char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(chegv)(CBLAS_INT *, char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zhegv)(CBLAS_INT *, char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(ssygvd)(CBLAS_INT *, char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsygvd)(CBLAS_INT *, char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(chegvd)(CBLAS_INT *, char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zhegvd)(CBLAS_INT *, char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(ssygvx)(CBLAS_INT *, char *, char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsygvx)(CBLAS_INT *, char *, char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(chegvx)(CBLAS_INT *, char *, char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zhegvx)(CBLAS_INT *, char *, char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        /* `d` and `e` are real even beside a complex `a`: they are the diagonal and
         * off-diagonal of a real symmetric tridiagonal form.  `tau` follows the flavor. */
        void BLAS_FUNC(ssytrd)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsytrd)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(chetrd)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zhetrd)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(ssygst)(CBLAS_INT *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsygst)(CBLAS_INT *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(chegst)(CBLAS_INT *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zhegst)(CBLAS_INT *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);
    }

    inline void gees(char jobvs, char sort, CBLAS_INT (*select)(f32 *, f32 *), CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *sdim, f32 *wr, f32 *wi, f32 *vs, CBLAS_INT ldvs, f32 *work, CBLAS_INT lwork, CBLAS_INT *bwork, CBLAS_INT *info)
        { BLAS_FUNC(sgees)(&jobvs, &sort, select, &n, a, &lda, sdim, wr, wi, vs, &ldvs, work, &lwork, bwork, info); }
    inline void gees(char jobvs, char sort, CBLAS_INT (*select)(f64 *, f64 *), CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *sdim, f64 *wr, f64 *wi, f64 *vs, CBLAS_INT ldvs, f64 *work, CBLAS_INT lwork, CBLAS_INT *bwork, CBLAS_INT *info)
        { BLAS_FUNC(dgees)(&jobvs, &sort, select, &n, a, &lda, sdim, wr, wi, vs, &ldvs, work, &lwork, bwork, info); }
    inline void gees(char jobvs, char sort, CBLAS_INT (*select)(c64 *), CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *sdim, c64 *w, c64 *vs, CBLAS_INT ldvs, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *bwork, CBLAS_INT *info)
        { BLAS_FUNC(cgees)(&jobvs, &sort, select, &n, a, &lda, sdim, w, vs, &ldvs, work, &lwork, rwork, bwork, info); }
    inline void gees(char jobvs, char sort, CBLAS_INT (*select)(c128 *), CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *sdim, c128 *w, c128 *vs, CBLAS_INT ldvs, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *bwork, CBLAS_INT *info)
        { BLAS_FUNC(zgees)(&jobvs, &sort, select, &n, a, &lda, sdim, w, vs, &ldvs, work, &lwork, rwork, bwork, info); }

    inline void gges(char jobvsl, char jobvsr, char sort, CBLAS_INT (*select)(f32 *, f32 *, f32 *), CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, CBLAS_INT *sdim, f32 *alphar, f32 *alphai, f32 *beta, f32 *vsl, CBLAS_INT ldvsl, f32 *vsr, CBLAS_INT ldvsr, f32 *work, CBLAS_INT lwork, CBLAS_INT *bwork, CBLAS_INT *info)
        { BLAS_FUNC(sgges)(&jobvsl, &jobvsr, &sort, select, &n, a, &lda, b, &ldb, sdim, alphar, alphai, beta, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, bwork, info); }
    inline void gges(char jobvsl, char jobvsr, char sort, CBLAS_INT (*select)(f64 *, f64 *, f64 *), CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, CBLAS_INT *sdim, f64 *alphar, f64 *alphai, f64 *beta, f64 *vsl, CBLAS_INT ldvsl, f64 *vsr, CBLAS_INT ldvsr, f64 *work, CBLAS_INT lwork, CBLAS_INT *bwork, CBLAS_INT *info)
        { BLAS_FUNC(dgges)(&jobvsl, &jobvsr, &sort, select, &n, a, &lda, b, &ldb, sdim, alphar, alphai, beta, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, bwork, info); }
    inline void gges(char jobvsl, char jobvsr, char sort, CBLAS_INT (*select)(c64 *, c64 *), CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, CBLAS_INT *sdim, c64 *alpha, c64 *beta, c64 *vsl, CBLAS_INT ldvsl, c64 *vsr, CBLAS_INT ldvsr, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *bwork, CBLAS_INT *info)
        { BLAS_FUNC(cgges)(&jobvsl, &jobvsr, &sort, select, &n, a, &lda, b, &ldb, sdim, alpha, beta, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, rwork, bwork, info); }
    inline void gges(char jobvsl, char jobvsr, char sort, CBLAS_INT (*select)(c128 *, c128 *), CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, CBLAS_INT *sdim, c128 *alpha, c128 *beta, c128 *vsl, CBLAS_INT ldvsl, c128 *vsr, CBLAS_INT ldvsr, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *bwork, CBLAS_INT *info)
        { BLAS_FUNC(zgges)(&jobvsl, &jobvsr, &sort, select, &n, a, &lda, b, &ldb, sdim, alpha, beta, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, rwork, bwork, info); }

    inline void gebal(char job, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *ilo, CBLAS_INT *ihi, f32 *scale, CBLAS_INT *info)
        { BLAS_FUNC(sgebal)(&job, &n, a, &lda, ilo, ihi, scale, info); }
    inline void gebal(char job, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *ilo, CBLAS_INT *ihi, f64 *scale, CBLAS_INT *info)
        { BLAS_FUNC(dgebal)(&job, &n, a, &lda, ilo, ihi, scale, info); }
    inline void gebal(char job, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *ilo, CBLAS_INT *ihi, f32 *scale, CBLAS_INT *info)
        { BLAS_FUNC(cgebal)(&job, &n, a, &lda, ilo, ihi, scale, info); }
    inline void gebal(char job, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *ilo, CBLAS_INT *ihi, f64 *scale, CBLAS_INT *info)
        { BLAS_FUNC(zgebal)(&job, &n, a, &lda, ilo, ihi, scale, info); }

    inline void gecon(char norm, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 anorm, f32 *rcond, f32 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(sgecon)(&norm, &n, a, &lda, &anorm, rcond, work, iwork, info); }
    inline void gecon(char norm, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 anorm, f64 *rcond, f64 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dgecon)(&norm, &n, a, &lda, &anorm, rcond, work, iwork, info); }
    inline void gecon(char norm, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 anorm, f32 *rcond, c64 *work, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cgecon)(&norm, &n, a, &lda, &anorm, rcond, work, rwork, info); }
    inline void gecon(char norm, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 anorm, f64 *rcond, c128 *work, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zgecon)(&norm, &n, a, &lda, &anorm, rcond, work, rwork, info); }

    inline void gesv(CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(sgesv)(&n, &nrhs, a, &lda, ipiv, b, &ldb, info); }
    inline void gesv(CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dgesv)(&n, &nrhs, a, &lda, ipiv, b, &ldb, info); }
    inline void gesv(CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cgesv)(&n, &nrhs, a, &lda, ipiv, b, &ldb, info); }
    inline void gesv(CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zgesv)(&n, &nrhs, a, &lda, ipiv, b, &ldb, info); }

    inline void getrf(CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(sgetrf)(&m, &n, a, &lda, ipiv, info); }
    inline void getrf(CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(dgetrf)(&m, &n, a, &lda, ipiv, info); }
    inline void getrf(CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(cgetrf)(&m, &n, a, &lda, ipiv, info); }
    inline void getrf(CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(zgetrf)(&m, &n, a, &lda, ipiv, info); }

    inline void getrs(char trans, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(sgetrs)(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, info); }
    inline void getrs(char trans, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dgetrs)(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, info); }
    inline void getrs(char trans, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cgetrs)(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, info); }
    inline void getrs(char trans, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zgetrs)(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, info); }

    inline void getc2(CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *ipiv, CBLAS_INT *jpiv, CBLAS_INT *info)
        { BLAS_FUNC(sgetc2)(&n, a, &lda, ipiv, jpiv, info); }
    inline void getc2(CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, CBLAS_INT *jpiv, CBLAS_INT *info)
        { BLAS_FUNC(dgetc2)(&n, a, &lda, ipiv, jpiv, info); }
    inline void getc2(CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, CBLAS_INT *jpiv, CBLAS_INT *info)
        { BLAS_FUNC(cgetc2)(&n, a, &lda, ipiv, jpiv, info); }
    inline void getc2(CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, CBLAS_INT *jpiv, CBLAS_INT *info)
        { BLAS_FUNC(zgetc2)(&n, a, &lda, ipiv, jpiv, info); }

    inline void gesc2(CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *rhs, CBLAS_INT *ipiv, CBLAS_INT *jpiv, f32 *scale)
        { BLAS_FUNC(sgesc2)(&n, a, &lda, rhs, ipiv, jpiv, scale); }
    inline void gesc2(CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *rhs, CBLAS_INT *ipiv, CBLAS_INT *jpiv, f64 *scale)
        { BLAS_FUNC(dgesc2)(&n, a, &lda, rhs, ipiv, jpiv, scale); }
    inline void gesc2(CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *rhs, CBLAS_INT *ipiv, CBLAS_INT *jpiv, f32 *scale)
        { BLAS_FUNC(cgesc2)(&n, a, &lda, rhs, ipiv, jpiv, scale); }
    inline void gesc2(CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *rhs, CBLAS_INT *ipiv, CBLAS_INT *jpiv, f64 *scale)
        { BLAS_FUNC(zgesc2)(&n, a, &lda, rhs, ipiv, jpiv, scale); }

    inline void getri(CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sgetri)(&n, a, &lda, ipiv, work, &lwork, info); }
    inline void getri(CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dgetri)(&n, a, &lda, ipiv, work, &lwork, info); }
    inline void getri(CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(cgetri)(&n, a, &lda, ipiv, work, &lwork, info); }
    inline void getri(CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zgetri)(&n, a, &lda, ipiv, work, &lwork, info); }

    inline void gesdd(char jobz, CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *s, f32 *u, CBLAS_INT ldu, f32 *vt, CBLAS_INT ldvt, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(sgesdd)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, info); }
    inline void gesdd(char jobz, CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *s, f64 *u, CBLAS_INT ldu, f64 *vt, CBLAS_INT ldvt, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dgesdd)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, info); }
    inline void gesdd(char jobz, CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *s, c64 *u, CBLAS_INT ldu, c64 *vt, CBLAS_INT ldvt, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(cgesdd)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, info); }
    inline void gesdd(char jobz, CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *s, c128 *u, CBLAS_INT ldu, c128 *vt, CBLAS_INT ldvt, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(zgesdd)(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, info); }

    inline void gesvd(char jobu, char jobvt, CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *s, f32 *u, CBLAS_INT ldu, f32 *vt, CBLAS_INT ldvt, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, info); }
    inline void gesvd(char jobu, char jobvt, CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *s, f64 *u, CBLAS_INT ldu, f64 *vt, CBLAS_INT ldvt, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, info); }
    inline void gesvd(char jobu, char jobvt, CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *s, c64 *u, CBLAS_INT ldu, c64 *vt, CBLAS_INT ldvt, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, info); }
    inline void gesvd(char jobu, char jobvt, CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *s, c128 *u, CBLAS_INT ldu, c128 *vt, CBLAS_INT ldvt, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, info); }

    inline void gels(char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sgels)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, info); }
    inline void gels(char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dgels)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, info); }
    inline void gels(char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(cgels)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, info); }
    inline void gels(char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zgels)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, info); }

    inline void geqp3(CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *jpvt, f32 *tau, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sgeqp3)(&m, &n, a, &lda, jpvt, tau, work, &lwork, info); }
    inline void geqp3(CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *jpvt, f64 *tau, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dgeqp3)(&m, &n, a, &lda, jpvt, tau, work, &lwork, info); }
    inline void geqp3(CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *jpvt, c64 *tau, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cgeqp3)(&m, &n, a, &lda, jpvt, tau, work, &lwork, rwork, info); }
    inline void geqp3(CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *jpvt, c128 *tau, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zgeqp3)(&m, &n, a, &lda, jpvt, tau, work, &lwork, rwork, info); }

    inline void gelsy(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, CBLAS_INT *jpvt, f32 cond, CBLAS_INT *rank, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sgelsy)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &cond, rank, work, &lwork, info); }
    inline void gelsy(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, CBLAS_INT *jpvt, f64 cond, CBLAS_INT *rank, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dgelsy)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &cond, rank, work, &lwork, info); }
    inline void gelsy(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, CBLAS_INT *jpvt, f32 cond, CBLAS_INT *rank, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cgelsy)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &cond, rank, work, &lwork, rwork, info); }
    inline void gelsy(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, CBLAS_INT *jpvt, f64 cond, CBLAS_INT *rank, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zgelsy)(&m, &n, &nrhs, a, &lda, b, &ldb, jpvt, &cond, rank, work, &lwork, rwork, info); }

    inline void gelsd(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *s, f32 cond, CBLAS_INT *rank, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(sgelsd)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &cond, rank, work, &lwork, iwork, info); }
    inline void gelsd(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *s, f64 cond, CBLAS_INT *rank, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dgelsd)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &cond, rank, work, &lwork, iwork, info); }
    inline void gelsd(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, f32 *s, f32 cond, CBLAS_INT *rank, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(cgelsd)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &cond, rank, work, &lwork, rwork, iwork, info); }
    inline void gelsd(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, f64 *s, f64 cond, CBLAS_INT *rank, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(zgelsd)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &cond, rank, work, &lwork, rwork, iwork, info); }

    inline void gelss(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *s, f32 cond, CBLAS_INT *rank, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sgelss)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &cond, rank, work, &lwork, info); }
    inline void gelss(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *s, f64 cond, CBLAS_INT *rank, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dgelss)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &cond, rank, work, &lwork, info); }
    inline void gelss(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, f32 *s, f32 cond, CBLAS_INT *rank, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cgelss)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &cond, rank, work, &lwork, rwork, info); }
    inline void gelss(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, f64 *s, f64 cond, CBLAS_INT *rank, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zgelss)(&m, &n, &nrhs, a, &lda, b, &ldb, s, &cond, rank, work, &lwork, rwork, info); }

    inline void gehrd(CBLAS_INT n, CBLAS_INT ilo, CBLAS_INT ihi, f32 *a, CBLAS_INT lda, f32 *tau, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sgehrd)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, info); }
    inline void gehrd(CBLAS_INT n, CBLAS_INT ilo, CBLAS_INT ihi, f64 *a, CBLAS_INT lda, f64 *tau, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dgehrd)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, info); }
    inline void gehrd(CBLAS_INT n, CBLAS_INT ilo, CBLAS_INT ihi, c64 *a, CBLAS_INT lda, c64 *tau, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(cgehrd)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, info); }
    inline void gehrd(CBLAS_INT n, CBLAS_INT ilo, CBLAS_INT ihi, c128 *a, CBLAS_INT lda, c128 *tau, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zgehrd)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, info); }

    inline void geqrf(CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *tau, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sgeqrf)(&m, &n, a, &lda, tau, work, &lwork, info); }
    inline void geqrf(CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *tau, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, info); }
    inline void geqrf(CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *tau, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(cgeqrf)(&m, &n, a, &lda, tau, work, &lwork, info); }
    inline void geqrf(CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *tau, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zgeqrf)(&m, &n, a, &lda, tau, work, &lwork, info); }

    inline void geqrfp(CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *tau, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sgeqrfp)(&m, &n, a, &lda, tau, work, &lwork, info); }
    inline void geqrfp(CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *tau, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dgeqrfp)(&m, &n, a, &lda, tau, work, &lwork, info); }
    inline void geqrfp(CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *tau, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(cgeqrfp)(&m, &n, a, &lda, tau, work, &lwork, info); }
    inline void geqrfp(CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *tau, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zgeqrfp)(&m, &n, a, &lda, tau, work, &lwork, info); }

    inline void gerqf(CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *tau, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sgerqf)(&m, &n, a, &lda, tau, work, &lwork, info); }
    inline void gerqf(CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *tau, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dgerqf)(&m, &n, a, &lda, tau, work, &lwork, info); }
    inline void gerqf(CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *tau, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(cgerqf)(&m, &n, a, &lda, tau, work, &lwork, info); }
    inline void gerqf(CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *tau, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zgerqf)(&m, &n, a, &lda, tau, work, &lwork, info); }

    inline void geequ(CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *r, f32 *c, f32 *rowcnd, f32 *colcnd, f32 *amax, CBLAS_INT *info)
        { BLAS_FUNC(sgeequ)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, info); }
    inline void geequ(CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *r, f64 *c, f64 *rowcnd, f64 *colcnd, f64 *amax, CBLAS_INT *info)
        { BLAS_FUNC(dgeequ)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, info); }
    inline void geequ(CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *r, f32 *c, f32 *rowcnd, f32 *colcnd, f32 *amax, CBLAS_INT *info)
        { BLAS_FUNC(cgeequ)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, info); }
    inline void geequ(CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *r, f64 *c, f64 *rowcnd, f64 *colcnd, f64 *amax, CBLAS_INT *info)
        { BLAS_FUNC(zgeequ)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, info); }

    inline void geequb(CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *r, f32 *c, f32 *rowcnd, f32 *colcnd, f32 *amax, CBLAS_INT *info)
        { BLAS_FUNC(sgeequb)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, info); }
    inline void geequb(CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *r, f64 *c, f64 *rowcnd, f64 *colcnd, f64 *amax, CBLAS_INT *info)
        { BLAS_FUNC(dgeequb)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, info); }
    inline void geequb(CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *r, f32 *c, f32 *rowcnd, f32 *colcnd, f32 *amax, CBLAS_INT *info)
        { BLAS_FUNC(cgeequb)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, info); }
    inline void geequb(CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *r, f64 *c, f64 *rowcnd, f64 *colcnd, f64 *amax, CBLAS_INT *info)
        { BLAS_FUNC(zgeequb)(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, info); }

    inline void geev(char jobvl, char jobvr, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *wr, f32 *wi, f32 *vl, CBLAS_INT ldvl, f32 *vr, CBLAS_INT ldvr, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sgeev)(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, info); }
    inline void geev(char jobvl, char jobvr, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *wr, f64 *wi, f64 *vl, CBLAS_INT ldvl, f64 *vr, CBLAS_INT ldvr, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dgeev)(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, info); }
    inline void geev(char jobvl, char jobvr, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *w, c64 *vl, CBLAS_INT ldvl, c64 *vr, CBLAS_INT ldvr, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cgeev)(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, info); }
    inline void geev(char jobvl, char jobvr, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *w, c128 *vl, CBLAS_INT ldvl, c128 *vr, CBLAS_INT ldvr, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zgeev)(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, info); }

    inline void ggev(char jobvl, char jobvr, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *alphar, f32 *alphai, f32 *beta, f32 *vl, CBLAS_INT ldvl, f32 *vr, CBLAS_INT ldvr, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sggev)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, work, &lwork, info); }
    inline void ggev(char jobvl, char jobvr, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *alphar, f64 *alphai, f64 *beta, f64 *vl, CBLAS_INT ldvl, f64 *vr, CBLAS_INT ldvr, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dggev)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta, vl, &ldvl, vr, &ldvr, work, &lwork, info); }
    inline void ggev(char jobvl, char jobvr, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, c64 *alpha, c64 *beta, c64 *vl, CBLAS_INT ldvl, c64 *vr, CBLAS_INT ldvr, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cggev)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, info); }
    inline void ggev(char jobvl, char jobvr, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, c128 *alpha, c128 *beta, c128 *vl, CBLAS_INT ldvl, c128 *vr, CBLAS_INT ldvr, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zggev)(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, info); }

    inline void gesvx(char fact, char trans, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, f32 *af, CBLAS_INT ldaf, CBLAS_INT *ipiv, char *equed, f32 *r, f32 *c, f32 *b, CBLAS_INT ldb, f32 *x, CBLAS_INT ldx, f32 *rcond, f32 *ferr, f32 *berr, f32 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(sgesvx)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldb, x, &ldx, rcond, ferr, berr, work, iwork, info); }
    inline void gesvx(char fact, char trans, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, f64 *af, CBLAS_INT ldaf, CBLAS_INT *ipiv, char *equed, f64 *r, f64 *c, f64 *b, CBLAS_INT ldb, f64 *x, CBLAS_INT ldx, f64 *rcond, f64 *ferr, f64 *berr, f64 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dgesvx)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldb, x, &ldx, rcond, ferr, berr, work, iwork, info); }
    inline void gesvx(char fact, char trans, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, c64 *af, CBLAS_INT ldaf, CBLAS_INT *ipiv, char *equed, f32 *r, f32 *c, c64 *b, CBLAS_INT ldb, c64 *x, CBLAS_INT ldx, f32 *rcond, f32 *ferr, f32 *berr, c64 *work, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cgesvx)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, info); }
    inline void gesvx(char fact, char trans, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, c128 *af, CBLAS_INT ldaf, CBLAS_INT *ipiv, char *equed, f64 *r, f64 *c, c128 *b, CBLAS_INT ldb, c128 *x, CBLAS_INT ldx, f64 *rcond, f64 *ferr, f64 *berr, c128 *work, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zgesvx)(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, equed, r, c, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, info); }

    inline void gtsv(CBLAS_INT n, CBLAS_INT nrhs, f32 *dl, f32 *d, f32 *du, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(sgtsv)(&n, &nrhs, dl, d, du, b, &ldb, info); }
    inline void gtsv(CBLAS_INT n, CBLAS_INT nrhs, f64 *dl, f64 *d, f64 *du, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dgtsv)(&n, &nrhs, dl, d, du, b, &ldb, info); }
    inline void gtsv(CBLAS_INT n, CBLAS_INT nrhs, c64 *dl, c64 *d, c64 *du, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cgtsv)(&n, &nrhs, dl, d, du, b, &ldb, info); }
    inline void gtsv(CBLAS_INT n, CBLAS_INT nrhs, c128 *dl, c128 *d, c128 *du, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zgtsv)(&n, &nrhs, dl, d, du, b, &ldb, info); }

    inline void gttrf(CBLAS_INT n, f32 *dl, f32 *d, f32 *du, f32 *du2, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(sgttrf)(&n, dl, d, du, du2, ipiv, info); }
    inline void gttrf(CBLAS_INT n, f64 *dl, f64 *d, f64 *du, f64 *du2, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(dgttrf)(&n, dl, d, du, du2, ipiv, info); }
    inline void gttrf(CBLAS_INT n, c64 *dl, c64 *d, c64 *du, c64 *du2, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(cgttrf)(&n, dl, d, du, du2, ipiv, info); }
    inline void gttrf(CBLAS_INT n, c128 *dl, c128 *d, c128 *du, c128 *du2, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(zgttrf)(&n, dl, d, du, du2, ipiv, info); }

    inline void gttrs(char trans, CBLAS_INT n, CBLAS_INT nrhs, f32 *dl, f32 *d, f32 *du, f32 *du2, CBLAS_INT *ipiv, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(sgttrs)(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, info); }
    inline void gttrs(char trans, CBLAS_INT n, CBLAS_INT nrhs, f64 *dl, f64 *d, f64 *du, f64 *du2, CBLAS_INT *ipiv, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dgttrs)(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, info); }
    inline void gttrs(char trans, CBLAS_INT n, CBLAS_INT nrhs, c64 *dl, c64 *d, c64 *du, c64 *du2, CBLAS_INT *ipiv, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cgttrs)(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, info); }
    inline void gttrs(char trans, CBLAS_INT n, CBLAS_INT nrhs, c128 *dl, c128 *d, c128 *du, c128 *du2, CBLAS_INT *ipiv, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zgttrs)(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, info); }

    inline void gtcon(char norm, CBLAS_INT n, f32 *dl, f32 *d, f32 *du, f32 *du2, CBLAS_INT *ipiv, f32 anorm, f32 *rcond, f32 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(sgtcon)(&norm, &n, dl, d, du, du2, ipiv, &anorm, rcond, work, iwork, info); }
    inline void gtcon(char norm, CBLAS_INT n, f64 *dl, f64 *d, f64 *du, f64 *du2, CBLAS_INT *ipiv, f64 anorm, f64 *rcond, f64 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dgtcon)(&norm, &n, dl, d, du, du2, ipiv, &anorm, rcond, work, iwork, info); }
    inline void gtcon(char norm, CBLAS_INT n, c64 *dl, c64 *d, c64 *du, c64 *du2, CBLAS_INT *ipiv, f32 anorm, f32 *rcond, c64 *work, CBLAS_INT *info)
        { BLAS_FUNC(cgtcon)(&norm, &n, dl, d, du, du2, ipiv, &anorm, rcond, work, info); }
    inline void gtcon(char norm, CBLAS_INT n, c128 *dl, c128 *d, c128 *du, c128 *du2, CBLAS_INT *ipiv, f64 anorm, f64 *rcond, c128 *work, CBLAS_INT *info)
        { BLAS_FUNC(zgtcon)(&norm, &n, dl, d, du, du2, ipiv, &anorm, rcond, work, info); }

    inline void gtsvx(char fact, char trans, CBLAS_INT n, CBLAS_INT nrhs, f32 *dl, f32 *d, f32 *du, f32 *dlf, f32 *df, f32 *duf, f32 *du2, CBLAS_INT *ipiv, f32 *b, CBLAS_INT ldb, f32 *x, CBLAS_INT ldx, f32 *rcond, f32 *ferr, f32 *berr, f32 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(sgtsvx)(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, iwork, info); }
    inline void gtsvx(char fact, char trans, CBLAS_INT n, CBLAS_INT nrhs, f64 *dl, f64 *d, f64 *du, f64 *dlf, f64 *df, f64 *duf, f64 *du2, CBLAS_INT *ipiv, f64 *b, CBLAS_INT ldb, f64 *x, CBLAS_INT ldx, f64 *rcond, f64 *ferr, f64 *berr, f64 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dgtsvx)(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, iwork, info); }
    inline void gtsvx(char fact, char trans, CBLAS_INT n, CBLAS_INT nrhs, c64 *dl, c64 *d, c64 *du, c64 *dlf, c64 *df, c64 *duf, c64 *du2, CBLAS_INT *ipiv, c64 *b, CBLAS_INT ldb, c64 *x, CBLAS_INT ldx, f32 *rcond, f32 *ferr, f32 *berr, c64 *work, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cgtsvx)(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, info); }
    inline void gtsvx(char fact, char trans, CBLAS_INT n, CBLAS_INT nrhs, c128 *dl, c128 *d, c128 *du, c128 *dlf, c128 *df, c128 *duf, c128 *du2, CBLAS_INT *ipiv, c128 *b, CBLAS_INT ldb, c128 *x, CBLAS_INT ldx, f64 *rcond, f64 *ferr, f64 *berr, c128 *work, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zgtsvx)(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, info); }

    inline void stev(char jobz, CBLAS_INT n, f32 *d, f32 *e, f32 *z, CBLAS_INT ldz, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(sstev)(&jobz, &n, d, e, z, &ldz, work, info); }
    inline void stev(char jobz, CBLAS_INT n, f64 *d, f64 *e, f64 *z, CBLAS_INT ldz, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dstev)(&jobz, &n, d, e, z, &ldz, work, info); }

    inline void stebz(char range, char order, CBLAS_INT n, f32 vl, f32 vu, CBLAS_INT il, CBLAS_INT iu, f32 abstol, f32 *d, f32 *e, CBLAS_INT *m, CBLAS_INT *nsplit, f32 *w, CBLAS_INT *iblock, CBLAS_INT *isplit, f32 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(sstebz)(&range, &order, &n, &vl, &vu, &il, &iu, &abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info); }
    inline void stebz(char range, char order, CBLAS_INT n, f64 vl, f64 vu, CBLAS_INT il, CBLAS_INT iu, f64 abstol, f64 *d, f64 *e, CBLAS_INT *m, CBLAS_INT *nsplit, f64 *w, CBLAS_INT *iblock, CBLAS_INT *isplit, f64 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dstebz)(&range, &order, &n, &vl, &vu, &il, &iu, &abstol, d, e, m, nsplit, w, iblock, isplit, work, iwork, info); }

    inline void sterf(CBLAS_INT n, f32 *d, f32 *e, CBLAS_INT *info)
        { BLAS_FUNC(ssterf)(&n, d, e, info); }
    inline void sterf(CBLAS_INT n, f64 *d, f64 *e, CBLAS_INT *info)
        { BLAS_FUNC(dsterf)(&n, d, e, info); }

    inline void stein(CBLAS_INT n, f32 *d, f32 *e, CBLAS_INT m, f32 *w, CBLAS_INT *iblock, CBLAS_INT *isplit, f32 *z, CBLAS_INT ldz, f32 *work, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(sstein)(&n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork, ifail, info); }
    inline void stein(CBLAS_INT n, f64 *d, f64 *e, CBLAS_INT m, f64 *w, CBLAS_INT *iblock, CBLAS_INT *isplit, f64 *z, CBLAS_INT ldz, f64 *work, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(dstein)(&n, d, e, &m, w, iblock, isplit, z, &ldz, work, iwork, ifail, info); }

    inline void stemr(char jobz, char range, CBLAS_INT n, f32 *d, f32 *e, f32 vl, f32 vu, CBLAS_INT il, CBLAS_INT iu, CBLAS_INT *m, f32 *w, f32 *z, CBLAS_INT ldz, CBLAS_INT nzc, CBLAS_INT *isuppz, CBLAS_INT *tryrac, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(sstemr)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, m, w, z, &ldz, &nzc, isuppz, tryrac, work, &lwork, iwork, &liwork, info); }
    inline void stemr(char jobz, char range, CBLAS_INT n, f64 *d, f64 *e, f64 vl, f64 vu, CBLAS_INT il, CBLAS_INT iu, CBLAS_INT *m, f64 *w, f64 *z, CBLAS_INT ldz, CBLAS_INT nzc, CBLAS_INT *isuppz, CBLAS_INT *tryrac, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(dstemr)(&jobz, &range, &n, d, e, &vl, &vu, &il, &iu, m, w, z, &ldz, &nzc, isuppz, tryrac, work, &lwork, iwork, &liwork, info); }

    inline void stevd(char jobz, CBLAS_INT n, f32 *d, f32 *e, f32 *z, CBLAS_INT ldz, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(sstevd)(&jobz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, info); }
    inline void stevd(char jobz, CBLAS_INT n, f64 *d, f64 *e, f64 *z, CBLAS_INT ldz, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(dstevd)(&jobz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, info); }

    inline void gbsv(CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, CBLAS_INT nrhs, f32 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(sgbsv)(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, info); }
    inline void gbsv(CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, CBLAS_INT nrhs, f64 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dgbsv)(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, info); }
    inline void gbsv(CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, CBLAS_INT nrhs, c64 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cgbsv)(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, info); }
    inline void gbsv(CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, CBLAS_INT nrhs, c128 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zgbsv)(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, info); }

    inline void gbtrf(CBLAS_INT m, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, f32 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(sgbtrf)(&m, &n, &kl, &ku, ab, &ldab, ipiv, info); }
    inline void gbtrf(CBLAS_INT m, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, f64 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(dgbtrf)(&m, &n, &kl, &ku, ab, &ldab, ipiv, info); }
    inline void gbtrf(CBLAS_INT m, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, c64 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(cgbtrf)(&m, &n, &kl, &ku, ab, &ldab, ipiv, info); }
    inline void gbtrf(CBLAS_INT m, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, c128 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(zgbtrf)(&m, &n, &kl, &ku, ab, &ldab, ipiv, info); }

    inline void gbtrs(char trans, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, CBLAS_INT nrhs, f32 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(sgbtrs)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, info); }
    inline void gbtrs(char trans, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, CBLAS_INT nrhs, f64 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dgbtrs)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, info); }
    inline void gbtrs(char trans, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, CBLAS_INT nrhs, c64 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cgbtrs)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, info); }
    inline void gbtrs(char trans, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, CBLAS_INT nrhs, c128 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zgbtrs)(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, info); }

    inline void gbcon(char norm, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, f32 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, f32 anorm, f32 *rcond, f32 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(sgbcon)(&norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, rcond, work, iwork, info); }
    inline void gbcon(char norm, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, f64 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, f64 anorm, f64 *rcond, f64 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dgbcon)(&norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, rcond, work, iwork, info); }
    inline void gbcon(char norm, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, c64 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, f32 anorm, f32 *rcond, c64 *work, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cgbcon)(&norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, rcond, work, rwork, info); }
    inline void gbcon(char norm, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, c128 *ab, CBLAS_INT ldab, CBLAS_INT *ipiv, f64 anorm, f64 *rcond, c128 *work, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zgbcon)(&norm, &n, &kl, &ku, ab, &ldab, ipiv, &anorm, rcond, work, rwork, info); }

    inline f32 langb(char norm, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, f32 *ab, CBLAS_INT ldab, f32 *work)
        { return BLAS_FUNC(slangb)(&norm, &n, &kl, &ku, ab, &ldab, work); }
    inline f64 langb(char norm, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, f64 *ab, CBLAS_INT ldab, f64 *work)
        { return BLAS_FUNC(dlangb)(&norm, &n, &kl, &ku, ab, &ldab, work); }
    inline f32 langb(char norm, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, c64 *ab, CBLAS_INT ldab, f32 *work)
        { return BLAS_FUNC(clangb)(&norm, &n, &kl, &ku, ab, &ldab, work); }
    inline f64 langb(char norm, CBLAS_INT n, CBLAS_INT kl, CBLAS_INT ku, c128 *ab, CBLAS_INT ldab, f64 *work)
        { return BLAS_FUNC(zlangb)(&norm, &n, &kl, &ku, ab, &ldab, work); }

    inline void pstrf(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *piv, CBLAS_INT *rank, f32 tol, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(spstrf)(&uplo, &n, a, &lda, piv, rank, &tol, work, info); }
    inline void pstrf(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *piv, CBLAS_INT *rank, f64 tol, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dpstrf)(&uplo, &n, a, &lda, piv, rank, &tol, work, info); }
    inline void pstrf(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *piv, CBLAS_INT *rank, f32 tol, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(cpstrf)(&uplo, &n, a, &lda, piv, rank, &tol, work, info); }
    inline void pstrf(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *piv, CBLAS_INT *rank, f64 tol, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(zpstrf)(&uplo, &n, a, &lda, piv, rank, &tol, work, info); }

    inline void pstf2(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *piv, CBLAS_INT *rank, f32 tol, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(spstf2)(&uplo, &n, a, &lda, piv, rank, &tol, work, info); }
    inline void pstf2(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *piv, CBLAS_INT *rank, f64 tol, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dpstf2)(&uplo, &n, a, &lda, piv, rank, &tol, work, info); }
    inline void pstf2(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *piv, CBLAS_INT *rank, f32 tol, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(cpstf2)(&uplo, &n, a, &lda, piv, rank, &tol, work, info); }
    inline void pstf2(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *piv, CBLAS_INT *rank, f64 tol, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(zpstf2)(&uplo, &n, a, &lda, piv, rank, &tol, work, info); }

    inline void posv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(sposv)(&uplo, &n, &nrhs, a, &lda, b, &ldb, info); }
    inline void posv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dposv)(&uplo, &n, &nrhs, a, &lda, b, &ldb, info); }
    inline void posv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cposv)(&uplo, &n, &nrhs, a, &lda, b, &ldb, info); }
    inline void posv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zposv)(&uplo, &n, &nrhs, a, &lda, b, &ldb, info); }

    inline void posvx(char fact, char uplo, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, f32 *af, CBLAS_INT ldaf, char *equed, f32 *s, f32 *b, CBLAS_INT ldb, f32 *x, CBLAS_INT ldx, f32 *rcond, f32 *ferr, f32 *berr, f32 *work, CBLAS_INT *irwork, CBLAS_INT *info)
        { BLAS_FUNC(sposvx)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr, work, irwork, info); }
    inline void posvx(char fact, char uplo, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, f64 *af, CBLAS_INT ldaf, char *equed, f64 *s, f64 *b, CBLAS_INT ldb, f64 *x, CBLAS_INT ldx, f64 *rcond, f64 *ferr, f64 *berr, f64 *work, CBLAS_INT *irwork, CBLAS_INT *info)
        { BLAS_FUNC(dposvx)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr, work, irwork, info); }
    inline void posvx(char fact, char uplo, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, c64 *af, CBLAS_INT ldaf, char *equed, f32 *s, c64 *b, CBLAS_INT ldb, c64 *x, CBLAS_INT ldx, f32 *rcond, f32 *ferr, f32 *berr, c64 *work, f32 *irwork, CBLAS_INT *info)
        { BLAS_FUNC(cposvx)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr, work, irwork, info); }
    inline void posvx(char fact, char uplo, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, c128 *af, CBLAS_INT ldaf, char *equed, f64 *s, c128 *b, CBLAS_INT ldb, c128 *x, CBLAS_INT ldx, f64 *rcond, f64 *ferr, f64 *berr, c128 *work, f64 *irwork, CBLAS_INT *info)
        { BLAS_FUNC(zposvx)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, equed, s, b, &ldb, x, &ldx, rcond, ferr, berr, work, irwork, info); }

    inline void pocon(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 anorm, f32 *rcond, f32 *work, CBLAS_INT *irwork, CBLAS_INT *info)
        { BLAS_FUNC(spocon)(&uplo, &n, a, &lda, &anorm, rcond, work, irwork, info); }
    inline void pocon(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 anorm, f64 *rcond, f64 *work, CBLAS_INT *irwork, CBLAS_INT *info)
        { BLAS_FUNC(dpocon)(&uplo, &n, a, &lda, &anorm, rcond, work, irwork, info); }
    inline void pocon(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 anorm, f32 *rcond, c64 *work, f32 *irwork, CBLAS_INT *info)
        { BLAS_FUNC(cpocon)(&uplo, &n, a, &lda, &anorm, rcond, work, irwork, info); }
    inline void pocon(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 anorm, f64 *rcond, c128 *work, f64 *irwork, CBLAS_INT *info)
        { BLAS_FUNC(zpocon)(&uplo, &n, a, &lda, &anorm, rcond, work, irwork, info); }

    inline void potrf(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *info)   { BLAS_FUNC(spotrf)(&uplo, &n, a, &lda, info); }
    inline void potrf(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *info)   { BLAS_FUNC(dpotrf)(&uplo, &n, a, &lda, info); }
    inline void potrf(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *info)   { BLAS_FUNC(cpotrf)(&uplo, &n, a, &lda, info); }
    inline void potrf(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *info)  { BLAS_FUNC(zpotrf)(&uplo, &n, a, &lda, info); }

    inline void potrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(spotrs)(&uplo, &n, &nrhs, a, &lda, b, &ldb, info); }
    inline void potrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dpotrs)(&uplo, &n, &nrhs, a, &lda, b, &ldb, info); }
    inline void potrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cpotrs)(&uplo, &n, &nrhs, a, &lda, b, &ldb, info); }
    inline void potrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zpotrs)(&uplo, &n, &nrhs, a, &lda, b, &ldb, info); }

    inline void potri(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *info)   { BLAS_FUNC(spotri)(&uplo, &n, a, &lda, info); }
    inline void potri(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *info)   { BLAS_FUNC(dpotri)(&uplo, &n, a, &lda, info); }
    inline void potri(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *info)   { BLAS_FUNC(cpotri)(&uplo, &n, a, &lda, info); }
    inline void potri(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *info)  { BLAS_FUNC(zpotri)(&uplo, &n, a, &lda, info); }

    inline void ptsv(CBLAS_INT n, CBLAS_INT nrhs, f32 *d, f32 *e, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(sptsv)(&n, &nrhs, d, e, b, &ldb, info); }
    inline void ptsv(CBLAS_INT n, CBLAS_INT nrhs, f64 *d, f64 *e, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dptsv)(&n, &nrhs, d, e, b, &ldb, info); }
    inline void ptsv(CBLAS_INT n, CBLAS_INT nrhs, f32 *d, c64 *e, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cptsv)(&n, &nrhs, d, e, b, &ldb, info); }
    inline void ptsv(CBLAS_INT n, CBLAS_INT nrhs, f64 *d, c128 *e, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zptsv)(&n, &nrhs, d, e, b, &ldb, info); }

    inline void pttrf(CBLAS_INT n, f32 *d, f32 *e, CBLAS_INT *info)   { BLAS_FUNC(spttrf)(&n, d, e, info); }
    inline void pttrf(CBLAS_INT n, f64 *d, f64 *e, CBLAS_INT *info)   { BLAS_FUNC(dpttrf)(&n, d, e, info); }
    inline void pttrf(CBLAS_INT n, f32 *d, c64 *e, CBLAS_INT *info)   { BLAS_FUNC(cpttrf)(&n, d, e, info); }
    inline void pttrf(CBLAS_INT n, f64 *d, c128 *e, CBLAS_INT *info)  { BLAS_FUNC(zpttrf)(&n, d, e, info); }

    /* The real overloads drop `uplo` rather than ignore it, so a wrapper that passes one to a
     * real flavor fails to compile instead of silently discarding it. */
    inline void pttrs(CBLAS_INT n, CBLAS_INT nrhs, f32 *d, f32 *e, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(spttrs)(&n, &nrhs, d, e, b, &ldb, info); }
    inline void pttrs(CBLAS_INT n, CBLAS_INT nrhs, f64 *d, f64 *e, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dpttrs)(&n, &nrhs, d, e, b, &ldb, info); }
    inline void pttrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f32 *d, c64 *e, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cpttrs)(&uplo, &n, &nrhs, d, e, b, &ldb, info); }
    inline void pttrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f64 *d, c128 *e, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zpttrs)(&uplo, &n, &nrhs, d, e, b, &ldb, info); }

    inline void pteqr(char compz, CBLAS_INT n, f32 *d, f32 *e, f32 *z, CBLAS_INT ldz, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(spteqr)(&compz, &n, d, e, z, &ldz, work, info); }
    inline void pteqr(char compz, CBLAS_INT n, f64 *d, f64 *e, f64 *z, CBLAS_INT ldz, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dpteqr)(&compz, &n, d, e, z, &ldz, work, info); }
    inline void pteqr(char compz, CBLAS_INT n, f32 *d, f32 *e, c64 *z, CBLAS_INT ldz, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(cpteqr)(&compz, &n, d, e, z, &ldz, work, info); }
    inline void pteqr(char compz, CBLAS_INT n, f64 *d, f64 *e, c128 *z, CBLAS_INT ldz, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(zpteqr)(&compz, &n, d, e, z, &ldz, work, info); }

    inline void ptsvx(char fact, CBLAS_INT n, CBLAS_INT nrhs, f32 *d, f32 *e, f32 *df, f32 *ef, f32 *b, CBLAS_INT ldb, f32 *x, CBLAS_INT ldx, f32 *rcond, f32 *ferr, f32 *berr, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(sptsvx)(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, rcond, ferr, berr, work, info); }
    inline void ptsvx(char fact, CBLAS_INT n, CBLAS_INT nrhs, f64 *d, f64 *e, f64 *df, f64 *ef, f64 *b, CBLAS_INT ldb, f64 *x, CBLAS_INT ldx, f64 *rcond, f64 *ferr, f64 *berr, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dptsvx)(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, rcond, ferr, berr, work, info); }
    inline void ptsvx(char fact, CBLAS_INT n, CBLAS_INT nrhs, f32 *d, c64 *e, f32 *df, c64 *ef, c64 *b, CBLAS_INT ldb, c64 *x, CBLAS_INT ldx, f32 *rcond, f32 *ferr, f32 *berr, c64 *work, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cptsvx)(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, info); }
    inline void ptsvx(char fact, CBLAS_INT n, CBLAS_INT nrhs, f64 *d, c128 *e, f64 *df, c128 *ef, c128 *b, CBLAS_INT ldb, c128 *x, CBLAS_INT ldx, f64 *rcond, f64 *ferr, f64 *berr, c128 *work, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zptsvx)(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, rcond, ferr, berr, work, rwork, info); }

    inline void sytrf(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(ssytrf)(&uplo, &n, a, &lda, ipiv, work, &lwork, info); }
    inline void sytrf(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dsytrf)(&uplo, &n, a, &lda, ipiv, work, &lwork, info); }
    inline void sytrf(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(csytrf)(&uplo, &n, a, &lda, ipiv, work, &lwork, info); }
    inline void sytrf(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zsytrf)(&uplo, &n, a, &lda, ipiv, work, &lwork, info); }
    inline void hetrf(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(chetrf)(&uplo, &n, a, &lda, ipiv, work, &lwork, info); }
    inline void hetrf(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zhetrf)(&uplo, &n, a, &lda, ipiv, work, &lwork, info); }

    inline void sytf2(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(ssytf2)(&uplo, &n, a, &lda, ipiv, info); }
    inline void sytf2(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(dsytf2)(&uplo, &n, a, &lda, ipiv, info); }
    inline void sytf2(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(csytf2)(&uplo, &n, a, &lda, ipiv, info); }
    inline void sytf2(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, CBLAS_INT *info)
        { BLAS_FUNC(zsytf2)(&uplo, &n, a, &lda, ipiv, info); }

    inline void sytrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(ssytrs)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, info); }
    inline void sytrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dsytrs)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, info); }
    inline void sytrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(csytrs)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, info); }
    inline void sytrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zsytrs)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, info); }
    inline void hetrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(chetrs)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, info); }
    inline void hetrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zhetrs)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, info); }

    inline void sytri(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(ssytri)(&uplo, &n, a, &lda, ipiv, work, info); }
    inline void sytri(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dsytri)(&uplo, &n, a, &lda, ipiv, work, info); }
    inline void sytri(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c64 *work, CBLAS_INT *info)
        { BLAS_FUNC(csytri)(&uplo, &n, a, &lda, ipiv, work, info); }
    inline void sytri(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c128 *work, CBLAS_INT *info)
        { BLAS_FUNC(zsytri)(&uplo, &n, a, &lda, ipiv, work, info); }
    inline void hetri(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c64 *work, CBLAS_INT *info)
        { BLAS_FUNC(chetri)(&uplo, &n, a, &lda, ipiv, work, info); }
    inline void hetri(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c128 *work, CBLAS_INT *info)
        { BLAS_FUNC(zhetri)(&uplo, &n, a, &lda, ipiv, work, info); }

    inline void syconv(char uplo, char way, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f32 *e, CBLAS_INT *info)
        { BLAS_FUNC(ssyconv)(&uplo, &way, &n, a, &lda, ipiv, e, info); }
    inline void syconv(char uplo, char way, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f64 *e, CBLAS_INT *info)
        { BLAS_FUNC(dsyconv)(&uplo, &way, &n, a, &lda, ipiv, e, info); }
    inline void syconv(char uplo, char way, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c64 *e, CBLAS_INT *info)
        { BLAS_FUNC(csyconv)(&uplo, &way, &n, a, &lda, ipiv, e, info); }
    inline void syconv(char uplo, char way, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c128 *e, CBLAS_INT *info)
        { BLAS_FUNC(zsyconv)(&uplo, &way, &n, a, &lda, ipiv, e, info); }

    inline void syequb(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *s, f32 *scond, f32 *amax, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(ssyequb)(&uplo, &n, a, &lda, s, scond, amax, work, info); }
    inline void syequb(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *s, f64 *scond, f64 *amax, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dsyequb)(&uplo, &n, a, &lda, s, scond, amax, work, info); }
    inline void syequb(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *s, f32 *scond, f32 *amax, c64 *work, CBLAS_INT *info)
        { BLAS_FUNC(csyequb)(&uplo, &n, a, &lda, s, scond, amax, work, info); }
    inline void syequb(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *s, f64 *scond, f64 *amax, c128 *work, CBLAS_INT *info)
        { BLAS_FUNC(zsyequb)(&uplo, &n, a, &lda, s, scond, amax, work, info); }
    inline void heequb(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *s, f32 *scond, f32 *amax, c64 *work, CBLAS_INT *info)
        { BLAS_FUNC(cheequb)(&uplo, &n, a, &lda, s, scond, amax, work, info); }
    inline void heequb(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *s, f64 *scond, f64 *amax, c128 *work, CBLAS_INT *info)
        { BLAS_FUNC(zheequb)(&uplo, &n, a, &lda, s, scond, amax, work, info); }

    inline void sycon(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f32 anorm, f32 *rcond, f32 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(ssycon)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, iwork, info); }
    inline void sycon(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f64 anorm, f64 *rcond, f64 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dsycon)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, iwork, info); }
    inline void sycon(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f32 anorm, f32 *rcond, c64 *work, CBLAS_INT *info)
        { BLAS_FUNC(csycon)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, info); }
    inline void sycon(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f64 anorm, f64 *rcond, c128 *work, CBLAS_INT *info)
        { BLAS_FUNC(zsycon)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, info); }
    inline void hecon(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f32 anorm, f32 *rcond, c64 *work, CBLAS_INT *info)
        { BLAS_FUNC(checon)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, info); }
    inline void hecon(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f64 anorm, f64 *rcond, c128 *work, CBLAS_INT *info)
        { BLAS_FUNC(zhecon)(&uplo, &n, a, &lda, ipiv, &anorm, rcond, work, info); }

    inline void sysv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f32 *b, CBLAS_INT ldb, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(ssysv)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, info); }
    inline void sysv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, f64 *b, CBLAS_INT ldb, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dsysv)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, info); }
    inline void sysv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c64 *b, CBLAS_INT ldb, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(csysv)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, info); }
    inline void sysv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c128 *b, CBLAS_INT ldb, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zsysv)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, info); }
    inline void hesv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c64 *b, CBLAS_INT ldb, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(chesv)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, info); }
    inline void hesv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, CBLAS_INT *ipiv, c128 *b, CBLAS_INT ldb, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zhesv)(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, info); }

    inline void sysvx(char fact, char uplo, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, f32 *af, CBLAS_INT ldaf, CBLAS_INT *ipiv, f32 *b, CBLAS_INT ldb, f32 *x, CBLAS_INT ldx, f32 *rcond, f32 *ferr, f32 *berr, f32 *work, CBLAS_INT lwork, CBLAS_INT *irwork, CBLAS_INT *info)
        { BLAS_FUNC(ssysvx)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, &lwork, irwork, info); }
    inline void sysvx(char fact, char uplo, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, f64 *af, CBLAS_INT ldaf, CBLAS_INT *ipiv, f64 *b, CBLAS_INT ldb, f64 *x, CBLAS_INT ldx, f64 *rcond, f64 *ferr, f64 *berr, f64 *work, CBLAS_INT lwork, CBLAS_INT *irwork, CBLAS_INT *info)
        { BLAS_FUNC(dsysvx)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, &lwork, irwork, info); }
    inline void sysvx(char fact, char uplo, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, c64 *af, CBLAS_INT ldaf, CBLAS_INT *ipiv, c64 *b, CBLAS_INT ldb, c64 *x, CBLAS_INT ldx, f32 *rcond, f32 *ferr, f32 *berr, c64 *work, CBLAS_INT lwork, f32 *irwork, CBLAS_INT *info)
        { BLAS_FUNC(csysvx)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, &lwork, irwork, info); }
    inline void sysvx(char fact, char uplo, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, c128 *af, CBLAS_INT ldaf, CBLAS_INT *ipiv, c128 *b, CBLAS_INT ldb, c128 *x, CBLAS_INT ldx, f64 *rcond, f64 *ferr, f64 *berr, c128 *work, CBLAS_INT lwork, f64 *irwork, CBLAS_INT *info)
        { BLAS_FUNC(zsysvx)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, &lwork, irwork, info); }
    inline void hesvx(char fact, char uplo, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, c64 *af, CBLAS_INT ldaf, CBLAS_INT *ipiv, c64 *b, CBLAS_INT ldb, c64 *x, CBLAS_INT ldx, f32 *rcond, f32 *ferr, f32 *berr, c64 *work, CBLAS_INT lwork, f32 *irwork, CBLAS_INT *info)
        { BLAS_FUNC(chesvx)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, &lwork, irwork, info); }
    inline void hesvx(char fact, char uplo, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, c128 *af, CBLAS_INT ldaf, CBLAS_INT *ipiv, c128 *b, CBLAS_INT ldb, c128 *x, CBLAS_INT ldx, f64 *rcond, f64 *ferr, f64 *berr, c128 *work, CBLAS_INT lwork, f64 *irwork, CBLAS_INT *info)
        { BLAS_FUNC(zhesvx)(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, rcond, ferr, berr, work, &lwork, irwork, info); }

    inline void syev(char jobz, char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *w, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(ssyev)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, info); }
    inline void syev(char jobz, char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *w, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dsyev)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, info); }
    inline void syev(char jobz, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *w, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cheev)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, info); }
    inline void syev(char jobz, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *w, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zheev)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, info); }

    inline void syevd(char jobz, char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *w, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(ssyevd)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, info); }
    inline void syevd(char jobz, char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *w, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(dsyevd)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, info); }
    inline void heevd(char jobz, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *w, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(cheevd)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, info); }
    inline void heevd(char jobz, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *w, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(zheevd)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, info); }

    inline void syevr(char jobz, char range, char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 vl, f32 vu, CBLAS_INT il, CBLAS_INT iu, f32 abstol, CBLAS_INT *m, f32 *w, f32 *z, CBLAS_INT ldz, CBLAS_INT *isuppz, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(ssyevr)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, info); }
    inline void syevr(char jobz, char range, char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 vl, f64 vu, CBLAS_INT il, CBLAS_INT iu, f64 abstol, CBLAS_INT *m, f64 *w, f64 *z, CBLAS_INT ldz, CBLAS_INT *isuppz, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(dsyevr)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, info); }
    inline void heevr(char jobz, char range, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 vl, f32 vu, CBLAS_INT il, CBLAS_INT iu, f32 abstol, CBLAS_INT *m, f32 *w, c64 *z, CBLAS_INT ldz, CBLAS_INT *isuppz, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(cheevr)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork, rwork, &lrwork, iwork, &liwork, info); }
    inline void heevr(char jobz, char range, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 vl, f64 vu, CBLAS_INT il, CBLAS_INT iu, f64 abstol, CBLAS_INT *m, f64 *w, c128 *z, CBLAS_INT ldz, CBLAS_INT *isuppz, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(zheevr)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork, rwork, &lrwork, iwork, &liwork, info); }

    inline void syevx(char jobz, char range, char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 vl, f32 vu, CBLAS_INT il, CBLAS_INT iu, f32 abstol, CBLAS_INT *m, f32 *w, f32 *z, CBLAS_INT ldz, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(ssyevx)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork, iwork, ifail, info); }
    inline void syevx(char jobz, char range, char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 vl, f64 vu, CBLAS_INT il, CBLAS_INT iu, f64 abstol, CBLAS_INT *m, f64 *w, f64 *z, CBLAS_INT ldz, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(dsyevx)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork, iwork, ifail, info); }
    inline void syevx(char jobz, char range, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 vl, f32 vu, CBLAS_INT il, CBLAS_INT iu, f32 abstol, CBLAS_INT *m, f32 *w, c64 *z, CBLAS_INT ldz, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(cheevx)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork, rwork, iwork, ifail, info); }
    inline void syevx(char jobz, char range, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 vl, f64 vu, CBLAS_INT il, CBLAS_INT iu, f64 abstol, CBLAS_INT *m, f64 *w, c128 *z, CBLAS_INT ldz, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(zheevx)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork, rwork, iwork, ifail, info); }

    inline void sygv(CBLAS_INT itype, char jobz, char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *w, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(ssygv)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, info); }
    inline void sygv(CBLAS_INT itype, char jobz, char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *w, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dsygv)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, info); }
    inline void sygv(CBLAS_INT itype, char jobz, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, f32 *w, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(chegv)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, info); }
    inline void sygv(CBLAS_INT itype, char jobz, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, f64 *w, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zhegv)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, info); }

    inline void sygvd(CBLAS_INT itype, char jobz, char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *w, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(ssygvd)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, iwork, &liwork, info); }
    inline void sygvd(CBLAS_INT itype, char jobz, char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *w, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(dsygvd)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, iwork, &liwork, info); }
    inline void hegvd(CBLAS_INT itype, char jobz, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, f32 *w, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(chegvd)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &lrwork, iwork, &liwork, info); }
    inline void hegvd(CBLAS_INT itype, char jobz, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, f64 *w, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(zhegvd)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &lrwork, iwork, &liwork, info); }

    inline void sygvx(CBLAS_INT itype, char jobz, char range, char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 vl, f32 vu, CBLAS_INT il, CBLAS_INT iu, f32 abstol, CBLAS_INT *m, f32 *w, f32 *z, CBLAS_INT ldz, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(ssygvx)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork, iwork, ifail, info); }
    inline void sygvx(CBLAS_INT itype, char jobz, char range, char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 vl, f64 vu, CBLAS_INT il, CBLAS_INT iu, f64 abstol, CBLAS_INT *m, f64 *w, f64 *z, CBLAS_INT ldz, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(dsygvx)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork, iwork, ifail, info); }
    inline void sygvx(CBLAS_INT itype, char jobz, char range, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, f32 vl, f32 vu, CBLAS_INT il, CBLAS_INT iu, f32 abstol, CBLAS_INT *m, f32 *w, c64 *z, CBLAS_INT ldz, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(chegvx)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork, rwork, iwork, ifail, info); }
    inline void sygvx(CBLAS_INT itype, char jobz, char range, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, f64 vl, f64 vu, CBLAS_INT il, CBLAS_INT iu, f64 abstol, CBLAS_INT *m, f64 *w, c128 *z, CBLAS_INT ldz, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(zhegvx)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, &lwork, rwork, iwork, ifail, info); }

    inline void sytrd(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *d, f32 *e, f32 *tau, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(ssytrd)(&uplo, &n, a, &lda, d, e, tau, work, &lwork, info); }
    inline void sytrd(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *d, f64 *e, f64 *tau, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dsytrd)(&uplo, &n, a, &lda, d, e, tau, work, &lwork, info); }
    inline void sytrd(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *d, f32 *e, c64 *tau, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(chetrd)(&uplo, &n, a, &lda, d, e, tau, work, &lwork, info); }
    inline void sytrd(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *d, f64 *e, c128 *tau, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zhetrd)(&uplo, &n, a, &lda, d, e, tau, work, &lwork, info); }

    inline void sygst(CBLAS_INT itype, char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(ssygst)(&itype, &uplo, &n, a, &lda, b, &ldb, info); }
    inline void sygst(CBLAS_INT itype, char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dsygst)(&itype, &uplo, &n, a, &lda, b, &ldb, info); }
    inline void sygst(CBLAS_INT itype, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(chegst)(&itype, &uplo, &n, a, &lda, b, &ldb, info); }
    inline void sygst(CBLAS_INT itype, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zhegst)(&itype, &uplo, &n, a, &lda, b, &ldb, info); }

}  // namespace lapack
