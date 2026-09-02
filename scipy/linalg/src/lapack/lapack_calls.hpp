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

        void BLAS_FUNC(stpttf)(char *, char *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dtpttf)(char *, char *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(ctpttf)(char *, char *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(ztpttf)(char *, char *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(stpttr)(char *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dtpttr)(char *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ctpttr)(char *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ztpttr)(char *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(stfttp)(char *, char *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dtfttp)(char *, char *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(ctfttp)(char *, char *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(ztfttp)(char *, char *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(stfttr)(char *, char *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dtfttr)(char *, char *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ctfttr)(char *, char *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ztfttr)(char *, char *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(strttf)(char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dtrttf)(char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(ctrttf)(char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(ztrttf)(char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(strttp)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dtrttp)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(ctrttp)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(ztrttp)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(stfsm)(char *, char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dtfsm)(char *, char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(ctfsm)(char *, char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, c64 *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(ztfsm)(char *, char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(sppcon)(char *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dppcon)(char *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cppcon)(char *, CBLAS_INT *, c64 *, f32 *, f32 *, c64 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(zppcon)(char *, CBLAS_INT *, c128 *, f64 *, f64 *, c128 *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(sppsv)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dppsv)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cppsv)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zppsv)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(spptrf)(char *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dpptrf)(char *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cpptrf)(char *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zpptrf)(char *, CBLAS_INT *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(spptri)(char *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dpptri)(char *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cpptri)(char *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zpptri)(char *, CBLAS_INT *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(spptrs)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dpptrs)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cpptrs)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zpptrs)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(ssfrk)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, f32 *);
        void BLAS_FUNC(dsfrk)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, f64 *);
        void BLAS_FUNC(chfrk)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, f32 *, c64 *);
        void BLAS_FUNC(zhfrk)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, f64 *, c128 *);

        void BLAS_FUNC(spftrf)(char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dpftrf)(char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cpftrf)(char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zpftrf)(char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(spftri)(char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dpftri)(char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cpftri)(char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zpftri)(char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(spftrs)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dpftrs)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cpftrs)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zpftrs)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(spbtrf)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dpbtrf)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cpbtrf)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zpbtrf)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(spbtrs)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dpbtrs)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cpbtrs)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zpbtrs)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(spbsv)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dpbsv)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cpbsv)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zpbsv)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(strtrs)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dtrtrs)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ctrtrs)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ztrtrs)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(strcon)(char *, char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dtrcon)(char *, char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ctrcon)(char *, char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, c64 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(ztrcon)(char *, char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, c128 *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(stbtrs)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dtbtrs)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ctbtrs)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ztbtrs)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(strtri)(char *, char *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dtrtri)(char *, char *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ctrtri)(char *, char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ztrtri)(char *, char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(slauum)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dlauum)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(clauum)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zlauum)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(slaswp)(CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dlaswp)(CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(claswp)(CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zlaswp)(CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(strexc)(char *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dtrexc)(char *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(ctrexc)(char *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ztrexc)(char *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(stgexc)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dtgexc)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ctgexc)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ztgexc)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(strsyl)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dtrsyl)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(ctrsyl)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(ztrsyl)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(strsen)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dtrsen)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ctrsen)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ztrsen)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(stgsen)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dtgsen)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ctgsen)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, c64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ztgsen)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, c128 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sorghr)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dorghr)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cunghr)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zunghr)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sorgqr)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dorgqr)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cungqr)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zungqr)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sorgrq)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dorgrq)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cungrq)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zungrq)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sormqr)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dormqr)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cunmqr)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zunmqr)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sormrz)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dormrz)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cunmrz)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zunmrz)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgeqrt)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dgeqrt)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cgeqrt)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zgeqrt)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(sgemqrt)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dgemqrt)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(cgemqrt)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(zgemqrt)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(stpqrt)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dtpqrt)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(ctpqrt)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(ztpqrt)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(stpmqrt)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dtpmqrt)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *);
        void BLAS_FUNC(ctpmqrt)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *);
        void BLAS_FUNC(ztpmqrt)(char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *);

        void BLAS_FUNC(stzrzf)(CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dtzrzf)(CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ctzrzf)(CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ztzrzf)(CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sorcsd)(char *, char *, char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dorcsd)(char *, char *, char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cuncsd)(char *, char *, char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zuncsd)(char *, char *, char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgejsv)(char *, char *, char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgejsv)(char *, char *, char *, char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(sgglse)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dgglse)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(cgglse)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, c64 *, c64 *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zgglse)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, c128 *, c128 *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(slasd4)(CBLAS_INT *, CBLAS_INT *, f32 *, f32 *, f32 *, f32 *, f32 *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dlasd4)(CBLAS_INT *, CBLAS_INT *, f64 *, f64 *, f64 *, f64 *, f64 *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(ssbev)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *);
        void BLAS_FUNC(dsbev)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *);

        void BLAS_FUNC(ssbevd)(char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsbevd)(char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(chbevd)(char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zhbevd)(char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        void BLAS_FUNC(ssbevx)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dsbevx)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(chbevx)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, c64 *, CBLAS_INT *, c64 *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(zhbevx)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, c128 *, CBLAS_INT *, c128 *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);

        f32 BLAS_FUNC(slamch)(char *);
        f64 BLAS_FUNC(dlamch)(char *);

        f32 BLAS_FUNC(slange)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *);
        f64 BLAS_FUNC(dlange)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *);
        f32 BLAS_FUNC(clange)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *);
        f64 BLAS_FUNC(zlange)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *);

        f32 BLAS_FUNC(slantr)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *);
        f64 BLAS_FUNC(dlantr)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *);
        f32 BLAS_FUNC(clantr)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *);
        f64 BLAS_FUNC(zlantr)(char *, char *, char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *);

        void BLAS_FUNC(slarfg)(CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *);
        void BLAS_FUNC(dlarfg)(CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *);
        void BLAS_FUNC(clarfg)(CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, c64 *);
        void BLAS_FUNC(zlarfg)(CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *);

        void BLAS_FUNC(slarf)(char *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, CBLAS_INT *, f32 *);
        void BLAS_FUNC(dlarf)(char *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, CBLAS_INT *, f64 *);
        void BLAS_FUNC(clarf)(char *, CBLAS_INT *, CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, c64 *, CBLAS_INT *, c64 *);
        void BLAS_FUNC(zlarf)(char *, CBLAS_INT *, CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, c128 *, CBLAS_INT *, c128 *);

        void BLAS_FUNC(slartg)(f32 *, f32 *, f32 *, f32 *, f32 *);
        void BLAS_FUNC(dlartg)(f64 *, f64 *, f64 *, f64 *, f64 *);
        void BLAS_FUNC(clartg)(c64 *, c64 *, f32 *, c64 *, c64 *);
        void BLAS_FUNC(zlartg)(c128 *, c128 *, f64 *, c128 *, c128 *);

        void BLAS_FUNC(crot)(CBLAS_INT *, c64 *, CBLAS_INT *, c64 *, CBLAS_INT *, f32 *, c64 *);
        void BLAS_FUNC(zrot)(CBLAS_INT *, c128 *, CBLAS_INT *, c128 *, CBLAS_INT *, f64 *, c128 *);


        void BLAS_FUNC(stgsyl)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, CBLAS_INT *, f32 *, f32 *, f32 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(dtgsyl)(char *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, CBLAS_INT *, f64 *, f64 *, f64 *, CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
        void BLAS_FUNC(ilaver)(CBLAS_INT *, CBLAS_INT *, CBLAS_INT *);
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
    inline void syevd(char jobz, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *w, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(cheevd)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, info); }
    inline void syevd(char jobz, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *w, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(zheevd)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, info); }

    inline void syevr(char jobz, char range, char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 vl, f32 vu, CBLAS_INT il, CBLAS_INT iu, f32 abstol, CBLAS_INT *m, f32 *w, f32 *z, CBLAS_INT ldz, CBLAS_INT *isuppz, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(ssyevr)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, info); }
    inline void syevr(char jobz, char range, char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 vl, f64 vu, CBLAS_INT il, CBLAS_INT iu, f64 abstol, CBLAS_INT *m, f64 *w, f64 *z, CBLAS_INT ldz, CBLAS_INT *isuppz, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(dsyevr)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, info); }
    inline void syevr(char jobz, char range, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 vl, f32 vu, CBLAS_INT il, CBLAS_INT iu, f32 abstol, CBLAS_INT *m, f32 *w, c64 *z, CBLAS_INT ldz, CBLAS_INT *isuppz, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(cheevr)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, isuppz, work, &lwork, rwork, &lrwork, iwork, &liwork, info); }
    inline void syevr(char jobz, char range, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 vl, f64 vu, CBLAS_INT il, CBLAS_INT iu, f64 abstol, CBLAS_INT *m, f64 *w, c128 *z, CBLAS_INT ldz, CBLAS_INT *isuppz, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
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
    inline void sygvd(CBLAS_INT itype, char jobz, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, f32 *w, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(chegvd)(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &lrwork, iwork, &liwork, info); }
    inline void sygvd(CBLAS_INT itype, char jobz, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, f64 *w, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
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

    /* The triangular storage conversions.  Each moves one triangle between the three layouts --
     * packed (`ap`, TP), rectangular full packed (`arf`, TF) and full (`a`, TR) -- so the pair of
     * letters in the name reads as "from, to", and only the destination is written. */
    inline void tpttf(char transr, char uplo, CBLAS_INT n, f32 *ap, f32 *arf, CBLAS_INT *info)
        { BLAS_FUNC(stpttf)(&transr, &uplo, &n, ap, arf, info); }
    inline void tpttf(char transr, char uplo, CBLAS_INT n, f64 *ap, f64 *arf, CBLAS_INT *info)
        { BLAS_FUNC(dtpttf)(&transr, &uplo, &n, ap, arf, info); }
    inline void tpttf(char transr, char uplo, CBLAS_INT n, c64 *ap, c64 *arf, CBLAS_INT *info)
        { BLAS_FUNC(ctpttf)(&transr, &uplo, &n, ap, arf, info); }
    inline void tpttf(char transr, char uplo, CBLAS_INT n, c128 *ap, c128 *arf, CBLAS_INT *info)
        { BLAS_FUNC(ztpttf)(&transr, &uplo, &n, ap, arf, info); }

    inline void tpttr(char uplo, CBLAS_INT n, f32 *ap, f32 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(stpttr)(&uplo, &n, ap, a, &lda, info); }
    inline void tpttr(char uplo, CBLAS_INT n, f64 *ap, f64 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(dtpttr)(&uplo, &n, ap, a, &lda, info); }
    inline void tpttr(char uplo, CBLAS_INT n, c64 *ap, c64 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(ctpttr)(&uplo, &n, ap, a, &lda, info); }
    inline void tpttr(char uplo, CBLAS_INT n, c128 *ap, c128 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(ztpttr)(&uplo, &n, ap, a, &lda, info); }

    inline void tfttp(char transr, char uplo, CBLAS_INT n, f32 *arf, f32 *ap, CBLAS_INT *info)
        { BLAS_FUNC(stfttp)(&transr, &uplo, &n, arf, ap, info); }
    inline void tfttp(char transr, char uplo, CBLAS_INT n, f64 *arf, f64 *ap, CBLAS_INT *info)
        { BLAS_FUNC(dtfttp)(&transr, &uplo, &n, arf, ap, info); }
    inline void tfttp(char transr, char uplo, CBLAS_INT n, c64 *arf, c64 *ap, CBLAS_INT *info)
        { BLAS_FUNC(ctfttp)(&transr, &uplo, &n, arf, ap, info); }
    inline void tfttp(char transr, char uplo, CBLAS_INT n, c128 *arf, c128 *ap, CBLAS_INT *info)
        { BLAS_FUNC(ztfttp)(&transr, &uplo, &n, arf, ap, info); }

    inline void tfttr(char transr, char uplo, CBLAS_INT n, f32 *arf, f32 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(stfttr)(&transr, &uplo, &n, arf, a, &lda, info); }
    inline void tfttr(char transr, char uplo, CBLAS_INT n, f64 *arf, f64 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(dtfttr)(&transr, &uplo, &n, arf, a, &lda, info); }
    inline void tfttr(char transr, char uplo, CBLAS_INT n, c64 *arf, c64 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(ctfttr)(&transr, &uplo, &n, arf, a, &lda, info); }
    inline void tfttr(char transr, char uplo, CBLAS_INT n, c128 *arf, c128 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(ztfttr)(&transr, &uplo, &n, arf, a, &lda, info); }

    inline void trttf(char transr, char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *arf, CBLAS_INT *info)
        { BLAS_FUNC(strttf)(&transr, &uplo, &n, a, &lda, arf, info); }
    inline void trttf(char transr, char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *arf, CBLAS_INT *info)
        { BLAS_FUNC(dtrttf)(&transr, &uplo, &n, a, &lda, arf, info); }
    inline void trttf(char transr, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *arf, CBLAS_INT *info)
        { BLAS_FUNC(ctrttf)(&transr, &uplo, &n, a, &lda, arf, info); }
    inline void trttf(char transr, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *arf, CBLAS_INT *info)
        { BLAS_FUNC(ztrttf)(&transr, &uplo, &n, a, &lda, arf, info); }

    inline void trttp(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *ap, CBLAS_INT *info)
        { BLAS_FUNC(strttp)(&uplo, &n, a, &lda, ap, info); }
    inline void trttp(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *ap, CBLAS_INT *info)
        { BLAS_FUNC(dtrttp)(&uplo, &n, a, &lda, ap, info); }
    inline void trttp(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *ap, CBLAS_INT *info)
        { BLAS_FUNC(ctrttp)(&uplo, &n, a, &lda, ap, info); }
    inline void trttp(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *ap, CBLAS_INT *info)
        { BLAS_FUNC(ztrttp)(&uplo, &n, a, &lda, ap, info); }

    /* `tfsm` is `trsm` for an RFP triangle: no `info`, and `alpha` is passed by value like the
     * BLAS scalars it mirrors. */
    inline void tfsm(char transr, char side, char uplo, char trans, char diag, CBLAS_INT m, CBLAS_INT n, f32 alpha, f32 *a, f32 *b, CBLAS_INT ldb)
        { BLAS_FUNC(stfsm)(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb); }
    inline void tfsm(char transr, char side, char uplo, char trans, char diag, CBLAS_INT m, CBLAS_INT n, f64 alpha, f64 *a, f64 *b, CBLAS_INT ldb)
        { BLAS_FUNC(dtfsm)(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb); }
    inline void tfsm(char transr, char side, char uplo, char trans, char diag, CBLAS_INT m, CBLAS_INT n, c64 alpha, c64 *a, c64 *b, CBLAS_INT ldb)
        { BLAS_FUNC(ctfsm)(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb); }
    inline void tfsm(char transr, char side, char uplo, char trans, char diag, CBLAS_INT m, CBLAS_INT n, c128 alpha, c128 *a, c128 *b, CBLAS_INT ldb)
        { BLAS_FUNC(ztfsm)(&transr, &side, &uplo, &trans, &diag, &m, &n, &alpha, a, b, &ldb); }

    /* The packed positive definite family: the same Cholesky routines as `po*`, over a matrix
     * held as the `n * (n + 1) / 2` elements of one triangle. */
    inline void ppcon(char uplo, CBLAS_INT n, f32 *ap, f32 anorm, f32 *rcond, f32 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(sppcon)(&uplo, &n, ap, &anorm, rcond, work, iwork, info); }
    inline void ppcon(char uplo, CBLAS_INT n, f64 *ap, f64 anorm, f64 *rcond, f64 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dppcon)(&uplo, &n, ap, &anorm, rcond, work, iwork, info); }
    inline void ppcon(char uplo, CBLAS_INT n, c64 *ap, f32 anorm, f32 *rcond, c64 *work, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(cppcon)(&uplo, &n, ap, &anorm, rcond, work, rwork, info); }
    inline void ppcon(char uplo, CBLAS_INT n, c128 *ap, f64 anorm, f64 *rcond, c128 *work, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(zppcon)(&uplo, &n, ap, &anorm, rcond, work, rwork, info); }

    inline void ppsv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f32 *ap, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(sppsv)(&uplo, &n, &nrhs, ap, b, &ldb, info); }
    inline void ppsv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f64 *ap, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dppsv)(&uplo, &n, &nrhs, ap, b, &ldb, info); }
    inline void ppsv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c64 *ap, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cppsv)(&uplo, &n, &nrhs, ap, b, &ldb, info); }
    inline void ppsv(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c128 *ap, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zppsv)(&uplo, &n, &nrhs, ap, b, &ldb, info); }

    inline void pptrf(char uplo, CBLAS_INT n, f32 *ap, CBLAS_INT *info)
        { BLAS_FUNC(spptrf)(&uplo, &n, ap, info); }
    inline void pptrf(char uplo, CBLAS_INT n, f64 *ap, CBLAS_INT *info)
        { BLAS_FUNC(dpptrf)(&uplo, &n, ap, info); }
    inline void pptrf(char uplo, CBLAS_INT n, c64 *ap, CBLAS_INT *info)
        { BLAS_FUNC(cpptrf)(&uplo, &n, ap, info); }
    inline void pptrf(char uplo, CBLAS_INT n, c128 *ap, CBLAS_INT *info)
        { BLAS_FUNC(zpptrf)(&uplo, &n, ap, info); }

    inline void pptri(char uplo, CBLAS_INT n, f32 *ap, CBLAS_INT *info)
        { BLAS_FUNC(spptri)(&uplo, &n, ap, info); }
    inline void pptri(char uplo, CBLAS_INT n, f64 *ap, CBLAS_INT *info)
        { BLAS_FUNC(dpptri)(&uplo, &n, ap, info); }
    inline void pptri(char uplo, CBLAS_INT n, c64 *ap, CBLAS_INT *info)
        { BLAS_FUNC(cpptri)(&uplo, &n, ap, info); }
    inline void pptri(char uplo, CBLAS_INT n, c128 *ap, CBLAS_INT *info)
        { BLAS_FUNC(zpptri)(&uplo, &n, ap, info); }

    inline void pptrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f32 *ap, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(spptrs)(&uplo, &n, &nrhs, ap, b, &ldb, info); }
    inline void pptrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, f64 *ap, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dpptrs)(&uplo, &n, &nrhs, ap, b, &ldb, info); }
    inline void pptrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c64 *ap, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cpptrs)(&uplo, &n, &nrhs, ap, b, &ldb, info); }
    inline void pptrs(char uplo, CBLAS_INT n, CBLAS_INT nrhs, c128 *ap, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zpptrs)(&uplo, &n, &nrhs, ap, b, &ldb, info); }

    /* `sfrk`/`hfrk` is `syrk`/`herk` writing its result into an RFP triangle.  `alpha` and
     * `beta` are real in every flavor, as they are for `herk`. */
    inline void sfrk(char transr, char uplo, char trans, CBLAS_INT n, CBLAS_INT k, f32 alpha, f32 *a, CBLAS_INT lda, f32 beta, f32 *c)
        { BLAS_FUNC(ssfrk)(&transr, &uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c); }
    inline void sfrk(char transr, char uplo, char trans, CBLAS_INT n, CBLAS_INT k, f64 alpha, f64 *a, CBLAS_INT lda, f64 beta, f64 *c)
        { BLAS_FUNC(dsfrk)(&transr, &uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c); }
    inline void sfrk(char transr, char uplo, char trans, CBLAS_INT n, CBLAS_INT k, f32 alpha, c64 *a, CBLAS_INT lda, f32 beta, c64 *c)
        { BLAS_FUNC(chfrk)(&transr, &uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c); }
    inline void sfrk(char transr, char uplo, char trans, CBLAS_INT n, CBLAS_INT k, f64 alpha, c128 *a, CBLAS_INT lda, f64 beta, c128 *c)
        { BLAS_FUNC(zhfrk)(&transr, &uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c); }

    /* The rectangular-full-packed Cholesky trio: `po*` over an RFP triangle, so `transr` and
     * `uplo` name the layout and there is no leading dimension to pass. */
    inline void pftrf(char transr, char uplo, CBLAS_INT n, f32 *a, CBLAS_INT *info)
        { BLAS_FUNC(spftrf)(&transr, &uplo, &n, a, info); }
    inline void pftrf(char transr, char uplo, CBLAS_INT n, f64 *a, CBLAS_INT *info)
        { BLAS_FUNC(dpftrf)(&transr, &uplo, &n, a, info); }
    inline void pftrf(char transr, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT *info)
        { BLAS_FUNC(cpftrf)(&transr, &uplo, &n, a, info); }
    inline void pftrf(char transr, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT *info)
        { BLAS_FUNC(zpftrf)(&transr, &uplo, &n, a, info); }

    inline void pftri(char transr, char uplo, CBLAS_INT n, f32 *a, CBLAS_INT *info)
        { BLAS_FUNC(spftri)(&transr, &uplo, &n, a, info); }
    inline void pftri(char transr, char uplo, CBLAS_INT n, f64 *a, CBLAS_INT *info)
        { BLAS_FUNC(dpftri)(&transr, &uplo, &n, a, info); }
    inline void pftri(char transr, char uplo, CBLAS_INT n, c64 *a, CBLAS_INT *info)
        { BLAS_FUNC(cpftri)(&transr, &uplo, &n, a, info); }
    inline void pftri(char transr, char uplo, CBLAS_INT n, c128 *a, CBLAS_INT *info)
        { BLAS_FUNC(zpftri)(&transr, &uplo, &n, a, info); }

    inline void pftrs(char transr, char uplo, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(spftrs)(&transr, &uplo, &n, &nrhs, a, b, &ldb, info); }
    inline void pftrs(char transr, char uplo, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dpftrs)(&transr, &uplo, &n, &nrhs, a, b, &ldb, info); }
    inline void pftrs(char transr, char uplo, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cpftrs)(&transr, &uplo, &n, &nrhs, a, b, &ldb, info); }
    inline void pftrs(char transr, char uplo, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zpftrs)(&transr, &uplo, &n, &nrhs, a, b, &ldb, info); }

    /* The banded Cholesky trio.  `kd` is the number of super- or subdiagonals, one less than
     * the row count of the band array. */
    inline void pbtrf(char uplo, CBLAS_INT n, CBLAS_INT kd, f32 *ab, CBLAS_INT ldab, CBLAS_INT *info)
        { BLAS_FUNC(spbtrf)(&uplo, &n, &kd, ab, &ldab, info); }
    inline void pbtrf(char uplo, CBLAS_INT n, CBLAS_INT kd, f64 *ab, CBLAS_INT ldab, CBLAS_INT *info)
        { BLAS_FUNC(dpbtrf)(&uplo, &n, &kd, ab, &ldab, info); }
    inline void pbtrf(char uplo, CBLAS_INT n, CBLAS_INT kd, c64 *ab, CBLAS_INT ldab, CBLAS_INT *info)
        { BLAS_FUNC(cpbtrf)(&uplo, &n, &kd, ab, &ldab, info); }
    inline void pbtrf(char uplo, CBLAS_INT n, CBLAS_INT kd, c128 *ab, CBLAS_INT ldab, CBLAS_INT *info)
        { BLAS_FUNC(zpbtrf)(&uplo, &n, &kd, ab, &ldab, info); }

    inline void pbtrs(char uplo, CBLAS_INT n, CBLAS_INT kd, CBLAS_INT nrhs, f32 *ab, CBLAS_INT ldab, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(spbtrs)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, info); }
    inline void pbtrs(char uplo, CBLAS_INT n, CBLAS_INT kd, CBLAS_INT nrhs, f64 *ab, CBLAS_INT ldab, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dpbtrs)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, info); }
    inline void pbtrs(char uplo, CBLAS_INT n, CBLAS_INT kd, CBLAS_INT nrhs, c64 *ab, CBLAS_INT ldab, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cpbtrs)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, info); }
    inline void pbtrs(char uplo, CBLAS_INT n, CBLAS_INT kd, CBLAS_INT nrhs, c128 *ab, CBLAS_INT ldab, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zpbtrs)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, info); }

    inline void pbsv(char uplo, CBLAS_INT n, CBLAS_INT kd, CBLAS_INT nrhs, f32 *ab, CBLAS_INT ldab, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(spbsv)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, info); }
    inline void pbsv(char uplo, CBLAS_INT n, CBLAS_INT kd, CBLAS_INT nrhs, f64 *ab, CBLAS_INT ldab, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dpbsv)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, info); }
    inline void pbsv(char uplo, CBLAS_INT n, CBLAS_INT kd, CBLAS_INT nrhs, c64 *ab, CBLAS_INT ldab, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(cpbsv)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, info); }
    inline void pbsv(char uplo, CBLAS_INT n, CBLAS_INT kd, CBLAS_INT nrhs, c128 *ab, CBLAS_INT ldab, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(zpbsv)(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, info); }

    /* Triangular solves, inverse and condition number, plus the two LU auxiliaries. */
    inline void trtrs(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT nrhs, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(strtrs)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, info); }
    inline void trtrs(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT nrhs, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dtrtrs)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, info); }
    inline void trtrs(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT nrhs, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(ctrtrs)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, info); }
    inline void trtrs(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT nrhs, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(ztrtrs)(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, info); }

    inline void trcon(char norm, char uplo, char diag, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *rcond, f32 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(strcon)(&norm, &uplo, &diag, &n, a, &lda, rcond, work, iwork, info); }
    inline void trcon(char norm, char uplo, char diag, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *rcond, f64 *work, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dtrcon)(&norm, &uplo, &diag, &n, a, &lda, rcond, work, iwork, info); }
    inline void trcon(char norm, char uplo, char diag, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *rcond, c64 *work, f32 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(ctrcon)(&norm, &uplo, &diag, &n, a, &lda, rcond, work, rwork, info); }
    inline void trcon(char norm, char uplo, char diag, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *rcond, c128 *work, f64 *rwork, CBLAS_INT *info)
        { BLAS_FUNC(ztrcon)(&norm, &uplo, &diag, &n, a, &lda, rcond, work, rwork, info); }

    inline void tbtrs(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT kd, CBLAS_INT nrhs, f32 *ab, CBLAS_INT ldab, f32 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(stbtrs)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, info); }
    inline void tbtrs(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT kd, CBLAS_INT nrhs, f64 *ab, CBLAS_INT ldab, f64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(dtbtrs)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, info); }
    inline void tbtrs(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT kd, CBLAS_INT nrhs, c64 *ab, CBLAS_INT ldab, c64 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(ctbtrs)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, info); }
    inline void tbtrs(char uplo, char trans, char diag, CBLAS_INT n, CBLAS_INT kd, CBLAS_INT nrhs, c128 *ab, CBLAS_INT ldab, c128 *b, CBLAS_INT ldb, CBLAS_INT *info)
        { BLAS_FUNC(ztbtrs)(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, info); }

    inline void trtri(char uplo, char diag, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(strtri)(&uplo, &diag, &n, a, &lda, info); }
    inline void trtri(char uplo, char diag, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(dtrtri)(&uplo, &diag, &n, a, &lda, info); }
    inline void trtri(char uplo, char diag, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(ctrtri)(&uplo, &diag, &n, a, &lda, info); }
    inline void trtri(char uplo, char diag, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(ztrtri)(&uplo, &diag, &n, a, &lda, info); }

    inline void lauum(char uplo, CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(slauum)(&uplo, &n, a, &lda, info); }
    inline void lauum(char uplo, CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(dlauum)(&uplo, &n, a, &lda, info); }
    inline void lauum(char uplo, CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(clauum)(&uplo, &n, a, &lda, info); }
    inline void lauum(char uplo, CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT *info)
        { BLAS_FUNC(zlauum)(&uplo, &n, a, &lda, info); }

    /* `laswp` numbers its pivots from 1 and its row range `k1`/`k2` from 1 as well. */
    inline void laswp(CBLAS_INT n, f32 *a, CBLAS_INT lda, CBLAS_INT k1, CBLAS_INT k2, CBLAS_INT *ipiv, CBLAS_INT incx)
        { BLAS_FUNC(slaswp)(&n, a, &lda, &k1, &k2, ipiv, &incx); }
    inline void laswp(CBLAS_INT n, f64 *a, CBLAS_INT lda, CBLAS_INT k1, CBLAS_INT k2, CBLAS_INT *ipiv, CBLAS_INT incx)
        { BLAS_FUNC(dlaswp)(&n, a, &lda, &k1, &k2, ipiv, &incx); }
    inline void laswp(CBLAS_INT n, c64 *a, CBLAS_INT lda, CBLAS_INT k1, CBLAS_INT k2, CBLAS_INT *ipiv, CBLAS_INT incx)
        { BLAS_FUNC(claswp)(&n, a, &lda, &k1, &k2, ipiv, &incx); }
    inline void laswp(CBLAS_INT n, c128 *a, CBLAS_INT lda, CBLAS_INT k1, CBLAS_INT k2, CBLAS_INT *ipiv, CBLAS_INT incx)
        { BLAS_FUNC(zlaswp)(&n, a, &lda, &k1, &k2, ipiv, &incx); }

    /* Schur reordering.  The real routines carry a workspace the complex ones have no use for,
     * so the two flavors take different argument lists rather than a `nullptr` placeholder --
     * overload resolution picks the right one from the element type alone.
     *
     * `tgexc` takes its `wantq`/`wantz` as Fortran LOGICALs, hence `CBLAS_INT` rather than the
     * option letter `trexc` uses for the same choice. */
    inline void trexc(char compq, CBLAS_INT n, f32 *t, CBLAS_INT ldt, f32 *q, CBLAS_INT ldq, CBLAS_INT *ifst, CBLAS_INT *ilst, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(strexc)(&compq, &n, t, &ldt, q, &ldq, ifst, ilst, work, info); }
    inline void trexc(char compq, CBLAS_INT n, f64 *t, CBLAS_INT ldt, f64 *q, CBLAS_INT ldq, CBLAS_INT *ifst, CBLAS_INT *ilst, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dtrexc)(&compq, &n, t, &ldt, q, &ldq, ifst, ilst, work, info); }
    inline void trexc(char compq, CBLAS_INT n, c64 *t, CBLAS_INT ldt, c64 *q, CBLAS_INT ldq, CBLAS_INT *ifst, CBLAS_INT *ilst, CBLAS_INT *info)
        { BLAS_FUNC(ctrexc)(&compq, &n, t, &ldt, q, &ldq, ifst, ilst, info); }
    inline void trexc(char compq, CBLAS_INT n, c128 *t, CBLAS_INT ldt, c128 *q, CBLAS_INT ldq, CBLAS_INT *ifst, CBLAS_INT *ilst, CBLAS_INT *info)
        { BLAS_FUNC(ztrexc)(&compq, &n, t, &ldt, q, &ldq, ifst, ilst, info); }

    inline void tgexc(CBLAS_INT wantq, CBLAS_INT wantz, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *q, CBLAS_INT ldq, f32 *z, CBLAS_INT ldz, CBLAS_INT *ifst, CBLAS_INT *ilst, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(stgexc)(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, ifst, ilst, work, &lwork, info); }
    inline void tgexc(CBLAS_INT wantq, CBLAS_INT wantz, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *q, CBLAS_INT ldq, f64 *z, CBLAS_INT ldz, CBLAS_INT *ifst, CBLAS_INT *ilst, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dtgexc)(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, ifst, ilst, work, &lwork, info); }
    inline void tgexc(CBLAS_INT wantq, CBLAS_INT wantz, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, c64 *q, CBLAS_INT ldq, c64 *z, CBLAS_INT ldz, CBLAS_INT *ifst, CBLAS_INT *ilst, CBLAS_INT *info)
        { BLAS_FUNC(ctgexc)(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, ifst, ilst, info); }
    inline void tgexc(CBLAS_INT wantq, CBLAS_INT wantz, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, c128 *q, CBLAS_INT ldq, c128 *z, CBLAS_INT ldz, CBLAS_INT *ifst, CBLAS_INT *ilst, CBLAS_INT *info)
        { BLAS_FUNC(ztgexc)(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, ifst, ilst, info); }

    /* `trsyl`'s `scale` is real in every flavor. */
    inline void trsyl(char trana, char tranb, CBLAS_INT isgn, CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *c, CBLAS_INT ldc, f32 *scale, CBLAS_INT *info)
        { BLAS_FUNC(strsyl)(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, info); }
    inline void trsyl(char trana, char tranb, CBLAS_INT isgn, CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *c, CBLAS_INT ldc, f64 *scale, CBLAS_INT *info)
        { BLAS_FUNC(dtrsyl)(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, info); }
    inline void trsyl(char trana, char tranb, CBLAS_INT isgn, CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, c64 *c, CBLAS_INT ldc, f32 *scale, CBLAS_INT *info)
        { BLAS_FUNC(ctrsyl)(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, info); }
    inline void trsyl(char trana, char tranb, CBLAS_INT isgn, CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, c128 *c, CBLAS_INT ldc, f64 *scale, CBLAS_INT *info)
        { BLAS_FUNC(ztrsyl)(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, info); }

    /* Schur cluster reordering.  The real and complex halves take genuinely different argument
     * lists, not just different types: the real `trsen` splits its eigenvalues into `wr`/`wi`
     * and carries an integer workspace, the complex one returns a single `w` and has neither.
     * `tgsen` splits the same way over `alphar`/`alphai` against `alpha`.
     *
     * `select`, `wantq` and `wantz` are Fortran LOGICALs, so `CBLAS_INT` throughout. */
    inline void trsen(char job, char compq, CBLAS_INT *select, CBLAS_INT n, f32 *t, CBLAS_INT ldt, f32 *q, CBLAS_INT ldq, f32 *wr, f32 *wi, CBLAS_INT *m, f32 *s, f32 *sep, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(strsen)(&job, &compq, select, &n, t, &ldt, q, &ldq, wr, wi, m, s, sep, work, &lwork, iwork, &liwork, info); }
    inline void trsen(char job, char compq, CBLAS_INT *select, CBLAS_INT n, f64 *t, CBLAS_INT ldt, f64 *q, CBLAS_INT ldq, f64 *wr, f64 *wi, CBLAS_INT *m, f64 *s, f64 *sep, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(dtrsen)(&job, &compq, select, &n, t, &ldt, q, &ldq, wr, wi, m, s, sep, work, &lwork, iwork, &liwork, info); }
    inline void trsen(char job, char compq, CBLAS_INT *select, CBLAS_INT n, c64 *t, CBLAS_INT ldt, c64 *q, CBLAS_INT ldq, c64 *w, CBLAS_INT *m, f32 *s, f32 *sep, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(ctrsen)(&job, &compq, select, &n, t, &ldt, q, &ldq, w, m, s, sep, work, &lwork, info); }
    inline void trsen(char job, char compq, CBLAS_INT *select, CBLAS_INT n, c128 *t, CBLAS_INT ldt, c128 *q, CBLAS_INT ldq, c128 *w, CBLAS_INT *m, f64 *s, f64 *sep, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(ztrsen)(&job, &compq, select, &n, t, &ldt, q, &ldq, w, m, s, sep, work, &lwork, info); }

    inline void tgsen(CBLAS_INT ijob, CBLAS_INT wantq, CBLAS_INT wantz, CBLAS_INT *select, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *alphar, f32 *alphai, f32 *beta, f32 *q, CBLAS_INT ldq, f32 *z, CBLAS_INT ldz, CBLAS_INT *m, f32 *pl, f32 *pr, f32 *dif, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(stgsen)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alphar, alphai, beta, q, &ldq, z, &ldz, m, pl, pr, dif, work, &lwork, iwork, &liwork, info); }
    inline void tgsen(CBLAS_INT ijob, CBLAS_INT wantq, CBLAS_INT wantz, CBLAS_INT *select, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *alphar, f64 *alphai, f64 *beta, f64 *q, CBLAS_INT ldq, f64 *z, CBLAS_INT ldz, CBLAS_INT *m, f64 *pl, f64 *pr, f64 *dif, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(dtgsen)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alphar, alphai, beta, q, &ldq, z, &ldz, m, pl, pr, dif, work, &lwork, iwork, &liwork, info); }
    inline void tgsen(CBLAS_INT ijob, CBLAS_INT wantq, CBLAS_INT wantz, CBLAS_INT *select, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, c64 *alpha, c64 *beta, c64 *q, CBLAS_INT ldq, c64 *z, CBLAS_INT ldz, CBLAS_INT *m, f32 *pl, f32 *pr, f32 *dif, c64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(ctgsen)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alpha, beta, q, &ldq, z, &ldz, m, pl, pr, dif, work, &lwork, iwork, &liwork, info); }
    inline void tgsen(CBLAS_INT ijob, CBLAS_INT wantq, CBLAS_INT wantz, CBLAS_INT *select, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, c128 *alpha, c128 *beta, c128 *q, CBLAS_INT ldq, c128 *z, CBLAS_INT ldz, CBLAS_INT *m, f64 *pl, f64 *pr, f64 *dif, c128 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(ztgsen)(&ijob, &wantq, &wantz, select, &n, a, &lda, b, &ldb, alpha, beta, q, &ldq, z, &ldz, m, pl, pr, dif, work, &lwork, iwork, &liwork, info); }

    /* Generating and applying the orthogonal/unitary factor.  LAPACK spells the same routine
     * `or*` for real and `un*` for complex, so these overload sets carry the real spelling and
     * the wrapper's Ctx picks the Python name. */
    inline void orghr(CBLAS_INT n, CBLAS_INT ilo, CBLAS_INT ihi, f32 *a, CBLAS_INT lda, f32 *tau, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sorghr)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, info); }
    inline void orghr(CBLAS_INT n, CBLAS_INT ilo, CBLAS_INT ihi, f64 *a, CBLAS_INT lda, f64 *tau, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dorghr)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, info); }
    inline void orghr(CBLAS_INT n, CBLAS_INT ilo, CBLAS_INT ihi, c64 *a, CBLAS_INT lda, c64 *tau, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(cunghr)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, info); }
    inline void orghr(CBLAS_INT n, CBLAS_INT ilo, CBLAS_INT ihi, c128 *a, CBLAS_INT lda, c128 *tau, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zunghr)(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, info); }

    inline void orgqr(CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, f32 *a, CBLAS_INT lda, f32 *tau, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sorgqr)(&m, &n, &k, a, &lda, tau, work, &lwork, info); }
    inline void orgqr(CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, f64 *a, CBLAS_INT lda, f64 *tau, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dorgqr)(&m, &n, &k, a, &lda, tau, work, &lwork, info); }
    inline void orgqr(CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, c64 *a, CBLAS_INT lda, c64 *tau, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(cungqr)(&m, &n, &k, a, &lda, tau, work, &lwork, info); }
    inline void orgqr(CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, c128 *a, CBLAS_INT lda, c128 *tau, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zungqr)(&m, &n, &k, a, &lda, tau, work, &lwork, info); }

    inline void orgrq(CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, f32 *a, CBLAS_INT lda, f32 *tau, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sorgrq)(&m, &n, &k, a, &lda, tau, work, &lwork, info); }
    inline void orgrq(CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, f64 *a, CBLAS_INT lda, f64 *tau, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dorgrq)(&m, &n, &k, a, &lda, tau, work, &lwork, info); }
    inline void orgrq(CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, c64 *a, CBLAS_INT lda, c64 *tau, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(cungrq)(&m, &n, &k, a, &lda, tau, work, &lwork, info); }
    inline void orgrq(CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, c128 *a, CBLAS_INT lda, c128 *tau, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zungrq)(&m, &n, &k, a, &lda, tau, work, &lwork, info); }

    inline void ormqr(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, f32 *a, CBLAS_INT lda, f32 *tau, f32 *c, CBLAS_INT ldc, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sormqr)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, info); }
    inline void ormqr(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, f64 *a, CBLAS_INT lda, f64 *tau, f64 *c, CBLAS_INT ldc, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dormqr)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, info); }
    inline void ormqr(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, c64 *a, CBLAS_INT lda, c64 *tau, c64 *c, CBLAS_INT ldc, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(cunmqr)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, info); }
    inline void ormqr(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, c128 *a, CBLAS_INT lda, c128 *tau, c128 *c, CBLAS_INT ldc, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zunmqr)(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, info); }

    inline void ormrz(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, CBLAS_INT l, f32 *a, CBLAS_INT lda, f32 *tau, f32 *c, CBLAS_INT ldc, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sormrz)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, info); }
    inline void ormrz(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, CBLAS_INT l, f64 *a, CBLAS_INT lda, f64 *tau, f64 *c, CBLAS_INT ldc, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dormrz)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, info); }
    inline void ormrz(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, CBLAS_INT l, c64 *a, CBLAS_INT lda, c64 *tau, c64 *c, CBLAS_INT ldc, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(cunmrz)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, info); }
    inline void ormrz(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, CBLAS_INT l, c128 *a, CBLAS_INT lda, c128 *tau, c128 *c, CBLAS_INT ldc, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zunmrz)(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, info); }

    /* Compact-WY and blocked QR.  These take a block size `nb` and hand back the triangular
     * block reflectors `t` instead of a `tau` vector, and their workspaces are sized from `nb`
     * rather than queried -- so none of them has an `lwork`. */
    inline void geqrt(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nb, f32 *a, CBLAS_INT lda, f32 *t, CBLAS_INT ldt, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(sgeqrt)(&m, &n, &nb, a, &lda, t, &ldt, work, info); }
    inline void geqrt(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nb, f64 *a, CBLAS_INT lda, f64 *t, CBLAS_INT ldt, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dgeqrt)(&m, &n, &nb, a, &lda, t, &ldt, work, info); }
    inline void geqrt(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nb, c64 *a, CBLAS_INT lda, c64 *t, CBLAS_INT ldt, c64 *work, CBLAS_INT *info)
        { BLAS_FUNC(cgeqrt)(&m, &n, &nb, a, &lda, t, &ldt, work, info); }
    inline void geqrt(CBLAS_INT m, CBLAS_INT n, CBLAS_INT nb, c128 *a, CBLAS_INT lda, c128 *t, CBLAS_INT ldt, c128 *work, CBLAS_INT *info)
        { BLAS_FUNC(zgeqrt)(&m, &n, &nb, a, &lda, t, &ldt, work, info); }

    inline void gemqrt(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, CBLAS_INT nb, f32 *v, CBLAS_INT ldv, f32 *t, CBLAS_INT ldt, f32 *c, CBLAS_INT ldc, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(sgemqrt)(&side, &trans, &m, &n, &k, &nb, v, &ldv, t, &ldt, c, &ldc, work, info); }
    inline void gemqrt(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, CBLAS_INT nb, f64 *v, CBLAS_INT ldv, f64 *t, CBLAS_INT ldt, f64 *c, CBLAS_INT ldc, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dgemqrt)(&side, &trans, &m, &n, &k, &nb, v, &ldv, t, &ldt, c, &ldc, work, info); }
    inline void gemqrt(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, CBLAS_INT nb, c64 *v, CBLAS_INT ldv, c64 *t, CBLAS_INT ldt, c64 *c, CBLAS_INT ldc, c64 *work, CBLAS_INT *info)
        { BLAS_FUNC(cgemqrt)(&side, &trans, &m, &n, &k, &nb, v, &ldv, t, &ldt, c, &ldc, work, info); }
    inline void gemqrt(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, CBLAS_INT nb, c128 *v, CBLAS_INT ldv, c128 *t, CBLAS_INT ldt, c128 *c, CBLAS_INT ldc, c128 *work, CBLAS_INT *info)
        { BLAS_FUNC(zgemqrt)(&side, &trans, &m, &n, &k, &nb, v, &ldv, t, &ldt, c, &ldc, work, info); }

    inline void tpqrt(CBLAS_INT m, CBLAS_INT n, CBLAS_INT l, CBLAS_INT nb, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *t, CBLAS_INT ldt, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(stpqrt)(&m, &n, &l, &nb, a, &lda, b, &ldb, t, &ldt, work, info); }
    inline void tpqrt(CBLAS_INT m, CBLAS_INT n, CBLAS_INT l, CBLAS_INT nb, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *t, CBLAS_INT ldt, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dtpqrt)(&m, &n, &l, &nb, a, &lda, b, &ldb, t, &ldt, work, info); }
    inline void tpqrt(CBLAS_INT m, CBLAS_INT n, CBLAS_INT l, CBLAS_INT nb, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, c64 *t, CBLAS_INT ldt, c64 *work, CBLAS_INT *info)
        { BLAS_FUNC(ctpqrt)(&m, &n, &l, &nb, a, &lda, b, &ldb, t, &ldt, work, info); }
    inline void tpqrt(CBLAS_INT m, CBLAS_INT n, CBLAS_INT l, CBLAS_INT nb, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, c128 *t, CBLAS_INT ldt, c128 *work, CBLAS_INT *info)
        { BLAS_FUNC(ztpqrt)(&m, &n, &l, &nb, a, &lda, b, &ldb, t, &ldt, work, info); }

    inline void tpmqrt(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, CBLAS_INT l, CBLAS_INT nb, f32 *v, CBLAS_INT ldv, f32 *t, CBLAS_INT ldt, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(stpmqrt)(&side, &trans, &m, &n, &k, &l, &nb, v, &ldv, t, &ldt, a, &lda, b, &ldb, work, info); }
    inline void tpmqrt(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, CBLAS_INT l, CBLAS_INT nb, f64 *v, CBLAS_INT ldv, f64 *t, CBLAS_INT ldt, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dtpmqrt)(&side, &trans, &m, &n, &k, &l, &nb, v, &ldv, t, &ldt, a, &lda, b, &ldb, work, info); }
    inline void tpmqrt(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, CBLAS_INT l, CBLAS_INT nb, c64 *v, CBLAS_INT ldv, c64 *t, CBLAS_INT ldt, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, c64 *work, CBLAS_INT *info)
        { BLAS_FUNC(ctpmqrt)(&side, &trans, &m, &n, &k, &l, &nb, v, &ldv, t, &ldt, a, &lda, b, &ldb, work, info); }
    inline void tpmqrt(char side, char trans, CBLAS_INT m, CBLAS_INT n, CBLAS_INT k, CBLAS_INT l, CBLAS_INT nb, c128 *v, CBLAS_INT ldv, c128 *t, CBLAS_INT ldt, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, c128 *work, CBLAS_INT *info)
        { BLAS_FUNC(ztpmqrt)(&side, &trans, &m, &n, &k, &l, &nb, v, &ldv, t, &ldt, a, &lda, b, &ldb, work, info); }

    inline void tzrzf(CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *tau, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(stzrzf)(&m, &n, a, &lda, tau, work, &lwork, info); }
    inline void tzrzf(CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *tau, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dtzrzf)(&m, &n, a, &lda, tau, work, &lwork, info); }
    inline void tzrzf(CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, c64 *tau, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(ctzrzf)(&m, &n, a, &lda, tau, work, &lwork, info); }
    inline void tzrzf(CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, c128 *tau, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(ztzrzf)(&m, &n, a, &lda, tau, work, &lwork, info); }

    /* The CS decomposition.  The complex routines carry a second, real workspace the real ones
     * have no use for, so the two take different argument lists; `theta` is real in both. */
    inline void orcsd(char jobu1, char jobu2, char jobv1t, char jobv2t, char trans, char signs, CBLAS_INT m, CBLAS_INT p, CBLAS_INT q, f32 *x11, CBLAS_INT ldx11, f32 *x12, CBLAS_INT ldx12, f32 *x21, CBLAS_INT ldx21, f32 *x22, CBLAS_INT ldx22, f32 *theta, f32 *u1, CBLAS_INT ldu1, f32 *u2, CBLAS_INT ldu2, f32 *v1t, CBLAS_INT ldv1t, f32 *v2t, CBLAS_INT ldv2t, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(sorcsd)(&jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q, x11, &ldx11, x12, &ldx12, x21, &ldx21, x22, &ldx22, theta, u1, &ldu1, u2, &ldu2, v1t, &ldv1t, v2t, &ldv2t, work, &lwork, iwork, info); }
    inline void orcsd(char jobu1, char jobu2, char jobv1t, char jobv2t, char trans, char signs, CBLAS_INT m, CBLAS_INT p, CBLAS_INT q, f64 *x11, CBLAS_INT ldx11, f64 *x12, CBLAS_INT ldx12, f64 *x21, CBLAS_INT ldx21, f64 *x22, CBLAS_INT ldx22, f64 *theta, f64 *u1, CBLAS_INT ldu1, f64 *u2, CBLAS_INT ldu2, f64 *v1t, CBLAS_INT ldv1t, f64 *v2t, CBLAS_INT ldv2t, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dorcsd)(&jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q, x11, &ldx11, x12, &ldx12, x21, &ldx21, x22, &ldx22, theta, u1, &ldu1, u2, &ldu2, v1t, &ldv1t, v2t, &ldv2t, work, &lwork, iwork, info); }
    inline void orcsd(char jobu1, char jobu2, char jobv1t, char jobv2t, char trans, char signs, CBLAS_INT m, CBLAS_INT p, CBLAS_INT q, c64 *x11, CBLAS_INT ldx11, c64 *x12, CBLAS_INT ldx12, c64 *x21, CBLAS_INT ldx21, c64 *x22, CBLAS_INT ldx22, f32 *theta, c64 *u1, CBLAS_INT ldu1, c64 *u2, CBLAS_INT ldu2, c64 *v1t, CBLAS_INT ldv1t, c64 *v2t, CBLAS_INT ldv2t, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(cuncsd)(&jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q, x11, &ldx11, x12, &ldx12, x21, &ldx21, x22, &ldx22, theta, u1, &ldu1, u2, &ldu2, v1t, &ldv1t, v2t, &ldv2t, work, &lwork, rwork, &lrwork, iwork, info); }
    inline void orcsd(char jobu1, char jobu2, char jobv1t, char jobv2t, char trans, char signs, CBLAS_INT m, CBLAS_INT p, CBLAS_INT q, c128 *x11, CBLAS_INT ldx11, c128 *x12, CBLAS_INT ldx12, c128 *x21, CBLAS_INT ldx21, c128 *x22, CBLAS_INT ldx22, f64 *theta, c128 *u1, CBLAS_INT ldu1, c128 *u2, CBLAS_INT ldu2, c128 *v1t, CBLAS_INT ldv1t, c128 *v2t, CBLAS_INT ldv2t, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(zuncsd)(&jobu1, &jobu2, &jobv1t, &jobv2t, &trans, &signs, &m, &p, &q, x11, &ldx11, x12, &ldx12, x21, &ldx21, x22, &ldx22, theta, u1, &ldu1, u2, &ldu2, v1t, &ldv1t, v2t, &ldv2t, work, &lwork, rwork, &lrwork, iwork, info); }

    /* `gejsv` and `lasd4` are real-only: LAPACK ships no complex counterpart for either. */
    inline void gejsv(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *sva, f32 *u, CBLAS_INT ldu, f32 *v, CBLAS_INT ldv, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(sgejsv)(&joba, &jobu, &jobv, &jobr, &jobt, &jobp, &m, &n, a, &lda, sva, u, &ldu, v, &ldv, work, &lwork, iwork, info); }
    inline void gejsv(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *sva, f64 *u, CBLAS_INT ldu, f64 *v, CBLAS_INT ldv, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dgejsv)(&joba, &jobu, &jobv, &jobr, &jobt, &jobp, &m, &n, a, &lda, sva, u, &ldu, v, &ldv, work, &lwork, iwork, info); }

    inline void gglse(CBLAS_INT m, CBLAS_INT n, CBLAS_INT p, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *c, f32 *d, f32 *x, f32 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(sgglse)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, info); }
    inline void gglse(CBLAS_INT m, CBLAS_INT n, CBLAS_INT p, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *c, f64 *d, f64 *x, f64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(dgglse)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, info); }
    inline void gglse(CBLAS_INT m, CBLAS_INT n, CBLAS_INT p, c64 *a, CBLAS_INT lda, c64 *b, CBLAS_INT ldb, c64 *c, c64 *d, c64 *x, c64 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(cgglse)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, info); }
    inline void gglse(CBLAS_INT m, CBLAS_INT n, CBLAS_INT p, c128 *a, CBLAS_INT lda, c128 *b, CBLAS_INT ldb, c128 *c, c128 *d, c128 *x, c128 *work, CBLAS_INT lwork, CBLAS_INT *info)
        { BLAS_FUNC(zgglse)(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, info); }

    inline void lasd4(CBLAS_INT n, CBLAS_INT i, f32 *d, f32 *z, f32 *delta, f32 rho, f32 *sigma, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(slasd4)(&n, &i, d, z, delta, &rho, sigma, work, info); }
    inline void lasd4(CBLAS_INT n, CBLAS_INT i, f64 *d, f64 *z, f64 *delta, f64 rho, f64 *sigma, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dlasd4)(&n, &i, d, z, delta, &rho, sigma, work, info); }

    /* The symmetric/Hermitian banded eigensolvers.  `sbev` is real-only; `sbevd`/`hbevd` and
     * `sbevx`/`hbevx` pair up, and there the complex flavors carry a real workspace the real
     * ones do not -- so, again, two argument lists rather than a placeholder. */
    inline void sbev(char jobz, char uplo, CBLAS_INT n, CBLAS_INT kd, f32 *ab, CBLAS_INT ldab, f32 *w, f32 *z, CBLAS_INT ldz, f32 *work, CBLAS_INT *info)
        { BLAS_FUNC(ssbev)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, info); }
    inline void sbev(char jobz, char uplo, CBLAS_INT n, CBLAS_INT kd, f64 *ab, CBLAS_INT ldab, f64 *w, f64 *z, CBLAS_INT ldz, f64 *work, CBLAS_INT *info)
        { BLAS_FUNC(dsbev)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, info); }

    inline void sbevd(char jobz, char uplo, CBLAS_INT n, CBLAS_INT kd, f32 *ab, CBLAS_INT ldab, f32 *w, f32 *z, CBLAS_INT ldz, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(ssbevd)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, iwork, &liwork, info); }
    inline void sbevd(char jobz, char uplo, CBLAS_INT n, CBLAS_INT kd, f64 *ab, CBLAS_INT ldab, f64 *w, f64 *z, CBLAS_INT ldz, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(dsbevd)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, iwork, &liwork, info); }
    inline void sbevd(char jobz, char uplo, CBLAS_INT n, CBLAS_INT kd, c64 *ab, CBLAS_INT ldab, f32 *w, c64 *z, CBLAS_INT ldz, c64 *work, CBLAS_INT lwork, f32 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(chbevd)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, info); }
    inline void sbevd(char jobz, char uplo, CBLAS_INT n, CBLAS_INT kd, c128 *ab, CBLAS_INT ldab, f64 *w, c128 *z, CBLAS_INT ldz, c128 *work, CBLAS_INT lwork, f64 *rwork, CBLAS_INT lrwork, CBLAS_INT *iwork, CBLAS_INT liwork, CBLAS_INT *info)
        { BLAS_FUNC(zhbevd)(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, info); }

    inline void sbevx(char jobz, char range, char uplo, CBLAS_INT n, CBLAS_INT kd, f32 *ab, CBLAS_INT ldab, f32 *q, CBLAS_INT ldq, f32 vl, f32 vu, CBLAS_INT il, CBLAS_INT iu, f32 abstol, CBLAS_INT *m, f32 *w, f32 *z, CBLAS_INT ldz, f32 *work, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(ssbevx)(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, info); }
    inline void sbevx(char jobz, char range, char uplo, CBLAS_INT n, CBLAS_INT kd, f64 *ab, CBLAS_INT ldab, f64 *q, CBLAS_INT ldq, f64 vl, f64 vu, CBLAS_INT il, CBLAS_INT iu, f64 abstol, CBLAS_INT *m, f64 *w, f64 *z, CBLAS_INT ldz, f64 *work, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(dsbevx)(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, iwork, ifail, info); }
    inline void sbevx(char jobz, char range, char uplo, CBLAS_INT n, CBLAS_INT kd, c64 *ab, CBLAS_INT ldab, c64 *q, CBLAS_INT ldq, f32 vl, f32 vu, CBLAS_INT il, CBLAS_INT iu, f32 abstol, CBLAS_INT *m, f32 *w, c64 *z, CBLAS_INT ldz, c64 *work, f32 *rwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(chbevx)(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, rwork, iwork, ifail, info); }
    inline void sbevx(char jobz, char range, char uplo, CBLAS_INT n, CBLAS_INT kd, c128 *ab, CBLAS_INT ldab, c128 *q, CBLAS_INT ldq, f64 vl, f64 vu, CBLAS_INT il, CBLAS_INT iu, f64 abstol, CBLAS_INT *m, f64 *w, c128 *z, CBLAS_INT ldz, c128 *work, f64 *rwork, CBLAS_INT *iwork, CBLAS_INT *ifail, CBLAS_INT *info)
        { BLAS_FUNC(zhbevx)(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, &vl, &vu, &il, &iu, &abstol, m, w, z, &ldz, work, rwork, iwork, ifail, info); }

    /* The auxiliaries.  These return values rather than filling an `info`, and `lamch` and
     * `ilaver` take no matrix at all -- `ilaver` takes nothing typed by a flavor, which is why
     * it is the one routine in this file with no template parameter. */
    inline f32 lamch(char cmach, f32 *) { return BLAS_FUNC(slamch)(&cmach); }
    inline f64 lamch(char cmach, f64 *) { return BLAS_FUNC(dlamch)(&cmach); }

    inline f32 lange(char norm, CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *work)
        { return BLAS_FUNC(slange)(&norm, &m, &n, a, &lda, work); }
    inline f64 lange(char norm, CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *work)
        { return BLAS_FUNC(dlange)(&norm, &m, &n, a, &lda, work); }
    inline f32 lange(char norm, CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *work)
        { return BLAS_FUNC(clange)(&norm, &m, &n, a, &lda, work); }
    inline f64 lange(char norm, CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *work)
        { return BLAS_FUNC(zlange)(&norm, &m, &n, a, &lda, work); }

    inline f32 lantr(char norm, char uplo, char diag, CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *work)
        { return BLAS_FUNC(slantr)(&norm, &uplo, &diag, &m, &n, a, &lda, work); }
    inline f64 lantr(char norm, char uplo, char diag, CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *work)
        { return BLAS_FUNC(dlantr)(&norm, &uplo, &diag, &m, &n, a, &lda, work); }
    inline f32 lantr(char norm, char uplo, char diag, CBLAS_INT m, CBLAS_INT n, c64 *a, CBLAS_INT lda, f32 *work)
        { return BLAS_FUNC(clantr)(&norm, &uplo, &diag, &m, &n, a, &lda, work); }
    inline f64 lantr(char norm, char uplo, char diag, CBLAS_INT m, CBLAS_INT n, c128 *a, CBLAS_INT lda, f64 *work)
        { return BLAS_FUNC(zlantr)(&norm, &uplo, &diag, &m, &n, a, &lda, work); }

    inline void larfg(CBLAS_INT n, f32 *alpha, f32 *x, CBLAS_INT incx, f32 *tau)
        { BLAS_FUNC(slarfg)(&n, alpha, x, &incx, tau); }
    inline void larfg(CBLAS_INT n, f64 *alpha, f64 *x, CBLAS_INT incx, f64 *tau)
        { BLAS_FUNC(dlarfg)(&n, alpha, x, &incx, tau); }
    inline void larfg(CBLAS_INT n, c64 *alpha, c64 *x, CBLAS_INT incx, c64 *tau)
        { BLAS_FUNC(clarfg)(&n, alpha, x, &incx, tau); }
    inline void larfg(CBLAS_INT n, c128 *alpha, c128 *x, CBLAS_INT incx, c128 *tau)
        { BLAS_FUNC(zlarfg)(&n, alpha, x, &incx, tau); }

    inline void larf(char side, CBLAS_INT m, CBLAS_INT n, f32 *v, CBLAS_INT incv, f32 tau, f32 *c, CBLAS_INT ldc, f32 *work)
        { BLAS_FUNC(slarf)(&side, &m, &n, v, &incv, &tau, c, &ldc, work); }
    inline void larf(char side, CBLAS_INT m, CBLAS_INT n, f64 *v, CBLAS_INT incv, f64 tau, f64 *c, CBLAS_INT ldc, f64 *work)
        { BLAS_FUNC(dlarf)(&side, &m, &n, v, &incv, &tau, c, &ldc, work); }
    inline void larf(char side, CBLAS_INT m, CBLAS_INT n, c64 *v, CBLAS_INT incv, c64 tau, c64 *c, CBLAS_INT ldc, c64 *work)
        { BLAS_FUNC(clarf)(&side, &m, &n, v, &incv, &tau, c, &ldc, work); }
    inline void larf(char side, CBLAS_INT m, CBLAS_INT n, c128 *v, CBLAS_INT incv, c128 tau, c128 *c, CBLAS_INT ldc, c128 *work)
        { BLAS_FUNC(zlarf)(&side, &m, &n, v, &incv, &tau, c, &ldc, work); }

    /* `lartg`'s cosine is real in every flavor; the sine follows the flavor. */
    inline void lartg(f32 f, f32 g, f32 *cs, f32 *sn, f32 *r)      { BLAS_FUNC(slartg)(&f, &g, cs, sn, r); }
    inline void lartg(f64 f, f64 g, f64 *cs, f64 *sn, f64 *r)      { BLAS_FUNC(dlartg)(&f, &g, cs, sn, r); }
    inline void lartg(c64 f, c64 g, f32 *cs, c64 *sn, c64 *r)      { BLAS_FUNC(clartg)(&f, &g, cs, sn, r); }
    inline void lartg(c128 f, c128 g, f64 *cs, c128 *sn, c128 *r)  { BLAS_FUNC(zlartg)(&f, &g, cs, sn, r); }

    /* `rot` exists only for the complex flavors: the real one is BLAS `srot`/`drot`. */
    inline void rot(CBLAS_INT n, c64 *x, CBLAS_INT incx, c64 *y, CBLAS_INT incy, f32 c, c64 s)
        { BLAS_FUNC(crot)(&n, x, &incx, y, &incy, &c, &s); }
    inline void rot(CBLAS_INT n, c128 *x, CBLAS_INT incx, c128 *y, CBLAS_INT incy, f64 c, c128 s)
        { BLAS_FUNC(zrot)(&n, x, &incx, y, &incy, &c, &s); }

    inline void ilaver(CBLAS_INT *major, CBLAS_INT *minor, CBLAS_INT *patch)
        { BLAS_FUNC(ilaver)(major, minor, patch); }

    /* `tgsyl` solves the generalized Sylvester pair.  The `.pyf` declares only the real
     * flavors, so only those are wrapped even though LAPACK ships complex ones too. */
    inline void tgsyl(char trans, CBLAS_INT ijob, CBLAS_INT m, CBLAS_INT n, f32 *a, CBLAS_INT lda, f32 *b, CBLAS_INT ldb, f32 *c, CBLAS_INT ldc, f32 *d, CBLAS_INT ldd, f32 *e, CBLAS_INT lde, f32 *f, CBLAS_INT ldf, f32 *scale, f32 *dif, f32 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(stgsyl)(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, work, &lwork, iwork, info); }
    inline void tgsyl(char trans, CBLAS_INT ijob, CBLAS_INT m, CBLAS_INT n, f64 *a, CBLAS_INT lda, f64 *b, CBLAS_INT ldb, f64 *c, CBLAS_INT ldc, f64 *d, CBLAS_INT ldd, f64 *e, CBLAS_INT lde, f64 *f, CBLAS_INT ldf, f64 *scale, f64 *dif, f64 *work, CBLAS_INT lwork, CBLAS_INT *iwork, CBLAS_INT *info)
        { BLAS_FUNC(dtgsyl)(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, work, &lwork, iwork, info); }

}  // namespace lapack
