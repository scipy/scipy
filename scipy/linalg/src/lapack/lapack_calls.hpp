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

}  // namespace lapack
