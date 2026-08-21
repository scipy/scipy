/**
 * @file
 * @brief Python wrappers for the LAPACK Symmetric/Hermitian routines.
 *
 */
#define PY_ARRAY_UNIQUE_SYMBOL scipy_lapack_ARRAY_API
#define NO_IMPORT_ARRAY
#include "lapack_helpers.hpp"
#include "lapack_calls.hpp"

namespace lapack {
    namespace capi {


        template <class T>
        static PyObject *sytrf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lower", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("sytrf", "O|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), info = 0;

            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(1, n));
            CHECK(lwork >= n || lwork == -1, lwork);

            ARRAY_OUT(CBLAS_INT, ipiv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));

            lapack::sytrf(lower ? 'L' : 'U', n, a.data<T>(), lda, ipiv.data<CBLAS_INT>(), work.data<T>(), lwork, &info);
            return make_result(a, ipiv, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *hetrf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lower", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("hetrf", "O|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), info = 0;

            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(1, n));
            CHECK(lwork >= n || lwork == -1, lwork);

            ARRAY_OUT(CBLAS_INT, ipiv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));

            lapack::hetrf(lower ? 'L' : 'U', n, a.data<T>(), lda, ipiv.data<CBLAS_INT>(), work.data<T>(), lwork, &info);
            return make_result(a, ipiv, static_cast<long long>(info));
        }


        /* The `_lwork` companions all take `(n, lower)` and hand back `work[0]`.  LAPACK's
         * query path writes only `WORK(1)` -- verified in `?sytrf`, `?sysv`, `?sysvx` and their
         * Hermitian twins -- so every array argument below is a nullptr. */

        template <class T>
        static PyObject *sytrf_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lower", nullptr};
            static constexpr Ctx<T> ctx("sytrf_lwork", "O|O", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            T work = 0;
            CBLAS_INT info = 0;
            lapack::sytrf(lower ? 'L' : 'U', n, nullptr, std::max<CBLAS_INT>(1, n), nullptr, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *hetrf_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lower", nullptr};
            static constexpr Ctx<T> ctx("hetrf_lwork", "O|O", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            T work = 0;
            CBLAS_INT info = 0;
            lapack::hetrf(lower ? 'L' : 'U', n, nullptr, std::max<CBLAS_INT>(1, n), nullptr, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sytf2(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lower", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("sytf2", "O|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), info = 0;

            ARRAY_OUT(CBLAS_INT, ipiv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));

            lapack::sytf2(lower ? 'L' : 'U', n, a.data<T>(), lda, ipiv.data<CBLAS_INT>(), &info);
            return make_result(a, ipiv, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sytrs(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "ipiv", "b", "lower", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("sytrs", "OOO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            /* `a` is `dimension(lda, n)` with `n` its *column* count, so a taller-than-wide
             * buffer is legal as long as `lda >= n`. */
            ARRAY_IN(T, a, 2);
            CBLAS_INT n = shape(a, 1), lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CHECKARRAY(lda >= n && n >= 0, a);

            ARRAY_IN(CBLAS_INT, ipiv, 1);
            CHECKARRAY(len(ipiv) == n, ipiv);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CBLAS_INT nrhs = shape(b, 1), ldb = std::max<CBLAS_INT>(1, shape(b, 0)), info = 0;
            CHECKARRAY(ldb >= n, b);

            lapack::sytrs(lower ? 'L' : 'U', n, nrhs, a.data<T>(), lda, ipiv.data<CBLAS_INT>(), b.data<T>(), ldb, &info);
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *hetrs(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "ipiv", "b", "lower", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("hetrs", "OOO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_IN(T, a, 2);
            CBLAS_INT n = shape(a, 1), lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CHECKARRAY(lda >= n && n >= 0, a);

            ARRAY_IN(CBLAS_INT, ipiv, 1);
            CHECKARRAY(len(ipiv) == n, ipiv);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CBLAS_INT nrhs = shape(b, 1), ldb = std::max<CBLAS_INT>(1, shape(b, 0)), info = 0;
            CHECKARRAY(ldb >= n, b);

            lapack::hetrs(lower ? 'L' : 'U', n, nrhs, a.data<T>(), lda, ipiv.data<CBLAS_INT>(), b.data<T>(), ldb, &info);
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sytri(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "ipiv", "lower", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("sytri", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT n = shape(a, 1), lda = std::max<CBLAS_INT>(1, shape(a, 0)), info = 0;

            ARRAY_IN(CBLAS_INT, ipiv, 1);
            CHECKARRAY(len(ipiv) == n, ipiv);

            /* `csytri`/`zsytri` want 2n where the real flavors and `?hetri` want n; the `.pyf`
             * says so by declaring `sytri` twice, once per prefix set. */
            CBLAS_INT work_len;
            if (!work_size(is_complex_v<T> ? 2LL * n : 1LL * n, &work_len)) { return nullptr; }
            ARRAY_HIDDEN(T, work, work_len);

            lapack::sytri(lower ? 'L' : 'U', n, a.data<T>(), lda, ipiv.data<CBLAS_INT>(), work.data<T>(), &info);
            return make_result(a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *hetri(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "ipiv", "lower", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("hetri", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT n = shape(a, 1), lda = std::max<CBLAS_INT>(1, shape(a, 0)), info = 0;

            ARRAY_IN(CBLAS_INT, ipiv, 1);
            CHECKARRAY(len(ipiv) == n, ipiv);
            ARRAY_HIDDEN(T, work, n);

            lapack::hetri(lower ? 'L' : 'U', n, a.data<T>(), lda, ipiv.data<CBLAS_INT>(), work.data<T>(), &info);
            return make_result(a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *syconv(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "ipiv", "lower", "way", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("syconv", "OO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);
            /* An integer here, not a letter: 0 converts, 1 reverts. */
            SCALAR_OPT(CBLAS_INT, way, 0);
            CHECK(way == 0 || way == 1, way);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), info = 0;

            ARRAY_IN(CBLAS_INT, ipiv, 1);
            CHECKARRAY(len(ipiv) == n, ipiv);
            ARRAY_OUT(T, e, 1, true, ctx.zeros(n));

            lapack::syconv(lower ? 'L' : 'U', way ? 'R' : 'C', n, a.data<T>(), lda, ipiv.data<CBLAS_INT>(), e.data<T>(), &info);
            return make_result(a, e, static_cast<long long>(info));
        }


        /* ======================= equilibration and conditioning ==================== */

        template <class T>
        static PyObject *syequb(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lower", nullptr};
            static constexpr Ctx<T> ctx("syequb", "O|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_IN(T, a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 1), lda = std::max<CBLAS_INT>(1, shape(a, 0)), work_len, info = 0;
            if (!work_size(3LL * n, &work_len)) { return nullptr; }

            ARRAY_OUT(R, s, 1, true, ctx.template zeros_as<R>(n));
            ARRAY_HIDDEN(T, work, work_len);
            R scond = 0, amax = 0;

            lapack::syequb(lower ? 'L' : 'U', n, a.data<T>(), lda, s.data<R>(), &scond, &amax, work.data<T>(), &info);
            return make_result(s, scond, amax, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *heequb(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lower", nullptr};
            static constexpr Ctx<T> ctx("heequb", "O|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_IN(T, a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 1), lda = std::max<CBLAS_INT>(1, shape(a, 0)), work_len, info = 0;
            if (!work_size(3LL * n, &work_len)) { return nullptr; }

            ARRAY_OUT(R, s, 1, true, ctx.template zeros_as<R>(n));
            ARRAY_HIDDEN(T, work, work_len);
            R scond = 0, amax = 0;

            lapack::heequb(lower ? 'L' : 'U', n, a.data<T>(), lda, s.data<R>(), &scond, &amax, work.data<T>(), &info);
            return make_result(s, scond, amax, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sycon(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "ipiv", "anorm", "lower", nullptr};
            static constexpr Ctx<T> ctx("sycon", "OOO|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_IN(T, a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), work_len, info = 0;

            ARRAY_IN(CBLAS_INT, ipiv, 1);
            CHECKARRAY(len(ipiv) == n, ipiv);
            SCALAR_REQ(R, anorm);

            if (!work_size(2LL * n, &work_len)) { return nullptr; }
            ARRAY_HIDDEN(T, work, work_len);
            R rcond = 0;

            /* The complex flavors carry no IWORK -- the estimator stays in the complex buffer,
             * as in `?gtcon` and `?gbcon`. */
            if constexpr (is_complex_v<T>) {
                lapack::sycon(lower ? 'L' : 'U', n, a.data<T>(), lda, ipiv.data<CBLAS_INT>(),
                              anorm, &rcond, work.data<T>(), &info);
            }
            else {
                ARRAY_HIDDEN(CBLAS_INT, iwork, n);
                lapack::sycon(lower ? 'L' : 'U', n, a.data<T>(), lda, ipiv.data<CBLAS_INT>(),
                              anorm, &rcond, work.data<T>(), iwork.data<CBLAS_INT>(), &info);
            }
            return make_result(rcond, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *hecon(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "ipiv", "anorm", "lower", nullptr};
            static constexpr Ctx<T> ctx("hecon", "OOO|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_IN(T, a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), work_len, info = 0;

            ARRAY_IN(CBLAS_INT, ipiv, 1);
            CHECKARRAY(len(ipiv) == n, ipiv);
            SCALAR_REQ(R, anorm);

            if (!work_size(2LL * n, &work_len)) { return nullptr; }
            ARRAY_HIDDEN(T, work, work_len);
            R rcond = 0;

            lapack::hecon(lower ? 'L' : 'U', n, a.data<T>(), lda, ipiv.data<CBLAS_INT>(), anorm, &rcond, work.data<T>(), &info);
            return make_result(rcond, static_cast<long long>(info));
        }


        /* ================================ simple drivers =========================== */

        template <class T>
        static PyObject *sysv(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "lwork", "lower", "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("sysv", "OO|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = shape(a, 0);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), ldb = shape(b, 0), info = 0;

            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(1, n));
            CHECK(lwork >= 1 || lwork == -1, lwork);

            ARRAY_OUT(CBLAS_INT, ipiv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));

            lapack::sysv(lower ? 'L' : 'U', n, nrhs, a.data<T>(), lda, ipiv.data<CBLAS_INT>(), b.data<T>(), ldb, work.data<T>(), lwork, &info);
            return make_result(a, ipiv, b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *hesv(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "lwork", "lower", "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("hesv", "OO|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = shape(a, 0);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), ldb = shape(b, 0), info = 0;

            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(1, n));
            CHECK(lwork >= 1 || lwork == -1, lwork);

            ARRAY_OUT(CBLAS_INT, ipiv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));

            lapack::hesv(lower ? 'L' : 'U', n, nrhs, a.data<T>(), lda, ipiv.data<CBLAS_INT>(),
                         b.data<T>(), ldb, work.data<T>(), lwork, &info);
            return make_result(a, ipiv, b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sysv_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lower", nullptr};
            static constexpr Ctx<T> ctx("sysv_lwork", "O|O", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            T work = 0;
            CBLAS_INT info = 0;
            lapack::sysv(lower ? 'L' : 'U', n, 1, nullptr, n, nullptr, nullptr, n, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *hesv_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lower", nullptr};
            static constexpr Ctx<T> ctx("hesv_lwork", "O|O", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            T work = 0;
            CBLAS_INT info = 0;
            lapack::hesv(lower ? 'L' : 'U', n, 1, nullptr, n, nullptr, nullptr, n, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        /* ================================ expert drivers =========================== */

        template <class T>
        static PyObject *sysvx(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "af", "ipiv", "lwork", "factored", "lower", "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("sysvx", "OO|OOOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            /* Integer flags where `gesvx` and `posvx` take option letters -- this file's own
             * spelling, kept as it is. */
            SCALAR_OPT(CBLAS_INT, factored, 0);
            CHECK(factored == 0 || factored == 1, factored);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = shape(a, 0);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), ldb = shape(b, 0), ldx = n, info = 0;

            ARRAY_OUT(T, af, 2, true, ctx.zeros(n, n));
            CHECKARRAY(shape(af, 0) == n && shape(af, 1) == n, af);
            CBLAS_INT ldaf = shape(af, 0);

            ARRAY_OUT(CBLAS_INT, ipiv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            CHECKARRAY(len(ipiv) == n, ipiv);

            CBLAS_INT minimum_lwork;
            if (!work_size(is_complex_v<T> ? 2LL * n : 3LL * n, &minimum_lwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(1, 3 * n));
            CHECK(lwork >= minimum_lwork || lwork == -1, lwork);

            ARRAY_OUT(T, x, 2, true, ctx.zeros(n, nrhs));
            ARRAY_OUT(R, ferr, 1, true, ctx.template zeros_as<R>(nrhs));
            ARRAY_OUT(R, berr, 1, true, ctx.template zeros_as<R>(nrhs));

            using W = std::conditional_t<is_complex_v<T>, R, CBLAS_INT>;
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            ARRAY_HIDDEN(W, irwork, n);
            R rcond = 0;

            lapack::sysvx(factored ? 'F' : 'N', lower ? 'L' : 'U', n, nrhs, a.data<T>(), lda,
                          af.data<T>(), ldaf, ipiv.data<CBLAS_INT>(), b.data<T>(), ldb,
                          x.data<T>(), ldx, &rcond, ferr.data<R>(), berr.data<R>(),
                          work.data<T>(), lwork, irwork.data<W>(), &info);
            return make_result(a, af, ipiv, b, x, rcond, ferr, berr,
                               static_cast<long long>(info));
        }


        template <class T>
        static PyObject *hesvx(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "af", "ipiv", "lwork", "factored", "lower", "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("hesvx", "OO|OOOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, factored, 0);
            CHECK(factored == 0 || factored == 1, factored);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            /* `a` and `b` are `intent(in,copy)` here, not `intent(in,copy,out)` as in `sysvx`:
             * the Hermitian driver hands back seven values, not nine. */
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = shape(a, 0);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), ldb = shape(b, 0), info = 0;

            ARRAY_OUT(T, af, 2, true, ctx.zeros(n, n));
            CHECKARRAY(shape(af, 0) == n && shape(af, 1) == n, af);
            CBLAS_INT ldaf = shape(af, 0);

            ARRAY_OUT(CBLAS_INT, ipiv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            CHECKARRAY(len(ipiv) == n, ipiv);

            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(1, 2 * n));
            CHECK(lwork >= 1 || lwork == -1, lwork);

            ARRAY_OUT(T, x, 2, true, ctx.zeros(n, nrhs));
            CBLAS_INT ldx = shape(x, 0);
            ARRAY_OUT(R, ferr, 1, true, ctx.template zeros_as<R>(nrhs));
            ARRAY_OUT(R, berr, 1, true, ctx.template zeros_as<R>(nrhs));

            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            ARRAY_HIDDEN(R, rwork, n);
            R rcond = 0;

            lapack::hesvx(factored ? 'F' : 'N', lower ? 'L' : 'U', n, nrhs, a.data<T>(), lda,
                          af.data<T>(), ldaf, ipiv.data<CBLAS_INT>(), b.data<T>(), ldb,
                          x.data<T>(), ldx, &rcond, ferr.data<R>(), berr.data<R>(),
                          work.data<T>(), lwork, rwork.data<R>(), &info);
            return make_result(af, ipiv, x, rcond, ferr, berr, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sysvx_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lower", nullptr};
            static constexpr Ctx<T> ctx("sysvx_lwork", "O|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            using W = std::conditional_t<is_complex_v<T>, R, CBLAS_INT>;
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            T work = 0;
            CBLAS_INT info = 0;
            lapack::sysvx('N', lower ? 'L' : 'U', n, 1, nullptr, n, nullptr, n, nullptr, nullptr, n, nullptr, n, nullptr, nullptr,
                          nullptr, &work, -1, static_cast<W *>(nullptr), &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *hesvx_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lower", nullptr};
            static constexpr Ctx<T> ctx("hesvx_lwork", "O|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            T work = 0;
            CBLAS_INT info = 0;
            lapack::hesvx('N', lower ? 'L' : 'U', n, 1, nullptr, n, nullptr, n, nullptr, nullptr, n, nullptr, n, nullptr,
                          nullptr, nullptr, &work, -1, static_cast<R *>(nullptr), &info);
            return make_result(work, static_cast<long long>(info));
        }


        /* ============================== eigenproblems ============================== */

        /* From here on a `sy`/`he` pair is *one* wrapper wherever the two `.pyf` blocks declare
         * the same Python signature: the complex flavor only adds a hidden `rwork` and uses a
         * different `lwork` formula.  The Python name comes from a compile-time ternary on the
         * `Ctx`, the way `gelsd` picks its kwlist.  Five pairs do not qualify -- `syevd`,
         * `syevr`, `sygvd` and two `_lwork` queries, where the complex side exposes an extra
         * `lrwork` argument or an extra return value -- and those stay as separate wrappers.
         *
         * The flag spelling is not uniform here and is kept as the `.pyf` has it: `syev` and
         * `syevd` take integer `compute_v`/`lower`, `syevr`/`syevx` add a `range` *letter*, and
         * the `sygv*` family takes `jobz`/`uplo` letters instead -- with `uplo` defaulting to
         * 'L' where `lower` defaults to 0, i.e. to the *other* triangle. */

        template <class T>
        static PyObject *syev(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "compute_v", "lower", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "heev" : "syev", "O|OOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, compute_v, 1);  CHECK(compute_v == 0 || compute_v == 1, compute_v);
            SCALAR_OPT(CBLAS_INT, lower, 0);      CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), info = 0;

            CBLAS_INT minimum_lwork;
            if (!work_size(is_complex_v<T> ? 2LL * n - 1 : 3LL * n - 1, &minimum_lwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, minimum_lwork);
            CHECK(lwork >= minimum_lwork, lwork);

            ARRAY_OUT(R, w, 1, true, ctx.template zeros_as<R>(n));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));

            if constexpr (is_complex_v<T>) {
                CBLAS_INT rwork_len;
                if (!work_size(3LL * n - 1, &rwork_len)) { return nullptr; }
                ARRAY_HIDDEN(R, rwork, rwork_len);
                lapack::syev(compute_v ? 'V' : 'N', lower ? 'L' : 'U', n, a.data<T>(), lda,
                             w.data<R>(), work.data<T>(), lwork, rwork.data<R>(), &info);
            }
            else {
                lapack::syev(compute_v ? 'V' : 'N', lower ? 'L' : 'U', n, a.data<T>(), lda,
                             w.data<R>(), work.data<T>(), lwork, &info);
            }
            return make_result(w, a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *syev_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lower", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "heev_lwork" : "syev_lwork", "O|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            T work = 0;
            CBLAS_INT info = 0, lda = std::max<CBLAS_INT>(1, n);
            if constexpr (is_complex_v<T>) {
                lapack::syev('N', lower ? 'L' : 'U', n, nullptr, lda, nullptr, &work, -1,
                             static_cast<R *>(nullptr), &info);
            }
            else {
                lapack::syev('N', lower ? 'L' : 'U', n, nullptr, lda, nullptr, &work, -1, &info);
            }
            return make_result(work, static_cast<long long>(info));
        }


        /* `syevd`/`heevd` split: the complex block adds an `lrwork` keyword and its own `rwork`
         * length, and its Fortran call interleaves them differently. */

        template <class T>
        static PyObject *syevd(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "compute_v", "lower", "lwork", "liwork",
                                           "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("syevd", "O|OOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, compute_v, 1);  CHECK(compute_v == 0 || compute_v == 1, compute_v);
            SCALAR_OPT(CBLAS_INT, lower, 0);      CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), info = 0;

            CBLAS_INT minimum_lwork, minimum_liwork;
            if (!work_size(compute_v ? 1LL + 6LL * n + 2LL * n * n : 2LL * n + 1, &minimum_lwork) ||
                !work_size(compute_v ? 3LL + 5LL * n : 1LL, &minimum_liwork)) {
                return nullptr;
            }
            SCALAR_OPT(CBLAS_INT, lwork, minimum_lwork);
            CHECK(lwork >= minimum_lwork || (n == 1 && lwork >= 1), lwork);
            SCALAR_OPT(CBLAS_INT, liwork, minimum_liwork);
            CHECK(liwork >= minimum_liwork || (n == 1 && liwork >= 1), liwork);

            ARRAY_OUT(T, w, 1, true, ctx.zeros(n));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            ARRAY_HIDDEN(CBLAS_INT, iwork, std::max<CBLAS_INT>(1, liwork));

            lapack::syevd(compute_v ? 'V' : 'N', lower ? 'L' : 'U', n, a.data<T>(), lda,
                          w.data<T>(), work.data<T>(), lwork, iwork.data<CBLAS_INT>(), liwork,
                          &info);
            return make_result(w, a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *heevd(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "compute_v", "lower", "lwork", "liwork", "lrwork",
                                           "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("heevd", "O|OOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, compute_v, 1);  CHECK(compute_v == 0 || compute_v == 1, compute_v);
            SCALAR_OPT(CBLAS_INT, lower, 0);      CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), info = 0;

            CBLAS_INT minimum_lwork, minimum_liwork, minimum_lrwork;
            if (!work_size(compute_v ? 2LL * n + 1LL * n * n : 1LL * n + 1, &minimum_lwork) ||
                !work_size(compute_v ? 3LL + 5LL * n : 1LL, &minimum_liwork) ||
                !work_size(compute_v ? 1LL + 5LL * n + 2LL * n * n : 1LL * n, &minimum_lrwork)) {
                return nullptr;
            }
            SCALAR_OPT(CBLAS_INT, lwork, minimum_lwork);
            CHECK(lwork >= minimum_lwork || (n == 1 && lwork >= 1), lwork);
            SCALAR_OPT(CBLAS_INT, liwork, minimum_liwork);
            CHECK(liwork >= minimum_liwork || (n == 1 && liwork >= 1), liwork);
            SCALAR_OPT(CBLAS_INT, lrwork, minimum_lrwork);
            CHECK(lrwork >= minimum_lrwork || (n == 1 && lrwork >= 1), lrwork);

            ARRAY_OUT(R, w, 1, true, ctx.template zeros_as<R>(n));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            ARRAY_HIDDEN(R, rwork, std::max<CBLAS_INT>(1, lrwork));
            ARRAY_HIDDEN(CBLAS_INT, iwork, std::max<CBLAS_INT>(1, liwork));

            lapack::heevd(compute_v ? 'V' : 'N', lower ? 'L' : 'U', n, a.data<T>(), lda,
                          w.data<R>(), work.data<T>(), lwork, rwork.data<R>(), lrwork,
                          iwork.data<CBLAS_INT>(), liwork, &info);
            return make_result(w, a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *syevd_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "compute_v", "lower", nullptr};
            static constexpr Ctx<T> ctx("syevd_lwork", "O|OO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, compute_v, 1);  CHECK(compute_v == 0 || compute_v == 1, compute_v);
            SCALAR_OPT(CBLAS_INT, lower, 0);      CHECK(lower == 0 || lower == 1, lower);

            T work = 0;
            CBLAS_INT iwork = 0, info = 0;
            lapack::syevd(compute_v ? 'V' : 'N', lower ? 'L' : 'U', n, nullptr,
                          std::max<CBLAS_INT>(1, n), nullptr, &work, -1, &iwork, -1, &info);
            return make_result(work, static_cast<long long>(iwork), static_cast<long long>(info));
        }


        template <class T>
        static PyObject *heevd_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "compute_v", "lower", nullptr};
            static constexpr Ctx<T> ctx("heevd_lwork", "O|OO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, compute_v, 1);  CHECK(compute_v == 0 || compute_v == 1, compute_v);
            SCALAR_OPT(CBLAS_INT, lower, 0);      CHECK(lower == 0 || lower == 1, lower);

            /* Four values back where the real query hands back three: the complex routine has a
             * third workspace to size. */
            T work = 0;
            R rwork = 0;
            CBLAS_INT iwork = 0, info = 0;
            lapack::heevd(compute_v ? 'V' : 'N', lower ? 'L' : 'U', n, nullptr,
                          std::max<CBLAS_INT>(1, n), nullptr, &work, -1, &rwork, -1, &iwork, -1,
                          &info);
            return make_result(work, static_cast<long long>(iwork), rwork,
                               static_cast<long long>(info));
        }


        /* `syevr`/`heevr` split: the complex block adds `lrwork`, and its `isuppz` is always
         * `2 * max(1, n)` rather than the conditional size -- an MKL work-around the `.pyf`
         * applies to the Hermitian half only.  Both are kept exactly as declared. */

        template <class T>
        static PyObject *syevr(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "compute_v", "range", "lower", "vl", "vu", "il",
                                           "iu", "abstol", "lwork", "liwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("syevr", "O|OOOOOOOOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, compute_v, 1);  CHECK(compute_v == 0 || compute_v == 1, compute_v);
            SCALAR_OPT(char, range, 'A');         CHECK(range == 'A' || range == 'V' || range == 'I', range);
            SCALAR_OPT(CBLAS_INT, lower, 0);      CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, n), m = 0, info = 0;

            SCALAR_OPT(CBLAS_INT, il, 1);   CHECK(il >= 1 && il <= n, il);
            SCALAR_OPT(CBLAS_INT, iu, n);   CHECK(n >= iu && iu >= il, iu);
            SCALAR_OPT(T, vl, 0.0);
            SCALAR_OPT(T, vu, 1.0);         CHECK(vu >= vl, vu);
            SCALAR_OPT(T, abstol, 0.0);

            CBLAS_INT minimum_lwork, minimum_liwork;
            if (!work_size(n <= 1 ? 1LL : 26LL * n, &minimum_lwork) ||
                !work_size(n <= 1 ? 1LL : 10LL * n, &minimum_liwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(1, 26 * n));
            CHECK(lwork >= minimum_lwork || lwork == -1, lwork);
            SCALAR_OPT(CBLAS_INT, liwork, std::max<CBLAS_INT>(1, 10 * n));
            CHECK(liwork >= minimum_liwork || liwork == -1, liwork);

            CBLAS_INT z_rows = compute_v ? std::max<CBLAS_INT>(0, n) : 0;
            CBLAS_INT z_cols = compute_v ? (range == 'I' ? iu - il + 1 : std::max<CBLAS_INT>(1, n)) : 0;
            /* LAPACK writes `isuppz` on the MRRR path -- `range` 'A', or 'I' spanning the
             * whole index range -- and also from the `n == 1` early return, which is guarded
             * by `jobz` alone.  The `.pyf` omits that last case, so `range='V'` at `n == 1`
             * handed LAPACK a zero-length buffer and it wrote two integers past the end. */
            CBLAS_INT isuppz_len =
                compute_v && (range == 'A' || (range == 'I' && iu - il + 1 == n) || n == 1)
                    ? 2 * n : 0;

            ARRAY_OUT(T, w, 1, true, ctx.zeros(n));
            ARRAY_OUT(T, z, 2, true, ctx.zeros(z_rows, z_cols));
            ARRAY_OUT(CBLAS_INT, isuppz, 1, true, ctx.template zeros_as<CBLAS_INT>(isuppz_len));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            ARRAY_HIDDEN(CBLAS_INT, iwork, std::max<CBLAS_INT>(1, liwork));
            CBLAS_INT ldz = std::max<CBLAS_INT>(1, shape(z, 0));

            lapack::syevr(compute_v ? 'V' : 'N', range, lower ? 'L' : 'U', n, a.data<T>(), lda,
                          vl, vu, il, iu, abstol, &m, w.data<T>(), z.data<T>(), ldz,
                          isuppz.data<CBLAS_INT>(), work.data<T>(), lwork,
                          iwork.data<CBLAS_INT>(), liwork, &info);
            return make_result(w, z, static_cast<long long>(m), isuppz,
                               static_cast<long long>(info));
        }


        template <class T>
        static PyObject *heevr(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "compute_v", "range", "lower", "vl", "vu", "il",
                                           "iu", "abstol", "lwork", "lrwork", "liwork",
                                           "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("heevr", "O|OOOOOOOOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, compute_v, 1);  CHECK(compute_v == 0 || compute_v == 1, compute_v);
            SCALAR_OPT(char, range, 'A');         CHECK(range == 'A' || range == 'V' || range == 'I', range);
            SCALAR_OPT(CBLAS_INT, lower, 0);      CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, n), m = 0, info = 0;

            SCALAR_OPT(CBLAS_INT, il, 1);   CHECK(il >= 1 && il <= n, il);
            SCALAR_OPT(CBLAS_INT, iu, n);   CHECK(n >= iu && iu >= il, iu);
            SCALAR_OPT(R, vl, 0.0);
            SCALAR_OPT(R, vu, 1.0);         CHECK(vu >= vl, vu);
            SCALAR_OPT(R, abstol, 0.0);

            CBLAS_INT minimum_lwork, minimum_lrwork, minimum_liwork;
            if (!work_size(n <= 1 ? 1LL : 2LL * n, &minimum_lwork) ||
                !work_size(n <= 1 ? 1LL : 24LL * n, &minimum_lrwork) ||
                !work_size(n <= 1 ? 1LL : 10LL * n, &minimum_liwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(1, 2 * n));
            CHECK(lwork >= minimum_lwork || lwork == -1, lwork);
            SCALAR_OPT(CBLAS_INT, lrwork, std::max<CBLAS_INT>(1, 24 * n));
            CHECK(lrwork >= minimum_lrwork || lrwork == -1, lrwork);
            SCALAR_OPT(CBLAS_INT, liwork, std::max<CBLAS_INT>(1, 10 * n));
            CHECK(liwork >= minimum_liwork || liwork == -1, liwork);

            CBLAS_INT z_rows = compute_v ? std::max<CBLAS_INT>(0, n) : 0;
            CBLAS_INT z_cols = compute_v ? (range == 'I' ? iu - il + 1 : std::max<CBLAS_INT>(1, n)) : 0;

            ARRAY_OUT(R, w, 1, true, ctx.template zeros_as<R>(n));
            ARRAY_OUT(T, z, 2, true, ctx.zeros(z_rows, z_cols));
            /* Same rule as `syevr`.  The `.pyf` allocated `2 * max(1, n)` unconditionally on
             * this half to work around an MKL bug that wrote `isuppz` where the reference did
             * not; its own comment expected a fix in MKL 2020 update 2, and MKL 2024.2.2 was
             * measured writing exactly the reference's entries, so the work-around is retired
             * and the two halves agree again. */
            CBLAS_INT isuppz_len =
                compute_v && (range == 'A' || (range == 'I' && iu - il + 1 == n) || n == 1)
                    ? 2 * n : 0;
            ARRAY_OUT(CBLAS_INT, isuppz, 1, true, ctx.template zeros_as<CBLAS_INT>(isuppz_len));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            ARRAY_HIDDEN(R, rwork, std::max<CBLAS_INT>(1, lrwork));
            ARRAY_HIDDEN(CBLAS_INT, iwork, std::max<CBLAS_INT>(1, liwork));
            CBLAS_INT ldz = std::max<CBLAS_INT>(1, shape(z, 0));

            lapack::heevr(compute_v ? 'V' : 'N', range, lower ? 'L' : 'U', n, a.data<T>(), lda,
                          vl, vu, il, iu, abstol, &m, w.data<R>(), z.data<T>(), ldz,
                          isuppz.data<CBLAS_INT>(), work.data<T>(), lwork, rwork.data<R>(),
                          lrwork, iwork.data<CBLAS_INT>(), liwork, &info);
            return make_result(w, z, static_cast<long long>(m), isuppz,
                               static_cast<long long>(info));
        }


        template <class T>
        static PyObject *syevr_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lower", nullptr};
            static constexpr Ctx<T> ctx("syevr_lwork", "O|O", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            T work = 0;
            CBLAS_INT iwork = 0, m = 0, info = 0, ld = std::max<CBLAS_INT>(1, n);
            lapack::syevr('N', 'A', lower ? 'L' : 'U', n, nullptr, ld, 0, 1, 1, 0, 0, &m, nullptr,
                          nullptr, ld, nullptr, &work, -1, &iwork, -1, &info);
            return make_result(work, static_cast<long long>(iwork), static_cast<long long>(info));
        }


        template <class T>
        static PyObject *heevr_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lower", nullptr};
            static constexpr Ctx<T> ctx("heevr_lwork", "O|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            /* `work, rwork, iwork` -- note the order differs from `heevd_lwork`, which reports
             * `work, iwork, rwork`.  Both follow their own `.pyf` argument list. */
            T work = 0;
            R rwork = 0;
            CBLAS_INT iwork = 0, m = 0, info = 0, ld = std::max<CBLAS_INT>(1, n);
            lapack::heevr('N', 'A', lower ? 'L' : 'U', n, nullptr, ld, 0, 1, 1, 0, 0, &m, nullptr,
                          nullptr, ld, nullptr, &work, -1, &rwork, -1, &iwork, -1, &info);
            return make_result(work, rwork, static_cast<long long>(iwork),
                               static_cast<long long>(info));
        }


        template <class T>
        static PyObject *syevx(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "compute_v", "range", "lower", "vl", "vu", "il",
                                           "iu", "abstol", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "heevx" : "syevx", "O|OOOOOOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, compute_v, 1);  CHECK(compute_v == 0 || compute_v == 1, compute_v);
            SCALAR_OPT(char, range, 'A');         CHECK(range == 'A' || range == 'V' || range == 'I', range);
            SCALAR_OPT(CBLAS_INT, lower, 0);      CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), m = 0, info = 0;

            SCALAR_OPT(CBLAS_INT, il, 1);   CHECK(il >= 1 && il <= n, il);
            SCALAR_OPT(CBLAS_INT, iu, n);   CHECK(n >= iu && iu >= il, iu);
            SCALAR_OPT(R, vl, 0.0);
            SCALAR_OPT(R, vu, 1.0);         CHECK(vu > vl, vu);
            SCALAR_OPT(R, abstol, 0.0);

            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(1, (is_complex_v<T> ? 2 : 8) * n));
            CHECK(lwork >= 1 || lwork == -1, lwork);

            CBLAS_INT z_rows = compute_v ? std::max<CBLAS_INT>(0, n) : 0;
            CBLAS_INT z_cols = compute_v ? (range == 'I' ? iu - il + 1 : std::max<CBLAS_INT>(1, n)) : 0;
            CBLAS_INT iwork_len;
            if (!work_size(5LL * n, &iwork_len)) { return nullptr; }

            ARRAY_OUT(R, w, 1, true, ctx.template zeros_as<R>(n));
            ARRAY_OUT(T, z, 2, true, ctx.zeros(z_rows, z_cols));
            ARRAY_OUT(CBLAS_INT, ifail, 1, true,
                      ctx.template zeros_as<CBLAS_INT>(compute_v ? n : 0));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            ARRAY_HIDDEN(CBLAS_INT, iwork, iwork_len);
            CBLAS_INT ldz = std::max<CBLAS_INT>(1, shape(z, 0));

            if constexpr (is_complex_v<T>) {
                CBLAS_INT rwork_len;
                if (!work_size(7LL * n, &rwork_len)) { return nullptr; }
                ARRAY_HIDDEN(R, rwork, rwork_len);
                lapack::syevx(compute_v ? 'V' : 'N', range, lower ? 'L' : 'U', n, a.data<T>(),
                              lda, vl, vu, il, iu, abstol, &m, w.data<R>(), z.data<T>(), ldz,
                              work.data<T>(), lwork, rwork.data<R>(), iwork.data<CBLAS_INT>(),
                              ifail.data<CBLAS_INT>(), &info);
            }
            else {
                lapack::syevx(compute_v ? 'V' : 'N', range, lower ? 'L' : 'U', n, a.data<T>(),
                              lda, vl, vu, il, iu, abstol, &m, w.data<R>(), z.data<T>(), ldz,
                              work.data<T>(), lwork, iwork.data<CBLAS_INT>(),
                              ifail.data<CBLAS_INT>(), &info);
            }
            return make_result(w, z, static_cast<long long>(m), ifail,
                               static_cast<long long>(info));
        }


        template <class T>
        static PyObject *syevx_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lower", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "heevx_lwork" : "syevx_lwork", "O|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            T work = 0;
            CBLAS_INT m = 0, info = 0, ld = std::max<CBLAS_INT>(1, n);
            if constexpr (is_complex_v<T>) {
                lapack::syevx('N', 'A', lower ? 'L' : 'U', n, nullptr, ld, 0, 1, 1, 0, 0, &m,
                              nullptr, nullptr, ld, &work, -1, static_cast<R *>(nullptr),
                              nullptr, nullptr, &info);
            }
            else {
                lapack::syevx('N', 'A', lower ? 'L' : 'U', n, nullptr, ld, 0, 1, 1, 0, 0, &m,
                              nullptr, nullptr, ld, &work, -1, nullptr, nullptr, &info);
            }
            return make_result(work, static_cast<long long>(info));
        }


        /* The generalized problems.  `itype` selects which of the three forms is solved; the
         * `.pyf` writes its guard as `check(itype > 0 || itype < 4)`, which is true for every
         * integer, so a bogus value reached LAPACK and came back as `info = -1` with an XERBLA
         * message on stderr.  Spelled here the way the same file writes it for `sygst`. */

        template <class T>
        static PyObject *sygv(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "itype", "jobz", "uplo", "lwork",
                                           "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "hegv" : "sygv", "OO|OOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, itype, 1);  CHECK(itype == 1 || itype == 2 || itype == 3, itype);
            SCALAR_OPT(char, jobz, 'V');      CHECK(jobz == 'N' || jobz == 'V', jobz);
            SCALAR_OPT(char, uplo, 'L');      CHECK(uplo == 'U' || uplo == 'L', uplo);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, n), info = 0;

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == shape(b, 1) && shape(b, 0) == n, b);
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0));

            CBLAS_INT default_lwork;
            if (!work_size(is_complex_v<T> ? 2LL * n - 1 : 3LL * n - 1, &default_lwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);
            CHECK(lwork > 0 || lwork == -1, lwork);

            ARRAY_OUT(R, w, 1, true, ctx.template zeros_as<R>(n));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));

            if constexpr (is_complex_v<T>) {
                CBLAS_INT rwork_len;
                if (!work_size(3LL * n - 2, &rwork_len)) { return nullptr; }
                ARRAY_HIDDEN(R, rwork, rwork_len);
                lapack::sygv(itype, jobz, uplo, n, a.data<T>(), lda, b.data<T>(), ldb,
                             w.data<R>(), work.data<T>(), lwork, rwork.data<R>(), &info);
            }
            else {
                lapack::sygv(itype, jobz, uplo, n, a.data<T>(), lda, b.data<T>(), ldb,
                             w.data<R>(), work.data<T>(), lwork, &info);
            }
            return make_result(w, a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sygv_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            /* `uplo` here where `syev_lwork` takes `lower` -- the `.pyf`'s own inconsistency. */
            static const char *kwlist[] = {"n", "uplo", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "hegv_lwork" : "sygv_lwork", "O|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(char, uplo, 'L');
            CHECK(uplo == 'U' || uplo == 'L', uplo);

            T work = 0;
            CBLAS_INT info = 0, ld = std::max<CBLAS_INT>(1, n);
            if constexpr (is_complex_v<T>) {
                lapack::sygv(1, 'N', uplo, n, nullptr, ld, nullptr, ld, nullptr, &work, -1,
                             static_cast<R *>(nullptr), &info);
            }
            else {
                lapack::sygv(1, 'N', uplo, n, nullptr, ld, nullptr, ld, nullptr, &work, -1, &info);
            }
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sygvd(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "itype", "jobz", "uplo", "lwork", "liwork",
                                           "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("sygvd", "OO|OOOOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, itype, 1);  CHECK(itype == 1 || itype == 2 || itype == 3, itype);
            SCALAR_OPT(char, jobz, 'V');      CHECK(jobz == 'N' || jobz == 'V', jobz);
            SCALAR_OPT(char, uplo, 'L');      CHECK(uplo == 'U' || uplo == 'L', uplo);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, n), info = 0;

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == shape(b, 1) && shape(b, 0) == n, b);
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0));

            CBLAS_INT default_lwork, default_liwork;
            if (!work_size(jobz == 'N' ? 2LL * n + 1 : 1LL + 6LL * n + 2LL * n * n, &default_lwork) ||
                !work_size(jobz == 'N' ? 1LL : 5LL * n + 3, &default_liwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);
            CHECK(lwork > 0 || lwork == -1, lwork);
            SCALAR_OPT(CBLAS_INT, liwork, default_liwork);
            CHECK(liwork > 0 || liwork == -1, liwork);

            ARRAY_OUT(T, w, 1, true, ctx.zeros(n));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            ARRAY_HIDDEN(CBLAS_INT, iwork, std::max<CBLAS_INT>(1, liwork));

            lapack::sygvd(itype, jobz, uplo, n, a.data<T>(), lda, b.data<T>(), ldb, w.data<T>(),
                          work.data<T>(), lwork, iwork.data<CBLAS_INT>(), liwork, &info);
            return make_result(w, a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *hegvd(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "itype", "jobz", "uplo", "lwork", "lrwork",
                                           "liwork", "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("hegvd", "OO|OOOOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, itype, 1);  CHECK(itype == 1 || itype == 2 || itype == 3, itype);
            SCALAR_OPT(char, jobz, 'V');      CHECK(jobz == 'N' || jobz == 'V', jobz);
            SCALAR_OPT(char, uplo, 'L');      CHECK(uplo == 'U' || uplo == 'L', uplo);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, n), info = 0;

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == shape(b, 1) && shape(b, 0) == n, b);
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0));

            CBLAS_INT default_lwork, default_lrwork, default_liwork;
            if (!work_size(jobz == 'N' ? 1LL * n + 1 : 1LL * n * (n + 2), &default_lwork) ||
                !work_size(jobz == 'N' ? 1LL * n : 2LL * n * n + 5LL * n + 1, &default_lrwork) ||
                !work_size(jobz == 'N' ? 1LL : 5LL * n + 3, &default_liwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);
            CHECK(lwork > 0 || lwork == -1, lwork);
            SCALAR_OPT(CBLAS_INT, lrwork, default_lrwork);
            CHECK(lrwork > 0 || lrwork == -1, lrwork);
            SCALAR_OPT(CBLAS_INT, liwork, default_liwork);
            CHECK(liwork > 0 || liwork == -1, liwork);

            ARRAY_OUT(R, w, 1, true, ctx.template zeros_as<R>(n));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            ARRAY_HIDDEN(R, rwork, std::max<CBLAS_INT>(1, lrwork));
            ARRAY_HIDDEN(CBLAS_INT, iwork, std::max<CBLAS_INT>(1, liwork));

            lapack::hegvd(itype, jobz, uplo, n, a.data<T>(), lda, b.data<T>(), ldb, w.data<R>(),
                          work.data<T>(), lwork, rwork.data<R>(), lrwork,
                          iwork.data<CBLAS_INT>(), liwork, &info);
            return make_result(w, a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sygvx(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "itype", "jobz", "range", "uplo", "vl", "vu",
                                           "il", "iu", "abstol", "lwork", "overwrite_a",
                                           "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "hegvx" : "sygvx", "OO|OOOOOOOOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, itype, 1);  CHECK(itype == 1 || itype == 2 || itype == 3, itype);
            SCALAR_OPT(char, jobz, 'V');      CHECK(jobz == 'N' || jobz == 'V', jobz);
            SCALAR_OPT(char, uplo, 'L');      CHECK(uplo == 'U' || uplo == 'L', uplo);
            SCALAR_OPT(char, range, 'A');     CHECK(range == 'A' || range == 'V' || range == 'I', range);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, n), m = 0, info = 0;

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == shape(b, 1) && shape(b, 0) == n, b);
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0));

            SCALAR_OPT(CBLAS_INT, il, 1);
            SCALAR_OPT(CBLAS_INT, iu, n);
            SCALAR_OPT(R, vl, 0.0);
            SCALAR_OPT(R, vu, 1.0);   CHECK(vu > vl, vu);
            SCALAR_OPT(R, abstol, 0.0);

            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(1, (is_complex_v<T> ? 2 : 8) * n));
            CHECK(lwork >= 1 || lwork == -1, lwork);

            /* `jobz`, not `compute_v`, drives the output shapes in this family. */
            CBLAS_INT z_rows = jobz == 'V' ? std::max<CBLAS_INT>(0, n) : 0;
            CBLAS_INT z_cols = jobz == 'V' ? (range == 'I' ? iu - il + 1 : std::max<CBLAS_INT>(1, n)) : 0;
            CBLAS_INT iwork_len;
            if (!work_size(5LL * n, &iwork_len)) { return nullptr; }

            ARRAY_OUT(R, w, 1, true, ctx.template zeros_as<R>(n));
            ARRAY_OUT(T, z, 2, true, ctx.zeros(z_rows, z_cols));
            ARRAY_OUT(CBLAS_INT, ifail, 1, true,
                      ctx.template zeros_as<CBLAS_INT>(jobz == 'N' ? 0 : n));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            ARRAY_HIDDEN(CBLAS_INT, iwork, iwork_len);
            CBLAS_INT ldz = std::max<CBLAS_INT>(1, shape(z, 0));

            if constexpr (is_complex_v<T>) {
                CBLAS_INT rwork_len;
                if (!work_size(7LL * n, &rwork_len)) { return nullptr; }
                ARRAY_HIDDEN(R, rwork, rwork_len);
                lapack::sygvx(itype, jobz, range, uplo, n, a.data<T>(), lda, b.data<T>(), ldb,
                              vl, vu, il, iu, abstol, &m, w.data<R>(), z.data<T>(), ldz,
                              work.data<T>(), lwork, rwork.data<R>(), iwork.data<CBLAS_INT>(),
                              ifail.data<CBLAS_INT>(), &info);
            }
            else {
                lapack::sygvx(itype, jobz, range, uplo, n, a.data<T>(), lda, b.data<T>(), ldb,
                              vl, vu, il, iu, abstol, &m, w.data<R>(), z.data<T>(), ldz,
                              work.data<T>(), lwork, iwork.data<CBLAS_INT>(),
                              ifail.data<CBLAS_INT>(), &info);
            }
            return make_result(w, z, static_cast<long long>(m), ifail,
                               static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sygvx_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "uplo", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "hegvx_lwork" : "sygvx_lwork", "O|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(char, uplo, 'L');
            CHECK(uplo == 'U' || uplo == 'L', uplo);

            T work = 0;
            CBLAS_INT m = 0, info = 0, ld = std::max<CBLAS_INT>(1, n);
            if constexpr (is_complex_v<T>) {
                lapack::sygvx(1, 'N', 'A', uplo, n, nullptr, ld, nullptr, ld, 0, 1, 1, 0, 0, &m,
                              nullptr, nullptr, ld, &work, -1, static_cast<R *>(nullptr),
                              nullptr, nullptr, &info);
            }
            else {
                lapack::sygvx(1, 'N', 'A', uplo, n, nullptr, ld, nullptr, ld, 0, 1, 1, 0, 0, &m,
                              nullptr, nullptr, ld, &work, -1, nullptr, nullptr, &info);
            }
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sytrd(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lower", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "hetrd" : "sytrd", "O|OOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 1), lda = std::max<CBLAS_INT>(shape(a, 0), 1), info = 0;

            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(n, 1));
            CHECK(lwork >= 1 || lwork == -1, lwork);

            /* `d` and `e` are real for every flavor -- the tridiagonal form of a Hermitian
             * matrix is real symmetric -- while `tau` follows the flavor. */
            CBLAS_INT nm1 = n > 0 ? n - 1 : 0;
            ARRAY_OUT(R, d, 1, true, ctx.template zeros_as<R>(n));
            ARRAY_OUT(R, e, 1, true, ctx.template zeros_as<R>(nm1));
            ARRAY_OUT(T, tau, 1, true, ctx.zeros(nm1));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));

            lapack::sytrd(lower ? 'L' : 'U', n, a.data<T>(), lda, d.data<R>(), e.data<R>(),
                          tau.data<T>(), work.data<T>(), lwork, &info);
            return make_result(a, d, e, tau, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sytrd_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lower", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "hetrd_lwork" : "sytrd_lwork", "O|O", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            T work = 0;
            CBLAS_INT info = 0;
            lapack::sytrd(lower ? 'L' : 'U', n, nullptr, std::max<CBLAS_INT>(1, n), nullptr,
                          nullptr, nullptr, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sygst(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "itype", "lower", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "hegst" : "sygst", "OO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, itype, 1);  CHECK(itype == 1 || itype == 2 || itype == 3, itype);
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), info = 0;

            ARRAY_IN(T, b, 2);
            CHECKARRAY(shape(b, 0) == n && shape(b, 1) == n, b);
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0));

            lapack::sygst(itype, lower ? 'L' : 'U', n, a.data<T>(), lda, b.data<T>(), ldb, &info);
            return make_result(a, static_cast<long long>(info));
        }


        PyMethodDef sym_herm_methods[] = {
            FAMILY(sytrf),
            FAMILY(sytrf_lwork),
            FAMILY(sytf2),
            FAMILY(sytrs),
            FAMILY(sytri),
            FAMILY(syconv),
            FAMILY(syequb),
            FAMILY(sycon),
            FAMILY(sysv),
            FAMILY(sysv_lwork),
            FAMILY(sysvx),
            FAMILY(sysvx_lwork),
            ROW(chetrf, hetrf, c64),
            ROW(zhetrf, hetrf, c128),
            ROW(chetrf_lwork, hetrf_lwork, c64),
            ROW(zhetrf_lwork, hetrf_lwork, c128),
            ROW(chetrs, hetrs, c64),
            ROW(zhetrs, hetrs, c128),
            ROW(chetri, hetri, c64),
            ROW(zhetri, hetri, c128),
            ROW(cheequb, heequb, c64),
            ROW(zheequb, heequb, c128),
            ROW(checon, hecon, c64),
            ROW(zhecon, hecon, c128),
            ROW(chesv, hesv, c64),
            ROW(zhesv, hesv, c128),
            ROW(chesv_lwork, hesv_lwork, c64),
            ROW(zhesv_lwork, hesv_lwork, c128),
            ROW(chesvx, hesvx, c64),
            ROW(zhesvx, hesvx, c128),
            ROW(chesvx_lwork, hesvx_lwork, c64),
            ROW(zhesvx_lwork, hesvx_lwork, c128),

            /* The eigen half.  A merged `sy`/`he` pair registers all four flavors against the
             * one template; the five split pairs register two apiece. */
            ROW(ssyev, syev, f32),               ROW(dsyev, syev, f64),
            ROW(cheev, syev, c64),               ROW(zheev, syev, c128),
            ROW(ssyev_lwork, syev_lwork, f32),   ROW(dsyev_lwork, syev_lwork, f64),
            ROW(cheev_lwork, syev_lwork, c64),   ROW(zheev_lwork, syev_lwork, c128),
            ROW(ssyevd, syevd, f32),             ROW(dsyevd, syevd, f64),
            ROW(cheevd, heevd, c64),             ROW(zheevd, heevd, c128),
            ROW(ssyevd_lwork, syevd_lwork, f32), ROW(dsyevd_lwork, syevd_lwork, f64),
            ROW(cheevd_lwork, heevd_lwork, c64), ROW(zheevd_lwork, heevd_lwork, c128),
            ROW(ssyevr, syevr, f32),             ROW(dsyevr, syevr, f64),
            ROW(cheevr, heevr, c64),             ROW(zheevr, heevr, c128),
            ROW(ssyevr_lwork, syevr_lwork, f32), ROW(dsyevr_lwork, syevr_lwork, f64),
            ROW(cheevr_lwork, heevr_lwork, c64), ROW(zheevr_lwork, heevr_lwork, c128),
            ROW(ssyevx, syevx, f32),             ROW(dsyevx, syevx, f64),
            ROW(cheevx, syevx, c64),             ROW(zheevx, syevx, c128),
            ROW(ssyevx_lwork, syevx_lwork, f32), ROW(dsyevx_lwork, syevx_lwork, f64),
            ROW(cheevx_lwork, syevx_lwork, c64), ROW(zheevx_lwork, syevx_lwork, c128),
            ROW(ssygv, sygv, f32),               ROW(dsygv, sygv, f64),
            ROW(chegv, sygv, c64),               ROW(zhegv, sygv, c128),
            ROW(ssygv_lwork, sygv_lwork, f32),   ROW(dsygv_lwork, sygv_lwork, f64),
            ROW(chegv_lwork, sygv_lwork, c64),   ROW(zhegv_lwork, sygv_lwork, c128),
            ROW(ssygvd, sygvd, f32),             ROW(dsygvd, sygvd, f64),
            ROW(chegvd, hegvd, c64),             ROW(zhegvd, hegvd, c128),
            ROW(ssygvx, sygvx, f32),             ROW(dsygvx, sygvx, f64),
            ROW(chegvx, sygvx, c64),             ROW(zhegvx, sygvx, c128),
            ROW(ssygvx_lwork, sygvx_lwork, f32), ROW(dsygvx_lwork, sygvx_lwork, f64),
            ROW(chegvx_lwork, sygvx_lwork, c64), ROW(zhegvx_lwork, sygvx_lwork, c128),
            ROW(ssytrd, sytrd, f32),             ROW(dsytrd, sytrd, f64),
            ROW(chetrd, sytrd, c64),             ROW(zhetrd, sytrd, c128),
            ROW(ssytrd_lwork, sytrd_lwork, f32), ROW(dsytrd_lwork, sytrd_lwork, f64),
            ROW(chetrd_lwork, sytrd_lwork, c64), ROW(zhetrd_lwork, sytrd_lwork, c128),
            ROW(ssygst, sygst, f32),             ROW(dsygst, sygst, f64),
            ROW(chegst, sygst, c64),             ROW(zhegst, sygst, c128),
            {nullptr, nullptr, 0, nullptr},
        };

    }  // namespace capi
}  // namespace lapack
