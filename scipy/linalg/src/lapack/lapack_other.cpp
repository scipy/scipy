/**
 * @file
 * @brief Python wrappers for various LAPACK routines.
 */
#define PY_ARRAY_UNIQUE_SYMBOL scipy_lapack_ARRAY_API
#define NO_IMPORT_ARRAY
#include "lapack_helpers.hpp"
#include "lapack_calls.hpp"

namespace lapack {
    namespace capi {


        template <class T>
        static PyObject *tpttf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "ap", "transr", "uplo", nullptr};
            static constexpr Ctx<T> ctx("tpttf", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);       CHECK(n >= 0, n);
            SCALAR_OPT(char, transr, 'N');
            if constexpr (is_complex_v<T>) { CHECK(transr == 'N' || transr == 'C', transr); }
            else                           { CHECK(transr == 'N' || transr == 'T', transr); }
            SCALAR_OPT(char, uplo, 'U');    CHECK(uplo == 'U' || uplo == 'L', uplo);

            ARRAY_IN(T, ap, 1);
            CBLAS_INT nt = len(ap), info = 0;
            CHECKARRAY(nt == 1LL * n * (n + 1) / 2, ap);

            ARRAY_OUT(T, arf, 1, true, ctx.zeros(nt));

            lapack::tpttf(transr, uplo, n, ap.data<T>(), arf.data<T>(), &info);
            return make_result(arf, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *tpttr(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "ap", "uplo", nullptr};
            static constexpr Ctx<T> ctx("tpttr", "OO|O", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);     CHECK(n >= 0, n);
            SCALAR_OPT(char, uplo, 'U');  CHECK(uplo == 'U' || uplo == 'L', uplo);

            ARRAY_IN(T, ap, 1);
            CBLAS_INT nt = len(ap), info = 0;
            CHECKARRAY(nt == 1LL * n * (n + 1) / 2, ap);

            /* `a` is `(n, n)` while `lda` is clamped to 1, so the `n == 0` call gets a leading
             * dimension of 1 over an empty array -- what the `.pyf` computes. */
            ARRAY_OUT(T, a, 2, true, ctx.zeros(n, n));
            CBLAS_INT lda = std::max<CBLAS_INT>(n, 1);

            lapack::tpttr(uplo, n, ap.data<T>(), a.data<T>(), lda, &info);
            return make_result(a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *tfttp(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "arf", "transr", "uplo", nullptr};
            static constexpr Ctx<T> ctx("tfttp", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);       CHECK(n >= 0, n);
            SCALAR_OPT(char, transr, 'N');
            if constexpr (is_complex_v<T>) { CHECK(transr == 'N' || transr == 'C', transr); }
            else                           { CHECK(transr == 'N' || transr == 'T', transr); }
            SCALAR_OPT(char, uplo, 'U');    CHECK(uplo == 'U' || uplo == 'L', uplo);

            ARRAY_IN(T, arf, 1);
            CBLAS_INT nt = len(arf), info = 0;
            CHECKARRAY(nt == 1LL * n * (n + 1) / 2, arf);

            ARRAY_OUT(T, ap, 1, true, ctx.zeros(nt));

            lapack::tfttp(transr, uplo, n, arf.data<T>(), ap.data<T>(), &info);
            return make_result(ap, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *tfttr(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "arf", "transr", "uplo", nullptr};
            static constexpr Ctx<T> ctx("tfttr", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);       CHECK(n >= 0, n);
            SCALAR_OPT(char, transr, 'N');
            if constexpr (is_complex_v<T>) { CHECK(transr == 'N' || transr == 'C', transr); }
            else                           { CHECK(transr == 'N' || transr == 'T', transr); }
            SCALAR_OPT(char, uplo, 'U');    CHECK(uplo == 'U' || uplo == 'L', uplo);

            ARRAY_IN(T, arf, 1);
            CBLAS_INT nt = len(arf), info = 0;
            CHECKARRAY(nt == 1LL * n * (n + 1) / 2, arf);

            CBLAS_INT lda = std::max<CBLAS_INT>(n, 1);
            ARRAY_OUT(T, a, 2, true, ctx.zeros(lda, n));

            lapack::tfttr(transr, uplo, n, arf.data<T>(), a.data<T>(), lda, &info);
            return make_result(a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *trttf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "transr", "uplo", nullptr};
            static constexpr Ctx<T> ctx("trttf", "O|OO", kwlist);
            PARSE_ARGS();

            SCALAR_OPT(char, transr, 'N');
            if constexpr (is_complex_v<T>) { CHECK(transr == 'N' || transr == 'C', transr); }
            else                           { CHECK(transr == 'N' || transr == 'T', transr); }
            SCALAR_OPT(char, uplo, 'U');    CHECK(uplo == 'U' || uplo == 'L', uplo);

            ARRAY_IN(T, a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 1), lda = std::max<CBLAS_INT>(shape(a, 0), 1), info = 0;
            npy_intp nt = 1LL * n * (n + 1) / 2;
            ARRAY_OUT(T, arf, 1, true, ctx.zeros(nt));

            lapack::trttf(transr, uplo, n, a.data<T>(), lda, arf.data<T>(), &info);
            return make_result(arf, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *trttp(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "uplo", nullptr};
            static constexpr Ctx<T> ctx("trttp", "O|O", kwlist);
            PARSE_ARGS();

            SCALAR_OPT(char, uplo, 'U');  CHECK(uplo == 'U' || uplo == 'L', uplo);

            ARRAY_IN(T, a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 1), lda = std::max<CBLAS_INT>(shape(a, 0), 1), info = 0;

            npy_intp nt = 1LL * n * (n + 1) / 2;
            ARRAY_OUT(T, ap, 1, true, ctx.zeros(nt));

            lapack::trttp(uplo, n, a.data<T>(), lda, ap.data<T>(), &info);
            return make_result(ap, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *tfsm(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"alpha", "a", "b", "transr", "side", "uplo", "trans",
                                           "diag", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("tfsm", "OOO|OOOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            SCALAR_REQ(T, alpha);
            SCALAR_OPT(char, transr, 'N');
            if constexpr (is_complex_v<T>) { CHECK(transr == 'N' || transr == 'C', transr); }
            else                           { CHECK(transr == 'N' || transr == 'T', transr); }
            SCALAR_OPT(char, side, 'L');    CHECK(side == 'L' || side == 'R', side);
            SCALAR_OPT(char, uplo, 'U');    CHECK(uplo == 'U' || uplo == 'L', uplo);
            SCALAR_OPT(char, trans, 'N');
            if constexpr (is_complex_v<T>) { CHECK(trans == 'N' || trans == 'C', trans); }
            else                           { CHECK(trans == 'N' || trans == 'T', trans); }
            SCALAR_OPT(char, diag, 'N');    CHECK(diag == 'U' || diag == 'N', diag);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CBLAS_INT m = shape(b, 0), n = shape(b, 1), ldb = std::max<CBLAS_INT>(m, 1);

            ARRAY_IN(T, a, 1);
            CBLAS_INT order = side == 'L' ? m : n;
            CHECKARRAY(len(a) == 1LL * order * (order + 1) / 2, a);

            lapack::tfsm(transr, side, uplo, trans, diag, m, n, alpha, a.data<T>(), b.data<T>(), ldb);
            return make_result(b);
        }


        template <class T>
        static PyObject *ppcon(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "ap", "anorm", "lower", nullptr};
            static constexpr Ctx<T> ctx("ppcon", "OOO|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);
            SCALAR_REQ(CBLAS_INT, n);         CHECK(n >= 0, n);

            ARRAY_IN(T, ap, 1);
            CHECKARRAY(len(ap) >= 1LL * n * (n + 1) / 2, ap);
            SCALAR_REQ(R, anorm);

            CBLAS_INT work_len, info = 0;
            if (!work_size(is_complex_v<T> ? 2LL * n : 3LL * n, &work_len)) { return nullptr; }
            using W = std::conditional_t<is_complex_v<T>, R, CBLAS_INT>;
            ARRAY_HIDDEN(T, work, work_len);
            ARRAY_HIDDEN(W, irwork, n);
            R rcond = 0;

            lapack::ppcon(lower ? 'L' : 'U', n, ap.data<T>(), anorm, &rcond, work.data<T>(), irwork.data<W>(), &info);
            return make_result(rcond, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *ppsv(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "ap", "b", "lower", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("ppsv", "OOO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);
            SCALAR_REQ(CBLAS_INT, n);         CHECK(n >= 0, n);

            ARRAY_IN(T, ap, 1);
            CHECKARRAY(len(ap) >= 1LL * n * (n + 1) / 2, ap);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0)), nrhs = shape(b, 1), info = 0;

            lapack::ppsv(lower ? 'L' : 'U', n, nrhs, ap.data<T>(), b.data<T>(), ldb, &info);
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pptrf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "ap", "lower", "overwrite_ap", nullptr};
            static constexpr Ctx<T> ctx("pptrf", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_ap);
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);
            SCALAR_REQ(CBLAS_INT, n);         CHECK(n >= 0, n);

            ARRAY_INOUT(T, ap, 1, overwrite_ap != 0);
            CHECKARRAY(len(ap) >= 1LL * n * (n + 1) / 2, ap);
            CBLAS_INT info = 0;

            lapack::pptrf(lower ? 'L' : 'U', n, ap.data<T>(), &info);
            return make_result(ap, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pptri(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "ap", "lower", "overwrite_ap", nullptr};
            static constexpr Ctx<T> ctx("pptri", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_ap);
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);
            SCALAR_REQ(CBLAS_INT, n);         CHECK(n >= 0, n);

            ARRAY_INOUT(T, ap, 1, overwrite_ap != 0);
            CHECKARRAY(len(ap) >= 1LL * n * (n + 1) / 2, ap);
            CBLAS_INT info = 0;

            lapack::pptri(lower ? 'L' : 'U', n, ap.data<T>(), &info);
            return make_result(ap, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pptrs(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "ap", "b", "lower", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("pptrs", "OOO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);
            SCALAR_REQ(CBLAS_INT, n);         CHECK(n >= 0, n);

            ARRAY_IN(T, ap, 1);
            CHECKARRAY(len(ap) >= 1LL * n * (n + 1) / 2, ap);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0)), nrhs = shape(b, 1), info = 0;

            lapack::pptrs(lower ? 'L' : 'U', n, nrhs, ap.data<T>(), b.data<T>(), ldb, &info);
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sfrk(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "k", "alpha", "a", "beta", "c", "transr", "uplo",
                                           "trans", "overwrite_c", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "hfrk" : "sfrk", "OOOOOO|OOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_c);
            SCALAR_OPT(char, transr, 'N');
            if constexpr (is_complex_v<T>) { CHECK(transr == 'N' || transr == 'C', transr); }
            else                           { CHECK(transr == 'N' || transr == 'T', transr); }
            SCALAR_OPT(char, uplo, 'U');   CHECK(uplo == 'U' || uplo == 'L', uplo);
            SCALAR_OPT(char, trans, 'N');
            if constexpr (is_complex_v<T>) { CHECK(trans == 'N' || trans == 'C', trans); }
            else                           { CHECK(trans == 'N' || trans == 'T', trans); }
            SCALAR_REQ(CBLAS_INT, n);      CHECK(n >= 0, n);
            SCALAR_REQ(CBLAS_INT, k);      CHECK(k >= 0, k);
            SCALAR_REQ(R, alpha);
            SCALAR_REQ(R, beta);

            ARRAY_IN(T, a, 2);
            CBLAS_INT lda = std::max<CBLAS_INT>(trans == 'N' ? n : k, 1);
            CHECKARRAY(shape(a, 1) == (trans == 'N' ? k : n), a);
            CHECKARRAY(shape(a, 0) == lda, a);

            ARRAY_INOUT(T, c, 1, overwrite_c != 0);
            CHECKARRAY(len(c) == 1LL * n * (n + 1) / 2, c);

            lapack::sfrk(transr, uplo, trans, n, k, alpha, a.data<T>(), lda, beta, c.data<T>());
            return make_result(c);
        }


        template <class T>
        static PyObject *pftrf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "a", "transr", "uplo", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("pftrf", "OO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(char, transr, 'N');
            if constexpr (is_complex_v<T>) { CHECK(transr == 'N' || transr == 'C', transr); }
            else                           { CHECK(transr == 'N' || transr == 'T', transr); }
            SCALAR_OPT(char, uplo, 'U');  CHECK(uplo == 'U' || uplo == 'L', uplo);
            SCALAR_REQ(CBLAS_INT, n);     CHECK(n >= 0, n);

            ARRAY_INOUT(T, a, 1, overwrite_a != 0);
            CHECKARRAY(len(a) == 1LL * n * (n + 1) / 2, a);
            CBLAS_INT info = 0;

            lapack::pftrf(transr, uplo, n, a.data<T>(), &info);
            return make_result(a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pftri(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "a", "transr", "uplo", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("pftri", "OO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(char, transr, 'N');
            if constexpr (is_complex_v<T>) { CHECK(transr == 'N' || transr == 'C', transr); }
            else                           { CHECK(transr == 'N' || transr == 'T', transr); }
            SCALAR_OPT(char, uplo, 'U');  CHECK(uplo == 'U' || uplo == 'L', uplo);
            SCALAR_REQ(CBLAS_INT, n);     CHECK(n >= 0, n);

            ARRAY_INOUT(T, a, 1, overwrite_a != 0);
            CHECKARRAY(len(a) == 1LL * n * (n + 1) / 2, a);
            CBLAS_INT info = 0;

            lapack::pftri(transr, uplo, n, a.data<T>(), &info);
            return make_result(a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pftrs(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "a", "b", "transr", "uplo", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("pftrs", "OOO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(char, transr, 'N');
            if constexpr (is_complex_v<T>) { CHECK(transr == 'N' || transr == 'C', transr); }
            else                           { CHECK(transr == 'N' || transr == 'T', transr); }
            SCALAR_OPT(char, uplo, 'U');  CHECK(uplo == 'U' || uplo == 'L', uplo);
            SCALAR_REQ(CBLAS_INT, n);     CHECK(n >= 0, n);

            ARRAY_IN(T, a, 1);
            CHECKARRAY(len(a) == 1LL * n * (n + 1) / 2, a);

            /* `b` may have more rows than the system needs; only the leading `n` are read. */
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) >= n, b);
            CBLAS_INT nrhs = shape(b, 1), ldb = std::max<CBLAS_INT>(shape(b, 0), 1), info = 0;

            lapack::pftrs(transr, uplo, n, nrhs, a.data<T>(), b.data<T>(), ldb, &info);
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pbtrf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"ab", "lower", "ldab", "overwrite_ab", nullptr};
            static constexpr Ctx<T> ctx("pbtrf", "O|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_ab);
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, ab, 2, overwrite_ab != 0);
            SCALAR_OPT(CBLAS_INT, ldab, shape(ab, 0));
            CHECK(shape(ab, 0) == ldab, ldab);

            CBLAS_INT n = shape(ab, 1), kd = shape(ab, 0) - 1, info = 0;

            lapack::pbtrf(lower ? 'L' : 'U', n, kd, ab.data<T>(), ldab, &info);
            return make_result(ab, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pbtrs(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"ab", "b", "lower", "ldab", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("pbtrs", "OO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);

            ARRAY_IN(T, ab, 2);
            SCALAR_OPT(CBLAS_INT, ldab, shape(ab, 0));
            CHECK(shape(ab, 0) == ldab, ldab);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CBLAS_INT n = shape(ab, 1), kd = shape(ab, 0) - 1;
            CBLAS_INT ldb = shape(b, 0), nrhs = shape(b, 1), info = 0;

            lapack::pbtrs(lower ? 'L' : 'U', n, kd, nrhs, ab.data<T>(), ldab, b.data<T>(), ldb, &info);
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pbsv(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"ab", "b", "lower", "ldab", "overwrite_ab",
                                           "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("pbsv", "OO|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_ab);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, ab, 2, overwrite_ab != 0);
            SCALAR_OPT(CBLAS_INT, ldab, shape(ab, 0));
            CHECK(shape(ab, 0) == ldab, ldab);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CBLAS_INT n = shape(ab, 1), kd = shape(ab, 0) - 1;
            CBLAS_INT ldb = shape(b, 0), nrhs = shape(b, 1), info = 0;

            lapack::pbsv(lower ? 'L' : 'U', n, kd, nrhs, ab.data<T>(), ldab, b.data<T>(), ldb, &info);
            return make_result(ab, b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *trtrs(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "lower", "trans", "unitdiag", "lda",
                                           "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("trtrs", "OO|OOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, lower, 0);     CHECK(lower == 0 || lower == 1, lower);
            SCALAR_OPT(CBLAS_INT, trans, 0);     CHECK(trans >= 0 && trans <= 2, trans);
            SCALAR_OPT(CBLAS_INT, unitdiag, 0);  CHECK(unitdiag == 0 || unitdiag == 1, unitdiag);

            ARRAY_IN(T, a, 2);
            SCALAR_OPT(CBLAS_INT, lda, shape(a, 0));
            CHECK(shape(a, 0) == lda, lda);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CBLAS_INT n = shape(a, 1), ldb = shape(b, 0), nrhs = shape(b, 1), info = 0;

            /* LAPACK requires both leading dimensions to be at least `MAX(1, n)`, so an empty
             * matrix -- whose shape gives 0 -- has to be floored on the way in.  The `.pyf`
             * passed the raw shapes and an empty triangular system came back `info == -7`
             * instead of trivially solved.  `lda` itself keeps the caller's value: it is a
             * checked argument and the check is against `a`, not against this floor. */
            lapack::trtrs(lower ? 'L' : 'U', trans ? (trans == 2 ? 'C' : 'T') : 'N',
                          unitdiag ? 'U' : 'N', n, nrhs, a.data<T>(),
                          std::max<CBLAS_INT>(1, lda), b.data<T>(),
                          std::max<CBLAS_INT>(1, ldb), &info);
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *trcon(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "norm", "uplo", "diag", nullptr};
            static constexpr Ctx<T> ctx("trcon", "O|OOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(char, norm, '1');  CHECK(norm == '1' || norm == 'I' || norm == 'O', norm);
            SCALAR_OPT(char, uplo, 'U');  CHECK(uplo == 'U' || uplo == 'L', uplo);
            SCALAR_OPT(char, diag, 'N');  CHECK(diag == 'N' || diag == 'U', diag);

            ARRAY_IN(T, a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 1), lda = std::max<CBLAS_INT>(1, n), work_len, info = 0;
            if (!work_size(3LL * n, &work_len)) { return nullptr; }
            using W = std::conditional_t<is_complex_v<T>, R, CBLAS_INT>;
            ARRAY_HIDDEN(T, work, work_len);
            ARRAY_HIDDEN(W, irwork, n);
            R rcond = 0;

            lapack::trcon(norm, uplo, diag, n, a.data<T>(), lda, &rcond, work.data<T>(), irwork.data<W>(), &info);
            return make_result(rcond, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *tbtrs(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"ab", "b", "uplo", "trans", "diag", "overwrite_b",
                                           nullptr};
            static constexpr Ctx<T> ctx("tbtrs", "OO|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(char, uplo, 'U');   CHECK(uplo == 'U' || uplo == 'L', uplo);
            /* Unlike the RFP routines, `tbtrs` takes all three letters in every flavor. */
            SCALAR_OPT(char, trans, 'N');  CHECK(trans == 'N' || trans == 'T' || trans == 'C', trans);
            SCALAR_OPT(char, diag, 'N');   CHECK(diag == 'N' || diag == 'U', diag);

            ARRAY_IN(T, ab, 2);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);

            CBLAS_INT ldab = std::max<CBLAS_INT>(1, shape(ab, 0));
            CBLAS_INT n = std::max<CBLAS_INT>(1, shape(ab, 1));
            CBLAS_INT kd = ldab - 1;
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0));
            CBLAS_INT nrhs = std::max<CBLAS_INT>(1, shape(b, 1)), info = 0;
            CHECK(ldb >= n, ldb);

            lapack::tbtrs(uplo, trans, diag, n, kd, nrhs, ab.data<T>(), ldab, b.data<T>(), ldb, &info);
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *trtri(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"c", "lower", "unitdiag", "overwrite_c", nullptr};
            static constexpr Ctx<T> ctx("trtri", "O|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);
            SCALAR_OPT(CBLAS_INT, lower, 0);     CHECK(lower == 0 || lower == 1, lower);
            SCALAR_OPT(CBLAS_INT, unitdiag, 0);  CHECK(unitdiag == 0 || unitdiag == 1, unitdiag);

            ARRAY_INOUT(T, c, 2, overwrite_c != 0);
            CHECKARRAY(shape(c, 0) == shape(c, 1), c);
            CBLAS_INT n = shape(c, 0), info = 0;

            /* The `.pyf` passed `n` as the leading dimension too, so a zero-order matrix reached
             * LAPACK with `lda == 0` and was refused with `info == -5`.  LAPACK's own bound is
             * `MAX(1, n)`, so the floor goes here and inverting an empty triangle succeeds. */
            lapack::trtri(lower ? 'L' : 'U', unitdiag ? 'U' : 'N', n, c.data<T>(),
                          std::max<CBLAS_INT>(1, n), &info);
            return make_result(c, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *lauum(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"c", "lower", "overwrite_c", nullptr};
            static constexpr Ctx<T> ctx("lauum", "O|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, c, 2, overwrite_c != 0);
            CHECKARRAY(shape(c, 0) == shape(c, 1), c);
            CBLAS_INT n = shape(c, 0), info = 0;

            /* Same zero-order floor as `trtri`; the `.pyf` refused an empty triangle with
             * `info == -4`. */
            lapack::lauum(lower ? 'L' : 'U', n, c.data<T>(), std::max<CBLAS_INT>(1, n), &info);
            return make_result(c, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *laswp(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "piv", "k1", "k2", "off", "inc", "overwrite_a",
                                           nullptr};
            static constexpr Ctx<T> ctx("laswp", "OO|OOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            ARRAY_IN(CBLAS_INT, piv, 1);

            CBLAS_INT nrows = shape(a, 0), n = shape(a, 1), npiv = len(piv);
            CHECKARRAY(npiv <= nrows, piv);

            SCALAR_OPT(CBLAS_INT, k1, 0);         CHECK(0 <= k1, k1);
            SCALAR_OPT(CBLAS_INT, inc, 1);        CHECK(inc > 0 || inc < 0, inc);
            SCALAR_OPT(CBLAS_INT, off, 0);        CHECK(off >= 0 && off < npiv, off);
            SCALAR_OPT(CBLAS_INT, k2, npiv - 1);  CHECK(k1 <= k2 && k2 < npiv - off, k2);

            npy_intp step = abs(inc);
            CHECK(1LL * k1 + 1 + 1LL * (k2 - k1) * step <= npiv - off, inc);
            ARRAY_HIDDEN(CBLAS_INT, pivots, npiv);
            const CBLAS_INT *supplied = piv.data<CBLAS_INT>();
            CBLAS_INT *shifted = pivots.data<CBLAS_INT>();
            for (CBLAS_INT i = 0; i < npiv; ++i) { shifted[i] = supplied[i] + 1; }

            lapack::laswp(n, a.data<T>(), nrows, k1 + 1, k2 + 1, shifted + off, inc);
            return make_result(a);
        }


        template <class T>
        static PyObject *trexc(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "q", "ifst", "ilst", "wantq", "overwrite_a",
                                           "overwrite_q", nullptr};
            static constexpr Ctx<T> ctx("trexc", "OOOO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_q);
            SCALAR_OPT(CBLAS_INT, wantq, 1);  CHECK(wantq == 0 || wantq == 1, wantq);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            ARRAY_INOUT(T, q, 2, overwrite_q != 0);
            SCALAR_REQ(CBLAS_INT, ifst);
            SCALAR_REQ(CBLAS_INT, ilst);

            CBLAS_INT n = shape(a, 1), info = 0;
            CBLAS_INT lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CBLAS_INT ldq = std::max<CBLAS_INT>(1, shape(q, 0));

            if constexpr (is_complex_v<T>) {
                lapack::trexc(wantq ? 'V' : 'N', n, a.data<T>(), lda, q.data<T>(), ldq, &ifst, &ilst, &info);
            }
            else {
                ARRAY_HIDDEN(T, work, n);
                lapack::trexc(wantq ? 'V' : 'N', n, a.data<T>(), lda, q.data<T>(), ldq, &ifst, &ilst, work.data<T>(), &info);
            }
            return make_result(a, q, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *tgexc(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            /* The complex flavors have no workspace, so they take no `lwork` either. */
            static const char *kwlist_real[] = {"a", "b", "q", "z", "ifst", "ilst", "wantq",
                                                "wantz", "lwork", "overwrite_a", "overwrite_b",
                                                "overwrite_q", "overwrite_z", nullptr};
            static const char *kwlist_complex[] = {"a", "b", "q", "z", "ifst", "ilst", "wantq",
                                                   "wantz", "overwrite_a", "overwrite_b",
                                                   "overwrite_q", "overwrite_z", nullptr};
            static constexpr Ctx<T> ctx("tgexc",
                                        is_complex_v<T> ? "OOOOOO|OOOOOO" : "OOOOOO|OOOOOOO",
                                        is_complex_v<T> ? kwlist_complex : kwlist_real);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_FLAG(overwrite_q);
            SCALAR_FLAG(overwrite_z);
            /* Fortran LOGICALs here, not the option letter `trexc` takes for the same choice. */
            SCALAR_OPT(CBLAS_INT, wantq, 1);  CHECK(wantq == 0 || wantq == 1, wantq);
            SCALAR_OPT(CBLAS_INT, wantz, 1);  CHECK(wantz == 0 || wantz == 1, wantz);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            ARRAY_INOUT(T, q, 2, overwrite_q != 0);
            ARRAY_INOUT(T, z, 2, overwrite_z != 0);
            SCALAR_REQ(CBLAS_INT, ifst);
            SCALAR_REQ(CBLAS_INT, ilst);

            CBLAS_INT n = shape(a, 1), info = 0;
            CBLAS_INT lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0));
            CBLAS_INT ldq = std::max<CBLAS_INT>(1, shape(q, 0));
            CBLAS_INT ldz = std::max<CBLAS_INT>(1, shape(z, 0));

            if constexpr (is_complex_v<T>) {
                lapack::tgexc(wantq, wantz, n, a.data<T>(), lda, b.data<T>(), ldb, q.data<T>(),
                              ldq, z.data<T>(), ldz, &ifst, &ilst, &info);
                return make_result(a, b, q, z, static_cast<long long>(info));
            }
            else {
                CBLAS_INT deflt;
                if (!work_size(4LL * n + 16, &deflt)) { return nullptr; }
                SCALAR_OPT(CBLAS_INT, lwork, deflt);
                CHECK(lwork == -1 || lwork >= 4 * n + 16, lwork);

                ARRAY_OUT(T, work, 1, true, ctx.zeros(std::max<CBLAS_INT>(lwork, 1)));
                lapack::tgexc(wantq, wantz, n, a.data<T>(), lda, b.data<T>(), ldb, q.data<T>(),
                              ldq, z.data<T>(), ldz, &ifst, &ilst, work.data<T>(), lwork, &info);
                return make_result(a, b, q, z, work, static_cast<long long>(info));
            }
        }


        template <class T>
        static PyObject *trsyl(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "c", "trana", "tranb", "isgn",
                                           "overwrite_c", nullptr};
            static constexpr Ctx<T> ctx("trsyl", "OOO|OOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_c);
            SCALAR_OPT(char, trana, 'N');
            if constexpr (is_complex_v<T>) { CHECK(trana == 'N' || trana == 'C', trana); }
            else { CHECK(trana == 'N' || trana == 'T' || trana == 'C', trana); }
            SCALAR_OPT(char, tranb, 'N');
            if constexpr (is_complex_v<T>) { CHECK(tranb == 'N' || tranb == 'C', tranb); }
            else { CHECK(tranb == 'N' || tranb == 'T' || tranb == 'C', tranb); }
            SCALAR_OPT(CBLAS_INT, isgn, 1);  CHECK(isgn == 1 || isgn == -1, isgn);

            ARRAY_IN(T, a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            ARRAY_IN(T, b, 2);
            CHECKARRAY(shape(b, 0) == shape(b, 1), b);
            ARRAY_INOUT(T, c, 2, overwrite_c != 0);

            CBLAS_INT m = shape(a, 0), n = shape(b, 0), info = 0;
            CBLAS_INT lda = shape(a, 0), ldb = shape(b, 0), ldc = shape(c, 0);
            R scale = 0;

            lapack::trsyl(trana, tranb, isgn, m, n, a.data<T>(), lda, b.data<T>(), ldb,
                          c.data<T>(), ldc, &scale, &info);
            return make_result(c, scale, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *trsen(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            /* The complex flavors have no integer workspace, so they take no `liwork`. */
            static const char *kwlist_real[] = {"select", "t", "q", "job", "wantq", "lwork",
                                                "liwork", "overwrite_t", "overwrite_q", nullptr};
            static const char *kwlist_complex[] = {"select", "t", "q", "job", "wantq", "lwork",
                                                   "overwrite_t", "overwrite_q", nullptr};
            static constexpr Ctx<T> ctx("trsen",
                                        is_complex_v<T> ? "OOO|OOOOO" : "OOO|OOOOOO",
                                        is_complex_v<T> ? kwlist_complex : kwlist_real);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_t);
            SCALAR_FLAG(overwrite_q);
            SCALAR_OPT(char, job, 'B');
            CHECK(job == 'N' || job == 'E' || job == 'V' || job == 'B', job);
            SCALAR_OPT(CBLAS_INT, wantq, 1);  CHECK(wantq == 0 || wantq == 1, wantq);

            /* A Fortran LOGICAL array, so `CBLAS_INT` -- a bool array converts on the way in. */
            ARRAY_IN(CBLAS_INT, select, 1);
            ARRAY_INOUT(T, t, 2, overwrite_t != 0);
            ARRAY_INOUT(T, q, 2, overwrite_q != 0);

            CBLAS_INT n = shape(t, 0), info = 0, m = 0;
            CBLAS_INT ldt = std::max<CBLAS_INT>(1, shape(t, 0));
            CBLAS_INT ldq = std::max<CBLAS_INT>(1, shape(q, 0));
            CHECKARRAY(len(select) == n, select);

            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(1, n));
            CHECK(lwork == -1 || lwork >= 1, lwork);
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(lwork, 1));
            R s_out = 0, sep = 0;

            if constexpr (is_complex_v<T>) {
                ARRAY_OUT(T, w, 1, true, ctx.zeros(n));
                lapack::trsen(job, wantq ? 'V' : 'N', select.data<CBLAS_INT>(), n, t.data<T>(),
                              ldt, q.data<T>(), ldq, w.data<T>(), &m, &s_out, &sep,
                              work.data<T>(), lwork, &info);
                return make_result(t, q, w, static_cast<long long>(m), s_out, sep,
                                   static_cast<long long>(info));
            }
            else {
                SCALAR_OPT(CBLAS_INT, liwork, 1);
                CHECK(liwork == -1 || liwork >= 1, liwork);
                ARRAY_HIDDEN(CBLAS_INT, iwork, std::max<CBLAS_INT>(1, liwork));
                ARRAY_OUT(T, wr, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, wi, 1, true, ctx.zeros(n));
                lapack::trsen(job, wantq ? 'V' : 'N', select.data<CBLAS_INT>(), n, t.data<T>(),
                              ldt, q.data<T>(), ldq, wr.data<T>(), wi.data<T>(), &m, &s_out, &sep,
                              work.data<T>(), lwork, iwork.data<CBLAS_INT>(), liwork, &info);
                return make_result(t, q, wr, wi, static_cast<long long>(m), s_out, sep,
                                   static_cast<long long>(info));
            }
        }


        template <class T>
        static PyObject *trsen_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"select", "t", "job", nullptr};
            static constexpr Ctx<T> ctx("trsen_lwork", "OO|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(char, job, 'B');
            CHECK(job == 'N' || job == 'E' || job == 'V' || job == 'B', job);

            ARRAY_IN(CBLAS_INT, select, 1);
            ARRAY_IN(T, t, 2);
            CBLAS_INT n = shape(t, 0), ld = std::max<CBLAS_INT>(1, n), info = 0, m = 0;
            CHECKARRAY(len(select) == n, select);

            T work{}, q{};
            R s_out = 0, sep = 0;

            if constexpr (is_complex_v<T>) {
                T w{};
                lapack::trsen(job, 'N', select.data<CBLAS_INT>(), n, t.data<T>(), ld, &q, ld,
                              &w, &m, &s_out, &sep, &work, -1, &info);
                return make_result(work, static_cast<long long>(info));
            }
            else {
                T wr{}, wi{};
                CBLAS_INT iwork = 0;
                lapack::trsen(job, 'N', select.data<CBLAS_INT>(), n, t.data<T>(), ld, &q, ld,
                              &wr, &wi, &m, &s_out, &sep, &work, -1, &iwork, -1, &info);
                return make_result(work, static_cast<long long>(iwork),
                                   static_cast<long long>(info));
            }
        }


        template <class T>
        static PyObject *tgsen(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"select", "a", "b", "q", "z", "ijob", "wantq", "wantz",
                                           "lwork", "liwork", "overwrite_a", "overwrite_b",
                                           "overwrite_q", "overwrite_z", nullptr};
            static constexpr Ctx<T> ctx("tgsen", "OOOOO|OOOOOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_FLAG(overwrite_q);
            SCALAR_FLAG(overwrite_z);
            SCALAR_OPT(CBLAS_INT, ijob, 4);   CHECK(ijob >= 0 && ijob <= 5, ijob);
            SCALAR_OPT(CBLAS_INT, wantq, 1);  CHECK(wantq == 0 || wantq == 1, wantq);
            SCALAR_OPT(CBLAS_INT, wantz, 1);  CHECK(wantz == 0 || wantz == 1, wantz);

            ARRAY_IN(CBLAS_INT, select, 1);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            ARRAY_INOUT(T, q, 2, overwrite_q != 0);
            ARRAY_INOUT(T, z, 2, overwrite_z != 0);

            CBLAS_INT n = shape(a, 0), info = 0, m = 0;
            CBLAS_INT lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0));
            CBLAS_INT ldq = std::max<CBLAS_INT>(1, shape(q, 0));
            CBLAS_INT ldz = std::max<CBLAS_INT>(1, shape(z, 0));
            CHECKARRAY(len(select) == n, select);

            CBLAS_INT lwork_deflt, liwork_deflt;
            if constexpr (is_complex_v<T>) {
                lwork_deflt = liwork_deflt = ijob == 0 ? 1 : n + 2;
            }
            else {
                if (!work_size(4LL * n + 16, &lwork_deflt)) { return nullptr; }
                if (!work_size(1LL * n + 6, &liwork_deflt)) { return nullptr; }
            }
            SCALAR_OPT(CBLAS_INT, lwork, lwork_deflt);
            CHECK(lwork == -1 || lwork >= 1, lwork);
            SCALAR_OPT(CBLAS_INT, liwork, liwork_deflt);
            CHECK(liwork == -1 || liwork >= 1, liwork);

            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(lwork, 1));
            ARRAY_HIDDEN(CBLAS_INT, iwork,
                         is_complex_v<T> ? liwork : std::max<CBLAS_INT>(1, liwork));
            ARRAY_OUT(R, dif, 1, true, ctx.template zeros_as<R>(2));
            R pl = 0, pr = 0;

            if constexpr (is_complex_v<T>) {
                ARRAY_OUT(T, alpha, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, beta, 1, true, ctx.zeros(n));
                lapack::tgsen(ijob, wantq, wantz, select.data<CBLAS_INT>(), n, a.data<T>(), lda,
                              b.data<T>(), ldb, alpha.data<T>(), beta.data<T>(), q.data<T>(), ldq,
                              z.data<T>(), ldz, &m, &pl, &pr, dif.data<R>(), work.data<T>(),
                              lwork, iwork.data<CBLAS_INT>(), liwork, &info);
                return make_result(a, b, alpha, beta, q, z, static_cast<long long>(m), pl, pr,
                                   dif, static_cast<long long>(info));
            }
            else {
                ARRAY_OUT(T, alphar, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, alphai, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, beta, 1, true, ctx.zeros(n));
                lapack::tgsen(ijob, wantq, wantz, select.data<CBLAS_INT>(), n, a.data<T>(), lda,
                              b.data<T>(), ldb, alphar.data<T>(), alphai.data<T>(),
                              beta.data<T>(), q.data<T>(), ldq, z.data<T>(), ldz, &m, &pl, &pr,
                              dif.data<R>(), work.data<T>(), lwork, iwork.data<CBLAS_INT>(),
                              liwork, &info);
                return make_result(a, b, alphar, alphai, beta, q, z, static_cast<long long>(m),
                                   pl, pr, dif, static_cast<long long>(info));
            }
        }


        template <class T>
        static PyObject *tgsen_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist_real[] = {"select", "a", "ijob", nullptr};
            static const char *kwlist_complex[] = {"select", "a", "b", "ijob", nullptr};
            static constexpr Ctx<T> ctx("tgsen_lwork",
                                        is_complex_v<T> ? "OOO|O" : "OO|O",
                                        is_complex_v<T> ? kwlist_complex : kwlist_real);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(CBLAS_INT, ijob, 4);  CHECK(ijob >= 0 && ijob <= 5, ijob);

            ARRAY_IN(CBLAS_INT, select, 1);
            ARRAY_IN(T, a, 2);
            CBLAS_INT n = shape(a, 0), ld = std::max<CBLAS_INT>(1, n), info = 0, m = 0;
            CHECKARRAY(len(select) == n, select);

            T work{}, q{}, z{};
            R pl = 0, pr = 0, dif = 0;
            CBLAS_INT iwork = 0;

            if constexpr (is_complex_v<T>) {
                ARRAY_IN(T, b, 2);
                /* Unlike the real query, this one hands over real `alpha`/`beta` arrays. */
                ARRAY_HIDDEN(T, alpha, n);
                ARRAY_HIDDEN(T, beta, n);
                lapack::tgsen(ijob, 0, 0, select.data<CBLAS_INT>(), n, a.data<T>(), ld,
                              b.data<T>(), ld, alpha.data<T>(), beta.data<T>(), &q, ld, &z, ld,
                              &m, &pl, &pr, &dif, &work, -1, &iwork, -1, &info);
            }
            else {
                T b{}, alphar{}, alphai{}, beta{};
                lapack::tgsen(ijob, 0, 0, select.data<CBLAS_INT>(), n, a.data<T>(), ld, &b, ld,
                              &alphar, &alphai, &beta, &q, ld, &z, ld, &m, &pl, &pr, &dif,
                              &work, -1, &iwork, -1, &info);
            }
            return make_result(work, static_cast<long long>(iwork), static_cast<long long>(info));
        }


        template <class T>
        static PyObject *orghr(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "tau", "lo", "hi", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "unghr" : "orghr", "OO|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), info = 0;

            ARRAY_IN(T, tau, 1);
            CHECKARRAY(len(tau) == (n > 0 ? n - 1 : 0), tau);

            SCALAR_OPT(CBLAS_INT, lo, 0);
            SCALAR_OPT(CBLAS_INT, hi, n - 1);
            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(hi - lo, 1));
            CHECK(lwork >= hi - lo, lwork);

            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(lwork, 1));
            lapack::orghr(n, lo + 1, hi + 1, a.data<T>(), n, tau.data<T>(), work.data<T>(),
                          lwork, &info);
            return make_result(a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *orghr_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lo", "hi", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "unghr_lwork" : "orghr_lwork",
                                        "O|OO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lo, 0);
            SCALAR_OPT(CBLAS_INT, hi, n - 1);

            T work{}, a{}, tau{};
            CBLAS_INT info = 0;
            lapack::orghr(n, lo + 1, hi + 1, &a, n, &tau, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *orgqr(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "tau", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "ungqr" : "orgqr", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            ARRAY_IN(T, tau, 1);

            CBLAS_INT m = shape(a, 0), n = shape(a, 1), k = len(tau), info = 0;
            CBLAS_INT deflt;
            if (!work_size(3LL * n, &deflt)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, deflt);
            CHECK(lwork >= n || lwork == -1, lwork);

            ARRAY_OUT(T, work, 1, true, ctx.zeros(std::max<CBLAS_INT>(lwork, 1)));
            lapack::orgqr(m, n, k, a.data<T>(), m, tau.data<T>(), work.data<T>(), lwork, &info);
            return make_result(a, work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *orgrq(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "tau", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "ungrq" : "orgrq", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            ARRAY_IN(T, tau, 1);

            CBLAS_INT m = shape(a, 0), n = shape(a, 1), k = len(tau), info = 0;
            CBLAS_INT deflt;
            if (!work_size(3LL * m, &deflt)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, deflt);
            CHECK(lwork >= m || lwork == -1, lwork);

            ARRAY_OUT(T, work, 1, true, ctx.zeros(std::max<CBLAS_INT>(lwork, 1)));
            lapack::orgrq(m, n, k, a.data<T>(), m, tau.data<T>(), work.data<T>(), lwork, &info);
            return make_result(a, work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *ormqr(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"side", "trans", "a", "tau", "c", "lwork",
                                           "overwrite_c", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "unmqr" : "ormqr", "OOOOOO|O", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);
            SCALAR_REQ(char, side);   CHECK(side == 'L' || side == 'R', side);
            SCALAR_REQ(char, trans);
            if constexpr (is_complex_v<T>) { CHECK(trans == 'N' || trans == 'C', trans); }
            else                           { CHECK(trans == 'N' || trans == 'T', trans); }

            ARRAY_IN(T, a, 2);
            ARRAY_IN(T, tau, 1);
            ARRAY_INOUT(T, c, 2, overwrite_c != 0);
            SCALAR_REQ(CBLAS_INT, lwork);

            CBLAS_INT m = shape(c, 0), n = shape(c, 1), k = shape(a, 1), info = 0;
            CBLAS_INT lda = shape(a, 0), ldc = shape(c, 0);
            CHECKARRAY(len(tau) == k, tau);

            ARRAY_OUT(T, work, 1, true, ctx.zeros(std::max<CBLAS_INT>(lwork, 1)));
            lapack::ormqr(side, trans, m, n, k, a.data<T>(), lda, tau.data<T>(), c.data<T>(),
                          ldc, work.data<T>(), lwork, &info);
            return make_result(c, work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *ormrz(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "tau", "c", "side", "trans", "lwork",
                                           "overwrite_c", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "unmrz" : "ormrz", "OOO|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);
            SCALAR_OPT(char, side, 'L');   CHECK(side == 'L' || side == 'R', side);
            SCALAR_OPT(char, trans, 'N');
            if constexpr (is_complex_v<T>) { CHECK(trans == 'N' || trans == 'C', trans); }
            else                           { CHECK(trans == 'N' || trans == 'T', trans); }

            /* `a` holds the `k`-by-`nt` factor `tzrzf` produced: `k` reflectors over a row of
             * length `nt`, of which the trailing `l = nt - k` columns carry the Z part. */
            ARRAY_IN(T, a, 2);
            CHECKARRAY(shape(a, 1) >= shape(a, 0), a);
            ARRAY_IN(T, tau, 1);
            ARRAY_INOUT(T, c, 2, overwrite_c != 0);

            CBLAS_INT m = shape(c, 0), n = shape(c, 1), info = 0;
            CBLAS_INT k = shape(a, 0), nt = shape(a, 1), l = nt - k;
            CBLAS_INT lda = std::max<CBLAS_INT>(shape(a, 0), 1);
            CBLAS_INT ldc = std::max<CBLAS_INT>(shape(c, 0), 1);
            CHECKARRAY((side == 'L' ? m : n) == nt, a);
            CHECKARRAY(len(tau) == k, tau);

            CBLAS_INT order = side == 'L' ? n : m;
            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(order, 1));
            CHECK(lwork >= order || lwork == -1, lwork);

            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(lwork, 1));
            lapack::ormrz(side, trans, m, n, k, l, a.data<T>(), lda, tau.data<T>(), c.data<T>(),
                          ldc, work.data<T>(), lwork, &info);
            return make_result(c, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *ormrz_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"m", "n", "side", "trans", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "unmrz_lwork" : "ormrz_lwork",
                                        "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_OPT(char, side, 'L');   CHECK(side == 'L' || side == 'R', side);
            SCALAR_OPT(char, trans, 'N');
            if constexpr (is_complex_v<T>) { CHECK(trans == 'N' || trans == 'C', trans); }
            else                           { CHECK(trans == 'N' || trans == 'T', trans); }
            SCALAR_REQ(CBLAS_INT, m);
            SCALAR_REQ(CBLAS_INT, n);

            CBLAS_INT k = side == 'L' ? m : n, l = 0, info = 0;
            T work{}, a{}, tau{}, c{};
            lapack::ormrz(side, trans, m, n, k, l, &a, k, &tau, &c, m, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *geqrt(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"nb", "a", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("geqrt", "OO|O", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_REQ(CBLAS_INT, nb);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);

            CBLAS_INT m = shape(a, 0), n = shape(a, 1), info = 0;
            CBLAS_INT lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CHECK(std::min(m, n) >= nb && nb >= 1, nb);

            ARRAY_OUT(T, t, 2, true, ctx.zeros(nb, std::min(m, n)));
            CBLAS_INT ldt = std::max<CBLAS_INT>(1, shape(t, 0));
            ARRAY_HIDDEN(T, work, nb, n);

            lapack::geqrt(m, n, nb, a.data<T>(), lda, t.data<T>(), ldt, work.data<T>(), &info);
            return make_result(a, t, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gemqrt(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            /* `side` and `trans` are optional here, so they trail the operands rather than
             * leading them as they do in `ormqr`. */
            static const char *kwlist[] = {"v", "t", "c", "side", "trans", "overwrite_c", nullptr};
            static constexpr Ctx<T> ctx("gemqrt", "OOO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);
            SCALAR_OPT(char, side, 'L');   CHECK(side == 'L' || side == 'R', side);
            SCALAR_OPT(char, trans, 'N');
            if constexpr (is_complex_v<T>) { CHECK(trans == 'N' || trans == 'C', trans); }
            else                           { CHECK(trans == 'N' || trans == 'T', trans); }

            ARRAY_IN(T, v, 2);
            ARRAY_IN(T, t, 2);
            ARRAY_INOUT(T, c, 2, overwrite_c != 0);

            CBLAS_INT m = shape(c, 0), n = shape(c, 1), info = 0;
            CBLAS_INT k = shape(v, 1), nb = shape(t, 0);
            CBLAS_INT ldv = std::max<CBLAS_INT>(1, shape(v, 0));
            CBLAS_INT ldt = std::max<CBLAS_INT>(1, shape(t, 0));
            CBLAS_INT ldc = std::max<CBLAS_INT>(1, shape(c, 0));
            CBLAS_INT order = side == 'L' ? m : n;
            CHECK(order >= k && k >= 0, k);
            CHECK(k >= nb && nb >= 1, nb);
            CHECKARRAY(shape(v, 0) == order, v);
            CHECKARRAY(shape(t, 1) == k, t);

            ARRAY_HIDDEN(T, work, 1LL * (side == 'L' ? n : m) * nb);
            lapack::gemqrt(side, trans, m, n, k, nb, v.data<T>(), ldv, t.data<T>(), ldt,
                           c.data<T>(), ldc, work.data<T>(), &info);
            return make_result(c, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *tpqrt(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"l", "nb", "a", "b", "overwrite_a", "overwrite_b",
                                           nullptr};
            static constexpr Ctx<T> ctx("tpqrt", "OOOO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_REQ(CBLAS_INT, l);
            SCALAR_REQ(CBLAS_INT, nb);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);

            /* The pentagonal block `b` gives the shape; `a` is the square upper triangle on top
             * of it, so it is `n`-by-`n`. */
            CBLAS_INT m = shape(b, 0), n = shape(b, 1), info = 0;
            CBLAS_INT lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0));
            CHECK(std::min(m, n) >= l && l >= 0, l);
            CHECK(n >= nb && nb >= 1, nb);
            CHECKARRAY(shape(a, 0) == n && shape(a, 1) == n, a);

            ARRAY_OUT(T, t, 2, true, ctx.zeros(nb, n));
            CBLAS_INT ldt = std::max<CBLAS_INT>(1, shape(t, 0));
            ARRAY_HIDDEN(T, work, nb, n);

            lapack::tpqrt(m, n, l, nb, a.data<T>(), lda, b.data<T>(), ldb, t.data<T>(), ldt,
                          work.data<T>(), &info);
            return make_result(a, b, t, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *tpmqrt(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"l", "v", "t", "a", "b", "side", "trans",
                                           "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("tpmqrt", "OOOOO|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(char, side, 'L');   CHECK(side == 'L' || side == 'R', side);
            SCALAR_OPT(char, trans, 'N');
            if constexpr (is_complex_v<T>) { CHECK(trans == 'N' || trans == 'C', trans); }
            else                           { CHECK(trans == 'N' || trans == 'T', trans); }
            SCALAR_REQ(CBLAS_INT, l);

            ARRAY_IN(T, v, 2);
            ARRAY_IN(T, t, 2);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);

            /* `c` is the pair `(a, b)`; `b` carries the shape and `a` the leading block, whose
             * orientation follows `side`. */
            CBLAS_INT m = shape(b, 0), n = shape(b, 1), info = 0;
            CBLAS_INT k = shape(t, 1), nb = shape(t, 0);
            CBLAS_INT ldv = std::max<CBLAS_INT>(1, shape(v, 0));
            CBLAS_INT ldt = std::max<CBLAS_INT>(1, shape(t, 0));
            CBLAS_INT lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0));
            CHECK(k >= l && l >= 0, l);
            CHECK(k >= nb && nb >= 1, nb);
            CHECKARRAY(shape(v, 0) == (side == 'L' ? m : n) && shape(v, 1) == k, v);
            CHECKARRAY(shape(a, 0) == (side == 'L' ? k : m), a);
            CHECKARRAY(shape(a, 1) == (side == 'L' ? n : k), a);

            ARRAY_HIDDEN(T, work, 1LL * (side == 'L' ? n : m) * nb);
            lapack::tpmqrt(side, trans, m, n, k, l, nb, v.data<T>(), ldv, t.data<T>(), ldt,
                           a.data<T>(), lda, b.data<T>(), ldb, work.data<T>(), &info);
            return make_result(a, b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *tzrzf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("tzrzf", "O|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            /* An upper trapezoidal matrix: never taller than it is wide. */
            CHECKARRAY(shape(a, 1) >= shape(a, 0), a);

            CBLAS_INT m = shape(a, 0), n = shape(a, 1), info = 0;
            CBLAS_INT lda = std::max<CBLAS_INT>(shape(a, 0), 1);
            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(m, 1));
            CHECK(lwork >= m, lwork);

            ARRAY_OUT(T, tau, 1, true, ctx.zeros(m));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(lwork, 1));

            lapack::tzrzf(m, n, a.data<T>(), lda, tau.data<T>(), work.data<T>(), lwork, &info);
            return make_result(a, tau, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *tzrzf_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"m", "n", nullptr};
            static constexpr Ctx<T> ctx("tzrzf_lwork", "OO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, m);
            SCALAR_REQ(CBLAS_INT, n);

            T work{}, a{}, tau{};
            CBLAS_INT info = 0;
            lapack::tzrzf(m, n, &a, std::max<CBLAS_INT>(m, 1), &tau, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *orcsd(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist_real[] = {"x11", "x12", "x21", "x22", "compute_u1",
                                                "compute_u2", "compute_v1t", "compute_v2t",
                                                "trans", "signs", "lwork", "overwrite_x11",
                                                "overwrite_x12", "overwrite_x21",
                                                "overwrite_x22", nullptr};
            static const char *kwlist_complex[] = {"x11", "x12", "x21", "x22", "compute_u1",
                                                   "compute_u2", "compute_v1t", "compute_v2t",
                                                   "trans", "signs", "lwork", "lrwork",
                                                   "overwrite_x11", "overwrite_x12",
                                                   "overwrite_x21", "overwrite_x22", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "uncsd" : "orcsd",
                                        is_complex_v<T> ? "OOOO|OOOOOOOOOOOO" : "OOOO|OOOOOOOOOOO",
                                        is_complex_v<T> ? kwlist_complex : kwlist_real);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_x11);
            SCALAR_FLAG(overwrite_x12);
            SCALAR_FLAG(overwrite_x21);
            SCALAR_FLAG(overwrite_x22);
            SCALAR_OPT(CBLAS_INT, compute_u1, 1);   CHECK(compute_u1 == 0 || compute_u1 == 1, compute_u1);
            SCALAR_OPT(CBLAS_INT, compute_u2, 1);   CHECK(compute_u2 == 0 || compute_u2 == 1, compute_u2);
            SCALAR_OPT(CBLAS_INT, compute_v1t, 1);  CHECK(compute_v1t == 0 || compute_v1t == 1, compute_v1t);
            SCALAR_OPT(CBLAS_INT, compute_v2t, 1);  CHECK(compute_v2t == 0 || compute_v2t == 1, compute_v2t);
            SCALAR_OPT(CBLAS_INT, trans, 0);        CHECK(trans == 0 || trans == 1, trans);
            SCALAR_OPT(CBLAS_INT, signs, 0);        CHECK(signs == 0 || signs == 1, signs);

            ARRAY_INOUT(T, x11, 2, overwrite_x11 != 0);
            ARRAY_INOUT(T, x12, 2, overwrite_x12 != 0);
            ARRAY_INOUT(T, x21, 2, overwrite_x21 != 0);
            ARRAY_INOUT(T, x22, 2, overwrite_x22 != 0);

            CBLAS_INT p = shape(x11, 0), q = shape(x11, 1);
            CBLAS_INT mmp = shape(x22, 0), mmq = shape(x22, 1);
            CHECK(p + mmp == q + mmq, p);
            CBLAS_INT m = p + mmp, info = 0;
            CHECKARRAY(shape(x12, 1) == mmq || shape(x12, 0) == p, x12);
            CHECKARRAY(shape(x21, 0) == mmp || shape(x21, 1) == q, x21);

            CBLAS_INT ldx11 = std::max<CBLAS_INT>(1, shape(x11, 0));
            CBLAS_INT ldx12 = std::max<CBLAS_INT>(1, shape(x12, 0));
            CBLAS_INT ldx21 = std::max<CBLAS_INT>(1, shape(x21, 0));
            CBLAS_INT ldx22 = std::max<CBLAS_INT>(1, shape(x22, 0));
            CBLAS_INT ldu1 = std::max<CBLAS_INT>(1, p), ldu2 = std::max<CBLAS_INT>(1, mmp);
            CBLAS_INT ldv1t = std::max<CBLAS_INT>(1, q), ldv2t = std::max<CBLAS_INT>(1, mmq);

            CBLAS_INT r = std::min(std::min(p, mmp), std::min(q, mmq));
            ARRAY_OUT(R, theta, 1, true, ctx.template zeros_as<R>(r));
            ARRAY_OUT(T, u1, 2, true, ctx.zeros(compute_u1 ? p : 0, compute_u1 ? p : 0));
            ARRAY_OUT(T, u2, 2, true, ctx.zeros(compute_u2 ? mmp : 0, compute_u2 ? mmp : 0));
            ARRAY_OUT(T, v1t, 2, true, ctx.zeros(compute_v1t ? q : 0, compute_v1t ? q : 0));
            ARRAY_OUT(T, v2t, 2, true, ctx.zeros(compute_v2t ? mmq : 0, compute_v2t ? mmq : 0));
            ARRAY_HIDDEN(CBLAS_INT, iwork, p + mmp - r);

            const char u1c = compute_u1 ? 'Y' : 'N', u2c = compute_u2 ? 'Y' : 'N';
            const char v1c = compute_v1t ? 'Y' : 'N', v2c = compute_v2t ? 'Y' : 'N';
            const char tc = trans ? 'T' : 'N', sc = signs ? 'O' : 'D';

            if constexpr (is_complex_v<T>) {
                CBLAS_INT ldeflt, lrdeflt;
                if (!work_size(2LL * m + std::max<CBLAS_INT>(1, std::max(mmp, mmq)) + 1, &ldeflt)) { return nullptr; }
                if (!work_size(5LL * std::max<CBLAS_INT>(1, q - 1) + 4LL * std::max<CBLAS_INT>(1, q) + 8LL * q + 1, &lrdeflt)) { return nullptr; }
                SCALAR_OPT(CBLAS_INT, lwork, ldeflt);    CHECK(lwork == -1 || lwork > 0, lwork);
                SCALAR_OPT(CBLAS_INT, lrwork, lrdeflt);  CHECK(lrwork == -1 || lrwork > 0, lrwork);
                ARRAY_HIDDEN(T, work, lwork);
                ARRAY_HIDDEN(R, rwork, lrwork);
                lapack::orcsd(u1c, u2c, v1c, v2c, tc, sc, m, p, q, x11.data<T>(), ldx11,
                              x12.data<T>(), ldx12, x21.data<T>(), ldx21, x22.data<T>(), ldx22,
                              theta.data<R>(), u1.data<T>(), ldu1, u2.data<T>(), ldu2,
                              v1t.data<T>(), ldv1t, v2t.data<T>(), ldv2t, work.data<T>(), lwork,
                              rwork.data<R>(), lrwork, iwork.data<CBLAS_INT>(), &info);
            }
            else {
                CBLAS_INT ldeflt;
                if (!work_size(2LL + 2LL * m + 5LL * std::max<CBLAS_INT>(1, q - 1) + 4LL * std::max<CBLAS_INT>(1, q) + 8LL * q, &ldeflt)) { return nullptr; }
                SCALAR_OPT(CBLAS_INT, lwork, ldeflt);  CHECK(lwork == -1 || lwork > 0, lwork);
                ARRAY_HIDDEN(T, work, lwork);
                lapack::orcsd(u1c, u2c, v1c, v2c, tc, sc, m, p, q, x11.data<T>(), ldx11,
                              x12.data<T>(), ldx12, x21.data<T>(), ldx21, x22.data<T>(), ldx22,
                              theta.data<R>(), u1.data<T>(), ldu1, u2.data<T>(), ldu2,
                              v1t.data<T>(), ldv1t, v2t.data<T>(), ldv2t, work.data<T>(), lwork,
                              iwork.data<CBLAS_INT>(), &info);
            }
            return make_result(x11, x12, x21, x22, theta, u1, u2, v1t, v2t,
                               static_cast<long long>(info));
        }


        template <class T>
        static PyObject *orcsd_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"m", "p", "q", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "uncsd_lwork" : "orcsd_lwork",
                                        "OOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(CBLAS_INT, m);
            SCALAR_REQ(CBLAS_INT, p);
            SCALAR_REQ(CBLAS_INT, q);

            CBLAS_INT mmp = m - p, mmq = m - q, iwork = 0, info = 0;
            CBLAS_INT ldx11 = std::max<CBLAS_INT>(1, p), ldx12 = std::max<CBLAS_INT>(1, p);
            CBLAS_INT ldx21 = std::max<CBLAS_INT>(1, mmp), ldx22 = std::max<CBLAS_INT>(1, mmp);
            CBLAS_INT ldu1 = std::max<CBLAS_INT>(1, p), ldu2 = std::max<CBLAS_INT>(1, mmp);
            CBLAS_INT ldv1t = std::max<CBLAS_INT>(1, q), ldv2t = std::max<CBLAS_INT>(1, mmq);
            T work{}, x11{}, x12{}, x21{}, x22{}, u1{}, u2{}, v1t{}, v2t{};
            R theta{};

            if constexpr (is_complex_v<T>) {
                R rwork{};
                lapack::orcsd('Y', 'Y', 'Y', 'Y', 'N', 'D', m, p, q, &x11, ldx11, &x12, ldx12,
                              &x21, ldx21, &x22, ldx22, &theta, &u1, ldu1, &u2, ldu2, &v1t,
                              ldv1t, &v2t, ldv2t, &work, -1, &rwork, -1, &iwork, &info);
                return make_result(work, rwork, static_cast<long long>(info));
            }
            else {
                lapack::orcsd('Y', 'Y', 'Y', 'Y', 'N', 'D', m, p, q, &x11, ldx11, &x12, ldx12,
                              &x21, ldx21, &x22, ldx22, &theta, &u1, ldu1, &u2, ldu2, &v1t,
                              ldv1t, &v2t, ldv2t, &work, -1, &iwork, &info);
                return make_result(work, static_cast<long long>(info));
            }
        }


        template <class T>
        static PyObject *gejsv(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "joba", "jobu", "jobv", "jobr", "jobt", "jobp",
                                           "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("gejsv", "O|OOOOOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, joba, 4);  CHECK(0 <= joba && joba < 6, joba);
            SCALAR_OPT(CBLAS_INT, jobu, 0);  CHECK(0 <= jobu && jobu < 4, jobu);
            SCALAR_OPT(CBLAS_INT, jobv, 0);
            CHECK(0 <= jobv && jobv < 4 && (jobv != 1 || jobu < 2), jobv);
            SCALAR_OPT(CBLAS_INT, jobr, 1);  CHECK(jobr == 0 || jobr == 1, jobr);
            SCALAR_OPT(CBLAS_INT, jobt, 0);  CHECK(jobt == 0 || jobt == 1, jobt);
            SCALAR_OPT(CBLAS_INT, jobp, 1);  CHECK(jobp == 0 || jobp == 1, jobp);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1), info = 0;
            CHECK(m >= n, m);
            CBLAS_INT lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CBLAS_INT ldu = std::max<CBLAS_INT>(1, jobu < 3 ? m : 1);
            CBLAS_INT ldv = std::max<CBLAS_INT>(1, jobv < 3 ? n : 1);

            const bool no_u = jobt == 0 && jobu == 3, no_v = jobt == 0 && jobv == 3;
            ARRAY_OUT(T, sva, 1, true, ctx.zeros(n));
            ARRAY_OUT(T, u, 2, true, ctx.zeros(no_u ? 0 : m, no_u ? 0 : (jobu == 1 ? m : n)));
            ARRAY_OUT(T, v, 2, true, ctx.zeros(no_v ? 0 : ldv, no_v ? 0 : n));

            long long floor = std::max({6LL * n + 2LL * n * n, 2LL * m + n, 4LL * n + 1LL * n * n,
                                        2LL * n + 1LL * n * n + 6, 7LL});
            CBLAS_INT deflt;
            if (!work_size(floor, &deflt)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, deflt);
            CHECK(lwork >= 7, lwork);

            ARRAY_HIDDEN(T, work, lwork);
            ARRAY_HIDDEN(CBLAS_INT, iwork, std::max<CBLAS_INT>(3, m + 3 * n));

            lapack::gejsv("CEFGAR"[joba], "UFWN"[jobu], "VJWN"[jobv], jobr ? 'R' : 'N',
                          jobt ? 'T' : 'N', jobp ? 'P' : 'N', m, n, a.data<T>(), lda,
                          sva.data<T>(), u.data<T>(), ldu, v.data<T>(), ldv, work.data<T>(),
                          lwork, iwork.data<CBLAS_INT>(), &info);

            /* LAPACK reports its scaling and rank diagnostics in the first few workspace slots;
             * the `.pyf` copies them out rather than handing back the whole buffer. */
            ARRAY_OUT(T, workout, 1, true, ctx.zeros(7));
            ARRAY_OUT(CBLAS_INT, iworkout, 1, true, ctx.template zeros_as<CBLAS_INT>(3));
            for (CBLAS_INT i = 0; i < 7; ++i) { workout.data<T>()[i] = work.data<T>()[i]; }
            for (CBLAS_INT i = 0; i < 3; ++i) { iworkout.data<CBLAS_INT>()[i] = iwork.data<CBLAS_INT>()[i]; }
            return make_result(sva, u, v, workout, iworkout, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gglse(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "c", "d", "lwork", "overwrite_a",
                                           "overwrite_b", "overwrite_c", "overwrite_d", nullptr};
            static constexpr Ctx<T> ctx("gglse", "OOOO|OOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_FLAG(overwrite_c);
            SCALAR_FLAG(overwrite_d);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1), p = shape(b, 0), info = 0;
            CHECK(m >= 0, m);
            CHECK(n >= 0, n);
            CHECK(p >= n - m && p >= 0, p);
            CHECKARRAY(shape(b, 1) == n, b);
            CBLAS_INT lda = std::max<CBLAS_INT>(shape(a, 0), 1);
            CBLAS_INT ldb = std::max<CBLAS_INT>(shape(b, 0), 1);

            ARRAY_INOUT(T, c, 1, overwrite_c != 0);
            CHECKARRAY(len(c) == m, c);
            /* `d` is copied but not handed back: the `.pyf` marks it `intent(in,copy)` with no
             * `out`, so the constraint vector the caller passed is left untouched. */
            ARRAY_INOUT(T, d, 1, overwrite_d != 0);
            CHECKARRAY(len(d) == p, d);

            ARRAY_OUT(T, x, 1, true, ctx.zeros(n));
            CBLAS_INT deflt;
            if (!work_size(1LL * m + n + p, &deflt)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, deflt);
            CHECK(lwork == -1 || lwork >= 1, lwork);
            ARRAY_HIDDEN(T, work, lwork);

            lapack::gglse(m, n, p, a.data<T>(), lda, b.data<T>(), ldb, c.data<T>(), d.data<T>(),
                          x.data<T>(), work.data<T>(), lwork, &info);
            return make_result(a, b, c, x, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gglse_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"m", "n", "p", nullptr};
            static constexpr Ctx<T> ctx("gglse_lwork", "OOO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, m);  CHECK(m >= 0, m);
            SCALAR_REQ(CBLAS_INT, n);  CHECK(n >= 0, n);
            SCALAR_REQ(CBLAS_INT, p);  CHECK(p >= n - m && p >= 0 && p <= n, p);

            T work{}, a{}, b{}, c{}, d{}, x{};
            CBLAS_INT info = 0;
            lapack::gglse(m, n, p, &a, std::max<CBLAS_INT>(1, m), &b, std::max<CBLAS_INT>(1, p),
                          &c, &d, &x, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *lasd4(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"i", "d", "z", "rho", nullptr};
            static constexpr Ctx<T> ctx("lasd4", "OOO|O", kwlist);
            PARSE_ARGS();

            ARRAY_IN(T, d, 1);
            CBLAS_INT n = len(d), info = 0;
            ARRAY_IN(T, z, 1);
            CHECKARRAY(len(z) == n, z);

            /* `i` selects which root to compute and is 0-based here, 1-based inside LAPACK. */
            SCALAR_REQ(CBLAS_INT, i);  CHECK(i >= 0 && i <= n - 1, i);
            SCALAR_OPT(T, rho, 1.0);

            ARRAY_OUT(T, delta, 1, true, ctx.zeros(n));
            ARRAY_OUT(T, work, 1, true, ctx.zeros(n));
            T sigma = 0;

            lapack::lasd4(n, i + 1, d.data<T>(), z.data<T>(), delta.data<T>(), rho, &sigma,
                          work.data<T>(), &info);
            return make_result(delta, sigma, work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sbev(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"ab", "compute_v", "lower", "ldab", "overwrite_ab",
                                           nullptr};
            static constexpr Ctx<T> ctx("sbev", "O|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_ab);
            if (P.raw("overwrite_ab") == nullptr) { overwrite_ab = 1; }
            SCALAR_OPT(CBLAS_INT, compute_v, 1);  CHECK(compute_v == 1 || compute_v == 0, compute_v);
            SCALAR_OPT(CBLAS_INT, lower, 0);      CHECK(lower == 0 || lower == 1, lower);

            /* `ab` is `intent(in, overwrite)`: reused in place when the flag allows, but never
             * handed back -- the band is consumed, not returned. */
            ARRAY_INOUT(T, ab, 2, overwrite_ab != 0);
            SCALAR_OPT(CBLAS_INT, ldab, shape(ab, 0));
            CHECK(shape(ab, 0) == ldab, ldab);

            CBLAS_INT n = shape(ab, 1), kd = shape(ab, 0) - 1, info = 0;
            CBLAS_INT ldz = compute_v ? n : 1, work_len;
            if (!work_size(std::max(1LL, 3LL * n - 1), &work_len)) { return nullptr; }

            ARRAY_OUT(T, w, 1, true, ctx.zeros(n));
            ARRAY_OUT(T, z, 2, true, ctx.zeros(ldz, ldz));
            ARRAY_HIDDEN(T, work, work_len);

            lapack::sbev(compute_v ? 'V' : 'N', lower ? 'L' : 'U', n, kd, ab.data<T>(), ldab,
                         w.data<T>(), z.data<T>(), ldz, work.data<T>(), &info);
            return make_result(w, z, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sbevd(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist_real[] = {"ab", "compute_v", "lower", "ldab", "liwork", "overwrite_ab", nullptr};
            static const char *kwlist_complex[] = {"ab", "compute_v", "lower", "ldab", "lrwork", "liwork", "overwrite_ab", nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "hbevd" : "sbevd",
                                        is_complex_v<T> ? "O|OOOOOO" : "O|OOOOO",
                                        is_complex_v<T> ? kwlist_complex : kwlist_real);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_ab);
            if (P.raw("overwrite_ab") == nullptr) { overwrite_ab = 1; }
            SCALAR_OPT(CBLAS_INT, compute_v, 1);  CHECK(compute_v == 1 || compute_v == 0, compute_v);
            SCALAR_OPT(CBLAS_INT, lower, 0);      CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, ab, 2, overwrite_ab != 0);
            SCALAR_OPT(CBLAS_INT, ldab, shape(ab, 0));
            CHECK(shape(ab, 0) == ldab, ldab);

            CBLAS_INT n = shape(ab, 1), kd = shape(ab, 0) - 1, info = 0;
            /* The complex `.pyf` forbids `n == 0` outright: its workspace formulas do not cover
             * the empty case.  The real one has no such check, and neither does this. */
            if constexpr (is_complex_v<T>) { CHECK(n > 0, n); }
            CBLAS_INT ldz = compute_v ? n : 1;

            CBLAS_INT liwork_deflt = compute_v ? 3 + 5 * n : 1;
            SCALAR_OPT(CBLAS_INT, liwork, liwork_deflt);
            CHECK(liwork >= (compute_v ? 3 + 5 * n : 1), liwork);

            ARRAY_OUT(R, w, 1, true, ctx.template zeros_as<R>(n));
            ARRAY_OUT(T, z, 2, true, ctx.zeros(ldz, ldz));
            ARRAY_HIDDEN(CBLAS_INT, iwork, liwork);

            CBLAS_INT lwork;
            if constexpr (is_complex_v<T>) {
                if (!work_size(std::max(1LL, compute_v ? 2LL * n * n : 1LL * n), &lwork)) { return nullptr; }
                CBLAS_INT lrwork_deflt;
                if (!work_size(compute_v ? 1LL + 5 * n + 2LL * n * n : 1LL * n, &lrwork_deflt)) { return nullptr; }
                SCALAR_OPT(CBLAS_INT, lrwork, lrwork_deflt);
                CHECK(lrwork >= lrwork_deflt, lrwork);
                ARRAY_HIDDEN(T, work, lwork);
                ARRAY_HIDDEN(R, rwork, lrwork);
                lapack::sbevd(compute_v ? 'V' : 'N', lower ? 'L' : 'U', n, kd, ab.data<T>(), ldab,
                              w.data<R>(), z.data<T>(), ldz, work.data<T>(), lwork,
                              rwork.data<R>(), lrwork, iwork.data<CBLAS_INT>(), liwork, &info);
            }
            else {
                if (!work_size(std::max(1LL, compute_v ? 1LL + 5 * n + 2LL * n * n : 2LL * n), &lwork)) { return nullptr; }
                ARRAY_HIDDEN(T, work, lwork);
                lapack::sbevd(compute_v ? 'V' : 'N', lower ? 'L' : 'U', n, kd, ab.data<T>(), ldab,
                              w.data<R>(), z.data<T>(), ldz, work.data<T>(), lwork,
                              iwork.data<CBLAS_INT>(), liwork, &info);
            }
            return make_result(w, z, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sbevx(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"ab", "vl", "vu", "il", "iu", "ldab", "compute_v",
                                           "range", "lower", "abstol", "mmax", "overwrite_ab",
                                           nullptr};
            static constexpr Ctx<T> ctx(is_complex_v<T> ? "hbevx" : "sbevx", "OOOOO|OOOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_ab);
            if (P.raw("overwrite_ab") == nullptr) { overwrite_ab = 1; }
            SCALAR_OPT(CBLAS_INT, compute_v, 1);  CHECK(compute_v == 1 || compute_v == 0, compute_v);
            SCALAR_OPT(CBLAS_INT, lower, 0);      CHECK(lower == 0 || lower == 1, lower);
            SCALAR_OPT(CBLAS_INT, range, 0);
            CHECK(range == 2 || range == 1 || range == 0, range);

            ARRAY_INOUT(T, ab, 2, overwrite_ab != 0);
            SCALAR_OPT(CBLAS_INT, ldab, shape(ab, 0));
            CHECK(shape(ab, 0) == ldab, ldab);

            CBLAS_INT n = shape(ab, 1), kd = shape(ab, 0) - 1, info = 0, m = 0;
            SCALAR_REQ(R, vl);
            SCALAR_REQ(R, vu);
            SCALAR_REQ(CBLAS_INT, il);  CHECK(il >= 1 && il <= n, il);
            SCALAR_REQ(CBLAS_INT, iu);  CHECK(iu >= 1 && iu <= n && iu >= il, iu);
            SCALAR_OPT(R, abstol, 0.0);

            CBLAS_INT ldq = compute_v ? n : 1, ldz = compute_v ? n : 1;
            /* `mmax` caps how many eigenvectors `z` has room for.  With `range == 1` the count
             * is not known until the call returns, so the default is the whole spectrum and a
             * caller who knows better should say so. */
            SCALAR_OPT(CBLAS_INT, mmax, compute_v ? (range == 2 ? iu - il + 1 : n) : 1);

            ARRAY_OUT(R, w, 1, true, ctx.template zeros_as<R>(n));
            ARRAY_OUT(T, z, 2, true, ctx.zeros(ldz, mmax));
            ARRAY_OUT(CBLAS_INT, ifail, 1, true,
                      ctx.template zeros_as<CBLAS_INT>(compute_v ? n : 1));
            ARRAY_HIDDEN(T, q, ldq, ldq);
            ARRAY_HIDDEN(CBLAS_INT, iwork, 5LL * n);

            /* LAPACK writes one column of `z` per eigenvalue it finds, and how many that is is
             * not known until it returns -- so an `mmax` below `n` is a buffer it can overrun.
             * f2py handed `z` straight over and a short `mmax` corrupted the heap. The call
             * gets a full-size staging buffer instead and the leading `mmax` columns are copied
             * back out; the result is what the caller asked for, and the common `mmax >= n` path
             * still writes into `z` directly. */
            const bool stage = compute_v != 0 && mmax < n;
            ARRAY_HIDDEN(T, staging, ldz, stage ? n : 0);
            T *zdata = stage ? staging.data<T>() : z.data<T>();

            const char jobz = compute_v ? 'V' : 'N';
            const char rangec = range > 0 ? (range == 1 ? 'V' : 'I') : 'A';
            const char uplo = lower ? 'L' : 'U';

            if constexpr (is_complex_v<T>) {
                ARRAY_HIDDEN(T, work, n);
                ARRAY_HIDDEN(R, rwork, 7LL * n);
                lapack::sbevx(jobz, rangec, uplo, n, kd, ab.data<T>(), ldab, q.data<T>(), ldq,
                              vl, vu, il, iu, abstol, &m, w.data<R>(), zdata, ldz,
                              work.data<T>(), rwork.data<R>(), iwork.data<CBLAS_INT>(),
                              ifail.data<CBLAS_INT>(), &info);
            }
            else {
                ARRAY_HIDDEN(T, work, 7LL * n);
                lapack::sbevx(jobz, rangec, uplo, n, kd, ab.data<T>(), ldab, q.data<T>(), ldq,
                              vl, vu, il, iu, abstol, &m, w.data<R>(), zdata, ldz,
                              work.data<T>(), iwork.data<CBLAS_INT>(), ifail.data<CBLAS_INT>(),
                              &info);
            }
            if (stage) {
                const T *from = staging.data<T>();
                T *to = z.data<T>();
                for (npy_intp i = 0; i < 1LL * ldz * mmax; ++i) { to[i] = from[i]; }
            }
            return make_result(w, z, static_cast<long long>(m), ifail,
                               static_cast<long long>(info));
        }


        template <class T>
        static PyObject *lamch(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"cmach", nullptr};
            static constexpr Ctx<T> ctx("lamch", "O", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(char, cmach);
            /* Every letter LAPACK defines is accepted and anything else returns zero, which is
             * what `?LAMCH` itself does -- the `.pyf` adds no check and neither does this. */
            return make_result(lapack::lamch(cmach, static_cast<T *>(nullptr)));
        }


        template <class T>
        static PyObject *lange(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"norm", "a", nullptr};
            static constexpr Ctx<T> ctx("lange", "OO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(char, norm);
            CHECK(norm == 'M' || norm == 'm' || norm == '1' || norm == 'O' || norm == 'o' ||
                  norm == 'I' || norm == 'i' || norm == 'F' || norm == 'f' || norm == 'E' ||
                  norm == 'e', norm);

            ARRAY_IN(T, a, 2);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT lda = std::max<CBLAS_INT>(1, m);
            /* The infinity norm is the only one that needs scratch, and it needs `m` of it; the
             * `.pyf` allocates `m + 1` for every norm and this keeps that. */
            ARRAY_HIDDEN(R, work, m + 1);

            return make_result(lapack::lange(norm, m, n, a.data<T>(), lda, work.data<R>()));
        }


        template <class T>
        static PyObject *lantr(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"norm", "a", "uplo", "diag", nullptr};
            static constexpr Ctx<T> ctx("lantr", "OO|OO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(char, norm);
            CHECK(norm == 'M' || norm == 'm' || norm == '1' || norm == 'O' || norm == 'o' ||
                  norm == 'I' || norm == 'i' || norm == 'F' || norm == 'f' || norm == 'E' ||
                  norm == 'e', norm);
            SCALAR_OPT(char, uplo, 'U');  CHECK(uplo == 'U' || uplo == 'L', uplo);
            SCALAR_OPT(char, diag, 'N');  CHECK(diag == 'N' || diag == 'U', diag);

            ARRAY_IN(T, a, 2);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT lda = std::max<CBLAS_INT>(1, m);
            ARRAY_HIDDEN(R, work, lda);

            return make_result(lapack::lantr(norm, uplo, diag, m, n, a.data<T>(), lda, work.data<R>()));
        }


        template <class T>
        static PyObject *larfg(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "alpha", "x", "incx", "overwrite_x", nullptr};
            static constexpr Ctx<T> ctx("larfg", "OOO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_x);
            SCALAR_REQ(CBLAS_INT, n);  CHECK(n >= 1, n);
            SCALAR_REQ(T, alpha);
            ARRAY_INOUT(T, x, 1, overwrite_x != 0);
            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx > 0 || incx < 0, incx);

            /* The reflector is built from `alpha` and the `n - 1` further entries `incx` steps
             * apart, so `x` holds `1 + (n - 2) * |incx|` of them. */
            npy_intp step = abs(incx);
            CHECKARRAY(len(x) == 1 + (n - 2) * step, x);

            T tau = 0;
            lapack::larfg(n, &alpha, x.data<T>(), incx, &tau);
            return make_result(alpha, x, tau);
        }


        template <class T>
        static PyObject *larf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"v", "tau", "c", "work", "side", "incv",
                                           "overwrite_c", nullptr};
            static constexpr Ctx<T> ctx("larf", "OOOO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);
            ARRAY_IN(T, v, 1);
            SCALAR_REQ(T, tau);
            ARRAY_INOUT(T, c, 2, overwrite_c != 0);
            /* The `.pyf` carries `! FIXME: work should not have been an input argument but kept
             * here for backwards compatibility!`.  It is kept: `work` stays a required argument
             * and is handed to LAPACK as the real workspace, exactly as f2py did. */
            ARRAY_IN(T, work, 1);
            SCALAR_OPT(char, side, 'L');     CHECK(side == 'L' || side == 'R', side);
            SCALAR_OPT(CBLAS_INT, incv, 1);  CHECK(incv > 0 || incv < 0, incv);

            CBLAS_INT m = shape(c, 0), n = shape(c, 1);
            CBLAS_INT ldc = std::max<CBLAS_INT>(1, shape(c, 0));
            npy_intp step = abs(incv);
            CHECKARRAY(len(v) == 1 + ((side == 'L' ? m : n) - 1) * step, v);
            CHECKARRAY(len(work) == (side == 'L' ? n : m), work);

            lapack::larf(side, m, n, v.data<T>(), incv, tau, c.data<T>(), ldc, work.data<T>());
            return make_result(c);
        }


        template <class T>
        static PyObject *lartg(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"f", "g", nullptr};
            static constexpr Ctx<T> ctx("lartg", "OO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(T, f);
            SCALAR_REQ(T, g);

            /* The cosine of a Givens rotation is real whatever the flavor; the sine is not. */
            R cs = 0;
            T sn = 0, r = 0;
            lapack::lartg(f, g, &cs, &sn, &r);
            return make_result(cs, sn, r);
        }


        template <class T>
        static PyObject *rot(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"x", "y", "c", "s", "n", "offx", "incx", "offy",
                                           "incy", "overwrite_x", "overwrite_y", nullptr};
            static constexpr Ctx<T> ctx("rot", "OOOO|OOOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_x);
            SCALAR_FLAG(overwrite_y);
            ARRAY_INOUT(T, x, 1, overwrite_x != 0);
            ARRAY_INOUT(T, y, 1, overwrite_y != 0);
            SCALAR_REQ(R, c);
            SCALAR_REQ(T, s);

            CBLAS_INT lx = len(x), ly = len(y);
            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx > 0 || incx < 0, incx);
            SCALAR_OPT(CBLAS_INT, incy, 1);  CHECK(incy > 0 || incy < 0, incy);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < lx, offx);
            SCALAR_OPT(CBLAS_INT, offy, 0);  CHECK(offy >= 0 && offy < ly, offy);

            npy_intp stepx = abs(incx), stepy = abs(incy);
            SCALAR_OPT(CBLAS_INT, n, (lx - 1 - offx) / stepx + 1);
            CHECK(lx - offx > (n - 1) * stepx, n);
            CHECK(ly - offy > (n - 1) * stepy, n);

            lapack::rot(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy, c, s);
            return make_result(x, y);
        }


        template <class T>
        static PyObject *tgsyl(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "c", "d", "e", "f", "trans", "ijob",
                                           "lwork", "overwrite_c", "overwrite_f", nullptr};
            static constexpr Ctx<T> ctx("tgsyl", "OOOOOO|OOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);
            SCALAR_FLAG(overwrite_f);
            SCALAR_OPT(char, trans, 'N');   CHECK(trans == 'N' || trans == 'T', trans);
            SCALAR_OPT(CBLAS_INT, ijob, 0); CHECK(ijob >= 0 && ijob <= 4, ijob);

            ARRAY_IN(T, a, 2);
            ARRAY_IN(T, b, 2);
            ARRAY_INOUT(T, c, 2, overwrite_c != 0);
            ARRAY_IN(T, d, 2);
            ARRAY_IN(T, e, 2);
            ARRAY_INOUT(T, f, 2, overwrite_f != 0);

            /* The `.pyf` takes the orders from the *column* counts of `a` and `b`. */
            CBLAS_INT m = shape(a, 1), n = shape(b, 1), info = 0;
            CBLAS_INT lda = std::max<CBLAS_INT>(shape(a, 0), 1);
            CBLAS_INT ldb = std::max<CBLAS_INT>(shape(b, 0), 1);
            CBLAS_INT ldc = std::max<CBLAS_INT>(shape(c, 0), 1);
            CBLAS_INT ldd = std::max<CBLAS_INT>(shape(d, 0), 1);
            CBLAS_INT lde = std::max<CBLAS_INT>(shape(e, 0), 1);
            CBLAS_INT ldf = std::max<CBLAS_INT>(shape(f, 0), 1);
            CHECKARRAY(shape(a, 0) == m, a);
            CHECKARRAY(shape(b, 0) == n, b);
            CHECKARRAY(shape(c, 0) == m && shape(c, 1) == n, c);
            CHECKARRAY(shape(d, 0) == m && shape(d, 1) == m, d);
            CHECKARRAY(shape(e, 0) == n && shape(e, 1) == n, e);
            CHECKARRAY(shape(f, 0) == m && shape(f, 1) == n, f);

            CBLAS_INT deflt;
            if (!work_size(2LL * m * n, &deflt)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, deflt);
            /* Only the two `ijob` values that estimate `dif` with a transposed system need the
             * full `2 * m * n`; the rest get by with `n`. */
            CHECK((ijob == 1 || ijob == 2) && trans == 'N' ? lwork >= std::max<CBLAS_INT>(1, deflt)
                                                           : (lwork >= n || lwork == -1), lwork);

            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            ARRAY_HIDDEN(CBLAS_INT, iwork, 1LL * m + n + 6);
            T scale = 0, dif = 0;

            lapack::tgsyl(trans, ijob, m, n, a.data<T>(), lda, b.data<T>(), ldb, c.data<T>(),
                          ldc, d.data<T>(), ldd, e.data<T>(), lde, f.data<T>(), ldf, &scale,
                          &dif, work.data<T>(), lwork, iwork.data<CBLAS_INT>(), &info);
            return make_result(c, f, scale, dif, static_cast<long long>(info));
        }


        /** @brief `ilaver` reports the LAPACK version and has no flavor at all, so it is a
         *         plain function rather than a template and gets a bare method-table row. */
        static PyObject *ilaver(PyObject *Py_UNUSED(self), PyObject *Py_UNUSED(args)) noexcept
        {
            CBLAS_INT major = 0, minor = 0, patch = 0;
            lapack::ilaver(&major, &minor, &patch);
            return make_result(static_cast<long long>(major), static_cast<long long>(minor),
                               static_cast<long long>(patch));
        }


        /** @brief Emit the four method-table rows of one orthogonal/unitary family:
         *         `s/dor<fam>` and `c/zun<fam>`, all bound to the one `or<fam>` template.
         *
         * The nine families below are the only place in the port where the real and complex
         * spellings differ by a fixed prefix while sharing a single wrapper, so this stays here
         * beside the table rather than joining FAMILY/ROW in `wrapper_helpers.hpp`.  It is not a
         * shape `lapack_sym_herm.cpp` could reuse: its `sy`/`he` pairs look similar but split
         * across two templates, so no four-row macro fits them.
         */
        #define FAMILY_ORUN(fam)                                                                  \
            ROW(sor##fam, or##fam, f32),  ROW(dor##fam, or##fam, f64),                            \
            ROW(cun##fam, or##fam, c64),  ROW(zun##fam, or##fam, c128)

        PyMethodDef other_methods[] = {
            FAMILY(tpttf),
            FAMILY(tpttr),
            FAMILY(tfttp),
            FAMILY(tfttr),
            FAMILY(trttf),
            FAMILY(trttp),
            FAMILY(tfsm),
            FAMILY(ppcon),
            FAMILY(ppsv),
            FAMILY(pptrf),
            FAMILY(pptri),
            FAMILY(pptrs),
            FAMILY(pftrf),
            FAMILY(pftri),
            FAMILY(pftrs),
            FAMILY(pbtrf),
            FAMILY(pbtrs),
            FAMILY(pbsv),
            FAMILY(trtrs),
            FAMILY(trcon),
            FAMILY(tbtrs),
            FAMILY(trtri),
            FAMILY(lauum),
            FAMILY(laswp),
            FAMILY(trexc),
            FAMILY(tgexc),
            FAMILY(trsyl),
            FAMILY(trsen),
            FAMILY(trsen_lwork),
            FAMILY(tgsen),
            FAMILY(tgsen_lwork),
            FAMILY(geqrt),
            FAMILY(gemqrt),
            FAMILY(tpqrt),
            FAMILY(tpmqrt),
            FAMILY(tzrzf),
            FAMILY(tzrzf_lwork),
            FAMILY(gglse),
            FAMILY(gglse_lwork),
            FAMILY(lange),
            FAMILY(lantr),
            FAMILY(larfg),
            FAMILY(larf),
            FAMILY(lartg),
            ROW(ssfrk, sfrk, f32),
            ROW(dsfrk, sfrk, f64),
            ROW(chfrk, sfrk, c64),
            ROW(zhfrk, sfrk, c128),
            ROW(stgsyl, tgsyl, f32),
            ROW(dtgsyl, tgsyl, f64),
            FAMILY_ORUN(ghr),
            FAMILY_ORUN(ghr_lwork),
            FAMILY_ORUN(gqr),
            FAMILY_ORUN(grq),
            FAMILY_ORUN(mqr),
            FAMILY_ORUN(mrz),
            FAMILY_ORUN(mrz_lwork),
            FAMILY_ORUN(csd),
            FAMILY_ORUN(csd_lwork),
            ROW(sgejsv, gejsv, f32),
            ROW(dgejsv, gejsv, f64),
            ROW(slasd4, lasd4, f32),
            ROW(dlasd4, lasd4, f64),
            ROW(ssbev, sbev, f32),
            ROW(dsbev, sbev, f64),
            ROW(ssbevd, sbevd, f32),
            ROW(dsbevd, sbevd, f64),
            ROW(chbevd, sbevd, c64),
            ROW(zhbevd, sbevd, c128),
            ROW(ssbevx, sbevx, f32),
            ROW(dsbevx, sbevx, f64),
            ROW(chbevx, sbevx, c64),
            ROW(zhbevx, sbevx, c128),
            ROW(slamch, lamch, f32),
            ROW(dlamch, lamch, f64),
            ROW(crot, rot, c64),
            ROW(zrot, rot, c128),
            {"ilaver", (PyCFunction)ilaver, METH_NOARGS, nullptr},
            {nullptr, nullptr, 0, nullptr},
        };

        #undef FAMILY_ORUN

    }  // namespace capi
}  // namespace lapack
