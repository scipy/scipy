/**
 * @file
 * @brief Python wrappers for the positive definite LAPACK routines.
 *
 */
#define PY_ARRAY_UNIQUE_SYMBOL scipy_lapack_ARRAY_API
#define NO_IMPORT_ARRAY
#include "lapack_helpers.hpp"
#include "lapack_calls.hpp"

namespace lapack {
    namespace capi {


        template <class T>
        static PyObject *pstrf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "tol", "lower", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("pstrf", "O|OOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(R, tol, -1.0);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CBLAS_INT rank_c = 0, info = 0, work_len;
            if (!work_size(2LL * n, &work_len)) { return nullptr; }

            ARRAY_OUT(CBLAS_INT, piv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            ARRAY_HIDDEN(R, work, work_len);

            lapack::pstrf(lower ? 'L' : 'U', n, a.data<T>(), lda, piv.data<CBLAS_INT>(), &rank_c, tol, work.data<R>(), &info);
            return make_result(a, piv, static_cast<long long>(rank_c), static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pstf2(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "tol", "lower", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("pstf2", "O|OOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(R, tol, -1.0);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CBLAS_INT rank_c = 0, info = 0, work_len;
            if (!work_size(2LL * n, &work_len)) { return nullptr; }

            ARRAY_OUT(CBLAS_INT, piv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            ARRAY_HIDDEN(R, work, work_len);

            lapack::pstf2(lower ? 'L' : 'U', n, a.data<T>(), lda, piv.data<CBLAS_INT>(), &rank_c, tol, work.data<R>(), &info);
            return make_result(a, piv, static_cast<long long>(rank_c), static_cast<long long>(info));
        }


        template <class T>
        static PyObject *posv(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "lower", "overwrite_a", "overwrite_b",
                                           nullptr};
            static constexpr Ctx<T> ctx("posv", "OO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), info = 0;

            lapack::posv(lower ? 'L' : 'U', n, nrhs, a.data<T>(), n, b.data<T>(), n, &info);
            return make_result(a, b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *posvx(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "fact", "af", "equed", "s", "lower",
                                           "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("posvx", "OO|OOOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(char, fact, 'E');
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);
            /* An input when `fact` is 'F', an output otherwise */
            SCALAR_OPT(char, equed, 'Y');

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = shape(a, 0);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), ldb = shape(b, 0);

            ARRAY_OUT(T, af, 2, true, ctx.zeros(n, n));
            CHECKARRAY(shape(af, 0) == n && shape(af, 1) == n, af);
            CBLAS_INT ldaf = shape(af, 0);

            ARRAY_OUT(R, s, 1, true, ctx.template zeros_as<R>(n));
            CHECKARRAY(len(s) == n, s);

            ARRAY_OUT(T, x, 2, true, ctx.zeros(n, nrhs));
            ARRAY_OUT(R, ferr, 1, true, ctx.template zeros_as<R>(nrhs));
            ARRAY_OUT(R, berr, 1, true, ctx.template zeros_as<R>(nrhs));

            CBLAS_INT ldx = shape(x, 0), work_len, info = 0;
            if (!work_size(is_complex_v<T> ? 2LL * n : 3LL * n, &work_len)) { return nullptr; }
            using W = std::conditional_t<is_complex_v<T>, R, CBLAS_INT>;
            ARRAY_HIDDEN(T, work, work_len);
            ARRAY_HIDDEN(W, irwork, n);
            R rcond = 0;

            lapack::posvx(fact, lower ? 'L' : 'U', n, nrhs, a.data<T>(), lda, af.data<T>(), ldaf,
                          &equed, s.data<R>(), b.data<T>(), ldb, x.data<T>(), ldx, &rcond,
                          ferr.data<R>(), berr.data<R>(), work.data<T>(), irwork.data<W>(), &info);
            return make_result(a, af, equed, s, b, x, rcond, ferr, berr, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pocon(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "anorm", "uplo", nullptr};
            static constexpr Ctx<T> ctx("pocon", "OO|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(char, uplo, 'U');

            ARRAY_IN(T, a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            SCALAR_REQ(R, anorm);

            CBLAS_INT n = shape(a, 0), lda = shape(a, 0), work_len, info = 0;
            if (!work_size(is_complex_v<T> ? 2LL * n : 3LL * n, &work_len)) { return nullptr; }
            using W = std::conditional_t<is_complex_v<T>, R, CBLAS_INT>;
            ARRAY_HIDDEN(T, work, work_len);
            ARRAY_HIDDEN(W, irwork, n);
            R rcond = 0;

            lapack::pocon(uplo, n, a.data<T>(), lda, anorm, &rcond, work.data<T>(), irwork.data<W>(), &info);
            return make_result(rcond, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *potrf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lower", "clean", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("potrf", "O|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);
            SCALAR_OPT(CBLAS_INT, clean, 1);  CHECK(clean == 0 || clean == 1, clean);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, n), info = 0;

            lapack::potrf(lower ? 'L' : 'U', n, a.data<T>(), lda, &info);

            if (clean) {
                T *c = a.data<T>();
                for (CBLAS_INT i = 0; i < n; ++i) {
                    for (CBLAS_INT j = i + 1; j < n; ++j) {
                        c[lower ? j * n + i : i * n + j] = T{};
                    }
                }
            }
            return make_result(a, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *potrs(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"c", "b", "lower", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("potrs", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_IN(T, c, 2);
            CHECKARRAY(shape(c, 0) == shape(c, 1), c);
            CBLAS_INT n = shape(c, 0);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), info = 0;

            lapack::potrs(lower ? 'L' : 'U', n, nrhs, c.data<T>(), n, b.data<T>(), n, &info);
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *potri(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"c", "lower", "overwrite_c", nullptr};
            static constexpr Ctx<T> ctx("potri", "O|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);
            SCALAR_OPT(CBLAS_INT, lower, 0);
            CHECK(lower == 0 || lower == 1, lower);

            ARRAY_INOUT(T, c, 2, overwrite_c != 0);
            CHECKARRAY(shape(c, 0) == shape(c, 1), c);
            CBLAS_INT n = shape(c, 0), info = 0;

            lapack::potri(lower ? 'L' : 'U', n, c.data<T>(), n, &info);
            return make_result(c, static_cast<long long>(info));
        }


        PyMethodDef pos_def_methods[] = {
            FAMILY(pstrf),
            FAMILY(pstf2),
            FAMILY(posv),
            FAMILY(posvx),
            FAMILY(pocon),
            FAMILY(potrf),
            FAMILY(potrs),
            FAMILY(potri),
            {nullptr, nullptr, 0, nullptr},
        };

    }  // namespace capi
}  // namespace lapack
