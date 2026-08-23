/**
 * @file
 * @brief Python wrappers for the LAPACK routines declared in `flapack_gen.pyf.src`.
 */
#define PY_ARRAY_UNIQUE_SYMBOL scipy_lapack_ARRAY_API
#define NO_IMPORT_ARRAY
#include "lapack_helpers.hpp"
#include "lapack_callback.hpp"
#include "lapack_calls.hpp"

namespace lapack {
    namespace capi {


        template <class T>
        static PyObject *gees(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {gees_select_traits<T>::kwname, "a", "compute_v", "sort_t",
                                        "lwork", gees_select_traits<T>::extra_kwname, "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("gees", "OO|OOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT n = shape(a, 1), lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CHECKARRAY(shape(a, 0) == n, a);

            SCALAR_OPT(CBLAS_INT, compute_v, 1);  CHECK(compute_v == 0 || compute_v == 1, compute_v);
            SCALAR_OPT(CBLAS_INT, sort_t, 0);     CHECK(sort_t == 0 || sort_t == 1, sort_t);

            CBLAS_INT default_lwork, minimum_lwork;
            if (!work_size(3LL * n, &default_lwork) || !work_size(2LL * n, &minimum_lwork)) {
                return nullptr;
            }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);
            CHECK(lwork == -1 || lwork >= minimum_lwork, lwork);

            CALLABLE_SELECT(gees, select);
            CBLAS_INT vs_rows = compute_v ? n : 1, ldvs = std::max<CBLAS_INT>(1, vs_rows);
            CBLAS_INT sdim = 0, info = 0;
            ARRAY_HIDDEN(CBLAS_INT, bwork, n);

            if constexpr (is_complex_v<T>) {
                using R = real_of_t<T>;
                ARRAY_OUT(T, w, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, vs, 2, true, ctx.zeros(vs_rows, n));
                ARRAY_OUT(T, work, 1, true, ctx.zeros(lwork > 1 ? lwork : 1));
                ARRAY_HIDDEN(R, rwork, n);

                CALLABLE_CALL(select, lapack::gees(compute_v ? 'V' : 'N', sort_t ? 'S' : 'N', select, n,
                                                    a.data<T>(), lda, &sdim, w.data<T>(), vs.data<T>(), ldvs,
                                                    work.data<T>(), lwork, rwork.data<R>(), bwork.data<CBLAS_INT>(), &info));
                return make_result(a, static_cast<long long>(sdim), w, vs, work,
                                static_cast<long long>(info));
            }
            else {
                ARRAY_OUT(T, wr, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, wi, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, vs, 2, true, ctx.zeros(vs_rows, n));
                ARRAY_OUT(T, work, 1, true, ctx.zeros(lwork > 1 ? lwork : 1));

                CALLABLE_CALL(select, lapack::gees(compute_v ? 'V' : 'N', sort_t ? 'S' : 'N', select, n,
                                                    a.data<T>(), lda, &sdim, wr.data<T>(), wi.data<T>(), vs.data<T>(), ldvs,
                                                    work.data<T>(), lwork, bwork.data<CBLAS_INT>(), &info));
                return make_result(a, static_cast<long long>(sdim), wr, wi, vs, work,
                                static_cast<long long>(info));
            }
        }


        template <class T>
        static PyObject *gges(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {gges_select_traits<T>::kwname, "a", "b", "jobvsl", "jobvsr",
                                           "sort_t", "ldvsl", "ldvsr", "lwork",
                                           gges_select_traits<T>::extra_kwname, "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("gges", "OOO|OOOOOOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CBLAS_INT n = shape(a, 1);
            CBLAS_INT lda = std::max<CBLAS_INT>(1, shape(a, 0)), ldb = std::max<CBLAS_INT>(1, shape(b, 0));
            CHECKARRAY(shape(a, 0) >= n, a);
            CHECKARRAY(shape(b, 1) == n && shape(b, 0) >= n, b);

            SCALAR_OPT(CBLAS_INT, jobvsl, 1);  CHECK(jobvsl == 0 || jobvsl == 1, jobvsl);
            SCALAR_OPT(CBLAS_INT, jobvsr, 1);  CHECK(jobvsr == 0 || jobvsr == 1, jobvsr);
            SCALAR_OPT(CBLAS_INT, sort_t, 0);  CHECK(sort_t == 0 || sort_t == 1, sort_t);
            SCALAR_OPT(CBLAS_INT, ldvsl, std::max<CBLAS_INT>(1, jobvsl ? n : 1));
            CHECK(ldvsl >= std::max<CBLAS_INT>(1, jobvsl ? n : 1), ldvsl);
            SCALAR_OPT(CBLAS_INT, ldvsr, std::max<CBLAS_INT>(1, jobvsr ? n : 1));
            CHECK(ldvsr >= std::max<CBLAS_INT>(1, jobvsr ? n : 1), ldvsr);

            CBLAS_INT default_lwork, minimum_lwork;
            if constexpr (is_complex_v<T>) {
                if (!work_size(2LL * n, &default_lwork) || !work_size(2LL * n, &minimum_lwork)) {
                    return nullptr;
                }
            }
            else {
                CBLAS_INT alternate_minimum;
                if (!work_size(8LL * n + 16, &default_lwork) || !work_size(8LL * n, &minimum_lwork) ||
                    !work_size(6LL * n + 16, &alternate_minimum)) {
                    return nullptr;
                }
                minimum_lwork = minimum_lwork > alternate_minimum ? minimum_lwork : alternate_minimum;
            }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);
            CHECK(lwork == -1 || lwork >= minimum_lwork, lwork);

            CALLABLE_SELECT(gges, select);
            CBLAS_INT sdim = 0, info = 0;
            ARRAY_HIDDEN(CBLAS_INT, bwork, n);

            if constexpr (is_complex_v<T>) {
                using R = real_of_t<T>;
                CBLAS_INT rwork_size;
                if (!work_size(8LL * n, &rwork_size)) { return nullptr; }
                ARRAY_OUT(T, alpha, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, beta, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, vsl, 2, true, ctx.zeros(ldvsl, n));
                ARRAY_OUT(T, vsr, 2, true, ctx.zeros(ldvsr, n));
                ARRAY_OUT(T, work, 1, true, ctx.zeros(lwork > 1 ? lwork : 1));
                ARRAY_HIDDEN(R, rwork, rwork_size);

                CALLABLE_CALL(select, lapack::gges(jobvsl ? 'V' : 'N', jobvsr ? 'V' : 'N', sort_t ? 'S' : 'N',
                                                    select, n, a.data<T>(), lda, b.data<T>(), ldb, &sdim, alpha.data<T>(),
                                                    beta.data<T>(), vsl.data<T>(), ldvsl, vsr.data<T>(), ldvsr, work.data<T>(),
                                                    lwork, rwork.data<R>(), bwork.data<CBLAS_INT>(), &info));
                return make_result(a, b, static_cast<long long>(sdim), alpha, beta, vsl, vsr, work,
                                static_cast<long long>(info));
            }
            else {
                ARRAY_OUT(T, alphar, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, alphai, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, beta, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, vsl, 2, true, ctx.zeros(ldvsl, n));
                ARRAY_OUT(T, vsr, 2, true, ctx.zeros(ldvsr, n));
                ARRAY_OUT(T, work, 1, true, ctx.zeros(lwork > 1 ? lwork : 1));

                CALLABLE_CALL(select, lapack::gges(jobvsl ? 'V' : 'N', jobvsr ? 'V' : 'N', sort_t ? 'S' : 'N',
                                                    select, n, a.data<T>(), lda, b.data<T>(), ldb, &sdim, alphar.data<T>(),
                                                    alphai.data<T>(), beta.data<T>(), vsl.data<T>(), ldvsl, vsr.data<T>(), ldvsr,
                                                    work.data<T>(), lwork, bwork.data<CBLAS_INT>(), &info));
                return make_result(a, b, static_cast<long long>(sdim), alphar, alphai, beta, vsl, vsr,
                                work, static_cast<long long>(info));
            }
        }


        /**
         * `lo`/`hi` are 0-based on the Python side and 1-based in Fortran, so they are shifted
         * across the call.  Here they are outputs, so the shift is a decrement afterwards.
         */
        template <class T>
        static PyObject *gebal(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "scale", "permute", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("gebal", "O|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, permute, 0);
            SCALAR_OPT(CBLAS_INT, scale, 0);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) >= shape(a, 1), a);
            CBLAS_INT n = shape(a, 1), lda = std::max<CBLAS_INT>(1, shape(a, 0));

            using R = real_of_t<T>;
            ARRAY_OUT(R, pivscale, 1, true, ctx.template zeros_as<R>(n));
            CBLAS_INT lo = 0, hi = 0, info = 0;

            lapack::gebal(permute ? (scale ? 'B' : 'P') : (scale ? 'S' : 'N'), n, a.data<T>(), lda, &lo, &hi, pivscale.data<R>(), &info);
            lo -= 1;
            hi -= 1;
            return make_result(a, static_cast<long long>(lo), static_cast<long long>(hi), pivscale, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gehrd(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lo", "hi", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("gehrd", "O|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0));

            SCALAR_OPT(CBLAS_INT, lo, 0);
            SCALAR_OPT(CBLAS_INT, hi, n - 1);
            ARRAY_OUT(T, tau, 1, true, ctx.zeros(n > 0 ? n - 1 : 0));

            CBLAS_INT minimum_lwork = n > 1 ? n : 1;
            SCALAR_OPT(CBLAS_INT, lwork, minimum_lwork);
            CHECK(lwork >= minimum_lwork, lwork);
            ARRAY_HIDDEN(T, work, lwork);

            CBLAS_INT info = 0;
            lapack::gehrd(n, lo + 1, hi + 1, a.data<T>(), lda, tau.data<T>(), work.data<T>(), lwork, &info);
            return make_result(a, tau, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gehrd_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "lo", "hi", nullptr};
            static constexpr Ctx<T> ctx("gehrd_lwork", "O|OO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, lo, 0);
            SCALAR_OPT(CBLAS_INT, hi, n - 1);

            T work{};
            CBLAS_INT info = 0;
            lapack::gehrd(n, lo + 1, hi + 1, nullptr, std::max<CBLAS_INT>(1, n), nullptr, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gesv(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("gesv", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0); CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0); CHECKARRAY(shape(a, 0) == shape(b, 0), b);

            CBLAS_INT n = shape(a, 0), nrhs = shape(b, 1), info = 0;
            CBLAS_INT lda = std::max<CBLAS_INT>(1, shape(a, 0)), ldb = std::max<CBLAS_INT>(1, shape(b, 0));
            ARRAY_OUT(CBLAS_INT, piv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            CBLAS_INT *pivots = piv.data<CBLAS_INT>();

            lapack::gesv(n, nrhs, a.data<T>(), lda, pivots, b.data<T>(), ldb, &info);
            for (CBLAS_INT i = 0; i < n; ++i) { pivots[i] -= 1; }
            return make_result(a, piv, b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gecon(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "anorm", "norm", nullptr};
            static constexpr Ctx<T> ctx("gecon", "OO|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(char, norm, '1');
            ARRAY_IN(T, a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            SCALAR_REQ(R, anorm);

            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), work_len, irwork_len;
            if (!work_size(is_complex_v<T> ? 2LL * n : 4LL * n, &work_len) ||
                !work_size(is_complex_v<T> ? 2LL * n : 1LL * n, &irwork_len)) {
                return nullptr;
            }
            using W = std::conditional_t<is_complex_v<T>, R, CBLAS_INT>;
            ARRAY_HIDDEN(T, work, work_len);
            ARRAY_HIDDEN(W, irwork, irwork_len);

            R rcond = 0;
            CBLAS_INT info = 0;
            lapack::gecon(norm, n, a.data<T>(), lda, anorm, &rcond, work.data<T>(), irwork.data<W>(), &info);
            return make_result(rcond, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *getrf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("getrf", "O|O", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);

            CBLAS_INT m = shape(a, 0), n = shape(a, 1), info = 0;
            CBLAS_INT npiv = m < n ? m : n, lda = std::max<CBLAS_INT>(1, shape(a, 0));
            ARRAY_OUT(CBLAS_INT, piv, 1, true, ctx.template zeros_as<CBLAS_INT>(npiv));
            CBLAS_INT *pivots = piv.data<CBLAS_INT>();

            lapack::getrf(m, n, a.data<T>(), lda, pivots, &info);
            for (CBLAS_INT i = 0; i < npiv; ++i) { pivots[i] -= 1; }
            return make_result(a, piv, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *getrs(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"lu", "piv", "b", "trans", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("getrs", "OOO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, trans, 0);
            CHECK(trans >= 0 && trans <= 2, trans);

            ARRAY_IN(T, lu, 2);
            CHECKARRAY(shape(lu, 0) == shape(lu, 1), lu);
            CBLAS_INT n = shape(lu, 0), ldlu = std::max<CBLAS_INT>(1, shape(lu, 0));

            ARRAY_IN(CBLAS_INT, piv, 1);
            CHECKARRAY(len(piv) == n, piv);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), ldb = std::max<CBLAS_INT>(1, shape(b, 0)), info = 0;
            ARRAY_HIDDEN(CBLAS_INT, pivots, n);
            const CBLAS_INT *supplied = piv.data<CBLAS_INT>();
            CBLAS_INT *shifted = pivots.data<CBLAS_INT>();
            for (CBLAS_INT i = 0; i < n; ++i) { shifted[i] = supplied[i] + 1; }

            lapack::getrs(trans ? (trans == 2 ? 'C' : 'T') : 'N', n, nrhs, lu.data<T>(), ldlu, shifted, b.data<T>(), ldb, &info);
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *getc2(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("getc2", "O|O", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0)), info = 0;

            ARRAY_OUT(CBLAS_INT, ipiv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            ARRAY_OUT(CBLAS_INT, jpiv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            CBLAS_INT *row_pivots = ipiv.data<CBLAS_INT>(), *col_pivots = jpiv.data<CBLAS_INT>();

            lapack::getc2(n, a.data<T>(), lda, row_pivots, col_pivots, &info);
            for (CBLAS_INT i = 0; i < n; ++i) { row_pivots[i] -= 1; col_pivots[i] -= 1; }
            return make_result(a, ipiv, jpiv, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gesc2(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"lu", "rhs", "ipiv", "jpiv", "overwrite_rhs", nullptr};
            static constexpr Ctx<T> ctx("gesc2", "OOOO|O", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_rhs);
            ARRAY_IN(T, lu, 2);
            CHECKARRAY(shape(lu, 0) == shape(lu, 1), lu);
            CHECKARRAY(shape(lu, 0) > 0, lu);
            CBLAS_INT n = shape(lu, 0), lda = std::max<CBLAS_INT>(1, shape(lu, 0));

            ARRAY_INOUT(T, rhs, 1, overwrite_rhs != 0);
            CHECKARRAY(len(rhs) == n, rhs);
            ARRAY_IN(CBLAS_INT, ipiv, 1);
            CHECKARRAY(len(ipiv) == n, ipiv);
            ARRAY_IN(CBLAS_INT, jpiv, 1);
            CHECKARRAY(len(jpiv) == n, jpiv);

            ARRAY_HIDDEN(CBLAS_INT, rows, n);
            ARRAY_HIDDEN(CBLAS_INT, cols, n);
            const CBLAS_INT *supplied_row = ipiv.data<CBLAS_INT>(), *supplied_col = jpiv.data<CBLAS_INT>();
            CBLAS_INT *row_pivots = rows.data<CBLAS_INT>(), *col_pivots = cols.data<CBLAS_INT>();
            for (CBLAS_INT i = 0; i < n; ++i) { row_pivots[i] = supplied_row[i] + 1; col_pivots[i] = supplied_col[i] + 1; }
            real_of_t<T> scale = 0;

            lapack::gesc2(n, lu.data<T>(), lda, rhs.data<T>(), row_pivots, col_pivots, &scale);
            return make_result(rhs, scale);
        }


        template <class T>
        static PyObject *getri(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"lu", "piv", "lwork", "overwrite_lu", nullptr};
            static constexpr Ctx<T> ctx("getri", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_lu);
            ARRAY_INOUT(T, lu, 2, overwrite_lu != 0);
            CHECKARRAY(shape(lu, 0) == shape(lu, 1), lu);
            CBLAS_INT n = shape(lu, 0), lda = std::max<CBLAS_INT>(1, shape(lu, 0)), info = 0;

            ARRAY_IN(CBLAS_INT, piv, 1);
            CHECKARRAY(len(piv) == n, piv);

            CBLAS_INT default_lwork;
            if (!work_size(3LL * n, &default_lwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);
            CHECK(lwork >= n, lwork);
            ARRAY_HIDDEN(T, work, lwork);

            ARRAY_HIDDEN(CBLAS_INT, pivots, n);
            const CBLAS_INT *supplied = piv.data<CBLAS_INT>();
            CBLAS_INT *shifted = pivots.data<CBLAS_INT>();
            for (CBLAS_INT i = 0; i < n; ++i) { shifted[i] = supplied[i] + 1; }

            lapack::getri(n, lu.data<T>(), lda, shifted, work.data<T>(), lwork, &info);
            return make_result(lu, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *getri_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", nullptr};
            static constexpr Ctx<T> ctx("getri_lwork", "O", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);

            T work{};
            CBLAS_INT info = 0;
            lapack::getri(n, nullptr, std::max<CBLAS_INT>(1, n), nullptr, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gesdd(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "compute_uv", "full_matrices", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("gesdd", "O|OOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, compute_uv, 1);     CHECK(compute_uv == 0 || compute_uv == 1, compute_uv);
            SCALAR_OPT(CBLAS_INT, full_matrices, 1);  CHECK(full_matrices == 0 || full_matrices == 1, full_matrices);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT minmn = m < n ? m : n, maxmn = m > n ? m : n;
            CBLAS_INT u0  = compute_uv ? m : 1;
            CBLAS_INT u1  = compute_uv ? (full_matrices ? m : minmn) : 1;
            CBLAS_INT vt0 = compute_uv ? (full_matrices ? n : minmn) : 1;
            CBLAS_INT vt1 = compute_uv ? n : 1;
            CBLAS_INT lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CBLAS_INT ldu = std::max<CBLAS_INT>(1, u0), ldvt = std::max<CBLAS_INT>(1, vt0);

            CBLAS_INT default_lwork;
            if constexpr (is_complex_v<T>) {
                if (!work_size(compute_uv ? 2LL * minmn * minmn + maxmn + 2LL * minmn : 2LL * minmn + maxmn, &default_lwork)) {
                    return nullptr;
                }
            }
            else {
                /* 10*minmn + 2 + 25*(25+8) is the pyf's own spelling of the no-vectors bound. */
                long long bound = 14LL * minmn + 4;
                if (10LL * minmn + 827 > bound) { bound = 10LL * minmn + 827; }
                if (!work_size(compute_uv ? 4LL * minmn * minmn + maxmn + 9LL * minmn : bound + maxmn, &default_lwork)) {
                    return nullptr;
                }
            }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);

            CBLAS_INT iwork_len;
            if (!work_size(8LL * minmn, &iwork_len)) { return nullptr; }
            ARRAY_OUT(R, s, 1, true, ctx.template zeros_as<R>(minmn));
            ARRAY_OUT(T, u, 2, true, ctx.zeros(u0, u1));
            ARRAY_OUT(T, vt, 2, true, ctx.zeros(vt0, vt1));
            ARRAY_HIDDEN(T, work, lwork);
            ARRAY_HIDDEN(CBLAS_INT, iwork, iwork_len);

            char job = compute_uv ? (full_matrices ? 'A' : 'S') : 'N';
            CBLAS_INT info = 0;

            if constexpr (is_complex_v<T>) {
                long long span = 5LL * minmn + 7;
                if (2LL * maxmn + 2LL * minmn + 1 > span) { span = 2LL * maxmn + 2LL * minmn + 1; }
                CBLAS_INT rwork_len;
                if (!work_size(compute_uv ? minmn * span : 7LL * minmn, &rwork_len)) { return nullptr; }
                ARRAY_HIDDEN(R, rwork, rwork_len);
                lapack::gesdd(job, m, n, a.data<T>(), lda, s.data<R>(), u.data<T>(), ldu, vt.data<T>(), ldvt,
                              work.data<T>(), lwork, rwork.data<R>(), iwork.data<CBLAS_INT>(), &info);
            }
            else {
                lapack::gesdd(job, m, n, a.data<T>(), lda, s.data<R>(), u.data<T>(), ldu, vt.data<T>(), ldvt,
                              work.data<T>(), lwork, iwork.data<CBLAS_INT>(), &info);
            }
            return make_result(u, s, vt, static_cast<long long>(info));
        }


        /**
         * With `lwork = -1` gesdd only writes `work[0]`, so every array it is handed is a null
         * pointer.  `&work` is what picks the flavor's overload.
         */
        template <class T>
        static PyObject *gesdd_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"m", "n", "compute_uv", "full_matrices", nullptr};
            static constexpr Ctx<T> ctx("gesdd_lwork", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, m);
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, compute_uv, 1);     CHECK(compute_uv == 0 || compute_uv == 1, compute_uv);
            SCALAR_OPT(CBLAS_INT, full_matrices, 1);  CHECK(full_matrices == 0 || full_matrices == 1, full_matrices);

            CBLAS_INT minmn = m < n ? m : n;
            CBLAS_INT u0  = compute_uv ? m : 1;
            CBLAS_INT vt0 = compute_uv ? (full_matrices ? n : minmn) : 1;

            char job = compute_uv ? (full_matrices ? 'A' : 'S') : 'N';
            CBLAS_INT lda = std::max<CBLAS_INT>(1, m);
            CBLAS_INT ldu = std::max<CBLAS_INT>(1, u0), ldvt = std::max<CBLAS_INT>(1, vt0);
            T work{};
            CBLAS_INT info = 0;

            if constexpr (is_complex_v<T>) {
                lapack::gesdd(job, m, n, nullptr, lda, nullptr, nullptr, ldu, nullptr, ldvt, &work, -1, nullptr, nullptr, &info);
            }
            else {
                lapack::gesdd(job, m, n, nullptr, lda, nullptr, nullptr, ldu, nullptr, ldvt, &work, -1, nullptr, &info);
            }
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gesvd(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "compute_uv", "full_matrices", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("gesvd", "O|OOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, compute_uv, 1);     CHECK(compute_uv == 0 || compute_uv == 1, compute_uv);
            SCALAR_OPT(CBLAS_INT, full_matrices, 1);  CHECK(full_matrices == 0 || full_matrices == 1, full_matrices);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT minmn = m < n ? m : n, maxmn = m > n ? m : n;
            CBLAS_INT u0  = compute_uv ? m : 1;
            CBLAS_INT u1  = compute_uv ? (full_matrices ? m : minmn) : 1;
            CBLAS_INT vt0 = compute_uv ? (full_matrices ? n : minmn) : 1;
            CBLAS_INT vt1 = compute_uv ? n : 1;
            CBLAS_INT lda = std::max<CBLAS_INT>(1, shape(a, 0));
            CBLAS_INT ldu = std::max<CBLAS_INT>(1, u0), ldvt = std::max<CBLAS_INT>(1, vt0);

            CBLAS_INT default_lwork;
            if constexpr (is_complex_v<T>) {
                if (!work_size(2LL * minmn + maxmn, &default_lwork)) { return nullptr; }
            }
            else {
                long long bound = 3LL * minmn + maxmn;
                if (5LL * minmn > bound) { bound = 5LL * minmn; }
                if (!work_size(bound, &default_lwork)) { return nullptr; }
            }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);

            ARRAY_OUT(R, s, 1, true, ctx.template zeros_as<R>(minmn));
            ARRAY_OUT(T, u, 2, true, ctx.zeros(u0, u1));
            ARRAY_OUT(T, vt, 2, true, ctx.zeros(vt0, vt1));
            ARRAY_HIDDEN(T, work, lwork);

            char job = compute_uv ? (full_matrices ? 'A' : 'S') : 'N';
            CBLAS_INT info = 0;

            if constexpr (is_complex_v<T>) {
                CBLAS_INT rwork_len;
                if (!work_size(5LL * minmn, &rwork_len)) { return nullptr; }
                ARRAY_HIDDEN(R, rwork, rwork_len);
                lapack::gesvd(job, job, m, n, a.data<T>(), lda, s.data<R>(), u.data<T>(), ldu, vt.data<T>(), ldvt, work.data<T>(), lwork, rwork.data<R>(), &info);
            }
            else {
                lapack::gesvd(job, job, m, n, a.data<T>(), lda, s.data<R>(), u.data<T>(), ldu, vt.data<T>(), ldvt, work.data<T>(), lwork, &info);
            }
            return make_result(u, s, vt, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gesvd_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"m", "n", "compute_uv", "full_matrices", nullptr};
            static constexpr Ctx<T> ctx("gesvd_lwork", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, m);
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, compute_uv, 1);     CHECK(compute_uv == 0 || compute_uv == 1, compute_uv);
            SCALAR_OPT(CBLAS_INT, full_matrices, 1);  CHECK(full_matrices == 0 || full_matrices == 1, full_matrices);

            CBLAS_INT minmn = m < n ? m : n;
            CBLAS_INT u0  = compute_uv ? m : 1;
            CBLAS_INT vt0 = compute_uv ? (full_matrices ? n : minmn) : 1;

            char job = compute_uv ? (full_matrices ? 'A' : 'S') : 'N';
            CBLAS_INT lda = std::max<CBLAS_INT>(1, m);
            CBLAS_INT ldu = std::max<CBLAS_INT>(1, u0), ldvt = std::max<CBLAS_INT>(1, vt0);
            T work{};
            CBLAS_INT info = 0;

            if constexpr (is_complex_v<T>) {
                lapack::gesvd(job, job, m, n, nullptr, lda, nullptr, nullptr, ldu, nullptr, ldvt, &work, -1, nullptr, &info);
            }
            else {
                lapack::gesvd(job, job, m, n, nullptr, lda, nullptr, nullptr, ldu, nullptr, ldvt, &work, -1, &info);
            }
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gels(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "trans", "lwork", "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("gels", "OO|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(char, trans, 'N');

            if constexpr (is_complex_v<T>) { CHECK(trans == 'N' || trans == 'C', trans); }
            else                           { CHECK(trans == 'N' || trans == 'T', trans); }

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT minmn = m < n ? m : n, maxmn = m > n ? m : n;

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == maxmn, b);
            CBLAS_INT nrhs = shape(b, 1);
            CBLAS_INT lda = std::max<CBLAS_INT>(1, m), ldb = std::max<CBLAS_INT>(1, maxmn);

            CBLAS_INT default_lwork;
            if (!work_size(1LL * minmn + (minmn > nrhs ? minmn : nrhs), &default_lwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);
            CHECK(lwork >= 1 || lwork == -1, lwork);

            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));

            CBLAS_INT info = 0;
            lapack::gels(trans, m, n, nrhs, a.data<T>(), lda, b.data<T>(), ldb, work.data<T>(), lwork, &info);
            return make_result(a, b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gels_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"m", "n", "nrhs", "trans", nullptr};
            static constexpr Ctx<T> ctx("gels_lwork", "OOO|O", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, m);     CHECK(m >= 0, m);
            SCALAR_REQ(CBLAS_INT, n);     CHECK(n >= 0, n);
            SCALAR_REQ(CBLAS_INT, nrhs);  CHECK(nrhs >= 0, nrhs);
            SCALAR_OPT(char, trans, 'N');
            /* Two checks rather than one with a ternary: CHECK stringizes its expression into the
            * message, and `is_complex_v<T> ? 'C' : 'T'` tells the caller nothing about which
            * letter to pass. */
            if constexpr (is_complex_v<T>) { CHECK(trans == 'N' || trans == 'C', trans); }
            else                           { CHECK(trans == 'N' || trans == 'T', trans); }

            CBLAS_INT maxmn = m > n ? m : n;
            CBLAS_INT lda = std::max<CBLAS_INT>(1, m), ldb = std::max<CBLAS_INT>(1, maxmn);

            T work{};
            CBLAS_INT info = 0;
            lapack::gels(trans, m, n, nrhs, nullptr, lda, nullptr, ldb, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        /**
         * Least squares by SVD.  The singular values and `cond` are real for every flavor, and the
         * complex routines take an extra `rwork` after `lwork`.  `work` is `intent(out)` here rather
         * than hidden, so it comes back in the tuple.
         */
        template <class T>
        static PyObject *gelss(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "cond", "lwork", "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("gelss", "OO|OOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT minmn = m < n ? m : n, maxmn = m > n ? m : n;

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == maxmn, b);
            CBLAS_INT nrhs = shape(b, 1);
            CBLAS_INT lda = std::max<CBLAS_INT>(1, m), ldb = std::max<CBLAS_INT>(1, maxmn);

            SCALAR_OPT(R, cond, -1.0);

            CBLAS_INT default_lwork;
            if constexpr (is_complex_v<T>) {
                if (!work_size(2LL * minmn + (maxmn > nrhs ? maxmn : nrhs), &default_lwork)) {
                    return nullptr;
                }
            }
            else {
                long long inner = 2LL * minmn;
                long long outer = maxmn > nrhs ? maxmn : nrhs;
                if (outer > inner) { inner = outer; }
                if (!work_size(3LL * minmn + inner, &default_lwork)) { return nullptr; }
            }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);
            CHECK(lwork >= 1 || lwork == -1, lwork);

            ARRAY_OUT(R, s, 1, true, ctx.template zeros_as<R>(minmn));
            ARRAY_OUT(T, work, 1, true, ctx.zeros(std::max<CBLAS_INT>(1, lwork)));
            CBLAS_INT rank = 0, info = 0;

            if constexpr (is_complex_v<T>) {
                CBLAS_INT rwork_len;
                if (!work_size(5LL * minmn, &rwork_len)) { return nullptr; }
                ARRAY_HIDDEN(R, rwork, rwork_len);
                lapack::gelss(m, n, nrhs, a.data<T>(), lda, b.data<T>(), ldb, s.data<R>(), cond, &rank, work.data<T>(), lwork, rwork.data<R>(), &info);
            }
            else {
                lapack::gelss(m, n, nrhs, a.data<T>(), lda, b.data<T>(), ldb, s.data<R>(), cond, &rank, work.data<T>(), lwork, &info);
            }
            return make_result(a, b, s, static_cast<long long>(rank), work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gelss_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"m", "n", "nrhs", "cond", "lwork", nullptr};
            static constexpr Ctx<T> ctx("gelss_lwork", "OOO|OO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(CBLAS_INT, m);
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_REQ(CBLAS_INT, nrhs);
            SCALAR_OPT(R, cond, -1.0);
            SCALAR_OPT(CBLAS_INT, lwork, -1);

            CBLAS_INT maxmn = m > n ? m : n;
            CBLAS_INT lda = std::max<CBLAS_INT>(1, m), ldb = std::max<CBLAS_INT>(1, maxmn);
            T work{};
            CBLAS_INT rank = 0, info = 0;

            if constexpr (is_complex_v<T>) {
                lapack::gelss(m, n, nrhs, nullptr, lda, nullptr, ldb, nullptr, cond, &rank, &work, lwork, nullptr, &info);
            }
            else {
                lapack::gelss(m, n, nrhs, nullptr, lda, nullptr, ldb, nullptr, cond, &rank, &work, lwork, &info);
            }
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gelsy(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "jptv", "cond", "lwork", "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("gelsy", "OOOOO|OO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT minmn = m < n ? m : n, maxmn = m > n ? m : n;

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == maxmn, b);
            CBLAS_INT nrhs = shape(b, 1);
            CBLAS_INT lda = std::max<CBLAS_INT>(1, m), ldb = std::max<CBLAS_INT>(1, maxmn);

            ARRAY_INOUT(CBLAS_INT, jptv, 1, true);
            CHECKARRAY(len(jptv) == n, jptv);

            SCALAR_REQ(R, cond);
            SCALAR_REQ(CBLAS_INT, lwork);

            CBLAS_INT minimum_lwork;
            if constexpr (is_complex_v<T>) {
                long long inner = 2LL * minmn;
                if (n + 1LL > inner) { inner = n + 1LL; }
                if (minmn + 1LL * nrhs > inner) { inner = minmn + 1LL * nrhs; }
                if (!work_size(minmn + inner, &minimum_lwork)) { return nullptr; }
            }
            else {
                long long inner = minmn + 3LL * n + 1;
                long long other = 2LL * minmn + nrhs;
                if (other > inner) { inner = other; }
                if (!work_size(inner, &minimum_lwork)) { return nullptr; }
            }
            CHECK(lwork >= minimum_lwork, lwork);

            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            CBLAS_INT rank = 0, info = 0;

            if constexpr (is_complex_v<T>) {
                CBLAS_INT rwork_len;
                if (!work_size(2LL * n, &rwork_len)) { return nullptr; }
                ARRAY_HIDDEN(R, rwork, rwork_len);
                lapack::gelsy(m, n, nrhs, a.data<T>(), lda, b.data<T>(), ldb, jptv.data<CBLAS_INT>(), cond, &rank, work.data<T>(), lwork, rwork.data<R>(), &info);
            }
            else {
                lapack::gelsy(m, n, nrhs, a.data<T>(), lda, b.data<T>(), ldb, jptv.data<CBLAS_INT>(), cond, &rank, work.data<T>(), lwork, &info);
            }
            return make_result(a, b, jptv, static_cast<long long>(rank), static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gelsy_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"m", "n", "nrhs", "cond", "lwork", nullptr};
            static constexpr Ctx<T> ctx("gelsy_lwork", "OOOO|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(CBLAS_INT, m);
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_REQ(CBLAS_INT, nrhs);
            SCALAR_REQ(R, cond);
            SCALAR_OPT(CBLAS_INT, lwork, -1);

            CBLAS_INT maxmn = m > n ? m : n;
            CBLAS_INT lda = std::max<CBLAS_INT>(1, m), ldb = std::max<CBLAS_INT>(1, maxmn);
            T work{};
            CBLAS_INT rank = 0, info = 0;

            if constexpr (is_complex_v<T>) {
                lapack::gelsy(m, n, nrhs, nullptr, lda, nullptr, ldb, nullptr, cond, &rank, &work, lwork, nullptr, &info);
            }
            else {
                lapack::gelsy(m, n, nrhs, nullptr, lda, nullptr, ldb, nullptr, cond, &rank, &work, lwork, &info);
            }
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gelsd(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist_real[] = {"a", "b", "lwork", "size_iwork", "cond",
                                                "overwrite_a", "overwrite_b", nullptr};
            static const char *kwlist_complex[] = {"a", "b", "lwork", "size_rwork", "size_iwork", "cond",
                                                "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("gelsd", is_complex_v<T> ? "OOOOO|OOO" : "OOOO|OOO",
                                        is_complex_v<T> ? kwlist_complex : kwlist_real);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT minmn = m < n ? m : n, maxmn = m > n ? m : n;

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == maxmn, b);
            CBLAS_INT nrhs = shape(b, 1);
            CBLAS_INT lda = std::max<CBLAS_INT>(1, m), ldb = std::max<CBLAS_INT>(1, maxmn);

            SCALAR_REQ(CBLAS_INT, lwork);
            if constexpr (is_complex_v<T>) { CHECK(lwork >= 1 || lwork == -1, lwork); }
            else                           { CHECK(lwork >= 1, lwork); }
            SCALAR_REQ(CBLAS_INT, size_iwork);
            SCALAR_OPT(R, cond, -1.0);

            ARRAY_OUT(R, s, 1, true, ctx.template zeros_as<R>(minmn));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            ARRAY_HIDDEN(CBLAS_INT, iwork, std::max<CBLAS_INT>(1, size_iwork));
            CBLAS_INT rank = 0, info = 0;

            if constexpr (is_complex_v<T>) {
                SCALAR_REQ(CBLAS_INT, size_rwork);
                ARRAY_HIDDEN(R, rwork, std::max<CBLAS_INT>(1, size_rwork));
                lapack::gelsd(m, n, nrhs, a.data<T>(), lda, b.data<T>(), ldb, s.data<R>(), cond, &rank, work.data<T>(), lwork, rwork.data<R>(), iwork.data<CBLAS_INT>(), &info);
            }
            else {
                lapack::gelsd(m, n, nrhs, a.data<T>(), lda, b.data<T>(), ldb, s.data<R>(), cond, &rank, work.data<T>(), lwork, iwork.data<CBLAS_INT>(), &info);
            }
            return make_result(b, s, static_cast<long long>(rank), static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gelsd_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"m", "n", "nrhs", "cond", "lwork", nullptr};
            static constexpr Ctx<T> ctx("gelsd_lwork", "OOO|OO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(CBLAS_INT, m);
            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_REQ(CBLAS_INT, nrhs);
            SCALAR_OPT(R, cond, -1.0);
            SCALAR_OPT(CBLAS_INT, lwork, -1);

            CBLAS_INT maxmn = m > n ? m : n;
            CBLAS_INT lda = std::max<CBLAS_INT>(1, m), ldb = std::max<CBLAS_INT>(1, maxmn);
            T work{};
            CBLAS_INT iwork = 0, rank = 0, info = 0;

            if constexpr (is_complex_v<T>) {
                R rwork{};
                lapack::gelsd(m, n, nrhs, nullptr, lda, nullptr, ldb, nullptr, cond, &rank, &work, lwork, &rwork, &iwork, &info);
                return make_result(work, rwork, static_cast<long long>(iwork), static_cast<long long>(info));
            }
            else {
                lapack::gelsd(m, n, nrhs, nullptr, lda, nullptr, ldb, nullptr, cond, &rank, &work, lwork, &iwork, &info);
                return make_result(work, static_cast<long long>(iwork), static_cast<long long>(info));
            }
        }


        template <class T>
        static PyObject *geqp3(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("geqp3", "O|OO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT minmn = m < n ? m : n, lda = std::max<CBLAS_INT>(1, shape(a, 0));

            CBLAS_INT default_lwork;
            if (!work_size(3LL * (n + 1), &default_lwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);
            CHECK(lwork >= n || lwork == -1, lwork);

            /* The column pivots are LAPACK's own 1-based numbering, unshifted, as in gelsy. */
            ARRAY_OUT(CBLAS_INT, jpvt, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            ARRAY_OUT(T, tau, 1, true, ctx.zeros(minmn));
            ARRAY_OUT(T, work, 1, true, ctx.zeros(std::max<CBLAS_INT>(1, lwork)));
            CBLAS_INT info = 0;

            if constexpr (is_complex_v<T>) {
                CBLAS_INT rwork_len;
                if (!work_size(2LL * n, &rwork_len)) { return nullptr; }
                ARRAY_HIDDEN(R, rwork, rwork_len);
                lapack::geqp3(m, n, a.data<T>(), lda, jpvt.data<CBLAS_INT>(), tau.data<T>(), work.data<T>(), lwork, rwork.data<R>(), &info);
            }
            else {
                lapack::geqp3(m, n, a.data<T>(), lda, jpvt.data<CBLAS_INT>(), tau.data<T>(), work.data<T>(), lwork, &info);
            }
            return make_result(a, jpvt, tau, work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *geqrf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("geqrf", "O|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT minmn = m < n ? m : n, lda = std::max<CBLAS_INT>(1, shape(a, 0));

            CBLAS_INT default_lwork;
            if (!work_size(3LL * n, &default_lwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);
            CHECK(lwork >= n || lwork == -1, lwork);

            ARRAY_OUT(T, tau, 1, true, ctx.zeros(minmn));
            ARRAY_OUT(T, work, 1, true, ctx.zeros(std::max<CBLAS_INT>(1, lwork)));
            CBLAS_INT info = 0;

            lapack::geqrf(m, n, a.data<T>(), lda, tau.data<T>(), work.data<T>(), lwork, &info);
            return make_result(a, tau, work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *geqrf_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"m", "n", nullptr};
            static constexpr Ctx<T> ctx("geqrf_lwork", "OO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, m);  CHECK(m > 0, m);
            SCALAR_REQ(CBLAS_INT, n);  CHECK(n > 0, n);

            T work{};
            CBLAS_INT info = 0;
            lapack::geqrf(m, n, nullptr, std::max<CBLAS_INT>(1, m), nullptr, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *geqrfp(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("geqrfp", "O|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) > 0, a);
            CHECKARRAY(shape(a, 1) > 0, a);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT minmn = m < n ? m : n, lda = std::max<CBLAS_INT>(1, shape(a, 0));

            SCALAR_OPT(CBLAS_INT, lwork, std::max<CBLAS_INT>(1, n));
            CHECK(lwork >= n || lwork == -1, lwork);

            ARRAY_OUT(T, tau, 1, true, ctx.zeros(minmn));
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            CBLAS_INT info = 0;

            lapack::geqrfp(m, n, a.data<T>(), lda, tau.data<T>(), work.data<T>(), lwork, &info);
            return make_result(a, tau, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *geqrfp_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"m", "n", nullptr};
            static constexpr Ctx<T> ctx("geqrfp_lwork", "OO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, m);  CHECK(m > 0, m);
            SCALAR_REQ(CBLAS_INT, n);  CHECK(n > 0, n);

            T work{};
            CBLAS_INT info = 0;
            lapack::geqrfp(m, n, nullptr, std::max<CBLAS_INT>(1, m), nullptr, &work, -1, &info);
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gerqf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("gerqf", "O|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT minmn = m < n ? m : n, lda = std::max<CBLAS_INT>(1, shape(a, 0));

            CBLAS_INT default_lwork;
            if (!work_size(3LL * m, &default_lwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);
            CHECK(lwork >= m || lwork == -1, lwork);

            ARRAY_OUT(T, tau, 1, true, ctx.zeros(minmn));
            ARRAY_OUT(T, work, 1, true, ctx.zeros(std::max<CBLAS_INT>(1, lwork)));
            CBLAS_INT info = 0;

            lapack::gerqf(m, n, a.data<T>(), lda, tau.data<T>(), work.data<T>(), lwork, &info);
            return make_result(a, tau, work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *geequ(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", nullptr};
            static constexpr Ctx<T> ctx("geequ", "O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            ARRAY_IN(T, a, 2);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1), lda = std::max<CBLAS_INT>(1, shape(a, 0));

            ARRAY_OUT(R, r, 1, true, ctx.template zeros_as<R>(m));
            ARRAY_OUT(R, c, 1, true, ctx.template zeros_as<R>(n));
            R rowcnd = 0, colcnd = 0, amax = 0;
            CBLAS_INT info = 0;

            lapack::geequ(m, n, a.data<T>(), lda, r.data<R>(), c.data<R>(), &rowcnd, &colcnd, &amax, &info);
            return make_result(r, c, rowcnd, colcnd, amax, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *geequb(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", nullptr};
            static constexpr Ctx<T> ctx("geequb", "O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            ARRAY_IN(T, a, 2);
            CBLAS_INT m = shape(a, 0), n = shape(a, 1), lda = std::max<CBLAS_INT>(1, shape(a, 0));

            ARRAY_OUT(R, r, 1, true, ctx.template zeros_as<R>(m));
            ARRAY_OUT(R, c, 1, true, ctx.template zeros_as<R>(n));
            R rowcnd = 0, colcnd = 0, amax = 0;
            CBLAS_INT info = 0;

            lapack::geequb(m, n, a.data<T>(), lda, r.data<R>(), c.data<R>(), &rowcnd, &colcnd, &amax, &info);
            return make_result(r, c, rowcnd, colcnd, amax, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *geev(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "compute_vl", "compute_vr", "lwork", "overwrite_a", nullptr};
            static constexpr Ctx<T> ctx("geev", "O|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_OPT(CBLAS_INT, compute_vl, 1);  CHECK(compute_vl == 1 || compute_vl == 0, compute_vl);
            SCALAR_OPT(CBLAS_INT, compute_vr, 1);  CHECK(compute_vr == 1 || compute_vr == 0, compute_vr);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0));

            CBLAS_INT default_lwork, minimum_lwork;
            if constexpr (is_complex_v<T>) {
                if (!work_size(2LL * n, &default_lwork) || !work_size(2LL * n, &minimum_lwork)) {
                    return nullptr;
                }
            }
            else {
                if (!work_size(4LL * n, &default_lwork) ||
                    !work_size((compute_vl || compute_vr) ? 4LL * n : 3LL * n, &minimum_lwork)) {
                    return nullptr;
                }
            }
            SCALAR_OPT(CBLAS_INT, lwork, default_lwork);
            CHECK(lwork >= minimum_lwork, lwork);

            CBLAS_INT vl_rows = compute_vl ? n : 1, ldvl = std::max<CBLAS_INT>(1, vl_rows);
            CBLAS_INT vr_rows = compute_vr ? n : 1, ldvr = std::max<CBLAS_INT>(1, vr_rows);
            ARRAY_HIDDEN(T, work, std::max<CBLAS_INT>(1, lwork));
            CBLAS_INT info = 0;

            if constexpr (is_complex_v<T>) {
                using R = real_of_t<T>;
                CBLAS_INT rwork_len;
                if (!work_size(2LL * n, &rwork_len)) { return nullptr; }
                ARRAY_OUT(T, w, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, vl, 2, true, ctx.zeros(vl_rows, n));
                ARRAY_OUT(T, vr, 2, true, ctx.zeros(vr_rows, n));
                ARRAY_HIDDEN(R, rwork, rwork_len);

                lapack::geev(compute_vl ? 'V' : 'N', compute_vr ? 'V' : 'N', n, a.data<T>(), lda,
                             w.data<T>(), vl.data<T>(), ldvl, vr.data<T>(), ldvr, work.data<T>(),
                             lwork, rwork.data<R>(), &info);
                return make_result(w, vl, vr, static_cast<long long>(info));
            }
            else {
                ARRAY_OUT(T, wr, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, wi, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, vl, 2, true, ctx.zeros(vl_rows, n));
                ARRAY_OUT(T, vr, 2, true, ctx.zeros(vr_rows, n));

                lapack::geev(compute_vl ? 'V' : 'N', compute_vr ? 'V' : 'N', n, a.data<T>(), lda,
                             wr.data<T>(), wi.data<T>(), vl.data<T>(), ldvl, vr.data<T>(), ldvr,
                             work.data<T>(), lwork, &info);
                return make_result(wr, wi, vl, vr, static_cast<long long>(info));
            }
        }


        template <class T>
        static PyObject *geev_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"n", "compute_vl", "compute_vr", nullptr};
            static constexpr Ctx<T> ctx("geev_lwork", "O|OO", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(CBLAS_INT, n);
            SCALAR_OPT(CBLAS_INT, compute_vl, 1);  CHECK(compute_vl == 1 || compute_vl == 0, compute_vl);
            SCALAR_OPT(CBLAS_INT, compute_vr, 1);  CHECK(compute_vr == 1 || compute_vr == 0, compute_vr);

            CBLAS_INT lda = std::max<CBLAS_INT>(1, n);
            CBLAS_INT ldvl = std::max<CBLAS_INT>(1, compute_vl ? n : 1);
            CBLAS_INT ldvr = std::max<CBLAS_INT>(1, compute_vr ? n : 1);
            T work{};
            CBLAS_INT info = 0;

            if constexpr (is_complex_v<T>) {
                real_of_t<T> rwork{};
                lapack::geev(compute_vl ? 'V' : 'N', compute_vr ? 'V' : 'N', n, nullptr, lda, nullptr, nullptr, ldvl, nullptr, ldvr, &work, -1, &rwork, &info);
            }
            else {
                lapack::geev(compute_vl ? 'V' : 'N', compute_vr ? 'V' : 'N', n, nullptr, lda, nullptr, nullptr, nullptr, ldvl, nullptr, ldvr, &work, -1, &info);
            }
            return make_result(work, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *ggev(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "compute_vl", "compute_vr", "lwork",
                                           "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("ggev", "OO|OOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(CBLAS_INT, compute_vl, 1);  CHECK(compute_vl == 1 || compute_vl == 0, compute_vl);
            SCALAR_OPT(CBLAS_INT, compute_vr, 1);  CHECK(compute_vr == 1 || compute_vr == 0, compute_vr);

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0));

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == shape(b, 1), b);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT ldb = std::max<CBLAS_INT>(1, shape(b, 0));

            CBLAS_INT minimum_lwork;
            if (!work_size(is_complex_v<T> ? 2LL * n : 8LL * n, &minimum_lwork)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, lwork, minimum_lwork);
            CHECK(lwork == -1 || lwork >= minimum_lwork, lwork);

            CBLAS_INT vl_rows = compute_vl ? n : 1, ldvl = std::max<CBLAS_INT>(1, vl_rows);
            CBLAS_INT vr_rows = compute_vr ? n : 1, ldvr = std::max<CBLAS_INT>(1, vr_rows);
            ARRAY_OUT(T, vl, 2, true, ctx.zeros(vl_rows, n));
            ARRAY_OUT(T, vr, 2, true, ctx.zeros(vr_rows, n));
            ARRAY_OUT(T, work, 1, true, ctx.zeros(std::max<CBLAS_INT>(1, lwork)));
            ARRAY_OUT(T, beta, 1, true, ctx.zeros(n));
            CBLAS_INT info = 0;

            if constexpr (is_complex_v<T>) {
                using R = real_of_t<T>;
                CBLAS_INT rwork_len;
                if (!work_size(8LL * n, &rwork_len)) { return nullptr; }
                ARRAY_OUT(T, alpha, 1, true, ctx.zeros(n));
                ARRAY_HIDDEN(R, rwork, rwork_len);

                lapack::ggev(compute_vl ? 'V' : 'N', compute_vr ? 'V' : 'N', n, a.data<T>(), lda,
                             b.data<T>(), ldb, alpha.data<T>(), beta.data<T>(), vl.data<T>(), ldvl,
                             vr.data<T>(), ldvr, work.data<T>(), lwork, rwork.data<R>(), &info);
                return make_result(alpha, beta, vl, vr, work, static_cast<long long>(info));
            }
            else {
                ARRAY_OUT(T, alphar, 1, true, ctx.zeros(n));
                ARRAY_OUT(T, alphai, 1, true, ctx.zeros(n));

                lapack::ggev(compute_vl ? 'V' : 'N', compute_vr ? 'V' : 'N', n, a.data<T>(), lda,
                             b.data<T>(), ldb, alphar.data<T>(), alphai.data<T>(), beta.data<T>(),
                             vl.data<T>(), ldvl, vr.data<T>(), ldvr, work.data<T>(), lwork, &info);
                return make_result(alphar, alphai, beta, vl, vr, work,
                                   static_cast<long long>(info));
            }
        }


        template <class T>
        static PyObject *gesvx(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", "fact", "trans", "af", "ipiv", "equed",
                                           "r", "c", "overwrite_a", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("gesvx", "OO|OOOOOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_a);
            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(char, fact, 'E');
            SCALAR_OPT(char, trans, 'N');
            /* The one option letter LAPACK writes back: an input when `fact` is 'F', an output
             * otherwise.  A plain mutable local covers both -- `&equed` goes to LAPACK and the
             * same variable goes to `make_result`, which renders a `char` as one-byte bytes the
             * way f2py's `Py_BuildValue("c")` did. */
            SCALAR_OPT(char, equed, 'B');

            ARRAY_INOUT(T, a, 2, overwrite_a != 0);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0), lda = std::max<CBLAS_INT>(1, shape(a, 0));

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), ldb = std::max<CBLAS_INT>(1, shape(b, 0));

            ARRAY_OUT(T, af, 2, true, ctx.zeros(n, n));
            CHECKARRAY(shape(af, 0) == n && shape(af, 1) == n, af);
            CBLAS_INT ldaf = std::max<CBLAS_INT>(1, shape(af, 0));

            ARRAY_OUT(CBLAS_INT, ipiv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            CHECKARRAY(len(ipiv) == n, ipiv);
            ARRAY_OUT(R, r, 1, true, ctx.template zeros_as<R>(n));
            CHECKARRAY(len(r) == n, r);
            ARRAY_OUT(R, c, 1, true, ctx.template zeros_as<R>(n));
            CHECKARRAY(len(c) == n, c);

            ARRAY_OUT(T, x, 2, true, ctx.zeros(n, nrhs));
            ARRAY_OUT(R, ferr, 1, true, ctx.template zeros_as<R>(nrhs));
            ARRAY_OUT(R, berr, 1, true, ctx.template zeros_as<R>(nrhs));

            CBLAS_INT work_len;
            if (!work_size(is_complex_v<T> ? 2LL * n : 4LL * n, &work_len)) { return nullptr; }
            ARRAY_HIDDEN(T, work, work_len);

            CBLAS_INT ldx = std::max<CBLAS_INT>(1, n);
            ARRAY_HIDDEN(CBLAS_INT, pivots, n);
            CBLAS_INT *shifted = pivots.data<CBLAS_INT>();
            CBLAS_INT *supplied = ipiv.data<CBLAS_INT>();
            R rcond = 0;
            CBLAS_INT info = 0;

            /* Shifted both ways, as `getrs` and `getri` do.  The `.pyf` only decrements after the
             * call, which leaves a supplied `ipiv` one too low for LAPACK -- a `fact='F'` re-solve
             * then returns a wrong answer with `info == 0` -- and decrements the caller's own
             * array on every call.  What is wrong is the asymmetry, not a missing shift:
             * `gtsvx` and `hesvx` are `intent(in,out)` pivots that shift at neither end and
             * are self-consistently 1-based, while the seven call sites that shift on the way
             * out all shift on the way in too.  `gesvx` alone does one without the other. */
            for (CBLAS_INT i = 0; i < n; ++i) { shifted[i] = supplied[i] + 1; }
            if constexpr (is_complex_v<T>) {
                CBLAS_INT rwork_len;
                if (!work_size(2LL * n, &rwork_len)) { return nullptr; }
                ARRAY_HIDDEN(R, rwork, rwork_len);
                lapack::gesvx(fact, trans, n, nrhs, a.data<T>(), lda, af.data<T>(), ldaf, shifted,
                              &equed, r.data<R>(), c.data<R>(), b.data<T>(), ldb, x.data<T>(), ldx,
                              &rcond, ferr.data<R>(), berr.data<R>(), work.data<T>(),
                              rwork.data<R>(), &info);
            }
            else {
                ARRAY_HIDDEN(CBLAS_INT, iwork, n);
                lapack::gesvx(fact, trans, n, nrhs, a.data<T>(), lda, af.data<T>(), ldaf, shifted,
                              &equed, r.data<R>(), c.data<R>(), b.data<T>(), ldb, x.data<T>(), ldx,
                              &rcond, ferr.data<R>(), berr.data<R>(), work.data<T>(),
                              iwork.data<CBLAS_INT>(), &info);
            }
            for (CBLAS_INT i = 0; i < n; ++i) { supplied[i] = shifted[i] - 1; }

            return make_result(a, af, ipiv, equed, r, c, b, x, rcond, ferr, berr,
                               static_cast<long long>(info));
        }


        PyMethodDef gen_methods[] = {
            FAMILY(gees),
            FAMILY(gges),
            FAMILY(gebal),
            FAMILY(gehrd),
            FAMILY(gehrd_lwork),
            FAMILY(gesv),
            FAMILY(gecon),
            FAMILY(getrf),
            FAMILY(getrs),
            FAMILY(getc2),
            FAMILY(gesc2),
            FAMILY(getri),
            FAMILY(getri_lwork),
            FAMILY(gesdd),
            FAMILY(gesdd_lwork),
            FAMILY(gesvd),
            FAMILY(gesvd_lwork),
            FAMILY(gels),
            FAMILY(gels_lwork),
            FAMILY(gelss),
            FAMILY(gelss_lwork),
            FAMILY(gelsy),
            FAMILY(gelsy_lwork),
            FAMILY(gelsd),
            FAMILY(gelsd_lwork),
            FAMILY(geqp3),
            FAMILY(geqrf),
            FAMILY(geqrf_lwork),
            FAMILY(geqrfp),
            FAMILY(geqrfp_lwork),
            FAMILY(gerqf),
            FAMILY(geequ),
            FAMILY(geequb),
            FAMILY(geev),
            FAMILY(geev_lwork),
            FAMILY(ggev),
            FAMILY(gesvx),
            {nullptr, nullptr, 0, nullptr},
        };

    }  // namespace capi
}  // namespace lapack
