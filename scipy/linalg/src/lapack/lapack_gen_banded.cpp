/**
 * @file
 * @brief Python wrappers for the general banded LAPACK routines.
 *
 */
#define PY_ARRAY_UNIQUE_SYMBOL scipy_lapack_ARRAY_API
#define NO_IMPORT_ARRAY
#include "lapack_helpers.hpp"
#include "lapack_calls.hpp"

namespace lapack {
    namespace capi {


        template <class T>
        static PyObject *gbsv(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"kl", "ku", "ab", "b", "overwrite_ab", "overwrite_b",
                                           nullptr};
            static constexpr Ctx<T> ctx("gbsv", "OOOO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_ab);
            SCALAR_FLAG(overwrite_b);
            SCALAR_REQ(CBLAS_INT, kl);
            SCALAR_REQ(CBLAS_INT, ku);

            ARRAY_INOUT(T, ab, 2, overwrite_ab != 0);
            /* Compared in `long long` so a wide `kl` cannot wrap into a match.  `shape(ab, 0)` is
             * already certified to fit CBLAS_INT, so `ldab` taken from it needs no narrowing. */
            CHECKARRAY(2LL * kl + ku + 1 == shape(ab, 0), ab);
            CBLAS_INT n = shape(ab, 1), ldab = shape(ab, 0);

            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), info = 0;

            ARRAY_OUT(CBLAS_INT, piv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            CBLAS_INT *pivots = piv.data<CBLAS_INT>();

            lapack::gbsv(n, kl, ku, nrhs, ab.data<T>(), ldab, pivots, b.data<T>(), n, &info);
            for (CBLAS_INT i = 0; i < n; ++i) { pivots[i] -= 1; }
            return make_result(ab, piv, b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gbtrf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"ab", "kl", "ku", "m", "n", "ldab", "overwrite_ab",
                                           nullptr};
            static constexpr Ctx<T> ctx("gbtrf", "OOO|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_ab);
            ARRAY_INOUT(T, ab, 2, overwrite_ab != 0);
            SCALAR_REQ(CBLAS_INT, kl);
            SCALAR_REQ(CBLAS_INT, ku);

            SCALAR_OPT(CBLAS_INT, m, shape(ab, 1));
            SCALAR_OPT(CBLAS_INT, n, shape(ab, 1));
            SCALAR_OPT(CBLAS_INT, ldab, std::max<CBLAS_INT>(shape(ab, 0), 1));
            CHECK(shape(ab, 0) == ldab, ldab);
            CHECKARRAY(shape(ab, 1) == n, ab);

            /* `m` is the one dimension that is genuinely free: a wide `ab` can be factored as a
             * shorter matrix, and the pivot vector follows `min(m, n)` rather than `n`. */
            CBLAS_INT mn = std::min(m, n), info = 0;
            ARRAY_OUT(CBLAS_INT, ipiv, 1, true, ctx.template zeros_as<CBLAS_INT>(mn));
            CBLAS_INT *pivots = ipiv.data<CBLAS_INT>();

            lapack::gbtrf(m, n, kl, ku, ab.data<T>(), ldab, pivots, &info);
            for (CBLAS_INT i = 0; i < mn; ++i) { pivots[i] -= 1; }
            return make_result(ab, ipiv, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gbtrs(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"ab", "kl", "ku", "b", "ipiv", "trans", "n", "ldab",
                                           "ldb", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("gbtrs", "OOOOO|OOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            ARRAY_IN(T, ab, 2);
            SCALAR_REQ(CBLAS_INT, kl);
            SCALAR_REQ(CBLAS_INT, ku);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);
            ARRAY_IN(CBLAS_INT, ipiv, 1);

            /* Unlike `getrs`, this `.pyf` puts no bound on `trans`: the ternary below maps every
             * value, anything above 1 landing on 'C'. */
            SCALAR_OPT(CBLAS_INT, trans, 0);
            SCALAR_OPT(CBLAS_INT, n, shape(ab, 1));
            SCALAR_OPT(CBLAS_INT, ldab, shape(ab, 0));
            CHECK(shape(ab, 0) == ldab, ldab);
            SCALAR_OPT(CBLAS_INT, ldb, shape(b, 0));
            CHECK(shape(b, 0) == ldb, ldb);
            CHECKARRAY(shape(ab, 1) == n, ab);
            CHECKARRAY(len(ipiv) == n, ipiv);

            CBLAS_INT nrhs = shape(b, 1), info = 0;
            ARRAY_HIDDEN(CBLAS_INT, pivots, n);
            const CBLAS_INT *supplied = ipiv.data<CBLAS_INT>();
            CBLAS_INT *shifted = pivots.data<CBLAS_INT>();
            for (CBLAS_INT i = 0; i < n; ++i) { shifted[i] = supplied[i] + 1; }

            lapack::gbtrs(trans ? (trans == 1 ? 'T' : 'C') : 'N', n, kl, ku, nrhs, ab.data<T>(),
                          ldab, shifted, b.data<T>(), ldb, &info);
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gbcon(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"kl", "ku", "ab", "ipiv", "anorm", "norm", "ldab",
                                           nullptr};
            static constexpr Ctx<T> ctx("gbcon", "OOOOO|OO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(char, norm, '1');
            SCALAR_REQ(CBLAS_INT, kl);  CHECK(kl >= 0, kl);
            SCALAR_REQ(CBLAS_INT, ku);  CHECK(ku >= 0, ku);

            ARRAY_IN(T, ab, 2);
            ARRAY_IN(CBLAS_INT, ipiv, 1);
            SCALAR_REQ(R, anorm);

            /* The factored band needs `kl` rows of fill-in above the matrix's own `kl + ku + 1`.
             * work_size does the checked narrowing; its `>= 1` floor is not reachable here, both
             * band widths having just been checked non-negative. */
            CBLAS_INT band_rows;
            if (!work_size(2LL * kl + ku + 1, &band_rows)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, ldab, band_rows);
            CHECK(ldab >= band_rows, ldab);

            CBLAS_INT n = shape(ab, 1);
            CHECKARRAY(shape(ab, 0) == ldab, ab);
            CHECKARRAY(len(ipiv) == n, ipiv);

            CBLAS_INT work_len;
            if (!work_size(is_complex_v<T> ? 2LL * n : 3LL * n, &work_len)) { return nullptr; }
            ARRAY_HIDDEN(T, work, work_len);

            R rcond = 0;
            CBLAS_INT info = 0;
            ARRAY_HIDDEN(CBLAS_INT, pivots, n);
            const CBLAS_INT *supplied = ipiv.data<CBLAS_INT>();
            CBLAS_INT *shifted = pivots.data<CBLAS_INT>();
            for (CBLAS_INT i = 0; i < n; ++i) { shifted[i] = supplied[i] + 1; }

            if constexpr (is_complex_v<T>) {
                ARRAY_HIDDEN(R, rwork, n);
                lapack::gbcon(norm, n, kl, ku, ab.data<T>(), ldab, shifted, anorm, &rcond,
                              work.data<T>(), rwork.data<R>(), &info);
            }
            else {
                ARRAY_HIDDEN(CBLAS_INT, iwork, n);
                lapack::gbcon(norm, n, kl, ku, ab.data<T>(), ldab, shifted, anorm, &rcond,
                              work.data<T>(), iwork.data<CBLAS_INT>(), &info);
            }
            return make_result(rcond, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *langb(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"norm", "kl", "ku", "ab", "ldab", nullptr};
            static constexpr Ctx<T> ctx("langb", "OOOO|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_REQ(char, norm);
            CHECK(norm == 'M' || norm == 'm' || norm == '1' || norm == 'O' || norm == 'o' ||
                  norm == 'I' || norm == 'i' || norm == 'F' || norm == 'f' || norm == 'E' ||
                  norm == 'e', norm);
            SCALAR_REQ(CBLAS_INT, kl);  CHECK(kl >= 0, kl);
            SCALAR_REQ(CBLAS_INT, ku);  CHECK(ku >= 0, ku);

            ARRAY_IN(T, ab, 2);

            /* No fill-in room here: `langb` only reads the matrix, so the band is `kl + ku + 1`
             * rows rather than the `2*kl + ku + 1` the factoring routines want. */
            CBLAS_INT band_rows;
            if (!work_size(1LL * kl + ku + 1, &band_rows)) { return nullptr; }
            SCALAR_OPT(CBLAS_INT, ldab, band_rows);
            CHECK(ldab >= band_rows, ldab);

            CBLAS_INT n = shape(ab, 1), work_len;
            CHECKARRAY(shape(ab, 0) == ldab, ab);
            if (!work_size(1LL * n + 1, &work_len)) { return nullptr; }
            ARRAY_HIDDEN(R, work, work_len);

            R n2 = lapack::langb(norm, n, kl, ku, ab.data<T>(), ldab, work.data<R>());
            RETURN(n2);
        }


        PyMethodDef gen_banded_methods[] = {
            FAMILY(gbsv),
            FAMILY(gbtrf),
            FAMILY(gbtrs),
            FAMILY(gbcon),
            FAMILY(langb),
            {nullptr, nullptr, 0, nullptr},
        };

    }  // namespace capi
}  // namespace lapack
