/**
 * @file
 * @brief Python wrappers for the positive definite tridiagonal LAPACK routines.
 *
 */
#define PY_ARRAY_UNIQUE_SYMBOL scipy_lapack_ARRAY_API
#define NO_IMPORT_ARRAY
#include "lapack_helpers.hpp"
#include "lapack_calls.hpp"

namespace lapack {
    namespace capi {


        template <class T>
        static PyObject *ptsv(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"d", "e", "b", "overwrite_d", "overwrite_e",
                                           "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("ptsv", "OOO|OOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_d);
            SCALAR_FLAG(overwrite_e);
            SCALAR_FLAG(overwrite_b);

            ARRAY_INOUT(R, d, 1, overwrite_d != 0);
            ARRAY_INOUT(T, e, 1, overwrite_e != 0);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);

            CBLAS_INT n = len(d), nm1 = n > 0 ? n - 1 : 0;
            CHECKARRAY(len(e) == nm1, e);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), info = 0;

            lapack::ptsv(n, nrhs, d.data<R>(), e.data<T>(), b.data<T>(), n, &info);
            /* `e` comes back as `du`: LAPACK overwrites it with the superdiagonal of U. */
            return make_result(d, e, b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pttrf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"d", "e", "overwrite_d", "overwrite_e", nullptr};
            static constexpr Ctx<T> ctx("pttrf", "OO|OO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_d);
            SCALAR_FLAG(overwrite_e);

            ARRAY_INOUT(R, d, 1, overwrite_d != 0);
            ARRAY_INOUT(T, e, 1, overwrite_e != 0);

            CBLAS_INT n = len(d), nm1 = n > 0 ? n - 1 : 0, info = 0;
            CHECKARRAY(len(e) == nm1, e);

            lapack::pttrf(n, d.data<R>(), e.data<T>(), &info);
            return make_result(d, e, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pttrs(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            /* The complex flavors take `lower`, saying whether `e` is the sub- or the
             * superdiagonal of the Hermitian factor; for a symmetric one there is nothing to
             * distinguish, so the real flavors have no such argument at all. */
            static const char *kwlist_real[] = {"d", "e", "b", "overwrite_b", nullptr};
            static const char *kwlist_complex[] = {"d", "e", "b", "lower", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("pttrs", "OOO|OO", is_complex_v<T> ? kwlist_complex : kwlist_real);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_b);

            ARRAY_IN(R, d, 1);
            ARRAY_IN(T, e, 1);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);

            CBLAS_INT n = len(d), nm1 = n > 0 ? n - 1 : 0;
            CHECKARRAY(len(e) == nm1, e);
            CBLAS_INT nrhs = shape(b, 1), ldb = std::max<CBLAS_INT>(1, shape(b, 0)), info = 0;

            if constexpr (is_complex_v<T>) {
                SCALAR_OPT(CBLAS_INT, lower, 0);
                CHECK(lower == 0 || lower == 1, lower);
                lapack::pttrs(lower ? 'L' : 'U', n, nrhs, d.data<R>(), e.data<T>(), b.data<T>(), ldb, &info);
            }
            else {
                lapack::pttrs(n, nrhs, d.data<R>(), e.data<T>(), b.data<T>(), ldb, &info);
            }
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *pteqr(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"d", "e", "z", "compute_z", "overwrite_d", "overwrite_e", "overwrite_z", nullptr};
            static constexpr Ctx<T> ctx("pteqr", "OOO|OOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_FLAG(overwrite_d);
            SCALAR_FLAG(overwrite_e);
            SCALAR_FLAG(overwrite_z);
            SCALAR_OPT(CBLAS_INT, compute_z, 0);
            CHECK(compute_z >= 0 && compute_z <= 2, compute_z);

            ARRAY_INOUT(R, d, 1, overwrite_d != 0);
            ARRAY_INOUT(R, e, 1, overwrite_e != 0);
            ARRAY_INOUT(T, z, 2, overwrite_z != 0);

            CBLAS_INT n = len(d), nm1 = n > 0 ? n - 1 : 0, info = 0, work_len;
            CHECKARRAY(len(e) == nm1, e);
            /* `z` is required even when it is unused: with `compute_z == 0` LAPACK never reads
             * it and whatever shape the caller passed stands, so only the other two modes
             * constrain it to `(max(1, n), n)`. */
            if (compute_z != 0) {
                CHECKARRAY(shape(z, 0) == std::max<CBLAS_INT>(1, n) && shape(z, 1) == n, z);
            }
            CBLAS_INT ldz = compute_z == 0 ? 1 : std::max<CBLAS_INT>(1, shape(z, 0));
            if (!work_size(4LL * n, &work_len)) { return nullptr; }
            ARRAY_HIDDEN(R, work, work_len);

            lapack::pteqr(compute_z ? (compute_z == 2 ? 'I' : 'V') : 'N', n, d.data<R>(), e.data<R>(), z.data<T>(), ldz, work.data<R>(), &info);
            return make_result(d, e, z, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *ptsvx(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"d", "e", "b", "fact", "df", "ef", nullptr};
            static constexpr Ctx<T> ctx("ptsvx", "OOO|OOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(char, fact, 'N');

            ARRAY_IN(R, d, 1);
            ARRAY_IN(T, e, 1);
            CBLAS_INT n = len(d), nm1 = n > 0 ? n - 1 : 0;
            CHECKARRAY(len(e) == nm1, e);

            ARRAY_IN(T, b, 2);
            CHECKARRAY(shape(b, 0) >= n, b);
            CBLAS_INT nrhs = shape(b, 1), ldb = std::max<CBLAS_INT>(1, shape(b, 0)), ldx = n;

            /* `intent(in,out)` without `copy`: supplied and written in place when `fact` is
             * 'F', filled in and handed back otherwise. */
            ARRAY_OUT(R, df, 1, true, ctx.template zeros_as<R>(n));   CHECKARRAY(len(df) == n, df);
            ARRAY_OUT(T, ef, 1, true, ctx.zeros(nm1));                CHECKARRAY(len(ef) == nm1, ef);

            ARRAY_OUT(T, x, 2, true, ctx.zeros(ldx, nrhs));
            ARRAY_OUT(R, ferr, 1, true, ctx.template zeros_as<R>(nrhs));
            ARRAY_OUT(R, berr, 1, true, ctx.template zeros_as<R>(nrhs));

            CBLAS_INT work_len, info = 0;
            if (!work_size(is_complex_v<T> ? 1LL * n : 2LL * n, &work_len)) { return nullptr; }
            ARRAY_HIDDEN(T, work, work_len);
            R rcond = 0;

            if constexpr (is_complex_v<T>) {
                ARRAY_HIDDEN(R, rwork, n);
                lapack::ptsvx(fact, n, nrhs, d.data<R>(), e.data<T>(), df.data<R>(), ef.data<T>(),
                              b.data<T>(), ldb, x.data<T>(), ldx, &rcond, ferr.data<R>(),
                              berr.data<R>(), work.data<T>(), rwork.data<R>(), &info);
            }
            else {
                lapack::ptsvx(fact, n, nrhs, d.data<R>(), e.data<T>(), df.data<R>(), ef.data<T>(),
                              b.data<T>(), ldb, x.data<T>(), ldx, &rcond, ferr.data<R>(),
                              berr.data<R>(), work.data<T>(), &info);
            }
            return make_result(df, ef, x, rcond, ferr, berr, static_cast<long long>(info));
        }


        PyMethodDef pos_def_tri_methods[] = {
            FAMILY(ptsv),
            FAMILY(pttrf),
            FAMILY(pttrs),
            FAMILY(pteqr),
            FAMILY(ptsvx),
            {nullptr, nullptr, 0, nullptr},
        };

    }  // namespace capi
}  // namespace lapack
