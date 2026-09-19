/**
 * @file
 * @brief Python wrappers for the general tridiagonal LAPACK routines.
 *
 */
#define PY_ARRAY_UNIQUE_SYMBOL scipy_lapack_ARRAY_API
#define NO_IMPORT_ARRAY
#include "lapack_helpers.hpp"
#include "lapack_calls.hpp"

namespace lapack {
    namespace capi {


        template <class T>
        static PyObject *gtsv(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"dl", "d", "du", "b", "overwrite_dl", "overwrite_d",
                                           "overwrite_du", "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("gtsv", "OOOO|OOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_dl);
            SCALAR_FLAG(overwrite_d);
            SCALAR_FLAG(overwrite_du);
            SCALAR_FLAG(overwrite_b);

            ARRAY_INOUT(T, dl, 1, overwrite_dl != 0);
            ARRAY_INOUT(T, d, 1, overwrite_d != 0);
            ARRAY_INOUT(T, du, 1, overwrite_du != 0);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);

            CBLAS_INT n = len(d), nm1 = n > 1 ? n - 1 : 0;
            CHECKARRAY(len(dl) == nm1, dl);
            CHECKARRAY(len(du) == nm1, du);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), ldb = std::max<CBLAS_INT>(1, n), info = 0;

            lapack::gtsv(n, nrhs, dl.data<T>(), d.data<T>(), du.data<T>(), b.data<T>(), ldb, &info);
            /* `dl` comes back holding `du2`, the second superdiagonal of U -- LAPACK overwrites it
             * in place, and the `.pyf` renames it on the way out with `out=du2`. */
            return make_result(dl, d, du, b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gttrf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"dl", "d", "du", "overwrite_dl", "overwrite_d",
                                           "overwrite_du", nullptr};
            static constexpr Ctx<T> ctx("gttrf", "OOO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_dl);
            SCALAR_FLAG(overwrite_d);
            SCALAR_FLAG(overwrite_du);

            ARRAY_INOUT(T, dl, 1, overwrite_dl != 0);
            ARRAY_INOUT(T, d, 1, overwrite_d != 0);
            ARRAY_INOUT(T, du, 1, overwrite_du != 0);

            CBLAS_INT n = len(d), nm1 = n > 1 ? n - 1 : 0, nm2 = n > 2 ? n - 2 : 0, info = 0;
            CHECKARRAY(len(dl) == nm1, dl);
            CHECKARRAY(len(du) == nm1, du);

            ARRAY_OUT(T, du2, 1, true, ctx.zeros(nm2));
            /* Handed straight to `gttrs`/`gtsvx`, so it stays in LAPACK's 1-based form.  The
             * `getrf` family 0-bases its pivots; this one never did, and the two consumers of
             * this vector are wrapped here to match. */
            ARRAY_OUT(CBLAS_INT, ipiv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));

            lapack::gttrf(n, dl.data<T>(), d.data<T>(), du.data<T>(), du2.data<T>(),
                          ipiv.data<CBLAS_INT>(), &info);
            return make_result(dl, d, du, du2, ipiv, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gttrs(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"dl", "d", "du", "du2", "ipiv", "b", "trans",
                                           "overwrite_b", nullptr};
            static constexpr Ctx<T> ctx("gttrs", "OOOOOO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);
            SCALAR_OPT(char, trans, 'N');
            CHECK(trans == 'N' || trans == 'T' || trans == 'C', trans);

            ARRAY_IN(T, dl, 1);
            ARRAY_IN(T, d, 1);
            ARRAY_IN(T, du, 1);
            ARRAY_IN(T, du2, 1);
            ARRAY_IN(CBLAS_INT, ipiv, 1);
            ARRAY_INOUT(T, b, 2, overwrite_b != 0);

            CBLAS_INT n = len(d), nm1 = n > 1 ? n - 1 : 0, nm2 = n > 2 ? n - 2 : 0;
            CHECKARRAY(len(dl) == nm1, dl);
            CHECKARRAY(len(du) == nm1, du);
            CHECKARRAY(len(du2) == nm2, du2);
            CHECKARRAY(len(ipiv) == n, ipiv);
            CHECKARRAY(shape(b, 0) == n, b);
            CBLAS_INT nrhs = shape(b, 1), ldb = std::max<CBLAS_INT>(1, shape(b, 0)), info = 0;

            lapack::gttrs(trans, n, nrhs, dl.data<T>(), d.data<T>(), du.data<T>(), du2.data<T>(),
                          ipiv.data<CBLAS_INT>(), b.data<T>(), ldb, &info);
            return make_result(b, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gtcon(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"dl", "d", "du", "du2", "ipiv", "anorm", "norm", nullptr};
            static constexpr Ctx<T> ctx("gtcon", "OOOOOO|O", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(char, norm, '1');

            ARRAY_IN(T, dl, 1);
            ARRAY_IN(T, d, 1);
            ARRAY_IN(T, du, 1);
            ARRAY_IN(T, du2, 1);
            ARRAY_IN(CBLAS_INT, ipiv, 1);
            SCALAR_REQ(R, anorm);

            CBLAS_INT n = len(d), nm1 = n > 1 ? n - 1 : 0, nm2 = n > 2 ? n - 2 : 0;
            CHECKARRAY(len(dl) == nm1, dl);
            CHECKARRAY(len(du) == nm1, du);
            CHECKARRAY(len(du2) == nm2, du2);
            CHECKARRAY(len(ipiv) == n, ipiv);

            CBLAS_INT work_len;
            if (!work_size(2LL * n, &work_len)) { return nullptr; }
            ARRAY_HIDDEN(T, work, work_len);

            R rcond = 0;
            CBLAS_INT info = 0;
            if constexpr (is_complex_v<T>) {
                lapack::gtcon(norm, n, dl.data<T>(), d.data<T>(), du.data<T>(), du2.data<T>(),
                              ipiv.data<CBLAS_INT>(), anorm, &rcond, work.data<T>(), &info);
            }
            else {
                ARRAY_HIDDEN(CBLAS_INT, iwork, n);
                lapack::gtcon(norm, n, dl.data<T>(), d.data<T>(), du.data<T>(), du2.data<T>(),
                              ipiv.data<CBLAS_INT>(), anorm, &rcond, work.data<T>(),
                              iwork.data<CBLAS_INT>(), &info);
            }
            return make_result(rcond, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *gtsvx(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"dl", "d", "du", "b", "fact", "trans", "dlf", "df",
                                           "duf", "du2", "ipiv", nullptr};
            static constexpr Ctx<T> ctx("gtsvx", "OOOO|OOOOOOO", kwlist);
            PARSE_ARGS();

            using R = real_of_t<T>;
            SCALAR_OPT(char, fact, 'N');   CHECK(fact == 'F' || fact == 'N', fact);
            SCALAR_OPT(char, trans, 'N');  CHECK(trans == 'N' || trans == 'C' || trans == 'T', trans);

            ARRAY_IN(T, dl, 1);
            ARRAY_IN(T, d, 1);
            ARRAY_IN(T, du, 1);
            CBLAS_INT n = len(d), nm1 = n > 1 ? n - 1 : 0, nm2 = n > 2 ? n - 2 : 0;
            CHECKARRAY(len(dl) == nm1, dl);
            CHECKARRAY(len(du) == nm1, du);

            ARRAY_IN(T, b, 2);
            CHECKARRAY(shape(b, 0) >= n, b);
            CBLAS_INT nrhs = shape(b, 1), ldb = std::max<CBLAS_INT>(1, shape(b, 0)), ldx = ldb;

            /* The five factorization buffers are `intent(in,out)` without `copy`: supplied when
             * `fact` is 'F' and written in place, filled in and handed back otherwise. */
            ARRAY_OUT(T, dlf, 1, true, ctx.zeros(nm1));   CHECKARRAY(len(dlf) == nm1, dlf);
            ARRAY_OUT(T, df, 1, true, ctx.zeros(n));      CHECKARRAY(len(df) == n, df);
            ARRAY_OUT(T, duf, 1, true, ctx.zeros(nm1));   CHECKARRAY(len(duf) == nm1, duf);
            ARRAY_OUT(T, du2, 1, true, ctx.zeros(nm2));   CHECKARRAY(len(du2) == nm2, du2);
            ARRAY_OUT(CBLAS_INT, ipiv, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            CHECKARRAY(len(ipiv) == n, ipiv);

            ARRAY_OUT(T, x, 2, true, ctx.zeros(ldx, nrhs));
            ARRAY_OUT(R, ferr, 1, true, ctx.template zeros_as<R>(nrhs));
            ARRAY_OUT(R, berr, 1, true, ctx.template zeros_as<R>(nrhs));

            CBLAS_INT work_len;
            if (!work_size(is_complex_v<T> ? 2LL * n : 3LL * n, &work_len)) { return nullptr; }
            ARRAY_HIDDEN(T, work, work_len);

            R rcond = 0;
            CBLAS_INT info = 0;
            if constexpr (is_complex_v<T>) {
                ARRAY_HIDDEN(R, rwork, n);
                lapack::gtsvx(fact, trans, n, nrhs, dl.data<T>(), d.data<T>(), du.data<T>(),
                              dlf.data<T>(), df.data<T>(), duf.data<T>(), du2.data<T>(),
                              ipiv.data<CBLAS_INT>(), b.data<T>(), ldb, x.data<T>(), ldx, &rcond,
                              ferr.data<R>(), berr.data<R>(), work.data<T>(), rwork.data<R>(), &info);
            }
            else {
                ARRAY_HIDDEN(CBLAS_INT, iwork, n);
                lapack::gtsvx(fact, trans, n, nrhs, dl.data<T>(), d.data<T>(), du.data<T>(),
                              dlf.data<T>(), df.data<T>(), duf.data<T>(), du2.data<T>(),
                              ipiv.data<CBLAS_INT>(), b.data<T>(), ldb, x.data<T>(), ldx, &rcond,
                              ferr.data<R>(), berr.data<R>(), work.data<T>(),
                              iwork.data<CBLAS_INT>(), &info);
            }
            return make_result(dlf, df, duf, du2, ipiv, x, rcond, ferr, berr,
                               static_cast<long long>(info));
        }


        template <class T>
        static PyObject *stev(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"d", "e", "compute_v", "overwrite_d", "overwrite_e",
                                           nullptr};
            static constexpr Ctx<T> ctx("stev", "OO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_d);
            SCALAR_FLAG(overwrite_e);
            SCALAR_OPT(CBLAS_INT, compute_v, 1);

            ARRAY_INOUT(T, d, 1, overwrite_d != 0);
            CBLAS_INT n = len(d);
            CHECKARRAY(n > 0, d);
            /* LAPACK destroys `e`; the caller never gets it back, so this is `intent(in,copy)` --
             * an ARRAY_INOUT that is absent from the result. */
            ARRAY_INOUT(T, e, 1, overwrite_e != 0);
            CHECKARRAY(len(e) == std::max<CBLAS_INT>(n - 1, 1), e);

            CBLAS_INT ldz = compute_v ? n : 1, work_len = 1, info = 0;
            if (compute_v && !work_size(2LL * n - 2, &work_len)) { return nullptr; }
            ARRAY_OUT(T, z, 2, true, ctx.zeros(ldz, ldz));
            ARRAY_HIDDEN(T, work, work_len);

            lapack::stev(compute_v ? 'V' : 'N', n, d.data<T>(), e.data<T>(), z.data<T>(), ldz,
                         work.data<T>(), &info);
            return make_result(d, z, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *stebz(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"d", "e", "range", "vl", "vu", "il", "iu", "tol",
                                           "order", nullptr};
            static constexpr Ctx<T> ctx("stebz", "OOOOOOOOO", kwlist);
            PARSE_ARGS();

            ARRAY_IN(T, d, 1);
            CBLAS_INT n = len(d);
            CHECKARRAY(n > 0, d);
            ARRAY_IN(T, e, 1);
            CHECKARRAY(len(e) == n - 1, e);

            SCALAR_REQ(CBLAS_INT, range);
            SCALAR_REQ(T, vl);
            SCALAR_REQ(T, vu);
            SCALAR_REQ(CBLAS_INT, il);
            SCALAR_REQ(CBLAS_INT, iu);
            SCALAR_REQ(T, tol);
            SCALAR_REQ(char, order);

            CBLAS_INT work_len, iwork_len;
            if (!work_size(4LL * n, &work_len) || !work_size(3LL * n, &iwork_len)) {
                return nullptr;
            }
            ARRAY_HIDDEN(T, work, work_len);
            ARRAY_HIDDEN(CBLAS_INT, iwork, iwork_len);

            ARRAY_OUT(T, w, 1, true, ctx.zeros(n));
            ARRAY_OUT(CBLAS_INT, iblock, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            ARRAY_OUT(CBLAS_INT, isplit, 1, true, ctx.template zeros_as<CBLAS_INT>(n));
            CBLAS_INT m = 0, nsplit = 0, info = 0;

            lapack::stebz(range > 0 ? (range == 1 ? 'V' : 'I') : 'A', order, n, vl, vu, il, iu, tol,
                          d.data<T>(), e.data<T>(), &m, &nsplit, w.data<T>(),
                          iblock.data<CBLAS_INT>(), isplit.data<CBLAS_INT>(), work.data<T>(),
                          iwork.data<CBLAS_INT>(), &info);
            return make_result(static_cast<long long>(m), w, iblock, isplit,
                               static_cast<long long>(info));
        }


        template <class T>
        static PyObject *sterf(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"d", "e", "overwrite_d", "overwrite_e", nullptr};
            static constexpr Ctx<T> ctx("sterf", "OO|OO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_d);
            SCALAR_FLAG(overwrite_e);

            ARRAY_INOUT(T, d, 1, overwrite_d != 0);
            CBLAS_INT n = len(d), nm1 = n > 1 ? n - 1 : 0, info = 0;
            ARRAY_INOUT(T, e, 1, overwrite_e != 0);
            CHECKARRAY(len(e) == nm1, e);

            lapack::sterf(n, d.data<T>(), e.data<T>(), &info);
            return make_result(d, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *stein(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"d", "e", "w", "iblock", "isplit", nullptr};
            static constexpr Ctx<T> ctx("stein", "OOOOO", kwlist);
            PARSE_ARGS();

            ARRAY_IN(T, d, 1);
            CBLAS_INT n = len(d);
            CHECKARRAY(n > 0, d);
            ARRAY_IN(T, e, 1);
            CHECKARRAY(len(e) == n - 1, e);
            ARRAY_IN(T, w, 1);
            ARRAY_IN(CBLAS_INT, iblock, 1);
            CHECKARRAY(len(iblock) == n, iblock);
            ARRAY_IN(CBLAS_INT, isplit, 1);
            CHECKARRAY(len(isplit) == n, isplit);

            CBLAS_INT m = len(w), ldz = n, work_len, info = 0;
            if (!work_size(5LL * n, &work_len)) { return nullptr; }
            ARRAY_OUT(T, z, 2, true, ctx.zeros(ldz, m));
            ARRAY_HIDDEN(T, work, work_len);
            ARRAY_HIDDEN(CBLAS_INT, iwork, n);
            ARRAY_HIDDEN(CBLAS_INT, ifail, m);

            lapack::stein(n, d.data<T>(), e.data<T>(), m, w.data<T>(), iblock.data<CBLAS_INT>(),
                          isplit.data<CBLAS_INT>(), z.data<T>(), ldz, work.data<T>(),
                          iwork.data<CBLAS_INT>(), ifail.data<CBLAS_INT>(), &info);
            return make_result(z, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *stemr(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"d", "e", "range", "vl", "vu", "il", "iu", "compute_v",
                                           "lwork", "liwork", "overwrite_d", "overwrite_e", nullptr};
            static constexpr Ctx<T> ctx("stemr", "OOOOOOO|OOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_d);
            SCALAR_FLAG(overwrite_e);
            SCALAR_OPT(CBLAS_INT, compute_v, 1);

            ARRAY_INOUT(T, d, 1, overwrite_d != 0);
            CBLAS_INT n = len(d);
            CHECKARRAY(n > 0, d);
            /* `intent(in,copy)`, as every other routine here that LAPACK documents as destroying
             * its `e` declares it (`stev`, `sterf`, `stevd`, and this routine's own `_lwork`
             * twin).  The `.pyf` says plain `intent(in)` on this one line alone, which hands
             * LAPACK the caller's buffer to overwrite -- ?STEMR takes `e` as [in,out] and uses
             * `e(n)` as scratch, which is why it is `n` long and not `n - 1`. */
            ARRAY_INOUT(T, e, 1, overwrite_e != 0);
            CHECKARRAY(len(e) == n, e);

            SCALAR_REQ(CBLAS_INT, range);
            SCALAR_REQ(T, vl);
            SCALAR_REQ(T, vu);
            SCALAR_REQ(CBLAS_INT, il);
            SCALAR_REQ(CBLAS_INT, iu);

            CBLAS_INT minimum_lwork, minimum_liwork, isuppz_len;
            if (!work_size(compute_v ? 18LL * n : 12LL * n, &minimum_lwork) ||
                !work_size(compute_v ? 10LL * n : 8LL * n, &minimum_liwork) ||
                !work_size(compute_v ? 2LL * n : 1LL, &isuppz_len)) {
                return nullptr;
            }
            SCALAR_OPT(CBLAS_INT, lwork, minimum_lwork);
            CHECK(lwork >= minimum_lwork, lwork);
            SCALAR_OPT(CBLAS_INT, liwork, minimum_liwork);
            CHECK(liwork >= minimum_liwork, liwork);

            ARRAY_OUT(T, w, 1, true, ctx.zeros(n));
            ARRAY_OUT(T, z, 2, true, ctx.zeros(n, n));
            ARRAY_HIDDEN(CBLAS_INT, isuppz, isuppz_len);
            ARRAY_HIDDEN(T, work, lwork);
            ARRAY_HIDDEN(CBLAS_INT, iwork, liwork);
            CBLAS_INT ldz = compute_v ? n : 1, m = 0, tryrac = 1, info = 0;

            lapack::stemr(compute_v ? 'V' : 'N', range > 0 ? (range == 1 ? 'V' : 'I') : 'A', n,
                          d.data<T>(), e.data<T>(), vl, vu, il, iu, &m, w.data<T>(), z.data<T>(),
                          ldz, n, isuppz.data<CBLAS_INT>(), &tryrac, work.data<T>(), lwork,
                          iwork.data<CBLAS_INT>(), liwork, &info);
            return make_result(static_cast<long long>(m), w, z, static_cast<long long>(info));
        }


        template <class T>
        static PyObject *stemr_lwork(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"d", "e", "range", "vl", "vu", "il", "iu", "compute_v",
                                           "overwrite_d", "overwrite_e", nullptr};
            static constexpr Ctx<T> ctx("stemr_lwork", "OOOOOOO|OOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_d);
            SCALAR_FLAG(overwrite_e);
            SCALAR_OPT(CBLAS_INT, compute_v, 1);

            ARRAY_INOUT(T, d, 1, overwrite_d != 0);
            CBLAS_INT n = len(d);
            CHECKARRAY(n > 0, d);
            ARRAY_INOUT(T, e, 1, overwrite_e != 0);
            CHECKARRAY(len(e) == n, e);

            SCALAR_REQ(CBLAS_INT, range);
            SCALAR_REQ(T, vl);
            SCALAR_REQ(T, vu);
            SCALAR_REQ(CBLAS_INT, il);
            SCALAR_REQ(CBLAS_INT, iu);

            /* ?STEMR returns from the query having written only WORK(1) and IWORK(1), and having
             * read D and E (through ?LARRC, when `range` selects an interval).  W and ISUPPZ it
             * never touches; Z it writes only under a separate NZC == -1 query, and NZC is `n`
             * here.  So those three are the nullptrs below. */
            T work = 0;
            CBLAS_INT iwork = 0, ldz = compute_v ? n : 1, m = 0, tryrac = 1, info = 0;

            lapack::stemr(compute_v ? 'V' : 'N', range > 0 ? (range == 1 ? 'V' : 'I') : 'A', n,
                          d.data<T>(), e.data<T>(), vl, vu, il, iu, &m, nullptr, nullptr, ldz, n,
                          nullptr, &tryrac, &work, -1, &iwork, -1, &info);
            return make_result(work, static_cast<long long>(iwork), static_cast<long long>(info));
        }


        template <class T>
        static PyObject *stevd(PyObject *Py_UNUSED(self), PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"d", "e", "compute_v", "lwork", "liwork", "overwrite_d",
                                           "overwrite_e", nullptr};
            static constexpr Ctx<T> ctx("stevd", "OO|OOOOO", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_d);
            SCALAR_FLAG(overwrite_e);
            SCALAR_OPT(CBLAS_INT, compute_v, 1);

            ARRAY_INOUT(T, d, 1, overwrite_d != 0);
            CBLAS_INT n = len(d);
            CHECKARRAY(n > 0, d);
            ARRAY_INOUT(T, e, 1, overwrite_e != 0);
            CHECKARRAY(len(e) == std::max<CBLAS_INT>(n - 1, 1), e);

            /* The divide-and-conquer workspaces are exact formulas, so `stevd` needs no `_lwork`
             * companion: 1 + 4n + n^2 reals and 3 + 5n integers with eigenvectors, 1 without. */
            CBLAS_INT minimum_lwork, minimum_liwork;
            if (!work_size(compute_v ? 1LL + 4LL * n + 1LL * n * n : 1LL, &minimum_lwork) ||
                !work_size(compute_v ? 3LL + 5LL * n : 1LL, &minimum_liwork)) {
                return nullptr;
            }
            SCALAR_OPT(CBLAS_INT, lwork, minimum_lwork);
            CHECK(lwork >= minimum_lwork, lwork);
            SCALAR_OPT(CBLAS_INT, liwork, minimum_liwork);
            CHECK(liwork >= minimum_liwork, liwork);

            CBLAS_INT ldz = compute_v ? n : 1, info = 0;
            ARRAY_OUT(T, z, 2, true, ctx.zeros(ldz, ldz));
            ARRAY_HIDDEN(T, work, lwork);
            ARRAY_HIDDEN(CBLAS_INT, iwork, liwork);

            lapack::stevd(compute_v ? 'V' : 'N', n, d.data<T>(), e.data<T>(), z.data<T>(), ldz,
                          work.data<T>(), lwork, iwork.data<CBLAS_INT>(), liwork, &info);
            return make_result(d, z, static_cast<long long>(info));
        }


        PyMethodDef gen_tri_methods[] = {
            FAMILY(gtsv),
            FAMILY(gttrf),
            FAMILY(gttrs),
            FAMILY(gtcon),
            FAMILY(gtsvx),
            ROW(sstev, stev, f32),
            ROW(dstev, stev, f64),
            ROW(sstebz, stebz, f32),
            ROW(dstebz, stebz, f64),
            ROW(ssterf, sterf, f32),
            ROW(dsterf, sterf, f64),
            ROW(sstein, stein, f32),
            ROW(dstein, stein, f64),
            ROW(sstemr, stemr, f32),
            ROW(dstemr, stemr, f64),
            ROW(sstemr_lwork, stemr_lwork, f32),
            ROW(dstemr_lwork, stemr_lwork, f64),
            ROW(sstevd, stevd, f32),
            ROW(dstevd, stevd, f64),
            {nullptr, nullptr, 0, nullptr},
        };

    }  // namespace capi
}  // namespace lapack
