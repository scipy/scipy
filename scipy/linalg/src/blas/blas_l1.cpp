/**
 * @file
 * @brief Level-1 BLAS wrappers.
 *
 * Wrapper conventions and the vocabulary macros (PARSE_ARGS, the SCALAR_ and ARRAY_ argument
 * declarations, CHECK, RETURN, ...) are documented in `blas_helpers.hpp`.  This file contributes
 * the method-table chunk `blas::capi::l1_methods`, merged into the module by `_blas_module.cpp`.
 *
 */
#define PY_ARRAY_UNIQUE_SYMBOL scipy_blas_ARRAY_API
#define NO_IMPORT_ARRAY
#include "blas_helpers.hpp"
#include "blas_calls.hpp"


// We use a nested namespace to avoid polluting the global namespace with the helper functions
// such as len(), shape() and abs() which are used for parameter checks.
namespace blas{
    namespace capi {


        template <class T>
        static PyObject *axpy(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"x", "y", "n", "a", "offx", "incx", "offy", "incy", nullptr};
            static constexpr Ctx<T> ctx("axpy", "OO|OOOOOO", kwlist);
            PARSE_ARGS();

            ARRAY_IN(x, 1);
            ARRAY_INOUT(y, 1, true);

            SCALAR_OPT(T, a, T(1));
            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx != 0, incx);
            SCALAR_OPT(CBLAS_INT, incy, 1);  CHECK(incy != 0, incy);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);
            SCALAR_OPT(CBLAS_INT, offy, 0);  CHECK(offy >= 0 && offy < len(y), offy);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - offx) / abs(incx));
            CHECK(len(y) - offy > (n - 1) * abs(incy), n);
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            blas::axpy(n, a, x.data<T>() + offx, incx, y.data<T>() + offy, incy);
            RETURN(y);   /* out=z */
        }


        template <class T>
        static PyObject *nrm2(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"x", "n", "offx", "incx", nullptr};
            /* tchar_fn prefix: snrm2/dnrm2/scnrm2/dznrm2. PARSE_ARGS uses ctx uniformly. */
            static constexpr Ctx<T> ctx(tchar_fn<T>(), "nrm2", "O|OOO", kwlist);
            PARSE_ARGS();

            ARRAY_IN(x, 1);

            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx > 0, incx);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - offx) / abs(incx));
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            real_of_t<T> n2 = blas::nrm2(n, x.data<T>() + offx, incx);
            RETURN(n2);   /* out=n2 */
        }


        template <class T>
        static PyObject *asum(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"x", "n", "offx", "incx", nullptr};
            static constexpr Ctx<T> ctx(tchar_fn<T>(), "asum", "O|OOO", kwlist);
            PARSE_ARGS();

            ARRAY_IN(x, 1);

            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx != 0, incx);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - offx) / abs(incx));
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            real_of_t<T> s = blas::asum(n, x.data<T>() + offx, incx);
            RETURN(s);   /* out=s */
        }


        template <class T>
        static PyObject *iamax(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"x", "n", "offx", "incx", nullptr};
            static constexpr Ctx<T> ctx(iflavor<T>(), "amax", "O|OOO", kwlist);
            PARSE_ARGS();

            ARRAY_IN(x, 1);

            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx != 0, incx);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - offx) / abs(incx));
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            /* Fortran's 1-based index is shifted to 0-based, as the .pyf callstatement did.
             * long long (not CBLAS_INT) so RETURN resolves to result_item(long long) -> PyLong. */
            long long idx = blas::iamax(n, x.data<T>() + offx, incx) - 1;
            RETURN(idx);   /* out=k */
        }


        template <class T>
        static PyObject *swap(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"x", "y", "n", "offx", "incx", "offy", "incy", nullptr};
            static constexpr Ctx<T> ctx("swap", "OO|OOOOO", kwlist);
            PARSE_ARGS();

            ARRAY_INOUT(x, 1, true);
            ARRAY_INOUT(y, 1, true);

            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx != 0, incx);
            SCALAR_OPT(CBLAS_INT, incy, 1);  CHECK(incy != 0, incy);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);
            SCALAR_OPT(CBLAS_INT, offy, 0);  CHECK(offy >= 0 && offy < len(y), offy);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - offx) / abs(incx));
            CHECK(len(y) - offy > (n - 1) * abs(incy), n);
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            blas::swap(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy);
            RETURN(x, y);   /* out=(x, y) */
        }


        template <class T>
        static PyObject *copy(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"x", "y", "n", "offx", "incx", "offy", "incy", nullptr};
            static constexpr Ctx<T> ctx("copy", "OO|OOOOO", kwlist);
            PARSE_ARGS();

            ARRAY_IN(x, 1);
            ARRAY_INOUT(y, 1, true);

            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx != 0, incx);
            SCALAR_OPT(CBLAS_INT, incy, 1);  CHECK(incy != 0, incy);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);
            SCALAR_OPT(CBLAS_INT, offy, 0);  CHECK(offy >= 0 && offy < len(y), offy);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - offx) / abs(incx));
            CHECK(len(y) - offy > (n - 1) * abs(incy), n);
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            blas::copy(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy);
            RETURN(y);   /* out=y */
        }


        template <class T>
        static PyObject *scal(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "x", "n", "offx", "incx", nullptr};
            static constexpr Ctx<T> ctx("scal", "OO|OOO", kwlist);
            PARSE_ARGS();

            /* f2py processes a before x: a bad scalar wins over a bad array */
            SCALAR_REQ(T, a);
            ARRAY_INOUT(x, 1, true);

            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx != 0, incx);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - offx) / abs(incx));
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            blas::scal(n, a, x.data<T>() + offx, incx);
            RETURN(x);   /* out=x */
        }


        /**
         * csscal/zdscal: scal for a *real* scale factor on complex data (a is the real type A,
         * not the data type T).  Unlike the regular scal, the .pyf declares x intent(in,out,copy),
         * so these carry an overwrite_x flag and copy by default.  Registered as
         * csscal = scal_real<c64, f32>, zdscal = scal_real<c128, f64>.
         * It is a distinct name rather than a second scal overload which otherwise works except
         * MSVC cannot resolve `scal<f32>` to a single function for the PyCFunction cast when
         * two function templates share the name.
         * */
        template <class T, class A>
        static PyObject *scal_real(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "x", "n", "offx", "incx", "overwrite_x", nullptr};
            static constexpr Ctx<T> ctx(tchar<T>(), "scal", "OO|OOOi", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_x);
            SCALAR_REQ(A, a);
            ARRAY_INOUT(x, 1, overwrite_x != 0);

            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx != 0, incx);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - offx) / abs(incx));
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            blas::scal(n, a, x.data<T>() + offx, incx);
            RETURN(x);   /* out=x */
        }


        template <class T>
        static PyObject *dot(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"x", "y", "n", "offx", "incx", "offy", "incy", nullptr};
            static constexpr Ctx<T> ctx("dot", "OO|OOOOO", kwlist);
            PARSE_ARGS();

            ARRAY_IN(x, 1);
            ARRAY_IN(y, 1);

            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx != 0, incx);
            SCALAR_OPT(CBLAS_INT, incy, 1);  CHECK(incy != 0, incy);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);
            SCALAR_OPT(CBLAS_INT, offy, 0);  CHECK(offy >= 0 && offy < len(y), offy);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - offx) / abs(incx));
            CHECK(len(y) - offy > (n - 1) * abs(incy), n);
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            real_of_t<T> xy = blas::dot(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy);
            RETURN(xy);   /* out=xy */
        }


        template <class T>
        static PyObject *dotu(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"x", "y", "n", "offx", "incx", "offy", "incy", nullptr};
            static constexpr Ctx<T> ctx("dotu", "OO|OOOOO", kwlist);
            PARSE_ARGS();

            ARRAY_IN(x, 1);
            ARRAY_IN(y, 1);

            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx != 0, incx);
            SCALAR_OPT(CBLAS_INT, incy, 1);  CHECK(incy != 0, incy);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);
            SCALAR_OPT(CBLAS_INT, offy, 0);  CHECK(offy >= 0 && offy < len(y), offy);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - offx) / abs(incx));
            CHECK(len(y) - offy > (n - 1) * abs(incy), n);
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            T xy = blas::dotu(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy);
            RETURN(xy);   /* out=xy */
        }


        template <class T>
        static PyObject *dotc(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"x", "y", "n", "offx", "incx", "offy", "incy", nullptr};
            static constexpr Ctx<T> ctx("dotc", "OO|OOOOO", kwlist);
            PARSE_ARGS();

            ARRAY_IN(x, 1);
            ARRAY_IN(y, 1);

            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx != 0, incx);
            SCALAR_OPT(CBLAS_INT, incy, 1);  CHECK(incy != 0, incy);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);
            SCALAR_OPT(CBLAS_INT, offy, 0);  CHECK(offy >= 0 && offy < len(y), offy);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - offx) / abs(incx));
            CHECK(len(y) - offy > (n - 1) * abs(incy), n);
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            T xy = blas::dotc(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy);
            RETURN(xy);   /* out=xy */
        }


        template <class T>
        static PyObject *rotg(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"a", "b", nullptr};
            static constexpr Ctx<T> ctx("rotg", "OO|", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(T, a);
            SCALAR_REQ(T, b);

            /* Though, c should have been REAL in crotg/zrotg, historically
             * _fblas assumed complex and returned garbage in c.imag. Now,
             * at least c.imag is 0. */
            real_of_t<T> c{};
            T s{};
            blas::rotg(a, b, c, s);

            T cc = static_cast<T>(c);   /* widen real c to the flavor type, imag = 0 */
            RETURN(cc, s);   /* out=(c, s) */
        }


        template <class T>
        static PyObject *rotmg(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"d1", "d2", "x1", "y1", nullptr};
            static constexpr Ctx<T> ctx("rotmg", "OOOO|", kwlist);
            PARSE_ARGS();

            SCALAR_REQ(T, d1);
            SCALAR_REQ(T, d2);
            SCALAR_REQ(T, x1);
            SCALAR_REQ(T, y1);

            py_ref param = ctx.zeros(5);
            if (!param) { return nullptr; }

            /* d1, d2, x1 are updated in place by the routine but exposed as intent(in), so
             * the updates are discarded, as in f2py. */
            blas::rotmg(d1, d2, x1, y1, param.data<T>());
            RETURN(param);   /* out=param */
        }


        template <class T>
        static PyObject *rot(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"x", "y", "c", "s", "n", "offx", "incx", "offy", "incy", "overwrite_x", "overwrite_y", nullptr};
            static constexpr Ctx<T> ctx(tchar<T>(), "rot", "OOOO|OOOOOii", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_x);
            SCALAR_FLAG(overwrite_y);
            ARRAY_INOUT(x, 1, overwrite_x != 0);
            ARRAY_INOUT(y, 1, overwrite_y != 0);

            /* c and s are real also for the complex flavors (csrot/zdrot) */
            SCALAR_REQ(real_of_t<T>, c);
            SCALAR_REQ(real_of_t<T>, s);

            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx != 0, incx);
            SCALAR_OPT(CBLAS_INT, incy, 1);  CHECK(incy != 0, incy);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);
            SCALAR_OPT(CBLAS_INT, offy, 0);  CHECK(offy >= 0 && offy < len(y), offy);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - 1 - offx) / abs(incx) + 1);
            CHECK(len(y) - offy > (n - 1) * abs(incy), n);
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            blas::rot(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy, c, s);
            RETURN(x, y);   /* out=(x, y) */
        }


        template <class T>
        static PyObject *rotm(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"x", "y", "param", "n", "offx", "incx", "offy", "incy", "overwrite_x", "overwrite_y", nullptr};
            static constexpr Ctx<T> ctx("rotm", "OOO|OOOOOii", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_x);
            SCALAR_FLAG(overwrite_y);
            ARRAY_INOUT(x, 1, overwrite_x != 0);
            ARRAY_INOUT(y, 1, overwrite_y != 0);
            ARRAY_IN(param, 1);
            /* f2py fixed the length at array creation ("0-th dimension must be fixed to 5") */
            CHECKARRAY(len(param) == 5, param);

            SCALAR_OPT(CBLAS_INT, incx, 1);  CHECK(incx != 0, incx);
            SCALAR_OPT(CBLAS_INT, incy, 1);  CHECK(incy != 0, incy);
            SCALAR_OPT(CBLAS_INT, offx, 0);  CHECK(offx >= 0 && offx < len(x), offx);
            SCALAR_OPT(CBLAS_INT, offy, 0);  CHECK(offy >= 0 && offy < len(y), offy);

            SCALAR_OPT(CBLAS_INT, n, (len(x) - offx) / abs(incx));
            CHECK(len(y) - offy > (n - 1) * abs(incy), n);
            CHECK(len(x) - offx > (n - 1) * abs(incx), n);

            blas::rotm(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy, param.data<T>());
            RETURN(x, y);   /* out=(x, y) */
        }


        PyMethodDef l1_methods[] = {
            BLAS_FAMILY(axpy),
            BLAS_FAMILY(copy),
            BLAS_FAMILY(rotg),
            BLAS_FAMILY(scal),
            BLAS_FAMILY(swap),
            /* Irregular function families are added individually */
            BLAS_ROW2(csscal, scal_real, c64, f32),
            BLAS_ROW2(zdscal, scal_real, c128, f64),
            BLAS_ROW(isamax, iamax, f32),
            BLAS_ROW(idamax, iamax, f64),
            BLAS_ROW(icamax, iamax, c64),
            BLAS_ROW(izamax, iamax, c128),
            BLAS_ROW(snrm2,  nrm2, f32),
            BLAS_ROW(dnrm2,  nrm2, f64),
            BLAS_ROW(scnrm2, nrm2, c64),
            BLAS_ROW(dznrm2, nrm2, c128),
            BLAS_ROW(sasum,  asum, f32),
            BLAS_ROW(dasum,  asum, f64),
            BLAS_ROW(scasum, asum, c64),
            BLAS_ROW(dzasum, asum, c128),
            BLAS_ROW(sdot,   dot,  f32),
            BLAS_ROW(ddot,   dot,  f64),
            BLAS_ROW(cdotu,  dotu, c64),
            BLAS_ROW(zdotu,  dotu, c128),
            BLAS_ROW(cdotc,  dotc, c64),
            BLAS_ROW(zdotc,  dotc, c128),
            BLAS_ROW(srot,   rot,  f32),
            BLAS_ROW(drot,   rot,  f64),
            BLAS_ROW(csrot,  rot,  c64),
            BLAS_ROW(zdrot,  rot,  c128),
            BLAS_ROW(srotm,  rotm,  f32),
            BLAS_ROW(drotm,  rotm,  f64),
            BLAS_ROW(srotmg, rotmg, f32),
            BLAS_ROW(drotmg, rotmg, f64),
            /* Sentinel */
            {nullptr, nullptr, 0, nullptr},
        };

    } // namespace blas::capi
}  // namespace blas
