/**
 * @file
 * @brief Level-1 BLAS wrappers.
 *
 * Wrapper conventions and the vocabulary macros (GETSCALAR, CHECKSCALAR, ...) are documented
 * in `blas_helpers.hpp`.  This file contributes the method-table chunk `blas::capi::l1_methods`,
 * merged into the module by `_blas_module.cpp`.
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
        static PyObject *axpy(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *x_obj, *y_obj, *n_obj = Py_None, *a_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None, *offy_obj = Py_None, *incy_obj = Py_None;
            static const char *kwlist[] = {"x", "y", "n", "a", "offx", "incx", "offy", "incy", nullptr};
            static const Ctx<T> ctx("axpy", "OO|OOOOOO", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &x_obj, &y_obj, &n_obj, &a_obj, &offx_obj, &incx_obj, &offy_obj, &incy_obj)) {
                return nullptr;
            }

            GETARRAY_IN(x, 1);
            GETARRAY_INOUT(y, 1, true);

            T a;             GETSCALAR(a, T(1));
            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;  GETSCALAR(incy, 1);  CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);
            CBLAS_INT offy;  GETSCALAR(offy, 0);  CHECKSCALAR(offy >= 0 && offy < len(y), offy);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - offx) / abs(incx));
            CHECKSCALAR(len(y) - offy > (n - 1) * abs(incy), n);
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            blas::axpy(n, a, x.data<T>() + offx, incx, y.data<T>() + offy, incy);
            return y.release();   /* out=z */
        }


        template <class T>
        static PyObject *nrm2(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *x_obj, *n_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None;
            static const char *kwlist[] = {"x", "n", "offx", "incx", nullptr};
            static const Ctx<T> ctx(tchar_fn<T>(), "nrm2", "O|OOO", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &x_obj, &n_obj, &offx_obj, &incx_obj)) {
                return nullptr;
            }

            GETARRAY_IN(x, 1);

            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx > 0, incx);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - offx) / abs(incx));
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            return PyFloat_FromDouble(blas::nrm2(n, x.data<T>() + offx, incx));   /* out=n2 */
        }


        template <class T>
        static PyObject *asum(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *x_obj, *n_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None;
            static const char *kwlist[] = {"x", "n", "offx", "incx", nullptr};
            static const Ctx<T> ctx(tchar_fn<T>(), "asum", "O|OOO", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &x_obj, &n_obj, &offx_obj, &incx_obj)) {
                return nullptr;
            }

            GETARRAY_IN(x, 1);

            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - offx) / abs(incx));
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            return PyFloat_FromDouble(blas::asum(n, x.data<T>() + offx, incx));   /* out=s */
        }


        template <class T>
        static PyObject *iamax(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *x_obj, *n_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None;
            static const char *kwlist[] = {"x", "n", "offx", "incx", nullptr};
            static const Ctx<T> ctx(iflavor<T>(), "amax", "O|OOO", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &x_obj, &n_obj, &offx_obj, &incx_obj)) {
                return nullptr;
            }

            GETARRAY_IN(x, 1);

            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - offx) / abs(incx));
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            /* Fortran's 1-based index is shifted to 0-based, as the .pyf callstatement did */
            return PyLong_FromLongLong(blas::iamax(n, x.data<T>() + offx, incx) - 1);   /* out=k */
        }


        template <class T>
        static PyObject *swap(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *x_obj, *y_obj, *n_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None, *offy_obj = Py_None, *incy_obj = Py_None;
            static const char *kwlist[] = {"x", "y", "n", "offx", "incx", "offy", "incy", nullptr};
            static const Ctx<T> ctx("swap", "OO|OOOOO", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &x_obj, &y_obj, &n_obj, &offx_obj, &incx_obj, &offy_obj, &incy_obj)) {
                return nullptr;
            }

            GETARRAY_INOUT(x, 1, true);
            GETARRAY_INOUT(y, 1, true);

            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;  GETSCALAR(incy, 1);  CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);
            CBLAS_INT offy;  GETSCALAR(offy, 0);  CHECKSCALAR(offy >= 0 && offy < len(y), offy);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - offx) / abs(incx));
            CHECKSCALAR(len(y) - offy > (n - 1) * abs(incy), n);
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            blas::swap(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy);
            return Py_BuildValue("NN", x.release(), y.release());   /* out=(x, y) */
        }


        template <class T>
        static PyObject *copy(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *x_obj, *y_obj, *n_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None, *offy_obj = Py_None, *incy_obj = Py_None;
            static const char *kwlist[] = {"x", "y", "n", "offx", "incx", "offy", "incy", nullptr};
            static const Ctx<T> ctx("copy", "OO|OOOOO", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &x_obj, &y_obj, &n_obj, &offx_obj, &incx_obj, &offy_obj, &incy_obj)) {
                return nullptr;
            }

            GETARRAY_IN(x, 1);
            GETARRAY_INOUT(y, 1, true);

            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;  GETSCALAR(incy, 1);  CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);
            CBLAS_INT offy;  GETSCALAR(offy, 0);  CHECKSCALAR(offy >= 0 && offy < len(y), offy);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - offx) / abs(incx));
            CHECKSCALAR(len(y) - offy > (n - 1) * abs(incy), n);
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            blas::copy(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy);
            return y.release();   /* out=y */
        }


        template <class T>
        static PyObject *scal(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *a_obj, *x_obj, *n_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None;
            static const char *kwlist[] = {"a", "x", "n", "offx", "incx", nullptr};
            static const Ctx<T> ctx("scal", "OO|OOO", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &a_obj, &x_obj, &n_obj, &offx_obj, &incx_obj)) {
                return nullptr;
            }

            /* f2py processes a before x: a bad scalar wins over a bad array */
            T a;  GETSCALAR_REQ(a);
            GETARRAY_INOUT(x, 1, true);

            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - offx) / abs(incx));
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            blas::scal(n, a, x.data<T>() + offx, incx);
            return x.release();   /* out=x */
        }


        /* csscal/zdscal: real scale factor on complex data; unlike the regular scal, the
         * .pyf declares x intent(in,out,copy), so these two carry an overwrite_x flag and
         * copy by default. */
        template <class T>
        static PyObject *scal_real(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *a_obj, *x_obj, *n_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None;
            int overwrite_x = 0;
            static const char *kwlist[] = {"a", "x", "n", "offx", "incx", "overwrite_x", nullptr};
            static const Ctx<T> ctx(tchar<T>(), "scal", "OO|OOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &a_obj, &x_obj, &n_obj, &offx_obj, &incx_obj, &overwrite_x)) {
                return nullptr;
            }

            real_of_t<T> a;  GETSCALAR_REQ(a);
            GETARRAY_INOUT(x, 1, overwrite_x != 0);

            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - offx) / abs(incx));
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            blas::scal(n, a, x.data<T>() + offx, incx);
            return x.release();   /* out=x */
        }


        template <class T>
        static PyObject *dot(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *x_obj, *y_obj, *n_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None, *offy_obj = Py_None, *incy_obj = Py_None;
            static const char *kwlist[] = {"x", "y", "n", "offx", "incx", "offy", "incy", nullptr};
            static const Ctx<T> ctx("dot", "OO|OOOOO", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &x_obj, &y_obj, &n_obj, &offx_obj, &incx_obj, &offy_obj, &incy_obj)) {
                return nullptr;
            }

            GETARRAY_IN(x, 1);
            GETARRAY_IN(y, 1);

            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;  GETSCALAR(incy, 1);  CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);
            CBLAS_INT offy;  GETSCALAR(offy, 0);  CHECKSCALAR(offy >= 0 && offy < len(y), offy);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - offx) / abs(incx));
            CHECKSCALAR(len(y) - offy > (n - 1) * abs(incy), n);
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            return PyFloat_FromDouble(blas::dot(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy));   /* out=xy */
        }


        template <class T>
        static PyObject *dotu(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *x_obj, *y_obj, *n_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None, *offy_obj = Py_None, *incy_obj = Py_None;
            static const char *kwlist[] = {"x", "y", "n", "offx", "incx", "offy", "incy", nullptr};
            static const Ctx<T> ctx("dotu", "OO|OOOOO", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &x_obj, &y_obj, &n_obj, &offx_obj, &incx_obj, &offy_obj, &incy_obj)) {
                return nullptr;
            }

            GETARRAY_IN(x, 1);
            GETARRAY_IN(y, 1);

            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;  GETSCALAR(incy, 1);  CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);
            CBLAS_INT offy;  GETSCALAR(offy, 0);  CHECKSCALAR(offy >= 0 && offy < len(y), offy);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - offx) / abs(incx));
            CHECKSCALAR(len(y) - offy > (n - 1) * abs(incy), n);
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            return to_pyobj(blas::dotu(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy));   /* out=xy */
        }


        template <class T>
        static PyObject *dotc(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *x_obj, *y_obj, *n_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None, *offy_obj = Py_None, *incy_obj = Py_None;
            static const char *kwlist[] = {"x", "y", "n", "offx", "incx", "offy", "incy", nullptr};
            static const Ctx<T> ctx("dotc", "OO|OOOOO", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &x_obj, &y_obj, &n_obj, &offx_obj, &incx_obj, &offy_obj, &incy_obj)) {
                return nullptr;
            }

            GETARRAY_IN(x, 1);
            GETARRAY_IN(y, 1);

            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;  GETSCALAR(incy, 1);  CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);
            CBLAS_INT offy;  GETSCALAR(offy, 0);  CHECKSCALAR(offy >= 0 && offy < len(y), offy);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - offx) / abs(incx));
            CHECKSCALAR(len(y) - offy > (n - 1) * abs(incy), n);
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            return to_pyobj(blas::dotc(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy));   /* out=xy */
        }


        template <class T>
        static PyObject *rotg(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *a_obj, *b_obj;
            static const char *kwlist[] = {"a", "b", nullptr};
            static const Ctx<T> ctx("rotg", "OO|", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(), &a_obj, &b_obj)) {
                return nullptr;
            }

            T a;  GETSCALAR_REQ(a);
            T b;  GETSCALAR_REQ(b);

            /* Though, c should have been REAL in crotg/zrotg, historically
             * _fblas assumed complex and returned garbage in c.imag. Now,
             * at least c.imag is 0. */
            real_of_t<T> c{};
            T s{};
            blas::rotg(a, b, c, s);

            return Py_BuildValue("NN", to_pyobj(static_cast<T>(c)), to_pyobj(s));   /* out=(c, s) */
        }


        template <class T>
        static PyObject *rotmg(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *d1_obj, *d2_obj, *x1_obj, *y1_obj;
            static const char *kwlist[] = {"d1", "d2", "x1", "y1", nullptr};
            static const Ctx<T> ctx("rotmg", "OOOO|", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &d1_obj, &d2_obj, &x1_obj, &y1_obj)) {
                return nullptr;
            }

            T d1;  GETSCALAR_REQ(d1);
            T d2;  GETSCALAR_REQ(d2);
            T x1;  GETSCALAR_REQ(x1);
            T y1;  GETSCALAR_REQ(y1);

            py_ref param = ctx.zeros(5);
            if (!param) { return nullptr; }

            /* d1, d2, x1 are updated in place by the routine but exposed as intent(in), so
             * the updates are discarded, as in f2py. */
            blas::rotmg(d1, d2, x1, y1, param.data<T>());
            return param.release();   /* out=param */
        }


        template <class T>
        static PyObject *rot(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *x_obj, *y_obj, *c_obj, *s_obj, *n_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None, *offy_obj = Py_None, *incy_obj = Py_None;
            int overwrite_x = 0, overwrite_y = 0;
            static const char *kwlist[] = {"x", "y", "c", "s", "n", "offx", "incx", "offy", "incy", "overwrite_x", "overwrite_y", nullptr};
            static const Ctx<T> ctx(tchar<T>(), "rot", "OOOO|OOOOOii", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &x_obj, &y_obj, &c_obj, &s_obj, &n_obj, &offx_obj, &incx_obj,
                                            &offy_obj, &incy_obj, &overwrite_x, &overwrite_y)) {
                return nullptr;
            }

            GETARRAY_INOUT(x, 1, overwrite_x != 0);
            GETARRAY_INOUT(y, 1, overwrite_y != 0);

            /* c and s are real also for the complex flavors (csrot/zdrot) */
            real_of_t<T> c;  GETSCALAR_REQ(c);
            real_of_t<T> s;  GETSCALAR_REQ(s);

            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;  GETSCALAR(incy, 1);  CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);
            CBLAS_INT offy;  GETSCALAR(offy, 0);  CHECKSCALAR(offy >= 0 && offy < len(y), offy);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - 1 - offx) / abs(incx) + 1);
            CHECKSCALAR(len(y) - offy > (n - 1) * abs(incy), n);
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            blas::rot(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy, c, s);
            return Py_BuildValue("NN", x.release(), y.release());   /* out=(x, y) */
        }


        template <class T>
        static PyObject *rotm(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *x_obj, *y_obj, *param_obj, *n_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None, *offy_obj = Py_None, *incy_obj = Py_None;
            int overwrite_x = 0, overwrite_y = 0;
            static const char *kwlist[] = {"x", "y", "param", "n", "offx", "incx", "offy", "incy", "overwrite_x", "overwrite_y", nullptr};
            static const Ctx<T> ctx("rotm", "OOO|OOOOOii", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &x_obj, &y_obj, &param_obj, &n_obj, &offx_obj, &incx_obj,
                                            &offy_obj, &incy_obj, &overwrite_x, &overwrite_y)) {
                return nullptr;
            }

            GETARRAY_INOUT(x, 1, overwrite_x != 0);
            GETARRAY_INOUT(y, 1, overwrite_y != 0);
            GETARRAY_IN(param, 1);
            /* f2py fixed the length at array creation ("0-th dimension must be fixed to 5") */
            CHECKARRAY(len(param) == 5, param);

            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;  GETSCALAR(incy, 1);  CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;  GETSCALAR(offx, 0);  CHECKSCALAR(offx >= 0 && offx < len(x), offx);
            CBLAS_INT offy;  GETSCALAR(offy, 0);  CHECKSCALAR(offy >= 0 && offy < len(y), offy);

            CBLAS_INT n;     GETSCALAR(n, (len(x) - offx) / abs(incx));
            CHECKSCALAR(len(y) - offy > (n - 1) * abs(incy), n);
            CHECKSCALAR(len(x) - offx > (n - 1) * abs(incx), n);

            blas::rotm(n, x.data<T>() + offx, incx, y.data<T>() + offy, incy, param.data<T>());
            return Py_BuildValue("NN", x.release(), y.release());   /* out=(x, y) */
        }


        PyMethodDef l1_methods[] = {
            BLAS_FAMILY(axpy),
            BLAS_FAMILY(copy),
            BLAS_FAMILY(rotg),
            BLAS_FAMILY(scal),
            BLAS_FAMILY(swap),
            /* Irregular function families are added individually */
            BLAS_ROW(csscal, scal_real, c64),
            BLAS_ROW(zdscal, scal_real, c128),
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
