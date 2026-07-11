/**
 * @file
 * @brief Level-2 BLAS wrappers.
 *
 * Wrapper conventions and the vocabulary macros (GETSCALAR, CHECKSCALAR, ...) are documented
 * in `blas_helpers.hpp`.  This file contributes the method-table chunk
 * `blas::capi::l2_methods`, merged into the module by `_blas_module.cpp`.
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
        static PyObject *gemv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *alpha_obj, *a_obj, *x_obj, *beta_obj = Py_None, *y_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None, *offy_obj = Py_None, *incy_obj = Py_None, *trans_obj = Py_None;
            int overwrite_y = 0;
            static const char *kwlist[] = {"alpha", "a", "x", "beta", "y", "offx", "incx", "offy", "incy", "trans", "overwrite_y", nullptr};
            static const Ctx<T> ctx("gemv", "OOO|OOOOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &alpha_obj, &a_obj, &x_obj, &beta_obj, &y_obj, &offx_obj, &incx_obj,
                                            &offy_obj, &incy_obj, &trans_obj, &overwrite_y)) {
                return nullptr;
            }

            CBLAS_INT trans;  GETSCALAR(trans, 0);  CHECKSCALAR(trans >= 0 && trans <= 2, trans);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;   GETSCALAR(incy, 1);   CHECKSCALAR(incy != 0, incy);
            T alpha;          GETSCALAR_REQ(alpha);
            T beta;           GETSCALAR(beta, T(0));

            GETARRAY_IN(a, 2);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT offy;   GETSCALAR(offy, 0);

            CBLAS_INT m = shape(a, 0), n = shape(a, 1);
            CBLAS_INT rows = trans ? n : m, cols = trans ? m : n;

            GETARRAY_IN(x, 1);
            CHECKARRAY(len(x) > offx + (cols - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            py_ref y = (y_obj == Py_None) ? ctx.zeros(1 + offy + (rows - 1) * abs(incy)) : ctx.inout(y_obj, 1, overwrite_y != 0, "y");
            if (!y) { return nullptr; }
            CHECKARRAY(len(y) > offy + (rows - 1) * abs(incy), y);
            CHECKARRAY(offy >= 0 && offy < len(y), y);

            blas::gemv(trans ? (trans == 2 ? 'C' : 'T') : 'N', m, n, alpha, a.data<T>(), m, x.data<T>() + offx, incx, beta, y.data<T>() + offy, incy);
            return y.release();
        }


        template <class T>
        static PyObject *gbmv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *m_obj, *n_obj, *kl_obj, *ku_obj, *alpha_obj, *a_obj, *x_obj, *incx_obj = Py_None, *offx_obj = Py_None, *beta_obj = Py_None, *y_obj = Py_None, *incy_obj = Py_None, *offy_obj = Py_None, *trans_obj = Py_None;
            int overwrite_y = 0;
            static const char *kwlist[] = {"m", "n", "kl", "ku", "alpha", "a", "x", "incx", "offx", "beta", "y", "incy", "offy", "trans", "overwrite_y", nullptr};
            static const Ctx<T> ctx("gbmv", "OOOOOOO|OOOOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &m_obj, &n_obj, &kl_obj, &ku_obj, &alpha_obj, &a_obj, &x_obj,
                                            &incx_obj, &offx_obj, &beta_obj, &y_obj, &incy_obj, &offy_obj,
                                            &trans_obj, &overwrite_y)) {
                return nullptr;
            }

            CBLAS_INT trans;  GETSCALAR(trans, 0);  CHECKSCALAR(trans >= 0 && trans <= 2, trans);
            CBLAS_INT kl;     GETSCALAR_REQ(kl);    CHECKSCALAR(kl >= 0, kl);
            CBLAS_INT ku;     GETSCALAR_REQ(ku);    CHECKSCALAR(ku >= 0, ku);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;   GETSCALAR(incy, 1);   CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT offy;   GETSCALAR(offy, 0);
            T alpha;          GETSCALAR_REQ(alpha);
            T beta;           GETSCALAR(beta, T(0));

            GETARRAY_IN(a, 2);
            CBLAS_INT m;      GETSCALAR_REQ(m);     CHECKSCALAR(m >= ku + kl + 1, m);
            CBLAS_INT n;      GETSCALAR_REQ(n);     CHECKSCALAR(n >= 0 && n == shape(a, 1), n);
            CBLAS_INT lda = shape(a, 0) > 1 ? shape(a, 0) : 1;

            /* the y/x lengths swap roles under transposition */
            py_ref y = (y_obj == Py_None) ? ctx.zeros(1 + offy + ((trans == 0 ? m : n) - 1) * abs(incy)) : ctx.inout(y_obj, 1, overwrite_y != 0, "y");
            if (!y) { return nullptr; }
            CHECKARRAY(len(y) > offy + ((trans == 0 ? m : n) - 1) * abs(incy), y);
            CHECKARRAY(offy >= 0 && offy < len(y), y);

            GETARRAY_IN(x, 1);
            CHECKARRAY(len(x) > offx + ((trans == 0 ? n : m) - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            blas::gbmv(trans ? (trans == 2 ? 'C' : 'T') : 'N', m, n, kl, ku, alpha,
                       a.data<T>(), lda, x.data<T>() + offx, incx, beta, y.data<T>() + offy, incy);
            return y.release();   /* out=yout */
        }


        /* The banded symmetric family has no complex-symmetric members in LAPACK (unlike
         * symv): sbmv is s/d only, the complex pair is the hermitian hbmv below. */
        template <class T>
        static PyObject *sbmv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *k_obj, *alpha_obj, *a_obj, *x_obj, *incx_obj = Py_None, *offx_obj = Py_None, *beta_obj = Py_None, *y_obj = Py_None, *incy_obj = Py_None, *offy_obj = Py_None, *lower_obj = Py_None;
            int overwrite_y = 0;
            static const char *kwlist[] = {"k", "alpha", "a", "x", "incx", "offx", "beta", "y", "incy", "offy", "lower", "overwrite_y", nullptr};
            static const Ctx<T> ctx("sbmv", "OOOO|OOOOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &k_obj, &alpha_obj, &a_obj, &x_obj, &incx_obj, &offx_obj,
                                            &beta_obj, &y_obj, &incy_obj, &offy_obj, &lower_obj,
                                            &overwrite_y)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;   GETSCALAR(incy, 1);   CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT offy;   GETSCALAR(offy, 0);
            T alpha;          GETSCALAR_REQ(alpha);
            T beta;           GETSCALAR(beta, T(0));

            GETARRAY_IN(a, 2);
            CBLAS_INT n = shape(a, 1);
            CBLAS_INT lda = shape(a, 0) > 1 ? shape(a, 0) : 1;
            CBLAS_INT k;      GETSCALAR_REQ(k);     CHECKSCALAR(k >= 0 && k <= lda - 1, k);

            py_ref y = (y_obj == Py_None) ? ctx.zeros(1 + offy + (n - 1) * abs(incy)) : ctx.inout(y_obj, 1, overwrite_y != 0, "y");
            if (!y) { return nullptr; }
            CHECKARRAY(len(y) > offy + (n - 1) * abs(incy), y);
            CHECKARRAY(offy >= 0 && offy < len(y), y);

            GETARRAY_IN(x, 1);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            blas::sbmv(lower ? 'L' : 'U', n, k, alpha, a.data<T>(), lda,
                       x.data<T>() + offx, incx, beta, y.data<T>() + offy, incy);
            return y.release();   /* out=yout */
        }


        template <class T>
        static PyObject *hbmv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *k_obj, *alpha_obj, *a_obj, *x_obj, *incx_obj = Py_None, *offx_obj = Py_None, *beta_obj = Py_None, *y_obj = Py_None, *incy_obj = Py_None, *offy_obj = Py_None, *lower_obj = Py_None;
            int overwrite_y = 0;
            static const char *kwlist[] = {"k", "alpha", "a", "x", "incx", "offx", "beta", "y", "incy", "offy", "lower", "overwrite_y", nullptr};
            static const Ctx<T> ctx("hbmv", "OOOO|OOOOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &k_obj, &alpha_obj, &a_obj, &x_obj, &incx_obj, &offx_obj,
                                            &beta_obj, &y_obj, &incy_obj, &offy_obj, &lower_obj,
                                            &overwrite_y)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;   GETSCALAR(incy, 1);   CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT offy;   GETSCALAR(offy, 0);
            T alpha;          GETSCALAR_REQ(alpha);
            T beta;           GETSCALAR(beta, T(0));

            GETARRAY_IN(a, 2);
            CBLAS_INT n = shape(a, 1);
            CBLAS_INT lda = shape(a, 0) > 1 ? shape(a, 0) : 1;
            CBLAS_INT k;      GETSCALAR_REQ(k);     CHECKSCALAR(k >= 0 && k <= lda - 1, k);

            py_ref y = (y_obj == Py_None) ? ctx.zeros(1 + offy + (n - 1) * abs(incy)) : ctx.inout(y_obj, 1, overwrite_y != 0, "y");
            if (!y) { return nullptr; }
            CHECKARRAY(len(y) > offy + (n - 1) * abs(incy), y);
            CHECKARRAY(offy >= 0 && offy < len(y), y);

            GETARRAY_IN(x, 1);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            blas::hbmv(lower ? 'L' : 'U', n, k, alpha, a.data<T>(), lda,
                       x.data<T>() + offx, incx, beta, y.data<T>() + offy, incy);
            return y.release();   /* out=yout */
        }


        template <class T>
        static PyObject *symv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *alpha_obj, *a_obj, *x_obj, *beta_obj = Py_None, *y_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None, *offy_obj = Py_None, *incy_obj = Py_None, *lower_obj = Py_None;
            int overwrite_y = 0;
            static const char *kwlist[] = {"alpha", "a", "x", "beta", "y", "offx", "incx", "offy", "incy", "lower", "overwrite_y", nullptr};
            static const Ctx<T> ctx("symv", "OOO|OOOOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &alpha_obj, &a_obj, &x_obj, &beta_obj, &y_obj, &offx_obj, &incx_obj,
                                            &offy_obj, &incy_obj, &lower_obj, &overwrite_y)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;   GETSCALAR(incy, 1);   CHECKSCALAR(incy != 0, incy);
            T alpha;          GETSCALAR_REQ(alpha);
            T beta;           GETSCALAR(beta, T(0));

            GETARRAY_IN(a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT offy;   GETSCALAR(offy, 0);
            CBLAS_INT n = shape(a, 0);

            GETARRAY_IN(x, 1);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            py_ref y = (y_obj == Py_None) ? ctx.zeros(1 + offy + (n - 1) * abs(incy)) : ctx.inout(y_obj, 1, overwrite_y != 0, "y");
            if (!y) { return nullptr; }
            CHECKARRAY(len(y) > offy + (n - 1) * abs(incy), y);
            CHECKARRAY(offy >= 0 && offy < len(y), y);

            blas::symv(lower ? 'L' : 'U', n, alpha, a.data<T>(), n, x.data<T>() + offx, incx, beta, y.data<T>() + offy, incy);
            return y.release();
        }


        template <class T>
        static PyObject *hemv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *alpha_obj, *a_obj, *x_obj, *beta_obj = Py_None, *y_obj = Py_None, *offx_obj = Py_None, *incx_obj = Py_None, *offy_obj = Py_None, *incy_obj = Py_None, *lower_obj = Py_None;
            int overwrite_y = 0;
            static const char *kwlist[] = {"alpha", "a", "x", "beta", "y", "offx", "incx", "offy", "incy", "lower", "overwrite_y", nullptr};
            static const Ctx<T> ctx("hemv", "OOO|OOOOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &alpha_obj, &a_obj, &x_obj, &beta_obj, &y_obj, &offx_obj, &incx_obj,
                                            &offy_obj, &incy_obj, &lower_obj, &overwrite_y)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;   GETSCALAR(incy, 1);   CHECKSCALAR(incy != 0, incy);
            T alpha;          GETSCALAR_REQ(alpha);
            T beta;           GETSCALAR(beta, T(0));

            GETARRAY_IN(a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT offy;   GETSCALAR(offy, 0);
            CBLAS_INT n = shape(a, 0);

            GETARRAY_IN(x, 1);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            py_ref y = (y_obj == Py_None) ? ctx.zeros(1 + offy + (n - 1) * abs(incy)) : ctx.inout(y_obj, 1, overwrite_y != 0, "y");
            if (!y) { return nullptr; }
            CHECKARRAY(len(y) > offy + (n - 1) * abs(incy), y);
            CHECKARRAY(offy >= 0 && offy < len(y), y);

            blas::hemv(lower ? 'L' : 'U', n, alpha, a.data<T>(), n, x.data<T>() + offx, incx, beta, y.data<T>() + offy, incy);
            return y.release();
        }


        /* Packed families: spmv/spr are six-flavored in netlib (cspmv/zspmv, cspr/zspr are
         * LAPACK auxiliary complex-symmetric routines, already exposed by the .pyf via
         * prefix6); the hermitian chp* members are separate wrappers.  hpr's alpha is
         * real; spr2/hpr2 has no complex-symmetric members. */
        template <class T>
        static PyObject *spmv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *n_obj, *alpha_obj, *ap_obj, *x_obj, *incx_obj = Py_None, *offx_obj = Py_None, *beta_obj = Py_None, *y_obj = Py_None, *incy_obj = Py_None, *offy_obj = Py_None, *lower_obj = Py_None;
            int overwrite_y = 0;
            static const char *kwlist[] = {"n", "alpha", "ap", "x", "incx", "offx", "beta", "y", "incy", "offy", "lower", "overwrite_y", nullptr};
            static const Ctx<T> ctx("spmv", "OOOO|OOOOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &n_obj, &alpha_obj, &ap_obj, &x_obj, &incx_obj, &offx_obj,
                                            &beta_obj, &y_obj, &incy_obj, &offy_obj, &lower_obj,
                                            &overwrite_y)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT n;      GETSCALAR_REQ(n);     CHECKSCALAR(n >= 0, n);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;   GETSCALAR(incy, 1);   CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT offy;   GETSCALAR(offy, 0);
            T alpha;          GETSCALAR_REQ(alpha);
            T beta;           GETSCALAR(beta, T(0));

            GETARRAY_IN(ap, 1);
            CHECKARRAY(len(ap) >= (npy_intp)n * (n + 1) / 2, ap);

            py_ref y = (y_obj == Py_None) ? ctx.zeros(1 + offy + (n - 1) * abs(incy)) : ctx.inout(y_obj, 1, overwrite_y != 0, "y");
            if (!y) { return nullptr; }
            CHECKARRAY(len(y) > offy + (n - 1) * abs(incy), y);
            CHECKARRAY(offy >= 0 && offy < len(y), y);

            GETARRAY_IN(x, 1);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            blas::spmv(lower ? 'L' : 'U', n, alpha, ap.data<T>(), x.data<T>() + offx, incx, beta, y.data<T>() + offy, incy);
            return y.release();   /* out=yout */
        }


        template <class T>
        static PyObject *hpmv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *n_obj, *alpha_obj, *ap_obj, *x_obj, *incx_obj = Py_None, *offx_obj = Py_None, *beta_obj = Py_None, *y_obj = Py_None, *incy_obj = Py_None, *offy_obj = Py_None, *lower_obj = Py_None;
            int overwrite_y = 0;
            static const char *kwlist[] = {"n", "alpha", "ap", "x", "incx", "offx", "beta", "y", "incy", "offy", "lower", "overwrite_y", nullptr};
            static const Ctx<T> ctx("hpmv", "OOOO|OOOOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &n_obj, &alpha_obj, &ap_obj, &x_obj, &incx_obj, &offx_obj,
                                            &beta_obj, &y_obj, &incy_obj, &offy_obj, &lower_obj,
                                            &overwrite_y)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT n;      GETSCALAR_REQ(n);     CHECKSCALAR(n >= 0, n);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;   GETSCALAR(incy, 1);   CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT offy;   GETSCALAR(offy, 0);
            T alpha;          GETSCALAR_REQ(alpha);
            T beta;           GETSCALAR(beta, T(0));

            GETARRAY_IN(ap, 1);
            CHECKARRAY(len(ap) >= (npy_intp)n * (n + 1) / 2, ap);

            py_ref y = (y_obj == Py_None) ? ctx.zeros(1 + offy + (n - 1) * abs(incy)) : ctx.inout(y_obj, 1, overwrite_y != 0, "y");
            if (!y) { return nullptr; }
            CHECKARRAY(len(y) > offy + (n - 1) * abs(incy), y);
            CHECKARRAY(offy >= 0 && offy < len(y), y);

            GETARRAY_IN(x, 1);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            blas::hpmv(lower ? 'L' : 'U', n, alpha, ap.data<T>(), x.data<T>() + offx, incx, beta, y.data<T>() + offy, incy);
            return y.release();   /* out=yout */
        }


        template <class T>
        static PyObject *spr(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *n_obj, *alpha_obj, *x_obj, *ap_obj, *incx_obj = Py_None, *offx_obj = Py_None, *lower_obj = Py_None;
            int overwrite_ap = 0;
            static const char *kwlist[] = {"n", "alpha", "x", "ap", "incx", "offx", "lower", "overwrite_ap", nullptr};
            static const Ctx<T> ctx("spr", "OOOO|OOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &n_obj, &alpha_obj, &x_obj, &ap_obj, &incx_obj, &offx_obj,
                                            &lower_obj, &overwrite_ap)) {
                return nullptr;
            }

            CBLAS_INT n;      GETSCALAR_REQ(n);     CHECKSCALAR(n >= 0, n);
            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            T alpha;          GETSCALAR_REQ(alpha);

            GETARRAY_IN(x, 1);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            GETARRAY_INOUT(ap, 1, overwrite_ap != 0);
            CHECKARRAY(len(ap) >= (npy_intp)n * (n + 1) / 2, ap);

            blas::spr(lower ? 'L' : 'U', n, alpha, x.data<T>() + offx, incx, ap.data<T>());
            return ap.release();   /* out=apu */
        }


        /* hpr's alpha is real (the diagonal of a hermitian update stays real) */
        template <class T>
        static PyObject *hpr(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *n_obj, *alpha_obj, *x_obj, *ap_obj, *incx_obj = Py_None, *offx_obj = Py_None, *lower_obj = Py_None;
            int overwrite_ap = 0;
            static const char *kwlist[] = {"n", "alpha", "x", "ap", "incx", "offx", "lower", "overwrite_ap", nullptr};
            static const Ctx<T> ctx("hpr", "OOOO|OOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &n_obj, &alpha_obj, &x_obj, &ap_obj, &incx_obj, &offx_obj,
                                            &lower_obj, &overwrite_ap)) {
                return nullptr;
            }

            CBLAS_INT n;      GETSCALAR_REQ(n);     CHECKSCALAR(n >= 0, n);
            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            real_of_t<T> alpha;  GETSCALAR_REQ(alpha);

            GETARRAY_IN(x, 1);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            GETARRAY_INOUT(ap, 1, overwrite_ap != 0);
            CHECKARRAY(len(ap) >= (npy_intp)n * (n + 1) / 2, ap);

            blas::hpr(lower ? 'L' : 'U', n, alpha, x.data<T>() + offx, incx, ap.data<T>());
            return ap.release();   /* out=apu */
        }


        template <class T>
        static PyObject *spr2(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *n_obj, *alpha_obj, *x_obj, *y_obj, *ap_obj, *incx_obj = Py_None, *offx_obj = Py_None, *incy_obj = Py_None, *offy_obj = Py_None, *lower_obj = Py_None;
            int overwrite_ap = 0;
            static const char *kwlist[] = {"n", "alpha", "x", "y", "ap", "incx", "offx", "incy", "offy", "lower", "overwrite_ap", nullptr};
            static const Ctx<T> ctx("spr2", "OOOOO|OOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &n_obj, &alpha_obj, &x_obj, &y_obj, &ap_obj, &incx_obj,
                                            &offx_obj, &incy_obj, &offy_obj, &lower_obj, &overwrite_ap)) {
                return nullptr;
            }

            CBLAS_INT n;      GETSCALAR_REQ(n);     CHECKSCALAR(n >= 0, n);
            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;   GETSCALAR(incy, 1);   CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT offy;   GETSCALAR(offy, 0);
            T alpha;          GETSCALAR_REQ(alpha);

            GETARRAY_IN(x, 1);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            GETARRAY_IN(y, 1);
            CHECKARRAY(len(y) > offy + (n - 1) * abs(incy), y);
            CHECKARRAY(offy >= 0 && offy < len(y), y);

            GETARRAY_INOUT(ap, 1, overwrite_ap != 0);
            CHECKARRAY(len(ap) >= (npy_intp)n * (n + 1) / 2, ap);

            blas::spr2(lower ? 'L' : 'U', n, alpha, x.data<T>() + offx, incx, y.data<T>() + offy, incy, ap.data<T>());
            return ap.release();   /* out=apu */
        }


        template <class T>
        static PyObject *hpr2(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *n_obj, *alpha_obj, *x_obj, *y_obj, *ap_obj, *incx_obj = Py_None, *offx_obj = Py_None, *incy_obj = Py_None, *offy_obj = Py_None, *lower_obj = Py_None;
            int overwrite_ap = 0;
            static const char *kwlist[] = {"n", "alpha", "x", "y", "ap", "incx", "offx", "incy", "offy", "lower", "overwrite_ap", nullptr};
            static const Ctx<T> ctx("hpr2", "OOOOO|OOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &n_obj, &alpha_obj, &x_obj, &y_obj, &ap_obj, &incx_obj,
                                            &offx_obj, &incy_obj, &offy_obj, &lower_obj, &overwrite_ap)) {
                return nullptr;
            }

            CBLAS_INT n;      GETSCALAR_REQ(n);     CHECKSCALAR(n >= 0, n);
            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT incy;   GETSCALAR(incy, 1);   CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT offy;   GETSCALAR(offy, 0);
            T alpha;          GETSCALAR_REQ(alpha);

            GETARRAY_IN(x, 1);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            GETARRAY_IN(y, 1);
            CHECKARRAY(len(y) > offy + (n - 1) * abs(incy), y);
            CHECKARRAY(offy >= 0 && offy < len(y), y);

            GETARRAY_INOUT(ap, 1, overwrite_ap != 0);
            CHECKARRAY(len(ap) >= (npy_intp)n * (n + 1) / 2, ap);

            blas::hpr2(lower ? 'L' : 'U', n, alpha, x.data<T>() + offx, incx, y.data<T>() + offy, incy, ap.data<T>());
            return ap.release();   /* out=apu */
        }


        template <class T>
        static PyObject *syr(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *alpha_obj, *x_obj, *lower_obj = Py_None, *incx_obj = Py_None, *offx_obj = Py_None, *n_obj = Py_None, *a_obj = Py_None;
            int overwrite_a = 0;
            static const char *kwlist[] = {"alpha", "x", "lower", "incx", "offx", "n", "a", "overwrite_a", nullptr};
            static const Ctx<T> ctx("syr", "OO|OOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &alpha_obj, &x_obj, &lower_obj, &incx_obj, &offx_obj,
                                            &n_obj, &a_obj, &overwrite_a)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            T alpha;          GETSCALAR_REQ(alpha);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);

            GETARRAY_IN(x, 1);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            CBLAS_INT n;      GETSCALAR(n, (len(x) - 1 - offx) / abs(incx) + 1);
            CHECKSCALAR(n <= (len(x) - 1 - offx) / abs(incx) + 1, n);
            CHECKSCALAR(n >= 0, n);

            py_ref a = (a_obj == Py_None) ? ctx.zeros(n, n) : ctx.inout(a_obj, 2, overwrite_a != 0, "a");
            if (!a) { return nullptr; }
            CHECKARRAY(shape(a, 0) == n && shape(a, 1) == n, a);

            blas::syr(lower ? 'L' : 'U', n, alpha, x.data<T>() + offx, incx, a.data<T>(), n);
            return a.release();   /* out=a */
        }


        /* her's alpha is real in the reference; the .pyf converted a full complex and
         * passed only its (leading) real component through the miscast function pointer.
         * Same conversion here -- complex alpha accepted, imaginary part ignored -- with
         * the honest real-pointer prototype. */
        template <class T>
        static PyObject *her(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *alpha_obj, *x_obj, *lower_obj = Py_None, *incx_obj = Py_None, *offx_obj = Py_None, *n_obj = Py_None, *a_obj = Py_None;
            int overwrite_a = 0;
            static const char *kwlist[] = {"alpha", "x", "lower", "incx", "offx", "n", "a", "overwrite_a", nullptr};
            static const Ctx<T> ctx("her", "OO|OOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &alpha_obj, &x_obj, &lower_obj, &incx_obj, &offx_obj,
                                            &n_obj, &a_obj, &overwrite_a)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            T alpha;          GETSCALAR_REQ(alpha);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);

            GETARRAY_IN(x, 1);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            CBLAS_INT n;      GETSCALAR(n, (len(x) - 1 - offx) / abs(incx) + 1);
            CHECKSCALAR(n <= (len(x) - 1 - offx) / abs(incx) + 1, n);
            CHECKSCALAR(n >= 0, n);

            py_ref a = (a_obj == Py_None) ? ctx.zeros(n, n) : ctx.inout(a_obj, 2, overwrite_a != 0, "a");
            if (!a) { return nullptr; }
            CHECKARRAY(shape(a, 0) == n && shape(a, 1) == n, a);

            blas::her(lower ? 'L' : 'U', n, alpha.real(), x.data<T>() + offx, incx, a.data<T>(), n);
            return a.release();   /* out=a */
        }


        template <class T>
        static PyObject *syr2(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *alpha_obj, *x_obj, *y_obj, *lower_obj = Py_None, *incx_obj = Py_None, *offx_obj = Py_None, *incy_obj = Py_None, *offy_obj = Py_None, *n_obj = Py_None, *a_obj = Py_None;
            int overwrite_a = 0;
            static const char *kwlist[] = {"alpha", "x", "y", "lower", "incx", "offx", "incy", "offy", "n", "a", "overwrite_a", nullptr};
            static const Ctx<T> ctx("syr2", "OOO|OOOOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &alpha_obj, &x_obj, &y_obj, &lower_obj, &incx_obj, &offx_obj,
                                            &incy_obj, &offy_obj, &n_obj, &a_obj, &overwrite_a)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            T alpha;          GETSCALAR_REQ(alpha);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT incy;   GETSCALAR(incy, 1);   CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offy;   GETSCALAR(offy, 0);

            GETARRAY_IN(x, 1);
            CHECKARRAY(offx >= 0 && offx < len(x), x);
            GETARRAY_IN(y, 1);
            CHECKARRAY(offy >= 0 && offy < len(y), y);

            npy_intp nx = (len(x) - 1 - offx) / abs(incx) + 1;
            npy_intp ny = (len(y) - 1 - offy) / abs(incy) + 1;
            CBLAS_INT n;      GETSCALAR(n, nx <= ny ? nx : ny);
            CHECKSCALAR(n <= (len(y) - 1 - offy) / abs(incy) + 1, n);
            CHECKSCALAR(n <= (len(x) - 1 - offx) / abs(incx) + 1, n);
            CHECKSCALAR(n >= 0, n);

            py_ref a = (a_obj == Py_None) ? ctx.zeros(n, n) : ctx.inout(a_obj, 2, overwrite_a != 0, "a");
            if (!a) { return nullptr; }
            CHECKARRAY(shape(a, 0) == n && shape(a, 1) == n, a);

            blas::syr2(lower ? 'L' : 'U', n, alpha, x.data<T>() + offx, incx, y.data<T>() + offy, incy, a.data<T>(), n);
            return a.release();   /* out=a */
        }


        template <class T>
        static PyObject *her2(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *alpha_obj, *x_obj, *y_obj, *lower_obj = Py_None, *incx_obj = Py_None, *offx_obj = Py_None, *incy_obj = Py_None, *offy_obj = Py_None, *n_obj = Py_None, *a_obj = Py_None;
            int overwrite_a = 0;
            static const char *kwlist[] = {"alpha", "x", "y", "lower", "incx", "offx", "incy", "offy", "n", "a", "overwrite_a", nullptr};
            static const Ctx<T> ctx("her2", "OOO|OOOOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &alpha_obj, &x_obj, &y_obj, &lower_obj, &incx_obj, &offx_obj,
                                            &incy_obj, &offy_obj, &n_obj, &a_obj, &overwrite_a)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            T alpha;          GETSCALAR_REQ(alpha);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;   GETSCALAR(offx, 0);
            CBLAS_INT incy;   GETSCALAR(incy, 1);   CHECKSCALAR(incy != 0, incy);
            CBLAS_INT offy;   GETSCALAR(offy, 0);

            GETARRAY_IN(x, 1);
            CHECKARRAY(offx >= 0 && offx < len(x), x);
            GETARRAY_IN(y, 1);
            CHECKARRAY(offy >= 0 && offy < len(y), y);

            npy_intp nx = (len(x) - 1 - offx) / abs(incx) + 1;
            npy_intp ny = (len(y) - 1 - offy) / abs(incy) + 1;
            CBLAS_INT n;      GETSCALAR(n, nx <= ny ? nx : ny);
            CHECKSCALAR(n <= (len(y) - 1 - offy) / abs(incy) + 1, n);
            CHECKSCALAR(n <= (len(x) - 1 - offx) / abs(incx) + 1, n);
            CHECKSCALAR(n >= 0, n);

            py_ref a = (a_obj == Py_None) ? ctx.zeros(n, n) : ctx.inout(a_obj, 2, overwrite_a != 0, "a");
            if (!a) { return nullptr; }
            CHECKARRAY(shape(a, 0) == n && shape(a, 1) == n, a);

            blas::her2(lower ? 'L' : 'U', n, alpha, x.data<T>() + offx, incx, y.data<T>() + offy, incy, a.data<T>(), n);
            return a.release();   /* out=a */
        }


        template <class T>
        static PyObject *ger(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *alpha_obj, *x_obj, *y_obj, *incx_obj = Py_None, *incy_obj = Py_None, *a_obj = Py_None;
            int overwrite_x = 0, overwrite_y = 0, overwrite_a = 0;
            static const char *kwlist[] = {"alpha", "x", "y", "incx", "incy", "a", "overwrite_x", "overwrite_y", "overwrite_a", nullptr};
            static const Ctx<T> ctx("ger", "OOO|OOOiii", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &alpha_obj, &x_obj, &y_obj, &incx_obj, &incy_obj, &a_obj,
                                            &overwrite_x, &overwrite_y, &overwrite_a)) {
                return nullptr;
            }

            T alpha;  GETSCALAR_REQ(alpha);
            GETARRAY_INOUT(x, 1, overwrite_x != 0);
            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx == 1 || incx == -1, incx);
            GETARRAY_INOUT(y, 1, overwrite_y != 0);
            CBLAS_INT incy;  GETSCALAR(incy, 1);  CHECKSCALAR(incy == 1 || incy == -1, incy);

            CBLAS_INT m = len(x), n = len(y);
            py_ref a = (a_obj == Py_None) ? ctx.zeros(m, n) : ctx.inout(a_obj, 2, overwrite_a != 0, "a");
            if (!a) { return nullptr; }
            /* f2py fixed a's shape at creation ("0-th dimension must be fixed to ..."); same ValueError, check-style message */
            CHECKARRAY(shape(a, 0) == m && shape(a, 1) == n, a);

            blas::ger(m, n, alpha, x.data<T>(), incx, y.data<T>(), incy, a.data<T>(), m);
            return a.release();
        }


        template <class T>
        static PyObject *geru(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *alpha_obj, *x_obj, *y_obj, *incx_obj = Py_None, *incy_obj = Py_None, *a_obj = Py_None;
            int overwrite_x = 0, overwrite_y = 0, overwrite_a = 0;
            static const char *kwlist[] = {"alpha", "x", "y", "incx", "incy", "a", "overwrite_x", "overwrite_y", "overwrite_a", nullptr};
            static const Ctx<T> ctx("geru", "OOO|OOOiii", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &alpha_obj, &x_obj, &y_obj, &incx_obj, &incy_obj, &a_obj,
                                            &overwrite_x, &overwrite_y, &overwrite_a)) {
                return nullptr;
            }

            T alpha;  GETSCALAR_REQ(alpha);
            GETARRAY_INOUT(x, 1, overwrite_x != 0);
            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx == 1 || incx == -1, incx);
            GETARRAY_INOUT(y, 1, overwrite_y != 0);
            CBLAS_INT incy;  GETSCALAR(incy, 1);  CHECKSCALAR(incy == 1 || incy == -1, incy);

            CBLAS_INT m = len(x), n = len(y);
            py_ref a = (a_obj == Py_None) ? ctx.zeros(m, n) : ctx.inout(a_obj, 2, overwrite_a != 0, "a");
            if (!a) { return nullptr; }
            CHECKARRAY(shape(a, 0) == m && shape(a, 1) == n, a);

            blas::geru(m, n, alpha, x.data<T>(), incx, y.data<T>(), incy, a.data<T>(), m);
            return a.release();
        }


        template <class T>
        static PyObject *gerc(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *alpha_obj, *x_obj, *y_obj, *incx_obj = Py_None, *incy_obj = Py_None, *a_obj = Py_None;
            int overwrite_x = 0, overwrite_y = 0, overwrite_a = 0;
            static const char *kwlist[] = {"alpha", "x", "y", "incx", "incy", "a", "overwrite_x", "overwrite_y", "overwrite_a", nullptr};
            static const Ctx<T> ctx("gerc", "OOO|OOOiii", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &alpha_obj, &x_obj, &y_obj, &incx_obj, &incy_obj, &a_obj,
                                            &overwrite_x, &overwrite_y, &overwrite_a)) {
                return nullptr;
            }

            T alpha;  GETSCALAR_REQ(alpha);
            GETARRAY_INOUT(x, 1, overwrite_x != 0);
            CBLAS_INT incx;  GETSCALAR(incx, 1);  CHECKSCALAR(incx == 1 || incx == -1, incx);
            GETARRAY_INOUT(y, 1, overwrite_y != 0);
            CBLAS_INT incy;  GETSCALAR(incy, 1);  CHECKSCALAR(incy == 1 || incy == -1, incy);

            CBLAS_INT m = len(x), n = len(y);
            py_ref a = (a_obj == Py_None) ? ctx.zeros(m, n) : ctx.inout(a_obj, 2, overwrite_a != 0, "a");
            if (!a) { return nullptr; }
            CHECKARRAY(shape(a, 0) == m && shape(a, 1) == n, a);

            blas::gerc(m, n, alpha, x.data<T>(), incx, y.data<T>(), incy, a.data<T>(), m);
            return a.release();
        }


        template <class T>
        static PyObject *trmv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *a_obj, *x_obj, *offx_obj = Py_None, *incx_obj = Py_None, *lower_obj = Py_None, *trans_obj = Py_None, *diag_obj = Py_None;
            int overwrite_x = 0;
            static const char *kwlist[] = {"a", "x", "offx", "incx", "lower", "trans", "diag", "overwrite_x", nullptr};
            static const Ctx<T> ctx("trmv", "OO|OOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &a_obj, &x_obj, &offx_obj, &incx_obj, &lower_obj, &trans_obj,
                                            &diag_obj, &overwrite_x)) {
                return nullptr;
            }

            CBLAS_INT trans;  GETSCALAR(trans, 0);  CHECKSCALAR(trans >= 0 && trans <= 2, trans);
            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT diag;   GETSCALAR(diag, 0);   CHECKSCALAR(diag == 0 || diag == 1, diag);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);

            GETARRAY_INOUT(x, 1, overwrite_x != 0);
            GETARRAY_IN(a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);

            CBLAS_INT offx;   GETSCALAR(offx, 0);   CHECKSCALAR(offx >= 0 && offx < len(x), offx);
            CBLAS_INT n = shape(a, 0);
            /* f2py reported this as a check on "hidden n"; an ordinary check on x now */
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);

            blas::trmv(lower ? 'L' : 'U', trans ? (trans == 2 ? 'C' : 'T') : 'N', diag ? 'U' : 'N',
                       n, a.data<T>(), n, x.data<T>() + offx, incx);
            return x.release();   /* out=x */
        }


        template <class T>
        static PyObject *trsv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *a_obj, *x_obj, *incx_obj = Py_None, *offx_obj = Py_None, *lower_obj = Py_None, *trans_obj = Py_None, *diag_obj = Py_None;
            int overwrite_x = 0;
            static const char *kwlist[] = {"a", "x", "incx", "offx", "lower", "trans", "diag", "overwrite_x", nullptr};
            static const Ctx<T> ctx("trsv", "OO|OOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &a_obj, &x_obj, &incx_obj, &offx_obj, &lower_obj, &trans_obj,
                                            &diag_obj, &overwrite_x)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT trans;  GETSCALAR(trans, 0);  CHECKSCALAR(trans >= 0 && trans <= 2, trans);
            CBLAS_INT diag;   GETSCALAR(diag, 0);   CHECKSCALAR(diag == 0 || diag == 1, diag);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;   GETSCALAR(offx, 0);

            GETARRAY_IN(a, 2);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);
            CBLAS_INT n = shape(a, 0);
            CBLAS_INT lda = shape(a, 0) > 1 ? shape(a, 0) : 1;

            GETARRAY_INOUT(x, 1, overwrite_x != 0);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            blas::trsv(lower ? 'L' : 'U', trans ? (trans == 2 ? 'C' : 'T') : 'N', diag ? 'U' : 'N',
                       n, a.data<T>(), lda, x.data<T>() + offx, incx);
            return x.release();   /* out=xout */
        }


        template <class T>
        static PyObject *tbmv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *k_obj, *a_obj, *x_obj, *incx_obj = Py_None, *offx_obj = Py_None, *lower_obj = Py_None, *trans_obj = Py_None, *diag_obj = Py_None;
            int overwrite_x = 0;
            static const char *kwlist[] = {"k", "a", "x", "incx", "offx", "lower", "trans", "diag", "overwrite_x", nullptr};
            static const Ctx<T> ctx("tbmv", "OOO|OOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &k_obj, &a_obj, &x_obj, &incx_obj, &offx_obj, &lower_obj,
                                            &trans_obj, &diag_obj, &overwrite_x)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT trans;  GETSCALAR(trans, 0);  CHECKSCALAR(trans >= 0 && trans <= 2, trans);
            CBLAS_INT diag;   GETSCALAR(diag, 0);   CHECKSCALAR(diag == 0 || diag == 1, diag);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;   GETSCALAR(offx, 0);

            GETARRAY_IN(a, 2);
            CBLAS_INT n = shape(a, 1);
            CBLAS_INT lda = shape(a, 0) > 1 ? shape(a, 0) : 1;
            CBLAS_INT k;      GETSCALAR_REQ(k);     CHECKSCALAR(k >= 0 && k <= lda - 1, k);

            GETARRAY_INOUT(x, 1, overwrite_x != 0);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            blas::tbmv(lower ? 'L' : 'U', trans ? (trans == 2 ? 'C' : 'T') : 'N', diag ? 'U' : 'N',
                       n, k, a.data<T>(), lda, x.data<T>() + offx, incx);
            return x.release();   /* out=xout */
        }


        template <class T>
        static PyObject *tbsv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *k_obj, *a_obj, *x_obj, *incx_obj = Py_None, *offx_obj = Py_None, *lower_obj = Py_None, *trans_obj = Py_None, *diag_obj = Py_None;
            int overwrite_x = 0;
            static const char *kwlist[] = {"k", "a", "x", "incx", "offx", "lower", "trans", "diag", "overwrite_x", nullptr};
            static const Ctx<T> ctx("tbsv", "OOO|OOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &k_obj, &a_obj, &x_obj, &incx_obj, &offx_obj, &lower_obj,
                                            &trans_obj, &diag_obj, &overwrite_x)) {
                return nullptr;
            }

            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT trans;  GETSCALAR(trans, 0);  CHECKSCALAR(trans >= 0 && trans <= 2, trans);
            CBLAS_INT diag;   GETSCALAR(diag, 0);   CHECKSCALAR(diag == 0 || diag == 1, diag);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;   GETSCALAR(offx, 0);

            GETARRAY_IN(a, 2);
            CBLAS_INT n = shape(a, 1);
            CBLAS_INT lda = shape(a, 0) > 1 ? shape(a, 0) : 1;
            CBLAS_INT k;      GETSCALAR_REQ(k);     CHECKSCALAR(k >= 0 && k <= lda - 1, k);

            GETARRAY_INOUT(x, 1, overwrite_x != 0);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            blas::tbsv(lower ? 'L' : 'U', trans ? (trans == 2 ? 'C' : 'T') : 'N', diag ? 'U' : 'N',
                       n, k, a.data<T>(), lda, x.data<T>() + offx, incx);
            return x.release();   /* out=xout */
        }


        template <class T>
        static PyObject *tpmv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *n_obj, *ap_obj, *x_obj, *incx_obj = Py_None, *offx_obj = Py_None, *lower_obj = Py_None, *trans_obj = Py_None, *diag_obj = Py_None;
            int overwrite_x = 0;
            static const char *kwlist[] = {"n", "ap", "x", "incx", "offx", "lower", "trans", "diag", "overwrite_x", nullptr};
            static const Ctx<T> ctx("tpmv", "OOO|OOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &n_obj, &ap_obj, &x_obj, &incx_obj, &offx_obj, &lower_obj,
                                            &trans_obj, &diag_obj, &overwrite_x)) {
                return nullptr;
            }

            CBLAS_INT n;      GETSCALAR_REQ(n);     CHECKSCALAR(n >= 0, n);
            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT trans;  GETSCALAR(trans, 0);  CHECKSCALAR(trans >= 0 && trans <= 2, trans);
            CBLAS_INT diag;   GETSCALAR(diag, 0);   CHECKSCALAR(diag == 0 || diag == 1, diag);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;   GETSCALAR(offx, 0);

            GETARRAY_IN(ap, 1);
            CHECKARRAY(len(ap) >= (npy_intp)n * (n + 1) / 2, ap);

            GETARRAY_INOUT(x, 1, overwrite_x != 0);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            blas::tpmv(lower ? 'L' : 'U', trans ? (trans == 2 ? 'C' : 'T') : 'N', diag ? 'U' : 'N',
                       n, ap.data<T>(), x.data<T>() + offx, incx);
            return x.release();   /* out=xout */
        }


        template <class T>
        static PyObject *tpsv(PyObject *, PyObject *args, PyObject *kwds)
        {
            PyObject *n_obj, *ap_obj, *x_obj, *incx_obj = Py_None, *offx_obj = Py_None, *lower_obj = Py_None, *trans_obj = Py_None, *diag_obj = Py_None;
            int overwrite_x = 0;
            static const char *kwlist[] = {"n", "ap", "x", "incx", "offx", "lower", "trans", "diag", "overwrite_x", nullptr};
            static const Ctx<T> ctx("tpsv", "OOO|OOOOOi", kwlist);

            if (!PyArg_ParseTupleAndKeywords(args, kwds, ctx.fmt(), ctx.kws(),
                                            &n_obj, &ap_obj, &x_obj, &incx_obj, &offx_obj, &lower_obj,
                                            &trans_obj, &diag_obj, &overwrite_x)) {
                return nullptr;
            }

            CBLAS_INT n;      GETSCALAR_REQ(n);     CHECKSCALAR(n >= 0, n);
            CBLAS_INT lower;  GETSCALAR(lower, 0);  CHECKSCALAR(lower == 0 || lower == 1, lower);
            CBLAS_INT trans;  GETSCALAR(trans, 0);  CHECKSCALAR(trans >= 0 && trans <= 2, trans);
            CBLAS_INT diag;   GETSCALAR(diag, 0);   CHECKSCALAR(diag == 0 || diag == 1, diag);
            CBLAS_INT incx;   GETSCALAR(incx, 1);   CHECKSCALAR(incx != 0, incx);
            CBLAS_INT offx;   GETSCALAR(offx, 0);

            GETARRAY_IN(ap, 1);
            CHECKARRAY(len(ap) >= (npy_intp)n * (n + 1) / 2, ap);

            GETARRAY_INOUT(x, 1, overwrite_x != 0);
            CHECKARRAY(len(x) > offx + (n - 1) * abs(incx), x);
            CHECKARRAY(offx >= 0 && offx < len(x), x);

            blas::tpsv(lower ? 'L' : 'U', trans ? (trans == 2 ? 'C' : 'T') : 'N', diag ? 'U' : 'N',
                       n, ap.data<T>(), x.data<T>() + offx, incx);
            return x.release();   /* out=xout */
        }


        PyMethodDef l2_methods[] = {
            BLAS_FAMILY(gbmv),
            BLAS_FAMILY(gemv),
            BLAS_FAMILY(tbmv),
            BLAS_FAMILY(tbsv),
            BLAS_FAMILY(tpmv),
            BLAS_FAMILY(tpsv),
            BLAS_FAMILY(trmv),
            BLAS_FAMILY(trsv),
            BLAS_FAMILY(spmv),
            BLAS_FAMILY(spr),
            BLAS_FAMILY(symv),
            BLAS_FAMILY(syr),
            /* Irregular function families are added individually */
            BLAS_ROW(chemv, hemv, c64),
            BLAS_ROW(zhemv, hemv, c128),
            BLAS_ROW(ssbmv, sbmv, f32),
            BLAS_ROW(dsbmv, sbmv, f64),
            BLAS_ROW(chbmv, hbmv, c64),
            BLAS_ROW(zhbmv, hbmv, c128),
            BLAS_ROW(chpmv, hpmv, c64),
            BLAS_ROW(zhpmv, hpmv, c128),
            BLAS_ROW(chpr,  hpr,  c64),
            BLAS_ROW(zhpr,  hpr,  c128),
            BLAS_ROW(sspr2, spr2, f32),
            BLAS_ROW(dspr2, spr2, f64),
            BLAS_ROW(chpr2, hpr2, c64),
            BLAS_ROW(zhpr2, hpr2, c128),
            BLAS_ROW(cher,  her,  c64),
            BLAS_ROW(zher,  her,  c128),
            BLAS_ROW(ssyr2, syr2, f32),
            BLAS_ROW(dsyr2, syr2, f64),
            BLAS_ROW(cher2, her2, c64),
            BLAS_ROW(zher2, her2, c128),
            BLAS_ROW(sger,  ger,  f32),
            BLAS_ROW(dger,  ger,  f64),
            BLAS_ROW(cgeru, geru, c64),
            BLAS_ROW(zgeru, geru, c128),
            BLAS_ROW(cgerc, gerc, c64),
            BLAS_ROW(zgerc, gerc, c128),
            /* Sentinel */
            {nullptr, nullptr, 0, nullptr},
        };

    } // namespace blas::capi
}  // namespace blas
