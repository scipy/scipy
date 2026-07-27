/**
 * @file
 * @brief Level-3 BLAS wrappers.
 *
 * Wrapper conventions and the vocabulary macros (PARSE_ARGS, the SCALAR_ and ARRAY_ argument
 * declarations, CHECK, RETURN, ...) are documented in `blas_helpers.hpp`.  This file contributes
 * the method-table chunk `blas::capi::l3_methods`, merged into the module by `_blas_module.cpp`.
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
        static PyObject *gemm(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"alpha", "a", "b", "beta", "c", "trans_a", "trans_b", "overwrite_c", nullptr};
            static const Ctx<T> ctx("gemm", "OOO|OOOOi", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);

            SCALAR_OPT(CBLAS_INT, trans_a, 0);  CHECK(trans_a >= 0 && trans_a <= 2, trans_a);
            SCALAR_OPT(CBLAS_INT, trans_b, 0);  CHECK(trans_b >= 0 && trans_b <= 2, trans_b);
            SCALAR_REQ(T, alpha);
            SCALAR_OPT(T, beta, T(0));

            ARRAY_IN(a, 2);
            ARRAY_IN(b, 2);

            CBLAS_INT m = trans_a ? shape(a, 1) : shape(a, 0);
            CBLAS_INT k = trans_a ? shape(a, 0) : shape(a, 1);
            /* f2py phrased the inner-dimension match as a check on "hidden n" */
            CHECKARRAY(trans_b ? shape(b, 1) == k : shape(b, 0) == k, b);
            CBLAS_INT n = trans_b ? shape(b, 0) : shape(b, 1);

            ARRAY_OUT(c, 2, overwrite_c != 0, ctx.zeros(m, n));
            CHECKARRAY(shape(c, 0) == m && shape(c, 1) == n, c);

            blas::gemm(trans_a ? (trans_a == 2 ? 'C' : 'T') : 'N', trans_b ? (trans_b == 2 ? 'C' : 'T') : 'N',
                       m, n, k, alpha, a.data<T>(), shape(a, 0), b.data<T>(), shape(b, 0),
                       beta, c.data<T>(), m);
            RETURN(c);   /* out=c */
        }


        /* symm is standard four-flavor BLAS at level 3 (csymm/zsymm are not auxiliaries);
         * hemm below carries the hermitian pair. */
        template <class T>
        static PyObject *symm(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"alpha", "a", "b", "beta", "c", "side", "lower", "overwrite_c", nullptr};
            static const Ctx<T> ctx("symm", "OOO|OOOOi", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);

            SCALAR_OPT(CBLAS_INT, side, 0);   CHECK(side == 0 || side == 1, side);
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);
            SCALAR_REQ(T, alpha);
            SCALAR_OPT(T, beta, T(0));

            ARRAY_IN(a, 2);
            ARRAY_IN(b, 2);

            CBLAS_INT m = side ? shape(b, 0) : shape(a, 0);
            CHECKARRAY(side ? shape(b, 1) == shape(a, 0) : shape(a, 1) == shape(b, 0), b);
            CBLAS_INT n = side ? shape(a, 1) : shape(b, 1);

            ARRAY_OUT(c, 2, overwrite_c != 0, ctx.zeros(m, n));
            CHECKARRAY(shape(c, 0) == m && shape(c, 1) == n, c);

            blas::symm(side ? 'R' : 'L', lower ? 'L' : 'U', m, n, alpha,
                       a.data<T>(), shape(a, 0), b.data<T>(), shape(b, 0), beta, c.data<T>(), m);
            RETURN(c);   /* out=c */
        }


        template <class T>
        static PyObject *hemm(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"alpha", "a", "b", "beta", "c", "side", "lower", "overwrite_c", nullptr};
            static const Ctx<T> ctx("hemm", "OOO|OOOOi", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);

            SCALAR_OPT(CBLAS_INT, side, 0);   CHECK(side == 0 || side == 1, side);
            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);
            SCALAR_REQ(T, alpha);
            SCALAR_OPT(T, beta, T(0));

            ARRAY_IN(a, 2);
            ARRAY_IN(b, 2);

            CBLAS_INT m = side ? shape(b, 0) : shape(a, 0);
            CHECKARRAY(side ? shape(b, 1) == shape(a, 0) : shape(a, 1) == shape(b, 0), b);
            CBLAS_INT n = side ? shape(a, 1) : shape(b, 1);

            ARRAY_OUT(c, 2, overwrite_c != 0, ctx.zeros(m, n));
            CHECKARRAY(shape(c, 0) == m && shape(c, 1) == n, c);

            blas::hemm(side ? 'R' : 'L', lower ? 'L' : 'U', m, n, alpha,
                       a.data<T>(), shape(a, 0), b.data<T>(), shape(b, 0), beta, c.data<T>(), m);
            RETURN(c);   /* out=c */
        }


        template <class T>
        static PyObject *syrk(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"alpha", "a", "beta", "c", "trans", "lower", "overwrite_c", nullptr};
            static const Ctx<T> ctx("syrk", "OO|OOOOi", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);

            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);
            SCALAR_OPT(CBLAS_INT, trans, 0);  CHECK(trans >= 0 && trans <= 2, trans);
            SCALAR_REQ(T, alpha);
            SCALAR_OPT(T, beta, T(0));

            ARRAY_IN(a, 2);
            CBLAS_INT n = trans ? shape(a, 1) : shape(a, 0);
            CBLAS_INT k = trans ? shape(a, 0) : shape(a, 1);

            ARRAY_OUT(c, 2, overwrite_c != 0, ctx.zeros(n, n));
            CHECKARRAY(shape(c, 0) == n && shape(c, 1) == n, c);

            blas::syrk(lower ? 'L' : 'U', trans ? (trans == 2 ? 'C' : 'T') : 'N',
                       n, k, alpha, a.data<T>(), shape(a, 0), beta, c.data<T>(), n);
            RETURN(c);   /* out=c */
        }


        /* herk's alpha AND beta are real in the reference; the .pyf declared them complex
         * and the miscast pointer passed the leading real component.  Same conversion,
         * honest prototype: complex accepted, imaginary parts ignored. */
        template <class T>
        static PyObject *herk(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"alpha", "a", "beta", "c", "trans", "lower", "overwrite_c", nullptr};
            static const Ctx<T> ctx("herk", "OO|OOOOi", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);

            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);
            SCALAR_OPT(CBLAS_INT, trans, 0);  CHECK(trans >= 0 && trans <= 2, trans);
            SCALAR_REQ(T, alpha);
            SCALAR_OPT(T, beta, T(0));

            ARRAY_IN(a, 2);
            CBLAS_INT n = trans ? shape(a, 1) : shape(a, 0);
            CBLAS_INT k = trans ? shape(a, 0) : shape(a, 1);

            ARRAY_OUT(c, 2, overwrite_c != 0, ctx.zeros(n, n));
            CHECKARRAY(shape(c, 0) == n && shape(c, 1) == n, c);

            blas::herk(lower ? 'L' : 'U', trans ? (trans == 2 ? 'C' : 'T') : 'N',
                       n, k, alpha.real(), a.data<T>(), shape(a, 0), beta.real(), c.data<T>(), n);
            RETURN(c);   /* out=c */
        }


        template <class T>
        static PyObject *syr2k(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"alpha", "a", "b", "beta", "c", "trans", "lower", "overwrite_c", nullptr};
            static const Ctx<T> ctx("syr2k", "OOO|OOOOi", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);

            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);
            SCALAR_OPT(CBLAS_INT, trans, 0);  CHECK(trans >= 0 && trans <= 2, trans);
            SCALAR_REQ(T, alpha);
            SCALAR_OPT(T, beta, T(0));

            ARRAY_IN(a, 2);
            ARRAY_IN(b, 2);

            CBLAS_INT n = trans ? shape(a, 1) : shape(a, 0);
            /* f2py phrased the a/b conformance as a check on "hidden k" */
            CHECKARRAY(trans ? shape(a, 0) == shape(b, 0) : shape(a, 1) == shape(b, 1), b);
            CBLAS_INT k = trans ? shape(a, 0) : shape(a, 1);

            ARRAY_OUT(c, 2, overwrite_c != 0, ctx.zeros(n, n));
            CHECKARRAY(shape(c, 0) == n && shape(c, 1) == n, c);

            blas::syr2k(lower ? 'L' : 'U', trans ? (trans == 2 ? 'C' : 'T') : 'N',
                        n, k, alpha, a.data<T>(), shape(a, 0), b.data<T>(), shape(b, 0),
                        beta, c.data<T>(), n);
            RETURN(c);   /* out=c */
        }


        /* her2k's beta is real in the reference (alpha stays complex); same arrangement as
         * herk. */
        template <class T>
        static PyObject *her2k(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"alpha", "a", "b", "beta", "c", "trans", "lower", "overwrite_c", nullptr};
            static const Ctx<T> ctx("her2k", "OOO|OOOOi", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_c);

            SCALAR_OPT(CBLAS_INT, lower, 0);  CHECK(lower == 0 || lower == 1, lower);
            SCALAR_OPT(CBLAS_INT, trans, 0);  CHECK(trans >= 0 && trans <= 2, trans);
            SCALAR_REQ(T, alpha);
            SCALAR_OPT(T, beta, T(0));

            ARRAY_IN(a, 2);
            ARRAY_IN(b, 2);

            CBLAS_INT n = trans ? shape(a, 1) : shape(a, 0);
            CHECKARRAY(trans ? shape(a, 0) == shape(b, 0) : shape(a, 1) == shape(b, 1), b);
            CBLAS_INT k = trans ? shape(a, 0) : shape(a, 1);

            ARRAY_OUT(c, 2, overwrite_c != 0, ctx.zeros(n, n));
            CHECKARRAY(shape(c, 0) == n && shape(c, 1) == n, c);

            blas::her2k(lower ? 'L' : 'U', trans ? (trans == 2 ? 'C' : 'T') : 'N',
                        n, k, alpha, a.data<T>(), shape(a, 0), b.data<T>(), shape(b, 0),
                        beta.real(), c.data<T>(), n);
            RETURN(c);   /* out=c */
        }


        template <class T>
        static PyObject *trmm(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"alpha", "a", "b", "side", "lower", "trans_a", "diag", "overwrite_b", nullptr};
            static const Ctx<T> ctx("trmm", "OOO|OOOOi", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);

            SCALAR_OPT(CBLAS_INT, side, 0);     CHECK(side == 0 || side == 1, side);
            SCALAR_OPT(CBLAS_INT, lower, 0);    CHECK(lower == 0 || lower == 1, lower);
            SCALAR_OPT(CBLAS_INT, trans_a, 0);  CHECK(trans_a >= 0 && trans_a <= 2, trans_a);
            SCALAR_OPT(CBLAS_INT, diag, 0);     CHECK(diag == 0 || diag == 1, diag);
            SCALAR_REQ(T, alpha);

            ARRAY_IN(a, 2);
            ARRAY_INOUT(b, 2, overwrite_b != 0);
            CBLAS_INT m = shape(b, 0), n = shape(b, 1);
            /* f2py's "hidden k" check: a must span the multiplied dimension */
            CHECKARRAY(shape(a, 1) >= (side ? n : m) && shape(a, 1) <= shape(a, 0), a);

            blas::trmm(side ? 'R' : 'L', lower ? 'L' : 'U', trans_a ? (trans_a == 2 ? 'C' : 'T') : 'N',
                       diag ? 'U' : 'N', m, n, alpha,
                       a.data<T>(), shape(a, 0) > 1 ? shape(a, 0) : 1,
                       b.data<T>(), m > 1 ? m : 1);
            RETURN(b);   /* out=b */
        }


        template <class T>
        static PyObject *trsm(PyObject *, PyObject *args, PyObject *kwds) noexcept
        {
            static const char *kwlist[] = {"alpha", "a", "b", "side", "lower", "trans_a", "diag", "overwrite_b", nullptr};
            static const Ctx<T> ctx("trsm", "OOO|OOOOi", kwlist);
            PARSE_ARGS();

            SCALAR_FLAG(overwrite_b);

            SCALAR_OPT(CBLAS_INT, side, 0);     CHECK(side == 0 || side == 1, side);
            SCALAR_OPT(CBLAS_INT, lower, 0);    CHECK(lower == 0 || lower == 1, lower);
            SCALAR_OPT(CBLAS_INT, trans_a, 0);  CHECK(trans_a >= 0 && trans_a <= 2, trans_a);
            SCALAR_OPT(CBLAS_INT, diag, 0);     CHECK(diag == 0 || diag == 1, diag);
            SCALAR_REQ(T, alpha);

            ARRAY_INOUT(b, 2, overwrite_b != 0);
            CBLAS_INT m = shape(b, 0), n = shape(b, 1);

            ARRAY_IN(a, 2);
            CHECKARRAY(shape(a, 0) == (side ? n : m), a);
            CHECKARRAY(shape(a, 0) == shape(a, 1), a);

            blas::trsm(side ? 'R' : 'L', lower ? 'L' : 'U', trans_a ? (trans_a == 2 ? 'C' : 'T') : 'N',
                    diag ? 'U' : 'N', m, n, alpha,
                    a.data<T>(), shape(a, 0) > 1 ? shape(a, 0) : 1,
                    b.data<T>(), m > 1 ? m : 1);
            RETURN(b);   /* out=x */
        }


        PyMethodDef l3_methods[] = {
            BLAS_FAMILY(gemm),
            BLAS_FAMILY(symm),
            BLAS_FAMILY(syr2k),
            BLAS_FAMILY(syrk),
            BLAS_FAMILY(trmm),
            BLAS_FAMILY(trsm),
            /* Irregular function families are added individually */
            BLAS_ROW(chemm,  hemm,  c64),
            BLAS_ROW(zhemm,  hemm,  c128),
            BLAS_ROW(cherk,  herk,  c64),
            BLAS_ROW(zherk,  herk,  c128),
            BLAS_ROW(cher2k, her2k, c64),
            BLAS_ROW(zher2k, her2k, c128),
            /* Sentinel */
            {nullptr, nullptr, 0, nullptr},
        };

    }  // namespace blas::capi
}  // namespace blas
