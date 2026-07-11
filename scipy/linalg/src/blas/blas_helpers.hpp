/**
 * @file
 * @brief Python/numpy <-> C++ translation layer for the BLAS wrappers.
 *
 * Three groups live here:
 *
 * 1. Scalar converters (`from_pyobj`): direct ports of the f2py-generated `int_from_pyobj` /
 *    `double_from_pyobj` / `complex_double_from_pyobj` (see `build/scipy/linalg/_fblasmodule.c`),
 *    deliberately preserving its (often overly) permissive coercions such as floats truncate
 *    to int, takes the real part of a complex, and element 0 of a sequence.
 *
 * 2. Array acquisition (`as_in` / `as_inout` / `zeros_vec`) built on the contemporary
 *    `PyArray_FromAny` instead of f2py's fortranobject machinery, plus the owning reference
 *    `py_ref`.  Observable copy/in-place behavior matches f2py.
 *
 * 3. `Ctx`, the per-routine context. f2py's generator derived, from each .pyf declaration,
 *    the flavored routine name, each argument's required-vs-keyword kind and ordinal
 *    ("4th keyword") error messages, and inlined the resulting message strings at every call
 *    site.
 *
 * @note The extension is built from several .cpp files, so numpy's C-API function table must
 *       be shared: every .cpp file defines `PY_ARRAY_UNIQUE_SYMBOL` to the same name *before*
 *       including this header, and every file except `_blas_module.cpp` (which calls
 *       `import_array()` once, in module init) also defines `NO_IMPORT_ARRAY`.  Omitting
 *       these gives each translation unit its own never-initialized API table and a segfault
 *       on the first `PyArray_*` call.
 */
#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <complex>
#include <type_traits>   /* std::is_same, used with real_of to pick symv-vs-hemv naming */
#include <cassert>       /* assert: buffer/kwlist invariants                            */
#include <cstdio>        /* snprintf: format names/messages into fixed buffers (no I/O) */
#include <cstring>       /* strchr/strcmp/strlen: parse the kwlist and PyArg format     */

#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#endif
#include "numpy/arrayobject.h"
#include "numpy/arrayscalars.h"
#include "numpy/npy_math.h"

#include "scipy_blas_defines.h"   /* CBLAS_INT, CBLAS_INT_MAX */

namespace blas {

    /**
     * Exception policy: f2py used to raise a module-specific `_fblas.error` which was
     * a bare Exception subclass, reachable only through the private dunder attribute
     * `_fblas.__fblas_error` which was not possible to catch by name. These wrappers
     * raise builtin types: ValueError for argument-check failures, TypeError for
     * conversion fallbacks, OverflowError where a range check fires.
     * */

    /** @brief Scalar type -> numpy dtype number (NPY_FLOAT/NPY_DOUBLE/NPY_CFLOAT/NPY_CDOUBLE). */
    template <class T> int npy_type();
    template <> inline int npy_type<float>()                { return NPY_FLOAT;   }
    template <> inline int npy_type<double>()               { return NPY_DOUBLE;  }
    template <> inline int npy_type<std::complex<float>>()  { return NPY_CFLOAT;  }
    template <> inline int npy_type<std::complex<double>>() { return NPY_CDOUBLE; }

    /**
     * @brief Scalar type -> routine-name prefix, regular naming (s/d/c/z): `saxpy`, `zgemv`, ...
     *
     * A few Level-1 routines mixing complex data with real scalars or results carry irregular
     * two-letter prefixes instead, in two opposite orders; see tchar() and tchar_fn() below.
     */
    template <class T> const char *flavor();
    template <> inline const char *flavor<float>()                { return "s"; }
    template <> inline const char *flavor<double>()               { return "d"; }
    template <> inline const char *flavor<std::complex<float>>()  { return "c"; }
    template <> inline const char *flavor<std::complex<double>>() { return "z"; }

    /**
     * @brief Irregular prefix, flavor letter first (s/d/cs/zd) -- the pyf `<tchar>` shorthand.
     *
     * Used by the routines whose *data* is complex but one *scalar argument* stays real:
     * `csrot`/`zdrot` (the cosine is real), `csscal`/`zdscal` (real scale factor).
     */
    template <class T> const char *tchar();
    template <> inline const char *tchar<float>()                { return "s"; }
    template <> inline const char *tchar<double>()               { return "d"; }
    template <> inline const char *tchar<std::complex<float>>()  { return "cs"; }
    template <> inline const char *tchar<std::complex<double>>() { return "zd"; }

    /**
     * @brief Irregular prefix, result letter first (s/d/sc/dz) -- pyf `<prefix3>`/`<prefix4>`.
     *
     * Used by the value-returning routines whose *result* is real even for complex data:
     * `scnrm2`/`dznrm2`, `scasum`/`dzasum`.  Note the opposite letter order from tchar().
     */
    template <class T> const char *tchar_fn();
    template <> inline const char *tchar_fn<float>()                { return "s"; }
    template <> inline const char *tchar_fn<double>()               { return "d"; }
    template <> inline const char *tchar_fn<std::complex<float>>()  { return "sc"; }
    template <> inline const char *tchar_fn<std::complex<double>>() { return "dz"; }

    /**
     * @brief Index-function prefix (is/id/ic/iz) -- the pyf `i<prefix>` pattern, `isamax`.
     *
     * The third naming irregularity: the `i` for "index-returning" precedes the flavor letter.
     */
    template <class T> const char *iflavor();
    template <> inline const char *iflavor<float>()                { return "is"; }
    template <> inline const char *iflavor<double>()               { return "id"; }
    template <> inline const char *iflavor<std::complex<float>>()  { return "ic"; }
    template <> inline const char *iflavor<std::complex<double>>() { return "iz"; }

    /**
     * @brief Real counterpart of a flavor: f32 -> f32, c64 -> f32, c128 -> f64.
     *
     * The scalar type of the real arguments in the irregular families (the pyf `<ftypereal>`):
     * `c`/`s` in `csrot`/`zdrot`, the scale factor in `csscal`/`zdscal`.
     */
    template <class T> struct real_of                  { using type = T; };
    template <class T> struct real_of<std::complex<T>> { using type = T; };
    template <class T> using real_of_t = typename real_of<T>::type;

    /** @brief Conversion-target names exactly as f2py spells them in "can't be converted to" messages. */
    template <class V> const char *conv_name();
    template <> inline const char *conv_name<CBLAS_INT>()            { return "int"; }
    template <> inline const char *conv_name<float>()                { return "float"; }
    template <> inline const char *conv_name<double>()               { return "double"; }
    template <> inline const char *conv_name<std::complex<float>>()  { return "complex_float"; }
    template <> inline const char *conv_name<std::complex<double>>() { return "complex_double"; }

    /**
     * @brief Scalar -> Python number, for value-returning routines and `intent(out)` scalars.
     *
     * Matches f2py's packing: `Py_BuildValue("f"/"d")` -> Python float,
     * `pyobj_from_complex_*` -> Python complex.
     */
    inline PyObject *to_pyobj(float v)                       { return PyFloat_FromDouble(v); }
    inline PyObject *to_pyobj(double v)                      { return PyFloat_FromDouble(v); }
    inline PyObject *to_pyobj(const std::complex<float> &v)  { return PyComplex_FromDoubles(v.real(), v.imag()); }
    inline PyObject *to_pyobj(const std::complex<double> &v) { return PyComplex_FromDoubles(v.real(), v.imag()); }


    /**
     * @brief Reshape @p *parr, in place of the pointer, to exactly @p ndim dimensions --
     *        f2py's `check_and_fix_dimensions`, specialized to all-free dimensions.
     *
     * f2py reinterpreted rank mismatches instead of rejecting them, and code in the wild
     * (including scipy's own test suite: `dgemm(3, [3], [-4])`) relies on it:
     * - fewer dims than declared: trailing unit axes are appended, `[1,2] -> [[1],[2]]`
     *   (a 1-D vector is a column);
     * - more dims than declared: unit axes are squeezed, and any excess non-unit axes
     *   collapse into the last dimension (a Fortran-order flatten, `[[1,2],[3,4]] -> rank-1
     *   [1,3,2,4]`).
     * The data is already F-contiguous here, so the reshape is a zero-copy view.  The total
     * size is preserved by construction (0-sized axes excepted; numpy then raises).
     *
     * The adjustment is *virtual*, as it was in f2py: the returned working view feeds the
     * size checks and the BLAS call, while the pre-reshape array -- handed back through
     * @p orig -- is what the wrapper returns to Python, preserving f2py's observable
     * behavior (`dgemm(3, [3], [-4], 3, [5])` returns shape `(1,)`, not `(1, 1)`).
     *
     * @param orig  Receives the owning reference to the pre-reshape array when a reshape
     *              happened, nullptr otherwise.
     * @return      New reference to the rank-@p ndim working view (the input itself when the
     *              rank already matches), or nullptr with an exception set and everything freed.
     */
    inline PyArrayObject *fix_rank(PyArrayObject *arr, int ndim, PyObject **orig)
    {
        *orig = nullptr;
        int nd = PyArray_NDIM(arr);
        if (nd == ndim) { return arr; }
        npy_intp newdims[NPY_MAXDIMS];
        if (nd < ndim) {
            for (int i = 0; i < nd; i++)     { npy_intp d = PyArray_DIM(arr, i); newdims[i] = d ? d : 1; }
            for (int i = nd; i < ndim; i++)  { newdims[i] = 1; }
        }
        else {
            int j = 0;
            for (int i = 0; i < ndim; i++) {
                while (j < nd && PyArray_DIM(arr, j) < 2) { j++; }
                newdims[i] = (j < nd) ? PyArray_DIM(arr, j++) : 1;
            }
            for (; j < nd; j++) {
                if (PyArray_DIM(arr, j) >= 2) { newdims[ndim - 1] *= PyArray_DIM(arr, j); }
            }
        }
        PyArray_Dims shape = {newdims, ndim};
        PyObject *view = PyArray_Newshape(arr, &shape, NPY_FORTRANORDER);
        if (view == nullptr) { Py_DECREF(arr); return nullptr; }
        *orig = (PyObject *)arr;
        return (PyArrayObject *)view;
    }

    /**
     * @brief Acquire a read-only input: a Fortran-ordered, aligned array of type T,
     *        rank-adjusted to exactly @p ndim (see fix_rank()).
     *
     * NumPy does not steal @p o; it returns a NEW reference we own -- the caller's array itself
     * when it already fits, or a converted copy.  Either way a single Py_DECREF releases it, so
     * the caller never needs to know whether a copy was made.
     *
     * @param o     Any Python object numpy can turn into an array (FORCECAST: dtype narrowing OK).
     * @param ndim  Declared rank; other ranks are reinterpreted as f2py did.
     * @return      New reference, or nullptr with (usually) an exception set.
     */
    template <class T>
    inline PyArrayObject *as_in(PyObject *o, int ndim, PyObject **orig)
    {
        int flags = NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_ALIGNED | NPY_ARRAY_FORCECAST;
        PyArrayObject *arr = (PyArrayObject *)PyArray_FromAny(o, PyArray_DescrFromType(npy_type<T>()), 0, 0, flags, nullptr);
        *orig = nullptr;
        return (arr == nullptr) ? nullptr : fix_rank(arr, ndim, orig);
    }

    /**
     * @brief Acquire a writable input and return to Python (f2py `intent(in,out)`).
     *
     * With @p overwrite, NumPy reuses the caller's array when it fits (the in-place path);
     * otherwise, or whenever the array doesn't fit, a fresh writable copy is made (f2py's
     * INTENT_COPY semantics). Rank mismatches are reinterpreted per fix_rank(); on the
     * in-place path the result is then a zero-copy *view* of the caller's buffer.
     *
     * @note We call PyArray_FromAny directly, not the PyArray_FROMANY macro: the macro ORs
     *       C_CONTIGUOUS into the flags whenever ENSURECOPY is set, which would conflict
     *       with the Fortran-order request for 2-D arrays.
     *
     * @param overwrite  true: reuse the caller's buffer when possible; false: always copy.
     * @return           New reference (possibly the caller's array, INCREF'd), or nullptr.
     */
    template <class T>
    inline PyArrayObject *as_inout(PyObject *o, int ndim, bool overwrite, PyObject **orig)
    {
        int flags = NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_ALIGNED | NPY_ARRAY_WRITEABLE | NPY_ARRAY_FORCECAST;
        if (!overwrite) { flags |= NPY_ARRAY_ENSURECOPY; }
        PyArrayObject *arr = (PyArrayObject *)PyArray_FromAny(o, PyArray_DescrFromType(npy_type<T>()), 0, 0, flags, nullptr);
        *orig = nullptr;
        return (arr == nullptr) ? nullptr : fix_rank(arr, ndim, orig);
    }

    /**
     * @brief A fresh zero-filled 1-D array of length @p n, for optional output vectors that were
     *        not supplied.
     */
    template <class T>
    inline PyArrayObject *zeros_vec(npy_intp n)
    {
        return (PyArrayObject *)PyArray_ZEROS(1, &n, npy_type<T>(), 0);
    }

    /** @brief Fortran-ordered zero-filled matrix, for optional `intent(in,out)` matrix arguments
     *         not supplied (e.g. `?ger`'s `a`, which f2py explicitly initialized to 0). */
    template <class T>
    inline PyArrayObject *zeros_mat(npy_intp m, npy_intp n)
    {
        npy_intp dims[2] = {m, n};
        return (PyArrayObject *)PyArray_ZEROS(2, dims, npy_type<T>(), 1);
    }


    /* ---------------- scalar converters, ported from the f2py-generated module ---------------- */

    /**
     * @brief Tri-state converter result, so failure messages are built only on failure.
     *
     * Converters report *how* they failed and `Ctx::scalar` formats lazily:
     *
     * - `ok`:       value written.
     * - `fail`:     a specific exception is already set and final (e.g. OverflowError with its
     *               own text); do not touch it.
     * - `fail_msg`: apply f2py's message policy -- keep an already-set exception's *type* but
     *               replace its message with the argument-specific text, or raise TypeError.
     */
    enum class conv { ok, fail, fail_msg };

    /**
     * @brief Port of f2py's `int_from_pyobj`, with the C-int target widened to CBLAS_INT.
     *
     * Accepts anything PyNumber_Long accepts (so 3.7 becomes 3), complex via its `.real`
     * attribute, and sequences via element 0; f2py's permissive coercions kept on purpose.
     *
     * @param v  Target; written only on success.
     */
    inline conv from_pyobj(CBLAS_INT *v, PyObject *obj)
    {
        PyObject *tmp = nullptr;

        if (PyLong_Check(obj)) {
            long long ll = PyLong_AsLongLong(obj);
            if (ll == -1 && PyErr_Occurred()) { return conv::fail; }
            if (ll > (long long)CBLAS_INT_MAX || ll < -(long long)CBLAS_INT_MAX - 1) {
                PyErr_SetString(PyExc_OverflowError, "Python int too large to convert to C int");
                return conv::fail;
            }
            *v = static_cast<CBLAS_INT>(ll);   /* lossless: range-checked above */
            return conv::ok;
        }

        tmp = PyNumber_Long(obj);
        if (tmp) {
            conv r = from_pyobj(v, tmp);
            Py_DECREF(tmp);
            return r;
        }

        if (PyComplex_Check(obj)) {
            PyErr_Clear();
            tmp = PyObject_GetAttrString(obj, "real");
        }
        else if (PyBytes_Check(obj) || PyUnicode_Check(obj)) {
            /* pass */
        }
        else if (PySequence_Check(obj)) {
            PyErr_Clear();
            tmp = PySequence_GetItem(obj, 0);
        }
        if (tmp) {
            if (from_pyobj(v, tmp) == conv::ok) { Py_DECREF(tmp); return conv::ok; }
            Py_DECREF(tmp);
        }
        return conv::fail_msg;
    }

    /** @brief Port of f2py's `double_from_pyobj`; same coercion rules as the int port. */
    inline conv from_pyobj(double *v, PyObject *obj)
    {
        PyObject *tmp = nullptr;

        if (PyFloat_Check(obj)) {
            *v = PyFloat_AsDouble(obj);
            return (*v == -1.0 && PyErr_Occurred()) ? conv::fail : conv::ok;
        }

        tmp = PyNumber_Float(obj);
        if (tmp) {
            *v = PyFloat_AsDouble(tmp);
            Py_DECREF(tmp);
            return (*v == -1.0 && PyErr_Occurred()) ? conv::fail : conv::ok;
        }

        if (PyComplex_Check(obj)) {
            PyErr_Clear();
            tmp = PyObject_GetAttrString(obj, "real");
        }
        else if (PyBytes_Check(obj) || PyUnicode_Check(obj)) {
            /* pass */
        }
        else if (PySequence_Check(obj)) {
            PyErr_Clear();
            tmp = PySequence_GetItem(obj, 0);
        }
        if (tmp) {
            if (from_pyobj(v, tmp) == conv::ok) { Py_DECREF(tmp); return conv::ok; }
            Py_DECREF(tmp);
        }
        return conv::fail_msg;
    }

    /** @brief Port of f2py's `float_from_pyobj`: converts through double, then narrows. */
    inline conv from_pyobj(float *v, PyObject *obj)
    {
        double d = 0.0;
        conv r = from_pyobj(&d, obj);
        if (r == conv::ok) { *v = (float)d; }
        return r;
    }

    /**
     * @brief Port of f2py's `complex_double_from_pyobj`.
     *
     * Handles Python complex, the numpy complex scalar types, 0-d arrays / array scalars, then
     * falls through the real cases by hand (Python has no PyNumber_Complex), and finally
     * sequences via element 0.
     */
    inline conv from_pyobj(std::complex<double> *v, PyObject *obj)
    {
        if (PyComplex_Check(obj)) {
            Py_complex c = PyComplex_AsCComplex(obj);
            *v = {c.real, c.imag};
            return conv::ok;
        }
        if (PyArray_IsScalar(obj, ComplexFloating)) {
            if (PyArray_IsScalar(obj, CFloat)) {
                npy_cfloat val;
                PyArray_ScalarAsCtype(obj, &val);
                *v = {(double)npy_crealf(val), (double)npy_cimagf(val)};
            }
            else if (PyArray_IsScalar(obj, CLongDouble)) {
                npy_clongdouble val;
                PyArray_ScalarAsCtype(obj, &val);
                *v = {(double)npy_creall(val), (double)npy_cimagl(val)};
            }
            else {  /* CDouble */
                npy_cdouble val;
                PyArray_ScalarAsCtype(obj, &val);
                *v = {npy_creal(val), npy_cimag(val)};
            }
            return conv::ok;
        }
        if (PyArray_CheckScalar(obj)) {  /* 0-dim array or other array scalar */
            PyArrayObject *arr;
            if (PyArray_Check(obj)) {
                arr = (PyArrayObject *)PyArray_Cast((PyArrayObject *)obj, NPY_CDOUBLE);
            }
            else {
                arr = (PyArrayObject *)PyArray_FromScalar(obj, PyArray_DescrFromType(NPY_CDOUBLE));
            }
            if (arr == nullptr) { return conv::fail; }
            npy_cdouble val = *(npy_cdouble *)PyArray_DATA(arr);
            *v = {npy_creal(val), npy_cimag(val)};
            Py_DECREF(arr);
            return conv::ok;
        }
        /* Python provides no PyNumber_Complex, so fall through the real cases by hand. */
        if (PyFloat_Check(obj)) {
            *v = {PyFloat_AsDouble(obj), 0.0};
            return (v->real() == -1.0 && PyErr_Occurred()) ? conv::fail : conv::ok;
        }
        if (PyLong_Check(obj)) {
            *v = {PyLong_AsDouble(obj), 0.0};
            return (v->real() == -1.0 && PyErr_Occurred()) ? conv::fail : conv::ok;
        }
        if (PySequence_Check(obj) && !(PyBytes_Check(obj) || PyUnicode_Check(obj))) {
            PyObject *tmp = PySequence_GetItem(obj, 0);
            if (tmp) {
                if (from_pyobj(v, tmp) == conv::ok) { Py_DECREF(tmp); return conv::ok; }
                Py_DECREF(tmp);
            }
        }
        return conv::fail_msg;
    }

    /** @brief Port of f2py's `complex_float_from_pyobj`: converts through complex_double, then narrows. */
    inline conv from_pyobj(std::complex<float> *v, PyObject *obj)
    {
        std::complex<double> cd;
        conv r = from_pyobj(&cd, obj);
        if (r == conv::ok) { *v = {(float)cd.real(), (float)cd.imag()}; }
        return r;
    }


    /**
     * @brief Owning reference to an acquired array; the reason the wrappers have no `goto fail`.
     *
     * Every acquisition returns exactly one new reference (see `as_in`/`as_inout`), so ownership
     * is uniform: hold it in a py_ref and every early `return nullptr` releases it automatically,
     * replacing f2py's nested-brace cleanup pyramid.  Move-only, like the reference it manages.
     */
    class py_ref {
    public:
        py_ref() = default;
        /** @param orig  The pre-reshape array when fix_rank() adjusted the rank (see there);
         *               it is what release() hands back to Python. */
        explicit py_ref(PyArrayObject *p, PyObject *orig = nullptr) : p_(p), orig_(orig) {}
        py_ref(py_ref &&o) noexcept : p_(o.p_), orig_(o.orig_) { o.p_ = nullptr; o.orig_ = nullptr; }
        py_ref(const py_ref &) = delete;
        py_ref &operator=(const py_ref &) = delete;
        py_ref &operator=(py_ref &&o) noexcept
        {
            if (this != &o) {
                Py_XDECREF(p_); Py_XDECREF(orig_);
                p_ = o.p_; orig_ = o.orig_;
                o.p_ = nullptr; o.orig_ = nullptr;
            }
            return *this;
        }
        ~py_ref() { Py_XDECREF(p_); Py_XDECREF(orig_); }

        /** @brief False when acquisition failed (held pointer is null). */
        explicit operator bool() const { return p_ != nullptr; }
        PyArrayObject *get() const { return p_; }
        /** @brief Data pointer as the wrapper's scalar type; the array's dtype guarantees the cast. */
        template <class T> T *data() const { return static_cast<T *>(PyArray_DATA(p_)); }
        /** @brief Hand the owned reference to Python -- the wrapper's return value.  When a
         *         rank adjustment made the working array a view, the caller-shaped original
         *         is returned instead (the data is shared, writes are visible in both). */
        PyObject *release()
        {
            PyObject *r = orig_ ? orig_ : (PyObject *)p_;
            if (orig_ != nullptr) { Py_XDECREF(p_); }
            p_ = nullptr; orig_ = nullptr;
            return r;
        }

    private:
        PyArrayObject *p_ = nullptr;
        PyObject *orig_ = nullptr;
    };

    /**
     * @brief Shorthand for the wrappers' size checks: `len(x)`, `shape(a, i)`, `abs(incx)`.
     *
     * The check expressions are stringized into error messages, so they should stay readable as
     * prose; these keep them close to the .pyf's `check(...)` vocabulary without macro tricks.
     *
     * `len`/`shape` return CBLAS_INT so dimensions can flow to the BLAS calls with no cast at
     * any call site.  The static_cast here is lossless: every py_ref the wrappers hold
     * came through `Ctx::checked` / `Ctx::zeros`, which reject any array with a dimension above
     * CBLAS_INT_MAX.
     *
     * @note `blas::abs` intentionally takes and returns npy_intp, so size products like
     *       `(n - 1) * abs(incx)` are computed 64-bit.  f2py evaluated them in C int, which can
     *       overflow before the comparison. Because the wrappers live in `namespace blas::capi`,
     *       unqualified `abs` finds this overload and never the C library's `::abs(int)`, that
     *       is the name lookup stops at `blas` before reaching the global scope.
     */
    inline CBLAS_INT len(const py_ref &x)          { return static_cast<CBLAS_INT>(PyArray_DIM(x.get(), 0)); }
    inline CBLAS_INT shape(const py_ref &x, int i) { return static_cast<CBLAS_INT>(PyArray_DIM(x.get(), i)); }
    inline npy_intp  abs(npy_intp v)               { return v < 0 ? -v : v; }


    /**
     * @brief Per-routine context: everything f2py's generator derived from the .pyf declaration.
     *
     * It derives the flavored routine name (`"saxpy"`), whether an argument is a required positional
     * or an optional keyword (its position relative to `'|'`), and its f2py ordinal (`"2nd argument"`
     * / `"4th keyword"`).
     *
     * @tparam T  The wrapper's scalar flavor (float/double/complex<float>/complex<double>).
     */
    template <class T>
    class Ctx {
    public:
        /**
         * @param prefix  Explicit routine-name prefix for the irregular L1 families -- pass the
         *                matching trait, e.g. `tchar_fn<T>()` so `nrm2` becomes `scnrm2` at
         *                `T = c64`.  Regular routines use the two-argument overload below.
         * @param name    Unprefixed routine name (`"axpy"`, `"nrm2"`).
         * @param pyfmt   PyArg_ParseTupleAndKeywords format without the `:name` suffix; the
         *                position of `'|'` defines how many leading arguments are required.
         * @param kwlist  Null-terminated argument-name list, in signature order.  Must outlive
         *                the Ctx (the wrappers use a function-local static).
         */
        Ctx(const char *prefix, const char *name, const char *pyfmt, const char *const *kwlist)
            : kwlist_(kwlist)
        {
            /* growth (LAPACK names, long kwlists) must fail loudly, not truncate silently */
            int r = snprintf(rout_, sizeof rout_, "%s%s", prefix, name);
            assert(r > 0 && static_cast<size_t>(r) < sizeof rout_ && "rout_ buffer too small");
            r = snprintf(fmt_, sizeof fmt_, "%s:_fblas.%s", pyfmt, rout_);
            assert(r > 0 && static_cast<size_t>(r) < sizeof fmt_ && "fmt_ buffer too small");
            (void)r;
            const char *bar = strchr(pyfmt, '|');
            nreq_ = (int)(bar ? bar - pyfmt : (std::ptrdiff_t)strlen(pyfmt));
        }

        /** @brief Regular naming: prepends the flavor letter, `"axpy"` -> `saxpy`/.../`zaxpy`. */
        Ctx(const char *name, const char *pyfmt, const char *const *kwlist)
            : Ctx(flavor<T>(), name, pyfmt, kwlist) {}

        /** @brief Full PyArg format, e.g. `"OO|OOOOOO:_fblas.saxpy"`. */
        const char *fmt() const { return fmt_; }
        /** @brief The kwlist as PyArg_ParseTupleAndKeywords wants it (non-const for the C API). */
        char **kws() const { return const_cast<char **>(kwlist_); }

        /**
         * @brief Position of @p kw in the kwlist (linear scan; the lists are < 12 entries).
         *
         * The GETSCALAR macro caches the result in a per-site `static const int`, so the scan
         * runs once per process, not once per call.
         */
        int index(const char *kw) const
        {
            int i = 0;
            while (kwlist_[i] && strcmp(kwlist_[i], kw) != 0) { i++; }
            /* A miss means the kwlist and a GETSCALAR/CHECK* name are out of sync -- a
            * transcription bug, not a runtime condition.  Fail loudly instead of letting the
            * end sentinel silently classify the argument as an optional keyword. */
            assert(kwlist_[i] != nullptr && "argument name missing from kwlist");
            return i;
        }

        /**
         * @brief Convert one scalar argument, unconditionally.
         *
         * The provided-or-default branch lives in the GETSCALAR macro (an omitted optional never
         * reaches this function); here the object is always converted.  On conversion failure
         * raises `_fblas.saxpy() 4th keyword (incx) can't be converted to int`.
         *
         * @param v    Conversion target; written only on success.
         * @param kw   Argument name as it appears in the kwlist.
         * @param idx  The kwlist position of @p kw, as returned by `index(kw)` (the GETSCALAR
         *             macros cache it per site); determines the "4th keyword" ordinal.
         * @return true on success; false with an exception set.
         */
        template <class V>
        bool scalar(V *v, PyObject *obj, const char *kw, int idx) const
        {
            conv r = from_pyobj(v, obj);
            if (r == conv::ok) { return true; }
            if (r == conv::fail_msg) {
                /* f2py's policy: keep an already-set exception's type, replace its message */
                PyObject *err = PyErr_Occurred();
                char pk[32], msg[128];
                poskind(pk, sizeof pk, idx);
                int n = snprintf(msg, sizeof msg, "_fblas.%s() %s (%s) can't be converted to %s",
                                rout_, pk, kw, conv_name<V>());
                assert(n > 0 && static_cast<size_t>(n) < sizeof msg && "msg buffer too small");
                (void)n;
                PyErr_SetString(err ? err : PyExc_TypeError, msg);
            }
            return false;
        }

        /**
         * @brief f2py CHECKSCALAR; no-op when @p ok.
         *
         * On failure raises `(incx>0||incx<0) failed for 4th keyword incx: saxpy:incx=0`.
         *
         * @param tcheck  The check expression's source text (the CHECKSCALAR macro stringizes it).
         * @param val     The offending value, printed after `name=`; widened to long long so
         *                ILP64 CBLAS_INT values print correctly on LLP64 platforms.
         */
        bool check_scalar(bool ok, const char *tcheck, const char *kw, long long val) const
        {
            if (ok) { return true; }
            char pk[32];
            poskind(pk, sizeof pk, index(kw));
            PyErr_Format(PyExc_ValueError, "(%s) failed for %s %s: %s:%s=%lld", tcheck, pk, kw, rout_, kw, val);
            return false;
        }

        /**
         * @brief f2py CHECKARRAY; no-op when @p ok.
         *
         * On failure raises `(shape(a,0)==shape(a,1)) failed for 2nd argument a`.
         */
        bool check_array(bool ok, const char *tcheck, const char *name) const
        {
            if (ok) { return true; }
            char pk[32];
            poskind(pk, sizeof pk, index(name));
            PyErr_Format(PyExc_ValueError, "(%s) failed for %s %s", tcheck, pk, name);
            return false;
        }

        /** @brief `as_in` (f2py `intent(in)`) with owning result and f2py's failure message. */
        py_ref in(PyObject *o, int ndim, const char *name) const
        {
            PyObject *orig = nullptr;
            /* two statements: orig must be written before it is read as an argument */
            PyArrayObject *arr = as_in<T>(o, ndim, &orig);
            return checked(arr, orig, name);
        }
        /** @brief `as_inout` (f2py `intent(in,out[,copy])`) with owning result and failure message. */
        py_ref inout(PyObject *o, int ndim, bool overwrite, const char *name) const
        {
            PyObject *orig = nullptr;
            PyArrayObject *arr = as_inout<T>(o, ndim, overwrite, &orig);
            return checked(arr, orig, name);
        }
        /** @brief `zeros_vec` with owning result, for optional output vectors not supplied.
         *         Guarded by the same CBLAS_INT_MAX certificate as acquired arrays, so
         *         `len()` stays lossless on internally allocated arrays too. */
        py_ref zeros(npy_intp n) const
        {
            if (n > (npy_intp)CBLAS_INT_MAX) {
                PyErr_Format(PyExc_OverflowError,
                            "_fblas.%s(): required output size %lld exceeds the BLAS integer"
                            " limit (%lld)",
                            rout_, static_cast<long long>(n), static_cast<long long>(CBLAS_INT_MAX));
                return py_ref(nullptr);
            }
            return py_ref(zeros_vec<T>(n));
        }
        /** @brief Two-dimensional variant of zeros(), Fortran-ordered, same certificate. */
        py_ref zeros(npy_intp m, npy_intp n) const
        {
            if (m > (npy_intp)CBLAS_INT_MAX || n > (npy_intp)CBLAS_INT_MAX) {
                PyErr_Format(PyExc_OverflowError,
                            "_fblas.%s(): required output shape (%lld, %lld) exceeds the BLAS"
                            " integer limit (%lld)",
                            rout_, static_cast<long long>(m), static_cast<long long>(n),
                            static_cast<long long>(CBLAS_INT_MAX));
                return py_ref(nullptr);
            }
            return py_ref(zeros_mat<T>(m, n));
        }

    private:
        /**
         * @brief Format the f2py position string, `"2nd argument"` / `"4th keyword"`.
         *
         * f2py numbers required positionals and optional keywords as two separate sequences
         * (gemv's `beta`, 4th in the signature, is the "1st keyword").  The generated `_fblas`
         * never exceeds "7th", so the simple 1st/2nd/3rd/Nth suffix rule is exact.
         */
        void poskind(char *buf, size_t size, int i) const
        {
            bool req = i < nreq_;
            int num = req ? i + 1 : i - nreq_ + 1;
            const char *suf = num == 1 ? "st" : num == 2 ? "nd" : num == 3 ? "rd" : "th";
            snprintf(buf, size, "%d%s %s", num, suf, req ? "argument" : "keyword");
        }

        /**
         * @brief Wrap a fresh array reference: creation-failure message plus the dimension guard.
         *
         * When creation failed without setting an exception, raises f2py's fallback message
         * `_fblas._fblas.saxpy: failed to create array from the 1st argument `x``
         * (doubled prefix kept -- it is what f2py emits).
         *
         * On success, this is the single point where npy_intp dimensions are certified to fit
         * CBLAS_INT: every dimension of every acquired array is checked against CBLAS_INT_MAX,
         * and an OverflowError mentions the routine, the argument, and the offending dimension.
         * Downstream, `len()`/`shape()` return CBLAS_INT relying on this certificate, so no
         * other narrowing exists anywhere in the wrappers.
         */
        py_ref checked(PyArrayObject *a, PyObject *orig, const char *name) const
        {
            if (a == nullptr) {
                Py_XDECREF(orig);
                if (!PyErr_Occurred()) {
                    char pk[32];
                    poskind(pk, sizeof pk, index(name));
                    PyErr_Format(PyExc_TypeError, "_fblas._fblas.%s: failed to create array from the %s `%s`",
                                rout_, pk, name);
                }
                return py_ref(nullptr);
            }
            for (int i = 0; i < PyArray_NDIM(a); i++) {
                if (PyArray_DIM(a, i) > (npy_intp)CBLAS_INT_MAX) {
                    char pk[32];
                    poskind(pk, sizeof pk, index(name));
                    PyErr_Format(PyExc_OverflowError,
                                "_fblas.%s() %s (%s): dimension %d with size %lld exceeds"
                                " the BLAS integer limit (%lld)",
                                rout_, pk, name, i,
                                static_cast<long long>(PyArray_DIM(a, i)),
                                static_cast<long long>(CBLAS_INT_MAX));
                    Py_DECREF(a);
                    Py_XDECREF(orig);
                    return py_ref(nullptr);
                }
            }
            return py_ref(a, orig);
        }

        char rout_[16];
        char fmt_[64];
        const char *const *kwlist_;
        int nreq_;
    };

}  // namespace blas


/**
 * ====================================================================================
 * Following are helper macros designed to mimic the legacy f2py-generated .pyf wrappers,
 * used until SciPy 2.0. They are not intended to be used in new code and only serve
 * to keep compatibility with the legacy, occasionally buggy, behavior.
 *
 * Though C++ frowns upon macros, these are the only way to concisely mimic f2py's
 * behavior without introducing a lot of boilerplate code. And they are strictly used
 * to perform stringization of the arguments and to mimic f2py's error messages.
 *
 * Conventions inside a wrapper (see the _blas_l*.cpp files):
 *
 * - `x_obj` is the raw Python argument exactly as `PyArg_ParseTupleAndKeywords` produced it.
 * - A bare array name (`x`, `y`, `a`, `b`) is the *acquired* array, held by an owning
 *   `py_ref` that releases its reference on every exit path.
 * - Bare lower-case names (`n`, `incx`, `alpha`, ...) are, if successful, converted C scalars.
 * - `ctx` is the per-routine context (flavored routine name, keyword list, required-argument
 *   count).  The vocabulary macros below expand against it by convention.
 * - Checks are stringified C-expressions evaluated as-is. Some helper macros stringify
 *   the variables into an error message, `len()`, `shape()` and `abs()` are actual functions
 *   (see above), so those expressions compile as written.
 * ====================================================================================
 */

/**
 * @brief Assign the optional scalar `<name>`: convert `<name>_obj` if the argument was
 *        provided, otherwise assign the .pyf default expression.
 *
 * `GETSCALAR(incx, 1)` is f2py's `integer optional :: incx = 1`, and expands to the branch
 * f2py generated: `if (incx_obj == Py_None) incx = 1; else <convert incx_obj into incx>`.
 * `None` -- omitted or explicitly passed -- selects the default, exactly as in f2py.  The
 * default may be the .pyf expression verbatim, over already-assigned values
 * (`GETSCALAR(n, (len(x) - offx) / abs(incx))`), evaluated only when needed.  Conversion
 * follows f2py's permissive rules (see `from_pyobj` above): floats truncate to int, complex
 * contributes its real part, a sequence contributes element 0.
 *
 * The default is cast to the variable's own type (`decltype(name)`) -- the single narrowing
 * point for derived defaults.  This is lossless by construction: defaults are transcribed
 * .pyf expressions over certified dimensions (see `checked()`), never user input; user input
 * takes the converter branch, which range-checks.
 *
 * @param name  Bare identifier: the target C variable (declared, possibly uninitialized --
 *              every path through the macro assigns it); also selects `<name>_obj` and the
 *              Python-level argument name.
 * @param def   The .pyf default value or expression, assigned when the argument is absent.
 * @note On conversion failure raises TypeError ("... can't be converted to int"; a more
 *       specific exception from the conversion attempt wins) and executes `return nullptr`
 *       in the enclosing wrapper; `py_ref` locals clean themselves up.
 */
#define GETSCALAR(name, def) \
    do { \
        if (name##_obj == Py_None) { name = static_cast<decltype(name)>(def); } \
        else { \
            static const int name##_kwidx = ctx.index(#name); \
            if (!ctx.scalar(&name, name##_obj, #name, name##_kwidx)) { return nullptr; } \
        } \
    } while (0)

/**
 * @brief Convert the required scalar `<name>` from `<name>_obj`, unconditionally.
 *
 * For required positionals (`alpha` in gemv/trsm) there is no default: whatever the caller
 * passed -- including an explicit `None` -- is converted, and failure raises f2py's
 * "can't be converted to" message.
 */
#define GETSCALAR_REQ(name) \
    do { \
        static const int name##_kwidx = ctx.index(#name); \
        if (!ctx.scalar(&name, name##_obj, #name, name##_kwidx)) { return nullptr; } \
    } while (0)

/**
 * @brief Validate a converted scalar against a `check(...)` C-expression.
 *
 * `CHECKSCALAR(incx != 0, incx)` stringizes the expression into f2py's CHECKSCALAR message
 * format, e.g. `(incx != 0) failed for 4th keyword incx: saxpy:incx=0`; the ordinal
 * ("4th keyword") comes from `ctx`.
 *
 * @param expr  The C expression for the check to be performed. Any bare comma needs to be parenthesized.
 * @param name  Bare identifier of the checked scalar; reported by name and value.
 *
 * @note On failure raises ValueError and executes `return nullptr` in the wrapper.
 */
#define CHECKSCALAR(expr, name) do { if (!ctx.check_scalar((expr), #expr, #name, static_cast<long long>(name))) { return nullptr; } } while (0)

/**
 * @brief Validate an acquired array against a `check(...)` C-expression.
 *
 * `CHECKARRAY(shape(a, 0) == shape(a, 1), a)` stringizes the expression into f2py's
 * CHECKARRAY message format, e.g. `(shape(a, 0) == shape(a, 1)) failed for 2nd argument a`.
 *
 * @param expr  The C expression for the check to be performed. Any bare comma needs to be parenthesized.
 * @param name  Bare identifier of the checked array (a `py_ref` local).
 *
 * @note On failure raises ValueError and executes `return nullptr` in the wrapper.
 */
#define CHECKARRAY(expr, name) do { if (!ctx.check_array((expr), #expr, #name)) { return nullptr; } } while (0)

/**
 * @brief Acquire `<name>_obj` as a Fortran-ordered array of the wrapper's flavor,
 *        declaring the owning local `py_ref <name>`.
 *
 * `GETARRAY_IN` is f2py `intent(in)`: the caller's array is used as-is when it already fits
 * (dtype, exact rank, F-contiguous, aligned), otherwise a converted copy is made.
 * `GETARRAY_INOUT` is `intent(in,out)`: with @p overwrite the caller's array is reused in
 * place when it fits; without it a fresh writable copy is always made (f2py `intent(copy)`).
 * The acquired array is what the wrapper returns to Python via `release()`.
 *
 * @param name       Bare identifier: declares `py_ref name` from `<name>_obj`.
 * @param ndim       Required rank, matched exactly.
 * @param overwrite  (INOUT only) `bool`: reuse the caller's buffer when possible.
 * @note Expands to two statements (declaration + failure return): use at block level only,
 *       never as the body of an unbraced `if`/`else`.  On failure an exception is already
 *       set (or f2py's array-creation message is raised) and the wrapper returns nullptr.
 */
#define GETARRAY_IN(name, ndim) \
    py_ref name = ctx.in(name##_obj, ndim, #name); \
    if (!name) { return nullptr; }
#define GETARRAY_INOUT(name, ndim, overwrite) \
    py_ref name = ctx.inout(name##_obj, ndim, overwrite, #name); \
    if (!name) { return nullptr; }

/**
 * @brief Emit the four method-table rows of one routine -- `s<name>`, `d<name>`, `c<name>`,
 *        `z<name>` -- each bound to the wrapper template instantiated at the matching scalar
 *        type.  Expand inside `namespace blas::capi`, where the wrappers and the width
 *        aliases resolve unqualified.
 */
#define BLAS_FAMILY(name) \
    {"s" #name, (PyCFunction)(void (*)())name<f32>,  METH_VARARGS | METH_KEYWORDS, nullptr}, \
    {"d" #name, (PyCFunction)(void (*)())name<f64>,  METH_VARARGS | METH_KEYWORDS, nullptr}, \
    {"c" #name, (PyCFunction)(void (*)())name<c64>,  METH_VARARGS | METH_KEYWORDS, nullptr}, \
    {"z" #name, (PyCFunction)(void (*)())name<c128>, METH_VARARGS | METH_KEYWORDS, nullptr}

/**
 * @brief One explicitly named method-table row, for the irregular L1 families where the
 *        Python name is not flavor letter + template name (`scnrm2` binds `nrm2<c64>`).
 *        The name spelled here must agree with the prefix trait passed to the wrapper's
 *        Ctx -- that is what the parity suite's message comparisons pin down.
 */
#define BLAS_ROW(pyname, name, T) \
    {#pyname, (PyCFunction)(void (*)())name<T>, METH_VARARGS | METH_KEYWORDS, nullptr}
