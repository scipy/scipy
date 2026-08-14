/**
 * @file
 * @brief Python/numpy <-> C++ translation layer shared by the BLAS and LAPACK wrappers.
 *
 * This header holds the machinery that used to be provided by f2py tool for wrapping BLAS and
 * LAPACK routines. Backwards compatibility is preserved as closely as possible to the point of
 * carrying over some design choices that are not ideal.
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
 * Nothing here names a Fortran library.  A wrapper module supplies its own identity through a
 * `module_id` trait (see below) and binds it once, so the same machinery serves `_fblas` and
 * `_flapack` while each keeps raising messages under its own name.
 *
 * @note The extension is built from several .cpp files, so numpy's C-API function table must
 *       be shared: every .cpp file defines `PY_ARRAY_UNIQUE_SYMBOL` to the same name *before*
 *       including this header, and every file except the module's `*_module.cpp` (which calls
 *       `import_array()` once, in module init) also defines `NO_IMPORT_ARRAY`.  Omitting
 *       these gives each translation unit its own never-initialized API table and a segfault
 *       on the first `PyArray_*` call.
 */
#pragma once

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <cassert>       /* assert: buffer/kwlist invariants                             */
#include <cstdio>        /* snprintf: format names/messages into fixed buffers (no I/O)  */
#include <cstring>       /* strchr/strcmp/strlen: parse the kwlist and PyArg format      */
#include <climits>       /* INT_MAX/INT_MIN: range check for the PyArg 'i'-style flags   */

#ifndef NPY_NO_DEPRECATED_API
#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#endif
#include "numpy/arrayobject.h"
#include "numpy/arrayscalars.h"
#include "numpy/npy_math.h"

#include "scipy_blas_defines.h"   /* CBLAS_INT, CBLAS_INT_MAX */
#include "wrapper_types.hpp"      /* f32..c128, real_of, is_complex, flavor */

namespace wrapper {

    /**
     * Exception policy: these wrappers raise builtin exception types.
     *
     * ValueError for argument-check failures, TypeError for conversion fallbacks, and
     * OverflowError for range-check failures. (f2py instead raised a custom `_fblas.error`
     * that could not be caught by name.)
     */

    /**
     * @brief What a wrapper module tells the machinery about itself.
     *
     * Both strings exist only to be printed.  `pyname` is the module every error message
     * attributes itself to: f2py used the `python module` name from the .pyf, so a module
     * replacing an f2py one must spell *that* name here (`_flapack`) rather than its own
     * while the two are built side by side.  `libname` names the Fortran library whose
     * `INTEGER` width bounds the array dimensions we accept.
     *
     * A module declares one of these and binds it once, which is the only place the library's
     * name is written down:
     *
     *     namespace blas {
     *         struct module {
     *             static constexpr const char *pyname  = "_fblas";
     *             static constexpr const char *libname = "BLAS";
     *         };
     *         template <class T> using Ctx = wrapper::Ctx<module, T>;
     *     }
     *
     * Passing it as a template parameter rather than storing it keeps `Ctx` a literal type
     * with no per-instance storage for either string, so the wrappers' `static constexpr Ctx`
     * still needs no run-time initialization.
     */

    /** @brief Scalar type -> numpy dtype number (NPY_FLOAT/NPY_DOUBLE/NPY_CFLOAT/NPY_CDOUBLE). */
    template <class T> int npy_type();
    template <> inline int npy_type<f32>()  { return NPY_FLOAT;   }
    template <> inline int npy_type<f64>()  { return NPY_DOUBLE;  }
    template <> inline int npy_type<c64>()  { return NPY_CFLOAT;  }
    template <> inline int npy_type<c128>() { return NPY_CDOUBLE; }
    /** LAPACK's integer-typed arrays (pivots, `bwork`) follow the ABI's integer width. */
    template <> inline int npy_type<CBLAS_INT>()            { return sizeof(CBLAS_INT) == sizeof(npy_int64) ? NPY_INT64 : NPY_INT32; }

    /** @brief Conversion-target names exactly as f2py spells them in "can't be converted to" messages. */
    template <class V> const char *conv_name();
    template <> inline const char *conv_name<CBLAS_INT>() { return "int"; }
    template <> inline const char *conv_name<f32>()       { return "float"; }
    template <> inline const char *conv_name<f64>()       { return "double"; }
    template <> inline const char *conv_name<c64>()       { return "complex_float"; }
    template <> inline const char *conv_name<c128>()      { return "complex_double"; }
    template <> inline const char *conv_name<char>()      { return "character"; }

    /**
     * @brief Scalar -> Python number, for value-returning routines and `intent(out)` scalars.
     *
     * Matches f2py's packing: `Py_BuildValue("f"/"d")` -> Python float,
     * `pyobj_from_complex_*` -> Python complex.
     */
    inline PyObject *to_pyobj(f32 v)         { return PyFloat_FromDouble(v); }
    inline PyObject *to_pyobj(f64 v)         { return PyFloat_FromDouble(v); }
    inline PyObject *to_pyobj(const c64 &v)  { return PyComplex_FromDoubles(v.real(), v.imag()); }
    inline PyObject *to_pyobj(const c128 &v) { return PyComplex_FromDoubles(v.real(), v.imag()); }


    /**
     * @brief Reshape @p *parr, in place of the pointer, to exactly @p ndim dimensions --
     *        f2py's `check_and_fix_dimensions`, specialized to all-free dimensions.
     *
     *  While trying to be as forgiving as possible, f2py reinterpreted rank mismatches
     *  instead of rejecting them, and code in the wild (including scipy's own test suite:
     *  `dgemm(3, [3], [-4])`) relies on it:
     *
     * - fewer dims than declared: trailing unit axes are appended, `[1,2] -> [[1],[2]]`
     *   (a 1-D vector is a column);
     * - more dims than declared: unit axes are squeezed, and any excess non-unit axes
     *   collapse into the last dimension (a Fortran-order flatten, `[[1,2],[3,4]] -> rank-1
     *   [1,3,2,4]`).
     *
     * The data is already F-contiguous here, so the reshape is a zero-copy view.  The total
     * size is preserved by construction (0-sized axes excepted; numpy then raises).
     *
     * The adjustment is *virtual*, as it was in f2py: the returned working view feeds the
     * size checks and the Fortran call, while the pre-reshape array -- handed back through
     * @p orig -- is what the wrapper returns to Python, preserving f2py's observable
     * behavior (`dgemm(3, [3], [-4], 3, [5])` returns shape `(1,)`, not `(1, 1)`).
     *
     * @param orig  Receives the owning reference to the pre-reshape array when a reshape
     *              happened, nullptr otherwise.
     * @return      New reference to the rank-@p ndim working view (the input itself when the
     *              rank already matches), or nullptr with an exception set and everything freed.
     */
    inline PyArrayObject *fix_rank(PyArrayObject *arr, int ndim, PyObject **orig) noexcept
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
     * @brief Coerce @p o to a @p descr -typed array with @p flags
     *
     * When @p o is already an ndarray, calls `PyArray_FromArray` directly, skipping
     * `PyArray_FromAny`'s dtype/shape discovery for non-array inputs (lists, scalars,
     * sequences, and other f2py's permissive cases) which take the `PyArray_FromAny`
     * path. Mimics f2py's `array_from_pyobj`.
     *
     * @param o     Reference to any Python object including PyArrayObject
     * @param descr Reference to the desired dtype descriptor
     * @param flags The resulting array requirements
     *
     * @return      New reference to the generated array or nullptr.
     *
     * @note Both steal the @p descr reference, so passing the fresh
     *       `PyArray_DescrFromType` result to either needs no extra INCREF/DECREF.
     */
    inline PyArrayObject *array_from_obj(PyObject *o, PyArray_Descr *descr, int flags) noexcept
    {
        if (PyArray_Check(o)) {
            return (PyArrayObject *)PyArray_FromArray((PyArrayObject *)o, descr, flags);
        }
        return (PyArrayObject *)PyArray_FromAny(o, descr, 0, 0, flags, nullptr);
    }

    /**
     * @brief Acquire a read-only input: a Fortran-ordered, aligned array of type T,
     *        rank-adjusted to exactly @p ndim (see fix_rank()).
     *
     * NumPy does not steal @p o; it returns a NEW reference we own -- the caller's array itself
     * when it already fits, or a converted copy.  Either way a single Py_DECREF releases it, so
     * the caller never needs to know whether a copy was made.
     *
     * @param o     Any Python object that NumPy can turn into an array (FORCECAST: dtype narrowing OK).
     * @param ndim  Declared rank; other ranks are reinterpreted as f2py did.
     * @return      New reference, or nullptr with an exception set.
     */
    template <class T>
    inline PyArrayObject *as_in(PyObject *o, int ndim, PyObject **orig) noexcept
    {
        int flags = NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_ALIGNED | NPY_ARRAY_FORCECAST;
        PyArrayObject *arr = array_from_obj(o, PyArray_DescrFromType(npy_type<T>()), flags);
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
    inline PyArrayObject *as_inout(PyObject *o, int ndim, bool overwrite, PyObject **orig) noexcept
    {
        int flags = NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_ALIGNED | NPY_ARRAY_WRITEABLE | NPY_ARRAY_FORCECAST;
        if (!overwrite) { flags |= NPY_ARRAY_ENSURECOPY; }
        PyArrayObject *arr = array_from_obj(o, PyArray_DescrFromType(npy_type<T>()), flags);
        *orig = nullptr;
        return (arr == nullptr) ? nullptr : fix_rank(arr, ndim, orig);
    }

    /**
     * @brief A fresh zero-filled 1-D array of length @p n, for optional output vectors that were
     *        not supplied.
     */
    template <class T>
    inline PyArrayObject *zeros_vec(npy_intp n) noexcept
    {
        return (PyArrayObject *)PyArray_ZEROS(1, &n, npy_type<T>(), 0);
    }

    /** @brief Fortran-ordered zero-filled matrix, for optional `intent(in,out)` matrix arguments
     *         not supplied (e.g. `?ger`'s `a`, which f2py explicitly initialized to 0). */
    template <class T>
    inline PyArrayObject *zeros_mat(npy_intp m, npy_intp n) noexcept
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
    inline conv from_pyobj(CBLAS_INT *v, PyObject *obj) noexcept
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
    inline conv from_pyobj(f64 *v, PyObject *obj) noexcept
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
    inline conv from_pyobj(f32 *v, PyObject *obj) noexcept
    {
        f64 d = 0.0;
        conv r = from_pyobj(&d, obj);
        if (r == conv::ok) { *v = (f32)d; }
        return r;
    }

    /**
     * @brief Port of f2py's `complex_double_from_pyobj`.
     *
     * Handles Python complex, the numpy complex scalar types, 0-d arrays / array scalars, then
     * falls through the real cases by hand (Python has no PyNumber_Complex), and finally
     * sequences via element 0.
     */
    inline conv from_pyobj(c128 *v, PyObject *obj) noexcept
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
                *v = {(f64)npy_crealf(val), (f64)npy_cimagf(val)};
            }
            else if (PyArray_IsScalar(obj, CLongDouble)) {
                npy_clongdouble val;
                PyArray_ScalarAsCtype(obj, &val);
                *v = {(f64)npy_creall(val), (f64)npy_cimagl(val)};
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
    inline conv from_pyobj(c64 *v, PyObject *obj) noexcept
    {
        c128 cd;
        conv r = from_pyobj(&cd, obj);
        if (r == conv::ok) { *v = {(f32)cd.real(), (f32)cd.imag()}; }
        return r;
    }

    /**
     * @brief Port of f2py's `character_from_pyobj`: the first byte of a string.
     *
     * LAPACK's option arguments are single letters ('N', 'T', 'U'), and callers spell them as
     * one-character strings. A byte string gives its first byte, a str its first byte once
     * encoded as ASCII, a string-typed array its first byte, and any other sequence its element
     * 0. Unlike the numeric converters this one takes no numbers: an int is not a letter.
     *
     * f2py also appended the offending object's type and size to the failure message; that is
     * dropped here, so a failure reads the same as every other converter's.
     *
     * @param v  Target; written only on success.
     */
    inline conv from_pyobj(char *v, PyObject *obj) noexcept
    {
        if (PyBytes_Check(obj)) {
            /* Even empty bytes carries its terminating null, so byte 0 is always readable. */
            *v = PyBytes_AS_STRING(obj)[0];
            return conv::ok;
        }
        if (PyUnicode_Check(obj)) {
            PyObject *tmp = PyUnicode_AsASCIIString(obj);
            if (tmp) {
                *v = PyBytes_AS_STRING(tmp)[0];
                Py_DECREF(tmp);
                return conv::ok;
            }
            return conv::fail_msg;   /* not ASCII, so not a LAPACK option letter */
        }
        if (PyArray_Check(obj)) {
            PyArrayObject *arr = reinterpret_cast<PyArrayObject *>(obj);
            if (PyArray_NBYTES(arr) < 1) { return conv::fail_msg; }
            int kind = PyDataType_TYPENUM(PyArray_DESCR(arr));
            if (kind == NPY_STRING || kind == NPY_UINT8) {
                *v = PyArray_BYTES(arr)[0];
                return conv::ok;
            }
            if (kind == NPY_UNICODE) {
                /* numpy holds UCS4 in native byte order, so decode instead of reading byte 0. */
                PyObject *tmp = PyUnicode_FromKindAndData(PyUnicode_4BYTE_KIND, PyArray_BYTES(arr), 1);
                if (tmp) {
                    conv r = from_pyobj(v, tmp);
                    Py_DECREF(tmp);
                    return r;
                }
            }
            return conv::fail_msg;
        }
        if (PySequence_Check(obj)) {
            PyObject *tmp = PySequence_GetItem(obj, 0);
            if (tmp) {
                conv r = from_pyobj(v, tmp);
                Py_DECREF(tmp);
                if (r == conv::ok) { return r; }
            }
            else { PyErr_Clear(); }
        }
        return conv::fail_msg;
    }

    /**
     * @brief Convert an integer flag with PyArg `i`-format (`__index__`) semantics.
     *
     * Uses PyNumber_Index: accepts an int, bool, or numpy integer; rejects a float, str, or None
     * with TypeError `'<type>' object cannot be interpreted as an integer`. Stricter than the
     * permissive from_pyobj(); kept separate so the synthetic `overwrite_*` flags reproduce
     * PyArg's message rather than f2py's "can't be converted to".
     *
     * @param out  Receives the converted int on success.
     * @param obj  Python object to convert.
     * @return     true on success (@p out written); false with an exception set.
     */
    inline bool from_pyobj_index(int *out, PyObject *obj) noexcept
    {
        PyObject *idx = PyNumber_Index(obj);
        if (idx == nullptr) { return false; }
        long v = PyLong_AsLong(idx);
        Py_DECREF(idx);
        if (v == -1 && PyErr_Occurred()) { return false; }
        if (v > INT_MAX) { PyErr_SetString(PyExc_OverflowError, "signed integer is greater than maximum"); return false; }
        if (v < INT_MIN) { PyErr_SetString(PyExc_OverflowError, "signed integer is less than minimum"); return false; }
        *out = static_cast<int>(v);
        return true;
    }


    /**
     * @brief Move-only RAII handle owning an acquired array; frees on scope exit.
     *
     * The acquisitions (as_in / as_inout) return owned references; storing them in a py_ref means
     * every early `return nullptr` releases them through the destructor, so error paths need no
     * manual cleanup. A handle owns one array normally, or two after a fix_rank() reshape (the
     * working array plus the caller-shaped original -- see private p_/orig_ below).
     */
    class py_ref {
    public:
        py_ref() = default;   // default constructor: empty handle, owns nothing

        /**
         * @param p     Acquired working array to own (null if acquisition failed).
         * @param orig  Set only when fix_rank() reshaped: the caller-shaped original that
         *              release() hands back. Defaults to null (no reshape).
         */
        explicit py_ref(PyArrayObject *p, PyObject *orig = nullptr) : p_(p), orig_(orig) {}

        // Move-only ownership: copies deleted, moves steal the references, dtor releases.
        py_ref(py_ref &&o) noexcept : p_(o.p_), orig_(o.orig_) { o.p_ = nullptr; o.orig_ = nullptr; } // move constructor
        py_ref(const py_ref &) = delete;                                                              // copy constructor (deleted)
        py_ref &operator=(const py_ref &) = delete;                                                   // copy assignment (deleted)
        py_ref &operator=(py_ref &&o) noexcept                                                        // move assignment
        {
            if (this != &o) {
                Py_XDECREF(p_); Py_XDECREF(orig_);
                p_ = o.p_; orig_ = o.orig_;
                o.p_ = nullptr; o.orig_ = nullptr;
            }
            return *this;
        }
        ~py_ref() { Py_XDECREF(p_); Py_XDECREF(orig_); }                                              // destructor: releases the references

        /** @brief True while an array is held; false after a failed acquisition or a release(). */
        explicit operator bool() const noexcept { return p_ != nullptr; }
        PyArrayObject *get() const noexcept { return p_; }
        /** @brief The working array's data as `T*`; the acquire-time dtype makes the cast sound. */
        template <class T> T *data() const noexcept { return static_cast<T *>(PyArray_DATA(p_)); }
        /**
         * @brief Transfer ownership and clean up.
         *
         * Returns the working array, or -- after a fix_rank() reshape -- the caller-shaped original
         * instead (it shares the working array's data, so the Fortran writes are visible in it). The
         * handle is left empty, so the destructor then does nothing.
         */
        PyObject *release()
        {
            PyObject *r = orig_ ? orig_ : (PyObject *)p_;
            if (orig_ != nullptr) { Py_XDECREF(p_); }
            p_ = nullptr; orig_ = nullptr;
            return r;
        }

    private:
        PyArrayObject *p_ = nullptr;   /**< Owned working array the Fortran call operates on. */
        PyObject *orig_ = nullptr;     /**< Non-null only after a fix_rank() reshape: the caller-shaped
                                            original, sharing p_'s data; handed back by release(). */
    };


    /**
     * @brief Shorthand for the wrappers' size checks: `len(x)`, `shape(a, i)`, `abs(incx)`.
     *
     * The check expressions are stringized into error messages, so they should stay readable as
     * prose; these keep them close to the .pyf's `check(...)` vocabulary without macro tricks.
     *
     * `len`/`shape` return CBLAS_INT so dimensions can flow to the Fortran calls with no cast at
     * any call site.  The static_cast here is lossless: every py_ref the wrappers hold
     * came through `Ctx::checked` / `Ctx::zeros`, which reject any array with a dimension above
     * CBLAS_INT_MAX.
     *
     * @note `abs` intentionally takes and returns npy_intp, so size products like
     *       `(n - 1) * abs(incx)` are computed 64-bit.  f2py evaluated them in C int, which can
     *       overflow before the comparison.  The wrappers live in `blas::capi` / `lapack::capi`
     *       and their enclosing namespace re-exports this overload with a using-declaration, so
     *       unqualified `abs` stops there and never reaches the C library's `::abs(int)`.
     */
    inline CBLAS_INT len(const py_ref &x)          { return static_cast<CBLAS_INT>(PyArray_DIM(x.get(), 0)); }
    inline CBLAS_INT shape(const py_ref &x, int i) { return static_cast<CBLAS_INT>(PyArray_DIM(x.get(), i)); }
    inline npy_intp  abs(npy_intp v)               { return v < 0 ? -v : v; }


    /**
     * @brief Uniform return packing for the wrappers (the RETURN macro below).
     *
     * Each result is either an acquired array (a `py_ref`, handed to Python via `release()`)
     * or a computed scalar (through `to_pyobj` / a Python int).  `make_result` maps each item
     * accordingly and, for more than one, packs them into a tuple in argument order --
     * replacing f2py's per-routine `Py_BuildValue("N"/"NN"/"NNNNNi", ...)` spellings with one
     * `RETURN(...)`.  Values are computed into locals first (no expressions tucked inside
     * RETURN); this only does the packing.
     */
    inline PyObject *result_item(py_ref &r)   { return r.release(); }
    inline PyObject *result_item(f32 v)       { return to_pyobj(v); }
    inline PyObject *result_item(f64 v)       { return to_pyobj(v); }
    inline PyObject *result_item(c64 v)       { return to_pyobj(v); }
    inline PyObject *result_item(c128 v)      { return to_pyobj(v); }
    inline PyObject *result_item(long long v) { return PyLong_FromLongLong(v); }
    /** A returned LAPACK option letter, as f2py's `Py_BuildValue("c")` gave it: one-byte bytes. */
    inline PyObject *result_item(char v)      { return PyBytes_FromStringAndSize(&v, 1); }

    template <class... A>
    inline PyObject *make_result(A &&...a) noexcept
    {
        if constexpr (sizeof...(A) == 1) {
            return result_item(a...);                       /* single value: no tuple */
        } else {
            PyObject *items[] = { result_item(a)... };      /* brace-init: left-to-right order */
            constexpr Py_ssize_t n = sizeof...(A);
            PyObject *t = PyTuple_New(n);
            if (t == nullptr) {
                for (Py_ssize_t i = 0; i < n; i++) { Py_XDECREF(items[i]); }
                return nullptr;
            }
            for (Py_ssize_t i = 0; i < n; i++) { PyTuple_SET_ITEM(t, i, items[i]); }
            return t;
        }
    }


    /**
     * @brief Per-routine context: everything f2py's generator derived from the .pyf declaration.
     *
     * It derives the flavored routine name (`"saxpy"`), whether an argument is a required positional
     * or an optional keyword (its position relative to `'|'`), and its f2py ordinal (`"2nd argument"`
     * / `"4th keyword"`).
     *
     * @tparam Mod  The module identity trait described at the top of this file; supplies the
     *              `pyname` every message is raised under and the `libname` bounding dimensions.
     * @tparam T    The wrapper's scalar flavor (f32/f64/c64/c128).
     */
    template <class Mod, class T>
    class Ctx {
    public:
        /**
         * @param prefix  Explicit routine-name prefix for the irregularly named families -- pass
         *                the matching trait, e.g. BLAS's `tchar_fn<T>()` so `nrm2` becomes
         *                `scnrm2` at `T = c64`.  Regular routines use the two-argument overload
         *                below.
         * @param name    Unprefixed routine name (`"axpy"`, `"nrm2"`).
         * @param pyfmt   PyArg_ParseTupleAndKeywords format without the `:name` suffix; the
         *                position of `'|'` defines how many leading arguments are required.
         * @param kwlist  Null-terminated argument-name list, in signature order.  Must outlive
         *                the Ctx (the wrappers use a function-local static).
         */
        constexpr Ctx(const char *prefix, const char *name, const char *pyfmt, const char *const *kwlist)
            : rout_{}, qual_{}, kwlist_(kwlist), nreq_(0)
        {
            /* Both names are laid down here so the ctor is constexpr and the wrappers' `static
             * constexpr Ctx` needs no run-time initialization.  @p pyfmt is *consumed* rather
             * than stored: all the machinery ever wants from it is where the '|' sits, so only
             * that count survives.  A name longer than its buffer overruns the array -> build
             * error. */
            std::size_t i = 0;
            for (const char *p = prefix; *p; ++p) { rout_[i++] = *p; }
            for (const char *p = name;   *p; ++p) { rout_[i++] = *p; }
            std::size_t j = 0;
            for (const char *p = Mod::pyname; *p; ++p) { qual_[j++] = *p; }
            qual_[j++] = '.';
            for (std::size_t k = 0; k < i;    ++k)     { qual_[j++] = rout_[k]; }
            const char *bar = pyfmt;
            while (*bar && *bar != '|') { ++bar; }
            nreq_ = static_cast<int>(bar - pyfmt);
        }

        /** @brief Regular naming: prepends the flavor letter, `"axpy"` -> `saxpy`/.../`zaxpy`. */
        constexpr Ctx(const char *name, const char *pyfmt, const char *const *kwlist) : Ctx(flavor<T>(), name, pyfmt, kwlist) {}

        /** @brief The keyword list, null-terminated and in signature order. */
        const char *const *kws() const noexcept { return kwlist_; }
        /** @brief Count of required (before-'|') arguments; drives ordinals and missing-arg checks. */
        int nreq() const noexcept { return nreq_; }
        /** @brief Module-qualified routine name every message is raised under, "_fblas.daxpy". */
        const char *qualname() const noexcept { return qual_; }

        /**
         * @brief Position of @p kw in the kwlist (linear scan; the lists are < 12 entries).
         *
         * Called only on cold paths -- to format the "Nth keyword" ordinal of an error message
         * (see scalar()/check_scalar()/checked()).
         */
        int index(const char *kw) const noexcept
        {
            int i = 0;
            while (kwlist_[i] && strcmp(kwlist_[i], kw) != 0) { i++; }
            /**
             * A miss means the kwlist and an argument name spelled in a macro are out of sync
             * Fail loudly instead of letting the end sentinel silently classify the argument
             * as an optional keyword.
             * */
            assert(kwlist_[i] != nullptr && "argument name missing from kwlist");
            return i;
        }

        /** @brief Like index() but returns -1 for an unknown name (used to flag unexpected kwargs). */
        int index_opt(const char *kw) const noexcept
        {
            int i = 0;
            while (kwlist_[i] && strcmp(kwlist_[i], kw) != 0) { i++; }
            return kwlist_[i] ? i : -1;
        }

        /**
         * @brief Convert one scalar argument, unconditionally.
         *
         * The provided-or-default branch lives in SCALAR_OPT (an omitted optional never reaches
         * this function); here the object is always converted.  On conversion failure raises
         * `_fblas.saxpy() 4th keyword (incx) can't be converted to int`.
         *
         * @param v    Conversion target; written only on success.
         * @param kw   Argument name as it appears in the kwlist; on failure its kwlist position
         *             (the "4th keyword" ordinal) is looked up lazily via `index(kw)`.
         * @return true on success; false with an exception set.
         */
        template <class V>
        bool scalar(V *v, PyObject *obj, const char *kw) const noexcept
        {
            conv r = from_pyobj(v, obj);
            if (r == conv::ok) { return true; }
            if (r == conv::fail_msg) {
                /* f2py's policy: keep an already-set exception's type, replace its message.
                 * index(kw) runs only on this cold failure path -- the ordinal exists purely for
                 * the message, so there is no reason to compute or cache it on the hot path. */
                PyObject *err = PyErr_Occurred();
                char pk[32], msg[128];
                poskind(pk, sizeof pk, index(kw));
                int n = snprintf(msg, sizeof msg, "%s() %s (%s) can't be converted to %s",
                                qualname(), pk, kw, conv_name<V>());
                assert(n > 0 && static_cast<size_t>(n) < sizeof msg && "msg buffer too small");
                (void)n;
                PyErr_SetString(err ? err : PyExc_TypeError, msg);
            }
            return false;
        }

        /**
         * @brief Raise a scalar `check(...)` failure as ValueError; no-op when @p ok.
         *
         * On failure raises `(incx>0||incx<0) failed for 4th keyword incx: saxpy:incx=0`, or
         * `(trans=='N'||trans=='T') failed for 1st keyword trans: gels:trans='X'` when the
         * argument is an option letter.  f2py printed those with `'%c'` rather than as a number,
         * and a letter reported as 88 helps nobody.
         *
         * @param tcheck  The check expression's source text (the CHECK macro stringizes it).
         * @param val     The offending value, printed after `name=`.  Integers widen to long
         *                long so ILP64 CBLAS_INT values print correctly on LLP64 platforms.
         */
        template <class V>
        bool check_scalar(bool ok, const char *tcheck, const char *kw, V val) const noexcept
        {
            if (ok) { return true; }
            char pk[32];
            poskind(pk, sizeof pk, index(kw));
            if constexpr (std::is_same_v<V, char>) {
                PyErr_Format(PyExc_ValueError, "(%s) failed for %s %s: %s:%s='%c'",
                             tcheck, pk, kw, rout_, kw, val);
            }
            else {
                PyErr_Format(PyExc_ValueError, "(%s) failed for %s %s: %s:%s=%lld",
                             tcheck, pk, kw, rout_, kw, static_cast<long long>(val));
            }
            return false;
        }

        /**
         * @brief Raise an array `check(...)` failure as ValueError; no-op when @p ok.
         *
         * On failure raises `(shape(a,0)==shape(a,1)) failed for 2nd argument a`.
         */
        bool check_array(bool ok, const char *tcheck, const char *name) const noexcept
        {
            if (ok) { return true; }
            char pk[32];
            poskind(pk, sizeof pk, index(name));
            PyErr_Format(PyExc_ValueError, "(%s) failed for %s %s", tcheck, pk, name);
            return false;
        }

        /** @brief `as_in` (f2py `intent(in)`) with owning result and f2py's failure message.
         *         @tparam V is the element type, defaulting to the wrapper's flavor; LAPACK also
         *         takes integer arrays alongside a float one (`getrs`'s pivots). */
        template <class V = T>
        py_ref in(PyObject *o, int ndim, const char *name) const noexcept
        {
            PyObject *orig = nullptr;
            /* two statements: orig must be written before it is read as an argument */
            PyArrayObject *arr = as_in<V>(o, ndim, &orig);
            return checked(arr, orig, name);
        }
        /** @brief `as_inout` (f2py `intent(in,out[,copy])`) with owning result and failure message.
         *         @tparam V is the element type, defaulting to the wrapper's flavor; LAPACK also
         *         hands back integer arrays (`gelsy`'s column pivots). */
        template <class V = T>
        py_ref inout(PyObject *o, int ndim, bool overwrite, const char *name) const noexcept
        {
            PyObject *orig = nullptr;
            PyArrayObject *arr = as_inout<V>(o, ndim, overwrite, &orig);
            return checked(arr, orig, name);
        }
        /** @brief `zeros_vec` with owning result, for optional output vectors not supplied.
         *         Guarded by the same CBLAS_INT_MAX certificate as acquired arrays, so
         *         `len()` stays lossless on internally allocated arrays too. */
        py_ref zeros(npy_intp n) const noexcept
        {
            return zeros_as<T>(n);
        }
        /** @brief zeros() with an explicit element type, for buffers whose type is not the
         *         wrapper's flavor (LAPACK's real `rwork`, integer `bwork`). */
        template <class V>
        py_ref zeros_as(npy_intp n) const noexcept
        {
            if (n > (npy_intp)CBLAS_INT_MAX) {
                PyErr_Format(PyExc_OverflowError,
                            "%s(): required output size %lld exceeds the %s integer"
                            " limit (%lld)",
                            qualname(), static_cast<long long>(n), Mod::libname,
                            static_cast<long long>(CBLAS_INT_MAX));
                return py_ref(nullptr);
            }
            return py_ref(zeros_vec<V>(n));
        }
        /** @brief Two-dimensional variant of zeros(), Fortran-ordered, same certificate. */
        py_ref zeros(npy_intp m, npy_intp n) const noexcept
        {
            return zeros_as<T>(m, n);
        }
        /** @brief Two-dimensional zeros_as(), Fortran-ordered, same certificate. */
        template <class V>
        py_ref zeros_as(npy_intp m, npy_intp n) const noexcept
        {
            if (m > (npy_intp)CBLAS_INT_MAX || n > (npy_intp)CBLAS_INT_MAX) {
                PyErr_Format(PyExc_OverflowError,
                            "%s(): required output shape (%lld, %lld) exceeds the %s"
                            " integer limit (%lld)",
                            qualname(), static_cast<long long>(m), static_cast<long long>(n),
                            Mod::libname, static_cast<long long>(CBLAS_INT_MAX));
                return py_ref(nullptr);
            }
            return py_ref(zeros_mat<V>(m, n));
        }

    private:
        /**
         * @brief Format the f2py position string, `"2nd argument"` / `"4th keyword"`.
         *
         * f2py numbers required positionals and optional keywords as two separate sequences
         * (gemv's `beta`, 4th in the signature, is the "1st keyword").  The generated modules
         * never exceed "7th", so the simple 1st/2nd/3rd/Nth suffix rule is exact.
         */
        void poskind(char *buf, size_t size, int i) const noexcept
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
        py_ref checked(PyArrayObject *a, PyObject *orig, const char *name) const noexcept
        {
            if (a == nullptr) {
                Py_XDECREF(orig);
                if (!PyErr_Occurred()) {
                    char pk[32];
                    poskind(pk, sizeof pk, index(name));
                    PyErr_Format(PyExc_TypeError, "%s.%s: failed to create array from the %s `%s`",
                                Mod::pyname, qualname(), pk, name);
                }
                return py_ref(nullptr);
            }
            for (int i = 0; i < PyArray_NDIM(a); i++) {
                if (PyArray_DIM(a, i) > (npy_intp)CBLAS_INT_MAX) {
                    char pk[32];
                    poskind(pk, sizeof pk, index(name));
                    PyErr_Format(PyExc_OverflowError,
                                "%s() %s (%s): dimension %d with size %lld exceeds"
                                " the %s integer limit (%lld)",
                                qualname(), pk, name, i,
                                static_cast<long long>(PyArray_DIM(a, i)),
                                Mod::libname, static_cast<long long>(CBLAS_INT_MAX));
                    Py_DECREF(a);
                    Py_XDECREF(orig);
                    return py_ref(nullptr);
                }
            }
            return py_ref(a, orig);
        }

        char rout_[16];               /**< Flavored routine name alone, "daxpy". */
        char qual_[32];               /**< Module-qualified, "_fblas.daxpy" -- what messages use. */
        const char *const *kwlist_;
        int nreq_;
    };


    /**
     * @brief Supplies each argument's raw Python object to the wrapper by name, so an argument is
     *        named exactly once (in the kwlist).  `raw("y")` returns the object passed for `y`
     *        positionally or by keyword, or nullptr when it was omitted.
     *
     * `parse()` performs the up-front structural validation, reproducing CPython's argument-parsing
     * error precedence: too-many-arguments, then missing-required, then given-by-name-and-position,
     * then unexpected-keyword (see parse() for why these must be separate passes).  Conversion,
     * size checks and array acquisition happen afterwards, in the SCALAR_ and ARRAY_ declaration
     * macros, which read from this cursor and `Ctx`.
     *
     * @tparam C  The wrapper's `Ctx` type, deduced by PARSE_ARGS() through class template
     *            argument deduction so no wrapper has to name it.
     */
    template <class C>
    class ArgCursor {
    public:
        ArgCursor(PyObject *args, PyObject *kwds, const C &ctx)
            : args_(args), kwds_(kwds), ctx_(ctx) {}

        /** @brief CPython-order structural validation; false with an exception set on failure. */
        bool parse() const noexcept
        {
            const char *const *kw = ctx_.kws();
            const char *fn = ctx_.qualname();
            Py_ssize_t npos = args_ ? PyTuple_GET_SIZE(args_) : 0;
            Py_ssize_t nkw = kwds_ ? PyDict_Size(kwds_) : 0;
            int nreq = ctx_.nreq();
            int nnames = 0;
            while (kw[nnames]) { nnames++; }

            /* CPython's PyArg validates in *separate passes*, not one interleaved loop -- so a
             * clash at an earlier index does not pre-empt a missing-required at a later one
             * (`ddot(x, x=x)` reports missing 'y', not the 'x' clash).  The order is:
             *   1. too many arguments;  2. missing required;
             *   3. argument given by name and position;  4. unexpected keyword. */

            /* 1. too many arguments: the count is positional + keyword *total*, so an extra,
             * duplicate, or unexpected keyword tips a fully-positional call over (e.g. rotg's
             * `OO|` reports "takes at most 2 (3 given)" for `drotg(a, b, foo=1)`, not an
             * unexpected-keyword error).  Keyword functions always word this "at most". */
            if (npos + nkw > nnames) {
                PyErr_Format(PyExc_TypeError, "%s() takes at most %d argument%s (%zd given)",
                             fn, nnames, nnames == 1 ? "" : "s", npos + nkw);
                return false;
            }

            /* 2. missing required arguments (signature order) */
            for (int i = 0; i < nreq; i++) {
                bool given = (i < npos) || (kwds_ && PyDict_GetItemString(kwds_, kw[i]) != nullptr);
                if (!given) {
                    PyErr_Format(PyExc_TypeError,
                                 "%s() missing required argument '%s' (pos %d)", fn, kw[i], i + 1);
                    return false;
                }
            }

            /* 3. arguments given by both name and position (signature order) */
            for (int i = 0; i < nnames && i < npos; i++) {
                if (kwds_ && PyDict_GetItemString(kwds_, kw[i]) != nullptr) {
                    PyErr_Format(PyExc_TypeError,
                                 "argument for %s() given by name ('%s') and position (%d)",
                                 fn, kw[i], i + 1);
                    return false;
                }
            }

            /* 4. unexpected keyword arguments */
            if (kwds_) {
                PyObject *key, *val;
                Py_ssize_t pos = 0;
                while (PyDict_Next(kwds_, &pos, &key, &val)) {
                    const char *k = PyUnicode_AsUTF8(key);
                    if (k == nullptr) { return false; }
                    if (ctx_.index_opt(k) < 0) {
                        PyErr_Format(PyExc_TypeError,
                                     "%s() got an unexpected keyword argument '%s'", fn, k);
                        return false;
                    }
                }
            }
            return true;
        }

        /**
         * @brief Raw (borrowed) object for @p name: its positional slot, else the keyword,
         *        else nullptr when the argument was omitted (an optional's default applies).
         *
         * Safe to call only after parse() succeeded (no position/keyword clash remains).
         */
        PyObject *raw(const char *name) const noexcept
        {
            return at(ctx_.index(name), name);
        }

        /**
         * @brief raw() for a result the caller is not offered a way to supply.
         *
         * Most routines produce their outputs outright: `geqp3` hands back `jpvt`, `tau` and
         * `work` without accepting any of them.  Those names are not in the kwlist, so the
         * supplied-or-not question has no argument to read and the answer is always "not
         * supplied".  Unlike raw() this is not a kwlist/macro-name typo, so it does not assert.
         */
        PyObject *raw_opt(const char *name) const noexcept
        {
            int idx = ctx_.index_opt(name);
            return idx < 0 ? nullptr : at(idx, name);
        }

    private:
        PyObject *at(int idx, const char *name) const noexcept
        {
            Py_ssize_t npos = args_ ? PyTuple_GET_SIZE(args_) : 0;
            if (idx < npos) { return PyTuple_GET_ITEM(args_, idx); }
            return kwds_ ? PyDict_GetItemString(kwds_, name) : nullptr;
        }

        PyObject *args_;
        PyObject *kwds_;
        const C &ctx_;
    };

}  // namespace wrapper


/**
 * ====================================================================================
 * Argument-handling vocabulary.  A wrapper parses its arguments with `PARSE_ARGS()` and then
 * declares each one, in processing order, with one macro line per argument.  Two names are
 * referenced by convention throughout: `P`, the ArgCursor that supplies each argument's raw
 * Python object by name, and `ctx`, the per-routine context (flavored routine name, keyword
 * list, required-argument count) that formats error messages.  Conventions:
 *
 * - A bare array name (`x`, `y`, `a`, `b`) is the *acquired* array, held by an owning `py_ref`
 *   that releases its reference on every exit path (so there is no manual cleanup).
 * - A bare lower-case scalar name (`n`, `incx`, `alpha`, ...) is the converted C scalar.
 * - A check is an ordinary C expression, stringized into its own error message; `len()`,
 *   `shape()` and `abs()` are real functions (see above) so the expressions compile as written.
 *
 * Naming: an argument-declaration macro is `<KIND>_<VARIANT>`.  Scalars: `SCALAR_REQ` (required),
 * `SCALAR_OPT` (optional, with a default), and `SCALAR_FLAG` (an optional `overwrite_*` flag using
 * strict `__index__` conversion).  Arrays: `ARRAY_IN` (read-only input), `ARRAY_INOUT` (read-write
 * in/out), `ARRAY_OUT` (a result the routine hands back), and `ARRAY_HIDDEN` (scratch that never
 * leaves).  Everything else is a verb for an action: `PARSE_ARGS`, `CHECK`, `CHECKARRAY`,
 * `RETURN`; LAPACK adds `CALLABLE_SELECT`.
 *
 * The array macros divide on **what happens to the array**, not on how it is acquired.  `IN` and
 * `INOUT` come from the caller; `OUT` goes to the caller; `HIDDEN` goes nowhere.  Whether the
 * caller may pre-supply an `ARRAY_OUT` buffer varies by routine (gemv's `y` yes, geqp3's `tau`
 * no) and is an implementation detail of that macro, not a reason to reclassify the array: the
 * name has to tell a reader reconstructing the Python signature which values come back.
 * `ARRAY_IN` and `ARRAY_INOUT` take the element type as their leading argument, as `ARRAY_OUT`'s
 * allocation expression and `ARRAY_HIDDEN` do, because LAPACK passes integer and real-counterpart
 * arrays alongside a routine's own flavor.
 * ====================================================================================
 */

/** @brief Declare `ArgCursor P` over `args`/`kwds` and run CPython-order structural parsing. */
#define PARSE_ARGS() \
    wrapper::ArgCursor P(args, kwds, ctx); \
    if (!P.parse()) { return nullptr; }

/**
 * @brief Acquire a required array argument, declaring the owning `py_ref <name>`.
 *
 * `ARRAY_IN` is a read-only input (`intent(in)`): the caller's array is used as-is when it
 * already fits (dtype, rank, Fortran-contiguous, aligned), otherwise a converted copy is made.
 * @p type is the element type, spelled out as `ARRAY_HIDDEN` spells it, because a routine can
 * take an integer array next to its float ones (`getrs`'s `piv` beside its `lu`).
 * `ARRAY_INOUT` is read-write (`intent(in,out)`): with @p overwrite the caller's buffer is
 * reused in place when it fits, otherwise a fresh writable copy is made.  The acquired array is
 * what the wrapper returns.  Two statements (declaration + failure return): block level only.
 */
#define ARRAY_IN(type, name, ndim) \
    py_ref name = ctx.template in<type>(P.raw(#name), ndim, #name); \
    if (!name) { return nullptr; }

#define ARRAY_INOUT(type, name, ndim, overwrite) \
    py_ref name = ctx.template inout<type>(P.raw(#name), ndim, overwrite, #name); \
    if (!name) { return nullptr; }

/**
 * @brief Acquire a result the wrapper produces and hands back to the caller.
 *
 * A supplied array is taken in/out and returned; an omitted one is created fresh, filled by the
 * Fortran call, and returned.  This is f2py's `intent(in,out)` optional output (gemv's `y`,
 * syr2/ger's `a`); the macro hides the supplied-vs-omitted branch so the wrapper body stays flat.
 *
 * Most LAPACK results cannot be supplied at all -- `geqp3` hands back `jpvt`, `tau` and `work`
 * without accepting any of them.  Those names are not among the wrapper's arguments, so the
 * omitted branch is the only reachable one and @p def is always what the caller gets back.
 *
 * Contrast ARRAY_HIDDEN, which is for buffers that never leave the wrapper.
 *
 * @param name       output array name; a `py_ref name` is declared holding the acquired-or-
 *                   allocated array -- the value the wrapper hands to RETURN.
 * @param ndim       required rank when the caller supplies the array.
 * @param overwrite  bool: reuse the caller's supplied array in place instead of copying it.
 * @param def        what the array defaults to when the caller omits it: a freshly allocated,
 *                   zero-filled result buffer sized for this routine's output -- `ctx.zeros(len)`
 *                   for a vector, `ctx.zeros(m, n)` for a matrix, `ctx.zeros_as<V>(...)` when the
 *                   element type is not the wrapper's flavor. To avoid garbage values, because
 *                   the routine can accumulate into it (gemv's `beta*y`, syr's
 *                   `a += alpha*x*x'` etc.), it is zero initialized.
 */
#define ARRAY_OUT(name, ndim, overwrite, def) \
    py_ref name = (P.raw_opt(#name) == nullptr) ? (def) : ctx.inout(P.raw_opt(#name), ndim, overwrite, #name); \
    if (!name) { return nullptr; }

/**
 * @brief Declare a zero-filled scratch buffer that never leaves the wrapper -- f2py's
 *        `intent(hide)` work arrays: `work`, `rwork`, `iwork`, `bwork`.
 *
 * `ARRAY_HIDDEN(T, work, lwork)` for a vector, `ARRAY_HIDDEN(R, rwork, 5 * minmn)` for a real
 * buffer beside a complex flavor.  @p type is explicit rather than defaulted to the wrapper's
 * flavor because these regularly use the real counterpart of a complex flavor (`rwork`) or the
 * ABI-width integer that backs Fortran `LOGICAL` (`bwork`).
 *
 * A result the routine hands back is ARRAY_OUT, even when the caller cannot supply it: the two
 * macros say whether the array escapes, which is what a reader reconstructing the wrapper's
 * Python signature needs to know.  Some routines return a buffer LAPACK calls `work` (`gees`,
 * `gges`, `gelss`, `geqp3` all report the optimal size in `work[0]`); those are ARRAY_OUT.
 *
 * Requires the enclosing wrapper's `ctx`; declares an owning `py_ref name`.
 */
#define ARRAY_HIDDEN(type, name, ...) \
    py_ref name = ctx.template zeros_as<type>(__VA_ARGS__); \
    if (!name) { return nullptr; }

/**
 * @brief Declare the scalar `<name>` of C type @p type and assign it: convert the caller's
 *        object when the argument is supplied, otherwise apply the default expression @p def.
 *
 * The default is cast to @p type.  Conversion follows the permissive `from_pyobj` rules (a
 * float truncates to int, a complex contributes its real part, a sequence its element 0); on
 * failure it raises the "can't be converted to" message and returns nullptr.
 */
#define SCALAR_OPT(type, name, def) \
    type name; \
    do { \
        PyObject *name##_raw = P.raw(#name); \
        if (name##_raw == nullptr) { name = static_cast<type>(def); } \
        else if (!ctx.scalar(&name, name##_raw, #name)) { return nullptr; } \
    } while (0)

/**
 * @brief Declare and convert a required scalar `<name>` of C type @p type (no default).
 *
 * For required positionals (`alpha`); parse() guarantees the object is present, so it is
 * converted unconditionally through the permissive `from_pyobj` path, raising the
 * "can't be converted to" message on failure.
 */
#define SCALAR_REQ(type, name) \
    type name; \
    do { \
        if (!ctx.scalar(&name, P.raw(#name), #name)) { return nullptr; } \
    } while (0)

/**
 * @brief Declare and assign an `int` flag `<name>` (default 0) -- the `overwrite_*` controls,
 *        which are not Fortran arguments but tell the wrapper whether it may write in place.
 *
 * A flag uses strict `__index__`-based conversion (see from_pyobj_index), which rejects the
 * float/str/None that SCALAR_OPT's permissive path would coerce.  Place SCALAR_FLAG immediately
 * after PARSE_ARGS(), ahead of the body scalars, so a bad flag is reported before a bad body
 * scalar (a fixed error precedence).
 */
#define SCALAR_FLAG(name) \
    int name = 0; \
    do { \
        PyObject *name##_raw = P.raw(#name); \
        if (name##_raw != nullptr && !wrapper::from_pyobj_index(&name, name##_raw)) { return nullptr; } \
    } while (0)

/**
 * @brief Validate a converted scalar against a `check(...)` C-expression.
 *
 * `CHECK(incx != 0, incx)` stringizes the expression into the message
 * `(incx != 0) failed for 4th keyword incx: daxpy:incx=0`.  @p name is the argument the check
 * is reported against (not necessarily a variable in @p expr, e.g. the n-bound checks report
 * `n`), and supplies the printed value.  On failure raises ValueError and returns nullptr.
 */
#define CHECK(expr, name) \
    do { if (!ctx.check_scalar((expr), #expr, #name, name)) { return nullptr; } } while (0)

/**
 * @brief Validate an acquired array against a `check(...)` C-expression.
 *
 * The array analog of CHECK: `CHECKARRAY(shape(a, 0) == shape(a, 1), a)` stringizes the
 * expression into `(shape(a, 0) == shape(a, 1)) failed for 2nd argument a`.  @p name is the
 * array the check is reported against (a `py_ref` local).  On failure raises ValueError and
 * returns nullptr.
 */
#define CHECKARRAY(expr, name) do { if (!ctx.check_array((expr), #expr, #name)) { return nullptr; } } while (0)

/** @brief Pack the wrapper's results and `return` them (see make_result): `RETURN(y)`,
 *         `RETURN(x, y)`, `RETURN(idx)`.  Compute values into locals first; no expressions here. */
#define RETURN(...) return wrapper::make_result(__VA_ARGS__)

/**
 * @brief Emit the four method-table rows of one routine -- `s<name>`, `d<name>`, `c<name>`,
 *        `z<name>` -- each bound to the wrapper template instantiated at the matching scalar
 *        type.  Expand inside `blas::capi` / `lapack::capi`, where the wrappers and the width
 *        aliases resolve unqualified.
 */
#define FAMILY(name) \
    {"s" #name, (PyCFunction)(void (*)())name<f32>,  METH_VARARGS | METH_KEYWORDS, nullptr}, \
    {"d" #name, (PyCFunction)(void (*)())name<f64>,  METH_VARARGS | METH_KEYWORDS, nullptr}, \
    {"c" #name, (PyCFunction)(void (*)())name<c64>,  METH_VARARGS | METH_KEYWORDS, nullptr}, \
    {"z" #name, (PyCFunction)(void (*)())name<c128>, METH_VARARGS | METH_KEYWORDS, nullptr}

/**
 * @brief One explicitly named method-table row, for the irregular families where the Python
 *        name is not flavor letter + template name (`scnrm2` binds `nrm2<c64>`).
 *        The name spelled here must agree with the prefix trait passed to the wrapper's
 *        Ctx -- that is what the parity suite's message comparisons pin down.
 */
#define ROW(pyname, name, T) \
    {#pyname, (PyCFunction)(void (*)())name<T>, METH_VARARGS | METH_KEYWORDS, nullptr}

/**
 * @brief Like ROW but for a wrapper template with two type parameters `name<T, A>`.
 *        Used by the real-scalar `scal` overload (`csscal` = `scal<c64, f32>`): the second
 *        parameter is the scalar type, which also distinguishes it from the regular `scal<T>`.
 */
#define ROW2(pyname, name, T, A) \
    {#pyname, (PyCFunction)(void (*)())name<T, A>, METH_VARARGS | METH_KEYWORDS, nullptr}
