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
 *    deliberately preserving its (often overly) permissive coercions.
 *
 * 2. Array acquisition (`as_in` / `as_inout` / `zeros_vec`) built on the contemporary
 *    `PyArray_FromAny` instead of f2py's fortranobject machinery, plus the owning reference
 *    `py_ref`, copy/in-place behavior matching f2py.
 *
 * 3. `Ctx`, the per-routine context; the flavored routine name, each argument's required vs.
 *     optional kind and ordinal ("4th keyword") error messages, and inlined the resulting
 *     message strings at every call site.
 *
 * @note The extension is built from several .cpp files, so numpy's C-API function table must
 *       be shared: every .cpp file defines `PY_ARRAY_UNIQUE_SYMBOL` to the same name *before*
 *       including this header, and every file except the module's `*_module.cpp` (which calls
 *       `import_array()` once, in module init) also defines `NO_IMPORT_ARRAY`.
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
#include "wrapper_types.hpp"      /* f32, f64, c64, c128, real_of, is_complex, flavor */

namespace wrapper {

    /**
     * Exception policy: these wrappers raise builtin exception types.
     *
     * ValueError for argument-check failures, TypeError for conversion fallbacks, and
     * OverflowError for range-check failures.
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
     * @brief Return @p arr reshaped to exactly @p ndim dimensions, as a view.
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
     * @param arr   Array to reshape.  Its reference is taken over on every path: handed back
     *              as the return value when the rank already matches, moved into @p orig when
     *              a view is made, and released when the reshape fails.
     * @param ndim  Rank the result must have.
     * @param orig  Receives the owning reference to the pre-reshape array when a reshape
     *              happened, nullptr otherwise.
     * @return      New reference to the rank-@p ndim working view (the input itself when the
     *              rank already matches), or nullptr.
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
     * @brief Acquire a writable input and return to Python (`intent(in,out)`).
     *
     * With @p overwrite, NumPy reuses the caller's array when it fits (the in-place path);
     * otherwise, or whenever the array doesn't fit, a fresh writable copy is made (f2py's
     * INTENT_COPY semantics). Rank mismatches are reinterpreted per fix_rank(); on the
     * in-place path the result is then a zero-copy *view* of the caller's buffer.
     *
     * @note We call `PyArray_FromAny` directly, not the `PyArray_FROMANY` macro: the macro ORs
     *       `C_CONTIGUOUS` into the flags whenever `ENSURECOPY` is set, which would conflict
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
     *         not supplied. */
    template <class T>
    inline PyArrayObject *zeros_mat(npy_intp m, npy_intp n) noexcept
    {
        npy_intp dims[2] = {m, n};
        return (PyArrayObject *)PyArray_ZEROS(2, dims, npy_type<T>(), 1);
    }


    /* ---------------- scalar converters, ported from the f2py-generated module ---------------- */

    /**
     * @brief Converter result enum to build error messages on failure avoiding hot path bloat.
     *
     * Converters report *how* they failed and `Ctx::scalar` formats lazily:
     *
     * - `ok`:       value written.
     * - `fail`:     a specific exception is already set and final (e.g. OverflowError with its
     *               own text); do not touch it.
     * - `fail_msg`: apply f2py's message policy, that is, keep an already-set exception's *type* but
     *               replace its message with the argument-specific text, or raise TypeError.
     */
    enum class conv { ok, fail, fail_msg };

    /**
     * @brief RAII `Py_EnterRecursiveCall`, for the converters' descent into element 0.
     *
     * The scalar converters below reproduce f2py's permissive coercion: an object that is not a
     * number but is a sequence contributes its element 0, tried again by calling the converter
     * recursively. Nesting is normally two or three deep, but a self-referential container,
     * `a = []; a.append(a)`, makes element 0 the sequence itself and the descent unbounded.
     * That recursion is pure C and pushes no Python frame and simply overflows. (f2py segfaults
     * on this input as well)
     *
     * `Py_EnterRecursiveCall` opts the descent into the interpreter's own depth accounting, so
     * the limit applies and `RecursionError` is raised instead. Pairing it with
     * `Py_LeaveRecursiveCall` by hand is error-prone here because the converters return from
     * several places, hence the destructor. Guard only the recursive step: the common,
     * non-recursive conversion stays free of the counter.
     */
    class recursion_guard {
        bool entered_;
    public:
        explicit recursion_guard(const char *where) noexcept : entered_(Py_EnterRecursiveCall(where) == 0) {}
        ~recursion_guard() { if (entered_) { Py_LeaveRecursiveCall(); } }
        recursion_guard(const recursion_guard &) = delete;
        recursion_guard &operator=(const recursion_guard &) = delete;
        /** @brief True when the descent may proceed; false with `RecursionError` set. */
        explicit operator bool() const noexcept { return entered_; }
    };

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
            recursion_guard descent(" while converting an argument");
            if (!descent) { Py_DECREF(tmp); return conv::fail; }   /* RecursionError is set */
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
            recursion_guard descent(" while converting an argument");
            if (!descent) { Py_DECREF(tmp); return conv::fail; }   /* RecursionError is set */
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
                recursion_guard descent(" while converting an argument");
                if (!descent) { Py_DECREF(tmp); return conv::fail; }   /* RecursionError is set */
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
            int kind = PyArray_TYPE(arr);
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
                recursion_guard descent(" while converting an argument");
                if (!descent) { Py_DECREF(tmp); return conv::fail; }   /* RecursionError is set */
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
     * working array plus the caller-shaped original. See private p_/orig_ below).
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

        // More move-only RAII stuff: copies are deleted, moves steal the references, destructor releases.
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
         * Returns the working array, or, after a fix_rank() reshape, the caller-shaped original
         * instead. The handle is left empty, so the destructor then does nothing.
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
     * The check expressions are stringized into error messages, so they should stay readable as prose
     *
     * `len`/`shape` return CBLAS_INT so dimensions can flow to the Fortran calls with no cast at any
     * call site. The static_cast here is lossless since every py_ref came through `Ctx::checked` /
     * `Ctx::zeros`, which reject any array with a dimension above CBLAS_INT_MAX.
     *
     * @note `abs` intentionally takes and returns `npy_intp`, so size products like
     *       `(n - 1) * abs(incx)` are computed 64-bit. f2py evaluated them in C int, which can
     *       overflow before the comparison. The wrappers live in `blas::capi` / `lapack::capi`
     *       and their enclosing namespace re-exports this overload with a using-declaration, so
     *       unqualified `abs` stops there and never reaches the C++ library's `::abs(int)`.
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
     * `RETURN(...)`. Values are computed into locals first (no expressions squeezed inside RETURN)
     */
    inline PyObject *result_item(py_ref &r)   { return r.release(); }
    inline PyObject *result_item(f32 v)       { return to_pyobj(v); }
    inline PyObject *result_item(f64 v)       { return to_pyobj(v); }
    inline PyObject *result_item(c64 v)       { return to_pyobj(v); }
    inline PyObject *result_item(c128 v)      { return to_pyobj(v); }
    inline PyObject *result_item(long long v) { return PyLong_FromLongLong(v); }

    /** A returned LAPACK option letter, as f2py's `Py_BuildValue("c")`. */
    inline PyObject *result_item(char v)      { return PyBytes_FromStringAndSize(&v, 1); }

    /**
     * @brief Steal and pack @p n owned references into a tuple.
     *
     * A null entry means the matching `result_item` ran out of memory; the tuple is abandoned in
     * that case.  Not a template, so it is emitted once rather than per distinct signature.
     *
     * @param items  @p n owned references, some possibly null.  All are consumed, on both the
     *               success and the failure path.
     * @param n      Number of entries in @p items.
     * @return       New reference to the tuple, or nullptr with an exception set.
     */
    inline PyObject *tuple_steal(PyObject **items, Py_ssize_t n) noexcept
    {
        bool complete = true;
        for (Py_ssize_t i = 0; i < n; i++) {
            if (items[i] == nullptr) { complete = false; }   /* a converter ran out of memory */
        }
        PyObject *t = complete ? PyTuple_New(n) : nullptr;
        if (t == nullptr) {
            for (Py_ssize_t i = 0; i < n; i++) { Py_XDECREF(items[i]); }
            return nullptr;
        }
        for (Py_ssize_t i = 0; i < n; i++) { PyTuple_SET_ITEM(t, i, items[i]); }
        return t;
    }

    /**
     * @brief Pack a wrapper's results: the bare value for one, a tuple for several.
     *
     * The braced-init-list ensures that it runs every `result_item` exactly once, left to right.
     * Function arguments have no such ordering, so the pack cannot be expanded into a call instead.
     *
     * @param a  The results in declaration order -- `py_ref` arrays, scalars, or `long long`.
     * @return   New reference to the value or tuple, or nullptr with an exception set.
     */
    template <class... A>
    inline PyObject *make_result(A &&...a) noexcept
    {
        if constexpr (sizeof...(A) == 1) {
            return result_item(a...);                       /* single value: no tuple */
        } else {
            PyObject *items[] = { result_item(a)... };      /* ordered, and all of them run */
            return tuple_steal(items, sizeof...(A));
        }
    }


    /**
     * @brief A routine's calling context: everything f2py's generator derived from the `.pyf`
     *        declaration -- its flavored name, keyword list and required-argument count, plus
     *        every message it can raise.
     *
     * It derives the flavored routine name (`"saxpy"`), whether an argument is a required
     * positional or an optional keyword (its position relative to `'|'`), and its f2py ordinal
     * (`"2nd argument"` / `"4th keyword"`).
     *
     * Here template specialization T is deliberately left out to reduce the binary bloat in
     * each template instantiation.
     *
     * `Mod` stays a template parameter and is instantiated exactly once anyway. Storing them
     * as members instead costs unnecessary size increase due to high number of instantiations.
     *
     * @tparam Mod  The module identity trait described at the top of this file; supplies the
     *              `pyname` every message is raised under and the `libname` bounding dimensions.
     */
    template <class Mod>
    class CtxBase {
    public:
        /**
         * @param prefix  Explicit routine-name prefix to accommodate for irregularly named families
         *                pass the matching trait, e.g. BLAS's `tchar_fn<T>()` so `nrm2` becomes
         *                `scnrm2` at `T = c64`.  Regular routines use the two-argument overload
         *                below.
         * @param name    Unprefixed routine name (`"axpy"`, `"nrm2"`).
         * @param pyfmt   `PyArg_ParseTupleAndKeywords` format without the `:name` suffix; the
         *                position of `'|'` defines kwarg separation.
         * @param kwlist  Null-terminated argument-name list, in signature order.
         */
        constexpr CtxBase(const char *prefix, const char *name, const char *pyfmt, const char *const *kwlist)
            : rout_{}, qual_{}, kwlist_(kwlist), nreq_(0)
        {
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

        /** @brief The keyword list, null-terminated and in signature order. */
        const char *const *kws() const noexcept { return kwlist_; }
        /** @brief Count of required (before-'|') arguments; drives ordinals and missing-arg checks. */
        int nreq() const noexcept { return nreq_; }
        /** @brief Module-qualified routine name every message is raised under, "_fblas.daxpy". */
        const char *qualname() const noexcept { return qual_; }

        /**
         * @brief Position of @p kw in the kwlist (linear scan; the lists are < 12 entries), or
         *        -1 when @p kw is not an argument at all.
         *
         * Called only on cold paths to format the "Nth keyword" ordinal of an error message.
         * A derived quantity has no ordinal -- `tbtrs` bounds `b` by the hidden `ldb` -- and
         * reports under its bare name; see poskind().
         */
        int index(const char *kw) const noexcept
        {
            int i = 0;
            while (kwlist_[i] && strcmp(kwlist_[i], kw) != 0) { i++; }
            return kwlist_[i] ? i : -1;
        }

        /**
         * @brief Convert one scalar argument, unconditionally.
         *
         * The provided-or-default branch lives in SCALAR_OPT. Here the object is always converted.
         * On conversion failure raises `_fblas.saxpy() 4th keyword (incx) can't be converted to int`.
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
                PyObject *err = PyErr_Occurred();
                char pk[32], msg[128];
                poskind(pk, sizeof pk, index(kw));
                int n = snprintf(msg, sizeof msg, "%s() %s(%s) can't be converted to %s", qualname(), pk, kw, conv_name<V>());
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
         * argument is an option letter.
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
                PyErr_Format(PyExc_ValueError, "(%s) failed for %s%s: %s:%s='%c'", tcheck, pk, kw, rout_, kw, val);
            }
            else {
                PyErr_Format(PyExc_ValueError, "(%s) failed for %s%s: %s:%s=%lld", tcheck, pk, kw, rout_, kw, static_cast<long long>(val));
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
            PyErr_Format(PyExc_ValueError, "(%s) failed for %s%s", tcheck, pk, name);
            return false;
        }

        /** @brief `as_in` (f2py `intent(in)`) with owning result and f2py's failure message.
         *         @tparam V is the element type, defaulting to the wrapper's flavor or int arrays. */
        template <class V>
        py_ref in(PyObject *o, int ndim, const char *name) const noexcept
        {
            PyObject *orig = nullptr;
            /* two statements: orig must be written before it is read as an argument */
            PyArrayObject *arr = as_in<V>(o, ndim, &orig);
            return checked(arr, orig, name);
        }
        /** @brief `as_inout` (f2py `intent(in,out[,copy])`) with owning result and failure message.
         *         @tparam V is the element type, defaulting to the wrapper's flavor or int arrays. */
        template <class V>
        py_ref inout(PyObject *o, int ndim, bool overwrite, const char *name) const noexcept
        {
            PyObject *orig = nullptr;
            PyArrayObject *arr = as_inout<V>(o, ndim, overwrite, &orig);
            return checked(arr, orig, name);
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
         * @brief Format the f2py position string, `"2nd argument "` / `"4th keyword "`, with the
         *        trailing space its callers splice in; empty for a name with no position (-1).
         */
        void poskind(char *buf, size_t size, int i) const noexcept
        {
            if (i < 0) { buf[0] = '\0'; return; }
            bool req = i < nreq_;
            int num = req ? i + 1 : i - nreq_ + 1;
            const char *suf = num == 1 ? "st" : num == 2 ? "nd" : num == 3 ? "rd" : "th";
            snprintf(buf, size, "%d%s %s ", num, suf, req ? "argument" : "keyword");
        }

        /**
         * @brief Take ownership of a freshly acquired array, rejecting any dimension too large
         *        for `CBLAS_INT`.
         *
         * Every acquired array passes through this, so this is the only `npy_intp` -> `CBLAS_INT`
         * narrowing in the wrappers; `len()` and `shape()` return `CBLAS_INT` on that guarantee.
         *
         * @param a     Array to acquire, or nullptr if acquisition failed.
         * @param orig  Pre-reshape array from `fix_rank`, or nullptr.
         * @param name  Argument name, for the message and its f2py ordinal.
         * @return      Owning reference to @p a, or a null `py_ref` with one of:
         *              - TypeError, when @p a is null and no exception is set.
         *              - OverflowError, naming the axis and its size, when a dimension exceeds
         *                `CBLAS_INT_MAX`.
         */
        py_ref checked(PyArrayObject *a, PyObject *orig, const char *name) const noexcept
        {
            if (a == nullptr) {
                Py_XDECREF(orig);
                if (!PyErr_Occurred()) {
                    char pk[32];
                    poskind(pk, sizeof pk, index(name));
                    PyErr_Format(PyExc_TypeError, "%s.%s: failed to create array from the %s`%s`",
                                Mod::pyname, qualname(), pk, name);
                }
                return py_ref(nullptr);
            }
            for (int i = 0; i < PyArray_NDIM(a); i++) {
                if (PyArray_DIM(a, i) > (npy_intp)CBLAS_INT_MAX) {
                    char pk[32];
                    poskind(pk, sizeof pk, index(name));
                    PyErr_Format(PyExc_OverflowError,
                                "%s() %s(%s): dimension %d with size %lld exceeds"
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
        const char *const *kwlist_;   /**< Argument names in signature order, null-terminated; borrowed so outlives the Ctx. */
        int nreq_;                    /**< Count of leading required arguments, taken from the `'|'` position. */
    };


    /**
     * @brief A CtxBase bound to one scalar flavor: `T` sets the routine-name prefix letter and
     *        the element type of `zeros()`.
     *
     * @tparam Mod  Module traits supplying `pyname` and `libname`.
     * @tparam T    Scalar flavor (f32/f64/c64/c128).
     */
    template <class Mod, class T>
    class Ctx : public CtxBase<Mod> {
    public:
        /** @brief Irregular naming: an explicit prefix, e.g. BLAS's `tchar_fn<T>()`. */
        constexpr Ctx(const char *prefix, const char *name, const char *pyfmt, const char *const *kwlist)
            : CtxBase<Mod>(prefix, name, pyfmt, kwlist) {}

        /** @brief Regular naming: prepends the flavor letter, `"axpy"` -> `saxpy`/.../`zaxpy`. */
        constexpr Ctx(const char *name, const char *pyfmt, const char *const *kwlist)
            : Ctx(flavor<T>(), name, pyfmt, kwlist) {}

        /** @brief `zeros_vec` with owning result, for optional output vectors not supplied.
         *         Guarded by the same CBLAS_INT_MAX certificate as acquired arrays. */
        py_ref zeros(npy_intp n) const noexcept { return this->template zeros_as<T>(n); }

        /** @brief Two-dimensional variant of zeros(), Fortran-ordered, same certificate. */
        py_ref zeros(npy_intp m, npy_intp n) const noexcept { return this->template zeros_as<T>(m, n); }
    };


    /**
     * @brief Validate a call's positional/keyword structure, reproducing CPython's checks,
     *        order and wording.
     *
     * Checks run as four separate passes, not one interleaved loop, so an earlier fault does not
     * pre-empt a later one: too many arguments, then missing required, then given by both name
     * and position, then unexpected keyword.
     *
     * Takes no scalar type, so it is a free function rather than a `Ctx` member; as a member it
     * would be emitted once per flavor.
     *
     * @param args  Positional-argument tuple, or nullptr.
     * @param kwds  Keyword-argument dict, or nullptr.
     * @param fn    Module-qualified routine name for messages ("_flapack.dgesvx").
     * @param kw    Null-terminated argument names, in signature order.
     * @param nreq  Number of leading required arguments, those before the format's `'|'`.
     * @return      true when the call is well-formed, false with a TypeError set otherwise.
     */
    inline bool parse_args(PyObject *args, PyObject *kwds, const char *fn, const char *const *kw, int nreq) noexcept
    {
        Py_ssize_t npos = args ? PyTuple_GET_SIZE(args) : 0;
        Py_ssize_t nkw = kwds ? PyDict_Size(kwds) : 0;
        int nnames = 0;
        while (kw[nnames]) { nnames++; }

        /** CPython's PyArg validates in *separate passes*, not one interleaved loop -- so a
         * clash at an earlier index does not pre-empt a missing-required at a later one
         * (`ddot(x, x=x)` reports missing 'y', not the 'x' clash).
         */

        /** 1. too many arguments: the count is positional + keyword *total*, so an extra,
         * duplicate, or unexpected keyword causes "takes at most 2 (3 given)". */
        if (npos + nkw > nnames) {
            PyErr_Format(PyExc_TypeError, "%s() takes at most %d argument%s (%zd given)", fn, nnames, nnames == 1 ? "" : "s", npos + nkw);
            return false;
        }

        /* 2. missing required arguments (signature order) */
        for (int i = 0; i < nreq; i++) {
            bool given = (i < npos) || (kwds && PyDict_GetItemString(kwds, kw[i]) != nullptr);
            if (!given) {PyErr_Format(PyExc_TypeError, "%s() missing required argument '%s' (pos %d)", fn, kw[i], i + 1);
                return false;
            }
        }

        /* 3. arguments given by both name and position (signature order) */
        for (int i = 0; i < nnames && i < npos; i++) {
            if (kwds && PyDict_GetItemString(kwds, kw[i]) != nullptr) {
                PyErr_Format(PyExc_TypeError, "argument for %s() given by name ('%s') and position (%d)", fn, kw[i], i + 1);
                return false;
            }
        }

        /** 4. unexpected keyword arguments.  The kwlist scan is `Ctx::index` spelled out,
         * so this stays free of the Ctx template. */
        if (kwds) {
            PyObject *key, *val;
            Py_ssize_t pos = 0;
            while (PyDict_Next(kwds, &pos, &key, &val)) {
                const char *k = PyUnicode_AsUTF8(key);
                if (k == nullptr) { return false; }
                int i = 0;
                while (kw[i] && strcmp(kw[i], k) != 0) { i++; }
                if (kw[i] == nullptr) {
                    PyErr_Format(PyExc_TypeError, "%s() got an unexpected keyword argument '%s'", fn, k);
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Looks each argument's raw Python object up by name, whether it arrived positionally
     *        or as a keyword.
     *
     * Holds `args`/`kwds` and resolves a name against the kwlist, so a wrapper spells each
     * argument once, in the kwlist. `parse()` validates the call's structure; the SCALAR_ and
     * ARRAY_ macros then convert and acquire, reading through it.
     *
     * Binds `CtxBase<Mod>` rather than `Ctx<Mod, T>`, so one instantiation serves every flavor.
     * An accompanying deduction guide maps a derived `Ctx<Mod, T>` argument back to
     * `ArgInspector<Mod>`; without it CTAD would deduce on the flavor and reinstate a copy per
     * flavor.
     *
     * @tparam Mod  Module traits supplying `pyname` and `libname`.
     */
    template <class Mod>
    class ArgInspector {
    public:
        ArgInspector(PyObject *args, PyObject *kwds, const CtxBase<Mod> &ctx) : args_(args), kwds_(kwds), ctx_(ctx) {}

        /** @brief CPython-order structural validation; false with an exception set on failure. */
        bool parse() const noexcept
        {
            return parse_args(args_, kwds_, ctx_.qualname(), ctx_.kws(), ctx_.nreq());
        }

        /**
         * @brief Borrowed object supplied for @p name, from its positional slot or the keywords.
         *
         * Call only after `parse()` has succeeded, which is what guarantees no name is filled
         * both positionally and by keyword.
         *
         * @param name  Argument name; must appear in the kwlist.
         * @return      Borrowed reference, or nullptr when the argument was omitted (the
         *              caller then applies the optional's default). Also nullptr with a
         *              SystemError set when @p name is not in the kwlist -- a macro or kwlist
         *              authoring fault rather than caller input, so it is raised in every build
         *              rather than left to an assert.
         */
        PyObject *raw(const char *name) const noexcept
        {
            int idx = ctx_.index(name);
            if (idx < 0) {
                PyErr_Format(PyExc_SystemError, "%s(): argument '%s' is missing from the kwlist", ctx_.qualname(), name);
                return nullptr;
            }
            return at(idx, name);
        }

        /**
         * @brief `raw()` for a name that need not be in the kwlist.
         *
         * Used for outputs a routine produces but does not accept -- `geqp3` returns `jpvt`,
         * `tau` and `work` without taking any of them -- where absence is the normal case and
         * not an authoring fault.
         *
         * @param name  Argument name; absence from the kwlist is not an error.
         * @return      Borrowed reference, or nullptr when @p name is absent from the kwlist or
         *              was omitted from the call.
         */
        PyObject *raw_opt(const char *name) const noexcept
        {
            int idx = ctx_.index(name);
            if (idx < 0) { return nullptr; }
            PyObject *o = at(idx, name);
            /* An explicit `None` means "not supplied", as it did for f2py: callers pass it to
             * select the omitted branch positionally, e.g. `gtsvx(..., dlf=None)` when `fact`
             * is 'N'.  Only the optional lookup does this; a required argument given as None
             * is a conversion error, not an omission. */
            return o == Py_None ? nullptr : o;
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
        const CtxBase<Mod> &ctx_;
    };

    template <class Mod, class T>
    ArgInspector(PyObject *, PyObject *, const Ctx<Mod, T> &) -> ArgInspector<Mod>;

}  // namespace wrapper


/* ====================================================================================
 * Argument-handling vocabulary
 *
 * A wrapper calls PARSE_ARGS(), then declares each argument with one macro line, in processing
 * order.  Two names are conventional throughout: `P`, the ArgInspector that supplies an
 * argument's raw Python object by name, and `ctx`, the per-routine context that formats
 * messages.
 *
 * Argument-declaration macros are named `<KIND>_<VARIANT>`:
 *
 *     SCALAR_REQ      required scalar
 *     SCALAR_OPT      optional scalar, with a default
 *     SCALAR_FLAG     optional `overwrite_*` flag, strict `__index__` conversion
 *     ARRAY_IN        read-only, supplied by the caller
 *     ARRAY_INOUT     read-write, supplied by the caller and returned
 *     ARRAY_OUT       result returned to the caller
 *     ARRAY_HIDDEN    scratch, never visible from Python
 *
 * Every other macro is a verb: PARSE_ARGS, CHECK, CHECKARRAY, RETURN, and LAPACK's
 * CALLABLE_SELECT.
 *
 * The array macros classify by destination, not by how the array is obtained: IN and INOUT
 * arrive from the caller, OUT is returned, HIDDEN is kept internal.  Whether a given ARRAY_OUT
 * accepts a caller-supplied buffer does not change its kind.  All four take the element type as
 * their leading argument, because LAPACK passes integer and real-counterpart arrays alongside a
 * routine's own flavor.
 *
 * Inside a wrapper body:
 *   - a bare array name (`x`, `y`, `a`, `b`) is the acquired array, owned by a `py_ref` that
 *     releases on every exit path;
 *   - a bare lower-case scalar name (`n`, `incx`, `alpha`) is the converted C scalar;
 *   - a check is an ordinary C expression, stringized into its own message; `len()`, `shape()`
 *     and `abs()` are functions, so checks compile as written.
 * ==================================================================================== */

/** @brief Declare `ArgInspector P` over `args`/`kwds` and run CPython-order structural parsing. */
#define PARSE_ARGS() \
    wrapper::ArgInspector P(args, kwds, ctx); \
    if (!P.parse()) { return nullptr; }

/**
 * @brief Acquire a read-only array argument (`intent(in)`), declaring the owning `py_ref <name>`.
 *
 * The caller's array is used as-is when it already satisfies dtype, rank, Fortran-contiguity and
 * alignment; otherwise a converted copy is used and the caller's array is left untouched.
 *
 * Expands to two statements, a declaration and a failure return, so it is valid at block level
 * only.
 *
 * @param type  Element type.  Given explicitly rather than taken from the wrapper's flavor,
 *              because a routine may take an integer array beside its float ones (`getrs`'s
 *              `piv` next to its `lu`).
 * @param name  Argument name; also the declared `py_ref` and the name used in messages.
 * @param ndim  Required rank; other ranks are reinterpreted as f2py did.
 */
#define ARRAY_IN(type, name, ndim) \
    py_ref name = ctx.template in<type>(P.raw(#name), ndim, #name); \
    if (!name) { return nullptr; }

/**
 * @brief Acquire a read-write array argument (`intent(in,out)`), declaring the owning
 *        `py_ref <name>`.
 *
 * The caller's buffer is written in place when @p overwrite is true and it already satisfies
 * dtype, rank, Fortran-contiguity and alignment; otherwise a writable copy is made.  Either way
 * the declared `py_ref` holds what the wrapper returns.
 *
 * Expands to two statements, a declaration and a failure return, so it is valid at block level
 * only.
 *
 * @param type       Element type, as for ARRAY_IN.
 * @param name       Argument name; also the declared `py_ref` and the name used in messages.
 * @param ndim       Required rank; other ranks are reinterpreted as f2py did.
 * @param overwrite  True permits writing through the caller's buffer.
 */
#define ARRAY_INOUT(type, name, ndim, overwrite) \
    py_ref name = ctx.template inout<type>(P.raw(#name), ndim, overwrite, #name); \
    if (!name) { return nullptr; }

/**
 * @brief Acquire a result the wrapper returns, declaring the owning `py_ref <name>`.
 *
 * When the caller supplied the array it is taken read-write and returned; when it is absent
 * @p def allocates it.  Many results cannot be supplied at all -- their names are not in the
 * kwlist -- in which case only the @p def branch is reachable.
 *
 * Use ARRAY_HIDDEN instead for a buffer that never reaches the caller.  The distinction is what
 * the array does, not whether the caller may supply one: a `work` buffer whose optimal size the
 * routine reports back (`gees`, `gges`, `gelss`, `geqp3`) is an ARRAY_OUT.
 *
 * @param type       Element type.  Needed by both branches: @p def allocates with it and a
 *                   supplied array is converted to it.  It is not always the wrapper's flavor --
 *                   `gesvx` returns `ipiv` (`CBLAS_INT`) and the scalings `r`, `c` (the real
 *                   counterpart) beside a flavor-typed `af`.
 * @param name       Result name; also the declared `py_ref` and the name used in messages.
 * @param ndim       Required rank when the caller supplies the array.
 * @param overwrite  True permits writing through a supplied array instead of copying it.
 * @param def        Allocation used when the array is absent: `ctx.zeros(len)` for a vector,
 *                   `ctx.zeros(m, n)` for a matrix, `ctx.zeros_as<V>(...)` for an element type
 *                   other than the wrapper's flavor.  Zero-filled, because a routine may
 *                   accumulate into the buffer rather than overwrite it (gemv's `beta * y`).
 */
#define ARRAY_OUT(type, name, ndim, overwrite, def) \
    PyObject *name##_raw = P.raw_opt(#name); \
    py_ref name = (name##_raw == nullptr) ? (def) \
                                          : ctx.template inout<type>(name##_raw, ndim, overwrite, #name); \
    if (!name) { return nullptr; }

/**
 * @brief Declare a zero-filled scratch buffer that never reaches the caller -- f2py's
 *        `intent(hide)` work arrays (`work`, `rwork`, `iwork`, `bwork`).
 *
 * Use ARRAY_OUT instead for any buffer the wrapper returns, including one the caller cannot
 * supply.
 *
 * @param type  Element type.  Given explicitly because these are regularly the real counterpart
 *              of a complex flavor (`rwork`) or the ABI-width integer backing Fortran `LOGICAL`
 *              (`bwork`).
 * @param name  Declared `py_ref`.
 * @param ...   Dimensions: one argument for a vector, two for a matrix.
 */
#define ARRAY_HIDDEN(type, name, ...) \
    py_ref name = ctx.template zeros_as<type>(__VA_ARGS__); \
    if (!name) { return nullptr; }

/**
 * @brief Declare an optional scalar `<name>` of C type @p type, converted from the caller's
 *        object or set from @p def when the argument is absent.
 *
 * Conversion is f2py's permissive `from_pyobj`: a float truncates to int, a complex contributes
 * its real part, a sequence its element 0.  On failure it raises TypeError and returns nullptr
 * from the wrapper.
 *
 * @param type  C type of the declared variable.
 * @param name  Argument name; also the declared variable.
 * @param def   Default expression, cast to @p type, used when the argument is absent.
 */
#define SCALAR_OPT(type, name, def) \
    type name; \
    do { \
        PyObject *name##_raw = P.raw(#name); \
        if (name##_raw == nullptr) { name = static_cast<type>(def); } \
        else if (!ctx.scalar(&name, name##_raw, #name)) { return nullptr; } \
    } while (0)

/**
 * @brief Declare a required scalar `<name>` of C type @p type, converted from the caller's
 *        object.
 *
 * `parse()` has already established that the argument is present, so the object is converted
 * unconditionally, by the same permissive rules as SCALAR_OPT.
 *
 * @param type  C type of the declared variable.
 * @param name  Argument name; also the declared variable.
 */
#define SCALAR_REQ(type, name) \
    type name; \
    do { \
        if (!ctx.scalar(&name, P.raw(#name), #name)) { return nullptr; } \
    } while (0)

/**
 * @brief Declare an `int` flag `<name>`, defaulting to 0 -- the `overwrite_*` controls, which
 *        are not Fortran arguments but decide whether the wrapper may write in place.
 *
 * Conversion follows Python's `__index__` protocol (`PyNumber_Index`), which accepts only
 * objects that are genuinely integers -- `int`, `bool`, numpy integer scalars.  A float, str or
 * None is rejected with a TypeError, where SCALAR_OPT's permissive path would have coerced it
 * and, for a flag, silently changed the meaning of the call.
 *
 * Place immediately after PARSE_ARGS(), ahead of the body scalars, so that a bad flag is
 * reported before a bad scalar; the resulting error precedence is fixed and tested.
 *
 * @param name  Argument name; also the declared variable.
 */
#define SCALAR_FLAG(name) \
    int name = 0; \
    do { \
        PyObject *name##_raw = P.raw(#name); \
        if (name##_raw != nullptr && !wrapper::from_pyobj_index(&name, name##_raw)) { return nullptr; } \
    } while (0)

/**
 * @brief Validate a converted scalar against a C expression, raising ValueError if it is false.
 *
 * The expression is stringized into the message, so `CHECK(incx != 0, incx)` produces
 * `(incx != 0) failed for 4th keyword incx: daxpy:incx=0`.
 *
 * @param expr  C expression that must hold.
 * @param name  Argument the failure is reported against, and whose value is printed.  It need
 *              not appear in @p expr -- bound checks on a dimension report `n` -- and it need
 *              not be an argument at all: a hidden dimension such as `tbtrs`'s `ldb` has no
 *              position to name, so it reports as `(ldb >= n) failed for ldb: dtbtrs:ldb=3`.
 */
#define CHECK(expr, name) \
    do { if (!ctx.check_scalar((expr), #expr, #name, name)) { return nullptr; } } while (0)

/**
 * @brief Validate an acquired array against a C expression, raising ValueError if it is false.
 *
 * The array counterpart of CHECK, and unlike CHECK it prints no value:
 * `CHECKARRAY(shape(a, 0) == shape(a, 1), a)` produces
 * `(shape(a, 0) == shape(a, 1)) failed for 2nd argument a`.
 *
 * @param expr  C expression that must hold.
 * @param name  Acquired array (a `py_ref` local) the failure is reported against.
 */
#define CHECKARRAY(expr, name) do { if (!ctx.check_array((expr), #expr, #name)) { return nullptr; } } while (0)

/**
 * @brief Pack the wrapper's results and return them: `RETURN(y)`, `RETURN(x, y)`.
 *
 * One argument returns that value alone, several return a tuple in the order given.  Compute
 * values into locals first; expressions do not belong here.
 *
 * @param ...  Results in declaration order -- `py_ref` arrays, scalars, or `long long`.
 */
#define RETURN(...) return wrapper::make_result(__VA_ARGS__)

/**
 * @brief Emit the four method-table rows of one routine -- `s<name>`, `d<name>`, `c<name>`,
 *        `z<name>` -- each bound to the wrapper template at the matching scalar type.
 *
 * Expand inside `blas::capi` or `lapack::capi`, where the wrapper templates and the width
 * aliases resolve unqualified.
 *
 * @param name  Wrapper template name, which is also the unprefixed routine name.
 */
#define FAMILY(name) \
    {"s" #name, (PyCFunction)(void (*)())name<f32>,  METH_VARARGS | METH_KEYWORDS, nullptr}, \
    {"d" #name, (PyCFunction)(void (*)())name<f64>,  METH_VARARGS | METH_KEYWORDS, nullptr}, \
    {"c" #name, (PyCFunction)(void (*)())name<c64>,  METH_VARARGS | METH_KEYWORDS, nullptr}, \
    {"z" #name, (PyCFunction)(void (*)())name<c128>, METH_VARARGS | METH_KEYWORDS, nullptr}

/**
 * @brief One explicitly named method-table row, for families whose Python name is not the
 *        flavor letter followed by the template name (`scnrm2` binds `nrm2<c64>`).
 *
 * @p pyname must match the prefix trait passed to the wrapper's `Ctx`, since that is the name
 * its error messages carry.
 *
 * @param pyname  Python-visible routine name.
 * @param name    Wrapper template name.
 * @param T       Scalar type to instantiate it at.
 */
#define ROW(pyname, name, T) \
    {#pyname, (PyCFunction)(void (*)())name<T>, METH_VARARGS | METH_KEYWORDS, nullptr}

/**
 * @brief ROW for a wrapper template taking two type parameters, `name<T, A>`.
 *
 * Used by the real-scalar `scal` overload, `csscal` = `scal<c64, f32>`, where the second
 * parameter is the scalar type and also what distinguishes it from the regular `scal<T>`.
 *
 * @param pyname  Python-visible routine name.
 * @param name    Wrapper template name.
 * @param T       First type parameter, the array element type.
 * @param A       Second type parameter, the scalar type.
 */
#define ROW2(pyname, name, T, A) \
    {#pyname, (PyCFunction)(void (*)())name<T, A>, METH_VARARGS | METH_KEYWORDS, nullptr}
