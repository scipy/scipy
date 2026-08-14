/**
 * @file
 * @brief What `namespace lapack` adds to the shared wrapper machinery.
 *
 * The Python/numpy <-> C++ translation layer itself lives in `wrapper_helpers.hpp` and knows
 * nothing about LAPACK.  This header binds it to this extension -- the module identity its
 * error messages are raised under -- and adds what is LAPACK's alone: workspace sizing and the
 * `intent(hide)` work-buffer declaration.
 *
 * Everything the machinery provides is re-exported into `namespace lapack` below, so a wrapper
 * in `lapack::capi` writes `Ctx<T>`, `py_ref`, `shape(a, 1)` unqualified and never spells
 * `wrapper::`.
 */
#pragma once

#include "wrapper_helpers.hpp"

/**
 * The name these wrappers raise their errors under.  It is deliberately *not* the extension's
 * own name: the module is built as `_flapack_cpp` only while it is being compared routine by
 * routine against the f2py-generated `_flapack` (see `lapack_module.cpp`), and every message
 * has to match what f2py emitted from `python module _flapack` in `flapack.pyf.src`.  When the
 * C++ module takes over the `_flapack` name, this stops being a distinction.
 */
#ifdef HAVE_BLAS_ILP64
#define FLAPACK_ERRNAME "_flapack_64"
#else
#define FLAPACK_ERRNAME "_flapack"
#endif

namespace lapack {

    /** @brief This extension's identity for the shared machinery (see `wrapper_helpers.hpp`). */
    struct module {
        static constexpr const char *pyname  = FLAPACK_ERRNAME;
        static constexpr const char *libname = "LAPACK";
    };

    /** @brief The per-routine context, bound to this module: wrappers just write `Ctx<T>`. */
    template <class T> using Ctx = wrapper::Ctx<module, T>;

    /* The rest of the machinery, unchanged and unqualified for the wrappers.  `abs` in
     * particular must be found here: it is the reason unqualified `abs(n)` in a wrapper
     * resolves to the npy_intp overload and never to the C library's `::abs(int)`. */
    using wrapper::py_ref;
    using wrapper::len;
    using wrapper::shape;
    using wrapper::abs;
    using wrapper::to_pyobj;
    using wrapper::make_result;
    using wrapper::npy_type;
    using wrapper::f32;
    using wrapper::f64;
    using wrapper::c64;
    using wrapper::c128;
    using wrapper::real_of;
    using wrapper::real_of_t;
    using wrapper::is_complex;
    using wrapper::is_complex_v;
    using wrapper::flavor;

    /**
     * @brief Narrow a workspace size to the LAPACK integer type, applying LAPACK's `>= 1` floor.
     *
     * Write @p size the way the `.pyf` spells it, with the leading factor as `LL` so the
     * arithmetic happens in `long long`: `work_size(8LL * n + 16, &lwork)`.  f2py evaluated the
     * same expressions in C `int` and could wrap silently; this raises instead.  The widened
     * arithmetic cannot itself wrap, because `n` is an array dimension and the array machinery
     * has already certified it against `CBLAS_INT_MAX`.
     *
     * @param out  Receives the size on success; untouched on failure.
     * @return     true on success; false with an OverflowError set.
     */
    inline bool work_size(long long size, CBLAS_INT *out) noexcept
    {
        if (size > static_cast<long long>(CBLAS_INT_MAX)) {
            PyErr_SetString(PyExc_OverflowError, "LAPACK workspace size exceeds the LAPACK integer limit");
            return false;
        }
        *out = size > 1 ? static_cast<CBLAS_INT>(size) : 1;
        return true;
    }

}  // namespace lapack
