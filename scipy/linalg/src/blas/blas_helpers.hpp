/**
 * @file
 * @brief What `namespace blas` adds to the shared wrapper machinery.
 *
 * The Python/numpy <-> C++ translation layer itself lives in `wrapper_helpers.hpp` and knows
 * nothing about BLAS.  This header binds it to this extension -- the module identity its error
 * messages are raised under -- and adds the naming irregularities that are BLAS's alone.
 *
 * Everything the machinery provides is re-exported into `namespace blas` below, so a wrapper in
 * `blas::capi` writes `Ctx<T>`, `py_ref`, `len(x)` unqualified and never spells `wrapper::`.
 *
 * @note Every .cpp of this extension defines `PY_ARRAY_UNIQUE_SYMBOL` to the same name before
 *       including this header, and all but `blas_module.cpp` also define `NO_IMPORT_ARRAY`;
 *       see `wrapper_helpers.hpp` for why.
 */
#pragma once

#include "wrapper_helpers.hpp"

/**
 * The LP64 build generates the module `_fblas` and the ILP64 build `_fblas_64`, mirroring the
 * legacy f2py naming.  The ILP64 dependency passes `-DHAVE_BLAS_ILP64`, which already selects
 * `CBLAS_INT = int64_t` and the ILP64 `BLAS_FUNC` symbol suffix in `scipy_blas_defines.h`, so
 * the name is the only remaining difference.  Both spellings the build needs are written out
 * here: the string that names the module and appears in every error message (`blas::module`
 * below), and the `PyInit_` entry point (`blas_module.cpp`).
 */
#ifdef HAVE_BLAS_ILP64
#define FBLAS_MODULE_STR "_fblas_64"
#define FBLAS_PYINIT     PyInit__fblas_64
#else
#define FBLAS_MODULE_STR "_fblas"
#define FBLAS_PYINIT     PyInit__fblas
#endif

namespace blas {

    /**
     * @brief This extension's identity, as the shared machinery needs it (see `wrapper_helpers.hpp`).
     *
     * `pyname` matches the `python module` name of the .pyf these wrappers replaced, so the
     * messages are byte-identical to the ones f2py used to raise -- including in the ILP64
     * build, whose f2py predecessor said `_fblas_64`.
     */
    struct module {
        static constexpr const char *pyname  = FBLAS_MODULE_STR;
        static constexpr const char *libname = "BLAS";
    };

    /** @brief The per-routine context, bound to this module: wrappers just write `Ctx<T>`. */
    template <class T> using Ctx = wrapper::Ctx<module, T>;

    /* The rest of the machinery, unchanged and unqualified for the wrappers.  `abs` in
     * particular must be found here: it is the reason unqualified `abs(incx)` in a wrapper
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
     * @brief Irregular prefix, flavor letter first (s/d/cs/zd) -- the pyf `<tchar>` shorthand.
     *
     * Used by the routines whose *data* is complex but one *scalar argument* stays real:
     * `csrot`/`zdrot` (the cosine is real), `csscal`/`zdscal` (real scale factor).
     */
    template <class T> constexpr const char *tchar();
    template <> inline constexpr const char *tchar<f32>()  { return "s"; }
    template <> inline constexpr const char *tchar<f64>()  { return "d"; }
    template <> inline constexpr const char *tchar<c64>()  { return "cs"; }
    template <> inline constexpr const char *tchar<c128>() { return "zd"; }

    /**
     * @brief Irregular prefix, result letter first (s/d/sc/dz) -- pyf `<prefix3>`/`<prefix4>`.
     *
     * Used by the value-returning routines whose *result* is real even for complex data:
     * `scnrm2`/`dznrm2`, `scasum`/`dzasum`.  Note the opposite letter order from tchar().
     */
    template <class T> constexpr const char *tchar_fn();
    template <> inline constexpr const char *tchar_fn<f32>()  { return "s"; }
    template <> inline constexpr const char *tchar_fn<f64>()  { return "d"; }
    template <> inline constexpr const char *tchar_fn<c64>()  { return "sc"; }
    template <> inline constexpr const char *tchar_fn<c128>() { return "dz"; }

    /**
     * @brief Index-function prefix (is/id/ic/iz) -- the pyf `i<prefix>` pattern, `isamax`.
     *
     * The third naming irregularity: the `i` for "index-returning" precedes the flavor letter.
     */
    template <class T> constexpr const char *iflavor();
    template <> inline constexpr const char *iflavor<f32>()  { return "is"; }
    template <> inline constexpr const char *iflavor<f64>()  { return "id"; }
    template <> inline constexpr const char *iflavor<c64>()  { return "ic"; }
    template <> inline constexpr const char *iflavor<c128>() { return "iz"; }

}  // namespace blas
