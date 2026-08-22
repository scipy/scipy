/**
 * @file
 * @brief Bridge between LAPACK's plain-function-pointer callback ABI (`gees`'s `SELECT`, `gges`'s
 *        `SELCTG`) and a Python callable.
 *
 * A wrapper installs a `CallbackFrame` for the duration of the Fortran call. The trampoline
 * finds it in thread-local storage, which the callback ABI forces since it carries no user-data
 * slot, and calls the Python callable with Python scalars. A failing callback aborts the Fortran
 * routine through `setjmp`/`longjmp` rather than letting it run to completion.
 *
 * The file holds the mechanism, then the trampolines, then the per-routine traits selecting one,
 * then the `CALLABLE_*` macros wrappers write.
 *
 * @note Wrappers must not release the interpreter around the Fortran call (no
 *       `Py_BEGIN_ALLOW_THREADS`) in any build: the trampoline re-enters Python synchronously on
 *       that same OS thread.
 */
#pragma once

#include <csetjmp>

#include "wrapper_helpers.hpp"   /* Python.h, CBLAS_INT, to_pyobj and the flavor aliases */

namespace lapack {

    /* The shared width aliases, used by the trampoline signatures below.  `lapack_calls.hpp` and
     * `lapack_helpers.hpp` re-export the same names, which is redundant but legal and lets any of
     * them be included on its own. */
    using wrapper::f32;
    using wrapper::f64;
    using wrapper::c64;
    using wrapper::c128;

    /**
     * @brief Per-call state for the one user callback active on this thread.
     *
     * One frame type serves every callback-taking routine; the trampoline decides how many scalars
     * to pass. Nested calls install nested frames, so a failure always returns to the innermost
     * Fortran call.
     */
    struct CallbackFrame {
        PyObject *callable;             /**< Borrowed; owned by the wrapper's arguments. */
        PyObject *extra_args;           /**< Borrowed tuple, appended after the LAPACK scalars. */
        Py_ssize_t argcount;            /**< Positional parameters the callable declares, or -1
                                             when it cannot be introspected (a builtin). */
        Py_ssize_t ndefaults;           /**< How many of those have defaults. */
        /* Zeroed here so CALLABLE_SELECT can brace-initialize the frame from the two members
         * it actually knows; `setjmp` fills this in immediately afterwards. */
        std::jmp_buf jmpbuf{};          /**< Set by setjmp in the wrapper; target of the abort jump. */
    };

    namespace detail { inline thread_local CallbackFrame *active_frame = nullptr; }

    /**
     * @brief Read how many positional arguments @p fn declares, and how many of them are optional.
     *
     * LAPACK offers a fixed number of scalars -- two for real `gees`, one for complex, three for
     * real `gges` -- but a caller may supply a callable that takes fewer, and f2py passed only as
     * many as it accepted.  `scipy.linalg.schur(sort=lambda x: x >= 0)` relies on that.
     *
     * @param fn         The callable.
     * @param argcount   Receives its positional count, or -1 when it declares none observably
     *                   (a builtin), meaning "pass whatever LAPACK offers".
     * @param ndefaults  Receives how many trailing parameters have defaults.
     */
    inline void callable_arity(PyObject *fn, Py_ssize_t *argcount, Py_ssize_t *ndefaults) noexcept
    {
        *argcount = -1;
        *ndefaults = 0;

        /* A bound method's `__code__` counts `self`, which the caller never supplies. */
        Py_ssize_t self_arg = PyFunction_Check(fn) ? 0 : 1;

        PyObject *code = PyObject_GetAttrString(fn, "__code__");
        if (code != nullptr) {
            PyObject *n = PyObject_GetAttrString(code, "co_argcount");
            Py_DECREF(code);
            if (n != nullptr) {
                Py_ssize_t v = PyLong_AsSsize_t(n);
                Py_DECREF(n);
                if (!(v == -1 && PyErr_Occurred())) { *argcount = v - self_arg; }
            }
        }
        PyObject *defaults = PyObject_GetAttrString(fn, "__defaults__");
        if (defaults != nullptr) {
            if (PyTuple_Check(defaults)) { *ndefaults = PyTuple_GET_SIZE(defaults); }
            Py_DECREF(defaults);
        }
        PyErr_Clear();   /* a missing attribute is not an error, it just means "unknown" */
    }

    /**
     * @brief Reject a callback that cannot be satisfied by @p maxnofargs scalars, before the
     *        Fortran call starts.
     *
     * f2py validated at call setup, so an unusable callback was refused even when the routine
     * would never have invoked it -- `sort_t=0` for `gees`.  Checking here rather than in the
     * trampoline keeps that, and keeps the failure on the calling thread with a plain `return`
     * instead of inside a callback that would have to `longjmp` out of Fortran.
     *
     * @return false with a TypeError set when the callable requires more positional arguments
     *         than are on offer.
     */
    inline bool check_callable_arity(Py_ssize_t argcount, Py_ssize_t ndefaults,
                                     Py_ssize_t maxnofargs, Py_ssize_t nextra,
                                     const char *kwname) noexcept
    {
        if (argcount < 0) { return true; }        /* not introspectable; f2py passed them all */
        Py_ssize_t offered = maxnofargs + nextra;
        Py_ssize_t required = argcount - ndefaults;
        if (required > offered) {
            PyErr_Format(PyExc_TypeError,
                         "%s requires %zd positional argument%s but at most %zd can be supplied",
                         kwname, required, required == 1 ? "" : "s", offered);
            return false;
        }
        return true;
    }

    /**
     * @brief RAII guard installing @p frame as this thread's active frame, restoring the previous
     *        one on scope exit.
     *
     * Must be constructed before the matching `setjmp`, so that the jump lands in a frame where
     * this object is already alive and is destroyed by the ordinary `return`.
     */
    class ScopedFrame {
    public:
        explicit ScopedFrame(CallbackFrame *frame) noexcept : prev_(detail::active_frame) { detail::active_frame = frame; }
        ~ScopedFrame() { detail::active_frame = prev_; }
        ScopedFrame(const ScopedFrame &) = delete;
        ScopedFrame &operator=(const ScopedFrame &) = delete;
    private:
        CallbackFrame *prev_;
    };

    /**
     * @brief Release callback arguments and jump back to the `setjmp` for @p f.
     *
     * @param argv Arguments to release; may be nullptr with @p argc 0 if already released.
     *
     * The Python error indicator must already be set, and is left set across the jump so the
     * wrapper can return nullptr directly. Keeping the exception out of an automatic object also
     * avoids reading a value modified after `setjmp`.
     */
    [[noreturn]] inline void discard_and_abort(CallbackFrame *f, PyObject *const *argv, Py_ssize_t argc) noexcept
    {
        for (Py_ssize_t i = 0; i < argc; i++) { Py_XDECREF(argv[i]); }
        std::longjmp(f->jmpbuf, 1);
    }

    /**
     * @brief Call the active frame's Python callable with @p argv and return the LAPACK LOGICAL
     *        value Fortran expects (0 or 1).
     *
     * @param argv LAPACK scalars to pass. Ownership transfers here: each entry is released exactly
     *             once on every exit path, so trampolines hand over fresh references and never
     *             decref them.
     *
     * Any failure -- a null scalar, a raising callable, a result whose `__bool__` raises, or the
     * recursion guard tripping -- leaves the exception set and jumps back to the wrapper, so LAPACK
     * stops instead of continuing with a synthetic "do not select".
     */
    inline CBLAS_INT invoke_or_abort(PyObject *const *argv, Py_ssize_t argc) noexcept
    {
        CallbackFrame *f = detail::active_frame;
        if (f == nullptr) {
            for (Py_ssize_t i = 0; i < argc; i++) { Py_XDECREF(argv[i]); }
            return 0;
        }

        for (Py_ssize_t i = 0; i < argc; i++) {
            if (argv[i] == nullptr) { discard_and_abort(f, argv, argc); }
        }

        Py_ssize_t nextra = f->extra_args == nullptr ? 0 : PyTuple_GET_SIZE(f->extra_args);
        if (nextra > PY_SSIZE_T_MAX - argc) {
            PyErr_NoMemory();
            discard_and_abort(f, argv, argc);
        }

        /* Pass only as many LAPACK scalars as the callable declares, f2py's rule: it offers
         * `argc + nextra` and the callable takes `argcount`, so the call gets the smaller of the
         * two.  A callable needing more than that was already refused by CALLABLE_SELECT. */
        Py_ssize_t declared = f->argcount < 0 ? argc + nextra : f->argcount;
        Py_ssize_t supplied = (argc + nextra) < declared ? (argc + nextra) : declared;
        Py_ssize_t nfixed = supplied - nextra < 0 ? 0 : supplied - nextra;
        /* The scalars beyond what is passed are still owned here and released with the rest. */
        for (Py_ssize_t i = nfixed; i < argc; i++) { Py_DECREF(argv[i]); }
        argc = nfixed;
        PyObject **call_argv = const_cast<PyObject **>(argv);
        PyObject **combined_argv = nullptr;
        if (nextra != 0) {
            combined_argv = static_cast<PyObject **>(PyMem_Malloc((argc + nextra) * sizeof(*combined_argv)));
            if (combined_argv == nullptr) {
                PyErr_NoMemory();
                discard_and_abort(f, argv, argc);
            }
            for (Py_ssize_t i = 0; i < argc; i++) { combined_argv[i] = argv[i]; }
            for (Py_ssize_t i = 0; i < nextra; i++) {
                combined_argv[argc + i] = PyTuple_GET_ITEM(f->extra_args, i);
            }
            call_argv = combined_argv;
        }

        /* Vectorcall leaves recursion accounting to its caller. The guard also bounds a callback
         * that recursively enters another callback-bearing LAPACK wrapper. */
        if (Py_EnterRecursiveCall(" while evaluating a LAPACK eigenvalue-sort callback")) {
            PyMem_Free(combined_argv);
            discard_and_abort(f, argv, argc);   /* RecursionError already set */
        }
        PyObject *result = PyObject_Vectorcall(f->callable, call_argv, (size_t)(argc + nextra), nullptr);
        Py_LeaveRecursiveCall();
        PyMem_Free(combined_argv);
        for (Py_ssize_t i = 0; i < argc; i++) { Py_DECREF(argv[i]); }   /* done with argv either way */
        if (result == nullptr) {
            discard_and_abort(f, nullptr, 0);
        }
        int truth = PyObject_IsTrue(result);   /* -1 if __bool__ raises */
        Py_DECREF(result);
        if (truth < 0) {
            discard_and_abort(f, nullptr, 0);
        }
        return truth;
    }

}  // namespace lapack


/**
 * @brief Define one SELECT trampoline in the 1-, 2- or 3-scalar shape its Fortran interface uses.
 *
 * `gees` passes a real eigenvalue as (real, imaginary) and a complex one whole; `gges` adds
 * `beta` to each. `inline` lets the definitions live in this header and still give Fortran one
 * address per trampoline.
 */
#define LAPACK_SELECT_TRAMPOLINE_1(fname, T) \
    extern "C" inline CBLAS_INT fname(T *a0) noexcept { \
        PyObject *argv[1] = { wrapper::to_pyobj(*a0) }; \
        return lapack::invoke_or_abort(argv, 1); \
    }

#define LAPACK_SELECT_TRAMPOLINE_2(fname, T) \
    extern "C" inline CBLAS_INT fname(T *a0, T *a1) noexcept { \
        PyObject *argv[2] = { wrapper::to_pyobj(*a0), wrapper::to_pyobj(*a1) }; \
        return lapack::invoke_or_abort(argv, 2); \
    }

#define LAPACK_SELECT_TRAMPOLINE_3(fname, T) \
    extern "C" inline CBLAS_INT fname(T *a0, T *a1, T *a2) noexcept { \
        PyObject *argv[3] = { wrapper::to_pyobj(*a0), wrapper::to_pyobj(*a1), wrapper::to_pyobj(*a2) }; \
        return lapack::invoke_or_abort(argv, 3); \
    }


namespace lapack {

    /* The trampolines themselves: `gees` is 2-scalar for a real flavor and 1-scalar for a complex
     * one, `gges` one wider in each case. */
    LAPACK_SELECT_TRAMPOLINE_2(scipy_lapack_sgees_select, f32)
    LAPACK_SELECT_TRAMPOLINE_2(scipy_lapack_dgees_select, f64)
    LAPACK_SELECT_TRAMPOLINE_1(scipy_lapack_cgees_select, c64)
    LAPACK_SELECT_TRAMPOLINE_1(scipy_lapack_zgees_select, c128)
    LAPACK_SELECT_TRAMPOLINE_3(scipy_lapack_sgges_select, f32)
    LAPACK_SELECT_TRAMPOLINE_3(scipy_lapack_dgges_select, f64)
    LAPACK_SELECT_TRAMPOLINE_2(scipy_lapack_cgges_select, c64)
    LAPACK_SELECT_TRAMPOLINE_2(scipy_lapack_zgges_select, c128)

    /**
     * @brief Per-flavor lookup for a routine's callback argument: its keyword names, the
     *        trampoline with the matching Fortran signature, and how many scalars that
     *        trampoline passes.
     *
     * The keyword names stay flavor-specific (`sselect`, `dselect`, ...) for compatibility with the
     * f2py-generated module; CALLABLE_SELECT reads them from here so no wrapper spells a flavor.
     */
    template <class T> struct gees_select_traits;
    template <> struct gees_select_traits<f32> {
        static constexpr const char *kwname = "sselect";
        static constexpr const char *extra_kwname = "sselect_extra_args";
        inline static constexpr auto fn = scipy_lapack_sgees_select;
        static constexpr Py_ssize_t maxnofargs = 2;
    };
    template <> struct gees_select_traits<f64> {
        static constexpr const char *kwname = "dselect";
        static constexpr const char *extra_kwname = "dselect_extra_args";
        inline static constexpr auto fn = scipy_lapack_dgees_select;
        static constexpr Py_ssize_t maxnofargs = 2;
    };
    template <> struct gees_select_traits<c64> {
        static constexpr const char *kwname = "cselect";
        static constexpr const char *extra_kwname = "cselect_extra_args";
        inline static constexpr auto fn = scipy_lapack_cgees_select;
        static constexpr Py_ssize_t maxnofargs = 1;
    };
    template <> struct gees_select_traits<c128> {
        static constexpr const char *kwname = "zselect";
        static constexpr const char *extra_kwname = "zselect_extra_args";
        inline static constexpr auto fn = scipy_lapack_zgees_select;
        static constexpr Py_ssize_t maxnofargs = 1;
    };

    template <class T> struct gges_select_traits;
    template <> struct gges_select_traits<f32> {
        static constexpr const char *kwname = "sselect";
        static constexpr const char *extra_kwname = "sselect_extra_args";
        inline static constexpr auto fn = scipy_lapack_sgges_select;
        static constexpr Py_ssize_t maxnofargs = 3;
    };
    template <> struct gges_select_traits<f64> {
        static constexpr const char *kwname = "dselect";
        static constexpr const char *extra_kwname = "dselect_extra_args";
        inline static constexpr auto fn = scipy_lapack_dgges_select;
        static constexpr Py_ssize_t maxnofargs = 3;
    };
    template <> struct gges_select_traits<c64> {
        static constexpr const char *kwname = "cselect";
        static constexpr const char *extra_kwname = "cselect_extra_args";
        inline static constexpr auto fn = scipy_lapack_cgges_select;
        static constexpr Py_ssize_t maxnofargs = 2;
    };
    template <> struct gges_select_traits<c128> {
        static constexpr const char *kwname = "zselect";
        static constexpr const char *extra_kwname = "zselect_extra_args";
        inline static constexpr auto fn = scipy_lapack_zgges_select;
        static constexpr Py_ssize_t maxnofargs = 2;
    };

}  // namespace lapack


/**
 * @brief Acquire a routine's Python callback argument and produce the function pointer to pass
 *        straight into the Fortran call -- the callback analog of ARRAY_IN.
 *
 *     CALLABLE_SELECT(gees, select);
 *     ...
 *     CALLABLE_CALL(select, lapack::gees(..., select, ...));
 *
 * Declares `<name>` (the trampoline for `T`), `<name>_obj`, `<name>_extra_args` and `<name>_frame`.
 * Requires `<routine>_select_traits<T>`, declared above.  The optional `*_extra_args` value must
 * be a tuple; it is appended to every callback call.
 */
#define CALLABLE_SELECT(routine, name) \
    PyObject *name##_obj = P.raw(lapack::routine##_select_traits<T>::kwname); \
    if (name##_obj == nullptr || !PyCallable_Check(name##_obj)) { \
        PyErr_Format(PyExc_TypeError, "%s must be callable", lapack::routine##_select_traits<T>::kwname); \
        return nullptr; \
    } \
    PyObject *name##_extra_args = P.raw(lapack::routine##_select_traits<T>::extra_kwname); \
    if (name##_extra_args != nullptr && !PyTuple_Check(name##_extra_args)) { \
        PyErr_Format(PyExc_TypeError, "%s must be a tuple", lapack::routine##_select_traits<T>::extra_kwname); \
        return nullptr; \
    } \
    Py_ssize_t name##_argcount, name##_ndefaults; \
    lapack::callable_arity(name##_obj, &name##_argcount, &name##_ndefaults); \
    if (!lapack::check_callable_arity(name##_argcount, name##_ndefaults, \
                                      lapack::routine##_select_traits<T>::maxnofargs, \
                                      name##_extra_args == nullptr ? 0 : PyTuple_GET_SIZE(name##_extra_args), \
                                      lapack::routine##_select_traits<T>::kwname)) { \
        return nullptr; \
    } \
    lapack::CallbackFrame name##_frame{ name##_obj, name##_extra_args, name##_argcount, name##_ndefaults }; \
    auto name = lapack::routine##_select_traits<T>::fn

/**
 * @brief Run one Fortran call with the callback declared by CALLABLE_SELECT installed.
 *
 * Place it where the call belongs, after all setup: a callback failure returns nullptr with its
 * exception already set. The call expression may contain ordinary commas.
 *
 * @note Anything the Fortran call writes to a local (`sdim`, `info`) is indeterminate on the
 *       failure path, per the usual `setjmp` rule; nothing may read those before checking.
 */
#define CALLABLE_CALL(name, call) \
    do { \
        lapack::ScopedFrame name##_scope(&name##_frame); \
        if (setjmp(name##_frame.jmpbuf) == 0) { call; } \
        else { return nullptr; } \
    } while (0)
