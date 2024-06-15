#include <Python.h>
#include <numpy/npy_math.h>

#include <stdarg.h>
#include <stdlib.h>

#include "sf_error.h"

#include "special/error.h"
// #include "sf_error.h"

const char *sf_error_messages[] = {
    "no error",
    "singularity",
    "underflow",
    "overflow",
    "too slow convergence",
    "loss of precision",
    "no result obtained",
    "domain error",
    "invalid input argument",
    "other error",
    NULL
};

extern "C" int wrap_PyUFunc_getfperr(void);

void sf_error_v(const char *func_name, sf_error_t code, const char *fmt, va_list ap) {
    /* Internal function which takes a va_list instead of variadic args.
     * Makes this easier to wrap in error handling used in special C++
     * namespace for special function kernels provided by SciPy. */
    PyGILState_STATE save;
    PyObject *scipy_special = NULL;
    char msg[2048], info[1024];
    static PyObject *py_SpecialFunctionWarning = NULL;
    sf_action_t action;

    if ((int) code < 0 || (int) code >= 10) {
        code = SF_ERROR_OTHER;
    }
    action = scipy_sf_error_get_action(code);
    if (action == SF_ERROR_IGNORE) {
        return;
    }

    if (func_name == NULL) {
        func_name = "?";
    }

    if (fmt != NULL && fmt[0] != '\0') {
        PyOS_vsnprintf(info, 1024, fmt, ap);
        PyOS_snprintf(msg, 2048, "scipy.special/%s: (%s) %s", func_name, sf_error_messages[(int) code], info);
    } else {
        PyOS_snprintf(msg, 2048, "scipy.special/%s: %s", func_name, sf_error_messages[(int) code]);
    }

#ifdef WITH_THREAD
    save = PyGILState_Ensure();
#endif

    if (PyErr_Occurred()) {
        goto skip_warn;
    }

    scipy_special = PyImport_ImportModule("scipy.special");
    if (scipy_special == NULL) {
        PyErr_Clear();
        goto skip_warn;
    }

    if (action == SF_ERROR_WARN) {
        py_SpecialFunctionWarning = PyObject_GetAttrString(scipy_special, "SpecialFunctionWarning");
    } else if (action == SF_ERROR_RAISE) {
        py_SpecialFunctionWarning = PyObject_GetAttrString(scipy_special, "SpecialFunctionError");
    } else {
        /* Sentinel, should never get here */
        py_SpecialFunctionWarning = NULL;
    }
    /* Done with scipy_special */
    Py_DECREF(scipy_special);

    if (py_SpecialFunctionWarning == NULL) {
        PyErr_Clear();
        goto skip_warn;
    }

    if (action == SF_ERROR_WARN) {
        PyErr_WarnEx(py_SpecialFunctionWarning, msg, 1);
        /*
         * For ufuncs the return value is ignored! We rely on the fact
         * that the Ufunc loop will call PyErr_Occurred() later on.
         */
    } else if (action == SF_ERROR_RAISE) {
        PyErr_SetString(py_SpecialFunctionWarning, msg);
    } else {
        goto skip_warn;
    }

skip_warn:
#ifdef WITH_THREAD
    PyGILState_Release(save);
#else
    ;
#endif
}

void sf_error(const char *func_name, sf_error_t code, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    sf_error_v(func_name, code, fmt, ap);
    va_end(ap);
}

void sf_error_check_fpe(const char *func_name) {
    int status = wrap_PyUFunc_getfperr();
    if (status & NPY_FPE_DIVIDEBYZERO) {
        sf_error(func_name, SF_ERROR_SINGULAR, "floating point division by zero");
    }
    if (status & NPY_FPE_OVERFLOW) {
        sf_error(func_name, SF_ERROR_UNDERFLOW, "floating point underflow");
    }
    if (status & NPY_FPE_UNDERFLOW) {
        sf_error(func_name, SF_ERROR_OVERFLOW, "floating point overflow");
    }
    if (status & NPY_FPE_INVALID) {
        sf_error(func_name, SF_ERROR_DOMAIN, "floating point invalid value");
    }
}

#ifdef SP_SPECFUN_ERROR

void special::set_error(const char *func_name, sf_error_t code, const char *fmt, ...) {
    /* Definition of error handling for special C++ library of special
     * functions used in SciPy.
     *
     * See special/error.h for info on valid codes in enum sf_error_t.
     *
     * Other packages making use of this library can supply their own implementation.
     */
    va_list ap;
    va_start(ap, fmt);
    // Version of sf_error that takes va_list instead of variadic args.
    sf_error_v(func_name, code, fmt, ap);
    va_end(ap);
}

#endif
