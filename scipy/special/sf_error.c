#include <Python.h>

#include <stdlib.h>
#include <stdarg.h>

#include "sf_error.h"

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

/* If this isn't volatile clang tries to optimize it away */
static volatile sf_action_t sf_error_actions[] = {
    SF_ERROR_IGNORE, /* SF_ERROR_OK */
    SF_ERROR_IGNORE, /* SF_ERROR_SINGULAR */
    SF_ERROR_IGNORE, /* SF_ERROR_UNDERFLOW */
    SF_ERROR_IGNORE, /* SF_ERROR_OVERFLOW */
    SF_ERROR_IGNORE, /* SF_ERROR_SLOW */
    SF_ERROR_IGNORE, /* SF_ERROR_LOSS */
    SF_ERROR_IGNORE, /* SF_ERROR_NO_RESULT */
    SF_ERROR_IGNORE, /* SF_ERROR_DOMAIN */
    SF_ERROR_IGNORE, /* SF_ERROR_ARG */
    SF_ERROR_IGNORE, /* SF_ERROR_OTHER */
    SF_ERROR_IGNORE  /* SF_ERROR__LAST */
};

extern int wrap_PyUFunc_getfperr(void);


void sf_error_set_action(sf_error_t code, sf_action_t action)
{
    sf_error_actions[(int)code] = action;
}


sf_action_t sf_error_get_action(sf_error_t code)
{
    return sf_error_actions[(int)code];
}


void sf_error_v(const char *func_name, sf_error_t code, const char *fmt, va_list ap)
{
    /* Internal function which takes a va_list instead of variadic args.
     * Makes this easier to wrap in error handling used in special C++
     * namespace for special function kernels provided by SciPy. */
    PyGILState_STATE save;
    PyObject *scipy_special = NULL;
    char msg[2048], info[1024];
    static PyObject *py_SpecialFunctionWarning = NULL;
    sf_action_t action;

    if ((int)code < 0 || (int)code >= 10) {
	code = SF_ERROR_OTHER;
    }
    action = sf_error_get_action(code);
    if (action == SF_ERROR_IGNORE) {
        return;
    }

    if (func_name == NULL) {
        func_name = "?";
    }

    if (fmt != NULL && fmt[0] != '\0') {
        PyOS_vsnprintf(info, 1024, fmt, ap);
        PyOS_snprintf(msg, 2048, "scipy.special/%s: (%s) %s",
                      func_name, sf_error_messages[(int)code], info);
    }
    else {
        PyOS_snprintf(msg, 2048, "scipy.special/%s: %s",
                      func_name, sf_error_messages[(int)code]);
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
	py_SpecialFunctionWarning =
	    PyObject_GetAttrString(scipy_special, "SpecialFunctionWarning");
    }
    else if (action == SF_ERROR_RAISE) {
	py_SpecialFunctionWarning =
	    PyObject_GetAttrString(scipy_special, "SpecialFunctionError");
    }
    else {
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
    }
    else if (action == SF_ERROR_RAISE) {
	PyErr_SetString(py_SpecialFunctionWarning, msg);
    }
    else {
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


#define UFUNC_FPE_DIVIDEBYZERO  1
#define UFUNC_FPE_OVERFLOW      2
#define UFUNC_FPE_UNDERFLOW     4
#define UFUNC_FPE_INVALID       8

void sf_error_check_fpe(const char *func_name)
{
    int status;
    status = wrap_PyUFunc_getfperr();
    if (status & UFUNC_FPE_DIVIDEBYZERO) {
        sf_error(func_name, SF_ERROR_SINGULAR, "floating point division by zero");
    }
    if (status & UFUNC_FPE_UNDERFLOW) {
        sf_error(func_name, SF_ERROR_UNDERFLOW, "floating point underflow");
    }
    if (status & UFUNC_FPE_OVERFLOW) {
        sf_error(func_name, SF_ERROR_OVERFLOW, "floating point overflow");
    }
    if (status & UFUNC_FPE_INVALID) {
        sf_error(func_name, SF_ERROR_DOMAIN, "floating point invalid value");
    }
}
