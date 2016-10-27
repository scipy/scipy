/*
 * Control error handling for special functions. The code is somewhat
 * convoluted since it has to support both special and cython_special,
 * which have different error handling models. In particular, special
 * invokes the usual Python warnings machinery whereas cython_special
 * uses a global errno style variable.
 */
#include <stdlib.h>
#include <stdarg.h>

#include <Python.h>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

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

extern int wrap_PyUFunc_getfperr(void);

#ifdef CYTHON_SPECIAL
/*
 * Global variable for tracking errors in cython_special. It needs to
 * be shared between the C functions and the C++ functions, so declare
 * it extern here and put the actual declaration in cython_special. If
 * we have openmp make it threadsafe for openmp threads. This includes
 * uses like cython.parallel.
 */
extern sf_error_t sf_errno;
#ifdef HAVE_OPENMP
#pragma omp threadprivate(sf_errno)
#endif

#endif

/*
 * Control whether warnings are printed in special
 */
static int print_error_messages = 0;


int sf_error_set_print(int flag)
{
#ifdef CYTHON_SPECIAL
    /* This is a noop in cython_special */
    return -1;
#else
    int old_flag = print_error_messages;
    print_error_messages = flag;
    return old_flag;
#endif
}


int sf_error_get_print()
{
#ifdef CYTHON_SPECIAL
    /* This is a noop in cython_special */
    return -1;
#else
    return print_error_messages;
#endif
}


sf_error_t sf_error_get_errno() {
#ifdef CYTHON_SPECIAL
    sf_error_t old_errno = sf_errno;
    sf_errno = SF_ERROR_OK;
    return old_errno;
#else
    /* This is a noop in special; return an invalid value */
    return SF_ERROR__LAST;
#endif
}


void sf_error(const char *func_name, sf_error_t code, const char *fmt, ...)
{
#ifdef CYTHON_SPECIAL
    if (code < 0 || code >= 10) {
	code = SF_ERROR_OTHER;
    }
    sf_errno = code;
#else
    char msg[2048], info[1024];
    static PyObject *py_SpecialFunctionWarning = NULL;
    va_list ap;

    if (!print_error_messages) {
        return;
    }

    if (func_name == NULL) {
        func_name = "?";
    }

    if ((int)code < 0 || (int)code >= 10) {
        code = SF_ERROR_OTHER;
    }

    if (fmt != NULL && fmt[0] != '\0') {
        va_start(ap, fmt);
        PyOS_vsnprintf(info, 1024, fmt, ap);
        va_end(ap);
        PyOS_snprintf(msg, 2048, "scipy.special/%s: (%s) %s",
                      func_name, sf_error_messages[(int)code], info);
    }
    else {
        PyOS_snprintf(msg, 2048, "scipy.special/%s: %s",
                      func_name, sf_error_messages[(int)code]);
    }

    {
#ifdef WITH_THREAD
        PyGILState_STATE save = PyGILState_Ensure();
#endif

        if (PyErr_Occurred())
            goto skip_warn;

        if (py_SpecialFunctionWarning == NULL) {
            PyObject *scipy_special = NULL;

            scipy_special = PyImport_ImportModule("scipy.special");
            if (scipy_special == NULL) {
                PyErr_Clear();
                goto skip_warn;
            }

            py_SpecialFunctionWarning = PyObject_GetAttrString(
                scipy_special, "SpecialFunctionWarning");
            if (py_SpecialFunctionWarning == NULL) {
                PyErr_Clear();
                goto skip_warn;
            }
        }

        if (py_SpecialFunctionWarning != NULL) {
            PyErr_WarnEx(py_SpecialFunctionWarning, msg, 1);
            /*
             * For ufuncs the return value is ignored! We rely on the
             * fact that the Ufunc loop will call PyErr_Occurred()
             * later on.
             */
        }

    skip_warn:
#ifdef WITH_THREAD
        PyGILState_Release(save);
#endif
    }
#endif /* CYTHON_SPECIAL */
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
