#include "Python.h"
#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"

static PyObject *scipy_special_SpecialFunctionWarning = NULL;

void scipy_special_raise_warning(char *fmt, ...)
{
    NPY_ALLOW_C_API_DEF
    char msg[1024];
    va_list ap;

    va_start(ap, fmt);
    PyOS_vsnprintf(msg, 1024, fmt, ap);
    va_end(ap);

    NPY_ALLOW_C_API
    PyErr_Warn(scipy_special_SpecialFunctionWarning, msg);
    NPY_DISABLE_C_API
}
