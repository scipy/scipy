#include "Python.h"
#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"


void scipy_special_raise_warning(char *fmt, ...)
{
  static PyObject *warning = NULL;

    NPY_ALLOW_C_API_DEF
    char msg[1024];
    va_list ap;

    va_start(ap, fmt);
    PyOS_vsnprintf(msg, 1024, fmt, ap);
    va_end(ap);

    NPY_ALLOW_C_API
      if (warning == NULL) {
        /* Get the PyObject pointer from scipy.special._cephes */
	PyObject* sp_basic = PyImport_ImportModule("scipy.special.basic");
	if(sp_basic != NULL) {
	  warning = PyObject_GetAttrString(sp_basic, "SpecialFunctionWarning");
	  Py_DECREF(sp_basic);
	  }
      }
    PyErr_Warn(warning, msg);
    NPY_DISABLE_C_API
}
