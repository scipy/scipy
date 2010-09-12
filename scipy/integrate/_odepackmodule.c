/*
    Multipack project.
 */
#include "multipack.h"
static PyObject *odepack_error;
#include "__odepack.h"
static struct PyMethodDef odepack_module_methods[] = {
{"odeint", (PyCFunction) odepack_odeint, METH_VARARGS|METH_KEYWORDS, doc_odeint},
{NULL,		NULL, 0, NULL}
};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_odepack",
    NULL,
    -1,
    odepack_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *PyInit__odepack(void)
{
    PyObject *m, *d, *s;

    m = PyModule_Create(&moduledef);
    import_array();
    d = PyModule_GetDict(m);

    s = PyUString_FromString(" 1.9 ");
    PyDict_SetItemString(d, "__version__", s);
    odepack_error = PyErr_NewException ("odpack.error", NULL, NULL);
    Py_DECREF(s);
    PyDict_SetItemString(d, "error", odepack_error);
    if (PyErr_Occurred()) {
        Py_FatalError("can't initialize module odepack");
    }
    return m;
}
#else
PyMODINIT_FUNC init_odepack(void) {
  PyObject *m, *d, *s;
  m = Py_InitModule("_odepack", odepack_module_methods);
  import_array();
  d = PyModule_GetDict(m);

  s = PyUString_FromString(" 1.9 ");
  PyDict_SetItemString(d, "__version__", s);
  odepack_error = PyErr_NewException ("odepack.error", NULL, NULL);
  Py_DECREF(s);
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module odepack");
}
#endif    
