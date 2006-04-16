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
void init_odepack() {
  PyObject *m, *d, *s;
  m = Py_InitModule("_odepack", odepack_module_methods);
  import_array();
  d = PyModule_GetDict(m);

  s = PyString_FromString(" 1.9 ");
  PyDict_SetItemString(d, "__version__", s);
  odepack_error = PyErr_NewException ("odepack.error", NULL, NULL);
  Py_DECREF(s);
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module odepack");
}
        
