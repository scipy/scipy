#ifdef __CPLUSPLUS__
extern "C" {
#endif
#include "Python.h"
static PyMethodDef module_methods[] = { {NULL,NULL} };
DL_EXPORT(void) initatlas_version(void) {
  PyObject *m = NULL;
#if defined(NO_ATLAS_INFO)
  printf("NO ATLAS INFO AVAILABLE\n");
#else
  void ATL_buildinfo(void);
  ATL_buildinfo();
#endif
  m = Py_InitModule("atlas_version", module_methods);
#if defined(ATLAS_INFO)
  {
    PyObject *d = PyModule_GetDict(m);
    PyDict_SetItemString(d,"ATLAS_VERSION",PyString_FromString(ATLAS_INFO));
  }
#endif
}
#ifdef __CPLUSCPLUS__
}
#endif
