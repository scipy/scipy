#ifdef __CPLUSPLUS__
extern "C" {
#endif
#include "Python.h"
static PyMethodDef module_methods[] = { {NULL,NULL} };
DL_EXPORT(void) initatlas_version(void) {
#if defined(NO_ATLAS_INFO)
  printf("NO ATLAS INFO AVAILABLE\n");
#else
#if defined(ATLAS_INFO)
  printf("ATLAS_VERSION=" ATLAS_INFO "\n");
#else
  void ATL_buildinfo(void);
  ATL_buildinfo();
#endif
#endif
  Py_InitModule("atlas_version", module_methods);
}
#ifdef __CPLUSCPLUS__
}
#endif
