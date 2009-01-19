#ifdef __CPLUSPLUS__
extern "C" {
#endif

#include "Python.h"

void ATL_buildinfo(void);

static char doc_version[] = "Print ATLAS build info.";

static  PyObject *version(PyObject* self, PyObject* args)
{
#if defined(NO_ATLAS_INFO)
        printf("NO ATLAS INFO AVAILABLE\n");
#else
        ATL_buildinfo();
#endif

        Py_INCREF(Py_None);
        return Py_None;
}

static PyMethodDef module_methods[] = { 
	{"version", version, METH_VARARGS, doc_version},
        {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initatlas_version(void) 
{
        PyObject *m = NULL;
        m = Py_InitModule("atlas_version", module_methods);
#if defined(ATLAS_INFO)
        {
                PyObject *d = PyModule_GetDict(m);
                PyDict_SetItemString(d, "ATLAS_VERSION",
                                     PyString_FromString(ATLAS_INFO));
        }
#endif
}

#ifdef __CPLUSCPLUS__
}
#endif
