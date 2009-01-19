#include "Python.h"

static PyObject* version(PyObject* self, PyObject* dummy)
{
#if defined(NO_ATLAS_INFO)
    printf("NO ATLAS INFO AVAILABLE\n");
#else
    void ATL_buildinfo(void);
    ATL_buildinfo();
#endif

    Py_INCREF(Py_None);
    return Py_None;
}

static char version_doc[] = "Print the build info from atlas.";

static PyMethodDef module_methods[] = {
    {"version", version, METH_VARARGS, version_doc},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initatlas_version(void)
{
    PyObject *m = NULL;
    m = Py_InitModule("atlas_version", module_methods);
#if defined(ATLAS_INFO)
    {
        PyObject *d = PyModule_GetDict(m);
        PyDict_SetItemString(d,"ATLAS_VERSION",PyString_FromString(ATLAS_INFO));
    }
#endif
}
