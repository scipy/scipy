#include "Python.h"
#include "numpy/npy_3kcompat.h"

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

#if PY_VERSION_HEX >= 0x03000000

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "atlas_version",
    NULL,
    -1,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *PyInit_atlas_version(void)
{
#define RETVAL m
    PyObject *m;
    m = PyModule_Create(&moduledef);
#else
#define RETVAL
PyMODINIT_FUNC initatlas_version(void)
{
    PyObject *m = NULL;
    m = Py_InitModule("atlas_version", module_methods);
#endif
    if (m == NULL) {
        return RETVAL;
    }
#if defined(ATLAS_INFO)
    {
        PyObject *d = PyModule_GetDict(m);
        PyDict_SetItemString(d,"ATLAS_VERSION",
                             PyUString_FromString(ATLAS_INFO));
    }
#endif
    return RETVAL;
}
