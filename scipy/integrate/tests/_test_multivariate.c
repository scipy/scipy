#include <Python.h>

#include "math.h"

const double PI = 3.141592653589793238462643383279502884;

static double
_multivariate_typical(int n, double *args)
{
    return cos(args[1] * args[0] - args[2] * sin(args[0])) / PI;
}

static double
_multivariate_indefinite(int n, double *args)
{
    return -exp(-args[0]) * log(args[0]);
}

static double
_multivariate_sin(int n, double *args)
{
    return sin(args[0]);
}

static double
_sin_0(double x, void *user_data)
{
    return sin(x);
}

static double
_sin_1(int ndim, double *x, void *user_data)
{
    return sin(x[0]);
}

static double
_sin_2(double x)
{
    return sin(x);
}

static double
_sin_3(int ndim, double *x)
{
    return sin(x[0]);
}


typedef struct {
    char *name;
    void *ptr;
} routine_t;


static const routine_t routines[] = {
    {"_multivariate_typical", &_multivariate_typical},
    {"_multivariate_indefinite", &_multivariate_indefinite},
    {"_multivariate_sin", &_multivariate_sin},
    {"_sin_0", &_sin_0},
    {"_sin_1", &_sin_1},
    {"_sin_2", &_sin_2},
    {"_sin_3", &_sin_3}
};


static int create_capsules(PyObject *module)
{
    PyObject *d, *capsule = NULL;
    int i;

    d = PyModule_GetDict(module);
    if (d == NULL) {
        goto fail;
    }

    for (i = 0; i < sizeof(routines) / sizeof(routine_t); ++i) {
        capsule = PyCapsule_New(routines[i].ptr, NULL, NULL);
        if (capsule == NULL) {
            goto fail;
        }

        if (PyDict_SetItemString(d, routines[i].name, capsule)) {
            goto fail;
        }

        Py_DECREF(capsule);
        capsule = NULL;
    }

    Py_XDECREF(capsule);
    return 0;

fail:
    Py_XDECREF(capsule);
    return -1;
}


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_test_multivariate",
    NULL,
    -1,
    NULL, /* Empty methods section */
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__test_multivariate(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (m == NULL) {
        return NULL;
    }
    if (create_capsules(m)) {
        Py_DECREF(m);
        return NULL;
    }
    return m;
}

#else

PyMODINIT_FUNC
init_test_multivariate(void)
{
    PyObject *m;
    m = Py_InitModule("_test_multivariate", NULL);
    if (m == NULL) {
        return;
    }
    create_capsules(m);
}
#endif
