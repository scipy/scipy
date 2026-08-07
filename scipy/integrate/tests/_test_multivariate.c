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


static int create_pointers(PyObject *module)
{
    PyObject *d, *obj = NULL;
    size_t i;

    d = PyModule_GetDict(module);
    if (d == NULL) {
        goto fail;
    }

    for (i = 0; i < sizeof(routines) / sizeof(routine_t); ++i) {
        obj = PyLong_FromVoidPtr(routines[i].ptr);
        if (obj == NULL) {
            goto fail;
        }

        if (PyDict_SetItemString(d, routines[i].name, obj)) {
            goto fail;
        }

        Py_DECREF(obj);
        obj = NULL;
    }

    Py_XDECREF(obj);
    return 0;

fail:
    Py_XDECREF(obj);
    return -1;
}


static struct PyModuleDef_Slot _test_multivariate_slots[] = {
    {Py_mod_exec, create_pointers},
    // signal that this module can be imported in isolated subinterpreters
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    // signal that this module supports running without an active GIL
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};


static struct PyModuleDef moduledef = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_test_multivariate",
    .m_size = 0,
    .m_methods = NULL,
    .m_slots = _test_multivariate_slots,
};

PyMODINIT_FUNC
PyInit__test_multivariate(void)
{
    return PyModuleDef_Init(&moduledef);
}
