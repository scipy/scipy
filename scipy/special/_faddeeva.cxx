#include <complex>

extern std::complex<double> Faddeeva_w(std::complex<double> z, double relerr);

extern "C" {

#include <Python.h>
#include <math.h>

#include "numpy/npy_math.h"
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"

#define RELERR 0    /* machine precision */

static void
wofz_loop(char **args, npy_intp *dimensions,
          npy_intp* steps, void* data)
{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *in = args[0], *out = args[1];
    npy_intp in_step = steps[0], out_step = steps[1];

    for (i = 0; i < n; i++) {
        std::complex<double> z(((npy_cdouble*)in)->real,
                               ((npy_cdouble*)in)->imag);
        std::complex<double> w = Faddeeva_w(z, RELERR);

        ((npy_cdouble*)out)->real = real(w);
        ((npy_cdouble*)out)->imag = imag(w);

        in += in_step;
        out += out_step;
    }
}

/*
 * Definitions for the ufuncs.
 */

static PyUFuncGenericFunction wofz_funcs[1] = {&wofz_loop};
static char types[2] = {NPY_CDOUBLE, NPY_CDOUBLE};
static void *data[1] = {NULL};

/* Module definition */

static PyMethodDef module_methods[] = {
    { NULL, NULL, 0, NULL }
};

static void _init_funcs(PyObject *m)
{
    PyObject *f, *d;

    d = PyModule_GetDict(m);

    f = PyUFunc_FromFuncAndData(wofz_funcs, data, types, 3, 1, 1,
                                PyUFunc_None, "wofz", NULL , 0);
    PyDict_SetItemString(d, "wofz", f);
    Py_DECREF(f);
}
    
#if PY_VERSION_HEX >= 0x03000000

static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_faddeeva",
    NULL,
    -1,
    module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__faddeeva()
{
    PyObject *m, *f, *d;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    import_array();
    import_umath();
    _init_funcs(m);
    return m;
}

#else

PyMODINIT_FUNC
init_faddeeva()
{
    PyObject *m, *f,  *d;
    m  = Py_InitModule("_faddeeva", module_methods);
    if (m == NULL) {
        return;
    }
    import_array();
    import_umath();
    _init_funcs(m);
}

#endif

}
