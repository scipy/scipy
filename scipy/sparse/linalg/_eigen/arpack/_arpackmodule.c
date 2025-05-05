/*
 * Python bindings for SciPy usage
 */

#ifndef __ARPACKMODULE_H
#define __ARPACKMODULE_H

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include "ARPACK/_arpack.h"

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); return NULL;}
static PyObject* arpack_error;

// ARPACK prototypes
struct ARPACK_arnoldi_update_vars_s;
struct ARPACK_arnoldi_update_vars_d;

void snaupd(struct ARPACK_arnoldi_update_vars_s *V, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void dnaupd(struct ARPACK_arnoldi_update_vars_d *V, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);


static PyObject*
snaupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    PyDictObject* input_dict = NULL;
    PyArrayObject* ap_resid = NULL;
    PyArrayObject* ap_v=NULL;
    PyArrayObject* ap_iparam=NULL;
    PyArrayObject* ap_ipntr=NULL;
    PyArrayObject* ap_workd=NULL;
    PyArrayObject* ap_workl=NULL;

    // Process input arguments
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!",
        &PyDict_Type, (PyObject **)&input_dict,  // O!
        &PyArray_Type, (PyObject **)&ap_resid,   // O!
        &PyArray_Type, (PyObject **)&ap_v,       // O!
        &PyArray_Type, (PyObject **)&ap_iparam,  // O!
        &PyArray_Type, (PyObject **)&ap_ipntr,   // O!
        &PyArray_Type, (PyObject **)&ap_workd,   // O!
        &PyArray_Type, (PyObject **)&ap_workl    // O!
        )
    )
    {
        return NULL;
    }

    // Check the input types
    if (!(PyArray_IS_C_CONTIGUOUS(ap_resid)) || (PyArray_TYPE(ap_resid) != NPY_FLOAT32))
    {
        PyErr_SetString(arpack_error, "resid must be a contiguous float32 array.");
        return NULL;
    }
    if (!(PyArray_IS_F_CONTIGUOUS(ap_v)) || (PyArray_TYPE(ap_v) != NPY_FLOAT32))
    {
        PyErr_SetString(arpack_error, "v must be a F-contiguous float32 array.");
        return NULL;
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_iparam)) || (PyArray_TYPE(ap_iparam) != NPY_INT32))
    {
        PyErr_SetString(arpack_error, "iparam must be a contiguous int32 array.");
        return NULL;
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_ipntr)) || (PyArray_TYPE(ap_ipntr) != NPY_INT32))
    {
        PyErr_SetString(arpack_error, "ipntr must be a contiguous int32 array.");
        return NULL;
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_workd)) || (PyArray_TYPE(ap_workd) != NPY_FLOAT32))
    {
        PyErr_SetString(arpack_error, "workd must be a contiguous float32 array.");
        return NULL;
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_workl)) || (PyArray_TYPE(ap_workl) != NPY_FLOAT32))
    {
        PyErr_SetString(arpack_error, "workl must be a contiguous float32 array.");
        return NULL;
    }

    if (PyArray_NDIM(ap_resid) != 1)
    {
        PyErr_SetString(arpack_error, "resid must be a 1D array.");
        return NULL;
    }
    if (PyArray_NDIM(ap_v) != 2)
    {
        PyErr_SetString(arpack_error, "v must be a 2D array.");
        return NULL;
    }
    if (PyArray_NDIM(ap_iparam) != 1)
    {
        PyErr_SetString(arpack_error, "iparam must be a 1D array.");
        return NULL;
    }
    if (PyArray_NDIM(ap_ipntr) != 1)
    {
        PyErr_SetString(arpack_error, "ipntr must be a 1D array.");
        return NULL;
    }
    if (PyArray_NDIM(ap_workd) != 1)
    {
        PyErr_SetString(arpack_error, "workd must be a 1D array.");
        return NULL;
    }
    if (PyArray_NDIM(ap_workl) != 1)
    {
        PyErr_SetString(arpack_error, "workl must be a 1D array.");
        return NULL;
    }

    int* iparam = PyArray_DATA(ap_iparam);
    int* ipntr = PyArray_DATA(ap_ipntr);
    float* resid = PyArray_DATA(ap_resid);
    float* v = PyArray_DATA(ap_v);
    float* workd = PyArray_DATA(ap_workd);
    float* workl = PyArray_DATA(ap_workl);

    // Map the input dict to the ARPACK structure

    #define STRUCT_FLOAT_FIELD_NAMES X(tol) X(getv0_rnorm0) X(aitr_betaj) X(aitr_rnorm1) X(aitr_wnorm) X(aup2_rnorm)
    #define STRUCT_INT_FIELD_NAMES X(ido) X(which) X(bmat) X(info) X(iter) X(maxiter) X(mode) \
                                   X(n) X(nconv) X(ncv) X(nev) X(np) X(numop) X(numpb) X(numreo) \
                                   X(shift) X(getv0_first) X(getv0_iter) X(getv0_itry) X(getv0_orth) \
                                   X(aitr_iter) X(aitr_j) X(aitr_orth1) X(aitr_orth2) X(aitr_restart) \
                                   X(aitr_step3) X(aitr_step4) X(aitr_ierr) X(aup2_initv) X(aup2_iter) \
                                   X(aup2_getv0) X(aup2_cnorm) X(aup2_kplusp) X(aup2_nev0) X(aup2_np0) \
                                   X(aup2_numcnv) X(aup2_update) X(aup2_ushift)
    #define STRUCT_FIELD_NAMES STRUCT_INT_FIELD_NAMES STRUCT_FLOAT_FIELD_NAMES




    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for floats and ints.

    // Declare and Initialize the ARPACK struct that will be populated from dict with zeros
    struct ARPACK_arnoldi_update_vars_s Vars = {0};

    #define X(name) Vars.name = 0;
    STRUCT_FIELD_NAMES
    #undef X

    // PyDict_GetItemString returns a borrowed reference.
    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = (float)PyFloat_AsDouble(name##_obj);
        STRUCT_FLOAT_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = (int)PyLong_AsLong(name##_obj);
        STRUCT_INT_FIELD_NAMES
    #undef X

    // Call ARPACK function
    snaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl);

    // Unpack the struct back to the dictionary
    #define X(name) do { \
            PyObject* tmp_##name = PyFloat_FromDouble((double)Vars.name); \
            if ((!tmp_##name) || (PyDict_SetItemString(input_dict, #name, tmp_##name) < 0)) { \
            Py_XDECREF(tmp_##name); \
            PYERR(arpack_error, "Setting '" #name "' failed."); \
            } \
            Py_DECREF(tmp_##name); \
        } while (0);
        STRUCT_FLOAT_FIELD_NAMES
    #undef X

    #define X(name) do { \
            PyObject* tmp_##name = PyLong_FromLong((long)Vars.name); \
            if ((!tmp_##name) || (PyDict_SetItemString(input_dict, #name, tmp_##name) < 0)) { \
                Py_XDECREF(tmp_##name); \
                PYERR(arpack_error, "Setting '" #name "' failed."); \
            } \
            Py_DECREF(tmp_##name); \
        } while (0);
        STRUCT_INT_FIELD_NAMES
    #undef X
    #undef STRUCT_FIELD_NAMES
    #undef STRUCT_INT_FIELD_NAMES
    #undef STRUCT_FLOAT_FIELD_NAMES


    Py_RETURN_NONE;
}



static PyObject*
dnaupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    PyDictObject* input_dict = NULL;
    PyArrayObject* ap_resid = NULL;
    PyArrayObject* ap_v=NULL;
    PyArrayObject* ap_iparam=NULL;
    PyArrayObject* ap_ipntr=NULL;
    PyArrayObject* ap_workd=NULL;
    PyArrayObject* ap_workl=NULL;

    // Process input arguments
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!",
        &PyDict_Type, (PyObject **)&input_dict,  // O!
        &PyArray_Type, (PyObject **)&ap_resid,   // O!
        &PyArray_Type, (PyObject **)&ap_v,       // O!
        &PyArray_Type, (PyObject **)&ap_iparam,  // O!
        &PyArray_Type, (PyObject **)&ap_ipntr,   // O!
        &PyArray_Type, (PyObject **)&ap_workd,   // O!
        &PyArray_Type, (PyObject **)&ap_workl    // O!
        )
    )
    {
        return NULL;
    }

    // Check the input types
    if (!(PyArray_IS_C_CONTIGUOUS(ap_resid)) || (PyArray_TYPE(ap_resid) != NPY_FLOAT64))
    {
        PyErr_SetString(arpack_error, "resid must be a contiguous float32 array.");
        return NULL;
    }
    if (!(PyArray_IS_F_CONTIGUOUS(ap_v)) || (PyArray_TYPE(ap_v) != NPY_FLOAT64))
    {
        PyErr_SetString(arpack_error, "v must be a F-contiguous float32 array.");
        return NULL;
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_iparam)) || (PyArray_TYPE(ap_iparam) != NPY_INT32))
    {
        PyErr_SetString(arpack_error, "iparam must be a contiguous int32 array.");
        return NULL;
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_ipntr)) || (PyArray_TYPE(ap_ipntr) != NPY_INT32))
    {
        PyErr_SetString(arpack_error, "ipntr must be a contiguous int32 array.");
        return NULL;
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_workd)) || (PyArray_TYPE(ap_workd) != NPY_FLOAT64))
    {
        PyErr_SetString(arpack_error, "workd must be a contiguous float32 array.");
        return NULL;
    }
    if (!(PyArray_IS_C_CONTIGUOUS(ap_workl)) || (PyArray_TYPE(ap_workl) != NPY_FLOAT64))
    {
        PyErr_SetString(arpack_error, "workl must be a contiguous float32 array.");
        return NULL;
    }

    if (PyArray_NDIM(ap_resid) != 1)
    {
        PyErr_SetString(arpack_error, "resid must be a 1D array.");
        return NULL;
    }
    if (PyArray_NDIM(ap_v) != 2)
    {
        PyErr_SetString(arpack_error, "v must be a 2D array.");
        return NULL;
    }
    if (PyArray_NDIM(ap_iparam) != 1)
    {
        PyErr_SetString(arpack_error, "iparam must be a 1D array.");
        return NULL;
    }
    if (PyArray_NDIM(ap_ipntr) != 1)
    {
        PyErr_SetString(arpack_error, "ipntr must be a 1D array.");
        return NULL;
    }
    if (PyArray_NDIM(ap_workd) != 1)
    {
        PyErr_SetString(arpack_error, "workd must be a 1D array.");
        return NULL;
    }
    if (PyArray_NDIM(ap_workl) != 1)
    {
        PyErr_SetString(arpack_error, "workl must be a 1D array.");
        return NULL;
    }

    int* iparam = PyArray_DATA(ap_iparam);
    int* ipntr = PyArray_DATA(ap_ipntr);
    double* resid = PyArray_DATA(ap_resid);
    double* v = PyArray_DATA(ap_v);
    double* workd = PyArray_DATA(ap_workd);
    double* workl = PyArray_DATA(ap_workl);

    // Map the input dict to the ARPACK structure

    #define STRUCT_DOUBLE_FIELD_NAMES X(tol) X(getv0_rnorm0) X(aitr_betaj) X(aitr_rnorm1) X(aitr_wnorm) X(aup2_rnorm)
    #define STRUCT_INT_FIELD_NAMES X(ido) X(which) X(bmat) X(info) X(iter) X(maxiter) X(mode) \
                                   X(n) X(nconv) X(ncv) X(nev) X(np) X(numop) X(numpb) X(numreo) \
                                   X(shift) X(getv0_first) X(getv0_iter) X(getv0_itry) X(getv0_orth) \
                                   X(aitr_iter) X(aitr_j) X(aitr_orth1) X(aitr_orth2) X(aitr_restart) \
                                   X(aitr_step3) X(aitr_step4) X(aitr_ierr) X(aup2_initv) X(aup2_iter) \
                                   X(aup2_getv0) X(aup2_cnorm) X(aup2_kplusp) X(aup2_nev0) X(aup2_np0) \
                                   X(aup2_numcnv) X(aup2_update) X(aup2_ushift)
    #define STRUCT_FIELD_NAMES STRUCT_INT_FIELD_NAMES STRUCT_DOUBLE_FIELD_NAMES




    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for floats and ints.

    // Declare and Initialize the ARPACK struct that will be populated from dict with zeros
    struct ARPACK_arnoldi_update_vars_s Vars = {0};

    #define X(name) Vars.name = 0;
    STRUCT_FIELD_NAMES
    #undef X

    // PyDict_GetItemString returns a borrowed reference.
    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = PyFloat_AsDouble(name##_obj);
        STRUCT_DOUBLE_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = (int)PyLong_AsLong(name##_obj);
        STRUCT_INT_FIELD_NAMES
    #undef X

    // Call ARPACK function
    snaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl);

    // Unpack the struct back to the dictionary
    #define X(name) do { \
            PyObject* tmp_##name = PyFloat_FromDouble(Vars.name); \
            if ((!tmp_##name) || (PyDict_SetItemString(input_dict, #name, tmp_##name) < 0)) { \
            Py_XDECREF(tmp_##name); \
            PYERR(arpack_error, "Setting '" #name "' failed."); \
            } \
            Py_DECREF(tmp_##name); \
        } while (0);
        STRUCT_DOUBLE_FIELD_NAMES
    #undef X

    #define X(name) do { \
            PyObject* tmp_##name = PyLong_FromLong((long)Vars.name); \
            if ((!tmp_##name) || (PyDict_SetItemString(input_dict, #name, tmp_##name) < 0)) { \
                Py_XDECREF(tmp_##name); \
                PYERR(arpack_error, "Setting '" #name "' failed."); \
            } \
            Py_DECREF(tmp_##name); \
        } while (0);
        STRUCT_INT_FIELD_NAMES
    #undef X
    #undef STRUCT_FIELD_NAMES
    #undef STRUCT_INT_FIELD_NAMES
    #undef STRUCT_DOUBLE_FIELD_NAMES


    Py_RETURN_NONE;
}




static char doc_snaupd[] = ("");
static char doc_dnaupd[] = ("");


// Sentinel terminated method list.
static struct PyMethodDef arpacklib_module_methods[] = {
    {"snaupd_wrap", snaupd_wrap, METH_VARARGS, doc_snaupd},
    {"dnaupd_wrap", dnaupd_wrap, METH_VARARGS, doc_dnaupd},
    {NULL, NULL, 0, NULL}
  };


struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_arpacklib",
    NULL,
    -1,
    arpacklib_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};


PyMODINIT_FUNC
PyInit__arpacklib(void)
{
    import_array();

    PyObject* module = PyModule_Create(&moduledef);
    if (module == NULL) { return NULL; }
    PyObject* mdict = PyModule_GetDict(module);
    if (mdict == NULL) { return NULL; }
    arpack_error = PyErr_NewException("_arpacklib.error", NULL, NULL);
    if (arpack_error == NULL) { return NULL; }
    if (PyDict_SetItemString(mdict, "error", arpack_error)) { return NULL; }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}



#endif // ifndef __ARPACKMODULE_H
