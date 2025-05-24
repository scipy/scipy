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

// The following macros are used to define the field names in the ARPACK struct.
#define STRUCT_INEXACT_FIELD_NAMES X(tol) X(getv0_rnorm0) X(aitr_betaj) X(aitr_rnorm1) X(aitr_wnorm) X(aup2_rnorm)
#define STRUCT_INT_FIELD_NAMES X(ido) X(which) X(bmat) X(info) X(iter) X(maxiter) X(mode) \
                               X(n) X(nconv) X(ncv) X(nev) X(np) X(numop) X(numpb) X(numreo) \
                               X(shift) X(getv0_first) X(getv0_iter) X(getv0_itry) X(getv0_orth) \
                               X(aitr_iter) X(aitr_j) X(aitr_orth1) X(aitr_orth2) X(aitr_restart) \
                               X(aitr_step3) X(aitr_step4) X(aitr_ierr) X(aup2_initv) X(aup2_iter) \
                               X(aup2_getv0) X(aup2_cnorm) X(aup2_kplusp) X(aup2_nev0) X(aup2_np0) \
                               X(aup2_numcnv) X(aup2_update) X(aup2_ushift)
#define STRUCT_FIELD_NAMES STRUCT_INT_FIELD_NAMES STRUCT_INEXACT_FIELD_NAMES


// ARPACK prototypes
struct ARPACK_arnoldi_update_vars_s;
struct ARPACK_arnoldi_update_vars_d;

void snaupd(struct ARPACK_arnoldi_update_vars_s *V, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void dnaupd(struct ARPACK_arnoldi_update_vars_d *V, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);
void cnaupd(struct ARPACK_arnoldi_update_vars_s *V, ARPACK_CPLXF_TYPE* resid, ARPACK_CPLXF_TYPE* v, int ldv, int* ipntr, ARPACK_CPLXF_TYPE* workd, ARPACK_CPLXF_TYPE* workl, float* rwork);
void znaupd(struct ARPACK_arnoldi_update_vars_d *V, ARPACK_CPLX_TYPE* resid, ARPACK_CPLX_TYPE* v, int ldv, int* ipntr, ARPACK_CPLX_TYPE* workd, ARPACK_CPLX_TYPE* workl, double* rwork);
void sneupd(struct ARPACK_arnoldi_update_vars_s *V, int rvec, int howmny, int* select, float* dr, float* di, float* z, int ldz, float sigmar, float sigmai, float* workev, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void dneupd(struct ARPACK_arnoldi_update_vars_d *V, int rvec, int howmny, int* select, double* dr, double* di, double* z, int ldz, double sigmar, double sigmai, double* workev, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);


static PyObject*
snaupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    PyObject* input_dict = NULL;
    PyArrayObject* ap_resid = NULL;
    PyArrayObject* ap_v=NULL;
    PyArrayObject* ap_ipntr=NULL;
    PyArrayObject* ap_workd=NULL;
    PyArrayObject* ap_workl=NULL;

    // Process input arguments
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!",
        &PyDict_Type, (PyObject **)&input_dict,  // O!
        &PyArray_Type, (PyObject **)&ap_resid,   // O!
        &PyArray_Type, (PyObject **)&ap_v,       // O!
        &PyArray_Type, (PyObject **)&ap_ipntr,   // O!
        &PyArray_Type, (PyObject **)&ap_workd,   // O!
        &PyArray_Type, (PyObject **)&ap_workl    // O!
        )
    )
    {
        return NULL;
    }

    int* ipntr = PyArray_DATA(ap_ipntr);
    float* resid = PyArray_DATA(ap_resid);
    float* v = PyArray_DATA(ap_v);
    float* workd = PyArray_DATA(ap_workd);
    float* workl = PyArray_DATA(ap_workl);

    // Map the input dict to the ARPACK structure

    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for floats and ints.

    // Declare and Initialize the ARPACK struct that will be populated from dict with zeros
    struct ARPACK_arnoldi_update_vars_s Vars = {0};

    #define X(name) Vars.name = 0;
    STRUCT_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = (float)PyFloat_AsDouble(name##_obj);
        STRUCT_INEXACT_FIELD_NAMES
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
        STRUCT_INEXACT_FIELD_NAMES
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

    Py_RETURN_NONE;
}


static PyObject*
dnaupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    PyObject* input_dict = NULL;
    PyArrayObject* ap_resid = NULL;
    PyArrayObject* ap_v=NULL;
    PyArrayObject* ap_ipntr=NULL;
    PyArrayObject* ap_workd=NULL;
    PyArrayObject* ap_workl=NULL;

    // Process input arguments
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!",
        &PyDict_Type, (PyObject **)&input_dict,  // O!
        &PyArray_Type, (PyObject **)&ap_resid,   // O!
        &PyArray_Type, (PyObject **)&ap_v,       // O!
        &PyArray_Type, (PyObject **)&ap_ipntr,   // O!
        &PyArray_Type, (PyObject **)&ap_workd,   // O!
        &PyArray_Type, (PyObject **)&ap_workl    // O!
        )
    )
    {
        return NULL;
    }

    int* ipntr = PyArray_DATA(ap_ipntr);
    double* resid = PyArray_DATA(ap_resid);
    double* v = PyArray_DATA(ap_v);
    double* workd = PyArray_DATA(ap_workd);
    double* workl = PyArray_DATA(ap_workl);

    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for floats and ints.

    // Declare and Initialize the ARPACK struct that will be populated from dict with zeros
    struct ARPACK_arnoldi_update_vars_d Vars = {0};

    #define X(name) Vars.name = 0;
    STRUCT_FIELD_NAMES
    #undef X

    #define X(name) \
    PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
    if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
    Vars.name = PyFloat_AsDouble(name##_obj);
    STRUCT_INEXACT_FIELD_NAMES
    #undef X

    #define X(name) \
    PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
    if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
    Vars.name = (int)PyLong_AsLong(name##_obj);
    STRUCT_INT_FIELD_NAMES
    #undef X

    // Call ARPACK function
    dnaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl);

    // Unpack the struct back to the dictionary
    #define X(name) do { \
            PyObject* tmp_##name = PyFloat_FromDouble(Vars.name); \
            if ((!tmp_##name) || (PyDict_SetItemString(input_dict, #name, tmp_##name) < 0)) { \
            Py_XDECREF(tmp_##name); \
            PYERR(arpack_error, "Setting '" #name "' failed."); \
            } \
            Py_DECREF(tmp_##name); \
        } while (0);
        STRUCT_INEXACT_FIELD_NAMES
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

    Py_RETURN_NONE;
}


static PyObject*
cnaupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    PyObject* input_dict = NULL;
    PyArrayObject* ap_resid = NULL;
    PyArrayObject* ap_v=NULL;
    PyArrayObject* ap_ipntr=NULL;
    PyArrayObject* ap_workd=NULL;
    PyArrayObject* ap_workl=NULL;
    PyArrayObject* ap_rwork=NULL;

    // Process input arguments
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!",
        &PyDict_Type, (PyObject **)&input_dict,  // O!
        &PyArray_Type, (PyObject **)&ap_resid,   // O!
        &PyArray_Type, (PyObject **)&ap_v,       // O!
        &PyArray_Type, (PyObject **)&ap_ipntr,   // O!
        &PyArray_Type, (PyObject **)&ap_workd,   // O!
        &PyArray_Type, (PyObject **)&ap_workl,   // O!
        &PyArray_Type, (PyObject **)&ap_rwork    // O!
        )
    )
    {
        return NULL;
    }

    int* ipntr = PyArray_DATA(ap_ipntr);
    ARPACK_CPLXF_TYPE* resid = PyArray_DATA(ap_resid);
    ARPACK_CPLXF_TYPE* v = PyArray_DATA(ap_v);
    ARPACK_CPLXF_TYPE* workd = PyArray_DATA(ap_workd);
    ARPACK_CPLXF_TYPE* workl = PyArray_DATA(ap_workl);
    float* rwork = PyArray_DATA(ap_rwork);

    // Map the input dict to the ARPACK structure

    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for floats and ints.

    // Declare and Initialize the ARPACK struct that will be populated from dict with zeros
    struct ARPACK_arnoldi_update_vars_s Vars = {0};

    #define X(name) Vars.name = 0;
    STRUCT_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = (float)PyFloat_AsDouble(name##_obj);
        STRUCT_INEXACT_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = (int)PyLong_AsLong(name##_obj);
        STRUCT_INT_FIELD_NAMES
    #undef X

    // Call ARPACK function
    cnaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl, rwork);

    // Unpack the struct back to the dictionary
    #define X(name) do { \
            PyObject* tmp_##name = PyFloat_FromDouble((double)Vars.name); \
            if ((!tmp_##name) || (PyDict_SetItemString(input_dict, #name, tmp_##name) < 0)) { \
            Py_XDECREF(tmp_##name); \
            PYERR(arpack_error, "Setting '" #name "' failed."); \
            } \
            Py_DECREF(tmp_##name); \
        } while (0);
        STRUCT_INEXACT_FIELD_NAMES
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

    Py_RETURN_NONE;
}


static PyObject*
znaupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    PyObject* input_dict = NULL;
    PyArrayObject* ap_resid = NULL;
    PyArrayObject* ap_v=NULL;
    PyArrayObject* ap_ipntr=NULL;
    PyArrayObject* ap_workd=NULL;
    PyArrayObject* ap_workl=NULL;
    PyArrayObject* ap_rwork=NULL;

    // Process input arguments
    if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!",
        &PyDict_Type, (PyObject **)&input_dict,  // O!
        &PyArray_Type, (PyObject **)&ap_resid,   // O!
        &PyArray_Type, (PyObject **)&ap_v,       // O!
        &PyArray_Type, (PyObject **)&ap_ipntr,   // O!
        &PyArray_Type, (PyObject **)&ap_workd,   // O!
        &PyArray_Type, (PyObject **)&ap_workl,   // O!
        &PyArray_Type, (PyObject **)&ap_rwork    // O!
        )
    )
    {
        return NULL;
    }

    int* ipntr = PyArray_DATA(ap_ipntr);
    ARPACK_CPLX_TYPE* resid = PyArray_DATA(ap_resid);
    ARPACK_CPLX_TYPE* v = PyArray_DATA(ap_v);
    ARPACK_CPLX_TYPE* workd = PyArray_DATA(ap_workd);
    ARPACK_CPLX_TYPE* workl = PyArray_DATA(ap_workl);
    double* rwork = PyArray_DATA(ap_rwork);

    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for floats and ints.

    // Declare and Initialize the ARPACK struct that will be populated from dict with zeros
    struct ARPACK_arnoldi_update_vars_d Vars = {0};

    #define X(name) Vars.name = 0;
    STRUCT_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = PyFloat_AsDouble(name##_obj);
        STRUCT_INEXACT_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = (int)PyLong_AsLong(name##_obj);
        STRUCT_INT_FIELD_NAMES
    #undef X

    // Call ARPACK function
    znaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl, rwork);

    // Unpack the struct back to the dictionary
    #define X(name) do { \
            PyObject* tmp_##name = PyFloat_FromDouble(Vars.name); \
            if ((!tmp_##name) || (PyDict_SetItemString(input_dict, #name, tmp_##name) < 0)) { \
            Py_XDECREF(tmp_##name); \
            PYERR(arpack_error, "Setting '" #name "' failed."); \
            } \
            Py_DECREF(tmp_##name); \
        } while (0);
        STRUCT_INEXACT_FIELD_NAMES
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

    Py_RETURN_NONE;
}


static PyObject*
sneupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{

    PyObject* input_dict = NULL;
    int want_ev = 0, howmny = 0, ldv = 0, ldz = 0;
    PyArrayObject* ap_select = NULL;
    float sigmar = 0.0;
    float sigmai = 0.0;
    PyArrayObject* ap_dr = NULL;
    PyArrayObject* ap_di = NULL;
    PyArrayObject* ap_v = NULL;
    PyArrayObject* ap_z = NULL;
    PyArrayObject* ap_workev = NULL;
    PyArrayObject* ap_resid = NULL;
    PyArrayObject* ap_ipntr = NULL;
    PyArrayObject* ap_workd = NULL;
    PyArrayObject* ap_workl = NULL;

    // Process input arguments
    if (!PyArg_ParseTuple(args, "O!iiO!O!O!O!ffO!O!O!O!O!O!",
        &PyDict_Type, (PyObject **)&input_dict,  // O!
        &want_ev,                                // i
        &howmny,                                 // i
        &PyArray_Type, (PyObject **)&ap_select,  // O!
        &PyArray_Type, (PyObject **)&ap_dr,      // O!
        &PyArray_Type, (PyObject **)&ap_di,      // O!
        &PyArray_Type, (PyObject **)&ap_z,       // O!
        &sigmar,                                 // f
        &sigmai,                                 // f
        &PyArray_Type, (PyObject **)&ap_workev,  // O!
        &PyArray_Type, (PyObject **)&ap_resid,   // O!
        &PyArray_Type, (PyObject **)&ap_v,       // O!
        &PyArray_Type, (PyObject **)&ap_ipntr,   // O!
        &PyArray_Type, (PyObject **)&ap_workd,   // O!
        &PyArray_Type, (PyObject **)&ap_workl    // O!
        )
    )
    {
        return NULL;
    }

    int* ipntr = PyArray_DATA(ap_ipntr);
    int* select = PyArray_DATA(ap_select);
    float* dr = PyArray_DATA(ap_dr);
    float* di = PyArray_DATA(ap_di);
    float* workev = PyArray_DATA(ap_workev);
    float* z = PyArray_DATA(ap_z);
    float* resid = PyArray_DATA(ap_resid);
    float* v = PyArray_DATA(ap_v);
    float* workd = PyArray_DATA(ap_workd);
    float* workl = PyArray_DATA(ap_workl);
    ldv = (int)PyArray_DIMS(ap_v)[0];
    ldz = (int)PyArray_DIMS(ap_z)[0];

    struct ARPACK_arnoldi_update_vars_s Vars = {0};

    #define X(name) Vars.name = 0;
    STRUCT_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = (float)PyFloat_AsDouble(name##_obj);
        STRUCT_INEXACT_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = (int)PyLong_AsLong(name##_obj);
        STRUCT_INT_FIELD_NAMES
    #undef X


    sneupd(&Vars, want_ev, howmny, select, dr, di, z, ldz, sigmar, sigmai, workev, resid, v, ldv, ipntr, workd, workl);

    Py_RETURN_NONE;

}

static PyObject*
dneupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{

    PyObject* input_dict = NULL;
    int want_ev = 0, howmny = 0, ldv = 0, ldz = 0;
    PyArrayObject* ap_select = NULL;
    double sigmar = 0.0;
    double sigmai = 0.0;
    PyArrayObject* ap_dr = NULL;
    PyArrayObject* ap_di = NULL;
    PyArrayObject* ap_v = NULL;
    PyArrayObject* ap_z = NULL;
    PyArrayObject* ap_workev = NULL;
    PyArrayObject* ap_resid = NULL;
    PyArrayObject* ap_ipntr = NULL;
    PyArrayObject* ap_workd = NULL;
    PyArrayObject* ap_workl = NULL;

    // Process input arguments
    if (!PyArg_ParseTuple(args, "O!iiO!O!O!O!ddO!O!O!O!O!O!",
        &PyDict_Type, (PyObject **)&input_dict,  // O!
        &want_ev,                                // i
        &howmny,                                 // i
        &PyArray_Type, (PyObject **)&ap_select,  // O!
        &PyArray_Type, (PyObject **)&ap_dr,      // O!
        &PyArray_Type, (PyObject **)&ap_di,      // O!
        &PyArray_Type, (PyObject **)&ap_z,       // O!
        &sigmar,                                 // d
        &sigmai,                                 // d
        &PyArray_Type, (PyObject **)&ap_workev,  // O!
        &PyArray_Type, (PyObject **)&ap_resid,   // O!
        &PyArray_Type, (PyObject **)&ap_v,       // O!
        &PyArray_Type, (PyObject **)&ap_ipntr,   // O!
        &PyArray_Type, (PyObject **)&ap_workd,   // O!
        &PyArray_Type, (PyObject **)&ap_workl    // O!
        )
    )
    {
        return NULL;
    }

    int* ipntr = PyArray_DATA(ap_ipntr);
    int* select = PyArray_DATA(ap_select);
    double* dr = PyArray_DATA(ap_dr);
    double* di = PyArray_DATA(ap_di);
    double* workev = PyArray_DATA(ap_workev);
    double* z = PyArray_DATA(ap_z);
    double* resid = PyArray_DATA(ap_resid);
    double* v = PyArray_DATA(ap_v);
    double* workd = PyArray_DATA(ap_workd);
    double* workl = PyArray_DATA(ap_workl);
    ldv = (int)PyArray_DIMS(ap_v)[0];
    ldz = (int)PyArray_DIMS(ap_z)[0];

    struct ARPACK_arnoldi_update_vars_d Vars = {0};

    #define X(name) Vars.name = 0;
    STRUCT_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = PyFloat_AsDouble(name##_obj);
        STRUCT_INEXACT_FIELD_NAMES
    #undef X

    #define X(name) \
        PyObject* name##_obj = PyDict_GetItemString(input_dict, #name); \
        if (!name##_obj) { PYERR(arpack_error, #name " not found in the dictionary."); } \
        Vars.name = (int)PyLong_AsLong(name##_obj);
        STRUCT_INT_FIELD_NAMES
    #undef X

    dneupd(&Vars, want_ev, howmny, select, dr, di, z, ldz, sigmar, sigmai, workev, resid, v, ldv, ipntr, workd, workl);

    Py_RETURN_NONE;

}


static PyObject*
cneupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    // This function is not implemented yet.
    PyErr_SetString(arpack_error, "cneupd_wrap is not implemented yet.");
    return NULL;
}


static PyObject*
zneupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    // This function is not implemented yet.
    PyErr_SetString(arpack_error, "zneupd_wrap is not implemented yet.");
    return NULL;
}


static PyObject*
ssaupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    // This function is not implemented yet.
    PyErr_SetString(arpack_error, "ssaupd_wrap is not implemented yet.");
    return NULL;
}

static PyObject*
dsaupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    // This function is not implemented yet.
    PyErr_SetString(arpack_error, "dsaupd_wrap is not implemented yet.");
    return NULL;
}

static PyObject*
sseupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    // This function is not implemented yet.
    PyErr_SetString(arpack_error, "sseupd_wrap is not implemented yet.");
    return NULL;
}


static PyObject*
dseupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    // This function is not implemented yet.
    PyErr_SetString(arpack_error, "dseupd_wrap is not implemented yet.");
    return NULL;
}


static char doc_snaupd[] = ("");
static char doc_dnaupd[] = ("");
static char doc_cnaupd[] = ("");
static char doc_znaupd[] = ("");
static char doc_sneupd[] = ("");
static char doc_dneupd[] = ("");
static char doc_cneupd[] = ("");
static char doc_zneupd[] = ("");
static char doc_ssaupd[] = ("");
static char doc_dsaupd[] = ("");
static char doc_sseupd[] = ("");
static char doc_dseupd[] = ("");


// Sentinel terminated method list.
static struct
PyMethodDef arpacklib_module_methods[] = {
    {"snaupd_wrap", snaupd_wrap, METH_VARARGS, doc_snaupd},
    {"dnaupd_wrap", dnaupd_wrap, METH_VARARGS, doc_dnaupd},
    {"cnaupd_wrap", cnaupd_wrap, METH_VARARGS, doc_cnaupd},
    {"znaupd_wrap", znaupd_wrap, METH_VARARGS, doc_znaupd},
    {"sneupd_wrap", sneupd_wrap, METH_VARARGS, doc_sneupd},
    {"dneupd_wrap", dneupd_wrap, METH_VARARGS, doc_dneupd},
    {"cneupd_wrap", cneupd_wrap, METH_VARARGS, doc_cneupd},
    {"zneupd_wrap", zneupd_wrap, METH_VARARGS, doc_zneupd},
    {"ssaupd_wrap", ssaupd_wrap, METH_VARARGS, doc_ssaupd},
    {"dsaupd_wrap", dsaupd_wrap, METH_VARARGS, doc_dsaupd},
    {"sseupd_wrap", sseupd_wrap, METH_VARARGS, doc_sseupd},
    {"dseupd_wrap", dseupd_wrap, METH_VARARGS, doc_dseupd},
    {NULL, NULL, 0, NULL}
};


static struct
PyModuleDef moduledef = {
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


#undef STRUCT_FIELD_NAMES
#undef STRUCT_INT_FIELD_NAMES
#undef STRUCT_INEXACT_FIELD_NAMES

#endif // ifndef __ARPACKMODULE_H
