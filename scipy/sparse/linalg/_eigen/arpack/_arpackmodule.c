/*
 * Python bindings for SciPy usage
 */

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include "arnaud/include/arnaud/arnaud.h"

#if defined(_MSC_VER)
    #define ARNAUD_cplx(real, imag) ((_Dcomplex){real, imag})
    #define ARNAUD_cplxf(real, imag) ((_Fcomplex){real, imag})
#else
    #define ARNAUD_cplx(real, imag) ((real) + (imag)*I)
    #define ARNAUD_cplxf(real, imag) ((real) + (imag)*I)
#endif

static PyObject* arpack_error_obj;

// The following macros are used to define the field names in the ARPACK struct.
#define STRUCT_INEXACT_FIELD_NAMES X(tol) X(getv0_rnorm0) X(aitr_betaj) X(aitr_rnorm1) X(aitr_wnorm) X(aup2_rnorm)
#define STRUCT_INT_FIELD_NAMES X(ido) X(which) X(bmat) X(info) X(iter) X(maxiter) X(mode) \
                               X(n) X(nconv) X(ncv) X(nev) X(np) \
                               X(shift) X(getv0_first) X(getv0_iter) X(getv0_itry) X(getv0_orth) \
                               X(aitr_iter) X(aitr_j) X(aitr_orth1) X(aitr_orth2) X(aitr_restart) \
                               X(aitr_step3) X(aitr_step4) X(aitr_ierr) X(aup2_initv) X(aup2_iter) \
                               X(aup2_getv0) X(aup2_cnorm) X(aup2_kplusp) X(aup2_nev0) X(aup2_np0) \
                               X(aup2_numcnv) X(aup2_update) X(aup2_ushift)
#define STRUCT_FIELD_NAMES STRUCT_INT_FIELD_NAMES STRUCT_INEXACT_FIELD_NAMES

// Error reporting macro
#define ARPACK_ERROR(func_name, field_name, operation) \
    do { \
        char error_msg[256]; \
        snprintf(error_msg, sizeof(error_msg), "%s: %s field '%s'", func_name, operation, field_name); \
        PyErr_SetString(arpack_error_obj, error_msg); \
    } while(0)

    
static int
pack_dict_to_state_s(PyObject* dict, struct ARNAUD_state_s* vars, const char* func_name)
{
    if ((!dict) || (!vars) || (!func_name)) {
        PyErr_SetString(arpack_error_obj, "Internal error: NULL pointer in pack function");
        return -1;
    }

    // Pack float fields
    #define X(name) \
        do { \
            PyObject* field_obj = PyDict_GetItemString(dict, #name); \
            if (!field_obj) { ARPACK_ERROR(func_name, #name, "Missing required"); return -1; } \
            vars->name = (float)PyFloat_AsDouble(field_obj); \
            if (PyErr_Occurred()) { ARPACK_ERROR(func_name, #name, "Failed to convert"); return -1; } \
        } while(0);
    STRUCT_INEXACT_FIELD_NAMES
    #undef X

    // Pack integer fields
    #define X(name) \
        do { \
            PyObject* field_obj = PyDict_GetItemString(dict, #name); \
            if (!field_obj) { ARPACK_ERROR(func_name, #name, "Missing required");  return -1; } \
            vars->name = (int)PyLong_AsLong(field_obj); \
            if (PyErr_Occurred()) { ARPACK_ERROR(func_name, #name, "Failed to convert"); return -1; } \
        } while(0);
    STRUCT_INT_FIELD_NAMES
    #undef X

    return 0;
}


static int
pack_dict_to_state_d(PyObject* dict, struct ARNAUD_state_d* vars, const char* func_name)
{
    if ((!dict) || (!vars) || (!func_name)) {
        PyErr_SetString(arpack_error_obj, "Internal error: NULL pointer in pack function");
        return -1;
    }

    // Pack double fields
    #define X(name) \
        do { \
            PyObject* field_obj = PyDict_GetItemString(dict, #name); \
            if (!field_obj) { ARPACK_ERROR(func_name, #name, "Missing required"); return -1; } \
            vars->name = PyFloat_AsDouble(field_obj); \
            if (PyErr_Occurred()) { ARPACK_ERROR(func_name, #name, "Failed to convert");  return -1; } \
        } while(0);
    STRUCT_INEXACT_FIELD_NAMES
    #undef X

    // Pack integer fields
    #define X(name) \
        do { \
            PyObject* field_obj = PyDict_GetItemString(dict, #name); \
            if (!field_obj) { ARPACK_ERROR(func_name, #name, "Missing required"); return -1; } \
            vars->name = (int)PyLong_AsLong(field_obj); \
            if (PyErr_Occurred()) { ARPACK_ERROR(func_name, #name, "Failed to convert"); return -1; } \
        } while(0);
    STRUCT_INT_FIELD_NAMES
    #undef X

    return 0;
}


static int
unpack_state_s_to_dict(const struct ARNAUD_state_s* vars, PyObject* dict, const char* func_name)
{
    if ((!dict) || (!vars) || (!func_name)) {
        PyErr_SetString(arpack_error_obj, "Internal error: NULL pointer in unpack function");
        return -1;
    }

    // Unpack float fields
    #define X(name) \
        do { \
            PyObject* tmp_obj = PyFloat_FromDouble((double)vars->name); \
            if ((!tmp_obj) || (PyDict_SetItemString(dict, #name, tmp_obj) < 0)) { \
                Py_XDECREF(tmp_obj); \
                ARPACK_ERROR(func_name, #name, "Failed to set"); \
                return -1; \
            } \
            Py_DECREF(tmp_obj); \
        } while(0);
    STRUCT_INEXACT_FIELD_NAMES
    #undef X

    // Unpack integer fields
    #define X(name) \
        do { \
            PyObject* tmp_obj = PyLong_FromLong((long)vars->name); \
            if ((!tmp_obj) || (PyDict_SetItemString(dict, #name, tmp_obj) < 0)) { \
                Py_XDECREF(tmp_obj); \
                ARPACK_ERROR(func_name, #name, "Failed to set"); \
                return -1; \
            } \
            Py_DECREF(tmp_obj); \
        } while(0);
    STRUCT_INT_FIELD_NAMES
    #undef X

    return 0;
}


static int 
unpack_state_d_to_dict(const struct ARNAUD_state_d* vars, PyObject* dict, const char* func_name)
{
    if ((!dict) || (!vars) || (!func_name)) {
        PyErr_SetString(arpack_error_obj, "Internal error: NULL pointer in unpack function");
        return -1;
    }

    // Unpack inexact (floating-point) fields
    #define X(name) \
        do { \
            PyObject* tmp_obj = PyFloat_FromDouble(vars->name); \
            if ((!tmp_obj) || (PyDict_SetItemString(dict, #name, tmp_obj) < 0)) { \
                Py_XDECREF(tmp_obj); \
                ARPACK_ERROR(func_name, #name, "Failed to set"); \
                return -1; \
            } \
            Py_DECREF(tmp_obj); \
        } while(0);
    STRUCT_INEXACT_FIELD_NAMES
    #undef X

    // Unpack integer fields
    #define X(name) \
        do { \
            PyObject* tmp_obj = PyLong_FromLong((long)vars->name); \
            if ((!tmp_obj) || (PyDict_SetItemString(dict, #name, tmp_obj) < 0)) { \
                Py_XDECREF(tmp_obj); \
                ARPACK_ERROR(func_name, #name, "Failed to set"); \
                return -1; \
            } \
            Py_DECREF(tmp_obj); \
        } while(0);
    STRUCT_INT_FIELD_NAMES
    #undef X

    return 0;
}


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

    struct ARNAUD_state_s Vars = {0};
    if (pack_dict_to_state_s(input_dict, &Vars, "snaupd_wrap") != 0) { return NULL; }

    // Call ARPACK function
    ARNAUD_snaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl);

    // Unpack the struct back to the dictionary
    if (unpack_state_s_to_dict(&Vars, input_dict, "snaupd_wrap") != 0) { return NULL; }

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

    int* ipntr = (int*)PyArray_DATA(ap_ipntr);
    double* resid = (double*)PyArray_DATA(ap_resid);
    double* v = (double*)PyArray_DATA(ap_v);
    double* workd = (double*)PyArray_DATA(ap_workd);
    double* workl = (double*)PyArray_DATA(ap_workl);

    struct ARNAUD_state_d Vars = {0};
    if (pack_dict_to_state_d(input_dict, &Vars, "dnaupd_wrap") != 0) { return NULL; }

    // Call ARPACK function
    ARNAUD_dnaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl);

    // Unpack the struct back to the dictionary
    if (unpack_state_d_to_dict(&Vars, input_dict, "dnaupd_wrap") != 0) { return NULL; }

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

    int* ipntr = (int*)PyArray_DATA(ap_ipntr);
    ARNAUD_CPLXF_TYPE* resid = (ARNAUD_CPLXF_TYPE*)PyArray_DATA(ap_resid);
    ARNAUD_CPLXF_TYPE* v = (ARNAUD_CPLXF_TYPE*)PyArray_DATA(ap_v);
    ARNAUD_CPLXF_TYPE* workd = (ARNAUD_CPLXF_TYPE*)PyArray_DATA(ap_workd);
    ARNAUD_CPLXF_TYPE* workl = (ARNAUD_CPLXF_TYPE*)PyArray_DATA(ap_workl);
    float* rwork = (float*)PyArray_DATA(ap_rwork);

    struct ARNAUD_state_s Vars = {0};
    if (pack_dict_to_state_s(input_dict, &Vars, "cnaupd_wrap") != 0) { return NULL; }

    // Call ARPACK function
    ARNAUD_cnaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl, rwork);

    // Unpack the struct back to the dictionary
    if (unpack_state_s_to_dict(&Vars, input_dict, "cnaupd_wrap") != 0) { return NULL; }

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

    int* ipntr = (int*)PyArray_DATA(ap_ipntr);
    ARNAUD_CPLX_TYPE* resid = (ARNAUD_CPLX_TYPE*)PyArray_DATA(ap_resid);
    ARNAUD_CPLX_TYPE* v = (ARNAUD_CPLX_TYPE*)PyArray_DATA(ap_v);
    ARNAUD_CPLX_TYPE* workd = (ARNAUD_CPLX_TYPE*)PyArray_DATA(ap_workd);
    ARNAUD_CPLX_TYPE* workl = (ARNAUD_CPLX_TYPE*)PyArray_DATA(ap_workl);
    double* rwork = (double*)PyArray_DATA(ap_rwork);

    struct ARNAUD_state_d Vars = {0};
    if (pack_dict_to_state_d(input_dict, &Vars, "znaupd_wrap") != 0) { return NULL; }

    // Call ARPACK function
    ARNAUD_znaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl, rwork);

    // Unpack the struct back to the dictionary
    if (unpack_state_d_to_dict(&Vars, input_dict, "znaupd_wrap") != 0) { return NULL; }

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

    int* ipntr = (int*)PyArray_DATA(ap_ipntr);
    int* select = (int*)PyArray_DATA(ap_select);
    float* dr = (float*)PyArray_DATA(ap_dr);
    float* di = (float*)PyArray_DATA(ap_di);
    float* workev = (float*)PyArray_DATA(ap_workev);
    float* z = (float*)PyArray_DATA(ap_z);
    float* resid = (float*)PyArray_DATA(ap_resid);
    float* v = (float*)PyArray_DATA(ap_v);
    float* workd = (float*)PyArray_DATA(ap_workd);
    float* workl = (float*)PyArray_DATA(ap_workl);
    ldv = (int)PyArray_DIMS(ap_v)[0];
    ldz = (int)PyArray_DIMS(ap_z)[0];

    struct ARNAUD_state_s Vars = {0};
    if (pack_dict_to_state_s(input_dict, &Vars, "sneupd_wrap") != 0) { return NULL; }

    ARNAUD_sneupd(&Vars, want_ev, howmny, select, dr, di, z, ldz, sigmar, sigmai, workev, resid, v, ldv, ipntr, workd, workl);

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

    int* ipntr = (int*)PyArray_DATA(ap_ipntr);
    int* select = (int*)PyArray_DATA(ap_select);
    double* dr = (double*)PyArray_DATA(ap_dr);
    double* di = (double*)PyArray_DATA(ap_di);
    double* workev = (double*)PyArray_DATA(ap_workev);
    double* z = (double*)PyArray_DATA(ap_z);
    double* resid = (double*)PyArray_DATA(ap_resid);
    double* v = (double*)PyArray_DATA(ap_v);
    double* workd = (double*)PyArray_DATA(ap_workd);
    double* workl = (double*)PyArray_DATA(ap_workl);
    ldv = (int)PyArray_DIMS(ap_v)[0];
    ldz = (int)PyArray_DIMS(ap_z)[0];

    struct ARNAUD_state_d Vars = {0};
    if (pack_dict_to_state_d(input_dict, &Vars, "dneupd_wrap") != 0) { return NULL; }

    ARNAUD_dneupd(&Vars, want_ev, howmny, select, dr, di, z, ldz, sigmar, sigmai, workev, resid, v, ldv, ipntr, workd, workl);

    Py_RETURN_NONE;

}


static PyObject*
cneupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    PyObject* input_dict = NULL;
    int want_ev = 0, howmny = 0, ldv = 0, ldz = 0;
    PyArrayObject* ap_select = NULL;
    Py_complex sigma = { .real = 0.0, .imag = 0.0 };
    PyArrayObject* ap_d = NULL;
    PyArrayObject* ap_v = NULL;
    PyArrayObject* ap_z = NULL;
    PyArrayObject* ap_workev = NULL;
    PyArrayObject* ap_resid = NULL;
    PyArrayObject* ap_ipntr = NULL;
    PyArrayObject* ap_workd = NULL;
    PyArrayObject* ap_workl = NULL;
    PyArrayObject* ap_rwork = NULL;

    // Process input arguments
    if (!PyArg_ParseTuple(args, "O!iiO!O!O!DO!O!O!O!O!O!O!",
        &PyDict_Type, (PyObject **)&input_dict,  // O!
        &want_ev,                                // i
        &howmny,                                 // i
        &PyArray_Type, (PyObject **)&ap_select,  // O!
        &PyArray_Type, (PyObject **)&ap_d,       // O!
        &PyArray_Type, (PyObject **)&ap_z,       // O!
        &sigma,                                  // D
        &PyArray_Type, (PyObject **)&ap_workev,  // O!
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

    int* ipntr = (int*)PyArray_DATA(ap_ipntr);
    int* select = (int*)PyArray_DATA(ap_select);
    ARNAUD_CPLXF_TYPE* d = (ARNAUD_CPLXF_TYPE*)PyArray_DATA(ap_d);
    ARNAUD_CPLXF_TYPE* workev = (ARNAUD_CPLXF_TYPE*)PyArray_DATA(ap_workev);
    ARNAUD_CPLXF_TYPE* z = (ARNAUD_CPLXF_TYPE*)PyArray_DATA(ap_z);
    ARNAUD_CPLXF_TYPE* resid = (ARNAUD_CPLXF_TYPE*)PyArray_DATA(ap_resid);
    ARNAUD_CPLXF_TYPE* v = (ARNAUD_CPLXF_TYPE*)PyArray_DATA(ap_v);
    ARNAUD_CPLXF_TYPE* workd = (ARNAUD_CPLXF_TYPE*)PyArray_DATA(ap_workd);
    ARNAUD_CPLXF_TYPE* workl = (ARNAUD_CPLXF_TYPE*)PyArray_DATA(ap_workl);
    float* rwork = PyArray_DATA(ap_rwork);
    ldv = (int)PyArray_DIMS(ap_v)[0];
    ldz = (int)PyArray_DIMS(ap_z)[0];
    ARNAUD_CPLXF_TYPE sigmaC = ARNAUD_cplxf((float)sigma.real, (float)sigma.imag);

    struct ARNAUD_state_s Vars = {0};
    if (pack_dict_to_state_s(input_dict, &Vars, "cneupd_wrap") != 0) { return NULL; }

    ARNAUD_cneupd(&Vars, want_ev, howmny, select, d, z, ldz, sigmaC, workev, resid, v, ldv, ipntr, workd, workl, rwork);

    Py_RETURN_NONE;
}


static PyObject*
zneupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    PyObject* input_dict = NULL;
    int want_ev = 0, howmny = 0, ldv = 0, ldz = 0;
    PyArrayObject* ap_select = NULL;
    Py_complex sigma = { .real = 0.0, .imag = 0.0 };
    PyArrayObject* ap_d = NULL;
    PyArrayObject* ap_v = NULL;
    PyArrayObject* ap_z = NULL;
    PyArrayObject* ap_workev = NULL;
    PyArrayObject* ap_resid = NULL;
    PyArrayObject* ap_ipntr = NULL;
    PyArrayObject* ap_workd = NULL;
    PyArrayObject* ap_workl = NULL;
    PyArrayObject* ap_rwork = NULL;

    // Process input arguments
    if (!PyArg_ParseTuple(args, "O!iiO!O!O!DO!O!O!O!O!O!O!",
        &PyDict_Type, (PyObject **)&input_dict,  // O!
        &want_ev,                                // i
        &howmny,                                 // i
        &PyArray_Type, (PyObject **)&ap_select,  // O!
        &PyArray_Type, (PyObject **)&ap_d,       // O!
        &PyArray_Type, (PyObject **)&ap_z,       // O!
        &sigma,                                  // D
        &PyArray_Type, (PyObject **)&ap_workev,  // O!
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

    int* ipntr = (int*)PyArray_DATA(ap_ipntr);
    int* select = (int*)PyArray_DATA(ap_select);
    ARNAUD_CPLX_TYPE* d = (ARNAUD_CPLX_TYPE*)PyArray_DATA(ap_d);
    ARNAUD_CPLX_TYPE* workev = (ARNAUD_CPLX_TYPE*)PyArray_DATA(ap_workev);
    ARNAUD_CPLX_TYPE* z = (ARNAUD_CPLX_TYPE*)PyArray_DATA(ap_z);
    ARNAUD_CPLX_TYPE* resid = (ARNAUD_CPLX_TYPE*)PyArray_DATA(ap_resid);
    ARNAUD_CPLX_TYPE* v = (ARNAUD_CPLX_TYPE*)PyArray_DATA(ap_v);
    ARNAUD_CPLX_TYPE* workd = (ARNAUD_CPLX_TYPE*)PyArray_DATA(ap_workd);
    ARNAUD_CPLX_TYPE* workl = (ARNAUD_CPLX_TYPE*)PyArray_DATA(ap_workl);
    double* rwork = PyArray_DATA(ap_rwork);
    ldv = (int)PyArray_DIMS(ap_v)[0];
    ldz = (int)PyArray_DIMS(ap_z)[0];
    ARNAUD_CPLX_TYPE sigmaC = ARNAUD_cplx(sigma.real, sigma.imag);

    struct ARNAUD_state_d Vars = {0};
    if (pack_dict_to_state_d(input_dict, &Vars, "zneupd_wrap") != 0) { return NULL; }

    ARNAUD_zneupd(&Vars, want_ev, howmny, select, d, z, ldz, sigmaC, workev, resid, v, ldv, ipntr, workd, workl, rwork);

    Py_RETURN_NONE;
}


static PyObject*
ssaupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
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

    struct ARNAUD_state_s Vars = {0};
    if (pack_dict_to_state_s(input_dict, &Vars, "ssaupd_wrap") != 0) { return NULL; }

    // Call ARPACK function
    ARNAUD_ssaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl);

    // Unpack the struct back to the dictionary
    if (unpack_state_s_to_dict(&Vars, input_dict, "ssaupd_wrap") != 0) { return NULL; }

    Py_RETURN_NONE;
}

static PyObject*
dsaupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
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

    int* ipntr = (int*)PyArray_DATA(ap_ipntr);
    double* resid = (double*)PyArray_DATA(ap_resid);
    double* v = (double*)PyArray_DATA(ap_v);
    double* workd = (double*)PyArray_DATA(ap_workd);
    double* workl = (double*)PyArray_DATA(ap_workl);

    struct ARNAUD_state_d Vars = {0};
    if (pack_dict_to_state_d(input_dict, &Vars, "dsaupd_wrap") != 0) { return NULL; }

    // Call ARPACK function
    ARNAUD_dsaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl);

    // Unpack the struct back to the dictionary
    if (unpack_state_d_to_dict(&Vars, input_dict, "dsaupd_wrap") != 0) { return NULL; }

    Py_RETURN_NONE;
}

static PyObject*
sseupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    PyObject* input_dict = NULL;
    int want_ev = 0, howmny = 0, ldv = 0, ldz = 0;
    PyArrayObject* ap_select = NULL;
    float sigma = 0.0;
    PyArrayObject* ap_d = NULL;
    PyArrayObject* ap_v = NULL;
    PyArrayObject* ap_z = NULL;
    PyArrayObject* ap_resid = NULL;
    PyArrayObject* ap_ipntr = NULL;
    PyArrayObject* ap_workd = NULL;
    PyArrayObject* ap_workl = NULL;

    // Process input arguments
    if (!PyArg_ParseTuple(args, "O!iiO!O!O!fO!O!O!O!O!",
        &PyDict_Type, (PyObject **)&input_dict,  // O!
        &want_ev,                                // i
        &howmny,                                 // i
        &PyArray_Type, (PyObject **)&ap_select,  // O!
        &PyArray_Type, (PyObject **)&ap_d,       // O!
        &PyArray_Type, (PyObject **)&ap_z,       // O!
        &sigma,                                  // f
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

    int* ipntr = (int*)PyArray_DATA(ap_ipntr);
    int* select = (int*)PyArray_DATA(ap_select);
    float* d = (float*)PyArray_DATA(ap_d);
    float* z = (float*)PyArray_DATA(ap_z);
    float* resid = (float*)PyArray_DATA(ap_resid);
    float* v = (float*)PyArray_DATA(ap_v);
    float* workd = (float*)PyArray_DATA(ap_workd);
    float* workl = (float*)PyArray_DATA(ap_workl);
    ldv = (int)PyArray_DIMS(ap_v)[0];
    ldz = (int)PyArray_DIMS(ap_z)[0];

    struct ARNAUD_state_s Vars = {0};
    if (pack_dict_to_state_s(input_dict, &Vars, "sseupd_wrap") != 0) { return NULL; }

    ARNAUD_sseupd(&Vars, want_ev, howmny, select, d, z, ldz, sigma, resid, v, ldv, ipntr, workd, workl);

    Py_RETURN_NONE;
}


static PyObject*
dseupd_wrap(PyObject* Py_UNUSED(dummy), PyObject* args)
{
    PyObject* input_dict = NULL;
    int want_ev = 0, howmny = 0, ldv = 0, ldz = 0;
    PyArrayObject* ap_select = NULL;
    double sigma = 0.0;
    PyArrayObject* ap_d = NULL;
    PyArrayObject* ap_v = NULL;
    PyArrayObject* ap_z = NULL;
    PyArrayObject* ap_resid = NULL;
    PyArrayObject* ap_ipntr = NULL;
    PyArrayObject* ap_workd = NULL;
    PyArrayObject* ap_workl = NULL;

    // Process input arguments
    if (!PyArg_ParseTuple(args, "O!iiO!O!O!dO!O!O!O!O!",
        &PyDict_Type, (PyObject **)&input_dict,  // O!
        &want_ev,                                // i
        &howmny,                                 // i
        &PyArray_Type, (PyObject **)&ap_select,  // O!
        &PyArray_Type, (PyObject **)&ap_d,       // O!
        &PyArray_Type, (PyObject **)&ap_z,       // O!
        &sigma,                                  // d
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

    int* ipntr = (int*)PyArray_DATA(ap_ipntr);
    int* select = (int*)PyArray_DATA(ap_select);
    double* d = (double*)PyArray_DATA(ap_d);
    double* z = (double*)PyArray_DATA(ap_z);
    double* resid = (double*)PyArray_DATA(ap_resid);
    double* v = (double*)PyArray_DATA(ap_v);
    double* workd = (double*)PyArray_DATA(ap_workd);
    double* workl = (double*)PyArray_DATA(ap_workl);
    ldv = (int)PyArray_DIMS(ap_v)[0];
    ldz = (int)PyArray_DIMS(ap_z)[0];

    struct ARNAUD_state_d Vars = {0};
    if (pack_dict_to_state_d(input_dict, &Vars, "dseupd_wrap") != 0) { return NULL; }

    ARNAUD_dseupd(&Vars, want_ev, howmny, select, d, z, ldz, sigma, resid, v, ldv, ipntr, workd, workl);

    Py_RETURN_NONE;
}


static char doc_snaupd[] = (
    "snaupd_wrap(state, resid, v, ipntr, workd, workl)\n\n"
    "Internal wrapper for the ARPACK ``snaupd`` routine.\n\n"
    "This function performs the reverse communication steps for the Implicitly Restarted Arnoldi Method to find eigenvalues "
    "and eigenvectors of a non-symmetric matrix.\n\n"
    "Parameters\n----------\n"
    "state : dict\n    A dictionary holding the state of the ARPACK solver. It is used for both input and output, and is modified in-place.\n"
    "resid : ndarray, float32\n    The residual vector.\n"
    "v : ndarray, float32\n    The Arnoldi basis vectors (workspace).\n"
    "ipntr : ndarray, int32\n    Pointer array for ``workd`` and ``v``.\n"
    "workd : ndarray, float32\n    Distributed work array.\n"
    "workl : ndarray, float32\n    Private work array.\n\n"
    "Returns\n-------\n"
    "None\n    The ``state`` dictionary is updated in-place to communicate results and the next required action back to the caller."
);

static char doc_dnaupd[] = (
    "dnaupd_wrap(state, resid, v, ipntr, workd, workl)\n\n"
    "Internal wrapper for the ARPACK ``dnaupd`` routine (double precision).\n\n"
    "This function performs the reverse communication steps for the Implicitly Restarted Arnoldi Method to find eigenvalues "
    "and eigenvectors of a non-symmetric matrix.\n\n"
    "Parameters\n----------\n"
    "state : dict\n    A dictionary holding the state of the ARPACK solver. It is used for both input and output, and is modified in-place.\n"
    "resid : ndarray, float64\n    The residual vector.\n"
    "v : ndarray, float64\n    The Arnoldi basis vectors (workspace).\n"
    "ipntr : ndarray, int32\n    Pointer array for ``workd`` and ``v``.\n"
    "workd : ndarray, float64\n    Distributed work array.\n"
    "workl : ndarray, float64\n    Private work array.\n\n"
    "Returns\n-------\n"
    "None\n    The ``state`` dictionary is updated in-place to communicate results and the next required action back to the caller."
);

static char doc_cnaupd[] = (
    "cnaupd_wrap(state, resid, v, ipntr, workd, workl, rwork)\n\n"
    "Internal wrapper for the ARPACK ``cnaupd`` routine (single complex).\n\n"
    "This function performs the reverse communication steps for the Implicitly Restarted Arnoldi Method to find eigenvalues "
    "and eigenvectors of a non-symmetric matrix.\n\n"
    "Parameters\n----------\n"
    "state : dict\n    A dictionary holding the state of the ARPACK solver. It is used for both input and output, and is modified in-place.\n"
    "resid : ndarray, complex64\n    The residual vector.\n"
    "v : ndarray, complex64\n    The Arnoldi basis vectors (workspace).\n"
    "ipntr : ndarray, int32\n    Pointer array for ``workd`` and ``v``.\n"
    "workd : ndarray, complex64\n    Distributed work array.\n"
    "workl : ndarray, complex64\n    Private work array.\n"
    "rwork : ndarray, float32\n    Real workspace array.\n\n"
    "Returns\n-------\n"
    "None\n    The ``state`` dictionary is updated in-place to communicate results and the next required action back to the caller."
);

static char doc_znaupd[] = (
    "znaupd_wrap(state, resid, v, ipntr, workd, workl, rwork)\n\n"
    "Internal wrapper for the ARPACK ``znaupd`` routine (double complex).\n\n"
    "This function performs the reverse communication steps for the Implicitly Restarted Arnoldi Method to find eigenvalues "
    "and eigenvectors of a non-symmetric matrix.\n\n"
    "Parameters\n----------\n"
    "state : dict\n    A dictionary holding the state of the ARPACK solver. It is used for both input and output, and is modified in-place.\n"
    "resid : ndarray, complex128\n    The residual vector.\n"
    "v : ndarray, complex128\n    The Arnoldi basis vectors (workspace).\n"
    "ipntr : ndarray, int32\n    Pointer array for ``workd`` and ``v``.\n"
    "workd : ndarray, complex128\n    Distributed work array.\n"
    "workl : ndarray, complex128\n    Private work array.\n"
    "rwork : ndarray, float64\n    Real workspace array.\n\n"
    "Returns\n-------\n"
    "None\n    The ``state`` dictionary is updated in-place to communicate results and the next required action back to the caller."
);

static char doc_sneupd[] = (
    "sneupd_wrap(state, rvec, howmny, select, dr, di, z, sigmar, sigmai, ...)\n\n"
    "Internal wrapper for the ARPACK ``sneupd`` routine.\n\n"
    "This function is called after ``snaupd`` converges to compute the final eigenvalues and, optionally, the eigenvectors.\n\n"
    "Parameters\n----------\n"
    "state : dict\n    The final state dictionary from a converged ``snaupd`` run.\n"
    "rvec : bool\n    If true, compute and return the eigenvectors in ``z``.\n"
    "howmny : str\n    Specifies which eigenvectors to compute ('A' for all, 'S' for selected).\n"
    "select : ndarray, int32\n    A logical array specifying which Ritz values to target.\n"
    "dr : ndarray, float32\n    On return, contains the real parts of the computed eigenvalues.\n"
    "di : ndarray, float32\n    On return, contains the imaginary parts of the computed eigenvalues.\n"
    "z : ndarray, float32\n    On return, contains the computed eigenvectors if ``rvec`` is true.\n"
    "sigmar : float\n    Real part of the shift used in shift-invert mode.\n"
    "sigmai : float\n    Imaginary part of the shift used in shift-invert mode.\n"
    "... : \n    (other array arguments: workev, resid, v, ipntr, workd, workl)\n\n"
    "Returns\n-------\n"
    "None\n    The result arrays ``dr``, ``di``, and ``z`` are updated in-place."
);

static char doc_dneupd[] = (
    "dneupd_wrap(state, rvec, howmny, select, dr, di, z, sigmar, sigmai, ...)\n\n"
    "Internal wrapper for the ARPACK ``dneupd`` routine (double precision).\n\n"
    "This function is called after ``dnaupd`` converges to compute the final eigenvalues and, optionally, the eigenvectors.\n\n"
    "Parameters\n----------\n"
    "state : dict\n    The final state dictionary from a converged ``dnaupd`` run.\n"
    "rvec : bool\n    If true, compute and return the eigenvectors in ``z``.\n"
    "howmny : str\n    Specifies which eigenvectors to compute ('A' for all, 'S' for selected).\n"
    "select : ndarray, int32\n    A logical array specifying which Ritz values to target.\n"
    "dr : ndarray, float64\n    On return, contains the real parts of the computed eigenvalues.\n"
    "di : ndarray, float64\n    On return, contains the imaginary parts of the computed eigenvalues.\n"
    "z : ndarray, float64\n    On return, contains the computed eigenvectors if ``rvec`` is true.\n"
    "sigmar : float\n    Real part of the shift used in shift-invert mode.\n"
    "sigmai : float\n    Imaginary part of the shift used in shift-invert mode.\n"
    "... : \n    (other array arguments: workev, resid, v, ipntr, workd, workl)\n\n"
    "Returns\n-------\n"
    "None\n    The result arrays ``dr``, ``di``, and ``z`` are updated in-place."
);

static char doc_cneupd[] = (
    "cneupd_wrap(state, rvec, howmny, select, d, z, sigma, ...)\n\n"
    "Internal wrapper for the ARPACK ``cneupd`` routine (single complex).\n\n"
    "This function is called after ``cnaupd`` converges to compute the final eigenvalues and, optionally, the eigenvectors.\n\n"
    "Parameters\n----------\n"
    "state : dict\n    The final state dictionary from a converged ``cnaupd`` run.\n"
    "rvec : bool\n    If true, compute and return the eigenvectors in ``z``.\n"
    "howmny : str\n    Specifies which eigenvectors to compute ('A' for all, 'S' for selected).\n"
    "select : ndarray, int32\n    A logical array specifying which Ritz values to target.\n"
    "d : ndarray, complex64\n    On return, contains the computed complex eigenvalues.\n"
    "z : ndarray, complex64\n    On return, contains the computed eigenvectors if ``rvec`` is true.\n"
    "sigma : complex\n    The complex shift used in shift-invert mode.\n"
    "... : \n    (other array arguments: workev, resid, v, ipntr, workd, workl, rwork)\n\n"
    "Returns\n-------\n"
    "None\n    The result arrays ``d`` and ``z`` are updated in-place."
);

static char doc_zneupd[] = (
    "zneupd_wrap(state, rvec, howmny, select, d, z, sigma, ...)\n\n"
    "Internal wrapper for the ARPACK ``zneupd`` routine (double complex).\n\n"
    "This function is called after ``znaupd`` converges to compute the final eigenvalues and, optionally, the eigenvectors.\n\n"
    "Parameters\n----------\n"
    "state : dict\n    The final state dictionary from a converged ``znaupd`` run.\n"
    "rvec : bool\n    If true, compute and return the eigenvectors in ``z``.\n"
    "howmny : str\n    Specifies which eigenvectors to compute ('A' for all, 'S' for selected).\n"
    "select : ndarray, int32\n    A logical array specifying which Ritz values to target.\n"
    "d : ndarray, complex128\n    On return, contains the computed complex eigenvalues.\n"
    "z : ndarray, complex128\n    On return, contains the computed eigenvectors if ``rvec`` is true.\n"
    "sigma : complex\n    The complex shift used in shift-invert mode.\n"
    "... : \n    (other array arguments: workev, resid, v, ipntr, workd, workl, rwork)\n\n"
    "Returns\n-------\n"
    "None\n    The result arrays ``d`` and ``z`` are updated in-place."
);

static char doc_ssaupd[] = (
    "ssaupd_wrap(state, resid, v, ipntr, workd, workl)\n\n"
    "Internal wrapper for the ARPACK ``ssaupd`` routine.\n\n"
    "This function performs the reverse communication steps for the Implicitly Restarted Lanczos Method to find eigenvalues "
    "and eigenvectors of a symmetric matrix.\n\n"
    "Parameters\n----------\n"
    "state : dict\n    A dictionary holding the state of the ARPACK solver. It is used for both input and output, and is modified in-place.\n"
    "resid : ndarray, float32\n    The residual vector.\n"
    "v : ndarray, float32\n    The Lanczos basis vectors (workspace).\n"
    "ipntr : ndarray, int32\n    Pointer array for ``workd`` and ``v``.\n"
    "workd : ndarray, float32\n    Distributed work array.\n"
    "workl : ndarray, float32\n    Private work array.\n\n"
    "Returns\n-------\n"
    "None\n    The ``state`` dictionary is updated in-place to communicate results and the next required action back to the caller."
);

static char doc_dsaupd[] = (
    "dsaupd_wrap(state, resid, v, ipntr, workd, workl)\n\n"
    "Internal wrapper for the ARPACK ``dsaupd`` routine (double precision).\n\n"
    "This function performs the reverse communication steps for the Implicitly Restarted Lanczos Method to find eigenvalues "
    "and eigenvectors of a symmetric matrix.\n\n"
    "Parameters\n----------\n"
    "state : dict\n    A dictionary holding the state of the ARPACK solver. It is used for both input and output, and is modified in-place.\n"
    "resid : ndarray, float64\n    The residual vector.\n"
    "v : ndarray, float64\n    The Lanczos basis vectors (workspace).\n"
    "ipntr : ndarray, int32\n    Pointer array for ``workd`` and ``v``.\n"
    "workd : ndarray, float64\n    Distributed work array.\n"
    "workl : ndarray, float64\n    Private work array.\n\n"
    "Returns\n-------\n"
    "None\n    The ``state`` dictionary is updated in-place to communicate results and the next required action back to the caller."
);

static char doc_sseupd[] = (
    "sseupd_wrap(state, rvec, howmny, select, d, z, sigma, ...)\n\n"
    "Internal wrapper for the ARPACK ``sseupd`` routine.\n\n"
    "This function is called after ``ssaupd`` converges to compute the final eigenvalues and, optionally, the eigenvectors.\n\n"
    "Parameters\n----------\n"
    "state : dict\n    The final state dictionary from a converged ``ssaupd`` run.\n"
    "rvec : bool\n    If true, compute and return the eigenvectors in ``z``.\n"
    "howmny : str\n    Specifies which eigenvectors to compute ('A' for all, 'S' for selected).\n"
    "select : ndarray, int32\n    A logical array specifying which Ritz values to target.\n"
    "d : ndarray, float32\n    On return, contains the computed eigenvalues.\n"
    "z : ndarray, float32\n    On return, contains the computed eigenvectors if ``rvec`` is true.\n"
    "sigma : float\n    The shift used in shift-invert mode.\n"
    "... : \n    (other array arguments: resid, v, ipntr, workd, workl)\n\n"
    "Returns\n-------\n"
    "None\n    The result arrays ``d`` and ``z`` are updated in-place."
);

static char doc_dseupd[] = (
    "dseupd_wrap(state, rvec, howmny, select, d, z, sigma, ...)\n\n"
    "Internal wrapper for the ARPACK ``dseupd`` routine (double precision).\n\n"
    "This function is called after ``dsaupd`` converges to compute the final eigenvalues and, optionally, the eigenvectors.\n\n"
    "Parameters\n----------\n"
    "state : dict\n    The final state dictionary from a converged ``dsaupd`` run.\n"
    "rvec : bool\n    If true, compute and return the eigenvectors in ``z``.\n"
    "howmny : str\n    Specifies which eigenvectors to compute ('A' for all, 'S' for selected).\n"
    "select : ndarray, int32\n    A logical array specifying which Ritz values to target.\n"
    "d : ndarray, float64\n    On return, contains the computed eigenvalues.\n"
    "z : ndarray, float64\n    On return, contains the computed eigenvectors if ``rvec`` is true.\n"
    "sigma : float\n    The shift used in shift-invert mode.\n"
    "... : \n    (other array arguments: resid, v, ipntr, workd, workl)\n\n"
    "Returns\n-------\n"
    "None\n    The result arrays ``d`` and ``z`` are updated in-place."
);


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


static struct PyModuleDef_Slot arpacklib_module_slots[] = {
#if PY_VERSION_HEX >= 0x030c00f0  // Python 3.12+
    // signal that this module can be imported in isolated subinterpreters
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#endif
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    // signal that this module supports running without an active GIL
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};


static struct
PyModuleDef moduledef = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_arpacklib",
    .m_size = 0,
    .m_methods = arpacklib_module_methods,
    .m_slots = arpacklib_module_slots,
};


PyMODINIT_FUNC
PyInit__arpacklib(void)
{
    import_array();
    return PyModuleDef_Init(&moduledef);
}


#undef STRUCT_FIELD_NAMES
#undef STRUCT_INT_FIELD_NAMES
#undef STRUCT_INEXACT_FIELD_NAMES
