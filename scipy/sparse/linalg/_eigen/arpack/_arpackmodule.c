/*
 * Python bindings for SciPy usage
 */

#include "stddef.h"

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include "ARPACK/_arpack.h"

#if defined(_MSC_VER)
    #define ARPACK_cplx(real, imag) ((_Dcomplex){real, imag})
    #define ARPACK_cplxf(real, imag) ((_Fcomplex){real, imag})
#else
    #define ARPACK_cplx(real, imag) ((real) + (imag)*I)
    #define ARPACK_cplxf(real, imag) ((real) + (imag)*I)
#endif

static PyObject* arpack_error;

typedef enum {
    FIELD_INT,
    FIELD_FLOATING,
} FieldType;

typedef struct {
    const char* name;
    size_t offset;
    FieldType type;
} GenericField;

// ARPACK has two structs that are part of the function signature and thus the API/ABI.
// They only differ in whether some floating fiealds are single or double precision, but
// because they're at the beginning of the struct (thus influencing the offsets of the
// fields for everything else), we cannot use a unified approach for both of them. see
// https://github.com/scipy/scipy/blob/main/scipy/sparse/linalg/_eigen/arpack/ARPACK/_arpack.h

// for offsetof, which needs a type
typedef struct ARPACK_arnoldi_update_vars_s ARPACK_arnoldi_update_vars_s;
typedef struct ARPACK_arnoldi_update_vars_d ARPACK_arnoldi_update_vars_d;

// to avoid duplication, use a macro so we define the fields only once; this needs to match the
// names (though not necessarily the order) of the fields in ARPACK_arnoldi_update_vars_{s,d}
#define DEFINE_ARPACK_FIELDS(suffix) \
static const GenericField arpack_fields_##suffix[] = { \
    {"tol",          offsetof(ARPACK_arnoldi_update_vars_##suffix, tol),           FIELD_FLOATING}, \
    {"getv0_rnorm0", offsetof(ARPACK_arnoldi_update_vars_##suffix, getv0_rnorm0),  FIELD_FLOATING}, \
    {"aitr_betaj",   offsetof(ARPACK_arnoldi_update_vars_##suffix, aitr_betaj),    FIELD_FLOATING}, \
    {"aitr_rnorm1",  offsetof(ARPACK_arnoldi_update_vars_##suffix, aitr_rnorm1),   FIELD_FLOATING}, \
    {"aitr_wnorm",   offsetof(ARPACK_arnoldi_update_vars_##suffix, aitr_wnorm),    FIELD_FLOATING}, \
    {"aup2_rnorm",   offsetof(ARPACK_arnoldi_update_vars_##suffix, aup2_rnorm),    FIELD_FLOATING}, \
    {"which",        offsetof(ARPACK_arnoldi_update_vars_##suffix, which),         FIELD_INT}, \
    {"ido",          offsetof(ARPACK_arnoldi_update_vars_##suffix, ido),           FIELD_INT}, \
    {"info",         offsetof(ARPACK_arnoldi_update_vars_##suffix, info),          FIELD_INT}, \
    {"bmat",         offsetof(ARPACK_arnoldi_update_vars_##suffix, bmat),          FIELD_INT}, \
    {"mode",         offsetof(ARPACK_arnoldi_update_vars_##suffix, mode),          FIELD_INT}, \
    {"n",            offsetof(ARPACK_arnoldi_update_vars_##suffix, n),             FIELD_INT}, \
    {"ncv",          offsetof(ARPACK_arnoldi_update_vars_##suffix, ncv),           FIELD_INT}, \
    {"nev",          offsetof(ARPACK_arnoldi_update_vars_##suffix, nev),           FIELD_INT}, \
    {"shift",        offsetof(ARPACK_arnoldi_update_vars_##suffix, shift),         FIELD_INT}, \
    {"maxiter",      offsetof(ARPACK_arnoldi_update_vars_##suffix, maxiter),       FIELD_INT}, \
    {"nconv",        offsetof(ARPACK_arnoldi_update_vars_##suffix, nconv),         FIELD_INT}, \
    {"iter",         offsetof(ARPACK_arnoldi_update_vars_##suffix, iter),          FIELD_INT}, \
    {"np",           offsetof(ARPACK_arnoldi_update_vars_##suffix, np),            FIELD_INT}, \
    {"getv0_first",  offsetof(ARPACK_arnoldi_update_vars_##suffix, getv0_first),   FIELD_INT}, \
    {"getv0_iter",   offsetof(ARPACK_arnoldi_update_vars_##suffix, getv0_iter),    FIELD_INT}, \
    {"getv0_itry",   offsetof(ARPACK_arnoldi_update_vars_##suffix, getv0_itry),    FIELD_INT}, \
    {"getv0_orth",   offsetof(ARPACK_arnoldi_update_vars_##suffix, getv0_orth),    FIELD_INT}, \
    {"aitr_iter",    offsetof(ARPACK_arnoldi_update_vars_##suffix, aitr_iter),     FIELD_INT}, \
    {"aitr_j",       offsetof(ARPACK_arnoldi_update_vars_##suffix, aitr_j),        FIELD_INT}, \
    {"aitr_orth1",   offsetof(ARPACK_arnoldi_update_vars_##suffix, aitr_orth1),    FIELD_INT}, \
    {"aitr_orth2",   offsetof(ARPACK_arnoldi_update_vars_##suffix, aitr_orth2),    FIELD_INT}, \
    {"aitr_restart", offsetof(ARPACK_arnoldi_update_vars_##suffix, aitr_restart),  FIELD_INT}, \
    {"aitr_step3",   offsetof(ARPACK_arnoldi_update_vars_##suffix, aitr_step3),    FIELD_INT}, \
    {"aitr_step4",   offsetof(ARPACK_arnoldi_update_vars_##suffix, aitr_step4),    FIELD_INT}, \
    {"aitr_ierr",    offsetof(ARPACK_arnoldi_update_vars_##suffix, aitr_ierr),     FIELD_INT}, \
    {"aup2_initv",   offsetof(ARPACK_arnoldi_update_vars_##suffix, aup2_initv),    FIELD_INT}, \
    {"aup2_iter",    offsetof(ARPACK_arnoldi_update_vars_##suffix, aup2_iter),     FIELD_INT}, \
    {"aup2_getv0",   offsetof(ARPACK_arnoldi_update_vars_##suffix, aup2_getv0),    FIELD_INT}, \
    {"aup2_cnorm",   offsetof(ARPACK_arnoldi_update_vars_##suffix, aup2_cnorm),    FIELD_INT}, \
    {"aup2_kplusp",  offsetof(ARPACK_arnoldi_update_vars_##suffix, aup2_kplusp),   FIELD_INT}, \
    {"aup2_nev0",    offsetof(ARPACK_arnoldi_update_vars_##suffix, aup2_nev0),     FIELD_INT}, \
    {"aup2_np0",     offsetof(ARPACK_arnoldi_update_vars_##suffix, aup2_np0),      FIELD_INT}, \
    {"aup2_numcnv",  offsetof(ARPACK_arnoldi_update_vars_##suffix, aup2_numcnv),   FIELD_INT}, \
    {"aup2_update",  offsetof(ARPACK_arnoldi_update_vars_##suffix, aup2_update),   FIELD_INT}, \
    {"aup2_ushift",  offsetof(ARPACK_arnoldi_update_vars_##suffix, aup2_ushift),   FIELD_INT}, \
};

// s = single precision; d = double precision
DEFINE_ARPACK_FIELDS(s)
DEFINE_ARPACK_FIELDS(d)
// by construction, same number of fields
#define ARPACK_FIELD_COUNT (sizeof(arpack_fields_s) / sizeof(arpack_fields_s[0]))


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

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_s[i].name;
        PyObject* obj = PyDict_GetItemString(input_dict, name);
        if (!obj) {
            PyErr_Format(arpack_error, "%s not found in the dictionary.", name);
        }

        // Use &Vars to get the struct's address; cast to char* so offset is in bytes,
        // not in units of the field's type (e.g. int or float).
        char* base = (char*)&Vars;
        void* ptr = base + arpack_fields_s[i].offset;

        switch (arpack_fields_s[i].type) {
            case FIELD_INT:
                *(int*)ptr = (int)PyLong_AsLong(obj);
                break;
            case FIELD_FLOATING:
                *(float*)ptr = (float)PyFloat_AsDouble(obj);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }
    }

    // Call ARPACK function
    ARPACK_snaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl);

    // Unpack the struct back to the dictionary
    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_s[i].name;
        // see above
        void* ptr = (char*)&Vars + arpack_fields_s[i].offset;

        PyObject* obj = NULL;
        switch (arpack_fields_s[i].type) {
            case FIELD_INT:
                obj = PyLong_FromLong((long)(*(int*)ptr));
                break;
            case FIELD_FLOATING:
                obj = PyFloat_FromDouble((double)(*(float*)ptr));
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }

        if (!obj || PyDict_SetItemString(input_dict, name, obj) < 0) {
            Py_XDECREF(obj);
            PyErr_Format(arpack_error, "Setting '%s' failed.", name);
        }
        Py_DECREF(obj);
    }

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

    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for floats and ints.

    // Declare and Initialize the ARPACK struct that will be populated from dict with zeros
    struct ARPACK_arnoldi_update_vars_d Vars = {0};

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_d[i].name;
        PyObject* obj = PyDict_GetItemString(input_dict, name);
        if (!obj) {
            PyErr_Format(arpack_error, "%s not found in the dictionary.", name);
        }
        void* ptr = (char*)&Vars + arpack_fields_d[i].offset;

        switch (arpack_fields_d[i].type) {
            case FIELD_INT:
                *(int*)ptr = (int)PyLong_AsLong(obj);
                break;
            case FIELD_FLOATING:
                *(double*)ptr = (double)PyFloat_AsDouble(obj);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }
    }

    // Call ARPACK function
    ARPACK_dnaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl);

    // Unpack the struct back to the dictionary
    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_d[i].name;
        void* ptr = (char*)&Vars + arpack_fields_d[i].offset;

        PyObject* obj = NULL;
        switch (arpack_fields_d[i].type) {
            case FIELD_INT:
                obj = PyLong_FromLong((long)(*(int*)ptr));
                break;
            case FIELD_FLOATING:
                obj = PyFloat_FromDouble(*(double*)ptr);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }

        if (!obj || PyDict_SetItemString(input_dict, name, obj) < 0) {
            Py_XDECREF(obj);
            PyErr_Format(arpack_error, "Setting '%s' failed.", name);
        }
        Py_DECREF(obj);
    }

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
    ARPACK_CPLXF_TYPE* resid = (ARPACK_CPLXF_TYPE*)PyArray_DATA(ap_resid);
    ARPACK_CPLXF_TYPE* v = (ARPACK_CPLXF_TYPE*)PyArray_DATA(ap_v);
    ARPACK_CPLXF_TYPE* workd = (ARPACK_CPLXF_TYPE*)PyArray_DATA(ap_workd);
    ARPACK_CPLXF_TYPE* workl = (ARPACK_CPLXF_TYPE*)PyArray_DATA(ap_workl);
    float* rwork = (float*)PyArray_DATA(ap_rwork);

    // Map the input dict to the ARPACK structure

    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for floats and ints.

    // Declare and Initialize the ARPACK struct that will be populated from dict with zeros
    struct ARPACK_arnoldi_update_vars_s Vars = {0};

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_s[i].name;
        PyObject* obj = PyDict_GetItemString(input_dict, name);
        if (!obj) {
            PyErr_Format(arpack_error, "%s not found in the dictionary.", name);
        }
        void* ptr = (char*)&Vars + arpack_fields_s[i].offset;

        switch (arpack_fields_s[i].type) {
            case FIELD_INT:
                *(int*)ptr = (int)PyLong_AsLong(obj);
                break;
            case FIELD_FLOATING:
                *(float*)ptr = (float)PyFloat_AsDouble(obj);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }
    }

    // Call ARPACK function
    ARPACK_cnaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl, rwork);

    // Unpack the struct back to the dictionary
    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_s[i].name;
        void* ptr = (char*)&Vars + arpack_fields_s[i].offset;

        PyObject* obj = NULL;
        switch (arpack_fields_s[i].type) {
            case FIELD_INT:
                obj = PyLong_FromLong((long)(*(int*)ptr));
                break;
            case FIELD_FLOATING:
                obj = PyFloat_FromDouble((double)(*(float*)ptr));
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }

        if (!obj || PyDict_SetItemString(input_dict, name, obj) < 0) {
            Py_XDECREF(obj);
            PyErr_Format(arpack_error, "Setting '%s' failed.", name);
        }
        Py_DECREF(obj);
    }

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
    ARPACK_CPLX_TYPE* resid = (ARPACK_CPLX_TYPE*)PyArray_DATA(ap_resid);
    ARPACK_CPLX_TYPE* v = (ARPACK_CPLX_TYPE*)PyArray_DATA(ap_v);
    ARPACK_CPLX_TYPE* workd = (ARPACK_CPLX_TYPE*)PyArray_DATA(ap_workd);
    ARPACK_CPLX_TYPE* workl = (ARPACK_CPLX_TYPE*)PyArray_DATA(ap_workl);
    double* rwork = (double*)PyArray_DATA(ap_rwork);

    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for floats and ints.

    // Declare and Initialize the ARPACK struct that will be populated from dict with zeros
    struct ARPACK_arnoldi_update_vars_d Vars = {0};

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_d[i].name;
        PyObject* obj = PyDict_GetItemString(input_dict, name);
        if (!obj) {
            PyErr_Format(arpack_error, "%s not found in the dictionary.", name);
        }
        void* ptr = (char*)&Vars + arpack_fields_d[i].offset;

        switch (arpack_fields_d[i].type) {
            case FIELD_INT:
                *(int*)ptr = (int)PyLong_AsLong(obj);
                break;
            case FIELD_FLOATING:
                *(double*)ptr = (double)PyFloat_AsDouble(obj);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }
    }

    // Call ARPACK function
    ARPACK_znaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl, rwork);

    // Unpack the struct back to the dictionary
    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_d[i].name;
        void* ptr = (char*)&Vars + arpack_fields_d[i].offset;

        PyObject* obj = NULL;
        switch (arpack_fields_d[i].type) {
            case FIELD_INT:
                obj = PyLong_FromLong((long)(*(int*)ptr));
                break;
            case FIELD_FLOATING:
                obj = PyFloat_FromDouble(*(double*)ptr);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }

        if (!obj || PyDict_SetItemString(input_dict, name, obj) < 0) {
            Py_XDECREF(obj);
            PyErr_Format(arpack_error, "Setting '%s' failed.", name);
        }
        Py_DECREF(obj);
    }

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

    struct ARPACK_arnoldi_update_vars_s Vars = {0};

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_s[i].name;
        PyObject* obj = PyDict_GetItemString(input_dict, name);
        if (!obj) {
            PyErr_Format(arpack_error, "%s not found in the dictionary.", name);
        }
        void* ptr = (char*)&Vars + arpack_fields_s[i].offset;

        switch (arpack_fields_s[i].type) {
            case FIELD_INT:
                *(int*)ptr = (int)PyLong_AsLong(obj);
                break;
            case FIELD_FLOATING:
                *(float*)ptr = (float)PyFloat_AsDouble(obj);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }
    }

    ARPACK_sneupd(&Vars, want_ev, howmny, select, dr, di, z, ldz, sigmar, sigmai, workev, resid, v, ldv, ipntr, workd, workl);

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

    struct ARPACK_arnoldi_update_vars_d Vars = {0};

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_d[i].name;
        PyObject* obj = PyDict_GetItemString(input_dict, name);
        if (!obj) {
            PyErr_Format(arpack_error, "%s not found in the dictionary.", name);
        }
        void* ptr = (char*)&Vars + arpack_fields_d[i].offset;

        switch (arpack_fields_d[i].type) {
            case FIELD_INT:
                *(int*)ptr = (int)PyLong_AsLong(obj);
                break;
            case FIELD_FLOATING:
                *(double*)ptr = (double)PyFloat_AsDouble(obj);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }
    }

    ARPACK_dneupd(&Vars, want_ev, howmny, select, dr, di, z, ldz, sigmar, sigmai, workev, resid, v, ldv, ipntr, workd, workl);

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
    ARPACK_CPLXF_TYPE* d = (ARPACK_CPLXF_TYPE*)PyArray_DATA(ap_d);
    ARPACK_CPLXF_TYPE* workev = (ARPACK_CPLXF_TYPE*)PyArray_DATA(ap_workev);
    ARPACK_CPLXF_TYPE* z = (ARPACK_CPLXF_TYPE*)PyArray_DATA(ap_z);
    ARPACK_CPLXF_TYPE* resid = (ARPACK_CPLXF_TYPE*)PyArray_DATA(ap_resid);
    ARPACK_CPLXF_TYPE* v = (ARPACK_CPLXF_TYPE*)PyArray_DATA(ap_v);
    ARPACK_CPLXF_TYPE* workd = (ARPACK_CPLXF_TYPE*)PyArray_DATA(ap_workd);
    ARPACK_CPLXF_TYPE* workl = (ARPACK_CPLXF_TYPE*)PyArray_DATA(ap_workl);
    float* rwork = PyArray_DATA(ap_rwork);
    ldv = (int)PyArray_DIMS(ap_v)[0];
    ldz = (int)PyArray_DIMS(ap_z)[0];
    ARPACK_CPLXF_TYPE sigmaC = ARPACK_cplxf((float)sigma.real, (float)sigma.imag);

    struct ARPACK_arnoldi_update_vars_s Vars = {0};

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_s[i].name;
        PyObject* obj = PyDict_GetItemString(input_dict, name);
        if (!obj) {
            PyErr_Format(arpack_error, "%s not found in the dictionary.", name);
        }
        void* ptr = (char*)&Vars + arpack_fields_s[i].offset;

        switch (arpack_fields_s[i].type) {
            case FIELD_INT:
                *(int*)ptr = (int)PyLong_AsLong(obj);
                break;
            case FIELD_FLOATING:
                *(float*)ptr = (float)PyFloat_AsDouble(obj);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }
    }

    ARPACK_cneupd(&Vars, want_ev, howmny, select, d, z, ldz, sigmaC, workev, resid, v, ldv, ipntr, workd, workl, rwork);

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
    ARPACK_CPLX_TYPE* d = (ARPACK_CPLX_TYPE*)PyArray_DATA(ap_d);
    ARPACK_CPLX_TYPE* workev = (ARPACK_CPLX_TYPE*)PyArray_DATA(ap_workev);
    ARPACK_CPLX_TYPE* z = (ARPACK_CPLX_TYPE*)PyArray_DATA(ap_z);
    ARPACK_CPLX_TYPE* resid = (ARPACK_CPLX_TYPE*)PyArray_DATA(ap_resid);
    ARPACK_CPLX_TYPE* v = (ARPACK_CPLX_TYPE*)PyArray_DATA(ap_v);
    ARPACK_CPLX_TYPE* workd = (ARPACK_CPLX_TYPE*)PyArray_DATA(ap_workd);
    ARPACK_CPLX_TYPE* workl = (ARPACK_CPLX_TYPE*)PyArray_DATA(ap_workl);
    double* rwork = PyArray_DATA(ap_rwork);
    ldv = (int)PyArray_DIMS(ap_v)[0];
    ldz = (int)PyArray_DIMS(ap_z)[0];
    ARPACK_CPLX_TYPE sigmaC = ARPACK_cplx(sigma.real, sigma.imag);

    struct ARPACK_arnoldi_update_vars_d Vars = {0};

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_d[i].name;
        PyObject* obj = PyDict_GetItemString(input_dict, name);
        if (!obj) {
            PyErr_Format(arpack_error, "%s not found in the dictionary.", name);
        }
        void* ptr = (char*)&Vars + arpack_fields_d[i].offset;

        switch (arpack_fields_d[i].type) {
            case FIELD_INT:
                *(int*)ptr = (int)PyLong_AsLong(obj);
                break;
            case FIELD_FLOATING:
                *(double*)ptr = (double)PyFloat_AsDouble(obj);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }
    }

    ARPACK_zneupd(&Vars, want_ev, howmny, select, d, z, ldz, sigmaC, workev, resid, v, ldv, ipntr, workd, workl, rwork);

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

    // Map the input dict to the ARPACK structure

    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for floats and ints.

    // Declare and Initialize the ARPACK struct that will be populated from dict with zeros
    struct ARPACK_arnoldi_update_vars_s Vars = {0};

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_s[i].name;
        PyObject* obj = PyDict_GetItemString(input_dict, name);
        if (!obj) {
            PyErr_Format(arpack_error, "%s not found in the dictionary.", name);
        }
        void* ptr = (char*)&Vars + arpack_fields_s[i].offset;

        switch (arpack_fields_s[i].type) {
            case FIELD_INT:
                *(int*)ptr = (int)PyLong_AsLong(obj);
                break;
            case FIELD_FLOATING:
                *(float*)ptr = (float)PyFloat_AsDouble(obj);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }
    }

    // Call ARPACK function
    ARPACK_ssaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl);

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_s[i].name;
        // see above
        void* ptr = (char*)&Vars + arpack_fields_s[i].offset;

        PyObject* obj = NULL;
        switch (arpack_fields_s[i].type) {
            case FIELD_INT:
                obj = PyLong_FromLong((long)(*(int*)ptr));
                break;
            case FIELD_FLOATING:
                obj = PyFloat_FromDouble((double)(*(float*)ptr));
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }

        if (!obj || PyDict_SetItemString(input_dict, name, obj) < 0) {
            Py_XDECREF(obj);
            PyErr_Format(arpack_error, "Setting '%s' failed.", name);
        }
        Py_DECREF(obj);
    }

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

    // Parse the dictionary, if the field is not found, raise an error.
    // Do it separately for floats and ints.

    // Declare and Initialize the ARPACK struct that will be populated from dict with zeros
    struct ARPACK_arnoldi_update_vars_d Vars = {0};

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_d[i].name;
        PyObject* obj = PyDict_GetItemString(input_dict, name);
        if (!obj) {
            PyErr_Format(arpack_error, "%s not found in the dictionary.", name);
        }
        void* ptr = (char*)&Vars + arpack_fields_d[i].offset;

        switch (arpack_fields_d[i].type) {
            case FIELD_INT:
                *(int*)ptr = (int)PyLong_AsLong(obj);
                break;
            case FIELD_FLOATING:
                *(double*)ptr = (double)PyFloat_AsDouble(obj);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }
    }

    // Call ARPACK function
    ARPACK_dsaupd(&Vars, resid, v, Vars.n, ipntr, workd, workl);

    // Unpack the struct back to the dictionary
    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_d[i].name;
        void* ptr = (char*)&Vars + arpack_fields_d[i].offset;

        PyObject* obj = NULL;
        switch (arpack_fields_d[i].type) {
            case FIELD_INT:
                obj = PyLong_FromLong((long)(*(int*)ptr));
                break;
            case FIELD_FLOATING:
                obj = PyFloat_FromDouble(*(double*)ptr);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }

        if (!obj || PyDict_SetItemString(input_dict, name, obj) < 0) {
            Py_XDECREF(obj);
            PyErr_Format(arpack_error, "Setting '%s' failed.", name);
        }
        Py_DECREF(obj);
    }

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

    struct ARPACK_arnoldi_update_vars_s Vars = {0};

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_s[i].name;
        PyObject* obj = PyDict_GetItemString(input_dict, name);
        if (!obj) {
            PyErr_Format(arpack_error, "%s not found in the dictionary.", name);
        }
        void* ptr = (char*)&Vars + arpack_fields_s[i].offset;

        switch (arpack_fields_s[i].type) {
            case FIELD_INT:
                *(int*)ptr = (int)PyLong_AsLong(obj);
                break;
            case FIELD_FLOATING:
                *(float*)ptr = (float)PyFloat_AsDouble(obj);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }
    }

    ARPACK_sseupd(&Vars, want_ev, howmny, select, d, z, ldz, sigma, resid, v, ldv, ipntr, workd, workl);

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

    struct ARPACK_arnoldi_update_vars_d Vars = {0};

    for (size_t i = 0; i < ARPACK_FIELD_COUNT; ++i) {
        const char* name = arpack_fields_d[i].name;
        PyObject* obj = PyDict_GetItemString(input_dict, name);
        if (!obj) {
            PyErr_Format(arpack_error, "%s not found in the dictionary.", name);
        }
        void* ptr = (char*)&Vars + arpack_fields_d[i].offset;

        switch (arpack_fields_d[i].type) {
            case FIELD_INT:
                *(int*)ptr = (int)PyLong_AsLong(obj);
                break;
            case FIELD_FLOATING:
                *(double*)ptr = (double)PyFloat_AsDouble(obj);
                break;
            default:
                PyErr_Format(arpack_error, "Unknown field type for %s", name);
        }
    }

    ARPACK_dseupd(&Vars, want_ev, howmny, select, d, z, ldz, sigma, resid, v, ldv, ipntr, workd, workl);

    Py_RETURN_NONE;
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
