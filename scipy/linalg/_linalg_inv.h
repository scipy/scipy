#ifndef _LINALG_INV_H
#define _LINALG_INV_H

#include "_common_array_utils.h"

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); return NULL;}
static PyObject* _linalg_inv_error;

void _inverse_s(const PyArrayObject* ap_Am, float* restrict ret_data, int* isIllconditioned, int* isSingular, int* info);
void _inverse_d(const PyArrayObject* ap_Am, double* restrict ret_data, int* isIllconditioned, int* isSingular, int* info);
void _inverse_c(const PyArrayObject* ap_Am, SCIPY_C* restrict ret_data, int* isIllconditioned, int* isSingular, int* info);
void _inverse_z(const PyArrayObject* ap_Am, SCIPY_Z* restrict ret_data, int* isIllconditioned, int* isSingular, int* info);


static PyObject*
_linalg_inv(PyObject* Py_UNUSED(dummy), PyObject* args) {

    PyArrayObject* ap_Am=NULL;
    PyArrayObject* ap_ret=NULL;
    void* mem_ret=NULL;
    void* work = NULL;
    int info = 0;
    int isIllconditioned = 0;
    int isSingular = 0;
    int info = 0;

    // Get the input array
    if (!PyArg_ParseTuple(args, ("O!"), &PyArray_Type, (PyObject **)&ap_Am)) {
        return NULL;
    }

    // Check for dtype compatibility
    if ((PyArray_TYPE(ap_Am) != NPY_FLOAT64) &&
         (PyArray_TYPE(ap_Am) != NPY_FLOAT32) &&
         (PyArray_TYPE(ap_Am) != NPY_COMPLEX64) &&
         (PyArray_TYPE(ap_Am) != NPY_COMPLEX128)
       ) {
        PYERR(_linalg_inv_error, "Input must be nD array of type "
              "float32, float64, complex64, or complex128.");
    }

    int ndim = PyArray_NDIM(ap_Am);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_Am);      // Array shape
    npy_intp n = shape[ndim - 1];                // Slice size
    if (n != shape[ndim - 2]) { PYERR(_linalg_inv_error, "Last two dimensions of the input must be the same."); }
    int input_type = PyArray_TYPE(ap_Am);        // Data type enum value
    npy_intp ret_dims = 1;
    for (int i = 0; i < ndim; i++) { ret_dims *= shape[i]; }  // Total size of the input array

    // Dispatch to the appropriate function based on the input type
    if (PyArray_TYPE(ap_Am) == NPY_FLOAT32)
    {
        mem_ret = malloc(ret_dims*sizeof(float));
        if (mem_ret == NULL) { PYERR(_linalg_inv_error, "Memory allocation failed in scipy.linalg.inv."); }
        _inverse_s(ap_Am, (float*)mem_ret, &isIllconditioned, &isSingular, &info);

    } else if (PyArray_TYPE(ap_Am) == NPY_FLOAT64) {
        mem_ret = malloc(ret_dims*sizeof(double));
        if (mem_ret == NULL) { PYERR(_linalg_inv_error, "Memory allocation failed in scipy.linalg.inv."); }
        _inverse_d(ap_Am, (double*)mem_ret, &isIllconditioned, &isSingular, &info);

    } else if (PyArray_TYPE(ap_Am) == NPY_COMPLEX64) {
        mem_ret = malloc(ret_dims*sizeof(SCIPY_C));
        if (mem_ret == NULL) { PYERR(_linalg_inv_error, "Memory allocation failed in scipy.linalg.inv."); }
        _inverse_c(ap_Am, (SCIPY_C*)mem_ret, &isIllconditioned, &isSingular, &info);

    } else if (PyArray_TYPE(ap_Am) == NPY_COMPLEX128) {
        mem_ret = malloc(ret_dims*sizeof(SCIPY_Z));
        if (mem_ret == NULL) { PYERR(_linalg_inv_error, "Memory allocation failed in scipy.linalg.inv."); }
        _inverse_z(ap_Am, (SCIPY_Z*)mem_ret, &isIllconditioned, &isSingular, &info);

    } else {
        PYERR(_linalg_inv_error, "Input must be nD array of type "
              "float32, float64, complex64, or complex128.");
    }

    // Return the result

    if (info < 0) {
        // Internal failure memory or LAPACK error, fail and return the error code in info
        free(mem_ret);
        PYERR(_linalg_inv_error, "Internal LAPACK failure in scipy.linalg.inv.");
    }
    ap_ret = (PyArrayObject*)PyArray_SimpleNewFromData(ndim, shape, PyArray_TYPE(ap_Am), mem_ret);
    return Py_BuildValue("Niii", PyArray_Return(ap_ret), isIllconditioned, isSingular, info);

}

static char doc_inv[] = ("Compute the matrix inverse.");

static struct PyMethodDef inv_module_methods[] = {
  {"_linalg_inv", _linalg_inv, METH_VARARGS, doc_inv},
  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_linalg_inv",
    NULL,
    -1,
    inv_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__linalg_inv(void)
{
    PyObject *module, *mdict;

    import_array();

    module = PyModule_Create(&moduledef);
    if (module == NULL) {
        return NULL;
    }

    mdict = PyModule_GetDict(module);
    if (mdict == NULL) {
        return NULL;
    }
    _linalg_inv_error = PyErr_NewException("_linalg_inv.error", NULL, NULL);
    if (_linalg_inv_error == NULL) {
        return NULL;
    }
    if (PyDict_SetItemString(mdict, "error", _linalg_inv_error)) {
        return NULL;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}



#endif // _LINALG_INV_H
