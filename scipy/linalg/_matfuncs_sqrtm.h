#ifndef _MATFUNCS_SQRTM_H
#define _MATFUNCS_SQRTM_H

#include "_common_array_utils.h"

void matrix_squareroot_s(const PyArrayObject* ap_Am, float* restrict ap_ret, int* isIllconditioned, int* isSingular, int* sq_info, int* view_as_complex);
void matrix_squareroot_d(const PyArrayObject* ap_Am, double* restrict ap_ret, int* isIllconditioned, int* isSingular, int* sq_info, int* view_as_complex);
void matrix_squareroot_c(const PyArrayObject* ap_Am, SCIPY_C* restrict ap_ret, int* isIllconditioned, int* isSingular, int* sq_info, int* unused         );
void matrix_squareroot_z(const PyArrayObject* ap_Am, SCIPY_Z* restrict ap_ret, int* isIllconditioned, int* isSingular, int* sq_info, int* unused         );


#define PYERR(errobj,message) {PyErr_SetString(errobj,message); return NULL;}

static PyObject* sqrtm_error;

/*
 * This function accepts an nD NumPy array of size (..., n, n) and, if possible,
 * computes the matrix square root for each (n, n) slice. The input must be at
 * least 2-dimensional.
 *
 * Real valued input has the possibility of resulting in a complex-valued result
 * and this is tracked internally and hence depending on the resulting flag the
 * caller should view the resulting single/double (n, 2n) array as (n, n)
 * single/double complex array.
 *
 * This function tests for the dtypes to be LAPACK compatible and for squareness.
 *
 * IMPORTANT: This function does not test for 0-dimensional, or nD arrays where
 * some dimensions are 0. It is expected to have this case caught by the caller.
*/
static PyObject*
recursive_schur_sqrtm(PyObject *dummy, PyObject *args) {
    int isComplex = 0, isIllconditioned = 0, isSingular = 0, info = 0;
    PyArrayObject *ap_Am=NULL;
    void* mem_ret = NULL;
    // Get the input array
    if (!PyArg_ParseTuple(args, ("O!"), &PyArray_Type, (PyObject **)&ap_Am))
    {
        return NULL;
    }

    // Check for dtype compatibility
    if ((PyArray_TYPE(ap_Am) != NPY_FLOAT64) &&
         (PyArray_TYPE(ap_Am) != NPY_FLOAT32) &&
         (PyArray_TYPE(ap_Am) != NPY_COMPLEX64) &&
         (PyArray_TYPE(ap_Am) != NPY_COMPLEX128)
       )
    {
        PYERR(sqrtm_error, "Input must be nD array of type "
              "float32, float64, complex64, or complex128.");
    }

    int ndim = PyArray_NDIM(ap_Am);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_Am);      // Array shape
    npy_intp n = shape[ndim - 1];                // Slice size
    int input_type = PyArray_TYPE(ap_Am);        // Data type enum value
    // Compare last two dimensions for squareness
    if (n != shape[ndim - 2])
    {
        PYERR(sqrtm_error, "Last two dimensions of the input must be the same.");
    }

    // Create the output array with the same shape as the input with twice the
    // number of entries to accomodate for potential complex-valued data from
    // real data which will be cast to complex e.g., "ret.asview(complex128)"
    // Example, (3, 3) -> (18), (4, 5, 5) -> (4, 50)
    npy_intp ret_dims = 1;
    for (int i = 0; i < ndim; i++) { ret_dims *= shape[i]; }
    if ((input_type == NPY_FLOAT32) || (input_type == NPY_FLOAT64))
    {
        ret_dims *= 2;
    }

    if (PyArray_TYPE(ap_Am) == NPY_FLOAT32)
    {
        mem_ret = malloc(ret_dims*sizeof(float));
        if (mem_ret == NULL) { PYERR(sqrtm_error, "Memory allocation failed in scipy.linalg.sqrtm."); }
        matrix_squareroot_s(ap_Am, (float*)mem_ret, &isIllconditioned, &isSingular, &info, &isComplex);

    } else if (PyArray_TYPE(ap_Am) == NPY_FLOAT64) {

        mem_ret = malloc(ret_dims*sizeof(double));
        if (mem_ret == NULL) { PYERR(sqrtm_error, "Memory allocation failed in scipy.linalg.sqrtm."); }
        matrix_squareroot_d(ap_Am, (double*)mem_ret, &isIllconditioned, &isSingular, &info, &isComplex);

    } else if (PyArray_TYPE(ap_Am) == NPY_COMPLEX64) {

        mem_ret = malloc(ret_dims*sizeof(SCIPY_C));
        if (mem_ret == NULL) { PYERR(sqrtm_error, "Memory allocation failed in scipy.linalg.sqrtm."); }
        matrix_squareroot_c(ap_Am, (SCIPY_C*)mem_ret, &isIllconditioned, &isSingular, &info, &isComplex);
        isComplex = 1;

    } else if (PyArray_TYPE(ap_Am) == NPY_COMPLEX128) {

        mem_ret = malloc(ret_dims*sizeof(SCIPY_Z));
        if (mem_ret == NULL) { PYERR(sqrtm_error, "Memory allocation failed in scipy.linalg.sqrtm."); }
        matrix_squareroot_z(ap_Am, (SCIPY_Z*)mem_ret, &isIllconditioned, &isSingular, &info, &isComplex);
        isComplex = 1;

    } else {
        PYERR(sqrtm_error, "Unsupported input data type to scipy.linalg.sqrtm C function.");
    }

    if (info < 0)
    {
        // Internal failure memory or LAPACK error, fail and return the error code in info
        free(mem_ret);
        Py_INCREF(Py_None);
        return Py_BuildValue("Niii", Py_None, isIllconditioned, isSingular, info);
    }

    if (!isComplex)
    {
        // Input and output is real, truncate the extra space at the end
        npy_intp new_size = (ret_dims/2)*(input_type == NPY_FLOAT32 ? sizeof(float) : sizeof(double));
        void* mem_ret_half = realloc(mem_ret, new_size);
        // Quite unlikely but still allowed to fail
        if (!mem_ret_half) {
            free(mem_ret);
            PYERR(sqrtm_error, "Memory reallocation failed.");
        }
        PyArrayObject* ap_ret = (PyArrayObject*)PyArray_SimpleNewFromData(ndim, shape, input_type, mem_ret_half);
        return Py_BuildValue("Niii",PyArray_Return(ap_ret), isIllconditioned, isSingular, info);

    } else if ((input_type == NPY_FLOAT32) || (input_type == NPY_FLOAT64)) {

        // Input was real, result is complex, then view the result as complex
        int new_type = (PyArray_TYPE(ap_Am) == NPY_FLOAT32 ? NPY_COMPLEX64 : NPY_COMPLEX128);
        PyArrayObject* ap_ret = (PyArrayObject*)PyArray_SimpleNewFromData(ndim, shape, new_type, mem_ret);
        return Py_BuildValue("Niii",PyArray_Return(ap_ret), isIllconditioned, isSingular, info);

    } else {

        // Input was complex, result is complex, only reshape
        PyArrayObject* ap_ret = (PyArrayObject*)PyArray_SimpleNewFromData(ndim, shape, PyArray_TYPE(ap_Am), mem_ret);
        // Return the result
        return Py_BuildValue("Niii",PyArray_Return(ap_ret), isIllconditioned, isSingular, info);
    }
}

static char doc_sqrtm[] = ("Compute the matrix square root by recursion.\n\n    "
                           "sqrtmA, isIllConditioned, isSingular, info = recursive_schur_sqrtm(A)\n\n");

static struct PyMethodDef sqrtm_module_methods[] = {
  {"recursive_schur_sqrtm", recursive_schur_sqrtm, METH_VARARGS, doc_sqrtm},
  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_matfuncs_schur_sqrtm",
    NULL,
    -1,
    sqrtm_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__matfuncs_schur_sqrtm(void)
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
    sqrtm_error = PyErr_NewException("_matfuncs_schur_sqrtm.error", NULL, NULL);
    if (sqrtm_error == NULL) {
        return NULL;
    }
    if (PyDict_SetItemString(mdict, "error", sqrtm_error)) {
        return NULL;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}

#endif // _MATFUNCS_SQRTM_H
