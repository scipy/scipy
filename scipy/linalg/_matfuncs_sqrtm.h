#ifndef _MATFUNCS_SQRTM_H
#define _MATFUNCS_SQRTM_H

#include "_common_array_utils.h"

void matrix_squareroot_s(const PyArrayObject* ap_Am, const PyArrayObject* ap_ret, int* isIllconditioned, int* isSingular, int* sq_info, int* view_as_complex);
void matrix_squareroot_d(const PyArrayObject* ap_Am, const PyArrayObject* ap_ret, int* isIllconditioned, int* isSingular, int* sq_info, int* view_as_complex);
void matrix_squareroot_c(const PyArrayObject* ap_Am, const PyArrayObject* ap_ret, int* isIllconditioned, int* isSingular, int* sq_info, int* unused         );
void matrix_squareroot_z(const PyArrayObject* ap_Am, const PyArrayObject* ap_ret, int* isIllconditioned, int* isSingular, int* sq_info, int* unused         );


#define PYERR(errobj,message) {PyErr_SetString(errobj,message); return NULL;}

static PyObject* sqrtm_error;

/*
This function accepts an nD NumPy array of size (..., n, n) and, if possible,
computes the matrix square root for each (n, n) slice. The input must be at
least 2-dimensional.

Real valued input has the possibility of resulting in a complex-valued result
and this is tracked internally and hence depending on the resulting flag the
caller should view the resulting single/double (n, 2n) array as (n, n)
single/double complex array.

This function tests for the dtypes to be LAPACK compatible and for squareness.

This function does not test for 0-dimensional, or nD arrays where some dimensions
are 0. It is expected to have this case caught by the caller.
*/
static PyObject*
recursive_schur_sqrtm(PyObject *dummy, PyObject *args) {
    int isComplex = 0, isIllconditioned = 0, isSingular = 0, info;
    PyArrayObject *ap_Am=NULL;
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
    // Compare last two dimensions for squareness
    if (n != shape[ndim - 2])
    {
        PYERR(sqrtm_error, "Last two dimensions of the input must be the same.");
    }

    npy_intp ret_dims[1] = {1};
    for (int i = 0; i < ndim; i++) { ret_dims[0] *= shape[i]; }

    // Create the output array with the same shape as the input with twice the
    // number of entries to accomodate for potential complex-valued data from
    // real data which will be cast to complex e.g., "ret.asview(complex128)"
    // Example, (3, 3) -> (18), (4, 5, 5) -> (4, 50)
    if ((PyArray_TYPE(ap_Am) == NPY_FLOAT32) || (PyArray_TYPE(ap_Am) == NPY_FLOAT64))
    {
        ret_dims[0] *= 2;
    }

    // Create the return array
    PyArrayObject* ap_ret = (PyArrayObject*)PyArray_ZEROS(1, ret_dims, PyArray_TYPE(ap_Am), 0);

    switch (PyArray_TYPE(ap_Am))
    {
        case (NPY_FLOAT32):
        {
            matrix_squareroot_s(ap_Am, ap_ret, &isIllconditioned, &isSingular, &info, &isComplex);
            break;
        }
        case (NPY_FLOAT64):
        {
            matrix_squareroot_d(ap_Am, ap_ret, &isIllconditioned, &isSingular, &info, &isComplex);
            break;
        }
        case (NPY_COMPLEX64):
        {
            matrix_squareroot_c(ap_Am, ap_ret, &isIllconditioned, &isSingular, &info, &isComplex);
            isComplex = 1;
            break;
        }
        case (NPY_COMPLEX128):
        {
            matrix_squareroot_z(ap_Am, ap_ret, &isIllconditioned, &isSingular, &info, &isComplex);
            isComplex = 1;
            break;
        }
        default:
        {
            PYERR(sqrtm_error, "Unsupported data type.");
        }
    }

    // Return the result
    return Py_BuildValue("Niiii",PyArray_Return(ap_ret), isComplex, isIllconditioned, isSingular, info);
}

static struct PyMethodDef sqrtm_module_methods[] = {
  {"recursive_schur_sqrtm", recursive_schur_sqrtm, METH_VARARGS, "Compute the matrix square root by recursion."},
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


#undef SQRTM_C
#undef SQRTM_Z

#endif // _MATFUNCS_SQRTM_H
