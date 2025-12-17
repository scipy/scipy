#include "src/_common_array_utils.h"

void matrix_exponential_s(PyArrayObject* a, float* restrict result, int* info);
void matrix_exponential_d(PyArrayObject* a, double* restrict result, int* info);
void matrix_exponential_c(PyArrayObject* a, SCIPY_C* restrict result, int* info);
void matrix_exponential_z(PyArrayObject* a, SCIPY_Z* restrict result, int* info);

void matrix_squareroot_s(const PyArrayObject* ap_Am, float* restrict ap_ret, int* isIllconditioned, int* isSingular, int* sq_info, int* view_as_complex);
void matrix_squareroot_d(const PyArrayObject* ap_Am, double* restrict ap_ret, int* isIllconditioned, int* isSingular, int* sq_info, int* view_as_complex);
void matrix_squareroot_c(const PyArrayObject* ap_Am, SCIPY_C* restrict ap_ret, int* isIllconditioned, int* isSingular, int* sq_info, int* unused         );
void matrix_squareroot_z(const PyArrayObject* ap_Am, SCIPY_Z* restrict ap_ret, int* isIllconditioned, int* isSingular, int* sq_info, int* unused         );

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); return NULL;}

static PyObject* matfuncs_error;


// A simple destructor for buffer attached to a NumPy array via a capsule.
static void
capsule_destructor(PyObject *capsule) {
    void *ptr = PyCapsule_GetPointer(capsule, NULL);
    free(ptr);
}



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
recursive_schur_sqrtm(PyObject* Py_UNUSED(dummy), PyObject *args) {
    int isComplex = 0, isIllconditioned = 0, isSingular = 0, info = 0;
    PyArrayObject *ap_Am = NULL;
    void* mem_ret = NULL;
    PyArrayObject* ap_ret = NULL;

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
        PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Input must be an nD-array"
                          "  of type float32, float64, complex64, or complex128.");
    }

    int ndim = PyArray_NDIM(ap_Am);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_Am);      // Array shape
    npy_intp n = shape[ndim - 1];                // Slice size
    int input_type = PyArray_TYPE(ap_Am);        // Data type enum value
    // Compare last two dimensions for squareness
    if (n != shape[ndim - 2])
    {
        PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Last two dimensions of the input must be the same.");
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
        if (mem_ret == NULL) { PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Memory allocation failed."); }
        matrix_squareroot_s(ap_Am, (float*)mem_ret, &isIllconditioned, &isSingular, &info, &isComplex);

    } else if (PyArray_TYPE(ap_Am) == NPY_FLOAT64) {

        mem_ret = malloc(ret_dims*sizeof(double));
        if (mem_ret == NULL) { PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Memory allocation failed."); }
        matrix_squareroot_d(ap_Am, (double*)mem_ret, &isIllconditioned, &isSingular, &info, &isComplex);

    } else if (PyArray_TYPE(ap_Am) == NPY_COMPLEX64) {

        mem_ret = malloc(ret_dims*sizeof(SCIPY_C));
        if (mem_ret == NULL) { PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Memory allocation failed."); }
        matrix_squareroot_c(ap_Am, (SCIPY_C*)mem_ret, &isIllconditioned, &isSingular, &info, &isComplex);
        isComplex = 1;

    } else if (PyArray_TYPE(ap_Am) == NPY_COMPLEX128) {

        mem_ret = malloc(ret_dims*sizeof(SCIPY_Z));
        if (mem_ret == NULL) { PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Memory allocation failed."); }
        matrix_squareroot_z(ap_Am, (SCIPY_Z*)mem_ret, &isIllconditioned, &isSingular, &info, &isComplex);
        isComplex = 1;

    } else {
        PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Unsupported input data type to C function.");
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
        void* new_ret = realloc(mem_ret, new_size);
        // Quite unlikely but still allowed to fail
        if (!new_ret) {
            free(mem_ret);
            PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Memory reallocation failed.");
        }
        mem_ret = new_ret;
        ap_ret = (PyArrayObject*)PyArray_SimpleNewFromData(ndim, shape, input_type, mem_ret);
        if (ap_ret == NULL) {
            free(mem_ret);
            PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Failed to create numpy array from data.");
        }
    } else if ((input_type == NPY_FLOAT32) || (input_type == NPY_FLOAT64)) {
        // Input was real, result is complex, then view the result as complex
        int new_type = (PyArray_TYPE(ap_Am) == NPY_FLOAT32 ? NPY_COMPLEX64 : NPY_COMPLEX128);
        ap_ret = (PyArrayObject*)PyArray_SimpleNewFromData(ndim, shape, new_type, mem_ret);
        if (ap_ret == NULL) {
            free(mem_ret);
            PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Failed to create numpy array from data.");
        }
    } else {
        // Input was complex, result is complex, only reshape
        ap_ret = (PyArrayObject*)PyArray_SimpleNewFromData(ndim, shape, PyArray_TYPE(ap_Am), mem_ret);
        if (ap_ret == NULL) {
            free(mem_ret);
            PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Failed to create numpy array from data.");
        }
    }

    // Attach the buffer to a capsule
    PyObject* capsule = PyCapsule_New(mem_ret, NULL, capsule_destructor);
    if (capsule == NULL) {
        Py_DECREF(ap_ret);
        free(mem_ret);
        PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Failed to create capsule.");
    }
    // For ref counting of the memory
    if (PyArray_SetBaseObject(ap_ret, capsule) == -1) {
        Py_DECREF(ap_ret);
        Py_DECREF(capsule);
        PYERR(matfuncs_error, "scipy.linalg._matfuncs:sqrtm: Failed to set array's base.");
    }

    // Return the result
    return Py_BuildValue("Niii", PyArray_Return(ap_ret), isIllconditioned, isSingular, info);
}


static PyObject*
matrix_exponential(PyObject* Py_UNUSED(dummy), PyObject *args) {
    int info = 0;
    PyArrayObject *ap_a=NULL;
    void* mem_ret = NULL;

    if (!PyArg_ParseTuple(args, ("O!"), &PyArray_Type, (PyObject **)&ap_a)) {
         return NULL;
    }

    if (
        ((PyArray_TYPE(ap_a) != NPY_FLOAT64) &&
         (PyArray_TYPE(ap_a) != NPY_FLOAT32) &&
         (PyArray_TYPE(ap_a) != NPY_COMPLEX64) &&
         (PyArray_TYPE(ap_a) != NPY_COMPLEX128)) ||
        (PyArray_NDIM(ap_a) < 2))
    {
        PYERR(matfuncs_error, "scipy.linalg._matfuncs:expm: Input must be at least a 2D-array"
                          "  of type float32, float64, complex64, or complex128.");
    }
    int ndim = PyArray_NDIM(ap_a);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_a);      // Array shape
    int input_type = PyArray_TYPE(ap_a);        // Data type enum value

    if (shape[ndim - 1] != shape[ndim - 2])
    {
        PYERR(matfuncs_error, "scipy.linalg._matfuncs:expm: Last two dimensions of the input must be the same.");
    }

    npy_intp ret_dims = 1;
    for (int i = 0; i < ndim; i++) { ret_dims *= shape[i]; }

    if (PyArray_TYPE(ap_a) == NPY_FLOAT32)
    {
        mem_ret = calloc(ret_dims, sizeof(float));
        if (mem_ret == NULL) { PYERR(matfuncs_error, "scipy.linalg._matfuncs:expm: Memory allocation failed."); }
        matrix_exponential_s(ap_a, (float*)mem_ret, &info);

    } else if (PyArray_TYPE(ap_a) == NPY_FLOAT64) {

        mem_ret = calloc(ret_dims, sizeof(double));
        if (mem_ret == NULL) { PYERR(matfuncs_error, "scipy.linalg._matfuncs:expm: Memory allocation failed."); }
        matrix_exponential_d(ap_a, (double*)mem_ret, &info);

    } else if (PyArray_TYPE(ap_a) == NPY_COMPLEX64) {

        mem_ret = calloc(ret_dims, sizeof(SCIPY_C));
        if (mem_ret == NULL) { PYERR(matfuncs_error, "scipy.linalg._matfuncs:expm: Memory allocation failed."); }
        matrix_exponential_c(ap_a, (SCIPY_C*)mem_ret, &info);

    } else if (PyArray_TYPE(ap_a) == NPY_COMPLEX128) {

        mem_ret = calloc(ret_dims, sizeof(SCIPY_Z));
        if (mem_ret == NULL) { PYERR(matfuncs_error, "scipy.linalg._matfuncs:expm: Memory allocation failed."); }
        matrix_exponential_z(ap_a, (SCIPY_Z*)mem_ret, &info);
    } else {
        PYERR(matfuncs_error, "scipy.linalg._matfuncs:expm: Unsupported input data type to C function.");
    }

    if (info < 0)
    {
        // Internal failure memory or LAPACK error, fail and return the error code in info
        free(mem_ret);
        Py_INCREF(Py_None);
        return Py_BuildValue("Ni", Py_None, info);
    }

    PyArrayObject* ap_ret = (PyArrayObject*)PyArray_SimpleNewFromData(ndim, shape, input_type, mem_ret);
    if (ap_ret == NULL) {
        free(mem_ret);
        PYERR(matfuncs_error, "scipy.linalg._matfuncs:expm: Failed to create numpy array from data.");
    }

    // For ref counting of the memory
    PyObject* capsule = PyCapsule_New(mem_ret, NULL, capsule_destructor);
    if (capsule == NULL) {
        Py_DECREF(ap_ret);
        free(mem_ret);
        PYERR(matfuncs_error, "scipy.linalg._matfuncs:expm: Failed to create capsule.");
    }
    // For ref counting of the memory
    if (PyArray_SetBaseObject(ap_ret, capsule) == -1) {
        Py_DECREF(ap_ret);
        Py_DECREF(capsule);
        PYERR(matfuncs_error, "scipy.linalg._matfuncs:expm: Failed to set array's base.");
    }

    // Return the result
    return Py_BuildValue("Ni", PyArray_Return(ap_ret), info);
}


static char doc_sqrtm[] = ("Compute the matrix square root by recursion.\n\n    "
                           "sqrtmA, isIllConditioned, isSingular, info = recursive_schur_sqrtm(A)\n\n");

static char doc_expm[] = ("Compute the matrix exponential by scaling-squaring.\n\n    "
                           "expmA, info = matrix_exponential(A)\n\n");


static struct PyMethodDef matfuncs_module_methods[] = {
  {"recursive_schur_sqrtm", recursive_schur_sqrtm, METH_VARARGS, doc_sqrtm},
  {"matrix_exponential", matrix_exponential, METH_VARARGS, doc_expm},
  {NULL, NULL, 0, NULL}
};


static int matfuncs_module_exec(PyObject *module) {
    if (_import_array() < 0) { return -1; }

    matfuncs_error = PyErr_NewException("scipy.linalg._matfuncs.error", NULL, NULL);
    if (matfuncs_error == NULL) {
        return -1;
    }
    Py_INCREF(matfuncs_error);
    if (PyModule_AddObject(module, "error", matfuncs_error) < 0) {
        Py_DECREF(matfuncs_error);
        return -1;
    }

    return 0;
}


static struct PyModuleDef_Slot matfuncs_module_slots[] = {
    {Py_mod_exec, matfuncs_module_exec},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};


static struct PyModuleDef moduledef = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_internal_matfuncs",
    .m_doc = "Algorithms for various matrix valued functions for SciPy",
    .m_size = 0,
    .m_methods = matfuncs_module_methods,
    .m_slots = matfuncs_module_slots
};


PyMODINIT_FUNC PyInit__internal_matfuncs(void) {
    return PyModuleDef_Init(&moduledef);
}
