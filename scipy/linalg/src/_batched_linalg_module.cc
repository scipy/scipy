#ifndef _LINALG_INV_H
#define _LINALG_INV_H

#include "_linalg_inv.hh"
#include "_linalg_solve.hh"
#include "_common_array_utils.hh"


#define PYERR(errobj,message) {PyErr_SetString(errobj,message); return NULL;}
static PyObject* _linalg_inv_error;


static PyObject*
_linalg_inv(PyObject* Py_UNUSED(dummy), PyObject* args) {

    PyArrayObject* ap_Am = NULL;
    PyArrayObject *ap_Ainv = NULL;
    CBLAS_INT info = 0;
    int isIllconditioned = 0;
    int isSingular = 0;
    St structure = St::NONE;
    int overwrite_a;

    // Get the input array
    if (!PyArg_ParseTuple(args, ("O!|np"), &PyArray_Type, (PyObject **)&ap_Am, &structure, &overwrite_a)) {
        return NULL;
    }

    // Check for dtype compatibility & array flags
    int typenum = PyArray_TYPE(ap_Am);
    bool dtype_ok = (typenum == NPY_FLOAT32)
                     || (typenum == NPY_FLOAT64)
                     || (typenum == NPY_COMPLEX64)
                     || (typenum == NPY_COMPLEX128);
    if(!dtype_ok || !PyArray_ISALIGNED(ap_Am)) {
        PyErr_SetString(PyExc_TypeError, "Expected a real or complex array.");
        return NULL;
    }

    int ndim = PyArray_NDIM(ap_Am);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_Am);      // Array shape
    npy_intp n = shape[ndim - 1];                // Slice size
    if (n != shape[ndim - 2]) {
        PYERR(PyExc_ValueError, "Last two dimensions of the input must be the same.")
    }

    overwrite_a = 0; // TODO: enable it

    if(!overwrite_a) {
        /* Allocate the output */
        ap_Ainv = (PyArrayObject *)PyArray_SimpleNew(ndim, shape, typenum);
        if(!ap_Ainv) {
            PyErr_NoMemory();
            return NULL;
        }
    }
    else {
        /* Reuse the memory buffer of the input array. */
        ap_Ainv = ap_Am;
        Py_INCREF(ap_Am);
    }

    void *buf = PyArray_DATA(ap_Ainv);
    switch(typenum) {
        case(NPY_FLOAT32):
            _inverse<float>(ap_Am, (float *)buf, structure, overwrite_a, &isIllconditioned, &isSingular, &info);
            break;
        case(NPY_FLOAT64):
            _inverse<double>(ap_Am, (double *)buf, structure, overwrite_a, &isIllconditioned, &isSingular, &info);
            break;
        case(NPY_COMPLEX64):
            _inverse<npy_complex64>(ap_Am, (npy_complex64 *)buf, structure, overwrite_a, &isIllconditioned, &isSingular, &info);
            break;
        case(NPY_COMPLEX128):
            _inverse<npy_complex128>(ap_Am, (npy_complex128 *)buf, structure, overwrite_a, &isIllconditioned, &isSingular, &info);
            break;
        default:
            PYERR(PyExc_RuntimeError, "Unknown array type.")
    }

    if(info < 0) {
        // Either OOM or internal LAPACK error.
        Py_DECREF(ap_Ainv);
        PYERR(PyExc_RuntimeError, "Internal LAPACK failure in scipy.linalg.inv.")
    }

    return Py_BuildValue("Niii", PyArray_Return(ap_Ainv), isIllconditioned, isSingular, info);
}


static PyObject*
_linalg_solve(PyObject* Py_UNUSED(dummy), PyObject* args) {

    PyArrayObject *ap_Am = NULL;
    PyArrayObject *ap_b = NULL;

    PyArrayObject *ap_x = NULL;
    int info = 0;
    int isIllconditioned = 0;
    int isSingular = 0;
    St structure = St::NONE;
    int overwrite_a = 0;
    int transposed = 0;
    int lower=0;

    // Get the input array
    if (!PyArg_ParseTuple(args, "O!O!|nppp", &PyArray_Type, (PyObject **)&ap_Am, &PyArray_Type, (PyObject **)&ap_b, &structure, &lower, &transposed, &overwrite_a)) {
        return NULL;
    }

    // Check for dtype compatibility & array flags
    int typenum = PyArray_TYPE(ap_Am);
    bool dtype_ok = (typenum == NPY_FLOAT32)
                     || (typenum == NPY_FLOAT64)
                     || (typenum == NPY_COMPLEX64)
                     || (typenum == NPY_COMPLEX128);
    if(!dtype_ok || !PyArray_ISALIGNED(ap_Am)) {
        PyErr_SetString(PyExc_TypeError, "Expected a real or complex array.");
        return NULL;
    }

    // Sanity check shapes
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp* shape = PyArray_SHAPE(ap_Am);
    if ((ndim < 2) || (shape[ndim - 1] != shape[ndim - 2])) {
        PYERR(PyExc_ValueError, "Last two dimensions of `a` must be the same.")
    }

    // At the python call site, 
    // 1) 1D `b` must have been converted in to 2D, and
    // 2) batch dimensions of `a` and `b` have been broadcast
    // Therefore, if `a.shape == (s, p, r, n, n)`, then `b.shape == (s, p, r, n, k)`
    // where `k` is the number of right-hand-sides.
    npy_intp ndim_b = PyArray_NDIM(ap_b);
    npy_intp *shape_b = PyArray_SHAPE(ap_b);

    bool dims_match = ndim_b == ndim;
    if (dims_match) {
        for (int i=0; i<ndim-1; i++) {
            dims_match = dims_match && (shape[i] == shape_b[i]);
        }
    }
    if (!dims_match){
        PyErr_SetString(PyExc_ValueError, "`a` and `b` shape mismatch.");
        return NULL;
    }

    // Allocate the output
    ap_x = (PyArrayObject *)PyArray_SimpleNew(ndim_b, shape_b, typenum);
    if(!ap_x) {
        PyErr_NoMemory();
        return NULL;
    }

    void *buf = PyArray_DATA(ap_x);
    switch(typenum) {
        case(NPY_FLOAT32):
            _solve<float>(ap_Am, ap_b, (float *)buf, structure, lower, transposed, overwrite_a, &isIllconditioned, &isSingular, &info);
            break;
        case(NPY_FLOAT64):
            _solve<double>(ap_Am, ap_b, (double *)buf, structure, lower, transposed, overwrite_a, &isIllconditioned, &isSingular, &info);
            break;
        case(NPY_COMPLEX64):
            _solve<npy_complex64>(ap_Am, ap_b, (npy_complex64 *)buf, structure, lower, transposed, overwrite_a, &isIllconditioned, &isSingular, &info);
            break;
        case(NPY_COMPLEX128):
            _solve<npy_complex128>(ap_Am, ap_b, (npy_complex128 *)buf, structure, lower, transposed, overwrite_a, &isIllconditioned, &isSingular, &info);
            break;
        default:
            PYERR(PyExc_RuntimeError, "Unknown array type.")
    }

    if(info < 0) {
        // Either OOM or internal LAPACK error.
        Py_DECREF(ap_x);
        PYERR(PyExc_RuntimeError, "Internal LAPACK failure in scipy.linalg.solve.")
    }


    return Py_BuildValue("Niii", PyArray_Return(ap_x), isIllconditioned, isSingular, info);
}


static char doc_inv[] = ("Compute the matrix inverse.");
static char doc_solve[] = ("Solve the linear system of equations.");

static struct PyMethodDef inv_module_methods[] = {
  {"_inv", _linalg_inv, METH_VARARGS, doc_inv},
  {"_solve", _linalg_solve, METH_VARARGS, doc_solve},
  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_batched_linalg",
    NULL,
    -1,
    inv_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__batched_linalg(void)
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
