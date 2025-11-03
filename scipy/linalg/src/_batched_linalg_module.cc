#ifndef _LINALG_INV_H
#define _LINALG_INV_H

#include <cstring>

#include "_linalg_inv.hh"
#include "_linalg_solve.hh"
#include "_linalg_svd.hh"
#include "_common_array_utils.hh"


static PyObject* _linalg_inv_error;


PyObject* 
convert_vec_status(SliceStatusVec& vec_status);


static PyObject*
_linalg_inv(PyObject* Py_UNUSED(dummy), PyObject* args) {

    PyArrayObject* ap_Am = NULL;
    PyArrayObject *ap_Ainv = NULL;
    int info = 0;
    SliceStatusVec vec_status;
    St structure = St::NONE;
    int overwrite_a;
    int lower;

    // Get the input array
    if (!PyArg_ParseTuple(args, ("O!|npp"), &PyArray_Type, (PyObject **)&ap_Am, &structure, &overwrite_a, &lower)) {
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
        PyErr_SetString(PyExc_ValueError, "Last two dimensions of the input must be the same.");
        return NULL;
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
            info = _inverse<float>(ap_Am, (float *)buf, structure, lower, overwrite_a, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _inverse<double>(ap_Am, (double *)buf, structure, lower, overwrite_a, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _inverse<npy_complex64>(ap_Am, (npy_complex64 *)buf, structure, lower, overwrite_a, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _inverse<npy_complex128>(ap_Am, (npy_complex128 *)buf, structure, lower, overwrite_a, vec_status);
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unknown array type.");
            return NULL;
    }

    if (info < 0) {
        // Either OOM or internal LAPACK error.
        Py_DECREF(ap_Ainv);
        PyErr_SetString(PyExc_RuntimeError, "Memory error in scipy.linalg.inv.");
        return NULL;
    }
    PyObject *ret_lst = convert_vec_status(vec_status);

    return Py_BuildValue("NN", PyArray_Return(ap_Ainv), ret_lst);
}


static PyObject*
_linalg_solve(PyObject* Py_UNUSED(dummy), PyObject* args) {

    PyArrayObject *ap_Am = NULL;
    PyArrayObject *ap_b = NULL;

    PyArrayObject *ap_x = NULL;
    int info = 0;
    SliceStatusVec vec_status;
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
        PyErr_SetString(PyExc_ValueError, "Last two dimensions of `a` must be the same.");
        return NULL;
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
            info = _solve<float>(ap_Am, ap_b, (float *)buf, structure, lower, transposed, overwrite_a, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _solve<double>(ap_Am, ap_b, (double *)buf, structure, lower, transposed, overwrite_a, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _solve<npy_complex64>(ap_Am, ap_b, (npy_complex64 *)buf, structure, lower, transposed, overwrite_a, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _solve<npy_complex128>(ap_Am, ap_b, (npy_complex128 *)buf, structure, lower, transposed, overwrite_a, vec_status);
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unknown array type.");
            return NULL;
    }

    if (info < 0) {
        // Either OOM error or requiested lwork too large.
        Py_DECREF(ap_x);
        PyErr_SetString(PyExc_MemoryError, "Memory error in scipy.linalg.solve.");
        return NULL;
    }
    PyObject *ret_lst = convert_vec_status(vec_status);

    return Py_BuildValue("NN", PyArray_Return(ap_x), ret_lst);
}


static PyObject*
_linalg_svd(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_Am = NULL;
    const char *lapack_driver = NULL;
    int compute_uv = 1;
    int full_matrices = 1;
    PyArrayObject *ap_S = NULL, *ap_U = NULL, *ap_Vh = NULL;

    int info = 0;
    SliceStatusVec vec_status;

    // Get the input array
    if (!PyArg_ParseTuple(args, "O!s|pp", &PyArray_Type, (PyObject **)&ap_Am,  &lapack_driver, &compute_uv, &full_matrices)) {
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

    // Basic checks of array dimensions
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    if (ndim < 2) {
        PyErr_SetString(PyExc_ValueError, "Expected at least a 2D array.");
    }

    npy_intp m = shape[ndim - 2];
    npy_intp n = shape[ndim - 1];
    npy_intp k = m < n ? m : n; 

    // Allocate the output(s)

    // S.dtype is real if A.dtype is complex
    npy_intp typenum_S = typenum;
    if (typenum_S == NPY_COMPLEX64) { typenum_S = NPY_FLOAT32; }
    else if (typenum_S == NPY_COMPLEX128) { typenum_S = NPY_FLOAT64; }

    // S.shape = (..., k)
    npy_intp shape_1[NPY_MAXDIMS];
    for(int i=0; i<PyArray_NDIM(ap_Am); i++) {
        shape_1[i] = PyArray_DIM(ap_Am, i);
    }
    shape_1[ndim - 2] = k;
    ap_S = (PyArrayObject *)PyArray_SimpleNew(ndim-1, shape_1, typenum_S);
    if(!ap_S) {
        PyErr_NoMemory();
        return NULL;
    }

    char jobz = compute_uv ? (full_matrices ? 'A' : 'S') : 'N';

    // U.shape = (..., m, m) or (..., m, k) or not referenced (for jobz='N')
    // Vh.shape = (..., n, n) or (..., k, n) or not referenced (for jobz='N')
    if (jobz != 'N') {
        npy_intp u_shape0, u_shape1, vh_shape0, vh_shape1;
        u_vh_shapes(m, n, jobz, &u_shape0, &u_shape1, &vh_shape0, &vh_shape1);

        shape_1[ndim-2] = u_shape0;
        shape_1[ndim-1] = u_shape1;

        ap_U = (PyArrayObject *)PyArray_SimpleNew(ndim, shape_1, typenum);
        if (!ap_U) {
            PyErr_NoMemory();
            return NULL;
        }

        shape_1[ndim-2] = vh_shape0;
        shape_1[ndim-1] = vh_shape1;

        ap_Vh = (PyArrayObject *)PyArray_SimpleNew(ndim, shape_1, typenum);
        if (!ap_Vh) {
            PyErr_NoMemory();
            return NULL;
        }
    }

    switch(typenum) {
        case(NPY_FLOAT32):
            info = _svd<float>(ap_Am, ap_U, ap_S, ap_Vh, jobz, lapack_driver, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _svd<double>(ap_Am, ap_U, ap_S, ap_Vh, jobz, lapack_driver, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _svd<npy_complex64>(ap_Am, ap_U, ap_S, ap_Vh, jobz, lapack_driver, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _svd<npy_complex128>(ap_Am, ap_U, ap_S, ap_Vh, jobz, lapack_driver, vec_status);
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unknown array type.");
    }

    if (info < 0) {
        // Either OOM or internal LAPACK error.
        Py_DECREF(ap_S);
        Py_DECREF(ap_U);
        Py_DECREF(ap_Vh);
        PyErr_SetString(PyExc_RuntimeError, "Memory error in scipy.linalg.svd.");
        return NULL;
    }

    PyObject *ret_lst = convert_vec_status(vec_status);

    if (compute_uv){
        return Py_BuildValue("NNNN", PyArray_Return(ap_U), PyArray_Return(ap_S), PyArray_Return(ap_Vh), ret_lst);
    } else {
        return Py_BuildValue("NN", PyArray_Return(ap_S), ret_lst);
    }

}


/*
 * Helper: convert a vector of slice error statuses to list of dicts
 */
PyObject* 
convert_vec_status(SliceStatusVec& vec_status) {
    PyObject *ret_dct = NULL;
    PyObject *ret_lst = NULL;

    if (vec_status.empty()) {
        ret_lst = PyList_New(0);
    } else {
        // Problems detected in some slices, report.

        ret_lst = PyList_New(0);
        for (size_t i=0; i<vec_status.size(); i++) {
            SliceStatus status = vec_status[i];
            ret_dct = Py_BuildValue(
                "{s:n,s:n,s:i,s:i,s:d,s:n}",
                "num", status.slice_num,
                "structure", status.structure,
                "is_singular", status.is_singular,
                "is_ill_conditioned", status.is_ill_conditioned,
                "rcond", status.rcond,
                "lapack_info", status.lapack_info
            );
            PyList_Append(ret_lst, ret_dct);
        }
    }
    return ret_lst;
}


static char doc_inv[] = ("Compute the matrix inverse.");
static char doc_solve[] = ("Solve the linear system of equations.");
static char doc_svd[] = ("SVD factorization.");

static struct PyMethodDef inv_module_methods[] = {
  {"_inv", _linalg_inv, METH_VARARGS, doc_inv},
  {"_solve", _linalg_solve, METH_VARARGS, doc_solve},
  {"_svd", _linalg_svd, METH_VARARGS, doc_svd},
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
