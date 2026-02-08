#include <cstring>
#include "_linalg_inv.hh"
#include "_linalg_solve.hh"
#include "_linalg_svd.hh"
#include "_linalg_lstsq.hh"
#include "_linalg_eig.hh"
#include "_common_array_utils.hh"


static PyObject* _linalg_inv_error;


PyObject*
convert_vec_status(SliceStatusVec& vec_status);

std::string
get_err_mesg(const std::string routine, const std::string func_name, int info);


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
    int overwrite_b = 0;
    int transposed = 0;
    int lower=0;

    // Get the input array
    if (!PyArg_ParseTuple(args, "O!O!|npppp", &PyArray_Type, (PyObject **)&ap_Am, &PyArray_Type, (PyObject **)&ap_b, &structure, &lower, &transposed, &overwrite_a, &overwrite_b)) {
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
_linalg_solve_banded(PyObject* Py_UNUSED(dummy), PyObject* args) {

    PyArrayObject *ap_Ab = NULL;
    PyArrayObject *ap_b = NULL;
    PyArrayObject *ap_kls = NULL; // lower bands
    PyArrayObject *ap_kus = NULL; // upper bands
    PyArrayObject *ap_x = NULL; // return object

    int overwrite_ab = 0;
    int overwrite_b = 0;
    int info;
    SliceStatusVec vec_status;

    // Get input data
    if (!PyArg_ParseTuple(args, "O!O!O!O!|pp", &PyArray_Type, (PyObject **)&ap_Ab, &PyArray_Type, (PyObject **)&ap_b, &PyArray_Type, (PyObject **)&ap_kls, &PyArray_Type, (PyObject **)&ap_kus, &overwrite_ab, &overwrite_b)) {
        PyErr_SetString(PyExc_ValueError, "Could not parse input.");
        return NULL;
    }

    // Check for dtype compatibility & array flags
    int typenum = PyArray_TYPE(ap_Ab);
    bool dtype_ok = (typenum == NPY_FLOAT32)
                     || (typenum == NPY_FLOAT64)
                     || (typenum == NPY_COMPLEX64)
                     || (typenum == NPY_COMPLEX128);
    if(!dtype_ok || !PyArray_ISALIGNED(ap_Ab)) {
        PyErr_SetString(PyExc_TypeError, "Expected a real or complex array.");
        return NULL;
    }

    if (!(PyArray_TYPE(ap_kls) == NPY_INT16) || !PyArray_ISALIGNED(ap_kls)) {
        PyErr_SetString(PyExc_TypeError, "Expected bounds to be integers.");
        return NULL;
    }

    // Sanity check shapes
    int ndim = PyArray_NDIM(ap_Ab);
    npy_intp* shape = PyArray_SHAPE(ap_Ab);
    if (ndim < 2) { // Can not perform a check on the two last dimensions: should be related to the banding
        PyErr_SetString(PyExc_ValueError, "Incorrect dimensions for ab.");
        return NULL;
    }

    // Allocate output object
    npy_intp ndim_b = PyArray_NDIM(ap_b);
    npy_intp *shape_b = PyArray_SHAPE(ap_b);

    bool dims_match = ndim_b == ndim;
    if (dims_match) {
        for (int i=0; i<ndim-2; i++) { // ndim - 2 can be different due to banded structure
            dims_match = dims_match && (shape[i] == shape_b[i]);
        }
    }
    if (!dims_match){
        PyErr_SetString(PyExc_ValueError, "`ab` and `b` shape mismatch.");
        return NULL;
    }

    // Allocate the output
    ap_x = (PyArrayObject *)PyArray_SimpleNew(ndim_b, shape_b, typenum);
    if(!ap_x) {
        PyErr_NoMemory();
        return NULL;
    }

    void *buf = PyArray_DATA(ap_x);
    switch (typenum) {
        case (NPY_FLOAT32):
            info = _solve_banded<float>(ap_Ab, ap_b, (float *)buf, ap_kls, ap_kus, overwrite_ab, overwrite_b, vec_status);
            break;

        case (NPY_FLOAT64):
            info = _solve_banded<double>(ap_Ab, ap_b, (double *)buf, ap_kls, ap_kus, overwrite_ab, overwrite_b, vec_status);
            break;

        case (NPY_COMPLEX64):
            info = _solve_banded<npy_complex64>(ap_Ab, ap_b, (npy_complex64 *)buf, ap_kls, ap_kus, overwrite_ab, overwrite_b, vec_status);
            break;

        case (NPY_COMPLEX128):
            info = _solve_banded<npy_complex128>(ap_Ab, ap_b, (npy_complex128 *)buf, ap_kls, ap_kus, overwrite_ab, overwrite_b, vec_status);
            break;

        default:
                PyErr_SetString(PyExc_RuntimeError, "Unknown array type.");
                return NULL;
    }

    if (info < 0) {
        // Either OOM error or requiested lwork too large.
        Py_DECREF(ap_x);
        PyErr_SetString(PyExc_MemoryError, "Memory error in scipy.linalg.solve_banded.");
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
    PyObject *ret_lst;

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
        return NULL;
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
            goto fail;
        }

        shape_1[ndim-2] = vh_shape0;
        shape_1[ndim-1] = vh_shape1;

        ap_Vh = (PyArrayObject *)PyArray_SimpleNew(ndim, shape_1, typenum);
        if (!ap_Vh) {
            PyErr_NoMemory();
            goto fail;
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
            goto fail;
    }

    if (info < 0) {
        // Either OOM or internal LAPACK error.
        PyErr_SetString(PyExc_RuntimeError, "Memory error in scipy.linalg.svd.");
        goto fail;
    }

    ret_lst = convert_vec_status(vec_status);

    if (compute_uv){
        return Py_BuildValue("NNNN", PyArray_Return(ap_U), PyArray_Return(ap_S), PyArray_Return(ap_Vh), ret_lst);
    } else {
        return Py_BuildValue("NN", PyArray_Return(ap_S), ret_lst);
    }

fail:
    Py_DECREF(ap_S);
    Py_XDECREF(ap_U);
    Py_XDECREF(ap_Vh);
    return NULL;
}



static PyObject*
_linalg_lstsq(PyObject* Py_UNUSED(dummy), PyObject* args) {

    PyArrayObject *ap_Am = NULL;
    PyArrayObject *ap_b = NULL;
    PyArrayObject *ap_S = NULL;
    PyArrayObject *ap_x = NULL;
    PyArrayObject *ap_rank = NULL;
    PyObject *ret_lst = NULL, *s_ret = NULL;
    double rcond;
    const char *lapack_driver = NULL;

    int info = 0;
    SliceStatusVec vec_status;

    // Get the input array
    if (!PyArg_ParseTuple(args, "O!O!ds", &PyArray_Type, (PyObject **)&ap_Am,  &PyArray_Type, (PyObject **)&ap_b, &rcond, &lapack_driver)) {
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

    if (rcond < 0) {
        PyErr_SetString(PyExc_ValueError, "Expected rcond >= 0.");
        return NULL;
    }

    // Sanity checks of array dimensions
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    if (ndim < 2) {
        PyErr_SetString(PyExc_ValueError, "Expected at least a 2D array.");
        return NULL;
    }

    // At the python call site,
    // 1) 1D `b` must have been converted into 2D, and
    // 2) batch dimensions of `a` and `b` have been broadcast
    // Therefore, if `a.shape == (s, p, r, m, n)`, then `b.shape == (s, p, r, m, nrhs)`
    // where `nrhs` is the number of right-hand-sides.
    npy_intp ndim_b = PyArray_NDIM(ap_b);
    npy_intp *shape_b = PyArray_SHAPE(ap_b);

    bool dims_match = ndim_b == ndim;
    if (dims_match) {
        for (int i=0; i<ndim-2; i++) {
            dims_match = dims_match && (shape[i] == shape_b[i]);
        }
    }
    if (!dims_match){
        PyErr_SetString(PyExc_ValueError, "`a` and `b` shape mismatch.");
        return NULL;
    }

    npy_intp m = shape[ndim - 2];
    npy_intp n = shape[ndim - 1];
    npy_intp min_mn = m < n ? m : n;
    npy_intp nrhs = PyArray_DIM(ap_b, ndim-1);

    // Allocate the output(s)
    npy_intp shape_1[NPY_MAXDIMS];
    for(int i=0; i<PyArray_NDIM(ap_Am); i++) {
        shape_1[i] = PyArray_DIM(ap_Am, i);
    }

    // x.shape = (..., N, NRHS)
    shape_1[ndim-2] = n;
    shape_1[ndim-1] = nrhs;
    ap_x = (PyArrayObject *)PyArray_SimpleNew(ndim, shape_1, typenum);
    if (!ap_x) {
        PyErr_NoMemory();
        goto fail;
    }

    // S array is not used by ?gelsy
    if (strcmp(lapack_driver, "gelsy") != 0) {
        // S.dtype is real if A.dtype is complex
        npy_intp typenum_S = typenum;
        if (typenum_S == NPY_COMPLEX64) { typenum_S = NPY_FLOAT32; }
        else if (typenum_S == NPY_COMPLEX128) { typenum_S = NPY_FLOAT64; }

        // S.shape = (..., min_mn)
        shape_1[ndim - 2] = min_mn;
        ap_S = (PyArrayObject *)PyArray_SimpleNew(ndim-1, shape_1, typenum_S);
        if (!ap_S) {
            PyErr_NoMemory();
            goto fail;
        }
    }

    // rank.shape = batch_shape
    ap_rank = (PyArrayObject *)PyArray_SimpleNew(ndim-2, shape_1, NPY_INT64);
    if (!ap_rank) {
        PyErr_NoMemory();
        goto fail;
    }

    switch(typenum) {
        case(NPY_FLOAT32):
            info = _lstsq<float>(ap_Am, ap_b, ap_S, ap_x, ap_rank, (float)rcond, lapack_driver, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _lstsq<double>(ap_Am, ap_b, ap_S, ap_x, ap_rank, rcond, lapack_driver, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _lstsq<npy_complex64>(ap_Am, ap_b, ap_S, ap_x, ap_rank, (float)rcond, lapack_driver, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _lstsq<npy_complex128>(ap_Am, ap_b, ap_S, ap_x, ap_rank, rcond, lapack_driver, vec_status);
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unknown array type.");
            goto fail;
    }

    if (info < 0) {
        // Some fatal error: OOM, LWORK query failed or lwork needed is too large
        PyErr_SetString(PyExc_MemoryError, get_err_mesg("lstsq", lapack_driver, info).c_str());
        goto fail;
    }
    ret_lst = convert_vec_status(vec_status);
    // with 'gelsy', we return None for `s`
    s_ret = ap_S != NULL ? PyArray_Return(ap_S) : Py_None;

    return Py_BuildValue("NNNN", PyArray_Return(ap_x), PyArray_Return(ap_rank), s_ret, ret_lst);

fail:
    Py_XDECREF(ap_x);
    Py_XDECREF(ap_rank);
    Py_XDECREF(ap_S);
    return NULL;
}



static PyObject*
_linalg_eig(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_Am = NULL;
    PyArrayObject *ap_Bm = NULL;
    PyArrayObject *ap_w = NULL;
    PyArrayObject *ap_beta = NULL;
    PyArrayObject *ap_vr = NULL;
    PyArrayObject *ap_vl = NULL;
    int compute_vl=0;
    int compute_vr=1;

    int info = 0;
    SliceStatusVec vec_status;

    // return values
    PyObject *ret_lst = NULL;
    PyObject *vl_ret = NULL;
    PyObject *vr_ret = NULL;
    PyObject *beta_ret = NULL;

    // Get the input array
    if (!PyArg_ParseTuple(args, "O!pp|O!",
            &PyArray_Type, (PyObject **)&ap_Am,
            &compute_vl, &compute_vr,
            &PyArray_Type, (PyObject **)&ap_Bm)
    ) {
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
        return NULL;
    }

    npy_intp n = shape[ndim - 1];

    if (PyArray_DIM(ap_Am, ndim-2) != n) {
        PyErr_SetString(PyExc_ValueError, "Expected a square matrix");
        return NULL;
    }

    // Allocate the output(s)
    // NB: we always allocate/return complex eigvalues/eigvecs even if
    // eigenvalues happen to be on the real axis (and then the python wrapper will
    // cast them to reals for backwards compat.

    npy_intp shape_1[NPY_MAXDIMS];
    for(int i=0; i<ndim; i++) {shape_1[i] = shape[i]; }

    // eigenvalues
    int w_typenum = typenum;
    if (typenum == NPY_FLOAT32) { w_typenum = NPY_COMPLEX64; }
    else if (typenum == NPY_FLOAT64) { w_typenum = NPY_COMPLEX128; }

    ap_w = (PyArrayObject *)PyArray_SimpleNew(ndim-1, shape_1, w_typenum);
    if (ap_w == NULL) {
        PyErr_NoMemory();
        return NULL;
    }

    if (ap_Bm != NULL) {
        ap_beta = (PyArrayObject *)PyArray_SimpleNew(ndim-1, shape_1, typenum);
        if (ap_beta == NULL) { PyErr_NoMemory(); goto fail; }
    }

    if (compute_vl) {
        ap_vl = (PyArrayObject *)PyArray_SimpleNew(ndim, shape, w_typenum);
        if (ap_vl == NULL) { PyErr_NoMemory(); goto fail; }
    }

    if (compute_vr) {
        ap_vr = (PyArrayObject *)PyArray_SimpleNew(ndim, shape, w_typenum);
        if (ap_vr == NULL) { PyErr_NoMemory(); goto fail; }
    }

    switch(typenum) {
        case(NPY_FLOAT32):
            info = _eig<float>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _eig<double>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _eig<npy_complex64>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _eig<npy_complex128>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, vec_status);
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unknown array type.");
            goto fail;
    }

    if (info < 0) {
        // Either OOM or internal LAPACK error.
        PyErr_SetString(PyExc_RuntimeError, "Memory error in scipy.linalg.eig.");
        goto fail;
    }

    // normal return
    ret_lst = convert_vec_status(vec_status);

    vl_ret = (ap_vl == NULL) ? Py_None : PyArray_Return(ap_vl);
    vr_ret = (ap_vr == NULL) ? Py_None : PyArray_Return(ap_vr);
    beta_ret = (ap_beta == NULL) ? Py_None : PyArray_Return(ap_beta);

    return Py_BuildValue("NNNNN", PyArray_Return(ap_w), beta_ret, vl_ret, vr_ret, ret_lst);

fail:
    Py_DECREF(ap_w);
    Py_XDECREF(ap_beta);
    Py_XDECREF(ap_vl);
    Py_XDECREF(ap_vr);
    return NULL;
}


/*
 * Helper: convert a vector of slice error statuses to list of dicts
 */
PyObject*
convert_vec_status(SliceStatusVec& vec_status) {
    PyObject *ret_dct = NULL;
    PyObject *ret_lst = PyList_New(0);

    if (vec_status.empty()) {
        return ret_lst;
    }

    // Problems detected in some slices, report.
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
        Py_DECREF(ret_dct);
    }

    return ret_lst;
}


/*
 * Helper for fatal errors (memory errors, mostly).
 */
std::string
get_err_mesg(const std::string routine, const std::string func_name, int info) {
    std::string mesg;
    mesg = "Memory error in scipy.linalg." + routine;
    mesg += " (Internal " + func_name + " returned " + std::to_string(info) + ").";
    return mesg;
}


static char doc_inv[] = ("Compute the matrix inverse.");
static char doc_solve[] = ("Solve the linear system of equations.");
static char doc_solve_banded[] = ("Solve the banded linear system of equations.");
static char doc_svd[] = ("SVD factorization.");
static char doc_lstsq[] = ("linear least squares.");
static char doc_eig[] = ("eigenvalue solver.");

static struct PyMethodDef inv_module_methods[] = {
  {"_inv", _linalg_inv, METH_VARARGS, doc_inv},
  {"_solve", _linalg_solve, METH_VARARGS, doc_solve},
  {"_solve_banded", _linalg_solve_banded, METH_VARARGS, doc_solve_banded},
  {"_svd", _linalg_svd, METH_VARARGS, doc_svd},
  {"_lstsq", _linalg_lstsq, METH_VARARGS, doc_lstsq},
  {"_eig", _linalg_eig, METH_VARARGS, doc_eig},
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
