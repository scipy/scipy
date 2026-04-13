#include <cstring>
#include "_linalg_inv.hh"
#include "_linalg_solve.hh"
#include "_linalg_svd.hh"
#include "_linalg_lstsq.hh"
#include "_linalg_eig.hh"
#include "_linalg_cholesky.hh"
#include "_linalg_qr.hh"
#include "_common_array_utils.hh"
#include "_linalg_lu_det.hh"


static PyObject* _linalg_inv_error;


PyObject*
convert_vec_status(SliceStatusVec& vec_status);

std::string
get_err_mesg(const std::string routine, const std::string func_name, int info);


/*
 * Minimal checks for an input array:
 *   - ndim >= 2
 *   - LAPACK-compatible dtype
 *   - aligned
 *
 * Sets an error if any of these is violated (users to bail out immediately).
 */
int _check_dtype_and_flags(PyArrayObject *ap_Am, const char *func) {
    int typenum = PyArray_TYPE(ap_Am);
    bool dtype_ok = (typenum == NPY_FLOAT32)
                     || (typenum == NPY_FLOAT64)
                     || (typenum == NPY_COMPLEX64)
                     || (typenum == NPY_COMPLEX128);

    if (!dtype_ok || !PyArray_ISALIGNED(ap_Am)) {
        PyErr_Format(PyExc_TypeError, "scipy.linag.%s : Expected a real or complex array.", func);
        return 0;
    }

    int ndim = PyArray_NDIM(ap_Am);
    if (ndim < 2) {
        PyErr_Format(PyExc_ValueError, "scipy.linalg.%s: Expected at least a 2D array.", func);
        return 0;
    }

    return 1;
}


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

    // Sanity check the input array
    if (!_check_dtype_and_flags(ap_Am, "inv")) {
        return NULL;
    }

    int typenum = PyArray_TYPE(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_Am);      // Array shape
    npy_intp n = shape[ndim - 1];                // Slice size
    if (n != shape[ndim - 2]) {
        PyErr_SetString(PyExc_ValueError, "Last two dimensions of the input must be the same.");
        return NULL;
    }

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
    if (!PyArg_ParseTuple(args, "O!O!|npppp",
        &PyArray_Type, (PyObject **)&ap_Am,
        &PyArray_Type, (PyObject **)&ap_b,
        &structure,
        &lower, &transposed,
        &overwrite_a, &overwrite_b)
    ){
        return NULL;
    }

    // Sanity check the input array
    if (!_check_dtype_and_flags(ap_Am, "solve")) {
        return NULL;
    }

    // Sanity check shapes
    int typenum = PyArray_TYPE(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp* shape = PyArray_SHAPE(ap_Am);
    if (shape[ndim - 1] != shape[ndim - 2]) {
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

    if (!overwrite_b) {
        /* Allocate the output */
        ap_x = (PyArrayObject *)PyArray_SimpleNew(ndim_b, shape_b, typenum);
        if(!ap_x) {
            PyErr_NoMemory();
            return NULL;
        }
    }
    else {
        /* Reuse the memory buffer of the input array. */
        ap_x = ap_b;
        Py_INCREF(ap_b);
    }

    void *buf = PyArray_DATA(ap_x);
    switch(typenum) {
        case(NPY_FLOAT32):
            info = _solve<float>(ap_Am, ap_b, (float *)buf, structure, lower, transposed, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _solve<double>(ap_Am, ap_b, (double *)buf, structure, lower, transposed, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _solve<npy_complex64>(ap_Am, ap_b, (npy_complex64 *)buf, structure, lower, transposed, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _solve<npy_complex128>(ap_Am, ap_b, (npy_complex128 *)buf, structure, lower, transposed, overwrite_a, overwrite_b, vec_status);
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unknown array type.");
            return NULL;
    }

    if (info < 0) {
        // Either OOM error or requested lwork too large.
        Py_DECREF(ap_x);
        PyErr_SetString(PyExc_MemoryError, "Memory error in scipy.linalg.solve.");
        return NULL;
    }
    PyObject *ret_lst = convert_vec_status(vec_status);

    return Py_BuildValue("NN", PyArray_Return(ap_x), ret_lst);
}



static PyObject*
_linalg_qr(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_A = NULL;

    int info = 0;
    SliceStatusVec vec_status;
    int overwrite_a = 0;
    QR_mode mode = QR_mode::FULL;
    int pivoting = 0;

    PyArrayObject *ap_Q = NULL, *ap_R = NULL, *ap_tau = NULL, *ap_jpvt = NULL;
    PyObject *ret_lst = NULL, *ret_Q = NULL, *ret_tau = NULL, *ret_jpvt = NULL;

    // Get the input array
    if (!PyArg_ParseTuple(args, "O!|pnp", &PyArray_Type, (PyObject **)&ap_A, &overwrite_a, &mode, &pivoting)) {
        return NULL;
    }

    // Sanity check the input array
    if (!_check_dtype_and_flags(ap_A, "qr")) {
        return NULL;
    }

    int typenum = PyArray_TYPE(ap_A);
    int ndim = PyArray_NDIM(ap_A);
    npy_intp* shape = PyArray_SHAPE(ap_A);

    // -------------------------------------------------------------------
    // Conditionally allocate return objects
    // -------------------------------------------------------------------
    npy_intp M = shape[ndim-2], N = shape[ndim-1];
    npy_intp K = std::min(M, N);

    npy_intp shape_Q[NPY_MAXDIMS];
    npy_intp shape_R[NPY_MAXDIMS];

    for (npy_intp i = 0; i < ndim; i++) {
        shape_Q[i] = shape[i];
        shape_R[i] = shape[i];
    }

    switch (mode) {
        case QR_mode::FULL:
        {
            shape_Q[ndim-1] = M;
            shape_R[ndim-1] = N;
            break;
        }

        case QR_mode::R:
        {
            // shape of `Q` irrelevant here
            shape_R[ndim-1] = N;
            break; // shape of `Q` irrelevant here
        }

        case QR_mode::RAW_MODE:
        {
            shape_Q[ndim-1] = N;
            shape_R[ndim-2] = K;
            shape_R[ndim-1] = N;
            break;
        }

        case QR_mode::ECONOMIC:
        {
            shape_Q[ndim-1] = K;
            shape_R[ndim-2] = K;
            shape_R[ndim-1] = N;
            break;
        }
    }

    /*
     * Allocation strategy:
     *
     * - Q: always needed except if `mode == r`. In the case that `mode == economic` or `mode == raw`
     *      the input buffer can be re-used if `overwrite_a` is set. (Albeit that it is possible that
     *      if M < N `Q` will have to be shrunk down at the python side to become `M x M`.)
     *
     *      An additional comment is that `Q` is returned with F-ordered slices for `mode == raw` as its
     *      output is destined to be fed into LAPACK, allowing fast-paths in memory copies.
     *
     * - R: always needed. If `mode == raw` or `mode == economic` a new array is always allocated for this.
     *      For the other two modes the input buffer can be re-used if `overwrite_a` is set.
     *
     * - tau: only allocated if `mode == raw`. Else a temporary buffer is allocated in the main loop as the
     *        reflectors shouldn't be stored anyways.
     *
     * - jpvt: only allocated if `pivoting` is set.
     */
    if (mode != QR_mode::R) {
        if (!overwrite_a || mode == QR_mode::FULL) {
            ap_Q = (PyArrayObject *)PyArray_SimpleNew(ndim, shape_Q, typenum);
            if (!ap_Q) {
                PyErr_NoMemory();
                goto fail_qr;
            }
        } else {
            Py_INCREF(ap_A);
            ap_Q = ap_A;
        }
    }

    if (mode == QR_mode::RAW_MODE || mode == QR_mode::ECONOMIC || (!overwrite_a && (mode == QR_mode::R || mode == QR_mode::FULL))) {
        // Return C-ordered object, hence `0` flag. Already pre-fill with zeros for efficiency.
        ap_R = (PyArrayObject *)PyArray_ZEROS(ndim, shape_R, typenum, 0);
        if (!ap_R) {
            PyErr_NoMemory();
            goto fail_qr;
        }
    } else {
        Py_INCREF(ap_A);
        ap_R = ap_A;
    }

    if (mode == QR_mode::RAW_MODE) {
        shape_Q[ndim-2] = K; // Just reuse `shape_Q`; not used any longer.
        ap_tau = (PyArrayObject *)PyArray_SimpleNew(ndim-1, shape_Q, typenum); // Just a vector, so `ndim-1`.
        if (!ap_tau) {
            PyErr_NoMemory();
            goto fail_qr;
        }

        /*
         * Since the point of the `Q` returned if `mode="raw"` is mostly to call other LAPACK routines (e.g. `qr_multiply`),
         * set the strides of `Q` to return F-ordered slices.
         *
         * Dimensions other than the last two are not affected so slice computation etc. stays intact.
         *
         * Since `overwrite_a` is only enabled when the input has F-ordered slices, this is not necessary in that case.
         */
        if (!overwrite_a) {
            npy_intp *strides_Q = PyArray_STRIDES(ap_Q);
            int sizeof_T = PyArray_ITEMSIZE(ap_Q);

            strides_Q[ndim-2] = sizeof_T;
            strides_Q[ndim-1] = M * sizeof_T;

            PyArray_UpdateFlags(ap_Q, NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_F_CONTIGUOUS);
        }
    }

    // Set all elements to 0 immediately as the pivots are also used as inputs in `geqp3`.
    // Allocate a C-ordered array, (hence the `0` magic number).
    if (pivoting) {
        shape_Q[ndim-2] = N;
        ap_jpvt = (PyArrayObject *)PyArray_ZEROS(ndim-1, shape_Q, sizeof(CBLAS_INT) == sizeof(NPY_INT32)? NPY_INT32 : NPY_INT64, 0);
        if (!ap_jpvt) {
            PyErr_NoMemory();
            goto fail_qr;
        }
    }

    switch(typenum) {
        case(NPY_FLOAT32):
            info = _qr<float>(ap_A, ap_Q, ap_R, ap_tau, ap_jpvt, overwrite_a, mode, pivoting, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _qr<double>(ap_A, ap_Q, ap_R, ap_tau, ap_jpvt, overwrite_a, mode, pivoting, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _qr<npy_complex64>(ap_A, ap_Q, ap_R, ap_tau, ap_jpvt, overwrite_a, mode, pivoting, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _qr<npy_complex128>(ap_A, ap_Q, ap_R, ap_tau, ap_jpvt, overwrite_a, mode, pivoting, vec_status);
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unknown array type.");
            return NULL;
    }

    if (info < 0) {
        // Either OOM error or requested lwork too large.
        PyErr_SetString(PyExc_MemoryError, "Memory error in scipy.linalg.qr.");
        goto fail_qr;
    }

    ret_lst = convert_vec_status(vec_status);

    ret_Q = (mode != QR_mode::R) ? PyArray_Return(ap_Q) : Py_None;
    ret_tau = (mode == QR_mode::RAW_MODE) ? PyArray_Return(ap_tau) : Py_None;
    ret_jpvt = (pivoting) ? PyArray_Return(ap_jpvt): Py_None;
    return Py_BuildValue("NNNNN", ret_Q, PyArray_Return(ap_R), ret_tau, ret_jpvt, ret_lst);

fail_qr:
    Py_XDECREF(ap_Q);
    Py_XDECREF(ap_R);
    Py_XDECREF(ap_tau);
    Py_XDECREF(ap_jpvt);
    return NULL;
}



static PyObject*
_linalg_svd(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_Am = NULL;
    const char *lapack_driver = NULL;
    int compute_uv = 1;
    int full_matrices = 1;
    int overwrite_a = 0;
    PyArrayObject *ap_S = NULL, *ap_U = NULL, *ap_Vh = NULL;
    PyObject *ret_lst;

    int info = 0;
    SliceStatusVec vec_status;

    // Get the input array
    if (!PyArg_ParseTuple(args, "O!s|ppp", &PyArray_Type, (PyObject **)&ap_Am,  &lapack_driver, &compute_uv, &full_matrices, &overwrite_a)) {
        return NULL;
    }

    // Sanity check the input array
    if (!_check_dtype_and_flags(ap_Am, "svd")) {
        return NULL;
    }

    int typenum = PyArray_TYPE(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);

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
            info = _svd<float>(ap_Am, ap_U, ap_S, ap_Vh, jobz, lapack_driver, overwrite_a, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _svd<double>(ap_Am, ap_U, ap_S, ap_Vh, jobz, lapack_driver, overwrite_a, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _svd<npy_complex64>(ap_Am, ap_U, ap_S, ap_Vh, jobz, lapack_driver, overwrite_a, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _svd<npy_complex128>(ap_Am, ap_U, ap_S, ap_Vh, jobz, lapack_driver, overwrite_a, vec_status);
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
    int overwrite_a = 0;
    int overwrite_b = 0;

    int info = 0;
    SliceStatusVec vec_status;

    // Get the input array
    if (!PyArg_ParseTuple(args, "O!O!ds|pp", &PyArray_Type, (PyObject **)&ap_Am,  &PyArray_Type, (PyObject **)&ap_b, &rcond, &lapack_driver, &overwrite_a, &overwrite_b)) {
        return NULL;
    }

    // Sanity check the input array
    if (!_check_dtype_and_flags(ap_Am, "lstsq")) {
        return NULL;
    }

    if (rcond < 0) {
        PyErr_SetString(PyExc_ValueError, "Expected rcond >= 0.");
        return NULL;
    }

    int typenum = PyArray_TYPE(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);

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
    CBLAS_INT max_mn = m > n ? m : n;
    CBLAS_INT ldb = max_mn > 1 ? max_mn : 1;
    npy_intp nrhs = PyArray_DIM(ap_b, ndim-1);

    // Allocate the output(s)
    npy_intp shape_1[NPY_MAXDIMS];
    for(int i=0; i<PyArray_NDIM(ap_Am); i++) {
        shape_1[i] = PyArray_DIM(ap_Am, i);
    }

    // x.shape = (..., N, NRHS)
    shape_1[ndim-2] = ldb;
    shape_1[ndim-1] = nrhs;
    if (!overwrite_b) { // Will only work if m > n, python side should have caught this.
        ap_x = (PyArrayObject *)PyArray_SimpleNew(ndim, shape_1, typenum);
        if (!ap_x) {
            PyErr_NoMemory();
            goto fail;
        }
    } else {
        Py_INCREF(ap_b);
        ap_x = ap_b;
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
            info = _lstsq<float>(ap_Am, ap_b, ap_S, ap_x, ap_rank, (float)rcond, lapack_driver, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _lstsq<double>(ap_Am, ap_b, ap_S, ap_x, ap_rank, rcond, lapack_driver, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _lstsq<npy_complex64>(ap_Am, ap_b, ap_S, ap_x, ap_rank, (float)rcond, lapack_driver, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _lstsq<npy_complex128>(ap_Am, ap_b, ap_S, ap_x, ap_rank, rcond, lapack_driver, overwrite_a, overwrite_b, vec_status);
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
    int overwrite_a=0;
    int overwrite_b=0;

    int info = 0;
    SliceStatusVec vec_status;

    // return values
    PyObject *ret_lst = NULL;
    PyObject *vl_ret = NULL;
    PyObject *vr_ret = NULL;
    PyObject *beta_ret = NULL;

    // Get the input array
    if (!PyArg_ParseTuple(args, "O!pppp|O!",
            &PyArray_Type, (PyObject **)&ap_Am,
            &compute_vl, &compute_vr,
            &overwrite_a, &overwrite_b,
            &PyArray_Type, (PyObject **)&ap_Bm)
    ) {
        return NULL;
    }

    // Sanity check the input array
    if (!_check_dtype_and_flags(ap_Am, "eig")) {
        return NULL;
    }

    int typenum = PyArray_TYPE(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);

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

        // Sanity check the input array
        if (!_check_dtype_and_flags(ap_Bm, "eig")) {
            return NULL;
        }

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
            info = _eig<float>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _eig<double>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _eig<npy_complex64>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _eig<npy_complex128>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, overwrite_a, overwrite_b, vec_status);
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



static PyObject*
_linalg_cholesky(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_Am = NULL;
    PyArrayObject *ap_Cm = NULL;
    PyObject *ret_lst = NULL;
    int lower=0, overwrite_a=0, clean=1;

    int info = 0;
    SliceStatusVec vec_status;

    // Get input
    if (!PyArg_ParseTuple(args, "O!|ppp", &PyArray_Type, (PyObject **)&ap_Am, &lower, &overwrite_a, &clean)) {
        return NULL;
    }

    // Sanity check the input array
    if (!_check_dtype_and_flags(ap_Am, "cholesky")) {
        return NULL;
    }

    int typenum = PyArray_TYPE(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    npy_intp n = shape[ndim - 1];
    if (PyArray_DIM(ap_Am, ndim-2) != n) {
        PyErr_SetString(PyExc_ValueError, "Expected a square matrix");
        return NULL;
    }

    // Allocate the output array if needed
    if (!overwrite_a) {
        ap_Cm = (PyArrayObject *)PyArray_ZEROS(ndim, shape, typenum, 0); // 0 to obtain C-ordered input
        if (ap_Cm == NULL) {
            PyErr_NoMemory();
            goto fail;
        }
    } else {
        Py_INCREF(ap_Am);
        ap_Cm = ap_Am;
    }

    switch(typenum) {
        case(NPY_FLOAT32):
            info = _cholesky<float>(ap_Am, ap_Cm, lower, overwrite_a, clean, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _cholesky<double>(ap_Am, ap_Cm, lower, overwrite_a, clean, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _cholesky<npy_complex64>(ap_Am, ap_Cm, lower, overwrite_a, clean, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _cholesky<npy_complex128>(ap_Am, ap_Cm, lower, overwrite_a, clean, vec_status);
            break;
    }

    if (info < 0) {
        // Out-of-memory or scipy internal error
        PyErr_SetString(PyExc_RuntimeError, "Memory error in scipy.linalg.cholesky.");
        goto fail;
    }

    // normal return
    ret_lst = convert_vec_status(vec_status);
    return Py_BuildValue("NN", PyArray_Return(ap_Cm), ret_lst);

fail:
    Py_XDECREF(ap_Cm);
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


static PyObject*
_linalg_lu(PyObject* Py_UNUSED(dummy), PyObject* args) {

    // Handle the input and metadata
    PyArrayObject *ap_a = NULL;
    int permute_l = 0;
    int overwrite_a = 0;

    if (!PyArg_ParseTuple(args, "O!pp", &PyArray_Type, (PyObject **)&ap_a, &permute_l, &overwrite_a)) {
        return NULL;
    }

    // Sanity check the input array
    if (!_check_dtype_and_flags(ap_a, "lu")) {
        return NULL;
    }

    int typenum = PyArray_TYPE(ap_a);
    int ndim = PyArray_NDIM(ap_a);
    if (ndim > LU_MAX_NDIM) {
        PyErr_SetString(PyExc_ValueError, "scipy.linalg.lu: Input array has too many dimensions.");
        return NULL;
    }

    PyArrayObject *ap_l = NULL;
    PyArrayObject *ap_u = NULL;
    PyArrayObject *ap_perm = NULL;
    void *scratch = NULL;
    CBLAS_INT *ipiv = NULL;
    int info = 0;
    CBLAS_INT *slice_info = NULL;
    LU_Context ctx = {};
    PyObject *result = NULL;

    npy_intp *shape = PyArray_SHAPE(ap_a);
    npy_intp *byte_strides = PyArray_STRIDES(ap_a);
    npy_intp elem_size = PyArray_ITEMSIZE(ap_a);
    CBLAS_INT m = (CBLAS_INT)shape[ndim - 2];
    CBLAS_INT n = (CBLAS_INT)shape[ndim - 1];
    CBLAS_INT minmn = m < n ? m : n;

    npy_intp num_of_slices = 1;
    for (int i = 0; i < ndim - 2; i++) { num_of_slices *= shape[i]; }

    int64_t ctx_shape[LU_MAX_NDIM];
    int64_t ctx_strides[LU_MAX_NDIM];
    for (int i = 0; i < ndim; i++) {
        ctx_shape[i] = (int64_t)shape[i];
        ctx_strides[i] = (int64_t)(byte_strides[i] / elem_size);
    }

    // Allocate perm (..., m)
    ap_perm = (PyArrayObject *)PyArray_SimpleNew(ndim - 1, shape, (sizeof(CBLAS_INT) == 8) ? NPY_INT64 : NPY_INT32);
    if (!ap_perm) { PyErr_NoMemory(); goto fail; }

    // Allocate scratch buffers (reused across slices)
    scratch = PyMem_Malloc(m * n * elem_size);
    ipiv = (CBLAS_INT*)PyMem_Malloc(minmn * sizeof(CBLAS_INT));
    slice_info = (CBLAS_INT*)PyMem_Calloc(num_of_slices, sizeof(CBLAS_INT));
    if (!scratch || !ipiv || !slice_info) { PyErr_NoMemory(); goto fail; }

    ctx.shape = ctx_shape;
    ctx.strides = ctx_strides;
    ctx.num_of_slices = num_of_slices;
    ctx.m = m;
    ctx.n = n;
    ctx.perm = (CBLAS_INT*)PyArray_DATA(ap_perm);
    ctx.ipiv = ipiv;
    ctx.ndim = ndim;
    ctx.permute_l = (bool)permute_l;
    ctx.overwrite_a = false;

    // Legacy overwrite_a: reuse input as the larger factor (2D + F-contiguous only).
    // Disabled when permute_l is true (permute_rows needs C-ordered data).
    if (overwrite_a && !permute_l && (ndim == 2) && PyArray_ISFARRAY_RO(ap_a) && PyArray_ISWRITEABLE(ap_a)) {

        ctx.overwrite_a = true;
        npy_intp other_shape[LU_MAX_NDIM];
        for (int i = 0; i < ndim; i++) { other_shape[i] = shape[i]; }

        if (m > n) {
            // Tall: L (m x n) is large — reuse input. Allocate U (n x n).
            ap_l = ap_a;
            Py_INCREF(ap_l);
            other_shape[ndim - 2] = minmn;
            ap_u = (PyArrayObject *)PyArray_ZEROS(ndim, other_shape, typenum, 0);
            if (!ap_u) { PyErr_NoMemory(); goto fail; }
        } else {
            // Fat/square: U (m x n) is large — reuse input. Allocate L (m x m).
            ap_u = ap_a;
            Py_INCREF(ap_u);
            other_shape[ndim - 1] = minmn;
            ap_l = (PyArrayObject *)PyArray_ZEROS(ndim, other_shape, typenum, 0);
            if (!ap_l) { PyErr_NoMemory(); goto fail; }
        }

    } else {

        // Allocate output arrays: L (..., m, minmn), U (..., minmn, n)
        shape[ndim - 1] = minmn;
        ap_l = (PyArrayObject *)PyArray_ZEROS(ndim, shape, typenum, 0);
        if (!ap_l) { PyErr_NoMemory(); shape[ndim - 1] = n; goto fail; }
        shape[ndim - 1] = n;

        shape[ndim - 2] = minmn;
        ap_u = (PyArrayObject *)PyArray_ZEROS(ndim, shape, typenum, 0);
        if (!ap_u) { PyErr_NoMemory(); shape[ndim - 2] = m; goto fail; }
        shape[ndim - 2] = m;
    }

    // Dispatch to templated C++ code
    switch (typenum) {
        case NPY_FLOAT32:
            info = lu_dispatch<float>(ctx, (float*)PyArray_DATA(ap_a), (float*)PyArray_DATA(ap_l), (float*)PyArray_DATA(ap_u), (float*)scratch, slice_info);
            break;
        case NPY_FLOAT64:
            info = lu_dispatch<double>(ctx, (double*)PyArray_DATA(ap_a), (double*)PyArray_DATA(ap_l), (double*)PyArray_DATA(ap_u), (double*)scratch, slice_info);
            break;
        case NPY_COMPLEX64:
            info = lu_dispatch<std::complex<float>>(ctx, (std::complex<float>*)PyArray_DATA(ap_a), (std::complex<float>*)PyArray_DATA(ap_l), (std::complex<float>*)PyArray_DATA(ap_u), (std::complex<float>*)scratch, slice_info);
            break;
        case NPY_COMPLEX128:
            info = lu_dispatch<std::complex<double>>(ctx, (std::complex<double>*)PyArray_DATA(ap_a), (std::complex<double>*)PyArray_DATA(ap_l), (std::complex<double>*)PyArray_DATA(ap_u), (std::complex<double>*)scratch, slice_info);
            break;
    }

    if (info < 0) {
        PyErr_SetString(PyExc_ValueError, get_err_mesg("lu", "getrf", info).c_str());
        goto fail;
    }

    // Free scratch buffers (no longer needed)
    PyMem_Free(scratch); scratch = NULL;
    PyMem_Free(ipiv); ipiv = NULL;

    // NOTE: slice_info contains per-slice LAPACK info (>0 means singular).
    // Not yet surfaced to Python — singular LU factors are still valid.
    // The intention is to eventually provide per-slice diagnostics.
    PyMem_Free(slice_info); slice_info = NULL;

    // Always return (P, L, U) with consistent types.
    // When permute_l is true, P is replaced with an empty array since
    // the permutation is already applied to L.
    if (permute_l) {
        Py_DECREF(ap_perm);
        // Return an empty int32 array instead of None so the return type is
        // always (ndarray, ndarray, ndarray) — avoids Optional in type stubs.
        // int32 since the array has zero elements; no overflow concern.
        npy_intp zero = 0;
        ap_perm = (PyArrayObject *)PyArray_EMPTY(1, &zero, NPY_INT32, 0);
        if (!ap_perm) { PyErr_NoMemory(); goto fail; }
    }
    result = Py_BuildValue("NNN", PyArray_Return(ap_perm), PyArray_Return(ap_l), PyArray_Return(ap_u));
    return result;

fail:
    PyMem_Free(scratch);
    PyMem_Free(ipiv);
    PyMem_Free(slice_info);
    Py_XDECREF(ap_l);
    Py_XDECREF(ap_u);
    Py_XDECREF(ap_perm);
    return NULL;
}


// ========================================================================
// Determinant
// ========================================================================

static PyObject*
_linalg_det(PyObject* Py_UNUSED(dummy), PyObject* args) {

    // Handle the input and metadata
    PyArrayObject *ap_a = NULL;
    int overwrite_a = 0;

    if (!PyArg_ParseTuple(args, "O!p", &PyArray_Type, (PyObject **)&ap_a, &overwrite_a)) {
        return NULL;
    }

    // Sanity check the input array
    if (!_check_dtype_and_flags(ap_a, "det")) {
        return NULL;
    }

    int typenum = PyArray_TYPE(ap_a);
    int ndim = PyArray_NDIM(ap_a);
    if (ndim > LU_MAX_NDIM) {
        PyErr_SetString(PyExc_ValueError, "scipy.linalg.det: Input array has too many dimensions.");
        return NULL;
    }

    npy_intp *shape = PyArray_SHAPE(ap_a);
    npy_intp *byte_strides = PyArray_STRIDES(ap_a);
    npy_intp elem_size = PyArray_ITEMSIZE(ap_a);
    CBLAS_INT n = (CBLAS_INT)shape[ndim - 1];

    if (shape[ndim - 2] != shape[ndim - 1]) {
        PyErr_SetString(PyExc_ValueError, "scipy.linalg.det: Last two dimensions must be square.");
        return NULL;
    }

    // All variables declared before any goto (C++ initialization rule)
    PyArrayObject *ap_det = NULL;
    void *scratch = NULL;
    CBLAS_INT *ipiv = NULL;
    CBLAS_INT *slice_info = NULL;
    LU_Context ctx = {};
    int info = 0;

    npy_intp num_of_slices = 1;
    for (int i = 0; i < ndim - 2; i++) { num_of_slices *= shape[i]; }

    int64_t ctx_shape[LU_MAX_NDIM];
    int64_t ctx_strides[LU_MAX_NDIM];
    for (int i = 0; i < ndim; i++) {
        ctx_shape[i] = (int64_t)shape[i];
        ctx_strides[i] = (int64_t)(byte_strides[i] / elem_size);
    }

    // Allocate output det array (...,) — batch dimensions only
    ap_det = (PyArrayObject *)PyArray_SimpleNew(ndim - 2, shape, typenum);
    if (!ap_det) { PyErr_NoMemory(); goto fail; }

    // Allocate scratch buffers
    scratch = PyMem_Malloc(n * n * elem_size);
    ipiv = (CBLAS_INT*)PyMem_Malloc(n * sizeof(CBLAS_INT));
    slice_info = (CBLAS_INT*)PyMem_Calloc(num_of_slices, sizeof(CBLAS_INT));
    if (!scratch || !ipiv || !slice_info) { PyErr_NoMemory(); goto fail; }

    ctx.shape = ctx_shape;
    ctx.strides = ctx_strides;
    ctx.num_of_slices = num_of_slices;
    ctx.m = n;
    ctx.n = n;
    ctx.perm = NULL;
    ctx.ipiv = ipiv;
    ctx.ndim = ndim;
    ctx.permute_l = false;
    ctx.overwrite_a = false;

    // Legacy overwrite_a: use input directly (2D + contiguous only).
    // For det, both C and F order work since det(A^T) = det(A).
    if (overwrite_a && (ndim == 2) && PyArray_ISWRITEABLE(ap_a)
        && (PyArray_ISFARRAY_RO(ap_a) || PyArray_ISCARRAY_RO(ap_a))) {
        ctx.overwrite_a = true;
    }

    // Dispatch to templated C++ code
    switch (typenum) {
        case NPY_FLOAT32:
            info = det_dispatch<float>(ctx, (float*)PyArray_DATA(ap_a), (float*)PyArray_DATA(ap_det), (float*)scratch, slice_info);
            break;
        case NPY_FLOAT64:
            info = det_dispatch<double>(ctx, (double*)PyArray_DATA(ap_a), (double*)PyArray_DATA(ap_det), (double*)scratch, slice_info);
            break;
        case NPY_COMPLEX64:
            info = det_dispatch<std::complex<float>>(ctx, (std::complex<float>*)PyArray_DATA(ap_a), (std::complex<float>*)PyArray_DATA(ap_det), (std::complex<float>*)scratch, slice_info);
            break;
        case NPY_COMPLEX128:
            info = det_dispatch<std::complex<double>>(ctx, (std::complex<double>*)PyArray_DATA(ap_a), (std::complex<double>*)PyArray_DATA(ap_det), (std::complex<double>*)scratch, slice_info);
            break;
    }

    if (info < 0) {
        PyErr_SetString(PyExc_ValueError, get_err_mesg("det", "getrf", info).c_str());
        goto fail;
    }

    PyMem_Free(scratch); scratch = NULL;
    PyMem_Free(ipiv); ipiv = NULL;

    // NOTE: slice_info[idx] > 0 means getrf detected exact singularity at
    // diagonal position info for that slice — the determinant is exactly zero,
    // not merely numerically small. This can be used to distinguish exact
    // singularity from near-singularity in the future.
    PyMem_Free(slice_info); slice_info = NULL;

    return (PyObject *)ap_det;

fail:
    PyMem_Free(scratch);
    PyMem_Free(ipiv);
    PyMem_Free(slice_info);
    Py_XDECREF(ap_det);
    return NULL;
}


template<typename T>
inline void
bandwidth_strided_scalar(const void *data, int64_t offset,
                         int64_t n, int64_t m, int64_t s0, int64_t s1,
                         int64_t *lower_band, int64_t *upper_band)
{
    const T *slice = (const T *)data + offset;
    const T zero = T(0);
    int64_t lb = 0, ub = 0;

    for (int64_t r = n - 1; r > 0; r--) {
        int64_t limit = r - lb;
        if (limit > m) { limit = m; }
        for (int64_t c = 0; c < limit; c++) {
            if (slice[r * s0 + c * s1] != zero) { lb = r - c; break; }
        }
        if (r <= lb) { break; }
    }

    for (int64_t r = 0; r < n - 1; r++) {
        for (int64_t c = m - 1; c > r + ub; c--) {
            if (slice[r * s0 + c * s1] != zero) { ub = c - r; break; }
        }
        if (r + ub + 1 > m) { break; }
    }

    *lower_band = lb;
    *upper_band = ub;
}


template<typename T>
inline void
bandwidth_contiguous_scalar(const void *data, int64_t offset,
                            int64_t n, int64_t m,
                            int64_t *lower_band, int64_t *upper_band)
{
    const T *slice = (const T *)data + offset;
    const T zero = T(0);
    int64_t lb = 0, ub = 0;

    for (int64_t r = n - 1; r > 0; r--) {
        int64_t limit = r - lb;
        if (limit > m) { limit = m; }
        for (int64_t c = 0; c < limit; c++) {
            if (slice[r * m + c] != zero) { lb = r - c; break; }
        }
        if (r <= lb) { break; }
    }

    for (int64_t r = 0; r < n - 1; r++) {
        for (int64_t c = m - 1; c > r + ub; c--) {
            if (slice[r * m + c] != zero) { ub = c - r; break; }
        }
        if (r + ub + 1 > m) { break; }
    }

    *lower_band = lb;
    *upper_band = ub;
}


static PyObject*
_linalg_bandwidth(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_a = NULL;

    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, (PyObject **)&ap_a)) {
        return NULL;
    }

    int ndim = PyArray_NDIM(ap_a);
    if (ndim < 2) {
        PyErr_SetString(PyExc_ValueError, "Input array must have at least 2 dimensions.");
        return NULL;
    }

    npy_intp *shape = PyArray_SHAPE(ap_a);
    int64_t n = shape[ndim - 2], m = shape[ndim - 1];

    int typenum = PyArray_TYPE(ap_a);
    npy_intp *byte_strides = PyArray_STRIDES(ap_a);
    npy_intp itemsize = PyArray_ITEMSIZE(ap_a);

    // longdouble/clongdouble are rejected by the Python wrapper in _misc.py
    bool has_contiguous = (typenum == NPY_FLOAT32) || (typenum == NPY_FLOAT64)
                          || (typenum == NPY_COMPLEX64) || (typenum == NPY_COMPLEX128);

    // Element strides for all dims
    int64_t elem_strides[LU_MAX_NDIM];
    for (int i = 0; i < ndim; i++) {
        elem_strides[i] = (int64_t)(byte_strides[i] / itemsize);
    }

    // Check if the last two dimensions are C- or F-contiguous
    // C-contiguous: inner slices are always C-contiguous (any ndim)
    // F-contiguous: inner slices are only F-contiguous for 2D
    bool inner_c_contig = PyArray_IS_C_CONTIGUOUS(ap_a);
    bool inner_f_contig = (ndim == 2) && PyArray_IS_F_CONTIGUOUS(ap_a);
    bool use_contiguous = has_contiguous && (inner_c_contig || inner_f_contig);

    int64_t n_eff = inner_f_contig ? m : n;
    int64_t m_eff = inner_f_contig ? n : m;

    int batch_ndim = ndim - 2;
    npy_intp num_slices = 1;
    for (int i = 0; i < batch_ndim; i++) { num_slices *= shape[i]; }

    PyArrayObject *ap_lb = (PyArrayObject *)PyArray_SimpleNew(batch_ndim, shape, NPY_INT64);
    PyArrayObject *ap_ub = (PyArrayObject *)PyArray_SimpleNew(batch_ndim, shape, NPY_INT64);
    if (!ap_lb || !ap_ub) {
        Py_XDECREF(ap_lb); Py_XDECREF(ap_ub);
        return PyErr_NoMemory();
    }

    int64_t *lb_data = (int64_t *)PyArray_DATA(ap_lb);
    int64_t *ub_data = (int64_t *)PyArray_DATA(ap_ub);
    const void *a_data = PyArray_DATA(ap_a);

    for (npy_intp idx = 0; idx < num_slices; idx++) {
        // Compute element offset to the idx-th 2D slice
        int64_t offset = 0;
        int64_t temp = idx;
        for (int i = batch_ndim - 1; i >= 0; i--) {
            offset += (temp % shape[i]) * elem_strides[i];
            temp /= shape[i];
        }

        if (use_contiguous) {
            switch (typenum) {
                case NPY_FLOAT32:    bandwidth_contiguous_scalar<float>                (a_data, offset, n_eff, m_eff, &lb_data[idx], &ub_data[idx]); break;
                case NPY_FLOAT64:    bandwidth_contiguous_scalar<double>               (a_data, offset, n_eff, m_eff, &lb_data[idx], &ub_data[idx]); break;
                case NPY_COMPLEX64:  bandwidth_contiguous_scalar<std::complex<float>>  (a_data, offset, n_eff, m_eff, &lb_data[idx], &ub_data[idx]); break;
                case NPY_COMPLEX128: bandwidth_contiguous_scalar<std::complex<double>> (a_data, offset, n_eff, m_eff, &lb_data[idx], &ub_data[idx]); break;
            }
            if (inner_f_contig) { std::swap(lb_data[idx], ub_data[idx]); }
        } else {
            int64_t s0 = elem_strides[ndim - 2];
            int64_t s1 = elem_strides[ndim - 1];

            // Normalize platform-dependent integer typenums to canonical
            // case labels. On Windows sizeof(int)==sizeof(long)==4 but
            // NPY_INT != NPY_LONG, and NPY_INT32 is only one of them.
            // Map all integer typenums to canonical NPY_INTxx by itemsize.
            int tn = typenum;
            if (PyTypeNum_ISINTEGER(tn) && !PyTypeNum_ISBOOL(tn)) {
                bool is_unsigned = PyTypeNum_ISUNSIGNED(tn);
                switch (itemsize) {
                    case 1: tn = is_unsigned ? NPY_UINT8  : NPY_INT8;  break;
                    case 2: tn = is_unsigned ? NPY_UINT16 : NPY_INT16; break;
                    case 4: tn = is_unsigned ? NPY_UINT32 : NPY_INT32; break;
                    case 8: tn = is_unsigned ? NPY_UINT64 : NPY_INT64; break;
                }
            }

            switch (tn) {
                case NPY_BOOL:        bandwidth_strided_scalar<npy_bool>                  (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                case NPY_INT8:        bandwidth_strided_scalar<npy_int8>                  (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                case NPY_INT16:       bandwidth_strided_scalar<npy_int16>                 (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                case NPY_INT32:       bandwidth_strided_scalar<npy_int32>                 (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                case NPY_INT64:       bandwidth_strided_scalar<npy_int64>                 (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                case NPY_UINT8:       bandwidth_strided_scalar<npy_uint8>                 (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                case NPY_UINT16:      bandwidth_strided_scalar<npy_uint16>                (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                case NPY_UINT32:      bandwidth_strided_scalar<npy_uint32>                (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                case NPY_UINT64:      bandwidth_strided_scalar<npy_uint64>                (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                case NPY_FLOAT:       bandwidth_strided_scalar<float>                     (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                case NPY_DOUBLE:      bandwidth_strided_scalar<double>                    (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                case NPY_COMPLEX64:   bandwidth_strided_scalar<std::complex<float>>       (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                case NPY_COMPLEX128:  bandwidth_strided_scalar<std::complex<double>>      (a_data, offset, n, m, s0, s1, &lb_data[idx], &ub_data[idx]); break;
                default:
                    Py_DECREF(ap_lb); Py_DECREF(ap_ub);
                    PyErr_SetString(PyExc_TypeError, "Unsupported dtype.");
                    return NULL;
            }
        }
    }

    // 2D: return tuple of np.int64 scalars; N-d: return tuple of int64 arrays
    if (ndim == 2) {
        PyArray_Descr *descr = PyArray_DescrFromType(NPY_INT64);
        PyObject *lb_obj = PyArray_Scalar(&lb_data[0], descr, NULL);
        PyObject *ub_obj = PyArray_Scalar(&ub_data[0], descr, NULL);
        Py_DECREF(descr);
        Py_DECREF(ap_lb); Py_DECREF(ap_ub);
        return Py_BuildValue("(NN)", lb_obj, ub_obj);
    }

    return Py_BuildValue("(NN)", (PyObject *)ap_lb, (PyObject *)ap_ub);
}


static char doc_det[] = (
    "_linalg_det(a, overwrite_a, /)\n"
    "\n"
    "Compute the determinant of a square matrix via LU factorization.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "a : (..., N, N) ndarray\n"
    "    Input array of type float32, float64, complex64, or complex128.\n"
    "overwrite_a : bool\n"
    "    If True and the input is 2D contiguous and writable, the input\n"
    "    buffer is used directly as the getrf workspace (destroyed on exit).\n"
    "\n"
    "Returns\n"
    "-------\n"
    "det : (...) ndarray\n"
    "    Determinant values with the same dtype as the input.\n"
);

static char doc_lu[] = (
    "_linalg_lu(a, permute_l, overwrite_a, /)\n"
    "\n"
    "LU factorization with partial pivoting.\n"
    "\n"
    "Computes P, L, U such that ``A = P @ L @ U`` where P is a permutation,\n"
    "L is unit lower triangular, and U is upper triangular.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "a : (..., M, N) ndarray\n"
    "    Input array of type float32, float64, complex64, or complex128.\n"
    "permute_l : bool\n"
    "    If True, L is returned already permuted (P @ L) and P is empty.\n"
    "overwrite_a : bool\n"
    "    If True and the input is 2D, F-contiguous, and writable, the input\n"
    "    buffer is reused as the larger factor (destroyed on exit).\n"
    "\n"
    "Returns\n"
    "-------\n"
    "p : (..., M) ndarray of int\n"
    "    Permutation indices. Empty array when permute_l is True.\n"
    "l : (..., M, K) ndarray\n"
    "    Lower triangular factor with unit diagonal, K = min(M, N).\n"
    "u : (..., K, N) ndarray\n"
    "    Upper triangular factor.\n"
);

static char doc_bandwidth[] = (
    "_bandwidth(a, /)\n"
    "\n"
    "Compute the lower and upper bandwidth of a 2D array.\n"
    "\n"
    "Scans the array for the outermost nonzero entries below and above the\n"
    "main diagonal. For C- or F-contiguous arrays an AVX2-accelerated path\n"
    "is used when available; non-contiguous arrays fall back to a strided\n"
    "scalar scan.\n"
    "\n"
    "Parameters\n"
    "----------\n"
    "a : (N, M) ndarray\n"
    "    Input array of type float32, float64, complex64, or complex128.\n"
    "\n"
    "Returns\n"
    "-------\n"
    "lower : int\n"
    "    Lower bandwidth. 0 means no subdiagonal entries, N-1 means full.\n"
    "upper : int\n"
    "    Upper bandwidth. 0 means no superdiagonal entries, M-1 means full.\n"
);
static char doc_inv[] = ("Compute the matrix inverse.");
static char doc_solve[] = ("Solve the linear system of equations.");
static char doc_svd[] = ("SVD factorization.");
static char doc_lstsq[] = ("linear least squares.");
static char doc_eig[] = ("eigenvalue solver.");
static char doc_cholesky[] = ("Cholesky factorization.");
static char doc_qr[] = ("Compute the qr decomposition.");

static struct PyMethodDef module_methods[] = {
  {"_bandwidth", _linalg_bandwidth, METH_VARARGS, doc_bandwidth},
  {"_det", _linalg_det, METH_VARARGS, doc_det},
  {"_lu", _linalg_lu, METH_VARARGS, doc_lu},
  {"_inv", _linalg_inv, METH_VARARGS, doc_inv},
  {"_solve", _linalg_solve, METH_VARARGS, doc_solve},
  {"_svd", _linalg_svd, METH_VARARGS, doc_svd},
  {"_lstsq", _linalg_lstsq, METH_VARARGS, doc_lstsq},
  {"_eig", _linalg_eig, METH_VARARGS, doc_eig},
  {"_cholesky", _linalg_cholesky, METH_VARARGS, doc_cholesky},
  {"_qr", _linalg_qr, METH_VARARGS, doc_qr},
  {NULL, NULL, 0, NULL}
};


static int module_exec(PyObject *module) {
    if (_import_array() < 0) { return -1; }

    _linalg_inv_error = PyErr_NewException("_batched_linalg.error", NULL, NULL);
    if (_linalg_inv_error == NULL) { return -1; }
    if (PyModule_AddObject(module, "error", _linalg_inv_error) < 0) {
        Py_DECREF(_linalg_inv_error);
        return -1;
    }

    return 0;
}

static PyModuleDef_Slot module_slots[] = {
    {Py_mod_exec, (void *)module_exec},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL}
};



// No designated initializers under /std:c++17.
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,                        /* m_base */
    "_batched_linalg",                            /* m_name */
    "Linear algebra extension module for SciPy",  /* m_doc */
    0,                                            /* m_size */
    module_methods,                               /* m_methods */
    module_slots,                                 /* m_slots */
    NULL,                                         /* m_traverse */
    NULL,                                         /* m_clear */
    NULL                                          /* m_free */
};

PyMODINIT_FUNC
PyInit__batched_linalg(void)
{
    return PyModuleDef_Init(&moduledef);
}
