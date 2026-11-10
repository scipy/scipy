#include "Python.h"

#include <cstring>
#include <complex>
#include <vector>
#include <tuple>
#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

#include "scipy_blas_defines.h"
#include "lapack_calls.hpp"
#include "wrapper_types.hpp"
#include "wrapper_helpers.hpp"   /* wrapper::py_ref */
#include "_common_array_utils.hh"

#include "_linalg_cholesky.hh"
#include "_linalg_eig.hh"
#include "_linalg_inv.hh"
#include "_linalg_lstsq.hh"
#include "_linalg_lu_det.hh"
#include "_linalg_qr.hh"
#include "_linalg_solve.hh"
#include "_linalg_svd.hh"


using namespace sp_linalg;

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


/**
 * @brief Raise the Python exception for a failed batched call.
 *
 * @param status  The failure; any value except BatchStatus::ok.
 * @param func    Public function name used in the message, e.g. "solve".
 * @return NULL, with the exception set.
 */
static PyObject*
set_batch_error(BatchStatus status, const char *func) {
    switch (status) {
        case BatchStatus::no_memory:
            return PyErr_Format(PyExc_MemoryError, "scipy.linalg.%s: could not allocate memory.", func);
        case BatchStatus::too_large:
            return PyErr_Format(PyExc_ValueError,
                "scipy.linalg.%s: the array is too large for the %d-bit integers of the LAPACK "
                "library; build SciPy with ILP64 support.", func, (int)(8 * sizeof(CBLAS_INT)));
        case BatchStatus::lapack_query_failed:
            return PyErr_Format(PyExc_RuntimeError, "scipy.linalg.%s: the LAPACK workspace query failed.", func);
        default:
            return PyErr_Format(PyExc_TypeError, "scipy.linalg.%s: unsupported array dtype.", func);
    }
}


static PyObject*
_linalg_inv(PyObject* Py_UNUSED(dummy), PyObject* args) {

    PyArrayObject* ap_Am = NULL;
    SliceStatusVec vec_status;
    St structure = St::NONE;
    int overwrite_a = 0;
    int lower = 0;

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

    // The output: a new array, or the input array itself when overwriting it
    wrapper::py_ref ap_Ainv;
    if (overwrite_a) {
        Py_INCREF(ap_Am);
        ap_Ainv = wrapper::py_ref(ap_Am);
    }
    else {
        ap_Ainv = wrapper::py_ref((PyArrayObject *)PyArray_SimpleNew(ndim, shape, typenum));
        if (!ap_Ainv) { return PyErr_NoMemory(); }
    }

    BatchStatus status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _inverse<T>(ap_Am, ap_Ainv.data<T>(), structure, lower, overwrite_a, vec_status);
    });
    if (status != BatchStatus::ok) { return set_batch_error(status, "inv"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    return Py_BuildValue("NN", PyArray_Return((PyArrayObject *)ap_Ainv.release()), ret_lst);
}



static PyObject*
_linalg_solve(PyObject* Py_UNUSED(dummy), PyObject* args) {

    PyArrayObject *ap_Am = NULL;
    PyArrayObject *ap_b = NULL;
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

    // The output: a new array, or the input array `b` itself when overwriting it
    wrapper::py_ref ap_x;
    if (overwrite_b) {
        Py_INCREF(ap_b);
        ap_x = wrapper::py_ref(ap_b);
    }
    else {
        ap_x = wrapper::py_ref((PyArrayObject *)PyArray_SimpleNew(ndim_b, shape_b, typenum));
        if (!ap_x) { return PyErr_NoMemory(); }
    }

    BatchStatus status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _solve<T>(ap_Am, ap_b, ap_x.data<T>(), structure, lower, transposed, overwrite_a, overwrite_b, vec_status);
    });
    if (status != BatchStatus::ok) { return set_batch_error(status, "solve"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    return Py_BuildValue("NN", PyArray_Return((PyArrayObject *)ap_x.release()), ret_lst);
}



static PyObject*
_linalg_qr(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_A = NULL;

    SliceStatusVec vec_status;
    int overwrite_a = 0;
    QR_mode mode = QR_mode::FULL;
    int pivoting = 0;

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
    wrapper::py_ref ap_Q, ap_R, ap_tau, ap_jpvt;

    if (mode != QR_mode::R) {
        if (!overwrite_a || mode == QR_mode::FULL) {
            ap_Q = wrapper::py_ref((PyArrayObject *)PyArray_SimpleNew(ndim, shape_Q, typenum));
            if (!ap_Q) { return PyErr_NoMemory(); }
        } else {
            Py_INCREF(ap_A);
            ap_Q = wrapper::py_ref(ap_A);
        }
    }

    if (mode == QR_mode::RAW_MODE || mode == QR_mode::ECONOMIC || (!overwrite_a && (mode == QR_mode::R || mode == QR_mode::FULL))) {
        // Return C-ordered object, hence `0` flag. Already pre-fill with zeros for efficiency.
        ap_R = wrapper::py_ref((PyArrayObject *)PyArray_ZEROS(ndim, shape_R, typenum, 0));
        if (!ap_R) { return PyErr_NoMemory(); }
    } else {
        Py_INCREF(ap_A);
        ap_R = wrapper::py_ref(ap_A);
    }

    if (mode == QR_mode::RAW_MODE) {
        shape_Q[ndim-2] = K; // Just reuse `shape_Q`; not used any longer.
        ap_tau = wrapper::py_ref((PyArrayObject *)PyArray_SimpleNew(ndim-1, shape_Q, typenum)); // Just a vector, so `ndim-1`.
        if (!ap_tau) { return PyErr_NoMemory(); }

        /*
         * Since the point of the `Q` returned if `mode="raw"` is mostly to call other LAPACK routines (e.g. `qr_multiply`),
         * set the strides of `Q` to return F-ordered slices.
         *
         * Dimensions other than the last two are not affected so slice computation etc. stays intact.
         *
         * Since `overwrite_a` is only enabled when the input has F-ordered slices, this is not necessary in that case.
         */
        if (!overwrite_a) {
            npy_intp *strides_Q = PyArray_STRIDES(ap_Q.get());
            int sizeof_T = PyArray_ITEMSIZE(ap_Q.get());

            strides_Q[ndim-2] = sizeof_T;
            strides_Q[ndim-1] = M * sizeof_T;

            PyArray_UpdateFlags(ap_Q.get(), NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_F_CONTIGUOUS);
        }
    }

    // Set all elements to 0 immediately as the pivots are also used as inputs in `geqp3`.
    // Allocate a C-ordered array, (hence the `0` magic number).
    if (pivoting) {
        shape_Q[ndim-2] = N;
        ap_jpvt = wrapper::py_ref((PyArrayObject *)PyArray_ZEROS(ndim-1, shape_Q, sizeof(CBLAS_INT) == sizeof(NPY_INT32)? NPY_INT32 : NPY_INT64, 0));
        if (!ap_jpvt) { return PyErr_NoMemory(); }
    }

    BatchStatus status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _qr<T>(ap_A, ap_Q.get(), ap_R.get(), ap_tau.get(), ap_jpvt.get(), overwrite_a, mode, pivoting, vec_status);
    });
    if (status != BatchStatus::ok) { return set_batch_error(status, "qr"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    PyObject *ret_Q = (mode != QR_mode::R) ? PyArray_Return((PyArrayObject *)ap_Q.release()) : Py_NewRef(Py_None);
    PyObject *ret_tau = (mode == QR_mode::RAW_MODE) ? PyArray_Return((PyArrayObject *)ap_tau.release()) : Py_NewRef(Py_None);
    PyObject *ret_jpvt = pivoting ? PyArray_Return((PyArrayObject *)ap_jpvt.release()) : Py_NewRef(Py_None);
    return Py_BuildValue("NNNNN", ret_Q, PyArray_Return((PyArrayObject *)ap_R.release()), ret_tau, ret_jpvt, ret_lst);
}



static PyObject*
_linalg_svd(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_Am = NULL;
    const char *lapack_driver = NULL;
    int compute_uv = 1;
    int full_matrices = 1;
    int overwrite_a = 0;
    SliceStatusVec vec_status;

    // Get the input array
    if (!PyArg_ParseTuple(args, "O!s|ppp", &PyArray_Type, (PyObject **)&ap_Am,  &lapack_driver, &compute_uv, &full_matrices, &overwrite_a)) {
        return NULL;
    }

    // Sanity check the input array
    if (!_check_dtype_and_flags(ap_Am, "svd")) {
        return NULL;
    }

    bool use_gesdd = (strcmp(lapack_driver, "gesdd") == 0);
    if (!use_gesdd && (strcmp(lapack_driver, "gesvd") != 0)) {
        return PyErr_Format(PyExc_ValueError, "scipy.linalg.svd: unknown lapack_driver '%s'.", lapack_driver);
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
    wrapper::py_ref ap_S((PyArrayObject *)PyArray_SimpleNew(ndim-1, shape_1, typenum_S));
    if (!ap_S) { return PyErr_NoMemory(); }

    char jobz = compute_uv ? (full_matrices ? 'A' : 'S') : 'N';

    // U.shape = (..., m, m) or (..., m, k) or not referenced (for jobz='N')
    // Vh.shape = (..., n, n) or (..., k, n) or not referenced (for jobz='N')
    wrapper::py_ref ap_U, ap_Vh;
    if (jobz != 'N') {
        npy_intp u_shape0, u_shape1, vh_shape0, vh_shape1;
        u_vh_shapes(m, n, jobz, &u_shape0, &u_shape1, &vh_shape0, &vh_shape1);

        shape_1[ndim-2] = u_shape0;
        shape_1[ndim-1] = u_shape1;
        ap_U = wrapper::py_ref((PyArrayObject *)PyArray_SimpleNew(ndim, shape_1, typenum));
        if (!ap_U) { return PyErr_NoMemory(); }

        shape_1[ndim-2] = vh_shape0;
        shape_1[ndim-1] = vh_shape1;
        ap_Vh = wrapper::py_ref((PyArrayObject *)PyArray_SimpleNew(ndim, shape_1, typenum));
        if (!ap_Vh) { return PyErr_NoMemory(); }
    }

    BatchStatus status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _svd<T>(ap_Am, ap_U.get(), ap_S.get(), ap_Vh.get(), jobz, use_gesdd, overwrite_a, vec_status);
    });
    if (status != BatchStatus::ok) { return set_batch_error(status, "svd"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    if (compute_uv) {
        return Py_BuildValue("NNNN", PyArray_Return((PyArrayObject *)ap_U.release()),
                                     PyArray_Return((PyArrayObject *)ap_S.release()),
                                     PyArray_Return((PyArrayObject *)ap_Vh.release()), ret_lst);
    }
    return Py_BuildValue("NN", PyArray_Return((PyArrayObject *)ap_S.release()), ret_lst);
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
            info = _lstsq<f32>(ap_Am, ap_b, ap_S, ap_x, ap_rank, (float)rcond, lapack_driver, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _lstsq<f64>(ap_Am, ap_b, ap_S, ap_x, ap_rank, rcond, lapack_driver, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _lstsq<c64>(ap_Am, ap_b, ap_S, ap_x, ap_rank, (float)rcond, lapack_driver, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _lstsq<c128>(ap_Am, ap_b, ap_S, ap_x, ap_rank, rcond, lapack_driver, overwrite_a, overwrite_b, vec_status);
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
            info = _eig<f32>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_FLOAT64):
            info = _eig<f64>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_COMPLEX64):
            info = _eig<c64>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, overwrite_a, overwrite_b, vec_status);
            break;
        case(NPY_COMPLEX128):
            info = _eig<c128>(ap_Am, ap_Bm, ap_w, ap_beta, ap_vl, ap_vr, overwrite_a, overwrite_b, vec_status);
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
_linalg_eigh(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_Am = NULL;
    PyArrayObject *ap_Bm = NULL;
    PyArrayObject *ap_w = NULL;
    PyArrayObject *ap_Z = NULL;
    int overwrite_a = 0;
    int overwrite_b = 0;
    int eigvals_only = 0;
    int lower = 0;
    int vals_range = 0;
    double vl = -std::numeric_limits<double>::infinity();
    double vu = std::numeric_limits<double>::infinity();
    int il = -1;
    int iu = -1;
    Eigh_driver lapack_driver;
    int itype = 1;

    int info = 0;
    SliceStatusVec vec_status;

    // Return objects
    PyObject *ret_w;
    PyObject *ret_Z;
    PyObject *ret_lst;

    if (!PyArg_ParseTuple(args, "O!ppipddiin|O!pi",
        &PyArray_Type, (PyObject **)&ap_Am, &overwrite_a,
        &eigvals_only, &vals_range, &lower,
        &vl, &vu, &il, &iu, &lapack_driver,
        &PyArray_Type, (PyArrayObject **)&ap_Bm, &overwrite_b, &itype)
    ) {
        return NULL;
    }

    /*
     * Convert int input to char for LAPACK
     *
     * - `uplo`: what triangular part of the matrix to use
     * - `jobz`: whether to compute eigenvectors, always applicable
     * - `range`: for drivers `evr`, `evx` and `gvx` what eigenvalues to compute:
     *      - 'A':  all eigenvalues, `N` in total, are returned unconditionally,
     *              `il`, `iu`, `vl` and `vu` are still sentinel values
     *      - 'V':  only eigenvalues `w_i` satisfying `vl < w_i <= vu` are returned, since the
     *              number of eigenvalues satisfying this condition is unknown beforehand,
     *              space for all `N` eigenvalues has to be allocated.
     *              `il` and `iu` remain sentinel values.
     *      - 'I':  when sorted in ascending order, return indices `il` to `iu` (inclusive),
     *              not yet accounted for Fortran indexing, the number of eigenvalues is known
     *              beforehand to be `iu - il + 1`, hence only that size should be allocated.
     *              `vl` and `vu` remain sentinels
     *
     * The number of eigenvalues/eigenvectors to be found is denoted with `M`. In order to be
     * able to return a properly sized array this is also returned to the python side.
     */
    char uplo = lower ? 'L' : 'U';
    char jobz = eigvals_only ? 'N' : 'V';
    char range = (vals_range == 0) ? 'A' : (vals_range == 1) ? 'V' : 'I';

    // Sanity check `a`
    if (!_check_dtype_and_flags(ap_Am, "eigh")) {
        return NULL;
    }


    /*
     * Allocation of return objects
     *
     * - `w` contains eigenvalues, always present, size (N,) per slice
     * - `Z` contains eigenvectors, only required when `jobz == 'V'`, size (N, M) per slice
     */
    int typenum = PyArray_TYPE(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    npy_intp N = shape[ndim - 1];
    if (shape[ndim - 2] != N) {
        PyErr_SetString(PyExc_ValueError, "Expected a square matrix");
        return NULL;
    }

    // `evr` does not return in place, so a new output array is always needed for `Z`
    int w_typenum = typenum;
    if (typenum == NPY_COMPLEX64) { w_typenum = NPY_FLOAT32; }
    if (typenum == NPY_COMPLEX128) { w_typenum = NPY_FLOAT64; }

    npy_intp shape_1[NPY_MAXDIMS];
    for (int i = 0; i < ndim ; i++) { shape_1[i] = shape[i]; }

    int M = (range == 'I') ? iu + 1 - il : N;
    shape_1[ndim - 2] = M;
    ap_w = (PyArrayObject *)PyArray_SimpleNew(ndim - 1, shape_1, w_typenum);
    if (ap_w == NULL) {
        PyErr_NoMemory();
        goto fail;
    }

    if (jobz != 'N') {
        if (overwrite_a) {
            Py_INCREF(ap_Am);
            ap_Z = ap_Am;
        } else {
            shape_1[ndim - 2] = N;
            shape_1[ndim - 1] = M;

            ap_Z = (PyArrayObject *)PyArray_SimpleNew(ndim, shape_1, typenum);
            if (ap_Z == NULL) {
                PyErr_NoMemory();
                goto fail;
            }
        }
    }

    // Dispatch to actual implementation
    // Pass in pointer to `M` to be able to return to python side
    switch (typenum) {
        case NPY_FLOAT32:
            info = _eigh<f32>(ap_Am, ap_Bm, ap_w, ap_Z, &M, overwrite_a, overwrite_b, itype, jobz, range, uplo, vl, vu, il, iu, lapack_driver, vec_status);
            break;
        case NPY_FLOAT64:
            info = _eigh<f64>(ap_Am, ap_Bm, ap_w, ap_Z, &M, overwrite_a, overwrite_b, itype, jobz, range, uplo, vl, vu, il, iu, lapack_driver, vec_status);
            break;
        case NPY_COMPLEX64:
            info = _eigh<c64>(ap_Am, ap_Bm, ap_w, ap_Z, &M, overwrite_a, overwrite_b, itype, jobz, range, uplo, vl, vu, il, iu, lapack_driver, vec_status);
            break;
        case NPY_COMPLEX128:
            info = _eigh<c128>(ap_Am, ap_Bm, ap_w, ap_Z, &M, overwrite_a, overwrite_b, itype, jobz, range, uplo, vl, vu, il, iu, lapack_driver, vec_status);
            break;
    }

    if (info < 0) {
        // Out-of-memory or scipy internal error
        PyErr_SetString(PyExc_RuntimeError, "Memory error in scipy.linalg.eigh.");
        goto fail;
    }

    // regular return
    ret_w = PyArray_Return(ap_w);
    ret_Z = (jobz == 'V') ? PyArray_Return(ap_Z) : Py_None;
    ret_lst = convert_vec_status(vec_status);

    return Py_BuildValue("NNnN", ret_w, ret_Z, (Py_ssize_t)M, ret_lst);

fail:
    Py_XDECREF(ap_w);
    Py_XDECREF(ap_Z);
    return NULL;
}



static PyObject*
_linalg_cholesky(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_Am = NULL;
    int lower=0, overwrite_a=0, clean=1;
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

    // The output: a new zero-filled C-ordered array, or the input array itself when overwriting it
    wrapper::py_ref ap_Cm;
    if (overwrite_a) {
        Py_INCREF(ap_Am);
        ap_Cm = wrapper::py_ref(ap_Am);
    }
    else {
        ap_Cm = wrapper::py_ref((PyArrayObject *)PyArray_ZEROS(ndim, shape, typenum, 0));
        if (!ap_Cm) { return PyErr_NoMemory(); }
    }

    BatchStatus status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _cholesky<T>(ap_Am, ap_Cm.get(), lower, overwrite_a, clean, vec_status);
    });
    if (status != BatchStatus::ok) { return set_batch_error(status, "cholesky"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    return Py_BuildValue("NN", PyArray_Return((PyArrayObject *)ap_Cm.release()), ret_lst);
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

    // Output shapes are derived from a copy; the input array's own shape is never modified
    npy_intp shape[LU_MAX_NDIM];
    std::copy(PyArray_SHAPE(ap_a), PyArray_SHAPE(ap_a) + ndim, shape);
    npy_intp *byte_strides = PyArray_STRIDES(ap_a);

    CBLAS_INT m, n;
    std::size_t scratch_count;
    if (!checked_narrow<CBLAS_INT>(shape[ndim - 2], &m) || !checked_narrow<CBLAS_INT>(shape[ndim - 1], &n) ||
        !checked_mul(m, n, &scratch_count)) {
        return set_batch_error(BatchStatus::too_large, "lu");
    }
    CBLAS_INT minmn = m < n ? m : n;

    npy_intp num_of_slices = 1;
    for (int i = 0; i < ndim - 2; i++) { num_of_slices *= shape[i]; }

    int64_t ctx_shape[LU_MAX_NDIM];
    int64_t ctx_strides[LU_MAX_NDIM];
    for (int i = 0; i < ndim; i++) {
        ctx_shape[i] = (int64_t)shape[i];
        ctx_strides[i] = (int64_t)byte_strides[i];
    }

    // perm (..., m)
    wrapper::py_ref ap_perm((PyArrayObject *)PyArray_SimpleNew(ndim - 1, shape, (sizeof(CBLAS_INT) == 8) ? NPY_INT64 : NPY_INT32));
    if (!ap_perm) { return PyErr_NoMemory(); }

    LU_Context ctx = {};
    ctx.shape = ctx_shape;
    ctx.strides = ctx_strides;
    ctx.num_of_slices = num_of_slices;
    ctx.m = m;
    ctx.n = n;
    ctx.perm = ap_perm.data<CBLAS_INT>();
    ctx.ndim = ndim;
    ctx.permute_l = (bool)permute_l;
    ctx.overwrite_a = false;

    wrapper::py_ref ap_l, ap_u;
    npy_intp l_shape[LU_MAX_NDIM], u_shape[LU_MAX_NDIM];
    std::copy(shape, shape + ndim, l_shape);
    std::copy(shape, shape + ndim, u_shape);

    // Legacy overwrite_a: reuse input as the larger factor (2D + F-contiguous only).
    // Disabled when permute_l is true (permute_rows needs C-ordered data).
    if (overwrite_a && !permute_l && (ndim == 2) && PyArray_ISFARRAY_RO(ap_a) && PyArray_ISWRITEABLE(ap_a)) {
        ctx.overwrite_a = true;
        Py_INCREF(ap_a);
        if (m > n) {
            // Tall: L (m x n) is large — reuse input. Allocate U (n x n).
            ap_l = wrapper::py_ref(ap_a);
            u_shape[ndim - 2] = minmn;
            ap_u = wrapper::py_ref((PyArrayObject *)PyArray_ZEROS(ndim, u_shape, typenum, 0));
        } else {
            // Fat/square: U (m x n) is large — reuse input. Allocate L (m x m).
            ap_u = wrapper::py_ref(ap_a);
            l_shape[ndim - 1] = minmn;
            ap_l = wrapper::py_ref((PyArrayObject *)PyArray_ZEROS(ndim, l_shape, typenum, 0));
        }
    } else {
        // Allocate output arrays: L (..., m, minmn), U (..., minmn, n)
        l_shape[ndim - 1] = minmn;
        u_shape[ndim - 2] = minmn;
        ap_l = wrapper::py_ref((PyArrayObject *)PyArray_ZEROS(ndim, l_shape, typenum, 0));
        ap_u = wrapper::py_ref((PyArrayObject *)PyArray_ZEROS(ndim, u_shape, typenum, 0));
    }
    if (!ap_l || !ap_u) { return PyErr_NoMemory(); }

    int info = 0;
    BatchStatus status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        // slice_info: per-slice getrf info (> 0: exactly singular); not yet surfaced to Python
        auto scratch = alloc_workspace<T>(ctx.overwrite_a ? 0 : scratch_count);
        auto ipiv = alloc_workspace<CBLAS_INT>(minmn);
        auto slice_info = alloc_workspace<CBLAS_INT>(num_of_slices);
        if (!scratch || !ipiv || !slice_info) { return BatchStatus::no_memory; }
        ctx.ipiv = ipiv.get();
        info = lu_dispatch<T>(ctx, static_cast<T *>(PyArray_DATA(ap_a)), ap_l.data<T>(), ap_u.data<T>(), scratch.get(), slice_info.get());
        return BatchStatus::ok;
    });
    if (status != BatchStatus::ok) { return set_batch_error(status, "lu"); }
    if (info < 0) {
        return PyErr_Format(PyExc_ValueError, "scipy.linalg.lu: internal error, getrf returned info = %d.", info);
    }

    // Always return (P, L, U) with consistent types.
    // When permute_l is true, P is replaced with an empty array since
    // the permutation is already applied to L.
    if (permute_l) {
        // Return an empty int32 array instead of None so the return type is
        // always (ndarray, ndarray, ndarray) — avoids Optional in type stubs.
        // int32 since the array has zero elements; no overflow concern.
        npy_intp zero = 0;
        ap_perm = wrapper::py_ref((PyArrayObject *)PyArray_EMPTY(1, &zero, NPY_INT32, 0));
        if (!ap_perm) { return PyErr_NoMemory(); }
    }
    return Py_BuildValue("NNN", PyArray_Return((PyArrayObject *)ap_perm.release()),
                                PyArray_Return((PyArrayObject *)ap_l.release()),
                                PyArray_Return((PyArrayObject *)ap_u.release()));
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

    if (shape[ndim - 2] != shape[ndim - 1]) {
        PyErr_SetString(PyExc_ValueError, "scipy.linalg.det: Last two dimensions must be square.");
        return NULL;
    }

    CBLAS_INT n;
    std::size_t scratch_count;
    if (!checked_narrow<CBLAS_INT>(shape[ndim - 1], &n) || !checked_mul(n, n, &scratch_count)) {
        return set_batch_error(BatchStatus::too_large, "det");
    }

    npy_intp num_of_slices = 1;
    for (int i = 0; i < ndim - 2; i++) { num_of_slices *= shape[i]; }

    int64_t ctx_shape[LU_MAX_NDIM];
    int64_t ctx_strides[LU_MAX_NDIM];
    for (int i = 0; i < ndim; i++) {
        ctx_shape[i] = (int64_t)shape[i];
        ctx_strides[i] = (int64_t)byte_strides[i];
    }

    // Allocate output det array (...,) — batch dimensions only
    wrapper::py_ref ap_det((PyArrayObject *)PyArray_SimpleNew(ndim - 2, shape, typenum));
    if (!ap_det) { return PyErr_NoMemory(); }

    LU_Context ctx = {};
    ctx.shape = ctx_shape;
    ctx.strides = ctx_strides;
    ctx.num_of_slices = num_of_slices;
    ctx.m = n;
    ctx.n = n;
    ctx.perm = NULL;
    ctx.ndim = ndim;
    ctx.permute_l = false;
    ctx.overwrite_a = false;

    // Legacy overwrite_a: use input directly (2D + contiguous only).
    // For det, both C and F order work since det(A^T) = det(A).
    if (overwrite_a && (ndim == 2) && PyArray_ISWRITEABLE(ap_a)
        && (PyArray_ISFARRAY_RO(ap_a) || PyArray_ISCARRAY_RO(ap_a))) {
        ctx.overwrite_a = true;
    }

    int info = 0;
    BatchStatus status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        // slice_info: per-slice getrf info (> 0: exactly singular, so the determinant is exactly
        // zero rather than numerically small); not yet surfaced to Python
        auto scratch = alloc_workspace<T>(ctx.overwrite_a ? 0 : scratch_count);
        auto ipiv = alloc_workspace<CBLAS_INT>(n);
        auto slice_info = alloc_workspace<CBLAS_INT>(num_of_slices);
        if (!scratch || !ipiv || !slice_info) { return BatchStatus::no_memory; }
        ctx.ipiv = ipiv.get();
        info = det_dispatch<T>(ctx, static_cast<T *>(PyArray_DATA(ap_a)), ap_det.data<T>(), scratch.get(), slice_info.get());
        return BatchStatus::ok;
    });
    if (status != BatchStatus::ok) { return set_batch_error(status, "det"); }
    if (info < 0) {
        return PyErr_Format(PyExc_ValueError, "scipy.linalg.det: internal error, getrf returned info = %d.", info);
    }

    return ap_det.release();
}


// ========================================================================
// Bandwidth
// ========================================================================

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
    npy_intp *strides = PyArray_STRIDES(ap_a);
    int64_t n = shape[ndim - 2], m = shape[ndim - 1];
    int64_t s0 = strides[ndim - 2], s1 = strides[ndim - 1];

    // Normalize platform-dependent integer typenums to canonical
    // case labels. On Windows sizeof(int)==sizeof(long)==4 but
    // NPY_INT != NPY_LONG, and NPY_INT32 is only one of them.
    // Map all integer typenums to canonical NPY_INTxx by itemsize.
    int typenum = PyArray_TYPE(ap_a);
    if (PyTypeNum_ISINTEGER(typenum) && !PyTypeNum_ISBOOL(typenum)) {
        bool is_unsigned = PyTypeNum_ISUNSIGNED(typenum);
        switch (PyArray_ITEMSIZE(ap_a)) {
            case 1: typenum = is_unsigned ? NPY_UINT8  : NPY_INT8;  break;
            case 2: typenum = is_unsigned ? NPY_UINT16 : NPY_INT16; break;
            case 4: typenum = is_unsigned ? NPY_UINT32 : NPY_INT32; break;
            case 8: typenum = is_unsigned ? NPY_UINT64 : NPY_INT64; break;
        }
    }

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
    const char *a_data = static_cast<const char *>(PyArray_DATA(ap_a));

    auto each_slice = [&](auto tag) {
        using T = decltype(tag);
        for (npy_intp idx = 0; idx < num_slices; idx++) {
            // byte offset of the idx-th 2D slice
            int64_t offset = 0, temp = idx;
            for (int i = batch_ndim - 1; i >= 0; i--) {
                offset += (temp % shape[i]) * strides[i];
                temp /= shape[i];
            }
            bandwidth_slice(reinterpret_cast<const T *>(a_data + offset), n, m, s0, s1, &lb_data[idx], &ub_data[idx]);
        }
    };

    // longdouble/clongdouble are rejected by the Python wrapper in _misc.py
    switch (typenum) {
        case NPY_BOOL:        each_slice(npy_bool{});   break;
        case NPY_INT8:        each_slice(npy_int8{});   break;
        case NPY_INT16:       each_slice(npy_int16{});  break;
        case NPY_INT32:       each_slice(npy_int32{});  break;
        case NPY_INT64:       each_slice(npy_int64{});  break;
        case NPY_UINT8:       each_slice(npy_uint8{});  break;
        case NPY_UINT16:      each_slice(npy_uint16{}); break;
        case NPY_UINT32:      each_slice(npy_uint32{}); break;
        case NPY_UINT64:      each_slice(npy_uint64{}); break;
        case NPY_FLOAT:       each_slice(f32{});        break;
        case NPY_DOUBLE:      each_slice(f64{});        break;
        case NPY_COMPLEX64:   each_slice(c64{});        break;
        case NPY_COMPLEX128:  each_slice(c128{});       break;
        default:
            Py_DECREF(ap_lb); Py_DECREF(ap_ub);
            PyErr_SetString(PyExc_TypeError, "Unsupported dtype.");
            return NULL;
    }

    // 2D: return tuple of int scalars; N-d: return tuple of int64 arrays
    if (ndim == 2) {
        // TODO: use PyLong_FromInt64 from py314+
        PyObject *lb_obj = PyLong_FromSsize_t((Py_ssize_t)lb_data[0]);
        PyObject *ub_obj = PyLong_FromSsize_t((Py_ssize_t)ub_data[0]);
        Py_DECREF(ap_lb); Py_DECREF(ap_ub);
        return Py_BuildValue("NN", lb_obj, ub_obj);
    }

    return Py_BuildValue("NN", (PyObject *)ap_lb, (PyObject *)ap_ub);
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
static char doc_eigh[] = ("eigenvalue solver for symmetric/hermitian matrices.");
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
  {"_eigh", _linalg_eigh, METH_VARARGS, doc_eigh},
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
