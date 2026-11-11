#include "Python.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <tuple>
#include <type_traits>
#include <vector>
#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

#include "scipy_blas_defines.h"
#include "lapack_calls.hpp"
#include "wrapper_types.hpp"
#include "_linalg_bandwidth.h"
#include "_common_array_utils.hpp"

#include "_linalg_cholesky.hpp"
#include "_linalg_eig.hpp"
#include "_linalg_inv.hpp"
#include "_linalg_lstsq.hpp"
#include "_linalg_lu_det.hpp"
#include "_linalg_qr.hpp"
#include "_linalg_solve.hpp"
#include "_linalg_svd.hpp"


using namespace sp_linalg;

static PyObject* _linalg_inv_error;


PyObject*
convert_vec_status(SliceStatusVec& vec_status);


/*
 * Minimal checks for an input array:
 *   - ndim >= 2
 *   - LAPACK-compatible dtype
 *   - aligned, in native byte order
 *
 * Sets an error if any of these is violated (users to bail out immediately).
 */
int _check_dtype_and_flags(PyArrayObject *ap_Am, const char *func) {
    int typenum = PyArray_TYPE(ap_Am);
    bool dtype_ok = (typenum == NPY_FLOAT32)
                     || (typenum == NPY_FLOAT64)
                     || (typenum == NPY_COMPLEX64)
                     || (typenum == NPY_COMPLEX128);

    if (!dtype_ok) {
        PyErr_Format(PyExc_TypeError, "scipy.linalg.%s: Expected an array of float32, float64, "
                     "complex64 or complex128.", func);
        return 0;
    }
    if (!PyArray_ISBEHAVED_RO(ap_Am)) {
        PyErr_Format(PyExc_TypeError, "scipy.linalg.%s: Expected an aligned array in native "
                     "byte order.", func);
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
                "scipy.linalg.%s: a dimension or workspace size exceeds the range of the %d-bit "
                "integers of the LAPACK library.", func, (int)(8 * sizeof(CBLAS_INT)));
        case BatchStatus::lapack_error:
            return PyErr_Format(PyExc_RuntimeError,
                "scipy.linalg.%s: a LAPACK routine rejected one of its arguments.", func);
        case BatchStatus::bad_dtype:
            return PyErr_Format(PyExc_TypeError, "scipy.linalg.%s: unsupported array dtype.", func);
        case BatchStatus::ok:
            break;
    }
    return PyErr_Format(PyExc_SystemError, "scipy.linalg.%s: no error to report.", func);
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
    OwnedArray ap_Ainv;
    if (overwrite_a) {
        Py_INCREF(ap_Am);
        ap_Ainv.reset(ap_Am);
    }
    else {
        ap_Ainv.reset((PyArrayObject *)PyArray_SimpleNew(ndim, shape, typenum));
        if (!ap_Ainv) { return NULL; }
    }

    BatchStatus status;
    Py_BEGIN_ALLOW_THREADS
    status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _inverse<T>(ap_Am, static_cast<T *>(PyArray_DATA(ap_Ainv.get())), structure, lower, overwrite_a, vec_status);
    });
    Py_END_ALLOW_THREADS
    if (status != BatchStatus::ok) { return set_batch_error(status, "inv"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    if (!ret_lst) { return NULL; }
    return Py_BuildValue("NN", PyArray_Return(ap_Ainv.release()), ret_lst);
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
    if (!_check_dtype_and_flags(ap_b, "solve")) { return NULL; }

    // Sanity check shapes
    int typenum = PyArray_TYPE(ap_Am);
    if (PyArray_TYPE(ap_b) != typenum) {
        PyErr_SetString(PyExc_ValueError, "`a` and `b` must have the same dtype.");
        return NULL;
    }
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
    OwnedArray ap_x;
    if (overwrite_b) {
        Py_INCREF(ap_b);
        ap_x.reset(ap_b);
    }
    else {
        ap_x.reset((PyArrayObject *)PyArray_SimpleNew(ndim_b, shape_b, typenum));
        if (!ap_x) { return NULL; }
    }

    BatchStatus status;
    Py_BEGIN_ALLOW_THREADS
    status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _solve<T>(ap_Am, ap_b, static_cast<T *>(PyArray_DATA(ap_x.get())), structure, lower, transposed, overwrite_a, overwrite_b, vec_status);
    });
    Py_END_ALLOW_THREADS
    if (status != BatchStatus::ok) { return set_batch_error(status, "solve"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    if (!ret_lst) { return NULL; }
    return Py_BuildValue("NN", PyArray_Return(ap_x.release()), ret_lst);
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
    OwnedArray ap_Q, ap_R, ap_tau, ap_jpvt;

    if (mode != QR_mode::R) {
        if (!overwrite_a || mode == QR_mode::FULL) {
            ap_Q.reset((PyArrayObject *)PyArray_SimpleNew(ndim, shape_Q, typenum));
            if (!ap_Q) { return NULL; }
        } else {
            Py_INCREF(ap_A);
            ap_Q.reset(ap_A);
        }
    }

    if (mode == QR_mode::RAW_MODE || mode == QR_mode::ECONOMIC || (!overwrite_a && (mode == QR_mode::R || mode == QR_mode::FULL))) {
        // Return C-ordered object, hence `0` flag. Already pre-fill with zeros for efficiency.
        ap_R.reset((PyArrayObject *)PyArray_ZEROS(ndim, shape_R, typenum, 0));
        if (!ap_R) { return NULL; }
    } else {
        Py_INCREF(ap_A);
        ap_R.reset(ap_A);
    }

    if (mode == QR_mode::RAW_MODE) {
        shape_Q[ndim-2] = K; // Just reuse `shape_Q`; not used any longer.
        ap_tau.reset((PyArrayObject *)PyArray_SimpleNew(ndim-1, shape_Q, typenum)); // Just a vector, so `ndim-1`.
        if (!ap_tau) { return NULL; }

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
        ap_jpvt.reset((PyArrayObject *)PyArray_ZEROS(ndim-1, shape_Q, (sizeof(CBLAS_INT) == 8) ? NPY_INT64 : NPY_INT32, 0));
        if (!ap_jpvt) { return NULL; }
    }

    BatchStatus status;
    Py_BEGIN_ALLOW_THREADS
    status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _qr<T>(ap_A, ap_Q.get(), ap_R.get(), ap_tau.get(), ap_jpvt.get(), overwrite_a, mode, pivoting, vec_status);
    });
    Py_END_ALLOW_THREADS
    if (status != BatchStatus::ok) { return set_batch_error(status, "qr"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    if (!ret_lst) { return NULL; }
    PyObject *ret_Q = (mode != QR_mode::R) ? PyArray_Return(ap_Q.release()) : Py_NewRef(Py_None);
    PyObject *ret_tau = (mode == QR_mode::RAW_MODE) ? PyArray_Return(ap_tau.release()) : Py_NewRef(Py_None);
    PyObject *ret_jpvt = pivoting ? PyArray_Return(ap_jpvt.release()) : Py_NewRef(Py_None);
    return Py_BuildValue("NNNNN", ret_Q, PyArray_Return(ap_R.release()), ret_tau, ret_jpvt, ret_lst);
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
    OwnedArray ap_S((PyArrayObject *)PyArray_SimpleNew(ndim-1, shape_1, typenum_S));
    if (!ap_S) { return NULL; }

    char jobz = compute_uv ? (full_matrices ? 'A' : 'S') : 'N';

    // U.shape = (..., m, m) or (..., m, k) or not referenced (for jobz='N')
    // Vh.shape = (..., n, n) or (..., k, n) or not referenced (for jobz='N')
    OwnedArray ap_U, ap_Vh;
    if (jobz != 'N') {
        npy_intp u_shape0, u_shape1, vh_shape0, vh_shape1;
        u_vh_shapes(m, n, jobz, &u_shape0, &u_shape1, &vh_shape0, &vh_shape1);

        shape_1[ndim-2] = u_shape0;
        shape_1[ndim-1] = u_shape1;
        ap_U.reset((PyArrayObject *)PyArray_SimpleNew(ndim, shape_1, typenum));
        if (!ap_U) { return NULL; }

        shape_1[ndim-2] = vh_shape0;
        shape_1[ndim-1] = vh_shape1;
        ap_Vh.reset((PyArrayObject *)PyArray_SimpleNew(ndim, shape_1, typenum));
        if (!ap_Vh) { return NULL; }
    }

    BatchStatus status;
    Py_BEGIN_ALLOW_THREADS
    status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _svd<T>(ap_Am, ap_U.get(), ap_S.get(), ap_Vh.get(), jobz, use_gesdd, overwrite_a, vec_status);
    });
    Py_END_ALLOW_THREADS
    if (status != BatchStatus::ok) { return set_batch_error(status, "svd"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    if (!ret_lst) { return NULL; }
    if (compute_uv) {
        return Py_BuildValue("NNNN", PyArray_Return(ap_U.release()),
                                     PyArray_Return(ap_S.release()),
                                     PyArray_Return(ap_Vh.release()), ret_lst);
    }
    return Py_BuildValue("NN", PyArray_Return(ap_S.release()), ret_lst);
}



static PyObject*
_linalg_lstsq(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_Am = NULL;
    PyArrayObject *ap_b = NULL;
    double rcond;
    const char *driver_name = NULL;
    int overwrite_a = 0;
    int overwrite_b = 0;
    SliceStatusVec vec_status;

    if (!PyArg_ParseTuple(args, "O!O!ds|pp", &PyArray_Type, (PyObject **)&ap_Am,  &PyArray_Type, (PyObject **)&ap_b, &rcond, &driver_name, &overwrite_a, &overwrite_b)) {
        return NULL;
    }
    if (!_check_dtype_and_flags(ap_Am, "lstsq")) { return NULL; }
    if (!_check_dtype_and_flags(ap_b, "lstsq")) { return NULL; }
    if (PyArray_TYPE(ap_b) != PyArray_TYPE(ap_Am)) {
        PyErr_SetString(PyExc_ValueError, "`a` and `b` must have the same dtype.");
        return NULL;
    }
    if (rcond < 0) {
        PyErr_SetString(PyExc_ValueError, "Expected rcond >= 0.");
        return NULL;
    }

    LstsqDriver driver;
    if (strcmp(driver_name, "gelsd") == 0) { driver = LstsqDriver::gelsd; }
    else if (strcmp(driver_name, "gelsy") == 0) { driver = LstsqDriver::gelsy; }
    else if (strcmp(driver_name, "gelss") == 0) { driver = LstsqDriver::gelss; }
    else { return PyErr_Format(PyExc_ValueError, "scipy.linalg.lstsq: unknown lapack_driver '%s'.", driver_name); }

    int typenum = PyArray_TYPE(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    npy_intp *shape_b = PyArray_SHAPE(ap_b);

    // The Python layer made b at least 2-D and broadcast the batch dimensions of a and b, so
    // a is (..., m, n) and b is (..., m, nrhs) with the same leading dimensions.
    if (PyArray_NDIM(ap_b) != ndim || !std::equal(shape, shape + ndim - 2, shape_b)) {
        PyErr_SetString(PyExc_ValueError, "`a` and `b` shape mismatch.");
        return NULL;
    }

    npy_intp m = shape[ndim - 2], n = shape[ndim - 1], nrhs = shape_b[ndim - 1];
    npy_intp min_mn = m < n ? m : n, max_mn = m < n ? n : m;

    // x is (..., max(1, m, n), nrhs): the solvers write x over the first n rows of b
    npy_intp shape_1[NPY_MAXDIMS];
    std::copy(shape, shape + ndim, shape_1);
    shape_1[ndim - 2] = max_mn > 1 ? max_mn : 1;
    shape_1[ndim - 1] = nrhs;
    OwnedArray ap_x;
    if (overwrite_b) {
        // b is 2-D, F-contiguous with m >= n: the solution is computed in place
        Py_INCREF(ap_b);
        ap_x.reset(ap_b);
    } else {
        ap_x.reset((PyArrayObject *)PyArray_SimpleNew(ndim, shape_1, typenum));
    }

    // S is real and (..., min(m, n)); gelsy computes no singular values
    OwnedArray ap_S;
    if (driver != LstsqDriver::gelsy) {
        int typenum_S = (typenum == NPY_COMPLEX64) ? NPY_FLOAT32 : (typenum == NPY_COMPLEX128) ? NPY_FLOAT64 : typenum;
        shape_1[ndim - 2] = min_mn;
        ap_S.reset((PyArrayObject *)PyArray_SimpleNew(ndim - 1, shape_1, typenum_S));
        if (!ap_S) { return NULL; }
    }

    // rank has the batch shape
    OwnedArray ap_rank((PyArrayObject *)PyArray_SimpleNew(ndim - 2, shape_1, NPY_INT64));
    if (!ap_x || !ap_rank) { return NULL; }

    BatchStatus status;
    Py_BEGIN_ALLOW_THREADS
    status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _lstsq<T>(ap_Am, ap_b, ap_S.get(), ap_x.get(), ap_rank.get(), rcond, driver, overwrite_a, overwrite_b, vec_status);
    });
    Py_END_ALLOW_THREADS
    if (status != BatchStatus::ok) { return set_batch_error(status, "lstsq"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    if (!ret_lst) { return NULL; }
    PyObject *ret_S = ap_S ? PyArray_Return(ap_S.release()) : Py_NewRef(Py_None);
    return Py_BuildValue("NNNN", PyArray_Return(ap_x.release()), PyArray_Return(ap_rank.release()), ret_S, ret_lst);
}



static PyObject*
_linalg_eig(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_Am = NULL;
    PyArrayObject *ap_Bm = NULL;
    int compute_vl = 0;
    int compute_vr = 1;
    int overwrite_a = 0;
    int overwrite_b = 0;
    SliceStatusVec vec_status;

    if (!PyArg_ParseTuple(args, "O!pppp|O!",
            &PyArray_Type, (PyObject **)&ap_Am,
            &compute_vl, &compute_vr,
            &overwrite_a, &overwrite_b,
            &PyArray_Type, (PyObject **)&ap_Bm)
    ) {
        return NULL;
    }
    if (!_check_dtype_and_flags(ap_Am, "eig")) { return NULL; }
    if (ap_Bm && !_check_dtype_and_flags(ap_Bm, "eig")) { return NULL; }

    int typenum = PyArray_TYPE(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    npy_intp n = shape[ndim - 1];
    if (shape[ndim - 2] != n) {
        PyErr_SetString(PyExc_ValueError, "Expected a square matrix");
        return NULL;
    }
    // the Python layer broadcast a and b to a common shape and dtype
    if (ap_Bm && (PyArray_TYPE(ap_Bm) != typenum || PyArray_NDIM(ap_Bm) != ndim
                  || !std::equal(shape, shape + ndim, PyArray_SHAPE(ap_Bm)))) {
        PyErr_SetString(PyExc_ValueError, "`a` and `b` must have the same shape and dtype.");
        return NULL;
    }

    // Eigenvalues and eigenvectors are complex even for real input; the Python wrapper casts
    // them back to real when every eigenvalue is real. w and beta are (..., n), vl and vr (..., n, n).
    int w_typenum = (typenum == NPY_FLOAT32) ? NPY_COMPLEX64 : (typenum == NPY_FLOAT64) ? NPY_COMPLEX128 : typenum;
    OwnedArray ap_w((PyArrayObject *)PyArray_SimpleNew(ndim - 1, shape, w_typenum));
    OwnedArray ap_beta, ap_vl, ap_vr;
    if (ap_Bm) { ap_beta.reset((PyArrayObject *)PyArray_SimpleNew(ndim - 1, shape, typenum)); }
    if (compute_vl) { ap_vl.reset((PyArrayObject *)PyArray_SimpleNew(ndim, shape, w_typenum)); }
    if (compute_vr) { ap_vr.reset((PyArrayObject *)PyArray_SimpleNew(ndim, shape, w_typenum)); }
    if (!ap_w || (ap_Bm && !ap_beta) || (compute_vl && !ap_vl) || (compute_vr && !ap_vr)) { return NULL; }

    BatchStatus status;
    Py_BEGIN_ALLOW_THREADS
    status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _eig<T>(ap_Am, ap_Bm, ap_w.get(), ap_beta.get(), ap_vl.get(), ap_vr.get(), overwrite_a, overwrite_b, vec_status);
    });
    Py_END_ALLOW_THREADS
    if (status != BatchStatus::ok) { return set_batch_error(status, "eig"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    if (!ret_lst) { return NULL; }
    PyObject *ret_beta = ap_beta ? PyArray_Return(ap_beta.release()) : Py_NewRef(Py_None);
    PyObject *ret_vl = ap_vl ? PyArray_Return(ap_vl.release()) : Py_NewRef(Py_None);
    PyObject *ret_vr = ap_vr ? PyArray_Return(ap_vr.release()) : Py_NewRef(Py_None);
    return Py_BuildValue("NNNNN", PyArray_Return(ap_w.release()), ret_beta, ret_vl, ret_vr, ret_lst);
}



static PyObject*
_linalg_eigh(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_Am = NULL;
    PyArrayObject *ap_Bm = NULL;
    int overwrite_a = 0;
    int overwrite_b = 0;
    int eigvals_only = 0;
    int lower = 0;
    int vals_range = 0;
    double vl = -std::numeric_limits<double>::infinity();
    double vu = std::numeric_limits<double>::infinity();
    int il = -1;
    int iu = -1;
    Py_ssize_t driver_index = 0;
    int itype = 1;
    SliceStatusVec vec_status;

    if (!PyArg_ParseTuple(args, "O!ppipddiin|O!pi",
        &PyArray_Type, (PyObject **)&ap_Am, &overwrite_a,
        &eigvals_only, &vals_range, &lower,
        &vl, &vu, &il, &iu, &driver_index,
        &PyArray_Type, (PyObject **)&ap_Bm, &overwrite_b, &itype)
    ) {
        return NULL;
    }

    Eigh_driver driver = (Eigh_driver)driver_index;
    bool generalized;
    switch (driver) {
        case Eigh_driver::EV: case Eigh_driver::EVD: case Eigh_driver::EVR: case Eigh_driver::EVX:
            generalized = false; break;
        case Eigh_driver::GV: case Eigh_driver::GVD: case Eigh_driver::GVX:
            generalized = true; break;
        default:
            return PyErr_Format(PyExc_ValueError, "scipy.linalg.eigh: unknown driver index %zd.", driver_index);
    }
    if (generalized != (ap_Bm != NULL)) {
        PyErr_SetString(PyExc_ValueError, "scipy.linalg.eigh: the gv, gvd and gvx drivers need `b` and the others take none.");
        return NULL;
    }

    /*
     * - `uplo`: the triangle of the matrix to use
     * - `jobz`: whether to compute eigenvectors
     * - `range`, for the evr, evx and gvx drivers:
     *      - 'A': all n eigenvalues; il, iu, vl and vu are sentinels
     *      - 'V': the eigenvalues w_i with vl < w_i <= vu; their number is unknown beforehand, so
     *             room for all n is allocated and the number found, m, is returned to Python
     *      - 'I': in ascending order, the 0-based indices il to iu inclusive, m = iu - il + 1
     */
    const char uplo = lower ? 'L' : 'U';
    const char jobz = eigvals_only ? 'N' : 'V';
    const char range = (vals_range == 0) ? 'A' : (vals_range == 1) ? 'V' : 'I';

    if (!_check_dtype_and_flags(ap_Am, "eigh")) { return NULL; }
    if (ap_Bm && !_check_dtype_and_flags(ap_Bm, "eigh")) { return NULL; }

    int typenum = PyArray_TYPE(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    npy_intp n = shape[ndim - 1];
    if (shape[ndim - 2] != n) {
        PyErr_SetString(PyExc_ValueError, "Expected a square matrix");
        return NULL;
    }
    // the Python layer broadcast a and b to a common shape and dtype
    if (ap_Bm && (PyArray_TYPE(ap_Bm) != typenum || PyArray_NDIM(ap_Bm) != ndim
                  || !std::equal(shape, shape + ndim, PyArray_SHAPE(ap_Bm)))) {
        PyErr_SetString(PyExc_ValueError, "`a` and `b` must have the same shape and dtype.");
        return NULL;
    }

    // w is real and (..., m); Z is (..., n, m), or the F-contiguous input itself when overwriting
    int w_typenum = (typenum == NPY_COMPLEX64) ? NPY_FLOAT32 : (typenum == NPY_COMPLEX128) ? NPY_FLOAT64 : typenum;
    npy_intp m = (range == 'I') ? iu + 1 - il : n;
    npy_intp shape_1[NPY_MAXDIMS];
    std::copy(shape, shape + ndim, shape_1);
    shape_1[ndim - 2] = m;
    OwnedArray ap_w((PyArrayObject *)PyArray_SimpleNew(ndim - 1, shape_1, w_typenum));
    OwnedArray ap_Z;
    if (jobz != 'N') {
        if (overwrite_a) {
            Py_INCREF(ap_Am);
            ap_Z.reset(ap_Am);
        } else {
            shape_1[ndim - 2] = n;
            shape_1[ndim - 1] = m;
            ap_Z.reset((PyArrayObject *)PyArray_SimpleNew(ndim, shape_1, typenum));
        }
    }
    if (!ap_w || (jobz != 'N' && !ap_Z)) { return NULL; }

    npy_intp m_found = m;
    BatchStatus status;
    Py_BEGIN_ALLOW_THREADS
    status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _eigh<T>(ap_Am, ap_Bm, ap_w.get(), ap_Z.get(), &m_found, overwrite_a, overwrite_b, itype,
                        jobz, range, uplo, vl, vu, il, iu, driver, vec_status);
    });
    Py_END_ALLOW_THREADS
    if (status != BatchStatus::ok) { return set_batch_error(status, "eigh"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    if (!ret_lst) { return NULL; }
    PyObject *ret_Z = ap_Z ? PyArray_Return(ap_Z.release()) : Py_NewRef(Py_None);
    return Py_BuildValue("NNnN", PyArray_Return(ap_w.release()), ret_Z, (Py_ssize_t)m_found, ret_lst);
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
    OwnedArray ap_Cm;
    if (overwrite_a) {
        Py_INCREF(ap_Am);
        ap_Cm.reset(ap_Am);
    }
    else {
        ap_Cm.reset((PyArrayObject *)PyArray_ZEROS(ndim, shape, typenum, 0));
        if (!ap_Cm) { return NULL; }
    }

    BatchStatus status;
    Py_BEGIN_ALLOW_THREADS
    status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _cholesky<T>(ap_Am, ap_Cm.get(), lower, overwrite_a, clean, vec_status);
    });
    Py_END_ALLOW_THREADS
    if (status != BatchStatus::ok) { return set_batch_error(status, "cholesky"); }

    PyObject *ret_lst = convert_vec_status(vec_status);
    if (!ret_lst) { return NULL; }
    return Py_BuildValue("NN", PyArray_Return(ap_Cm.release()), ret_lst);
}


/*
 * Helper: convert a vector of slice error statuses to list of dicts
 */
PyObject*
convert_vec_status(SliceStatusVec& vec_status) {
    PyObject *ret_lst = PyList_New(0);
    if (!ret_lst) { return NULL; }

    for (const SliceStatus &status : vec_status) {
        PyObject *ret_dct = Py_BuildValue(
            "{s:n,s:n,s:i,s:i,s:d,s:n}",
            "num", status.slice_num,
            "structure", status.structure,
            "is_singular", status.is_singular,
            "is_ill_conditioned", status.is_ill_conditioned,
            "rcond", status.rcond,
            "lapack_info", status.lapack_info
        );
        if (!ret_dct || PyList_Append(ret_lst, ret_dct) < 0) {
            Py_XDECREF(ret_dct);
            Py_DECREF(ret_lst);
            return NULL;
        }
        Py_DECREF(ret_dct);
    }
    return ret_lst;
}



static PyObject*
_linalg_lu(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_a = NULL;
    int permute_l = 0;
    int overwrite_a = 0;

    if (!PyArg_ParseTuple(args, "O!pp", &PyArray_Type, (PyObject **)&ap_a, &permute_l, &overwrite_a)) {
        return NULL;
    }
    if (!_check_dtype_and_flags(ap_a, "lu")) { return NULL; }

    int typenum = PyArray_TYPE(ap_a);
    int ndim = PyArray_NDIM(ap_a);
    npy_intp *shape = PyArray_SHAPE(ap_a);
    npy_intp m = shape[ndim - 2], n = shape[ndim - 1], minmn = m < n ? m : n;

    // perm (..., m); L (..., m, minmn); U (..., minmn, n)
    OwnedArray ap_perm((PyArrayObject *)PyArray_SimpleNew(ndim - 1, shape, (sizeof(CBLAS_INT) == 8) ? NPY_INT64 : NPY_INT32));
    npy_intp l_shape[NPY_MAXDIMS], u_shape[NPY_MAXDIMS];
    std::copy(shape, shape + ndim, l_shape);
    std::copy(shape, shape + ndim, u_shape);
    l_shape[ndim - 1] = minmn;
    u_shape[ndim - 2] = minmn;
    OwnedArray ap_l((PyArrayObject *)PyArray_ZEROS(ndim, l_shape, typenum, 0));
    OwnedArray ap_u((PyArrayObject *)PyArray_ZEROS(ndim, u_shape, typenum, 0));
    if (!ap_perm || !ap_l || !ap_u) { return NULL; }

    // A 2-D F-contiguous writeable input is the getrf workspace itself: no copy is made.
    bool in_place = overwrite_a && (ndim == 2) && PyArray_ISFARRAY_RO(ap_a) && PyArray_ISWRITEABLE(ap_a);

    BatchStatus status;
    Py_BEGIN_ALLOW_THREADS
    status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _lu<T>(ap_a, static_cast<T *>(PyArray_DATA(ap_l.get())), static_cast<T *>(PyArray_DATA(ap_u.get())),
                      static_cast<CBLAS_INT *>(PyArray_DATA(ap_perm.get())), (bool)permute_l, in_place);
    });
    Py_END_ALLOW_THREADS
    if (status != BatchStatus::ok) { return set_batch_error(status, "lu"); }

    if (permute_l) {
        // P is already applied to L. Return an empty int32 array rather than None so the
        // return type is always three arrays.
        npy_intp zero = 0;
        ap_perm.reset((PyArrayObject *)PyArray_EMPTY(1, &zero, NPY_INT32, 0));
        if (!ap_perm) { return NULL; }
    }
    return Py_BuildValue("NNN", PyArray_Return(ap_perm.release()),
                                PyArray_Return(ap_l.release()),
                                PyArray_Return(ap_u.release()));
}


// ========================================================================
// Determinant
// ========================================================================

static PyObject*
_linalg_det(PyObject* Py_UNUSED(dummy), PyObject* args) {
    PyArrayObject *ap_a = NULL;
    int overwrite_a = 0;

    if (!PyArg_ParseTuple(args, "O!p", &PyArray_Type, (PyObject **)&ap_a, &overwrite_a)) {
        return NULL;
    }
    if (!_check_dtype_and_flags(ap_a, "det")) { return NULL; }

    int typenum = PyArray_TYPE(ap_a);
    int ndim = PyArray_NDIM(ap_a);
    npy_intp *shape = PyArray_SHAPE(ap_a);
    if (shape[ndim - 2] != shape[ndim - 1]) {
        PyErr_SetString(PyExc_ValueError, "scipy.linalg.det: Last two dimensions must be square.");
        return NULL;
    }

    // one determinant per slice, so the batch dimensions only; double precision for single precision input
    int det_typenum = (typenum == NPY_FLOAT32) ? NPY_FLOAT64 : (typenum == NPY_COMPLEX64) ? NPY_COMPLEX128 : typenum;
    OwnedArray ap_det((PyArrayObject *)PyArray_SimpleNew(ndim - 2, shape, det_typenum));
    if (!ap_det) { return NULL; }

    // A 2-D contiguous writeable input is the getrf workspace itself; det(A^T) = det(A), so
    // either order works.
    bool in_place = overwrite_a && (ndim == 2) && PyArray_ISWRITEABLE(ap_a)
                    && (PyArray_ISFARRAY_RO(ap_a) || PyArray_ISCARRAY_RO(ap_a));

    BatchStatus status;
    Py_BEGIN_ALLOW_THREADS
    status = dispatch_dtype(typenum, [&](auto tag) {
        using T = decltype(tag);
        return _det<T>(ap_a, static_cast<det_type_t<T> *>(PyArray_DATA(ap_det.get())), in_place);
    });
    Py_END_ALLOW_THREADS
    if (status != BatchStatus::ok) { return set_batch_error(status, "det"); }

    return (PyObject *)ap_det.release();
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
    if (!PyArray_ISBEHAVED_RO(ap_a)) {
        PyErr_SetString(PyExc_TypeError, "Input array must be aligned and in native byte order.");
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

    OwnedArray ap_lb((PyArrayObject *)PyArray_SimpleNew(batch_ndim, shape, NPY_INT64));
    OwnedArray ap_ub((PyArrayObject *)PyArray_SimpleNew(batch_ndim, shape, NPY_INT64));
    if (!ap_lb || !ap_ub) { return NULL; }

    int64_t *lb_data = static_cast<int64_t *>(PyArray_DATA(ap_lb.get()));
    int64_t *ub_data = static_cast<int64_t *>(PyArray_DATA(ap_ub.get()));
    char *a_data = PyArray_BYTES(ap_a);

    auto each_slice = [&](auto tag) {
        using T = decltype(tag);
        for (npy_intp idx = 0; idx < num_slices; idx++) {
            const T *slice_ptr = compute_slice_ptr<const T>(idx, a_data, ndim, shape, strides);
            bandwidth_slice(slice_ptr, n, m, s0, s1, &lb_data[idx], &ub_data[idx]);
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
            PyErr_SetString(PyExc_TypeError, "Unsupported dtype.");
            return NULL;
    }

    // 2D: return tuple of int scalars; N-d: return tuple of int64 arrays
    if (ndim == 2) {
        // TODO: use PyLong_FromInt64 from py314+
        PyObject *lb_obj = PyLong_FromSsize_t((Py_ssize_t)lb_data[0]);
        PyObject *ub_obj = PyLong_FromSsize_t((Py_ssize_t)ub_data[0]);
        return Py_BuildValue("NN", lb_obj, ub_obj);
    }

    return Py_BuildValue("NN", (PyObject *)ap_lb.release(), (PyObject *)ap_ub.release());
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
