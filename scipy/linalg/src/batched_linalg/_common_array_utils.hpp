#pragma once

namespace sp_linalg {


using wrapper::f32;
using wrapper::f64;
using wrapper::c64;
using wrapper::c128;
using wrapper::real_of_t;
using wrapper::complex_of_t;
using wrapper::is_complex_v;


// `*con` routines need iwork for real flavors and rwork for complex ones.
template <class T> using irwork_t = typename std::conditional<is_complex_v<T>, real_of_t<T>, CBLAS_INT>::type;


/**
 * @brief Outcome of a whole batched call, as opposed to the per-slice SliceStatus.
 *
 * A failing slice does not fail the call: drivers record it in the SliceStatusVec and return ok.
 */
enum class BatchStatus {
    ok,
    no_memory,            ///< A workspace allocation failed.
    too_large,            ///< A dimension or workspace size passed to LAPACK does not fit CBLAS_INT.
    lapack_error,         ///< A LAPACK routine rejected an argument (info < 0).
    bad_dtype,            ///< The array dtype is not float32, float64, complex64 or complex128.
};


/**
 * @brief The deleter of Workspace: returns the memory with PyMem_RawFree.
 */
struct WorkspaceDeleter {
    void operator()(void *p) const noexcept { PyMem_RawFree(p); }
};


/**
 * @brief Owning pointer to an uninitialized workspace array obtained from PyMem_RawMalloc.
 */
template <class T>
using Workspace = std::unique_ptr<T[], WorkspaceDeleter>;


/**
 * @brief Allocate an uninitialized workspace array.
 *
 * @param count  Number of elements of type T; 0 gives a valid pointer to no elements.
 * @return An owning pointer to the array, never empty.
 * @throws std::bad_alloc if `count * sizeof(T)` exceeds SIZE_MAX or the allocation fails.
 */
template <class T>
Workspace<T> alloc_workspace(std::uint64_t count)
{
    if (count > SIZE_MAX / sizeof(T)) { throw std::bad_alloc(); }
    T *p = static_cast<T *>(PyMem_RawMalloc(static_cast<std::size_t>(count) * sizeof(T)));
    if (!p) { throw std::bad_alloc(); }
    return Workspace<T>(p);
}


/**
 * @brief The deleter of OwnedArray: releases the reference with Py_XDECREF.
 */
struct OwnedArrayDeleter {
    void operator()(PyArrayObject *p) const noexcept { Py_XDECREF(p); }
};


/**
 * @brief Owning reference to a NumPy array; `release()` hands the reference to the caller.
 */
using OwnedArray = std::unique_ptr<PyArrayObject, OwnedArrayDeleter>;


/**
 * @brief Convert a nonnegative size to a narrower integer type, detecting loss.
 *
 * @param v    The size; negative values are rejected.
 * @param out  Receives `v` converted to To; written only on success.
 * @return true if `v` is nonnegative and representable in To.
 */
template <class To, class From>
bool checked_narrow(From v, To *out) noexcept
{
    static_assert(std::is_integral_v<To> && std::is_integral_v<From>);
    if constexpr (std::is_signed_v<From>) {
        if (v < 0) { return false; }
    }
    using UFrom = std::make_unsigned_t<From>;
    using UTo = std::make_unsigned_t<To>;
    if (static_cast<UFrom>(v) > static_cast<UTo>(std::numeric_limits<To>::max())) { return false; }
    *out = static_cast<To>(v);
    return true;
}


/**
 * @brief Call a generic callable with a value of the C++ type matching a NumPy dtype.
 *
 * Callers run this with the thread state detached (Py_BEGIN_ALLOW_THREADS), so @p f may read
 * PyArrayObject fields (PyArray_DATA, PyArray_SHAPE, ...) but must not create objects, change
 * reference counts or set exceptions.
 *
 * @param typenum  NumPy type number; NPY_FLOAT32, NPY_FLOAT64, NPY_COMPLEX64 and NPY_COMPLEX128 are
 *                 dispatched.
 * @param f        Callable invoked as `f(f32{})`, `f(f64{})`, `f(c64{})` or `f(c128{})`, returning a
 *                 BatchStatus; it may throw std::bad_alloc, e.g. from alloc_workspace.
 * @return The result of @p f; BatchStatus::bad_dtype for any other @p typenum;
 *         BatchStatus::no_memory if @p f throws std::bad_alloc.
 */
template <class F>
BatchStatus dispatch_dtype(int typenum, F &&f) noexcept
{
    try {
        switch (typenum) {
            case NPY_FLOAT32:    return f(f32{});
            case NPY_FLOAT64:    return f(f64{});
            case NPY_COMPLEX64:  return f(c64{});
            case NPY_COMPLEX128: return f(c128{});
            default:             return BatchStatus::bad_dtype;
        }
    } catch (const std::bad_alloc &) {
        return BatchStatus::no_memory;
    }
}


// Structure tags; python side maps assume_a strings to these values
enum St : Py_ssize_t
{
    NONE = -1,
    GENERAL = 0,
    DIAGONAL = 11,
    TRIDIAGONAL = 31,
    BANDED = 41,
    UPPER_TRIANGULAR = 21,
    LOWER_TRIANGULAR = 22,
    POS_DEF = 101,
    SYM = 201,
    HER = 211
};

// QR mode tags; python side maps mode strings to these values
enum QR_mode : Py_ssize_t
{
    FULL = 1,
    R = 11,
    RAW_MODE = 21,
    ECONOMIC = 31
};

// Eigh driver tags; pythn side maps driver strings to these values
enum Eigh_driver : Py_ssize_t
{
    EV = 1,
    EVD = 2,
    EVR = 3,
    EVX = 4,
    GV = 10,
    GVD = 11,
    GVX = 12
};


/*
 * Rich return object
 */
struct SliceStatus {
    Py_ssize_t slice_num;
    Py_ssize_t structure;
    int is_singular;
    int is_ill_conditioned;
    double rcond;
    Py_ssize_t lapack_info;
};


inline void init_status(SliceStatus& slice_status, npy_intp idx, St slice_structure) {
    slice_status.slice_num = idx;
    slice_status.structure = (Py_ssize_t)slice_structure;
    slice_status.is_singular = 0;
    slice_status.is_ill_conditioned = 0;
    slice_status.rcond = 0;
    slice_status.lapack_info = 0;
}


typedef std::vector<SliceStatus> SliceStatusVec;


/*
 * When looping over slices, each `solve_slice_XYZ` return the `status` struct,
 * which records one of three possible outcomes:
 *
 *  - the slice is singular or LAPACK returned non-zero `info`
 *  - the slice is detected to be ill-conditioned
 *  - all is well
 *
 * In the first two cases, we record the `status` in `vec_status`.
 * For non-recoverable errors (singularity of non-zero info), we terminate the loop over the slices.
 */
inline int
_detect_problems(const SliceStatus& slice_status, SliceStatusVec& vec_status) {
    if ((slice_status.lapack_info < 0) || (slice_status.is_singular)) {
        vec_status.push_back(slice_status);
        return 1;
    }
    else if (slice_status.is_ill_conditioned) {
        vec_status.push_back(slice_status);
    }
    return 0;
}


/*
 * lwork defensive handler :
 *  cf https://github.com/scipy/scipy/blob/v1.15.2/scipy/linalg/lapack.py#L1004
 *
 *  Round floating-point lwork returned by lapack to integer.
 *
 *  Several LAPACK routines compute optimal values for LWORK, which
 *  they return in a floating-point variable. However, for large
 *  values of LWORK, single-precision floating point is not sufficient
 *  to hold the exact value --- some LAPACK versions (<= 3.5.0 at
 *  least) truncate the returned integer to single precision and in
 *  some cases this can be smaller than the required value.
 *
 *  The fudge_factor comes from
 *  https://github.com/scipy/scipy/blob/v1.15.2/scipy/linalg/_basic.py#L1154
 *
 *  A fudge factor of 1% was added in commit https://github.com/scipy/scipy/commit/dfb543c147c
 *  to avoid a "curious segfault with 500x500 matrices and OpenBLAS".
 */
template<typename T>
CBLAS_INT _calc_lwork(T _lwrk, f64 fudge_factor=1.0) {
    using real_type = real_of_t<T>;

    real_type value = std::real(_lwrk) * fudge_factor;
    if (std::is_same<real_type, f32>::value) {
        // Single-precision routine -- take next fp value to work
        // around possible truncation in LAPACK code
        value = std::nextafter(value, std::numeric_limits<real_type>::infinity());
    }

    CBLAS_INT lwork;
    // 2**digits is one above the largest CBLAS_INT and exact in either precision, so the values
    // strictly below it, and only those, convert to CBLAS_INT without overflow
    const real_type limit = std::ldexp(real_type(1), std::numeric_limits<CBLAS_INT>::digits);
    if ((value < 0) || !(value < limit)) {
        // Too large lwork required - Computation cannot be performed with standard LAPACK
        lwork = -1;
    }
    else {
        lwork = value > 0 ? (CBLAS_INT)value : 1;
    }
    return lwork;
}


/**
 * @brief Pointer to the first element of the idx-th 2-D slice of an array of ndim >= 2.
 *
 * @param idx      Slice index in C order over the leading ndim - 2 dimensions.
 * @param data     First element of the array; may point to const.
 * @param ndim     Number of dimensions.
 * @param shape    Shape of the array.
 * @param strides  Byte strides of the array.
 * @return @p data advanced by the byte offset of the slice.
 */
template<typename T>
T* compute_slice_ptr(npy_intp idx, T *data, npy_intp ndim, const npy_intp *shape, const npy_intp *strides) {
    npy_intp offset = 0;
    npy_intp temp_idx = idx;
    for (int i = ndim - 3; i >= 0; i--) {
        offset += (temp_idx % shape[i]) * strides[i];
        temp_idx /= shape[i];
    }
    using byte_t = std::conditional_t<std::is_const_v<T>, const char, char>;
    return reinterpret_cast<T *>(reinterpret_cast<byte_t *>(data) + offset);
}


/**
 * @brief Copy a strided m-by-n matrix into a column-major buffer.
 *
 * @param dst  Column-major destination; element (i, j) is written to `dst[i + j*ld]`.
 * @param src  First element of the source matrix.
 * @param m    Number of rows.
 * @param n    Number of columns.
 * @param s0   Byte stride between consecutive rows of @p src; may be negative.
 * @param s1   Byte stride between consecutive columns of @p src; may be negative.
 * @param ld   Leading dimension of @p dst, at least @p m.
 */
template<typename T>
void copy_strided_to_f(T *dst, const T *src, npy_intp m, npy_intp n, npy_intp s0, npy_intp s1, npy_intp ld) noexcept
{
    const char *base = reinterpret_cast<const char *>(src);
    // Walk the source along its smaller stride. The magnitudes are unsigned: NumPy accepts any
    // stride for a dimension of length 1, INT64_MIN included, and std::abs of it is undefined.
    std::uint64_t abs_s0 = s0 < 0 ? 0 - static_cast<std::uint64_t>(s0) : static_cast<std::uint64_t>(s0);
    std::uint64_t abs_s1 = s1 < 0 ? 0 - static_cast<std::uint64_t>(s1) : static_cast<std::uint64_t>(s1);
    if (abs_s1 <= abs_s0) {
        for (npy_intp i = 0; i < m; i++) {
            const char *src_row = base + i * s0;
            for (npy_intp j = 0; j < n; j++) {
                dst[i + j * ld] = *reinterpret_cast<const T *>(src_row + j * s1);
            }
        }
    }
    else {
        for (npy_intp j = 0; j < n; j++) {
            const char *src_col = base + j * s1;
            T *dst_col = dst + j * ld;
            for (npy_intp i = 0; i < m; i++) {
                dst_col[i] = *reinterpret_cast<const T *>(src_col + i * s0);
            }
        }
    }
}


/*
 * Copy n-by-m F-ordered `src` to C-ordered `dst`.
 *
 * `src` is ldb-by-m, F-ordered
 * `dst` is n-by-m, C-ordered
 *
 * The default is to have src and dst of the same size (ldb=-1 means ldb=n)
 */
template<typename T>
void copy_slice_F_to_C(T* dst, const T* src, const npy_intp n, const npy_intp m, npy_intp ldb=-1) {

    if (ldb == -1) {ldb = n;}

    for (npy_intp i = 0; i < n; i++) {
        for (npy_intp j = 0; j < m; j++) {
            dst[i*m + j] = src[i + j*ldb];
        }
    }
}


/*
 * Copy one triangle of a strided m-by-n `src` to C-ordered `dst`.
 * Note that by specifying m and/or n smaller than the input array `src`
 * it is possible to only copy a subset of the triangle.
 *
 * Only elements in the `uplo` triangle are copied;
 * `dst` is assumed zero-initialized for the other triangle.
 *
 * The function does not use the symmetry of `src` — it reads exactly the
 * triangle specified by `uplo`. The caller is responsible for ensuring
 * that the correct triangle is populated in `src`.
 *
 * `s0` and `s1` are the byte strides between rows and between columns of `src`.
 *
 * Examples (matrix dimension m x n, element size sz):
 *   C-contiguous src (s0=n*sz, s1=sz): src[i*n + j] -> dst[i*n + j]
 *      both sequential in the inner loop - effectively a partial memcpy.
 *   F-contiguous src (s0=sz, s1=m*sz): src[i + j*m] -> dst[i*n + j],
 *      column-sequential reads, row-sequential writes.
 */
template<typename T>
void copy_triangle_to_C(T *dst, const T *src, const npy_intp m, const npy_intp n, const npy_intp s0, const npy_intp s1, const char uplo) {
    const char *base = reinterpret_cast<const char *>(src);
    auto at = [&](npy_intp i, npy_intp j) { return *reinterpret_cast<const T *>(base + i * s0 + j * s1); };
    if (uplo == 'L') {
        for (npy_intp i = 0; i < m; i++) {
            npy_intp stop = std::min(i + 1, n);
            for (npy_intp j = 0; j < stop; j++) {
                dst[i * n + j] = at(i, j);
            }
        }
    } else {
        for (npy_intp i = 0; i < m; i++) {
            for (npy_intp j = i; j < n; j++) {
                dst[i * n + j] = at(i, j);
            }
        }
    }
}


/*
 * 1-norm of a matrix
 */

template<typename T>
real_of_t<T>
norm1_(T* A, const npy_intp n)
{
    using real_type = real_of_t<T>;

    real_type norm = 0.0;
    for (CBLAS_INT i = 0; i < n; i++) {
        real_type tmp = 0.0;
        for (CBLAS_INT j = 0; j < n; j++) {
            tmp += std::abs(A[i * n + j]);
        }

        if (tmp > norm) { norm = tmp; }
    }

    return norm;
}


template<typename T>
real_of_t<T>
norm1_sym_herm_upper(T* A, T* work, const npy_intp n)
{
    using real_type = real_of_t<T>;

    Py_ssize_t i, j;
    real_type temp = 0.0;
    real_type *rwork = (real_type *)work;

    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { rwork[i] = std::abs(A[i]);
     }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) {
        // only loop over the upper triangle
        rwork[i] += std::abs(A[i*n + i]);
        for (j = i+1; j < n; j++) {
            temp = std::abs(A[i*n + j]);
            rwork[j] += temp;
            rwork[i] += temp;
        }
    }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}


template<typename T>
real_of_t<T>
norm1_sym_herm_lower(T* A, T* work, const npy_intp n)
{
    using real_type = real_of_t<T>;

    Py_ssize_t i, j;
    real_type temp = 0.0;
    real_type *rwork = (real_type *)work;

    for (i = 0; i < n; i++) { rwork[i] = 0.0; }

    for (i=0; i < n; i++) {
        rwork[i] += std::abs(A[i*n + i]);
        for (j=0; j < i; j++) {
            temp = std::abs(A[i*n + j]);
            rwork[j] += temp;
            rwork[i] += temp;
        }
    }

    temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}


template<typename T>
real_of_t<T>
norm1_sym_herm(char uplo, T *A, T *work, const npy_intp n) {
    // NB: transpose for the F order; uplo is 'U' or 'L' by construction
    return (uplo == 'U') ? norm1_sym_herm_lower(A, work, n) : norm1_sym_herm_upper(A, work, n);
}


template<typename T>
real_of_t<T>
norm1_tridiag(T* dl, T *d, T *du, T *work, const npy_intp n) {
    using real_type = real_of_t<T>;

    real_type *rwork = (real_type *)work;

    npy_intp i;
    for (i=0; i<n; i++) {
        rwork[i] = std::abs(d[i]);
    }
    for (i=0; i<n-1; i++) {
        rwork[i] += std::abs(dl[i]);
    }
    for (i=1; i<n-1; i++) {
        rwork[i] += std::abs(du[i-1]);
    }

    real_type temp = 0.0;
    for (i = 0; i < n; i++) { if (rwork[i] > temp) { temp = rwork[i]; } }
    return temp;
}

/*
 * Compute the 1 norm of a matrix `A`, but assume it is already in its banded
 * form `ab` as constructed by `to_banded`. It is assumed that the size of `ab`
 * is always such that its number of rows is `2 * kl + ku + 1`.
 */
template <typename T>
real_of_t<T>
norm1_banded(T* ab, const npy_intp kl, const npy_intp ku, T* work, const npy_intp n) {
    using real_type = real_of_t<T>;

    real_type *rwork = (real_type *)work;

    npy_intp i, j;
    npy_intp ldab = 2 * kl + ku + 1;

    for (i = 0; i < n; i++) {
        rwork[i] = std::abs(ab[i * ldab + kl + ku]);
    }

    for (i = 0; i < kl; i++) { // run over lower bands
        for (j = 0; j < n - i - 1; j++) {
            rwork[j] += std::abs(ab[j * ldab + kl + ku + i + 1]);
        }
    }


    for (i = 0; i < ku; i++) { // run over upper bands
        for (j = i + 1; j < n; j++) {
            rwork[j] += std::abs(ab[j * ldab + kl + ku - i - 1]);
        }
    }

    real_type temp = 0.0;
    for (i = 0; i < n; i++) {if (rwork[i] > temp) {temp = rwork[i];} }
    return temp;
}


/***************************
 ***  Structure detection
 ***************************/

/**
 * @brief Lower and upper bandwidth of a strided n-by-m matrix of any numeric type.
 *
 * @param data   First element of the matrix.
 * @param n      Number of rows.
 * @param m      Number of columns.
 * @param s0     Byte stride between consecutive rows; may be negative.
 * @param s1     Byte stride between consecutive columns; may be negative.
 * @param lower  Receives the largest r - c over the nonzero entries (r, c), or 0.
 * @param upper  Receives the largest c - r over the nonzero entries (r, c), or 0.
 */
template<typename T>
void
bandwidth_strided(const T *data, int64_t n, int64_t m, int64_t s0, int64_t s1, int64_t *lower, int64_t *upper) noexcept
{
    const char *base = reinterpret_cast<const char *>(data);
    auto at = [&](int64_t r, int64_t c) { return *reinterpret_cast<const T *>(base + r * s0 + c * s1); };
    const T zero = T(0);
    int64_t lb = 0, ub = 0;

    for (int64_t r = n - 1; r > 0; r--) {
        int64_t ncols = r - lb;
        if (ncols > m) { ncols = m; }
        for (int64_t c = 0; c < ncols; c++) {
            if (at(r, c) != zero) { lb = r - c; break; }
        }
        if (r <= lb) { break; }
    }

    for (int64_t r = 0; r < n; r++) {
        for (int64_t c = m - 1; c > r + ub; c--) {
            if (at(r, c) != zero) { ub = c - r; break; }
        }
        if (r + ub + 1 > m) { break; }
    }

    *lower = lb;
    *upper = ub;
}


/**
 * @brief Lower and upper bandwidth of one n-by-m matrix slice.
 *
 * float32, float64, complex64 and complex128 slices that are contiguous in C or F order use the
 * bandwidth_{s,d,c,z} kernels of _linalg_bandwidth.h; every other slice uses bandwidth_strided.
 *
 * @param data   First element of the slice.
 * @param n      Number of rows.
 * @param m      Number of columns.
 * @param s0     Byte stride between consecutive rows; may be negative.
 * @param s1     Byte stride between consecutive columns; may be negative.
 * @param lower  Receives the largest r - c over the nonzero entries (r, c), or 0.
 * @param upper  Receives the largest c - r over the nonzero entries (r, c), or 0.
 */
template<typename T>
void
bandwidth_slice(const T *data, int64_t n, int64_t m, int64_t s0, int64_t s1, int64_t *lower, int64_t *upper) noexcept
{
    if constexpr (std::is_same_v<T, f32> || std::is_same_v<T, f64> ||
                  std::is_same_v<T, c64> || std::is_same_v<T, c128>) {
        constexpr int64_t size = sizeof(T);
        const bool c_order = (s1 == size) && (s0 == m * size);
        const bool f_order = (s0 == size) && (s1 == n * size);
        if (c_order || f_order) {
            // an F-ordered n-by-m slice is the C-ordered m-by-n transpose: swap the outputs
            const int64_t rows = c_order ? n : m, cols = c_order ? m : n;
            int64_t *lo = c_order ? lower : upper, *up = c_order ? upper : lower;
            if constexpr (std::is_same_v<T, f32>)      { bandwidth_s(data, rows, cols, lo, up); }
            else if constexpr (std::is_same_v<T, f64>) { bandwidth_d(data, rows, cols, lo, up); }
            else if constexpr (std::is_same_v<T, c64>) { bandwidth_c(data, rows, cols, lo, up); }
            else                                       { bandwidth_z(data, rows, cols, lo, up); }
            return;
        }
    }
    bandwidth_strided(data, n, m, s0, s1, lower, upper);
}


template<typename T>
void
detect_bandwidths(T* data, npy_intp ndim, npy_intp outer_size, npy_intp *shape, npy_intp *strides, int64_t *kl, int64_t *ku, int64_t *ldab_max) {
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        T* slice_ptr = compute_slice_ptr(idx, data, ndim, shape, strides);

        bandwidth_slice(slice_ptr, shape[ndim-2], shape[ndim-1], strides[ndim-2], strides[ndim-1], &kl[idx], &ku[idx]);
        if (2 * kl[idx] + ku[idx] + 1 > *ldab_max) {
            *ldab_max = 2 * kl[idx] + ku[idx] + 1;
        }
    }
}


template<typename T>
std::tuple<bool, bool>
is_sym_or_herm(const T *data, npy_intp n) {
    // Return a pair of (is_symmetric, is_hermitian)
    bool all_sym = true, all_herm = true;

    for (npy_intp i=0; i < n; i++) {
        for (npy_intp j=0; j < n; j++) {
            T elem1 = data[i*n + j];
            T elem2 = data[i + j*n];
            all_sym = all_sym && (elem1 == elem2);
            all_herm = all_herm && (elem1 == std::conj(elem2));
            if(!(all_sym || all_herm)) {
                // short-circuit : it's neither symmetric not hermitian
                return std::make_tuple(false, false);
            }
        }
    }
    return std::make_tuple(all_sym, all_herm);
}


/*
 * Common matrices
 */

// fill np.triu(a) from np.tril(a) or np.tril(a) from np.triu(a)
template<typename T>
inline void
fill_other_triangle(char uplo, T *data, npy_intp n) {
    if (uplo == 'U') {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=i+1; j<n; j++){
                if constexpr(is_complex_v<T>) {
                    data[j + i*n] = std::conj(data[i + j*n]);
                } else {
                    data[j + i*n] = data[i + j*n];
                }
            }
        }
    } else {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=0; j<i+1; j++){
                if constexpr(is_complex_v<T>) {
                    data[j + i*n] = std::conj(data[i + j*n]);
                } else {
                    data[j + i*n] = data[i + j*n];
                }
            }
        }
    }
}

// XXX deduplicate (conj or noconj)
template<typename T>
inline void
fill_other_triangle_noconj(char uplo, T *data, npy_intp n) {
    if (uplo == 'U') {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=i+1; j<n; j++){
                data[j + i*n] = data[i + j*n];
            }
        }
    } else {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=0; j<i+1; j++){
                data[j + i*n] = data[i + j*n];
            }
        }
    }
}


/*
 * Helper for converting an NxN matrix to tridiagonal form:
 * extract the diagonals of `data` (a full NxN matrix) into du, d, dl
 *
 * the upper and lower subdiagonals, `dl` and `du`, are have length N-1
 * the main diagonal, `d`, has length N
 *
 */
template<typename T>
inline void
to_tridiag(const T *data, npy_intp N, T *du, T *d, T *dl) {
    for (npy_intp i=0; i<N; i++) {
        d[i] = data[i + i*N];
    }
    for (npy_intp i=0; i<N-1; i++) {
        dl[i] = data[i*N + i + 1];
    }
    for (npy_intp i=0; i<N-1; i++) {
        du[i] = data[(i+1)*N + i];
    }
}


/*
 * Helper function for reshuffling a banded slice into the appropriate
 * structure for ?gbcon and ?gbtrf. `s1` and `s2` contain the strides in
 * the column and row direction (`ndim` - 2 and `ndim` - 1, respectively).
 * The result is stored in `ab` in Fortran order.
 *
 * It is assumed that `ab` provides at least `ldab` x `n` memory elements,
 * where ldab >= 2 * `kl` + `ku` + 1
 *
 * Reference: https://www.netlib.org/lapack/explore-html/df/dd6/group__gbtrf_ga682f53142f0398f83f5461c277d23ba2.html#ga682f53142f0398f83f5461c277d23ba2
 */
template<typename T>
inline void
to_banded(const T *data, npy_intp n, npy_intp kl, npy_intp ku, npy_intp ldab, T *ab, npy_intp s1, npy_intp s2) {
    // s1, s2 are byte strides: address elements through a char pointer
    const char *base = reinterpret_cast<const char *>(data);
    auto at = [&](npy_intp byte_offset) { return *reinterpret_cast<const T *>(base + byte_offset); };
    npy_intp i, j;

    // main diagonal
    for (i = 0; i < n; i++) {
        ab[(i + 1) * ldab - kl - 1] = at(i * (s1 + s2));
    }

    // lower bands
    for (i = 0; i < kl; i++) {
        for (j = 0; j < n - i - 1; j++) {
            ab[(j + 1) * ldab - kl + i] = at((i + 1) * s1 + j * (s1 + s2));
        }
    }

    // upper bands
    for (i = 0; i < ku; i++) {
        for (j = i + 1; j < n; j++) {
            ab[(j + 1) * ldab - kl - i - 2] = at((i + 1) * s2 + (j - i - 1) * (s1 + s2));
        }
    }
}


template<typename T>
inline void
zero_other_triangle(char uplo, T *data, const npy_intp m, npy_intp n = -1, npy_intp lda = -1) {
    if (n == -1) { n = m; }
    if (lda == -1) { lda = m; }

    if (uplo == 'U') {
        for (npy_intp i=0; i<n; i++) {
            for (npy_intp j=i+1; j<m; j++){
                data[j + i*lda] = 0.0;
            }
        }
    } else {
        for (npy_intp i=0; i<n; i++) {
            npy_intp stop = std::min(i, m);
            for (npy_intp j=0; j < stop; j++){
                data[j + i*lda] = 0.0;
            }
        }
    }
}

} // namespace sp_linalg
