/*
 * LAPACK declarations.
 */
#pragma once

namespace sp_linalg {

using namespace lapack;
using wrapper::f32;
using wrapper::f64;
using wrapper::c64;
using wrapper::c128;
using wrapper::real_of_t;
using wrapper::complex_of_t;
using wrapper::is_complex_v;


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


void init_status(SliceStatus& slice_status, npy_intp idx, St slice_structure) {
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
int
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
    if((std::is_same<real_type, f32>::value) ||
       (std::is_same<real_type, c64>::value)
    ) {
        // Single-precision routine -- take next fp value to work
        // around possible truncation in LAPACK code
        value = std::nextafter(value, std::numeric_limits<real_type>::infinity());
    }

    CBLAS_INT lwork;
    if ((value < 0) || !(value <= (real_type)std::numeric_limits<CBLAS_INT>::max())) {
        // Too large lwork required - Computation cannot be performed with standard LAPACK
        lwork = -1;
    }
    else {
        lwork = value > 0 ? (CBLAS_INT)value : 1;
    }
    return lwork;
}


/*
 * Given an array of ndim >= 2, compute a pointer to the start of the 2D "slice" idx
 */
template<typename T>
T* compute_slice_ptr(npy_intp idx, T *Am_data, npy_intp ndim, npy_intp *shape, npy_intp *strides) {
    npy_intp offset = 0;
    npy_intp temp_idx = idx;
    for (int i = ndim - 3; i >= 0; i--) {
        offset += (temp_idx % shape[i]) * strides[i];
        temp_idx /= shape[i];
    }
    T* slice_ptr = (T *)(Am_data + (offset/sizeof(T)));
    return slice_ptr;
}


/*
 * Copy n-by-m slice from slice_ptr to dst.
 */
template<typename T>
void copy_slice(T* dst, const T* slice_ptr, const npy_intp n, const npy_intp m, const npy_intp s2, const npy_intp s1) {

    for (npy_intp i = 0; i < n; i++) {
        for (npy_intp j = 0; j < m; j++) {
            dst[i * m + j] = *(slice_ptr + (i*s2/sizeof(T)) + (j*s1/sizeof(T)));
        }
    }
}


/*
 * Copy n-by-m C-order slice from slice_ptr to dst in F-order.
 *
 * `src` is n-by-m, strided
 * `dst` is ldb-by-m, F-ordered.
 *
 * The default is to have src and dst of the same size (ldb=-1 means ldb=n).
 */
template<typename T>
void copy_slice_F(T* dst, const T* slice_ptr, const npy_intp n, const npy_intp m, const npy_intp s2, const npy_intp s1, npy_intp ldb=-1) {

    if (ldb == -1) {ldb = n;}

    for (npy_intp i = 0; i < n; i++) {
        for (npy_intp j = 0; j < m; j++) {
            dst[i + j*ldb] = *(slice_ptr + (i*s2/sizeof(T)) + (j*s1/sizeof(T)));  // == src[i*m + j]
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
 * Examples (matrix dimension m x n):
 *   C-contiguous src (s0=n, s1=1): src[i*n + j] -> dst[i*n + j]
 *      both sequential in the inner loop - effectively a partial memcpy.
 *   F-contiguous src (s0=1, s2=m): src[i + j*m] -> dst[i*n + j],
 *      column-sequential reads, row-sequential writes.
 */
template<typename T>
void copy_triangle_to_C(T *dst, const T *src, const npy_intp m, const npy_intp n, const npy_intp s0, const npy_intp s1, const char uplo) {
    if (uplo == 'L') {
        for (npy_intp i = 0; i < m; i++) {
            npy_intp stop = std::min(i + 1, n);
            for (npy_intp j = 0; j < stop; j++) {
                dst[i * n + j] = *(src + i * s0 + j * s1);
            }
        }
    } else {
        for (npy_intp i = 0; i < m; i++) {
            for (npy_intp j = i; j < n; j++) {
                dst[i * n + j] = *(src + i * s0 + j * s1);
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
    // NB: transpose for the F order
    if (uplo == 'U') {return norm1_sym_herm_lower(A, work, n);}
    else if (uplo == 'L') {return norm1_sym_herm_upper(A, work, n);}
    else {throw std::runtime_error("uplo at norms");}
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

template<typename T>
void
bandwidth(T* data, npy_intp n, npy_intp m, npy_intp* lower_band, npy_intp* upper_band)
{
    T zero = T(0.);

    Py_ssize_t lb = 0, ub = 0;
    for (Py_ssize_t c = 0; c < m-1; c++)
    {
        for (Py_ssize_t r = n-1; r > c + lb; r--)
        {
            if (data[c*n + r] != zero) { lb = r - c; break; }
        }
        if (c + lb + 1 > m) { break; }
    }
    for (Py_ssize_t c = m-1; c > 0; c--)
    {
        for (Py_ssize_t r = 0; r < c - ub; r++)
        {
            if (data[c*n + r] != zero) { ub = c - r; break; }

        }
        if (c <= ub) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}


/*
 * Overload of the original `bandwidth` function that allows to take into
 * account the strides of the matrix to avoid having to explicitly set a
 * flag regarding the ordering of the matrix.
 *
 * The addressing is done using `npy_intp` instead of `Py_ssize_t` for
 * consistency.
 */
template<typename T>
void
bandwidth_strided(T* data, npy_intp n, npy_intp m, npy_intp s1, npy_intp s2, npy_intp *lower_band, npy_intp *upper_band)
{
    T zero = T(0.);

    s1 = s1 / sizeof(T);
    s2 = s2 / sizeof(T);
    npy_intp lb = 0, ub = 0;
    for (npy_intp c = 0; c < m-1; c++) {
        for (npy_intp r = n-1; r > c + lb; r--) {
            if (data[c * s2 + r * s1] != zero) { lb = r - c; break; }
        }
        if (c + lb + 1 > m) { break; }
    }
    for (npy_intp c = m-1; c > 0; c--) {
        for (npy_intp r = 0; r < c - ub; r++) {
            if (data[c * s2 + r * s1] != zero) { ub = c - r; break; }
        }
        if (c <= ub) { break; }
    }
    *lower_band = lb;
    *upper_band = ub;
}


template<typename T>
void
detect_bandwidths(T* data, npy_intp ndim, npy_intp outer_size, npy_intp *shape, npy_intp *strides, npy_intp *kl, npy_intp *ku, npy_intp *ldab_max) {
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        T* slice_ptr = compute_slice_ptr(idx, data, ndim, shape, strides);

        bandwidth_strided(slice_ptr, shape[ndim-2], shape[ndim-1], strides[ndim-2], strides[ndim-1], &kl[idx], &ku[idx]);
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


template<typename T>
inline void
swap_cf(T* src, T* dst, const Py_ssize_t r, const Py_ssize_t c, const Py_ssize_t n)
{
    Py_ssize_t i, j, ith_row, r2, c2;
    T *bb = dst;
    T *aa = src;
    if ((r < 16) && (c < 16)) {
        for (j = 0; j < c; j++)
        {
            ith_row = 0;
            for (i = 0; i < r; i++) {
                bb[ith_row] = aa[i];
                ith_row += n;
            }
            aa += n;
            bb += 1;
        }
    } else {
        // If tall
        if (r > c)
        {
            r2 = r/2;
            swap_cf(src, dst, r2, c, n);
            swap_cf(src + r2, dst+(r2)*n, r-r2, c, n);
        } else {  // Nope
            c2 = c/2;
            swap_cf(src, dst, r, c2, n);
            swap_cf(src+(c2)*n, dst+c2, r, c-c2, n);
        }
    }
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
    s1 = s1 / sizeof(T);
    s2 = s2 / sizeof(T);
    npy_intp i, j;

    // main diagonal
    for (i = 0; i < n; i++) {
        ab[(i + 1) * ldab - kl - 1] = data[i * (s1 + s2)];
    }

    // lower bands
    for (i = 0; i < kl; i++) {
        for (j = 0; j < n - i - 1; j++) {
            ab[(j + 1) * ldab - kl + i] = data[(i + 1) * s1 + j * (s1 + s2)];
        }
    }

    // upper bands
    for (i = 0; i < ku; i++) {
        for (j = i + 1; j < n; j++) {
            ab[(j + 1) * ldab - kl - i - 2] = data[(i + 1) * s2 + (j - i - 1) * (s1 + s2)];
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
