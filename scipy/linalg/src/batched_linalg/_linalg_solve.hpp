/*
 * Templated loops for `linalg.solve`
 */
#pragma once

namespace sp_linalg {

// Dense array solve with getrf, gecon and getrs
template<typename T>
inline void solve_slice_general(
    CBLAS_INT N, CBLAS_INT NRHS, T *data, CBLAS_INT *ipiv, T *b_data, char trans, irwork_t<T> *irwork, T *work,
    SliceStatus& status
) {
    using real_type = real_of_t<T>;

    CBLAS_INT info;
    char norm = '1';
    real_type rcond;
    real_type anorm = norm1_(data, (npy_intp)N);

    lapack::getrf(N, N, data, N, ipiv, &info);

    status.lapack_info = (Py_ssize_t)info;
    if (info == 0){
        // getrf success, check the condition number
        lapack::gecon(norm, N, data, N, anorm, &rcond, work, irwork, &info);

        status.rcond = (double)rcond;
        if (info >= 0) {
            status.is_ill_conditioned = (rcond != rcond) || (rcond < std::numeric_limits<real_type>::epsilon());

            // finally, solve
            lapack::getrs(trans, N, NRHS, data, N, ipiv, b_data, N, &info);
            status.is_singular = (info > 0);
        }
    }
    else if (info > 0) {
        // trf detected singularity
        status.is_singular = 1;
    }
}


// triangular solve with trtrs
template<typename T>
inline void solve_slice_triangular(
    char uplo, char diag, CBLAS_INT N, CBLAS_INT NRHS, T *data,  T *b_data, char trans, T *work, irwork_t<T> *irwork,
    SliceStatus& status
) {
    using real_type = real_of_t<T>;

    CBLAS_INT info;
    char norm = '1';
    real_type rcond;

    lapack::trtrs(uplo, trans, diag, N, NRHS, data, N, b_data, N, &info);

    status.lapack_info = (Py_ssize_t)info;
    status.is_singular  = (info > 0);
    if(info >= 0) {
        lapack::trcon(norm, uplo, diag, N, data, N, &rcond, work, irwork, &info);

        if (info >= 0) {
            status.is_ill_conditioned = (rcond != rcond) || (rcond < std::numeric_limits<real_type>::epsilon());
            status.rcond = (double)rcond;
        }
    }
}


// Cholesky solve with potrf, pocon and potrs
template<typename T>
inline void solve_slice_cholesky(
    char uplo, CBLAS_INT N, CBLAS_INT NRHS, T *data, T *b_data, T* work, irwork_t<T> *irwork,
    SliceStatus& status
) {
    using real_type = real_of_t<T>;

    CBLAS_INT info;
    real_type rcond;
    real_type anorm = norm1_sym_herm(uplo, data, work, (npy_intp)N);

    lapack::potrf(uplo, N, data, N, &info);

    status.lapack_info = (Py_ssize_t)info;
    if (info == 0) {
        // potrf success
        lapack::pocon(uplo, N, data, N, anorm, &rcond, work, irwork, &info);

        if (info >= 0) {
            status.rcond = (double)rcond;
            status.is_ill_conditioned = (rcond != rcond) || (rcond < std::numeric_limits<real_type>::epsilon());

            // finally, solve
            lapack::potrs(uplo, N, NRHS, data, N, b_data, N, &info);
            status.is_singular = (info > 0);
        }
    }
    else if (info > 0) {
        // trf detected singularity
        status.is_singular = 1;
    }
}


// Symmetric/hermitian solve with sytrf/hetrf and sytrs/hetrs
template<typename T>
void solve_slice_sym_herm(
    char uplo, CBLAS_INT N, CBLAS_INT NRHS, T *data, T *b_data, CBLAS_INT *ipiv, T *work, irwork_t<T> *irwork, CBLAS_INT lwork,
    bool is_symm_not_herm,
    SliceStatus& status
) {
    using real_type = real_of_t<T>;

    CBLAS_INT info;
    real_type rcond;
    real_type anorm = norm1_sym_herm(uplo, data, work, (npy_intp)N);

    if constexpr (!is_complex_v<T>) {
        lapack::sytrf(uplo, N, data, N, ipiv, work, lwork, &info);
    } else {
        if (is_symm_not_herm) {
            lapack::sytrf(uplo, N, data, N, ipiv, work, lwork, &info);
        } else {
            lapack::hetrf(uplo, N, data, N, ipiv, work, lwork, &info);
        }
    }

    status.lapack_info = (Py_ssize_t)info;
    if (info == 0) {
        // {sy,he}trf success
        if constexpr (!is_complex_v<T>) {
            lapack::sycon(uplo, N, data, N, ipiv, anorm, &rcond, work, irwork, &info);
        } else {
            if (is_symm_not_herm) {
                lapack::sycon(uplo, N, data, N, ipiv, anorm, &rcond, work, &info);
            } else {
                lapack::hecon(uplo, N, data, N, ipiv, anorm, &rcond, work, &info);
            }
        }

        if (info >= 0) {
            status.rcond = (double)rcond;
            status.is_ill_conditioned = (rcond != rcond) || (rcond < std::numeric_limits<real_type>::epsilon());

            // finally, solve
            if constexpr (!is_complex_v<T>) {
                lapack::sytrs(uplo, N, NRHS, data, N, ipiv, b_data, N, &info);
            } else {
                if (is_symm_not_herm) {
                    lapack::sytrs(uplo, N, NRHS, data, N, ipiv, b_data, N, &info);
                } else {
                    lapack::hetrs(uplo, N, NRHS, data, N, ipiv, b_data, N, &info);
                }
            }

            status.is_singular = (info > 0);
        }
    }
    else if (info > 0) {
        // trf detected singularity
        status.is_singular = 1;
    }
}


template<typename T>
void solve_slice_tridiag(
    CBLAS_INT N, CBLAS_INT NRHS, T *data, CBLAS_INT *ipiv, char trans,
    T *b_data,
    T *work, T *work2, irwork_t<T> *irwork,
    SliceStatus& status
) {
    using real_type = real_of_t<T>;
    // work is 4*n, is for dl, d, du, du2
    // work2 is 2*n, is for trcon's work array

    // copy the diagonals of `data`, use `work` for storage
    T *du2 = work;
    T *du = &work[N];
    T *d = &work[2*N];
    T *dl = &work[3*N];
    to_tridiag(data, N, du, d, dl);

    CBLAS_INT info;
    char norm = '1';
    real_type rcond;
    real_type anorm = norm1_tridiag(dl, d, du, work2, (npy_intp)N);

    lapack::gttrf(N, dl, d, du, du2, ipiv, &info);

    status.lapack_info = (Py_ssize_t)info;
    if (info == 0){
        // gttrf success, check the condition number
        if constexpr (!is_complex_v<T>) {
            lapack::gtcon(norm, N, dl, d, du, du2, ipiv, anorm, &rcond, work2, irwork, &info);
        } else {
            lapack::gtcon(norm, N, dl, d, du, du2, ipiv, anorm, &rcond, work2, &info);
        }

        status.rcond = (double)rcond;
        if (info >= 0) {
            status.is_ill_conditioned = (rcond != rcond) || (rcond < std::numeric_limits<real_type>::epsilon());

            // finally, solve
            lapack::gttrs(trans, N, NRHS, dl, d, du, du2, ipiv, b_data, N, &info);
            status.is_singular = (info > 0);
        }
    }
    else if (info > 0) {
        // trf detected singularity
        status.is_singular = 1;
    }
}

// Banded array solve, assumes the data is fed in in the banded structure `ab` already.
template<typename T>
inline void solve_slice_banded(
    char trans, CBLAS_INT N, CBLAS_INT NRHS, T *ab, CBLAS_INT *ipiv, T *b_data, T *work2, irwork_t<T> *irwork,
    CBLAS_INT kl, CBLAS_INT ku, SliceStatus &status
) {
    using real_type = real_of_t<T>;

    CBLAS_INT ldab = 2 * kl + ku + 1;

    CBLAS_INT info;
    char norm = '1';
    real_type rcond;
    real_type anorm = norm1_banded(ab, kl, ku, work2, N);

    lapack::gbtrf(N, N, kl, ku, ab, ldab, ipiv, &info);
    status.lapack_info = (Py_ssize_t)info;
    if (info == 0) {
        // gbtrf success, check condition number
        lapack::gbcon(norm, N, kl, ku, ab, ldab, ipiv, anorm, &rcond, work2, irwork, &info);

        status.rcond = (double)rcond;
        if (info >= 0) {
            status.is_ill_conditioned = (rcond != rcond) || (rcond < std::numeric_limits<real_type>::epsilon());

            // finally, solve
            lapack::gbtrs(trans, N, kl, ku, NRHS, ab, ldab, ipiv, b_data, N, &info);
            status.is_singular = (info > 0);
        }
    }
    else if (info > 0) {
        // trf detected singularity
        status.is_singular = 1;
    }
}


// Diagonal array solve
template<typename T>
inline void solve_slice_diagonal(
    CBLAS_INT N, CBLAS_INT NRHS, T *data, T *b_data, SliceStatus& status
) {
    using real_type = real_of_t<T>;

    T zero(0.), one(1.);
    real_type maxa(0.), maxinva(0.);

    for (npy_intp j=0; j<N; j++) {
        T ajj = data[j*N + j];

        status.is_singular = (ajj == zero);
        if (status.is_singular) {
            status.lapack_info = j;
            return;
        }

        T inv_ajj = one / ajj;
        for (npy_intp i=0; i<NRHS; i++) {
            b_data[j + i*N] *= inv_ajj;
        }

        // condition number
        real_type absa = std::abs(ajj), absinva = std::abs(inv_ajj);

        if(absa > maxa) {maxa = absa;}
        if(absinva > maxinva) {maxinva = absinva;}
    }
    double cond = (double)maxa * (double)maxinva;
    double rcond = 1.0 / cond;
    status.is_ill_conditioned = (rcond != rcond) || (rcond < std::numeric_limits<real_type>::epsilon());
    status.rcond = rcond;
}


// Specific implementation for the case of `solve(..., assume_a=banded)` to avoid
// cluttering the main loop in `_solve` by branching to this function early.
template<typename T>
BatchStatus
_solve_assume_banded(PyArrayObject *ap_Am, PyArrayObject *ap_b, T *ret_data, char trans, int overwrite_a, int overwrite_b, SliceStatus slice_status, SliceStatusVec &vec_status)
{
    // --------------------------------------------------------------------
    // Input Array Attributes, duplicate of `_solve()`
    // --------------------------------------------------------------------
    T* Am_data = (T *)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_Am);      // Array shape
    npy_intp n = shape[ndim - 1];                // Slice size
    npy_intp* strides = PyArray_STRIDES(ap_Am);
    // Get the number of slices to traverse if more than one; np.prod(shape[:-2])
    npy_intp outer_size = 1;
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    T *bm_data = (T *)PyArray_DATA(ap_b);
    npy_intp ndim_b = PyArray_NDIM(ap_b);
    npy_intp *shape_b = PyArray_SHAPE(ap_b);
    npy_intp *strides_b = PyArray_STRIDES(ap_b);
    npy_intp nrhs = PyArray_DIM(ap_b, ndim_b -1); // Number of right-hand-sides

    CBLAS_INT intn, int_nrhs;
    if (!checked_narrow<CBLAS_INT>(n, &intn) || !checked_narrow<CBLAS_INT>(nrhs, &int_nrhs)) {
        return BatchStatus::too_large;
    }

    // Bandwidth detection per slice. Required to first do a pass to find the
    // maximal `kl` and `ku` to find the minimal size the array will need to
    // have. To avoid having to call `bandwidth` twice per slice, the results
    // are stored in `ks`: lower bandwidths first, then upper bandwidths.
    // element counts are products of dimensions of existing arrays, so they fit size_t
    auto ks = alloc_workspace<int64_t>(2 * (std::size_t)outer_size);
    auto ipiv = alloc_workspace<CBLAS_INT>(n);
    auto irwork = alloc_workspace<irwork_t<T>>((is_complex_v<T> ? 3 : 1) * (std::size_t)n);
    if (!ks || !ipiv || !irwork) { return BatchStatus::no_memory; }

    int64_t ldab_max = 0;
    int64_t *kls = ks.get();
    int64_t *kus = ks.get() + outer_size;
    detect_bandwidths(Am_data, ndim, outer_size, shape, strides, kls, kus, &ldab_max);

    // every slice's ldab = 2*kl + ku + 1 <= ldab_max is passed to LAPACK
    if (ldab_max > CBLAS_INT_MAX) { return BatchStatus::too_large; }

    /*
     * Chop up buffer in parts:
     *
     *   ldab_max * n       3 * n       b_data_size
     * |---------------|-------------|---------------|
     * ^               ^             ^
     * ab              work          b_data
     *
     * - `ab` is a buffer holding the "padded" banded form of matrix `a`
     * - `work` is needed for `gbcon` (fixed size)
     * - `b_data` is a buffer for the rhs of the system, not needed if `overwrite_b` is set (size = 0 then)
     */
    // ldab_max <= 3*n, so these counts are bounded by the size of `a` as well
    const std::size_t ab_count = (std::size_t)ldab_max * n, work_count = 3 * (std::size_t)n;
    const std::size_t b_count = overwrite_b ? 0 : (std::size_t)n * nrhs;
    auto buffer = alloc_workspace<T>(ab_count + work_count + b_count);
    if (!buffer) { return BatchStatus::no_memory; }

    // Chop up buffer
    T *cursor = buffer.get();
    T *ab = cursor;   cursor += ab_count;
    T *work = cursor; cursor += work_count;
    T *b_data = overwrite_b ? bm_data : cursor;

    // Main loop traversal, taken from `_solve`
    // A failing slice is recorded in vec_status and ends the loop; the call itself still returns ok.
    for (npy_intp idx = 0; idx < outer_size; idx++) {

        T* slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);

        // Directly take into banded storage
        npy_intp ldab = 2 * kls[idx] + kus[idx] + 1;
        to_banded(slice_ptr, n, kls[idx], kus[idx], ldab, ab, strides[ndim-2], strides[ndim-1]);

        if (!overwrite_b) {
            T* slice_ptr_b = compute_slice_ptr(idx, bm_data, ndim, shape_b, strides_b);
            copy_strided_to_f(b_data, slice_ptr_b, n, nrhs, strides_b[ndim-2], strides_b[ndim-1], n);
        }
        // NB. `overwrite_b` is only set when the input is 2D F-ordered, so no copy needed then.
        // A generalization to ndim > 2 will also require computing the pointer here.

        // structure is known to be banded
        init_status(slice_status, idx, St::BANDED);
        // kl, ku <= ldab_max <= CBLAS_INT_MAX, checked above
        solve_slice_banded(trans, intn, int_nrhs, ab, ipiv.get(), b_data, work, irwork.get(),
                           (CBLAS_INT)kls[idx], (CBLAS_INT)kus[idx], slice_status);

        if (_detect_problems(slice_status, vec_status) != 0) { return BatchStatus::ok; }

        if (!overwrite_b) {
            // Put result in C-order in return buffer
            copy_slice_F_to_C(&ret_data[idx * n * nrhs], b_data, n, nrhs);
        } // If `overwrite_b` is true the result is already set.
    }

    return BatchStatus::ok;
}

template<typename T>
BatchStatus
_solve(PyArrayObject* ap_Am, PyArrayObject *ap_b, T* ret_data, St structure, int lower, int transposed, int overwrite_a, int overwrite_b, SliceStatusVec& vec_status)
{
    char trans = transposed ? 'T' : 'N';
    int64_t lower_band = 0, upper_band = 0;
    bool is_symm = false, is_herm = false;
    char uplo = lower ? 'L' : 'U';
    St slice_structure = St::NONE;
    bool posdef_fallback = true;
    SliceStatus slice_status;

    // branch early for `assume_a = banded`
    if (structure == St::BANDED) {
        return _solve_assume_banded(ap_Am, ap_b, ret_data, trans, overwrite_a, overwrite_b, slice_status, vec_status);
    }


    // --------------------------------------------------------------------
    // Input Array Attributes
    // --------------------------------------------------------------------
    T* Am_data = (T *)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_Am);      // Array shape
    npy_intp n = shape[ndim - 1];                // Slice size
    npy_intp* strides = PyArray_STRIDES(ap_Am);
    // Get the number of slices to traverse if more than one; np.prod(shape[:-2])
    npy_intp outer_size = 1;
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    T *bm_data = (T *)PyArray_DATA(ap_b);
    npy_intp ndim_b = PyArray_NDIM(ap_b);
    npy_intp *shape_b = PyArray_SHAPE(ap_b);
    npy_intp *strides_b = PyArray_STRIDES(ap_b);
    npy_intp nrhs = PyArray_DIM(ap_b, ndim_b -1); // Number of right-hand-sides


    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    CBLAS_INT intn, int_nrhs, info;
    if (!checked_narrow<CBLAS_INT>(n, &intn) || !checked_narrow<CBLAS_INT>(nrhs, &int_nrhs)) {
        return BatchStatus::too_large;
    }

    T tmp = 0.0;
    lapack::sytrf(uplo, intn, NULL, intn, NULL, &tmp, -1, &info);
    if (info != 0) { return BatchStatus::lapack_error; }

    // _calc_lwork returns -1 when the queried lwork does not fit CBLAS_INT
    CBLAS_INT lwork_sytrf = _calc_lwork(tmp);
    if (lwork_sytrf < 0) { return BatchStatus::too_large; }

    // gecon needs lwork of at least 4*n
    CBLAS_INT lwork;
    if (!checked_narrow<CBLAS_INT>(std::max(4 * (std::size_t)n, (std::size_t)lwork_sytrf), &lwork)) {
        return BatchStatus::too_large;
    }

    /*
     * Chop the buffer into parts:
     *
     *    nn       nb        2n      lwork
     * |--------|--------|--------|---------|
     * ^        ^        ^        ^
     * data     data_b   work2    work
     *
     * - data is for A (lhs), absent with overwrite_a
     * - data_b is for b (rhs), absent with overwrite_b
     * - work2 is for the tridiag solver, trcon's work array
     * - work is for all other LAPACK functions
     */
    // element counts are products of dimensions of existing arrays, so they fit size_t
    const std::size_t nn = (std::size_t)n * n, nb = (std::size_t)n * nrhs, n2 = 2 * (std::size_t)n;
    auto buffer = alloc_workspace<T>((overwrite_a ? 0 : nn) + (overwrite_b ? 0 : nb) + n2 + lwork);
    auto ipiv = alloc_workspace<CBLAS_INT>(n);
    auto irwork = alloc_workspace<irwork_t<T>>((is_complex_v<T> ? 3 : 1) * (std::size_t)n);   // {ge,po,tr}con need rwork or iwork
    if (!buffer || !ipiv || !irwork) { return BatchStatus::no_memory; }

    T *cursor = buffer.get();
    T *data = Am_data, *data_b = ret_data;
    if (!overwrite_a) { data = cursor; cursor += nn; }
    if (!overwrite_b) { data_b = cursor; cursor += nb; }
    T *work2 = cursor; cursor += n2;
    T *work = cursor;

    // normalize the structure detection inputs
    if (structure == St::POS_DEF) {
        posdef_fallback = false;
    }
    else if (structure == St::SYM) {
        is_symm = true;
    }
    else if (structure == St::HER) {
        is_herm = true;
    }
    if (structure == St::LOWER_TRIANGULAR) {
        uplo = 'L';
    }
    else if (structure == St::UPPER_TRIANGULAR) {
        uplo = 'U';
    }

    // Main loop to traverse the slices
    // A failing slice is recorded in vec_status and ends the loop; the call itself still returns ok.
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);

        if (!overwrite_a) {
            // copy the slice
            copy_strided_to_f(data, slice_ptr, n, n, strides[ndim-2], strides[ndim-1], n);
        }

        if (!overwrite_b) {
            // copy the r.h.s, too;
            T *slice_ptr_b = compute_slice_ptr(idx, bm_data, ndim_b, shape_b, strides_b);
            copy_strided_to_f(data_b, slice_ptr_b, n, nrhs, strides_b[ndim-2], strides_b[ndim-1], n);
        }

        // detect the structure if not given
        slice_structure = structure;

        if (slice_structure == St::NONE) {
            // Get the bandwidth of the slice
            // `data` is the F-ordered copy of the slice
            bandwidth_slice(data, n, n, (int64_t)sizeof(T), n * (int64_t)sizeof(T), &lower_band, &upper_band);

            if ((upper_band == 0) && (lower_band == 0)) {
                slice_structure = St::DIAGONAL;
            }
            else if ((upper_band == 1) && (lower_band == 1) && (n > 3)) {
                slice_structure = St::TRIDIAGONAL;
            }
            else if(lower_band == 0) {
                slice_structure = St::UPPER_TRIANGULAR;
                uplo = 'U';
            } else if (upper_band == 0) {
                slice_structure = St::LOWER_TRIANGULAR;
                uplo = 'L';
            } else {
                // Check if symmetric/hermitian
                std::tie(is_symm, is_herm) = is_sym_or_herm(data, n);
                if (is_herm || (is_symm && !is_complex_v<T>)) {
                    // either real symmetric or complex hermitian; try Cholesky first,
                    // fall back to sym/her if it fails
                    if (!overwrite_a) {
                        slice_structure = St::POS_DEF;
                    }
                    else {
                        // working in-place: cannot try Cholesky, have to use the
                        // right structure straight away
                        slice_structure = is_symm ? St::SYM : St::HER;
                    }
                }
                else if (is_symm && is_complex_v<T>) {
                    // complex symmetric, not hermitian
                    slice_structure = St::SYM;
                }
                else {
                    // give up auto-detection
                    slice_structure = St::GENERAL;
                }
            }
        }

        init_status(slice_status, idx, slice_structure);

        switch(slice_structure) {
            case St::DIAGONAL:
            {
                solve_slice_diagonal(intn, int_nrhs, data, data_b, slice_status);
                if (_detect_problems(slice_status, vec_status) != 0) {
                    // fail fast and loud
                    return BatchStatus::ok;
                }
                break;
            }
            case St::TRIDIAGONAL:
            {
                solve_slice_tridiag(intn, int_nrhs, data, ipiv.get(), trans, data_b, work, work2, irwork.get(), slice_status);
                if (_detect_problems(slice_status, vec_status) != 0) {
                    return BatchStatus::ok;
                }
                break;
            }
            case St::UPPER_TRIANGULAR:
            case St::LOWER_TRIANGULAR:
            {
                char diag = 'N';
                solve_slice_triangular(uplo, diag, intn, int_nrhs, data, data_b, trans, work, irwork.get(), slice_status);
                if (_detect_problems(slice_status, vec_status) != 0) {
                    return BatchStatus::ok;
                }
                break;
            }
            case St::POS_DEF:
            {
                solve_slice_cholesky(uplo, intn, int_nrhs, data, data_b, work, irwork.get(), slice_status);

                if ((slice_status.lapack_info == 0) || (!slice_status.is_singular) ) {
                    // success (maybe ill-conditioned)
                    if(slice_status.is_ill_conditioned) {
                        vec_status.push_back(slice_status);
                    }
                    break;
                }
                else { // potrf failed
                    if(posdef_fallback) {
                        // restore
                        copy_strided_to_f(data, slice_ptr, n, n, strides[ndim-2], strides[ndim-1], n);
                        init_status(slice_status, idx, slice_structure);

                        // no break: fall back to the general solver
                    }
                    else {
                        // potrf failed but no fallback
                        vec_status.push_back(slice_status);
                        break;
                    }
                }
            }
            case St::SYM:  // pos def fails, fall through to here
            case St::HER:
            {
                solve_slice_sym_herm(uplo, intn, int_nrhs, data, data_b, ipiv.get(), work, irwork.get(), lwork, (is_symm && !is_herm), slice_status);
                if (_detect_problems(slice_status, vec_status) != 0) {
                    return BatchStatus::ok;
                }
                break;
            }
            default:
            {
                // general matrix solve
                solve_slice_general(intn, int_nrhs, data, ipiv.get(), data_b, trans, irwork.get(), work, slice_status);
                if (_detect_problems(slice_status, vec_status) != 0) {
                    return BatchStatus::ok;
                }
            }
        }

        if (!overwrite_b) {
            // Swap back to the C order
            copy_slice_F_to_C(&ret_data[idx*n*nrhs], data_b, n, nrhs);
        }
    }

    return BatchStatus::ok;
}


} // namespace sp_linalg
