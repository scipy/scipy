/*
 * Templated loops for `linalg.inv`
 */
#pragma once

namespace sp_linalg {

// Dense array inversion with getrf, gecon and getri
template<typename T>
void invert_slice_general(
    CBLAS_INT N, T *data, CBLAS_INT *ipiv, irwork_t<T> *irwork, T *work, CBLAS_INT lwork,
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

            // finally, invert
            lapack::getri(N, data, N, ipiv, work, lwork, &info);
            status.is_singular = (info > 0);
        }
    }
    else if (info > 0) {
        // trf detected singularity
        status.is_singular = 1;
    }
}


///////////////////////////

// Symmetric/hermitian array inversion with potrf, pocon and potri
template<typename T>
void invert_slice_cholesky(
    char uplo, CBLAS_INT N, T *data, T* work, irwork_t<T> *irwork,
    SliceStatus& status
) {
    using real_type = real_of_t<T>;

    CBLAS_INT info;
    real_type anorm = norm1_sym_herm(uplo, data, work, (npy_intp)N);

    real_type rcond;

    lapack::potrf(uplo, N, data, N, &info);

    status.lapack_info = (Py_ssize_t)info;
    if (info == 0) {
        // potrf success
        lapack::pocon(uplo, N, data, N, anorm, &rcond, work, irwork, &info);

        if (info >= 0) {
            status.rcond = (double)rcond;
            status.is_ill_conditioned = (rcond != rcond) || (rcond < std::numeric_limits<real_type>::epsilon());

            // finally, invert
            lapack::potri(uplo, N, data, N, &info);
            status.is_singular = (info > 0);
        }
    }
    else if (info > 0) {
        // trf detected singularity
        status.is_singular = 1;
    }
}

// Symmetric/hermitian array inversion with sytrf/hetrf and sytri/hetri
template<typename T>
void invert_slice_sym_herm(
    char uplo, CBLAS_INT N, T *data, CBLAS_INT *ipiv, T *work, irwork_t<T> *irwork, CBLAS_INT lwork,
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
        if(is_symm_not_herm) {
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

            // finally, invert
            if constexpr (!is_complex_v<T>) {
                lapack::sytri(uplo, N, data, N, ipiv, work, &info);
            } else {
                if (is_symm_not_herm) {
                    lapack::sytri(uplo, N, data, N, ipiv, work, &info);
                } else {
                    lapack::hetri(uplo, N, data, N, ipiv, work, &info);
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


// triangular array inversion with trtri
template<typename T>
void invert_slice_triangular(
    char uplo, char diag, CBLAS_INT N, T *data, T *work, irwork_t<T> *irwork,
    SliceStatus& status
) {
    using real_type = real_of_t<T>;

    CBLAS_INT info;
    char norm = '1';
    real_type rcond;

    lapack::trtri(uplo, diag, N, data, N, &info);
    status.is_singular  = (info > 0);

    status.lapack_info = (Py_ssize_t)info;
    if(info >= 0) {
        lapack::trcon(norm, uplo, diag, N, data, N, &rcond, work, irwork, &info);

        if (info >= 0) {
            status.is_ill_conditioned = (rcond != rcond) || (rcond < std::numeric_limits<real_type>::epsilon());
            status.rcond = (double)rcond;
        }
    }
}


// Diagonal array inversion
template<typename T>
inline void invert_slice_diagonal(
    CBLAS_INT N, T *data, SliceStatus& status
) {
    using real_type = real_of_t<T>;

    T zero(0.), one(1.);
    real_type maxa(0.), maxinva(0.);

    for (npy_intp j=0; j<N; j++) {
        T ajj = data[j*N + j];

        status.is_singular = (ajj == zero);
        if (status.is_singular) {
            status.lapack_info  = j;
            return;
        }

        T inv_ajj = one / ajj;
        data[j*N + j] = inv_ajj;

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


template<typename T>
BatchStatus
_inverse(PyArrayObject* ap_Am, T* ret_data, St structure, int lower, int overwrite_a, SliceStatusVec& vec_status)
{
    int64_t lower_band = 0, upper_band = 0;
    bool is_symm = false, is_herm = false;
    char uplo = lower ? 'L' : 'U';
    St slice_structure = St::NONE;
    bool posdef_fallback = true;
    SliceStatus slice_status;

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

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    CBLAS_INT intn, info;
    if (!checked_narrow<CBLAS_INT>(n, &intn)) { return BatchStatus::too_large; }

    T tmp = 0.0;
    lapack::getri(intn, NULL, intn, NULL, &tmp, -1, &info);
    if (info != 0) { return BatchStatus::lapack_error; }
    CBLAS_INT lwork_getri = _calc_lwork(tmp, 1.01);

    lapack::sytrf(uplo, intn, NULL, intn, NULL, &tmp, -1, &info);
    if (info != 0) { return BatchStatus::lapack_error; }
    CBLAS_INT lwork_sytrf = _calc_lwork(tmp);

    // _calc_lwork returns -1 when the queried lwork does not fit CBLAS_INT
    if ((lwork_getri < 0) || (lwork_sytrf < 0)) { return BatchStatus::too_large; }

    // gecon needs lwork of at least 4*n
    CBLAS_INT lwork;
    if (!checked_narrow<CBLAS_INT>(std::max({4 * (std::uint64_t)n, (std::uint64_t)lwork_getri, (std::uint64_t)lwork_sytrf}), &lwork)) {
        return BatchStatus::too_large;
    }

    /*
     * Finally, we can start allocating memory.
     *
     * The key point is that LAPACK always operates on F-ordered arrays.
     * The memory strategy thus depends on the `overwrite_a` value.
     *
     * For `overwrite_a=False` (default), we:
     *   - allocate a temp buffer (`data` below) once
     *   - for each slice, we
     *       - copy-and-transpose the slice into the temp buffer,
     *       - feed the buffer to LAPACK
     *       - copy-and-transpose the result back to the C-ordered result array
     *
     * For `overwrite_a=True`, we assume that
     *   - `ret_data` may point to the same memory as `ap_Am` array
     *   - the caller had ensured that the input array is Fortran-ordered
     *   - the caller wants to get the result also Fortran ordered
     *
     * It's a caller's responsibility to make sure that these pre-conditions are met,
     * none of them is checked here.
     * Therefore, if `overwrite_a = True`, we skip the copy-and-transpose steps above,
     * and `ret_data` will simply contain the result from the LAPACK call.
     *
     */
    // each count is bounded by an existing array or by CBLAS_INT, so the 64-bit sums are exact
    const std::uint64_t nn = (std::uint64_t)n * n;
    auto buffer = alloc_workspace<T>((overwrite_a ? 0 : nn) + lwork);
    auto ipiv = alloc_workspace<CBLAS_INT>(n);
    auto irwork = alloc_workspace<irwork_t<T>>((is_complex_v<T> ? 3 : 1) * (std::uint64_t)n);   // {ge,po,tr}con need rwork or iwork

    // Chop buffer into parts: the slice (unless working in place), then LAPACK's work array
    T *cursor = buffer.get();
    T *data = ret_data;
    if (!overwrite_a) { data = cursor; cursor += nn; }
    T *work = cursor;

    /*
     * Normalize the structure detection inputs.
     */
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

    /*
     * Main loop to traverse the slices.
     * A failing slice is recorded in vec_status and ends the loop; the call itself still returns ok.
     */
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);

        if (!overwrite_a) {
            copy_strided_to_f(data, slice_ptr, n, n, strides[ndim-2], strides[ndim-1], n);
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
            else if(lower_band == 0) {
                slice_structure = St::UPPER_TRIANGULAR;
                uplo = 'U';
            } else if (upper_band == 0) {
                slice_structure = St::LOWER_TRIANGULAR;
                uplo = 'L';
            } else {
                // Check if symmetric/hermitian
                std::tie(is_symm, is_herm) = is_sym_or_herm(data, n);

                if constexpr (!is_complex_v<T>) {
                    // Real: is_symm and is_herm are always equal
                    if (is_symm) {
                        /*
                         * If working on a copy (overwrite_a is False):
                         *    try Cholesky first, fall back to sytrf if it fails
                         * If working in-place, do the inversion in one go,
                         *    (if Cholesky failed, it already destroyed the input)
                         */
                        slice_structure = overwrite_a ? St::SYM : St::POS_DEF ;
                    }
                    else {
                        slice_structure = St::GENERAL;
                    }
                }
                else {
                    // Complex
                    if (!is_symm && !is_herm) {
                        slice_structure = St::GENERAL;
                    }
                    else if (is_herm) {
                        // Hermitian (may also be symmetric if entries are real)
                        // try Cholesky first, fall back to hetrf if it fails
                        slice_structure = overwrite_a ? St::HER : St::POS_DEF ;
                    }
                    else {
                        // is_symm && !is_herm: complex symmetric, not hermitian
                        slice_structure = St::SYM;
                    }
                }
            }
        }

        init_status(slice_status, idx, slice_structure);

        // Use the appropriate LAPACK function for the `slice_structure`.
        switch(slice_structure) {
            case St::DIAGONAL:
            {
                invert_slice_diagonal(intn, data, slice_status);
                if (_detect_problems(slice_status, vec_status) != 0) {
                    // fail fast and loud
                    return BatchStatus::ok;
                }
                break;
            }
            case St::UPPER_TRIANGULAR:
            case St::LOWER_TRIANGULAR:
            {
                char diag = 'N';
                invert_slice_triangular(uplo, diag, intn, data, work, irwork.get(), slice_status);

                if (_detect_problems(slice_status, vec_status) != 0) {
                    return BatchStatus::ok;
                }
                zero_other_triangle(uplo, data, intn);
                break;
            }
            case St::POS_DEF:
            {
                invert_slice_cholesky(uplo, intn, data, work, irwork.get(), slice_status);

                if ((slice_status.lapack_info == 0) || (!slice_status.is_singular) ) {
                    // success (maybe ill-conditioned)
                    if(slice_status.is_ill_conditioned) {
                        vec_status.push_back(slice_status);
                    }
                    fill_other_triangle(uplo, data, intn);
                    break;
                }
                else { // potrf failed
                    if(posdef_fallback) {
                        // restore
                        copy_strided_to_f(data, slice_ptr, n, n, strides[ndim-2], strides[ndim-1], n);
                        init_status(slice_status, idx, slice_structure);

                        // no break: fall back to the symmetric solver
                    }
                    else {
                        // potrf failed but no fallback
                        vec_status.push_back(slice_status);
                        break;
                    }
                }
            }
            case St::SYM:     // NB: if POS_DEF failed, fall-through to here
            case St::HER:
            {
                if constexpr (!is_complex_v<T>) {
                    // Real: always use sytrf/sytri
                    invert_slice_sym_herm(uplo, intn, data, ipiv.get(), work, irwork.get(), lwork, true, slice_status);
                }
                else {
                    // Complex: use sytrf if symmetric-only, hetrf if hermitian
                    invert_slice_sym_herm(uplo, intn, data, ipiv.get(), work, irwork.get(), lwork, !is_herm, slice_status);
                }

                if (_detect_problems(slice_status, vec_status) != 0) {
                    return BatchStatus::ok;
                }

                if constexpr (!is_complex_v<T>) {
                    // Real symmetric
                    fill_other_triangle_noconj(uplo, data, intn);
                }
                else {
                    // Complex: depends on whether symmetric or hermitian
                    if (!is_herm) {
                        fill_other_triangle_noconj(uplo, data, intn);
                    }
                    else {
                        fill_other_triangle(uplo, data, intn);
                    }
                }
                break;
            }
            default:
            {
                // general matrix inverse
                invert_slice_general(intn, data, ipiv.get(), irwork.get(), work, lwork, slice_status);

                if (_detect_problems(slice_status, vec_status) != 0) {
                    return BatchStatus::ok;
                }
            }
        } // end of `switch(slice_structure)`

        if (!overwrite_a) {
            // transpose the F-ordered result into the C-ordered output slice
            copy_slice_F_to_C(&ret_data[idx*n*n], data, n, n, n);
        }
    } // end of `for(idx=...)`

    return BatchStatus::ok;
}

} // namespace sp_linalg
