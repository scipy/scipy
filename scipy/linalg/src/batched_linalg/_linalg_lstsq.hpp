/*
 * Templated loops for `linalg.lstsq`
 */
#pragma once

namespace sp_linalg {

enum class LstsqDriver { gelsd, gelsy, gelss };


/**
 * @brief Batched least squares, `linalg.lstsq`: minimise ||a @ x - b|| for every slice.
 *
 * @param[in]  ap_Am        Input array of shape (..., m, n), any strides.
 * @param[in]  ap_b         Right-hand sides of shape (..., m, nrhs), any strides, with the batch
 *                          shape of @p ap_Am.
 * @param[out] ap_S         Singular values, real (..., min(m, n)) array; NULL for gelsy.
 * @param[out] ap_x         C-contiguous (..., max(m, n), nrhs) buffer. The first n rows of each
 *                          slice receive x, the remaining rows the residual vectors. With
 *                          @p overwrite_b this is @p ap_b itself, F-contiguous with m >= n.
 * @param[out] ap_rank      C-contiguous (...) int64 effective ranks.
 * @param[in]  rcond        Singular values below rcond times the largest are treated as zero.
 * @param[in]  driver       The LAPACK routine.
 * @param[in]  overwrite_a  If true, @p ap_Am is 2-D and F-contiguous and is used as the
 *                          workspace for a.
 * @param[in]  overwrite_b  If true, @p ap_b is 2-D and F-contiguous with m >= n and receives
 *                          the solution in place.
 * @param[out] vec_status   Receives the status of the first failing slice, if any.
 * @return ok; too_large if a dimension or workspace size does not fit CBLAS_INT; no_memory;
 *         lapack_error if the workspace query fails.
 */
template<typename T>
BatchStatus
_lstsq(PyArrayObject *ap_Am, PyArrayObject *ap_b, PyArrayObject *ap_S, PyArrayObject *ap_x, PyArrayObject *ap_rank,
       f64 rcond, LstsqDriver driver, int overwrite_a, int overwrite_b, SliceStatusVec& vec_status)
{
    using real_type = real_of_t<T>; // f32 if T==c64 etc
    SliceStatus slice_status;

    // --------------------------------------------------------------------
    // Input array attributes
    // --------------------------------------------------------------------
    T *Am_data = (T *)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    npy_intp *strides = PyArray_STRIDES(ap_Am);
    npy_intp m = shape[ndim - 2], n = shape[ndim - 1];

    // the entry point checked that a and b share ndim and the batch shape
    T *bm_data = (T *)PyArray_DATA(ap_b);
    npy_intp *shape_b = PyArray_SHAPE(ap_b);
    npy_intp *strides_b = PyArray_STRIDES(ap_b);
    npy_intp nrhs = shape_b[ndim - 1];

    npy_intp outer_size = 1;
    for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i]; }

    const npy_intp min_mn = m < n ? m : n, max_mn = m < n ? n : m;

    CBLAS_INT intm, intn, int_nrhs, ldb;
    if (!checked_narrow<CBLAS_INT>(m, &intm) || !checked_narrow<CBLAS_INT>(n, &intn) ||
        !checked_narrow<CBLAS_INT>(nrhs, &int_nrhs) || !checked_narrow<CBLAS_INT>(max_mn > 1 ? max_mn : 1, &ldb)) {
        return BatchStatus::too_large;
    }
    const CBLAS_INT lda = intm > 1 ? intm : 1;
    const real_type r_rcond = (real_type)rcond;
    CBLAS_INT rank = 0, info = 0;

    // The LAPACK routines evaluate their workspace sizes in CBLAS_INT, already in the workspace
    // query, so the largest of those sizes is computed here first, after Reference-LAPACK.
    // It grows with min(m, n)*nrhs, ldb*nrhs (the comparisons in gelss), the block size times
    // nrhs, the m*m factor of the path that starts with an LQ factorization, and, in the real
    // workspace of the complex gelsd, 3*smlsiz*nrhs + n*(1 + nrhs). The sizes are counted in
    // uint64_t, which holds them exactly for every CBLAS_INT operand below 2**31.
    {
        const std::uint64_t um = m, un = n, unrhs = nrhs, uk = min_mn, umx = max_mn, uldb = ldb;
        // ILAENV(1, ...) block size of the routines each driver calls: 64 for the ?TZRZF and ?ORMRZ
        // of gelsy, 32 for the others
        const std::uint64_t nb = (driver == LstsqDriver::gelsy) ? 64 : 32;
        const std::uint64_t smlsiz = 25;    // ILAENV(9, ...), the size of the subproblems in gelsd
        // NLVL = MAX(INT(LOG2(min(m, n) / (smlsiz + 1))) + 1, 0), truncating towards zero
        const long lvl = (uk > 0) ? (long)std::log2(double(uk) / double(smlsiz + 1)) + 1 : 0;
        const std::uint64_t nlvl = lvl > 0 ? (std::uint64_t)lvl : 0;
        // n >= MNTHR = ILAENV(6, ...) = INT(1.6 * min(m, n)) selects the LQ path when m < n
        const bool lq_path = (m < n) && un >= (std::uint64_t)(float(um) * 1.6f);

        // the block-size terms: (m + n)*nb from the bidiagonal or QR stage, nrhs*nb from applying
        // the orthogonal factors to b
        std::uint64_t need = 3*umx + 2*uk + nb*std::max(um + un + 1, unrhs);
        if (driver == LstsqDriver::gelsd) {
            const std::uint64_t wlalsd = 9*uk + 2*uk*smlsiz + 8*uk*nlvl + uk*unrhs + (smlsiz + 1)*(smlsiz + 1);
            need = std::max(need, 3*uk + wlalsd);
            if (lq_path) { need = std::max(need, um*um + 4*um + std::max(wlalsd, um*unrhs + um)); }
            if constexpr (is_complex_v<T>) {
                const std::uint64_t lrwork = 10*uk + 2*uk*smlsiz + 8*uk*nlvl + 3*smlsiz*unrhs +
                                             std::max((smlsiz + 1)*(smlsiz + 1), un*(1 + unrhs) + 2*unrhs);
                need = std::max(need, lrwork);
            }
        } else if (driver == LstsqDriver::gelss) {
            need = std::max(need, uldb*unrhs + 4*umx + uk);
            if (lq_path) { need = std::max(need, um*um + 4*um + std::max({um, unrhs, un}) + um*unrhs); }
        }
        if (need > (std::uint64_t)CBLAS_INT_MAX) { return BatchStatus::too_large; }
    }

    // outputs; S is NULL for gelsy
    real_type *ptr_S = ap_S ? (real_type *)PyArray_DATA(ap_S) : NULL;
    T *ptr_x = (T *)PyArray_DATA(ap_x);
    npy_int64 *ptr_rank = (npy_int64 *)PyArray_DATA(ap_rank);

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    // query lwork; gelsd also reports the sizes of rwork (complex flavors) and iwork
    T tmp = 0.0;
    real_type tmp_rwork = 0;
    CBLAS_INT liwork = 0;
    switch (driver) {
        case LstsqDriver::gelsd:
            if constexpr (!is_complex_v<T>) {
                lapack::gelsd(intm, intn, int_nrhs, NULL, lda, NULL, ldb, NULL, r_rcond, &rank, &tmp, -1, &liwork, &info);
            } else {
                lapack::gelsd(intm, intn, int_nrhs, NULL, lda, NULL, ldb, NULL, r_rcond, &rank, &tmp, -1, &tmp_rwork, &liwork, &info);
            }
            break;
        case LstsqDriver::gelsy:
            if constexpr (!is_complex_v<T>) {
                lapack::gelsy(intm, intn, int_nrhs, NULL, lda, NULL, ldb, NULL, r_rcond, &rank, &tmp, -1, &info);
            } else {
                lapack::gelsy(intm, intn, int_nrhs, NULL, lda, NULL, ldb, NULL, r_rcond, &rank, &tmp, -1, NULL, &info);
            }
            break;
        case LstsqDriver::gelss:
            if constexpr (!is_complex_v<T>) {
                lapack::gelss(intm, intn, int_nrhs, NULL, lda, NULL, ldb, NULL, r_rcond, &rank, &tmp, -1, &info);
            } else {
                lapack::gelss(intm, intn, int_nrhs, NULL, lda, NULL, ldb, NULL, r_rcond, &rank, &tmp, -1, NULL, &info);
            }
            break;
    }
    if (info != 0) { return BatchStatus::lapack_error; }

    // _calc_lwork returns -1 when the queried size does not fit CBLAS_INT
    CBLAS_INT lwork = _calc_lwork(tmp);
    if (lwork < 0) { return BatchStatus::too_large; }

    // real work for the complex flavors: gelsd reports it, gelsy needs 2n and gelss 5 min(m, n)
    std::size_t rwork_size = 0;
    if constexpr (is_complex_v<T>) {
        if (driver == LstsqDriver::gelsd) {
            CBLAS_INT lrwork = _calc_lwork(tmp_rwork);
            if (lrwork < 0) { return BatchStatus::too_large; }
            rwork_size = lrwork;
        } else {
            rwork_size = (driver == LstsqDriver::gelsy) ? 2 * (std::size_t)n : 5 * (std::size_t)min_mn;
        }
    }
    // integer work: iwork for gelsd, jpvt for gelsy
    std::size_t iwork_size = 0;
    if (driver == LstsqDriver::gelsd) { iwork_size = liwork; }
    if (driver == LstsqDriver::gelsy) { iwork_size = n; }

    /*
     * Allocate memory and chop the buffer into parts
     *
     *    lwork     data_a_size   data_b_size
     * |----------|-------------|-------------|
     * ^          ^             ^
     * work       data_a        data_b
     *
     * data_a is m*n unless overwrite_a, data_b is ldb*nrhs unless overwrite_b (then both = 0)
     */
    // element counts are products of dimensions of existing arrays, so they fit size_t
    const std::size_t data_a_size = overwrite_a ? 0 : (std::size_t)m * n;
    const std::size_t data_b_size = overwrite_b ? 0 : (std::size_t)ldb * nrhs;

    auto buf = alloc_workspace<T>(lwork + data_a_size + data_b_size);
    auto rwork = alloc_workspace<real_type>(rwork_size);
    auto iwork = alloc_workspace<CBLAS_INT>(iwork_size);
    if (!buf || !rwork || !iwork) { return BatchStatus::no_memory; }

    // partition the workspace; overwrite_x inputs are 2-D F-contiguous and already in place
    T *cursor = buf.get();
    T *work = cursor; cursor += lwork;
    T *data_a = Am_data;
    if (!overwrite_a) { data_a = cursor; cursor += data_a_size; }
    T *data_b = bm_data;
    if (!overwrite_b) { data_b = cursor; }

    // --------------------------------------------------------------------
    // Main loop to traverse the slices
    // A failing slice is recorded in vec_status and ends the loop; the call itself still returns ok.
    // --------------------------------------------------------------------
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        init_status(slice_status, idx, St::GENERAL);

        if (!overwrite_a) {
            T *slice_ptr = compute_slice_ptr(idx, Am_data, ndim, shape, strides);
            copy_strided_to_f(data_a, slice_ptr, m, n, strides[ndim-2], strides[ndim-1], m);
        }
        if (!overwrite_b) {
            // the solvers need ldb = max(1, m, n): b fills the first m rows, x the first n
            T *slice_ptr_b = compute_slice_ptr(idx, bm_data, ndim, shape_b, strides_b);
            copy_strided_to_f(data_b, slice_ptr_b, m, nrhs, strides_b[ndim-2], strides_b[ndim-1], ldb);
        }

        switch (driver) {
            case LstsqDriver::gelsd:
                if constexpr (!is_complex_v<T>) {
                    lapack::gelsd(intm, intn, int_nrhs, data_a, lda, data_b, ldb, ptr_S, r_rcond, &rank, work, lwork, iwork.get(), &info);
                } else {
                    lapack::gelsd(intm, intn, int_nrhs, data_a, lda, data_b, ldb, ptr_S, r_rcond, &rank, work, lwork, rwork.get(), iwork.get(), &info);
                }
                break;
            case LstsqDriver::gelsy:
                std::fill_n(iwork.get(), n, CBLAS_INT(0));   // jpvt = 0: every column is free
                if constexpr (!is_complex_v<T>) {
                    lapack::gelsy(intm, intn, int_nrhs, data_a, lda, data_b, ldb, iwork.get(), r_rcond, &rank, work, lwork, &info);
                } else {
                    lapack::gelsy(intm, intn, int_nrhs, data_a, lda, data_b, ldb, iwork.get(), r_rcond, &rank, work, lwork, rwork.get(), &info);
                }
                break;
            case LstsqDriver::gelss:
                if constexpr (!is_complex_v<T>) {
                    lapack::gelss(intm, intn, int_nrhs, data_a, lda, data_b, ldb, ptr_S, r_rcond, &rank, work, lwork, &info);
                } else {
                    lapack::gelss(intm, intn, int_nrhs, data_a, lda, data_b, ldb, ptr_S, r_rcond, &rank, work, lwork, rwork.get(), &info);
                }
                break;
        }

        if (info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            return BatchStatus::ok;
        }

        // the F-ordered solution goes to the C-ordered output unless it was computed in place in b
        if (!overwrite_b) { copy_slice_F_to_C(ptr_x, data_b, ldb, nrhs, ldb); }
        *ptr_rank = (npy_int64)rank;

        // S is filled in place; advance the C-contiguous outputs
        if (ptr_S) { ptr_S += min_mn; }
        ptr_x += (npy_intp)ldb * nrhs;
        ptr_rank += 1;
    }

    return BatchStatus::ok;
}

} // namespace sp_linalg
