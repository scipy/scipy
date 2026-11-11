/*
 * Templated loops for `linalg.eig` and `linalg.eigh`
 */
#pragma once

namespace sp_linalg {

/**
 * @brief Expand the eigenvectors of a real ?geev/?ggev call into complex ones.
 *
 * Real LAPACK returns a real column for a real eigenvalue and, for a complex conjugate pair
 * (j, j + 1), the columns Re(u) and Im(u) with u(j) = Re + i Im and u(j + 1) = Re - i Im.
 *
 * @param[out] dst  C-ordered n-by-n complex eigenvectors.
 * @param[in]  v    F-ordered real eigenvectors with leading dimension @p ldv.
 * @param[in]  ldv  Leading dimension of @p v.
 * @param[in]  n    Matrix dimension.
 * @param[in]  wi   Imaginary parts of the eigenvalues; zero marks a real eigenvalue.
 */
template<typename real_type>
void
transform_eigvecs(complex_of_t<real_type> *dst, const real_type *v, npy_intp ldv, npy_intp n, const real_type *wi) {
    using cmplx_type = complex_of_t<real_type>;
    for (npy_intp j = 0; j < n; ) {
        if (wi[j] == real_type(0)) {
            for (npy_intp i = 0; i < n; i++) {
                dst[i*n + j] = cmplx_type(v[i + j*ldv], real_type(0));
            }
            j += 1;
        } else {
            for (npy_intp i = 0; i < n; i++) {
                real_type re = v[i + j*ldv], im = v[i + (j+1)*ldv];
                dst[i*n + j] = cmplx_type(re, im);
                dst[i*n + j + 1] = cmplx_type(re, -im);
            }
            j += 2;
        }
    }
}


/**
 * @brief Batched standard or generalized eigenproblem, `linalg.eig`.
 *
 * ?geev when @p ap_Bm is NULL, else ?ggev with alpha in @p ap_w and beta in @p ap_beta.
 *
 * @param[in]  ap_Am        Input array of shape (..., n, n), any strides.
 * @param[in]  ap_Bm        Matrices of the shape of @p ap_Am for the generalized problem; NULL
 *                          otherwise.
 * @param[out] ap_w         Complex C-contiguous (..., n) eigenvalues, or alpha.
 * @param[out] ap_beta      C-contiguous (..., n) beta of the generalized problem; NULL otherwise.
 * @param[out] ap_vl        Complex C-contiguous (..., n, n) left eigenvectors; NULL if not wanted.
 * @param[out] ap_vr        Complex C-contiguous (..., n, n) right eigenvectors; NULL if not wanted.
 * @param[in]  overwrite_a  If true, @p ap_Am is 2-D and F-contiguous and is used in place.
 * @param[in]  overwrite_b  If true, @p ap_Bm is 2-D and F-contiguous and is used in place.
 * @param[out] vec_status   Receives the status of the first failing slice, if any.
 * @return ok; too_large if n or the workspace does not fit CBLAS_INT; lapack_error if the
 *         workspace query fails.
 * @throws std::bad_alloc if a workspace cannot be allocated.
 */
template<typename T>
BatchStatus
_eig(PyArrayObject *ap_Am, PyArrayObject *ap_Bm, PyArrayObject *ap_w, PyArrayObject *ap_beta,
     PyArrayObject *ap_vl, PyArrayObject *ap_vr, int overwrite_a, int overwrite_b, SliceStatusVec& vec_status)
{
    using real_type = real_of_t<T>;   // f32 if T==c64 etc
    using cmplx_type = complex_of_t<T>;
    SliceStatus slice_status;

    // --------------------------------------------------------------------
    // Input array attributes
    // --------------------------------------------------------------------
    char *Am_data = PyArray_BYTES(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    npy_intp *strides = PyArray_STRIDES(ap_Am);
    npy_intp n = shape[ndim - 1];

    const bool generalized = (ap_Bm != NULL);
    char *Bm_data = generalized ? PyArray_BYTES(ap_Bm) : NULL;
    npy_intp *strides_B = generalized ? PyArray_STRIDES(ap_Bm) : NULL;

    npy_intp outer_size = 1;
    for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i]; }

    // outputs are complex and C-contiguous
    cmplx_type *ptr_W = (cmplx_type *)PyArray_DATA(ap_w);
    T *ptr_beta = generalized ? (T *)PyArray_DATA(ap_beta) : NULL;
    cmplx_type *ptr_vl = ap_vl ? (cmplx_type *)PyArray_DATA(ap_vl) : NULL;
    cmplx_type *ptr_vr = ap_vr ? (cmplx_type *)PyArray_DATA(ap_vr) : NULL;
    const char jobvl = ptr_vl ? 'V' : 'N', jobvr = ptr_vr ? 'V' : 'N';

    CBLAS_INT intn, info = 0;
    if (!checked_narrow<CBLAS_INT>(n, &intn)) { return BatchStatus::too_large; }

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    // the complex flavors read rwork already during the workspace query: 2n for geev, 8n for ggev
    std::uint64_t rwork_size = 0;
    if constexpr (is_complex_v<T>) { rwork_size = (generalized ? 8 : 2) * (std::uint64_t)n; }
    auto rwork = alloc_workspace<real_type>(rwork_size);

    // query lwork
    T tmp = 0.0;
    if (generalized) {
        if constexpr (!is_complex_v<T>) {
            lapack::ggev(jobvl, jobvr, intn, NULL, intn, NULL, intn, NULL, NULL, NULL, NULL, intn, NULL, intn, &tmp, -1, &info);
        } else {
            lapack::ggev(jobvl, jobvr, intn, NULL, intn, NULL, intn, NULL, NULL, NULL, intn, NULL, intn, &tmp, -1, rwork.get(), &info);
        }
    } else {
        if constexpr (!is_complex_v<T>) {
            lapack::geev(jobvl, jobvr, intn, NULL, intn, NULL, NULL, NULL, intn, NULL, intn, &tmp, -1, &info);
        } else {
            lapack::geev(jobvl, jobvr, intn, NULL, intn, NULL, NULL, intn, NULL, intn, &tmp, -1, rwork.get(), &info);
        }
    }
    if (info != 0) { return BatchStatus::lapack_error; }

    // _calc_lwork returns -1 when the queried lwork does not fit CBLAS_INT
    CBLAS_INT lwork = _calc_lwork(tmp);
    if (lwork < 0) { return BatchStatus::too_large; }

    /*
     * Allocate memory and chop the buffer into parts
     *
     *    lwork      n      n/0     n/0     n*n/0    n*n/0    n*n/0    n*n/0
     * |----------|------|-------|-------|--------|--------|--------|--------|
     * ^          ^      ^       ^       ^        ^        ^        ^
     * work       wr     wi      beta    data_A   data_B   buf_vl   buf_vr
     *
     * - wr holds the eigenvalues (alpha for the generalized problem); the real flavors return
     *   the imaginary parts separately in wi
     * - beta is written by the generalized problem only
     * - data_A and data_B are the F-ordered copies of the slices unless they are used in place
     * - buf_vl and buf_vr hold the eigenvectors when they are requested
     */
    // each count is bounded by an existing array or by CBLAS_INT, so the 64-bit sums are exact
    const std::uint64_t nn = (std::uint64_t)n * n;
    const std::uint64_t wi_size = is_complex_v<T> ? 0 : (std::uint64_t)n;
    const std::uint64_t beta_size = generalized ? (std::uint64_t)n : 0;
    const std::uint64_t a_size = overwrite_a ? 0 : nn;
    const std::uint64_t b_size = (generalized && !overwrite_b) ? nn : 0;
    const std::uint64_t vl_size = ptr_vl ? nn : 0, vr_size = ptr_vr ? nn : 0;

    auto buf = alloc_workspace<T>(lwork + n + wi_size + beta_size + a_size + b_size + vl_size + vr_size);

    // partition the workspace; overwrite_x inputs are 2-D F-contiguous and already in place
    T *cursor = buf.get();
    T *work = cursor; cursor += lwork;
    T *wr = cursor; cursor += n;
    T *wi = NULL;
    if (wi_size) { wi = cursor; cursor += wi_size; }
    T *beta = NULL;
    if (beta_size) { beta = cursor; cursor += beta_size; }
    T *data_A = static_cast<T *>(PyArray_DATA(ap_Am));
    if (!overwrite_a) { data_A = cursor; cursor += a_size; }
    T *data_B = generalized ? static_cast<T *>(PyArray_DATA(ap_Bm)) : NULL;
    if (b_size) { data_B = cursor; cursor += b_size; }
    T *buf_vl = NULL;
    if (vl_size) { buf_vl = cursor; cursor += vl_size; }
    T *buf_vr = NULL;
    if (vr_size) { buf_vr = cursor; }

    // --------------------------------------------------------------------
    // Main loop to traverse the slices
    // A failing slice is recorded in vec_status and ends the loop; the call itself still returns ok.
    // --------------------------------------------------------------------
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        init_status(slice_status, idx, St::GENERAL);

        if (!overwrite_a) {
            T *slice_ptr = compute_slice_ptr<T>(idx, Am_data, ndim, shape, strides);
            copy_strided_to_f(data_A, slice_ptr, n, n, strides[ndim-2], strides[ndim-1], n);
        }
        if (b_size) {
            T *slice_ptr = compute_slice_ptr<T>(idx, Bm_data, ndim, shape, strides_B);
            copy_strided_to_f(data_B, slice_ptr, n, n, strides_B[ndim-2], strides_B[ndim-1], n);
        }

        if (generalized) {
            if constexpr (!is_complex_v<T>) {
                lapack::ggev(jobvl, jobvr, intn, data_A, intn, data_B, intn, wr, wi, beta, buf_vl, intn, buf_vr, intn, work, lwork, &info);
            } else {
                lapack::ggev(jobvl, jobvr, intn, data_A, intn, data_B, intn, wr, beta, buf_vl, intn, buf_vr, intn, work, lwork, rwork.get(), &info);
            }
        } else {
            if constexpr (!is_complex_v<T>) {
                lapack::geev(jobvl, jobvr, intn, data_A, intn, wr, wi, buf_vl, intn, buf_vr, intn, work, lwork, &info);
            } else {
                lapack::geev(jobvl, jobvr, intn, data_A, intn, wr, buf_vl, intn, buf_vr, intn, work, lwork, rwork.get(), &info);
            }
        }

        if (info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            return BatchStatus::ok;
        }

        // eigenvalues, beta and the transposed eigenvectors go to the C-contiguous outputs
        if constexpr (is_complex_v<T>) {
            std::copy_n(wr, n, ptr_W);
            if (buf_vl) { copy_slice_F_to_C(ptr_vl, buf_vl, n, n, n); ptr_vl += nn; }
            if (buf_vr) { copy_slice_F_to_C(ptr_vr, buf_vr, n, n, n); ptr_vr += nn; }
        } else {
            // the real flavors return (wr, wi) and real eigenvector columns in the ?geev convention
            for (npy_intp i = 0; i < n; i++) { ptr_W[i] = cmplx_type(wr[i], wi[i]); }
            if (buf_vl) { transform_eigvecs(ptr_vl, buf_vl, n, n, wi); ptr_vl += nn; }
            if (buf_vr) { transform_eigvecs(ptr_vr, buf_vr, n, n, wi); ptr_vr += nn; }
        }
        ptr_W += n;
        if (beta) { std::copy_n(beta, n, ptr_beta); ptr_beta += n; }
    }

    return BatchStatus::ok;
}


/**
 * @brief Batched symmetric or Hermitian eigenproblem, `linalg.eigh`.
 *
 * @param[in]  ap_Am        Input array of shape (..., n, n), any strides; the @p uplo triangle is
 *                          read.
 * @param[in]  ap_Bm        Positive definite matrices of the shape of @p ap_Am for the generalized
 *                          drivers gv, gvd and gvx; NULL for the others.
 * @param[out] ap_w         Real C-contiguous (..., m) eigenvalues, m = iu - il + 1 for range 'I'
 *                          and n otherwise.
 * @param[out] ap_Z         C-contiguous (..., n, m) eigenvectors for @p jobz 'V', NULL for 'N'.
 *                          With @p overwrite_a this is @p ap_Am itself and receives them in F order.
 * @param[out] m_found      Receives the number of eigenvalues found, less than m only for range 'V'.
 * @param[in]  overwrite_a  If true, @p ap_Am is 2-D and F-contiguous and is used in place.
 * @param[in]  overwrite_b  If true, @p ap_Bm is 2-D and F-contiguous and is used in place.
 * @param[in]  itype        Generalized problem type 1, 2 or 3.
 * @param[in]  jobz         'N' for eigenvalues only, 'V' for eigenvectors too.
 * @param[in]  range        'A' all eigenvalues, 'V' those in (vl, vu], 'I' the il-th to iu-th.
 * @param[in]  uplo         'L' or 'U', the triangle that holds the data.
 * @param[in]  vl           Lower bound for range 'V'.
 * @param[in]  vu           Upper bound for range 'V'.
 * @param[in]  il           First 0-based index for range 'I'.
 * @param[in]  iu           Last 0-based index for range 'I'.
 * @param[in]  driver       The LAPACK routine.
 * @param[out] vec_status   Receives the status of the first failing slice, if any.
 * @return ok; too_large if n or a workspace size does not fit CBLAS_INT; lapack_error if the
 *         workspace query fails.
 * @throws std::bad_alloc if a workspace cannot be allocated.
 */
template<typename T>
BatchStatus
_eigh(PyArrayObject *ap_Am, PyArrayObject *ap_Bm, PyArrayObject *ap_w, PyArrayObject *ap_Z, npy_intp *m_found,
      int overwrite_a, int overwrite_b, int itype, char jobz, char range, char uplo,
      f64 vl, f64 vu, int il, int iu, Eigh_driver driver, SliceStatusVec& vec_status)
{
    using real_type = real_of_t<T>;
    SliceStatus slice_status;

    // --------------------------------------------------------------------
    // Input array attributes
    // --------------------------------------------------------------------
    char *A_data = PyArray_BYTES(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp *shape = PyArray_SHAPE(ap_Am);
    npy_intp *strides = PyArray_STRIDES(ap_Am);
    npy_intp n = shape[ndim - 1];
    const npy_intp m = (range == 'I') ? iu - il + 1 : n;   // eigenvalues per slice in the outputs

    const bool generalized = (ap_Bm != NULL);
    char *B_data = generalized ? PyArray_BYTES(ap_Bm) : NULL;
    npy_intp *strides_B = generalized ? PyArray_STRIDES(ap_Bm) : NULL;

    npy_intp outer_size = 1;
    for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i]; }

    // outputs; w is (..., m) and Z is (..., n, m), both C-contiguous, unless Z is A in place
    real_type *ptr_w = (real_type *)PyArray_DATA(ap_w);
    T *ptr_Z = (jobz != 'N') ? (T *)PyArray_DATA(ap_Z) : NULL;

    CBLAS_INT intn, intm, info = 0;
    if (!checked_narrow<CBLAS_INT>(n, &intn) || !checked_narrow<CBLAS_INT>(m, &intm)) { return BatchStatus::too_large; }
    const CBLAS_INT int_itype = itype, int_il = il + 1, int_iu = iu + 1;   // LAPACK indices are 1-based
    const real_type r_vl = (real_type)vl, r_vu = (real_type)vu, abstol = 0;

    // ev, evd, gv and gvd leave the eigenvectors in A; evr, evx and gvx write them to Z
    const bool vectors_in_a = (driver == Eigh_driver::EV || driver == Eigh_driver::EVD ||
                               driver == Eigh_driver::GV || driver == Eigh_driver::GVD);
    const bool expert = (driver == Eigh_driver::EVX || driver == Eigh_driver::GVX);

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    // With eigenvectors, ?syevd/?sygvd need lwork >= 1 + 6n + 2n^2, and ?heevd/?hegvd need
    // lwork >= 2n + n^2 and lrwork >= 1 + 5n + 2n^2. LAPACK evaluates these in CBLAS_INT; when
    // they overflow, the workspace query returns a value too small to run with instead of
    // failing, so they are checked here first.
    if ((driver == Eigh_driver::EVD || driver == Eigh_driver::GVD) && jobz == 'V') {
        const std::uint64_t un = n, lwmin = is_complex_v<T> ? un * un + 2 * un : 2 * un * un + 6 * un + 1;
        const std::uint64_t lrwmin = is_complex_v<T> ? 2 * un * un + 5 * un + 1 : 0;
        if (lwmin > (std::uint64_t)CBLAS_INT_MAX || lrwmin > (std::uint64_t)CBLAS_INT_MAX) { return BatchStatus::too_large; }
    }

    // query lwork; the divide-and-conquer and RRR drivers also report lrwork and liwork
    T tmp_work = 0.0;
    real_type tmp_rwork = 0.0;
    CBLAS_INT tmp_iwork = 0;
    switch (driver) {
        case Eigh_driver::EV:
            if constexpr (!is_complex_v<T>) {
                lapack::syev(jobz, uplo, intn, NULL, intn, NULL, &tmp_work, -1, &info);
            } else {
                lapack::heev(jobz, uplo, intn, NULL, intn, NULL, &tmp_work, -1, &tmp_rwork, &info);
            }
            break;
        case Eigh_driver::EVD:
            if constexpr (!is_complex_v<T>) {
                lapack::syevd(jobz, uplo, intn, NULL, intn, NULL, &tmp_work, -1, &tmp_iwork, -1, &info);
            } else {
                lapack::heevd(jobz, uplo, intn, NULL, intn, NULL, &tmp_work, -1, &tmp_rwork, -1, &tmp_iwork, -1, &info);
            }
            break;
        case Eigh_driver::EVR:
            if constexpr (!is_complex_v<T>) {
                lapack::syevr(jobz, range, uplo, intn, NULL, intn, r_vl, r_vu, int_il, int_iu, abstol, &intm, NULL, NULL, intn, NULL, &tmp_work, -1, &tmp_iwork, -1, &info);
            } else {
                lapack::heevr(jobz, range, uplo, intn, NULL, intn, r_vl, r_vu, int_il, int_iu, abstol, &intm, NULL, NULL, intn, NULL, &tmp_work, -1, &tmp_rwork, -1, &tmp_iwork, -1, &info);
            }
            break;
        case Eigh_driver::EVX:
            if constexpr (!is_complex_v<T>) {
                lapack::syevx(jobz, range, uplo, intn, NULL, intn, r_vl, r_vu, int_il, int_iu, abstol, &intm, NULL, NULL, intn, &tmp_work, -1, &tmp_iwork, NULL, &info);
            } else {
                lapack::heevx(jobz, range, uplo, intn, NULL, intn, r_vl, r_vu, int_il, int_iu, abstol, &intm, NULL, NULL, intn, &tmp_work, -1, &tmp_rwork, &tmp_iwork, NULL, &info);
            }
            break;
        case Eigh_driver::GV:
            if constexpr (!is_complex_v<T>) {
                lapack::sygv(int_itype, jobz, uplo, intn, NULL, intn, NULL, intn, NULL, &tmp_work, -1, &info);
            } else {
                lapack::hegv(int_itype, jobz, uplo, intn, NULL, intn, NULL, intn, NULL, &tmp_work, -1, &tmp_rwork, &info);
            }
            break;
        case Eigh_driver::GVD:
            if constexpr (!is_complex_v<T>) {
                lapack::sygvd(int_itype, jobz, uplo, intn, NULL, intn, NULL, intn, NULL, &tmp_work, -1, &tmp_iwork, -1, &info);
            } else {
                lapack::hegvd(int_itype, jobz, uplo, intn, NULL, intn, NULL, intn, NULL, &tmp_work, -1, &tmp_rwork, -1, &tmp_iwork, -1, &info);
            }
            break;
        case Eigh_driver::GVX:
            if constexpr (!is_complex_v<T>) {
                lapack::sygvx(int_itype, jobz, range, uplo, intn, NULL, intn, NULL, intn, r_vl, r_vu, int_il, int_iu, abstol, &intm, NULL, NULL, intn, &tmp_work, -1, &tmp_iwork, NULL, &info);
            } else {
                lapack::hegvx(int_itype, jobz, range, uplo, intn, NULL, intn, NULL, intn, r_vl, r_vu, int_il, int_iu, abstol, &intm, NULL, NULL, intn, &tmp_work, -1, &tmp_rwork, &tmp_iwork, NULL, &info);
            }
            break;
    }
    if (info != 0) { return BatchStatus::lapack_error; }

    // _calc_lwork returns -1 when the queried size does not fit CBLAS_INT
    CBLAS_INT lwork = _calc_lwork(tmp_work);
    if (lwork < 0) { return BatchStatus::too_large; }

    // real work, complex flavors only: fixed for ev/gv (3n - 2) and evx/gvx (7n), queried otherwise
    std::uint64_t lrwork = 0;
    if constexpr (is_complex_v<T>) {
        if (driver == Eigh_driver::EV || driver == Eigh_driver::GV) {
            lrwork = 3 * (std::uint64_t)n - 2;
        } else if (expert) {
            lrwork = 7 * (std::uint64_t)n;
        } else {
            CBLAS_INT queried = _calc_lwork(tmp_rwork);
            if (queried < 0) { return BatchStatus::too_large; }
            lrwork = queried;
        }
    }
    // integer work: none for ev/gv, 5n plus an n-long ifail for evx/gvx, queried otherwise; evr also
    // needs 2m for isuppz
    std::uint64_t liwork = 0, isuppz_size = 0, ifail_size = 0;
    if (driver == Eigh_driver::EVD || driver == Eigh_driver::GVD || driver == Eigh_driver::EVR) { liwork = tmp_iwork; }
    if (expert) { liwork = 5 * (std::uint64_t)n; ifail_size = n; }
    if (driver == Eigh_driver::EVR) { isuppz_size = 2 * (std::uint64_t)m; }

    CBLAS_INT int_lrwork, int_liwork;
    if (!checked_narrow<CBLAS_INT>(lrwork, &int_lrwork) || !checked_narrow<CBLAS_INT>(liwork, &int_liwork)) {
        return BatchStatus::too_large;
    }

    /*
     * Three buffers, chopped into parts:
     *
     *   T:          work (lwork) | Z (n*m, evr/evx/gvx) | A (n*n unless in place) | B (n*n, generalized unless in place)
     *   real_type:  w (n, when range 'I' wants fewer than n values: LAPACK still writes n) | rwork (lrwork)
     *   CBLAS_INT:  iwork (liwork) | isuppz (2m, evr) | ifail (n, evx/gvx)
     */
    // each count is bounded by an existing array or by CBLAS_INT, so the 64-bit sums are exact
    const std::uint64_t nn = (std::uint64_t)n * n;
    const std::uint64_t z_size = (ptr_Z && !vectors_in_a) ? (std::uint64_t)n * m : 0;
    const std::uint64_t a_size = overwrite_a ? 0 : nn;
    const std::uint64_t b_size = (generalized && !overwrite_b) ? nn : 0;
    const std::uint64_t w_size = (m < n) ? (std::uint64_t)n : 0;

    auto buf = alloc_workspace<T>(lwork + z_size + a_size + b_size);
    auto rbuf = alloc_workspace<real_type>(w_size + lrwork);
    auto ibuf = alloc_workspace<CBLAS_INT>(liwork + isuppz_size + ifail_size);

    // partition the workspaces; overwrite_x inputs are 2-D F-contiguous and already in place
    T *cursor = buf.get();
    T *work = cursor; cursor += lwork;
    T *buf_Z = NULL;
    if (z_size) { buf_Z = cursor; cursor += z_size; }
    T *buf_A = static_cast<T *>(PyArray_DATA(ap_Am));
    if (!overwrite_a) { buf_A = cursor; cursor += a_size; }
    T *buf_B = generalized ? static_cast<T *>(PyArray_DATA(ap_Bm)) : NULL;
    if (b_size) { buf_B = cursor; }

    real_type *rcursor = rbuf.get();
    real_type *buf_w = NULL;
    if (w_size) { buf_w = rcursor; rcursor += w_size; }
    real_type *rwork = rcursor;

    CBLAS_INT *icursor = ibuf.get();
    CBLAS_INT *iwork = icursor; icursor += liwork;
    CBLAS_INT *isuppz = icursor; icursor += isuppz_size;
    CBLAS_INT *ifail = icursor;

    // --------------------------------------------------------------------
    // Main loop to traverse the slices
    // A failing slice is recorded in vec_status and ends the loop; the call itself still returns ok.
    // --------------------------------------------------------------------
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        init_status(slice_status, idx, St::GENERAL);

        if (!overwrite_a) {
            T *slice_ptr = compute_slice_ptr<T>(idx, A_data, ndim, shape, strides);
            copy_strided_to_f(buf_A, slice_ptr, n, n, strides[ndim-2], strides[ndim-1], n);
        }
        if (b_size) {
            T *slice_ptr = compute_slice_ptr<T>(idx, B_data, ndim, shape, strides_B);
            copy_strided_to_f(buf_B, slice_ptr, n, n, strides_B[ndim-2], strides_B[ndim-1], n);
        }
        real_type *w = buf_w ? buf_w : ptr_w;   // LAPACK writes n eigenvalues; the output holds m

        switch (driver) {
            case Eigh_driver::EV:
                if constexpr (!is_complex_v<T>) {
                    lapack::syev(jobz, uplo, intn, buf_A, intn, w, work, lwork, &info);
                } else {
                    lapack::heev(jobz, uplo, intn, buf_A, intn, w, work, lwork, rwork, &info);
                }
                break;
            case Eigh_driver::EVD:
                if constexpr (!is_complex_v<T>) {
                    lapack::syevd(jobz, uplo, intn, buf_A, intn, w, work, lwork, iwork, int_liwork, &info);
                } else {
                    lapack::heevd(jobz, uplo, intn, buf_A, intn, w, work, lwork, rwork, int_lrwork, iwork, int_liwork, &info);
                }
                break;
            case Eigh_driver::EVR:
                if constexpr (!is_complex_v<T>) {
                    lapack::syevr(jobz, range, uplo, intn, buf_A, intn, r_vl, r_vu, int_il, int_iu, abstol, &intm, w, buf_Z, intn, isuppz, work, lwork, iwork, int_liwork, &info);
                } else {
                    lapack::heevr(jobz, range, uplo, intn, buf_A, intn, r_vl, r_vu, int_il, int_iu, abstol, &intm, w, buf_Z, intn, isuppz, work, lwork, rwork, int_lrwork, iwork, int_liwork, &info);
                }
                break;
            case Eigh_driver::EVX:
                if constexpr (!is_complex_v<T>) {
                    lapack::syevx(jobz, range, uplo, intn, buf_A, intn, r_vl, r_vu, int_il, int_iu, abstol, &intm, w, buf_Z, intn, work, lwork, iwork, ifail, &info);
                } else {
                    lapack::heevx(jobz, range, uplo, intn, buf_A, intn, r_vl, r_vu, int_il, int_iu, abstol, &intm, w, buf_Z, intn, work, lwork, rwork, iwork, ifail, &info);
                }
                break;
            case Eigh_driver::GV:
                if constexpr (!is_complex_v<T>) {
                    lapack::sygv(int_itype, jobz, uplo, intn, buf_A, intn, buf_B, intn, w, work, lwork, &info);
                } else {
                    lapack::hegv(int_itype, jobz, uplo, intn, buf_A, intn, buf_B, intn, w, work, lwork, rwork, &info);
                }
                break;
            case Eigh_driver::GVD:
                if constexpr (!is_complex_v<T>) {
                    lapack::sygvd(int_itype, jobz, uplo, intn, buf_A, intn, buf_B, intn, w, work, lwork, iwork, int_liwork, &info);
                } else {
                    lapack::hegvd(int_itype, jobz, uplo, intn, buf_A, intn, buf_B, intn, w, work, lwork, rwork, int_lrwork, iwork, int_liwork, &info);
                }
                break;
            case Eigh_driver::GVX:
                if constexpr (!is_complex_v<T>) {
                    lapack::sygvx(int_itype, jobz, range, uplo, intn, buf_A, intn, buf_B, intn, r_vl, r_vu, int_il, int_iu, abstol, &intm, w, buf_Z, intn, work, lwork, iwork, ifail, &info);
                } else {
                    lapack::hegvx(int_itype, jobz, range, uplo, intn, buf_A, intn, buf_B, intn, r_vl, r_vu, int_il, int_iu, abstol, &intm, w, buf_Z, intn, work, lwork, rwork, iwork, ifail, &info);
                }
                break;
        }

        if (info != 0) {
            slice_status.lapack_info = (Py_ssize_t)info;
            vec_status.push_back(slice_status);

            // cut it short on error in any slice
            return BatchStatus::ok;
        }

        if (buf_w) { std::copy_n(buf_w, m, ptr_w); }
        ptr_w += m;

        if (ptr_Z) {
            T *vectors = vectors_in_a ? buf_A : buf_Z;
            if (overwrite_a) {
                // Z is the F-contiguous input itself: it already holds the vectors, or receives
                // the first m columns from the Z buffer
                if (!vectors_in_a) { std::copy_n(buf_Z, z_size, ptr_Z); }
            } else {
                copy_slice_F_to_C(ptr_Z, vectors, n, m, n);
                ptr_Z += (npy_intp)n * m;
            }
        }

        // range 'V' finds an unknown number of eigenvalues; the Python layer slices with it
        *m_found = intm;
    }

    return BatchStatus::ok;
}

} // namespace sp_linalg
