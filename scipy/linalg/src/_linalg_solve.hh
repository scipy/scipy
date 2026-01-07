/*
 * Templated loops for `linalg.solve`
 */
#include "Python.h"
#include <iostream>
#include "cpython/code.h"
#include "numpy/arrayobject.h"
#include "numpy/ndarraytypes.h"
#include "numpy/npy_math.h"
#include "npy_cblas.h"
#include "_npymath.hh"
#include "_common_array_utils.hh"



// Dense array solve with getrf, gecon and getrs
template<typename T>
inline void solve_slice_general(
    CBLAS_INT N, CBLAS_INT NRHS, T *data, CBLAS_INT *ipiv, T *b_data, char trans, void *irwork, T *work,
    SliceStatus& status
) {
    using real_type = typename type_traits<T>::real_type;

    CBLAS_INT info;
    char norm = '1';
    real_type rcond;
    real_type anorm = norm1_(data, work, (npy_intp)N);

    getrf(&N, &N, data, &N, ipiv, &info);

    status.lapack_info = (Py_ssize_t)info;
    if (info == 0){
        // getrf success, check the condition number
        gecon(&norm, &N, data, &N, &anorm, &rcond, work, irwork, &info);

        status.rcond = (double)rcond;
        if (info >= 0) {
            status.is_ill_conditioned = (rcond != rcond) || (rcond < numeric_limits<real_type>::eps);

            // finally, solve
            getrs(&trans, &N, &NRHS, data, &N, ipiv, b_data, &N, &info);
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
    char uplo, char diag, CBLAS_INT N, CBLAS_INT NRHS, T *data,  T *b_data, char trans, T *work, void *irwork,
    SliceStatus& status
) {
    using real_type = typename type_traits<T>::real_type;

    CBLAS_INT info;
    char norm = '1';
    real_type rcond;

    trtrs(&uplo, &trans, &diag, &N, &NRHS, data, &N, b_data, &N, &info);

    status.lapack_info = (Py_ssize_t)info;
    status.is_singular  = (info > 0);
    if(info >= 0) {
        trcon(&norm, &uplo, &diag, &N, data, &N, &rcond, work, irwork, &info);
        if (info >= 0) {
            status.is_ill_conditioned = (rcond != rcond) || (rcond < numeric_limits<real_type>::eps);
            status.rcond = (double)rcond;
        }
    }
}


// Cholesky solve with potrf, pocon and potrs
template<typename T>
inline void solve_slice_cholesky(
    char uplo, CBLAS_INT N, CBLAS_INT NRHS, T *data, T *b_data, T* work, void *irwork,
    SliceStatus& status
) {
    using real_type = typename type_traits<T>::real_type;

    CBLAS_INT info;
    real_type rcond;
    real_type anorm = norm1_sym_herm(uplo, data, work, (npy_intp)N);

    potrf(&uplo, &N, data, &N, &info);

    status.lapack_info = (Py_ssize_t)info;
    if (info == 0) {
        // potrf success
        pocon(&uplo, &N, data, &N, &anorm, &rcond, work, irwork, &info);

        if (info >= 0) {
            status.rcond = (double)rcond;
            status.is_ill_conditioned = (rcond != rcond) || (rcond < numeric_limits<real_type>::eps);

            // finally, solve
            potrs(&uplo, &N, &NRHS, data, &N, b_data, &N, &info);
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
    char uplo, CBLAS_INT N, CBLAS_INT NRHS, T *data, T *b_data, CBLAS_INT *ipiv, T *work, void *irwork, CBLAS_INT lwork,
    bool is_symm_not_herm,
    SliceStatus& status
) {
    using real_type = typename type_traits<T>::real_type;

    CBLAS_INT info;
    real_type rcond;
    real_type anorm = norm1_sym_herm(uplo, data, work, (npy_intp)N);

    if(is_symm_not_herm) {
        sytrf(&uplo, &N, data, &N, ipiv, work, &lwork, &info);
    } else {
        hetrf(&uplo, &N, data, &N, ipiv, work, &lwork, &info);
    }

    status.lapack_info = (Py_ssize_t)info;
    if (info == 0) {
        // {sy,he}trf success
        if (is_symm_not_herm) {
            sycon(&uplo, &N, data, &N, ipiv, &anorm, &rcond, work, irwork, &info);
        } else {
            hecon(&uplo, &N, data, &N, ipiv, &anorm, &rcond, work, irwork, &info);
        }

        if (info >= 0) {
            status.rcond = (double)rcond;
            status.is_ill_conditioned = (rcond != rcond) || (rcond < numeric_limits<real_type>::eps);

            // finally, solve
            if (is_symm_not_herm) {
                sytrs(&uplo, &N, &NRHS, data, &N, ipiv, b_data, &N, &info);
            } else {
                hetrs(&uplo, &N, &NRHS, data, &N, ipiv, b_data, &N, &info);
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
    T *work, T *work2, void *irwork,
    SliceStatus& status
) {
    using real_type = typename type_traits<T>::real_type;
    // work is 4*n, is for dl, d, du, du2
    // work2 is 2*n, is for trcon's work array

    // copy the diagonals of `data`, use `work` for storage
    T *du2 = work;
    T *du = &work[N];
    T *d = &work[2*N];
    T *dl = &work[3*N];
    to_tridiag(data, N, du, d, dl);

    CBLAS_INT info;
    CBLAS_INT *iwork = (CBLAS_INT *)irwork;
    char norm = '1';
    real_type rcond;
    real_type anorm = norm1_tridiag(dl, d, du, work2, (npy_intp)N);

    gttrf(&N, dl, d, du, du2, ipiv, &info);

    status.lapack_info = (Py_ssize_t)info;
    if (info == 0){
        // gttrf success, check the condition number
        gtcon(&norm, &N, dl, d, du, du2, ipiv, &anorm, &rcond, work2, iwork, &info);

        status.rcond = (double)rcond;
        if (info >= 0) {
            status.is_ill_conditioned = (rcond != rcond) || (rcond < numeric_limits<real_type>::eps);

            // finally, solve
            gttrs(&trans, &N, &NRHS, dl, d, du, du2, ipiv, b_data, &N, &info);
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
    char trans, CBLAS_INT N, CBLAS_INT NRHS, T *ab, CBLAS_INT *ipiv, T *b_data, T *work2, void *irwork,
    CBLAS_INT kl, CBLAS_INT ku, SliceStatus &status
) {
    using real_type = typename type_traits<T>::real_type;

    CBLAS_INT ldab = 2 * kl + ku + 1;

    CBLAS_INT info;
    char norm = '1';
    real_type rcond;
    real_type anorm = norm1_banded(ab, kl, ku, work2, N);

    gbtrf(&N, &N, &kl, &ku, ab, &ldab, ipiv, &info);
    status.lapack_info = (Py_ssize_t)info;
    if (info == 0) {
        // gbtrf success, check condition number
        gbcon(&norm, &N, &kl, &ku, ab, &ldab, ipiv, &anorm, &rcond, work2, irwork, &info);

        status.rcond = (double)rcond;
        if (info >= 0) {
            status.is_ill_conditioned = (rcond != rcond) || (rcond < numeric_limits<real_type>::eps);

            // finally, solve
            gbtrs(&trans, &N, &kl, &ku, &NRHS, ab, &ldab, ipiv, b_data, &N, &info);
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
    using real_type = typename type_traits<T>::real_type;
    using value_type = typename type_traits<T>::value_type;
    value_type *pdata = reinterpret_cast<value_type *>(data);
    value_type *p_bdata = reinterpret_cast<value_type *>(b_data);

    value_type zero(0.), one(1.);
    real_type maxa(0.), maxinva(0.);

    for (CBLAS_INT j=0; j<N; j++) {
        value_type ajj = pdata[j*N + j];

        status.is_singular = (ajj == zero);
        if (status.is_singular) {
            status.lapack_info = j;
            return;
        }

        value_type inv_ajj = one / ajj;
        for (CBLAS_INT i=0; i<NRHS; i++) {
            p_bdata[j + i*N] *= inv_ajj;
        }

        // condition number
        real_type absa = std::abs(ajj), absinva = std::abs(inv_ajj);

        if(absa > maxa) {maxa = absa;}
        if(absinva > maxinva) {maxinva = absinva;}
    }
    status.is_ill_conditioned = maxa * maxinva > 1./ numeric_limits<real_type>::eps;
    status.rcond = maxa * maxinva;
}


// Specific implementation for the case of `solve(..., assume_a=banded)` to avoid
// cluttering the main loop in `_solve` by branching to this function early.
template<typename T>
int
_solve_assume_banded(PyArrayObject *ap_Am, PyArrayObject *ap_b, T *ret_data, char trans, int overwrite_a, SliceStatus slice_status, SliceStatusVec &vec_status)
{
    using real_type = typename type_traits<T>::real_type;

    CBLAS_INT info;
    npy_intp *ks = NULL; // For storage of the bandwidths
    T* buffer = NULL;

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

    CBLAS_INT intn = (CBLAS_INT)n, int_nrhs = (CBLAS_INT)nrhs;

    // General allocations
    CBLAS_INT *ipiv = (CBLAS_INT *)malloc(intn * sizeof(CBLAS_INT));
    if (ipiv == NULL) {
        info = -102;
        return int(info);
    }

    void *irwork;
    if (type_traits<T>::is_complex) {
        irwork = malloc(intn * sizeof(real_type));
    } else {
        irwork = malloc(intn * sizeof(CBLAS_INT));
    }

    if (irwork == NULL) {
        free(ipiv);
        info = -102;
        return int(info);
    }

    // Bandwidth detection per slice
    npy_intp kl_max = 0;
    npy_intp ku_max = 0;
    ks = (npy_intp *)malloc(2 * outer_size * sizeof(npy_intp));

    if (ks == NULL) {
        free(ipiv);
        free(irwork);
        free(ks);
        info = -102;
        return (int)info;
    }

    npy_intp *kls = &ks[0]; // Lower bandwidths
    npy_intp *kus = &ks[outer_size]; // Upper bandwidths
    detect_bandwidths(Am_data, ndim, outer_size, shape, strides, kls, kus, &kl_max, &ku_max);

    buffer = (T *)malloc((n * nrhs + 3 * n + (2 * kl_max + ku_max + 1) * n) * sizeof(T));

    if (buffer == NULL) {
        free(ipiv);
        free(irwork);
        free(ks);
        free(buffer);
        info = -102;
        return int(info);
    }

    // Chop up buffer
    T* b_data = &buffer[0];
    T* work2 = &buffer[n * nrhs]; // for `gbcon` call
    T *ab = &buffer[n * nrhs + 3 * n];

    // Main loop traversal, taken from `_solve`
    for (npy_intp idx = 0; idx < outer_size; idx++) {

        npy_intp offset = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset += (temp_idx % shape[i]) * strides[i];
            temp_idx /= shape[i];
        }
        T* slice_ptr = (T *)(Am_data + (offset / sizeof(T)));

        // Directly take into banded storage
        npy_intp ldab = 2 * kls[idx] + kus[idx] + 1;
        to_banded(slice_ptr, n, kls[idx], kus[idx], ldab, ab, strides[ndim-2], strides[ndim-1]);

        // XXX: dedupe (cfr. _solve)
        offset = 0;
        temp_idx = idx;
        for (int i = ndim_b - 3; i >= 0; i--) {
            offset += (temp_idx & shape_b[i]) * strides_b[i];
            temp_idx /= shape_b[i];
        }
        T* slice_ptr_b = (T *)(bm_data + (offset / sizeof(T)));
        copy_slice_F(b_data, slice_ptr_b, n, nrhs, strides_b[ndim-2], strides_b[ndim-1]);

        // structure is known to be banded
        init_status(slice_status, idx, St::BANDED);
        solve_slice_banded(trans, intn, int_nrhs, ab, ipiv, b_data, work2, irwork, kls[idx], kus[idx], slice_status);

        if ((slice_status.lapack_info < 0) || (slice_status.is_singular)) {
            vec_status.push_back(slice_status);
            goto free_exit_banded;
        } else if (slice_status.is_ill_conditioned) {
            vec_status.push_back(slice_status);
        }

        // Put result in C-order in return buffer
        copy_slice_F_to_C(&ret_data[idx * n * nrhs], b_data, n, nrhs);
    }

free_exit_banded:
    free(ipiv);
    free(irwork);
    free(ks);
    free(buffer);

    return 1;
}

template<typename T>
int
_solve(PyArrayObject* ap_Am, PyArrayObject *ap_b, T* ret_data, St structure, int lower, int transposed, int overwrite_a, SliceStatusVec& vec_status)
{
    using real_type = typename type_traits<T>::real_type; // float if T==npy_cfloat etc

    char trans = transposed ? 'T' : 'N';
    npy_intp lower_band = 0, upper_band = 0;
    bool is_symm_or_herm = false, is_symm_not_herm = false;
    char uplo = lower ? 'L' : 'U';
    St slice_structure = St::NONE;
    bool posdef_fallback = true;
    SliceStatus slice_status;

    // branch early for `assume_a = banded`
    if (structure == St::BANDED) {
        return _solve_assume_banded(ap_Am, ap_b, ret_data, trans, overwrite_a, slice_status, vec_status);
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
    CBLAS_INT intn = (CBLAS_INT)n, int_nrhs = (CBLAS_INT)nrhs, lwork=-1, info;

    T tmp = numeric_limits<T>::zero;
    sytrf(&uplo, &intn, NULL, &intn, NULL, &tmp, &lwork, &info);
    if (info != 0) { info = -100; return (int)info; }

    lwork = _calc_lwork(tmp);
    if ((lwork < 0) ||
        (n > std::numeric_limits<CBLAS_INT>::max() / 4)
    ) {
        // Too large lwork required - Computation cannot be performed.
        // if CBLAS_INT is 32-bit, need ILP64; if already using ILP64, we're out of luck.
        return -99;
    }

    // gecon needs lwork of at least 4*n
    lwork = (4*n > lwork ? 4*n : lwork);

    T* buffer = (T *)malloc((2*n*n + n*nrhs + 2*n + lwork)*sizeof(T));
    if (NULL == buffer) { info = -101; return (int)info; }

    // Chop the buffer into parts
    T* data = &buffer[0];
    T* scratch = &buffer[n*n];

    T *data_b = &buffer[2*n*n];
    T *work2 = &buffer[2*n*n + n*nrhs]; // 2*n for is for tridiag's trcon; XXX malloc it only if needed?
    T* work = &buffer[2*n*n + n*nrhs + 2*n];

    CBLAS_INT* ipiv = (CBLAS_INT *)malloc(n*sizeof(CBLAS_INT));
    if (ipiv == NULL) {
        free(buffer);
        info = -102;
        return (int)info;
    }

    // {ge,po,tr}con need rwork or iwork
    void *irwork;
    if (type_traits<T>::is_complex) {
        irwork = malloc(3*n*sizeof(real_type));   // {po,tr}con need at least 3*n
    } else {
        irwork = malloc(n*sizeof(CBLAS_INT));
    }
    if (irwork == NULL) {
        free(buffer);
        free(ipiv);
        info = -102;
        return (int)info;
    }

    // normalize the structure detection inputs
    if (structure == St::POS_DEF) {
        posdef_fallback = false;
    }
    else if (structure == St::SYM) {
        is_symm_not_herm = true;
    }
    else if (structure == St::HER) {
        is_symm_not_herm = false;
    }
    if (structure == St::LOWER_TRIANGULAR) {
        uplo = 'L';
    }
    else if (structure == St::UPPER_TRIANGULAR) {
        uplo = 'U';
    }

    // Main loop to traverse the slices
    for (npy_intp idx = 0; idx < outer_size; idx++) {

        npy_intp offset = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset += (temp_idx % shape[i]) * strides[i];
            temp_idx /= shape[i];
        }
        T* slice_ptr = (T *)(Am_data + (offset/sizeof(T)));
        copy_slice(scratch, slice_ptr, n, n, strides[ndim-2], strides[ndim-1]); // XXX: make it in one go
        swap_cf(scratch, data, n, n, n);

        // copy the r.h.s, too; XXX: dedupe
        offset = 0;
        temp_idx = idx;
        for (int i = ndim_b - 3; i >= 0; i--) {
            offset += (temp_idx % shape_b[i]) * strides_b[i];
            temp_idx /= shape_b[i];
        }
        T *slice_ptr_b = (T *)(bm_data + (offset/sizeof(T)));
        copy_slice_F(data_b, slice_ptr_b, n, nrhs, strides_b[ndim-2], strides_b[ndim-1]);

        // detect the structure if not given
        slice_structure = structure;

        if (slice_structure == St::NONE) {
            // Get the bandwidth of the slice
            bandwidth(data, n, n, &lower_band, &upper_band);

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
                std::tie(is_symm_or_herm, is_symm_not_herm) = is_sym_herm(data, n);
                if (is_symm_or_herm) {
                    slice_structure = St::POS_DEF;
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

                if ((slice_status.lapack_info < 0) || (slice_status.is_singular)) {
                    vec_status.push_back(slice_status);
                    goto free_exit;
                }
                else if (slice_status.is_ill_conditioned) {
                    vec_status.push_back(slice_status);
                }

                break;
            }
            case St::TRIDIAGONAL:
            {
                solve_slice_tridiag(intn, int_nrhs, data, ipiv, trans, data_b, work, work2, irwork, slice_status);

                if ((slice_status.lapack_info < 0) || (slice_status.is_singular)) {
                    vec_status.push_back(slice_status);
                    goto free_exit;
                }
                else if (slice_status.is_ill_conditioned) {
                    vec_status.push_back(slice_status);
                }

                break;
            }
            case St::UPPER_TRIANGULAR:
            case St::LOWER_TRIANGULAR:
            {
                char diag = 'N';
                solve_slice_triangular(uplo, diag, intn, int_nrhs, data, data_b, trans, work, irwork, slice_status);

                if ((slice_status.lapack_info < 0) || (slice_status.is_singular)) {
                    vec_status.push_back(slice_status);
                    goto free_exit;
                }
                else if (slice_status.is_ill_conditioned) {
                    vec_status.push_back(slice_status);
                }

                zero_other_triangle(uplo, data, intn);
                break;
            }
            case St::POS_DEF:
            {
                solve_slice_cholesky(uplo, intn, int_nrhs, data, data_b, work, irwork, slice_status);

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
                        copy_slice(scratch, slice_ptr, n, n, strides[ndim-2], strides[ndim-1]);
                        swap_cf(scratch, data, n, n, n);
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
                solve_slice_sym_herm(uplo, intn, int_nrhs, data, data_b, ipiv, work, irwork, lwork, is_symm_not_herm, slice_status);

                if ((slice_status.lapack_info < 0) || (slice_status.is_singular )) {
                    vec_status.push_back(slice_status);
                    goto free_exit;
                }
                else if (slice_status.is_ill_conditioned) {
                    vec_status.push_back(slice_status);
                }

                if (is_symm_not_herm) {
                    fill_other_triangle_noconj(uplo, data, intn);
                }
                else {
                    fill_other_triangle(uplo, data, intn);
                }
                break;
            }
            default:
            {
                // general matrix solve
                solve_slice_general(intn, int_nrhs, data, ipiv, data_b, trans, irwork, work, slice_status);

                if ((slice_status.lapack_info != 0) || slice_status.is_singular || slice_status.is_ill_conditioned) {
                    // some problem detected, store data to report
                    vec_status.push_back(slice_status);
                }

            }
        }

        if (slice_status.is_singular == 1) {
            // nan_matrix(data, n);
            goto free_exit;     // fail fast and loud
        }

        // Swap back to the C order
        copy_slice_F_to_C(&ret_data[idx*n*nrhs], data_b, n, nrhs);
    }

free_exit:
    free(buffer);
    free(irwork);
    free(ipiv);
    return 1;
}
