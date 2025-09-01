/*
 * Templated loops for `linalg.inv`
 */
#include "Python.h"
#include <iostream>
#include <vector>
#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"
#include "npy_cblas.h"
#include "_npymath.hh"
#include "_common_array_utils.hh"


// Dense array inversion with getrf, gecon and getri
template<typename T>
void invert_slice_general(
    CBLAS_INT N, T *data, CBLAS_INT *ipiv, void *irwork, T *work, CBLAS_INT lwork,
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

            // finally, invert
            getri(&N, data, &N, ipiv, work, &lwork, &info);
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
    char uplo, CBLAS_INT N, T *data, T* work, void *irwork,
    SliceStatus& status
) {
    using real_type = typename type_traits<T>::real_type;

    CBLAS_INT info;
    real_type anorm = norm1_sym_herm(uplo, data, work, (npy_intp)N);

    real_type rcond;

    potrf(&uplo, &N, data, &N, &info);

    status.lapack_info = (Py_ssize_t)info;
    if (info == 0) {
        // potrf success
        pocon(&uplo, &N, data, &N, &anorm, &rcond, work, irwork, &info);

        if (info >= 0) {
            status.rcond = (double)rcond;
            status.is_ill_conditioned = (rcond != rcond) || (rcond < numeric_limits<real_type>::eps);

            // finally, invert
            potri(&uplo, &N, data, &N, &info);
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
    char uplo, CBLAS_INT N, T *data, CBLAS_INT *ipiv, T *work, void *irwork, CBLAS_INT lwork,
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

            // finally, invert
            if (is_symm_not_herm) {
                sytri(&uplo, &N, data, &N, ipiv, work, &info);
            } else {
                hetri(&uplo, &N, data, &N, ipiv, work, &info);
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
    char uplo, char diag, CBLAS_INT N, T *data, T *work, void *irwork,
    SliceStatus& status
) {
    using real_type = typename type_traits<T>::real_type;

    CBLAS_INT info;
    char norm = '1';
    real_type rcond;

    trtri(&uplo, &diag, &N, data, &N, &info);
    status.is_singular  = (info > 0);

    status.lapack_info = (Py_ssize_t)info;
    if(info >= 0) {

        trcon(&norm, &uplo, &diag, &N, data, &N, &rcond, work, irwork, &info);
        if (info >= 0) {
            status.is_ill_conditioned = (rcond != rcond) || (rcond < numeric_limits<real_type>::eps);
            status.rcond = (double)rcond;
        }
    }
}


template<typename T>
int
_inverse(PyArrayObject* ap_Am, T* ret_data, St structure, int lower, int overwrite_a, SliceStatusVec& vec_status)
{
    using real_type = typename type_traits<T>::real_type; // float if T==npy_cfloat etc

    npy_intp lower_band = 0, upper_band = 0;
    bool is_symm_or_herm = false, is_symm_not_herm = false;
    char uplo;
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
    T tmp = numeric_limits<T>::zero;
    CBLAS_INT intn = (CBLAS_INT)n, lwork = -1, info;

    getri(&intn, NULL, &intn, NULL, &tmp, &lwork, &info);
    if (info != 0) { info = -100; return (int)info; }

    /*
     * The factor of 1.01 here mirrors
     * https://github.com/scipy/scipy/blob/v1.15.2/scipy/linalg/_basic.py#L1154
     *
     * It was added in commit
     * https://github.com/scipy/scipy/commit/dfb543c147c
     * to avoid a "curious segfault with 500x500 matrices and OpenBLAS".
     */
    lwork = (CBLAS_INT)(1.01 * real_part(tmp));

    lwork = (4*n > lwork ? 4*n : lwork); // gecon needs at least 4*n
    T* buffer = (T *)malloc((2*n*n + lwork)*sizeof(T));
    if (NULL == buffer) { info = -101; return (int)info; }

    // Chop buffer into parts, one for data and one for work
    T* data = &buffer[0];
    T* scratch = &buffer[n*n];
    T* work = &buffer[2*n*n];

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
    uplo = lower? 'L' : 'U';
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

        // detect the structure if not given
        slice_structure = structure;
        if (slice_structure == St::NONE) {
            // Get the bandwidth of the slice
            bandwidth(data, n, n, &lower_band, &upper_band);

            if(lower_band == 0) {
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
            case St::UPPER_TRIANGULAR:
            case St::LOWER_TRIANGULAR:
            {
                char diag = 'N';
                invert_slice_triangular(uplo, diag, intn, data, work, irwork, slice_status);

                if ((slice_status.lapack_info < 0) || (slice_status.is_singular )) {
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
                invert_slice_cholesky(uplo, intn, data, work, irwork, slice_status);

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
                invert_slice_sym_herm(uplo, intn, data, ipiv, work, irwork, lwork, is_symm_not_herm, slice_status);

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
                // general matrix inverse
                invert_slice_general(intn, data, ipiv, irwork, work, lwork, slice_status);

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

        // Swap back to original order
        swap_cf(data, &ret_data[idx*n*n], n, n, n);
    }

free_exit:
    free(buffer);
    free(irwork);
    free(ipiv);
    return 1;
}

