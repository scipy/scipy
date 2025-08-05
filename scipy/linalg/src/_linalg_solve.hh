/*
 * Templated loops for `linalg.solve`
 */
#include "Python.h"
#include <iostream>
#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"
#include "npy_cblas.h"
#include "_npymath.hh"
#include "_common_array_utils.hh"



// Dense array solve with getrf, gecon and getrs
template<typename T>
inline CBLAS_INT solve_slice_general(
    CBLAS_INT N, CBLAS_INT NRHS, T *data, CBLAS_INT *ipiv, T *b_data, char trans, void *irwork, T *work, CBLAS_INT lwork,
    int* isIllconditioned, int* isSingular
) {
    using real_type = typename type_traits<T>::real_type;

    CBLAS_INT info;
    char norm = '1';
    real_type rcond;
    real_type anorm = norm1_(data, work, (npy_intp)N);

    getrf(&N, &N, data, &N, ipiv, &info);

    if (info == 0){
        // getrf success, check the condition number
        gecon(&norm, &N, data, &N, &anorm, &rcond, work, irwork, &info);

        if (info >= 0) {
            *isIllconditioned = (rcond != rcond) || (rcond < numeric_limits<real_type>::eps);

            // finally, solve
            getrs(&trans, &N, &NRHS, data, &N, ipiv, b_data, &N, &info);
            *isSingular = (info > 0);
        }
    }
    else if (info > 0) {
        // trf detected singularity
        *isSingular = 1;
    }

    return info;
}


// triangular solve with trtrs
template<typename T>
inline CBLAS_INT solve_slice_triangular(
    char uplo, char diag, CBLAS_INT N, CBLAS_INT NRHS, T *data,  T *b_data, char trans, T *work, void *irwork,
    int* isIllconditioned, int* isSingular
) {
    using real_type = typename type_traits<T>::real_type;

    CBLAS_INT info;
    char norm = '1';
    real_type rcond;

    trtrs(&uplo, &trans, &diag, &N, &NRHS, data, &N, b_data, &N, &info);

    *isSingular  = (info > 0);

    if(info >= 0) {

        trcon(&norm, &uplo, &diag, &N, data, &N, &rcond, work, irwork, &info);
        if (info >= 0) {
            *isIllconditioned = (rcond != rcond) || (rcond < numeric_limits<real_type>::eps);
        }
    }
    return info;
}


// Cholesky solve with potrf, pocon and potrs
template<typename T>
inline CBLAS_INT solve_slice_cholesky(
    char uplo, CBLAS_INT N, CBLAS_INT NRHS, T *data, T *b_data, T* work, void *irwork,
    int* isIllconditioned, int* isSingular
) {
    using real_type = typename type_traits<T>::real_type;

    CBLAS_INT info;
    real_type anorm = norm1_sym_herm(uplo, data, work, (npy_intp)N);

    real_type rcond;

    potrf(&uplo, &N, data, &N, &info);
    if (info == 0) {
        // potrf success
        pocon(&uplo, &N, data, &N, &anorm, &rcond, work, irwork, &info);

        if (info >= 0) {
            *isIllconditioned = (rcond != rcond) || (rcond < numeric_limits<real_type>::eps);

            // finally, invert
            potrs(&uplo, &N, &NRHS, data, &N, b_data, &N, &info);
            *isSingular = (info > 0);
        }
    }
    else if (info > 0) {
        // trf detected singularity
        *isSingular = 1;
    }

    return info;
}


template<typename T>
void _solve(PyArrayObject* ap_Am, PyArrayObject *ap_b, T* ret_data, St structure, int lower, int transposed, int overwrite_a, int* isIllconditioned, int* isSingular, int* info)
{
    using real_type = typename type_traits<T>::real_type; // float if T==npy_cfloat etc

    *isIllconditioned = 0;
    *isSingular = 0;
    char trans = transposed ? 'T' : 'N'; 
    npy_intp lower_band = 0, upper_band = 0;
    bool is_symm = false;
    char uplo = 'X';    // sentinel
    St slice_structure = St::NONE;
    bool posdef_fallback = true;

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
    CBLAS_INT intn = (CBLAS_INT)n, int_nrhs = (CBLAS_INT)nrhs, lwork = -1;

    lwork = 4*n; // gecon needs at least 4*n
    T* buffer = (T *)malloc((2*n*n + n*nrhs + lwork)*sizeof(T));
    if (NULL == buffer) { *info = -101; return; }

    // Chop the buffer into parts, one for data and one for work
    T* data = &buffer[0];
    T* scratch = &buffer[n*n];

    T *data_b = &buffer[2*n*n];
    T* work = &buffer[2*n*n + n*nrhs];

    CBLAS_INT* ipiv = (CBLAS_INT *)malloc(n*sizeof(CBLAS_INT));
    if (ipiv == NULL) { free(ipiv); *info = -102; return; }

    // {ge,po,tr}con need rwork or iwork
    void *irwork;
    if (type_traits<T>::is_complex) {
        irwork = malloc(3*n*sizeof(real_type));   // {po,tr}con need at least 3*n
    } else {
        irwork = malloc(n*sizeof(CBLAS_INT));
    }
    if (irwork == NULL) { free(irwork); *info = -102; return; }

    // normalize the structure detection inputs
    if (structure == St::POS_DEF) {
        posdef_fallback = false;
        uplo = lower ? 'L' : 'U';
    }
    else {
        if (structure == St::POS_DEF_UPPER) {
            structure = St::POS_DEF;
            uplo = 'U';
            posdef_fallback = false;
        }
        else if (structure == St::POS_DEF_LOWER) {
            structure = St::POS_DEF;
            uplo = 'L';
            posdef_fallback = false;
        }
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

            if(lower_band == 0) {
                slice_structure = St::UPPER_TRIANGULAR;
                uplo = 'U';
            } else if (upper_band == 0) {
                slice_structure = St::LOWER_TRIANGULAR;
                uplo = 'L';
            } else {
                // Check if symmetric/hermitian
                is_symm = is_sym_herm(data, n);
                if (is_symm) {
                    slice_structure = St::POS_DEF;
                    uplo = 'U';
                }
                else {
                    // give up auto-detection
                    slice_structure = St::GENERAL;
                }
            }
        }

        switch(slice_structure) {
            case St::UPPER_TRIANGULAR:
            case St::LOWER_TRIANGULAR:
            {
                char diag = 'N';
                *info = solve_slice_triangular(uplo, diag, intn, int_nrhs, data, data_b, trans, work, irwork, isIllconditioned, isSingular);

                if ((*info < 0) || (*isSingular )) { goto free_exit;}
                zero_other_triangle(uplo, data, intn);
                break;
            }
            case St::POS_DEF:
            {
                *info = solve_slice_cholesky(uplo, intn, int_nrhs, data, data_b, work, irwork, isIllconditioned, isSingular);

                if ((*info == 0) || (*isSingular == 0) ) {
                    // success
                    fill_other_triangle(uplo, data, intn);
                    break;
                }
                else { // potrf failed
                    if(posdef_fallback) {
                        // restore
                        copy_slice(scratch, slice_ptr, n, n, strides[ndim-2], strides[ndim-1]);
                        swap_cf(scratch, data, n, n, n);

                        // no break: fall back to the general solver
                    }
                    else {
                        // potrf failed but no fallback
                        break;
                    }
                }
            }
            default:
            {
                // general matrix inverse
                *info = solve_slice_general(intn, int_nrhs, data, ipiv, data_b, trans, irwork, work, lwork, isIllconditioned, isSingular);
            }
        }

        if (*isSingular == 1) {
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
    return;
}
