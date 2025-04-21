#include "_linalg_inv.h"

float
snorm1(float* A, float* work, const Py_ssize_t n)
{
    Py_ssize_t i, j;
    float temp = 0.0;
    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { work[i] = fabsf(A[i]); }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) { for (j = 0; j < n; j++) { work[j] += fabsf(A[i*n + j]); } }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (work[i] > temp) { temp = work[i]; } }
    return temp;
}

double
dnorm1(double* A, double* work, const Py_ssize_t n)
{
    Py_ssize_t i, j;
    double temp = 0.0;
    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { work[i] = fabs(A[i]); }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) { for (j = 0; j < n; j++) { work[j] += fabs(A[i*n + j]); } }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (work[i] > temp) { temp = work[i]; } }
    return temp;
}

float
cnorm1(SCIPY_C* A, float* work, const Py_ssize_t n)
{
    Py_ssize_t i, j;
    float temp = 0.0;
    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { work[i] = cabsf(A[i]); }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) { for (j = 0; j < n; j++) { work[j] += cabsf(A[i*n + j]); } }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (work[i] > temp) { temp = work[i]; } }
    return temp;
}

double
znorm1(SCIPY_Z* A, double* work, const Py_ssize_t n)
{
    Py_ssize_t i, j;
    double temp = 0.0;
    // Write absolute values of first row of A to work
    for (i = 0; i < n; i++) { work[i] = cabs(A[i]); }
    // Add absolute values of remaining rows of A to work
    for (i = 1; i < n; i++) { for (j = 0; j < n; j++) { work[j] += cabs(A[i*n + j]); } }
    temp = 0.0;
    for (i = 0; i < n; i++) { if (work[i] > temp) { temp = work[i]; } }
    return temp;
}


void _inverse_s(const PyArrayObject* ap_Am, float* restrict ret_data, int* isIllconditioned, int* isSingular, int* info)
{
    *isIllconditioned = 0;
    *isSingular = 0;
    int lower_band = 0, upper_band = 0;

    // --------------------------------------------------------------------
    // Input Array Attributes
    // --------------------------------------------------------------------
    float* restrict Am_data = (float*)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_Am);      // Array shape
    npy_intp n = shape[ndim - 1];                // Slice size
    npy_intp* restrict strides = PyArray_STRIDES(ap_Am);
    // Get the number of slices to traverse if more than one; np.prod(shape[:-2])
    npy_intp outer_size = 1;
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    float tmp_float = 0.0f;
    int intn = (int)n, lwork = -1;
    sgetri_(&intn, NULL, NULL, NULL, &tmp_float, &lwork, &info);
    if (info != 0) { *info = -100; return; }
    lwork = (int)tmp_float;
    lwork = (4*n > lwork ? 4*n : lwork); // gecon needs at least 4*n
    float* buffer = malloc((n*n + lwork)*sizeof(float));
    if (NULL == buffer) { *info = -101; return; }
    // Chop buffer into two parts, one for data and one for work
    float* restrict data = &buffer[0];
    float* restrict work = &buffer[n*n];
    int* ipiv = malloc(n*sizeof(int));
    if (ipiv == NULL) { free(buffer); *info = -102; return; }

    // Main loop to traverse the slices
    for (npy_intp idx = 0; idx < outer_size; idx++) {

        npy_intp offset = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset += (temp_idx % shape[i]) * strides[i];
            temp_idx /= shape[i];
        }
        float* restrict slice_ptr = (float*)(Am_data + (offset/sizeof(float)));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                work[i * n + j] = *(slice_ptr + (i*strides[ndim - 2]/sizeof(float)) + (j*strides[ndim - 1]/sizeof(float)));
            }
        }
        swap_cf_s(work, data, n, n, n);

        // Get the bandwidth of the slice
        bandwidth_s(data, n, n, &lower_band, &upper_band);

        // TODO: Check if the array is symmetric...

        if ((lower_band == 0) || (upper_band == 0)) {
            char uplo = (lower_band == 0) ? 'U' : 'L';
            strtri(&uplo, "N", &intn, data, &intn, &info);

            // Don't forgive LAPACK errors
            if (info != 0) { *info = -103; goto free_exit; }

            // Check for singularity
            if (info > 0) { *isSingular = 1; }
            // TODO: Add trcon here to check for condition number and set isIllconditioned

        } else {
            float anorm = snorm1(data, work, n);
            // Dense array inversion with getrf, gecon and getri
            sgetrf_(&intn, &intn, data, &intn, ipiv, &info);
            if (info < 0) { *info = -104; goto free_exit; }
            if (info == 0)  // getrf success
            {
                sgecon_("1", &intn, data, &intn, &anorm, &tmp_float, work, &ipiv, &info);
                if (info < 0) { *info = -105; goto free_exit; }
                sgetri_(&intn, data, &intn, ipiv, work, &lwork, &info);
                if (info < 0) { *info = -106; goto free_exit; }
            } else {
                *isSingular = 1;
            }

            if (*isSingular == 1) {
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        data[i * n + j] = NAN;
                    }
                }
            } else {
                // If singular, set isIllconditioned to 1
                *isIllconditioned = 1;
            }
        }
        // Swap back to original order
        swap_cf_s(data, &ret_data[idx*n*n], n, n, n);
    }

free_exit:
    free(buffer);
    free(ipiv);
    return;
}


void _inverse_d(const PyArrayObject* ap_Am, double* restrict ret_data, int* isIllconditioned, int* isSingular, int* info)
{
    *isIllconditioned = 0;
    *isSingular = 0;
    int lower_band = 0, upper_band = 0;

    // --------------------------------------------------------------------
    // Input Array Attributes
    // --------------------------------------------------------------------
    double* restrict Am_data = (double*)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);              // Number of dimensions
    npy_intp* shape = PyArray_SHAPE(ap_Am);      // Array shape
    npy_intp n = shape[ndim - 1];                // Slice size
    npy_intp* restrict strides = PyArray_STRIDES(ap_Am);
    // Get the number of slices to traverse if more than one; np.prod(shape[:-2])
    npy_intp outer_size = 1;
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }

    // --------------------------------------------------------------------
    // Workspace computation and allocation
    // --------------------------------------------------------------------
    double tmp_float = 0.0;
    int intn = (int)n, lwork = -1;
    dgetri_(&intn, NULL, NULL, NULL, &tmp_float, &lwork, &info);
    if (info != 0) { *info = -100; return; }
    lwork = (int)tmp_float;
    lwork = (4*n > lwork ? 4*n : lwork); // gecon needs at least 4*n
    double* buffer = malloc((n*n + lwork)*sizeof(double));
    if (NULL == buffer) { *info = -101; return; }
    // Chop buffer into two parts, one for data and one for work
    double* restrict data = &buffer[0];
    double* restrict work = &buffer[n*n];
    int* ipiv = malloc(n*sizeof(int));
    if (ipiv == NULL) { free(buffer); *info = -102; return; }

    // Main loop to traverse the slices
    for (npy_intp idx = 0; idx < outer_size; idx++) {

        npy_intp offset = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset += (temp_idx % shape[i]) * strides[i];
            temp_idx /= shape[i];
        }
        double* restrict slice_ptr = (double*)(Am_data + (offset/sizeof(double)));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                work[i * n + j] = *(slice_ptr + (i*strides[ndim - 2]/sizeof(double)) + (j*strides[ndim - 1]/sizeof(double)));
            }
        }
        swap_cf_d(work, data, n, n, n);

        // Check if array is triangular
        bandwidth_s(data, n, n, &lower_band, &upper_band);
        if ((lower_band == 0) || (upper_band == 0)) {
            char uplo = (lower_band == 0) ? 'U' : 'L';
            dtrtri(&uplo, "N", &intn, data, &intn, &info);

            // Don't forgive LAPACK errors
            if (info != 0) { *info = -103; goto free_exit; }

            // Check for singularity
            if (info > 0) { *isSingular = 1; }
            // TODO: Add trcon here to check for condition number and set isIllconditioned

        } else {
            float anorm = snorm1(data, work, n);
            // Dense array inversion with getrf, gecon and getri
            dgetrf_(&intn, &intn, data, &intn, ipiv, &info);
            if (info < 0) { *info = -104; goto free_exit; }
            if (info == 0)  // getrf success
            {
                dgecon_("1", &intn, data, &intn, &anorm, &tmp_float, work, &ipiv, &info);
                if (info < 0) { *info = -105; goto free_exit; }
                dgetri_(&intn, data, &intn, ipiv, work, &lwork, &info);
                if (info < 0) { *info = -106; goto free_exit; }
            } else {
                *isSingular = 1;
            }

            if (*isSingular == 1) {
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        data[i * n + j] = NAN;
                    }
                }
            } else {
                // If singular, set isIllconditioned to 1
                *isIllconditioned = 1;
            }
        }
        // Swap back to original order
        swap_cf_d(data, &ret_data[idx*n*n], n, n, n);
    }

free_exit:
    free(buffer);
    free(ipiv);
    return;
}


void _inverse_c(const PyArrayObject* ap_Am, SCIPY_C* restrict ret_data, int* isIllconditioned, int* isSingular, int* info) {}
void _inverse_z(const PyArrayObject* ap_Am, SCIPY_Z* restrict ret_data, int* isIllconditioned, int* isSingular, int* info) {}

