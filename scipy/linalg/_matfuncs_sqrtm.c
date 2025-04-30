#include "_matfuncs_sqrtm.h"

static int sqrtm_recursion_s(float* T, npy_intp bign, npy_intp n);
static int sqrtm_recursion_d(double* T, npy_intp bign, npy_intp n);
static int sqrtm_recursion_c(SCIPY_C* T, npy_intp bign, npy_intp n);
static int sqrtm_recursion_z(SCIPY_Z* T, npy_intp bign, npy_intp n);
static inline void zebra_pattern_s(float* restrict data, const Py_ssize_t n) {for (Py_ssize_t i=n-1; i>=0; i--) { data[2*i] = data[i]; data[2*i + 1] = 0.0f; }}
static inline void zebra_pattern_d(double* restrict data, const Py_ssize_t n) {for (Py_ssize_t i=n-1; i>=0; i--) { data[2*i] = data[i]; data[2*i + 1] = 0.0; }}

void
matrix_squareroot_s(const PyArrayObject* ap_Am, float* restrict ret_data, int* isIllconditioned, int* isSingular, int* sq_info, int* view_as_complex)
{
    float aa, bb, cc, dd, cs, sn;
    // Setting indicators to False.
    // isComplex is set to True if the square-root turns out to be complex-valued (if exists)
    // isIllconditioned is set to True if the solution of the Sylvester equation is ill-conditioned
    // isSingular is set to True if the input is exactly singular
    *view_as_complex = 0;
    *isIllconditioned = 0;
    *isSingular = 0;
    int isComplex = 0;
    int upcasted_to_complex = 0;

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
    int info = 0, sdim = 0, lwork = -1, intn = (int)n;
    float tmp_float = 0.0f, one = 1.0f, zero = 0.0f;
    sgees_("V", "N", NULL, &intn, NULL, &intn, &sdim, NULL, NULL, NULL, &intn, &tmp_float, &lwork, NULL, &info);
    // Improbable to fail at lwork query but check anyway
    if (info != 0) { *sq_info = -100; return;}
    lwork = (int)tmp_float;
    // 2n*n + 2n*n for data and vs for potentially converting to complex if needed
    // n + n for wr, wi (needed for gees calls)
    // lwork for work
    size_t buffer_size = 4*n*n + 2*n + lwork;
    float* restrict buffer = malloc(buffer_size*sizeof(float));
    if (buffer == NULL) { *sq_info = -101; return; }

    // --------------------------------------------------------------------
    // Pointers to the variables inside the buffer for easy access
    // --------------------------------------------------------------------

    // Since we already allocated n*2n for the input, we can use the first and
    // second n*n for C to F conversion when needed.
    float* restrict data = &buffer[0];
    float* restrict data2 = &buffer[n*n];
    float* restrict vs = &buffer[2*n*n];
    float* restrict wr = &buffer[4*n*n];
    float* restrict wi = &buffer[4*n*n + n];
    float* restrict work = &buffer[4*n*n + 2*n];
    // Cover the first 2n*n and vs area with a complex pointer just in case
    SCIPY_C* complex_data = &((SCIPY_C*)buffer)[0];
    SCIPY_C* complex_vs = &((SCIPY_C*)buffer)[n*n];


    /*====================================================================
    |                    MAIN nxn SLICE LOOP                             |
    ====================================================================*/
    for (npy_intp idx = 0; idx < outer_size; idx++) {
        /*------------------------------------------------------------------
        Copy the current slice into buffer in F-layout for LAPACK calls. Note
        that the input can be in C or F layout and not necessarily contiguous.
        Since the input is predominantly C contiguous, we copy the data in
        C-layout first and then transpose it to F-layout.

        The output is always a C contiguous array hence we visit slices in
        C-order for easy insertion to the output array; by just offsetting
        idx*n*n * sizeof(dtype). For this reason, we  start the offset
        calculation from the inner-most outer dimension.
        ------------------------------------------------------------------*/
        npy_intp offset = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset += (temp_idx % shape[i]) * strides[i];
            temp_idx /= shape[i];
        }
        float* restrict slice_ptr = (float*)(Am_data + (offset/sizeof(float)));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                data2[i * n + j] = *(slice_ptr + (i*strides[ndim - 2]/sizeof(float)) + (j*strides[ndim - 1]/sizeof(float)));
            }
        }
        swap_cf_s(data2, data, n, n, n);

        isComplex = 0;
        // ------------------------------------------------------------------------
        // Check if array is (quasi)upper triangular.
        // ------------------------------------------------------------------------
        int isSchur = isschurf(data, n);

        if (!isSchur)
        {
            sgees_("V", "N", NULL, &intn, data, &intn, &sdim, wr, wi, vs, &intn, work, &lwork, NULL, &info);
            if (info != 0)
            {
                free(buffer);
                *sq_info = -102;
                return;
            }
        } else {
            // Manually get the eigenvalues with possible 2x2 blocks on the diagonal
            for (Py_ssize_t col = 0; col < n; col++)
            {
                if (col == n - 1)
                {
                    // Last column 1x1 block
                    wr[col] = data[col*n + col];
                    wi[col] = 0.0f;
                } else if (data[col*n + col + 1] == 0.0f) {
                    // Real eigenvalue
                    wr[col] = data[col*n + col];
                    wi[col] = 0.0f;
                } else {
                    // 2x2 block
                    aa = data[col*n + col];
                    bb = data[(col + 1)*n + col];
                    cc = data[col*n + col + 1];
                    dd = data[(col + 1)*n + col + 1];
                    slanv2_(&aa, &bb, &cc, &dd, &wr[col], &wi[col], &wr[col+1], &wi[col+1], &cs, &sn);
                    col++;
                }
            }
        }

        // Check for singularity and negative eigenvalues on the real axis
        for (Py_ssize_t i = 0; i < n; i++)
        {
            if (wi[i] == 0.0f) {
                if (wr[i] < 0.0f) { isComplex = 1; }
                if (wr[i] == 0.0f) { *isSingular = 1; }
            }
        }

        // Compute the square root of the Schur form
        if (isComplex)
        {
            // -- Convert real to complex
            // The squareroot is determined to be complex valued.
            // The data is already in F-layout so we need to insert a zero
            // between each float to make it complex-valued array with imaginary
            // part as zero. Same is true for the Schur vectors in vs.

            // Hence data n*n will become 2n*n with the zebra pattern
            // [real, real, real, ...] --> [real, 0, real, 0, ...]
            // We start from the end of the 2n*n array and insert zeros.
            zebra_pattern_s(buffer, n*n);
            // If needed, do the same for the Schur vectors
            if (!isSchur) { zebra_pattern_s(vs, n*n); }

            // Convert real schur form to complex schur form in float arithmetic
            for (Py_ssize_t col = 0; col < (n - 1); col++)
            {
                if (data[2*(col*n + col + 1)] == 0.0f)
                {
                    // Real eigenvalue
                    continue;
                } else {
                    // The only guarantee we have is c != 0, hence we use
                    // [lambda - d, c] eigenvector for triangularization as
                    // [b, \lambda - a] might become 0.

                    // 2x2-block eigenvalues are still in wr and wi. So we can
                    // avoid complex arithmetic.
                    float u2 = data[2*(col*n + col + 1)];                       // c
                    float d_term = data[2*((col+1)*n + col + 1)];               // d
                    float mag = hypotf(hypotf(wr[col] - d_term, wi[col]), u2);  // ||[lambda - d, c]||
                    float u1re = (wr[col] - d_term) / mag;                      // (lambda - d).real / ||lambda - d||
                    float u1im = wi[col] / mag;                                 // (lambda - d).imag / ||lambda - d||
                    u2 /= mag;                                                  // c / ||lambda - d||

                    // Apply conjugate from left
                    // [u1', u2] @ [temp1] = [temp1*u1' + temp2*u2]
                    // [-u2, u1]   [temp2]   [-temp1*u2 + temp2*u1]
                    for (Py_ssize_t i = col; i < n; i++)
                    {
                        float temp1 = data[2*i*n + 2*col];
                        float temp2 = data[2*i*n + 2*col + 2];
                        data[2*i*n + 2*col    ] =  u1re*temp1 + u2*temp2;
                        data[2*i*n + 2*col + 1] = -u1im*temp1;
                        data[2*i*n + 2*col + 2] = -u2*temp1   + u1re*temp2;
                        data[2*i*n + 2*col + 3] =  u1im*temp2;
                    }

                    // Apply from right both to data and vs, temps are complex valued
                    // [temp1, temp2] @ [u1, -u2] = [temp1*u1 + temp2*u2, -temp1*u2 + temp2*u1']
                    //                  [u2, u1']

                    // (x + yi) (w - zi) = (xw + yz) + (yw - xz)i
                    // (x + yi) (w + zi) = (xw - yz) + (yw + xz)i
                    for (Py_ssize_t i = 0; i <= col + 1; i++)
                    {
                        float temp1re = data[2*n*col       + 2*i];
                        float temp1im = data[2*n*col       + 2*i + 1];
                        float temp2re = data[2*n*(col + 1) + 2*i];
                        float temp2im = data[2*n*(col + 1) + 2*i + 1];
                        data[2*n*col       + 2*i    ] =  temp1re*u1re - temp1im*u1im + temp2re*u2;
                        data[2*n*col       + 2*i + 1] =  temp1re*u1im + temp1im*u1re + temp2im*u2;
                        data[2*n*(col + 1) + 2*i    ] = -temp1re*u2   + temp2re*u1re + temp2im*u1im;
                        data[2*n*(col + 1) + 2*i + 1] = -temp1im*u2   - temp2re*u1im + temp2im*u1re;
                    }

                    // Clean up the noise, zero out the subdiagonal
                    data[2*n*col + 2*col + 2] = 0.0f;
                    data[2*n*col + 2*col + 3] = 0.0f;

                    if (!isSchur)
                    {
                        for (Py_ssize_t i = 0; i < n; i++)
                        {
                            float temp1 = vs[2*n*col       + 2*i];
                            float temp2 = vs[2*n*(col + 1) + 2*i];
                            vs[2*n*col       + 2*i    ] =  u1re*temp1 + u2*temp2;
                            vs[2*n*col       + 2*i + 1] =  u1im*temp1;
                            vs[2*n*(col + 1) + 2*i    ] = -u2*temp1   + u1re*temp2;
                            vs[2*n*(col + 1) + 2*i + 1] = -u1im*temp2;
                        }
                    }
                    // Skip next column
                    col++;
                }
            }

            // Compute the complex square root
            info = sqrtm_recursion_c(complex_data, n, n);
            // Unapply the Schur decomposition, use the return array current
            // slice as scratch space for the matrix multiplications
            if (!isSchur)
            {
                // ret_data = vs * data
                // data = ret_data * vs^H
                SCIPY_C c_one = CPLX_C(1.0f, 0.0f);
                SCIPY_C c_zero = CPLX_C(0.0f, 0.0f);
                cgemm_("N", "N", &intn, &intn, &intn, &c_one, complex_vs, &intn, complex_data, &intn, &c_zero, &((SCIPY_C*)ret_data)[idx*n*n], &intn);
                cgemm_("N", "C", &intn, &intn, &intn, &c_one, &((SCIPY_C*)ret_data)[idx*n*n], &intn, complex_vs, &intn, &c_zero, complex_data, &intn);
            }

        } else {
            // Compute the square root
            info = sqrtm_recursion_s(data, n, n);
            if (!isSchur)
            {
                // Apply the Schur decomposition, use the return array, current
                // slice as scratch space for the matrix multiplication.
                int new_address = (upcasted_to_complex? 2*idx*n*n : idx*n*n);
                sgemm_("N", "N", &intn, &intn, &intn, &one, vs, &intn, data, &intn, &zero, &ret_data[new_address], &intn);
                sgemm_("N", "T", &intn, &intn, &intn, &one, &ret_data[new_address], &intn, vs, &intn, &zero, data, &intn);
            }
        }

        // Set the isIllconditioned flag if the Sylvester solver failed
        if (info != 0) { *isIllconditioned = 1; }
        // Now that we have the result, if this is the first result with complex
        // values, we need to upcast the return array too hence it gets the zebra
        // pattern too.
        if (isComplex)
        {
            if (!upcasted_to_complex)
            {
                // Signal the caller to view the return array as complex and
                // avoid doing this multiple times.
                *view_as_complex = 1;

                // Introduce zebra pattern until end of current slice, (idx+1)*n*n
                // by going backwards.
                zebra_pattern_s(ret_data, n*n*(idx + 1));
                // Avoid doing this multiple times
                upcasted_to_complex = 1;
            }
            // Copy the complex result to the return array in C layout
            swap_cf_c(complex_data, &((SCIPY_C*)ret_data)[idx*n*n], n, n, n);
        } else {

            // Copy the real result to the return array in C layout somewhere
            // depending on whether the return array was upcasted to complex
            // previously.
            swap_cf_s(data, &ret_data[(upcasted_to_complex? 2*idx*n*n : idx*n*n)], n, n, n);

            // If the return array was upcasted to complex previously, we
            // need to add the zebra pattern for the current real slice which
            // is now living at idx*(2*n*n).
            if (upcasted_to_complex) {
                zebra_pattern_s(&ret_data[idx*2*n*n], n*n);
            }
        }
    }
    /*====================================================================
    |                  END OF nxn SLICE LOOP                             |
    ====================================================================*/
    free(buffer);
    return;
}


void
matrix_squareroot_d(const PyArrayObject* ap_Am, double* restrict ret_data, int* isIllconditioned, int* isSingular, int* sq_info, int* view_as_complex)
{
    double aa, bb, cc, dd, cs, sn;
    *view_as_complex = 0;
    *isIllconditioned = 0;
    *isSingular = 0;
    int isComplex = 0;
    int upcasted_to_complex = 0;

    double* restrict Am_data = (double*)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp* shape = PyArray_SHAPE(ap_Am);
    npy_intp n = shape[ndim - 1];
    npy_intp* restrict strides = PyArray_STRIDES(ap_Am);

    npy_intp outer_size = 1;
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }
    int info = 0, sdim = 0, lwork = -1, intn = (int)n;
    double tmp_float = 0.0, one = 1.0, zero = 0.0;
    dgees_("V", "N", NULL, &intn, NULL, &intn, &sdim, NULL, NULL, NULL, &intn, &tmp_float, &lwork, NULL, &info);
    if (info != 0) { *sq_info = -100; return; }
    lwork = (int)tmp_float;
    size_t buffer_size = 4*n*n + 2*n + lwork;
    double* restrict buffer = malloc(buffer_size*sizeof(double));
    if (buffer == NULL) { *sq_info = -101; return; }

    double* restrict data = &buffer[0];
    double* restrict data2 = &buffer[n*n];
    double* restrict vs = &buffer[2*n*n];
    double* restrict wr = &buffer[4*n*n];
    double* restrict wi = &buffer[4*n*n + n];
    double* restrict work = &buffer[4*n*n + 2*n];
    SCIPY_Z* complex_data = &((SCIPY_Z*)buffer)[0];
    SCIPY_Z* complex_vs = &((SCIPY_Z*)buffer)[n*n];

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
                data2[i * n + j] = *(slice_ptr + (i*strides[ndim - 2]/sizeof(double)) + (j*strides[ndim - 1]/sizeof(double)));
            }
        }
        swap_cf_d(data2, data, n, n, n);

        isComplex = 0;
        int isSchur = isschur(data, n);

        if (!isSchur)
        {
            dgees_("V", "N", NULL, &intn, data, &intn, &sdim, wr, wi, vs, &intn, work, &lwork, NULL, &info);
            if (info != 0)
            {
                free(buffer);
                *sq_info = -102;
                return;
            }
        } else {

            for (Py_ssize_t col = 0; col < n; col++)
            {
                if (col == n - 1)
                {
                    wr[col] = data[col*n + col];
                    wi[col] = 0.0;
                } else if (data[col*n + col + 1] == 0.0) {
                    wr[col] = data[col*n + col];
                    wi[col] = 0.0;
                } else {
                    aa = data[col*n + col];
                    bb = data[(col + 1)*n + col];
                    cc = data[col*n + col + 1];
                    dd = data[(col + 1)*n + col + 1];
                    dlanv2_(&aa, &bb, &cc, &dd, &wr[col], &wi[col], &wr[col+1], &wi[col+1], &cs, &sn);
                    col++;
                }
            }
        }

        for (Py_ssize_t i = 0; i < n; i++)
        {
            if (wi[i] == 0.0) {
                if (wr[i] < 0.0) { isComplex = 1; }
                if (wr[i] == 0.0) { *isSingular = 1; }
            }
        }

        if (isComplex)
        {
            zebra_pattern_d(buffer, n*n);
            if (!isSchur) { zebra_pattern_d(vs, n*n); }

            for (Py_ssize_t col = 0; col < (n - 1); col++)
            {
                if (data[2*(col*n + col + 1)] == 0.0)
                {
                    continue;
                } else {
                    double u2 = data[2*(col*n + col + 1)];
                    double d_term = data[2*((col+1)*n + col + 1)];
                    double mag = hypot(hypot(wr[col] - d_term, wi[col]), u2);
                    double u1re = (wr[col] - d_term) / mag;
                    double u1im = wi[col] / mag;
                    u2 = u2 / mag;

                    for (Py_ssize_t i = col; i < n; i++)
                    {
                        double temp1 = data[2*i*n + 2*col];
                        double temp2 = data[2*i*n + 2*col + 2];
                        data[2*i*n + 2*col    ] =  u1re*temp1 + u2*temp2;
                        data[2*i*n + 2*col + 1] = -u1im*temp1;
                        data[2*i*n + 2*col + 2] = -u2*temp1   + u1re*temp2;
                        data[2*i*n + 2*col + 3] =  u1im*temp2;
                    }

                    for (Py_ssize_t i = 0; i <= col + 1; i++)
                    {
                        double temp1re = data[2*n*col       + 2*i];
                        double temp1im = data[2*n*col       + 2*i + 1];
                        double temp2re = data[2*n*(col + 1) + 2*i];
                        double temp2im = data[2*n*(col + 1) + 2*i + 1];
                        data[2*n*col       + 2*i    ] =  temp1re*u1re - temp1im*u1im + temp2re*u2;
                        data[2*n*col       + 2*i + 1] =  temp1re*u1im + temp1im*u1re + temp2im*u2;
                        data[2*n*(col + 1) + 2*i    ] = -temp1re*u2   + temp2re*u1re + temp2im*u1im;
                        data[2*n*(col + 1) + 2*i + 1] = -temp1im*u2   - temp2re*u1im + temp2im*u1re;
                    }

                    data[2*n*col + 2*col + 2] = 0.0;
                    data[2*n*col + 2*col + 3] = 0.0;

                    if (!isSchur)
                    {
                        for (Py_ssize_t i = 0; i < n; i++)
                        {
                            double temp1 = vs[2*n*col       + 2*i];
                            double temp2 = vs[2*n*(col + 1) + 2*i];
                            vs[2*n*col       + 2*i    ] =  u1re*temp1 + u2*temp2;
                            vs[2*n*col       + 2*i + 1] =  u1im*temp1;
                            vs[2*n*(col + 1) + 2*i    ] = -u2*temp1   + u1re*temp2;
                            vs[2*n*(col + 1) + 2*i + 1] = -u1im*temp2;
                        }
                    }
                    col++;
                }
            }
            info = sqrtm_recursion_z(complex_data, n, n);

            if (!isSchur)
            {
                SCIPY_Z c_one = CPLX_Z(1.0, 0.0);
                SCIPY_Z c_zero = CPLX_Z(0.0, 0.0);
                zgemm_("N", "N", &intn, &intn, &intn, &c_one, complex_vs, &intn, complex_data, &intn, &c_zero, &((SCIPY_Z*)ret_data)[idx*n*n], &intn);
                zgemm_("N", "C", &intn, &intn, &intn, &c_one, &((SCIPY_Z*)ret_data)[idx*n*n], &intn, complex_vs, &intn, &c_zero, complex_data, &intn);
            }
        } else {
            info = sqrtm_recursion_d(data, n, n);
            if (!isSchur)
            {
                int new_address = (upcasted_to_complex? 2*idx*n*n : idx*n*n);
                dgemm_("N", "N", &intn, &intn, &intn, &one, vs, &intn, data, &intn, &zero, &ret_data[new_address], &intn);
                dgemm_("N", "T", &intn, &intn, &intn, &one, &ret_data[new_address], &intn, vs, &intn, &zero, data, &intn);
            }
        }

        if (info != 0) { *isIllconditioned = 1; }
        if (isComplex)
        {
            if (!upcasted_to_complex)
            {
                *view_as_complex = 1;
                zebra_pattern_d(ret_data, n*n*(idx + 1));
                upcasted_to_complex = 1;
            }
            swap_cf_z(complex_data, &((SCIPY_Z*)ret_data)[idx*n*n], n, n, n);
        } else {
            swap_cf_d(data, &ret_data[(upcasted_to_complex? 2*idx*n*n : idx*n*n)], n, n, n);
            if (upcasted_to_complex) { zebra_pattern_d(&ret_data[2*idx*n*n], n*n); }
        }
    }
    free(buffer);
    return;
}


void
matrix_squareroot_c(const PyArrayObject* ap_Am, SCIPY_C* restrict ret_data, int* isIllconditioned, int* isSingular, int* sq_info, int* unused)
{

    *isIllconditioned = 0;
    *isSingular = 0;

    SCIPY_C* restrict Am_data = (SCIPY_C*)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp* shape = PyArray_SHAPE(ap_Am);
    npy_intp n = shape[ndim - 1];
    npy_intp* restrict strides = PyArray_STRIDES(ap_Am);

    npy_intp outer_size = 1;
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }
    int info = 0, sdim = 0, lwork = -1, intn = (int)n;
    SCIPY_C tmp_float = CPLX_C(0.0f, 0.0f), cone = CPLX_C(1.0f, 0.0f), czero = CPLX_C(0.0f, 0.0f);
    cgees_("V", "N", NULL, &intn, NULL, &intn, &sdim, NULL, NULL, &intn, &tmp_float, &lwork, NULL, NULL, &info);
    if (info != 0) { *sq_info = -100; return; }

    lwork = (int)crealf(tmp_float);
    size_t buffer_size = 2*n*n + 2*n + lwork;
    SCIPY_C* buffer = malloc(buffer_size*sizeof(SCIPY_C));
    if (buffer == NULL) { *sq_info = -101; return; }


    SCIPY_C* restrict data = &buffer[0];
    SCIPY_C* restrict vs = &buffer[n*n];
    SCIPY_C* restrict w = &buffer[2*n*n];
    float* restrict rwork = &((float*)buffer)[2*n*n + n];
    SCIPY_C* restrict work = &buffer[2*n*n + 2*n];

    for (npy_intp idx = 0; idx < outer_size; idx++) {
        npy_intp offset = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset += (temp_idx % shape[i]) * strides[i];
            temp_idx /= shape[i];
        }
        SCIPY_C* restrict slice_ptr = (SCIPY_C*)(Am_data + (offset/sizeof(SCIPY_C)));
        for (Py_ssize_t i = 0; i < n; i++) {
            for (Py_ssize_t j = 0; j < n; j++) {
                vs[i * n + j] = *(slice_ptr + (i*strides[ndim - 2]/sizeof(SCIPY_C)) + (j*strides[ndim - 1]/sizeof(SCIPY_C)));
            }
        }

        swap_cf_c(vs, data, n, n, n);

        // Check if array is upper triangular
        int isSchur = 1;
        for (Py_ssize_t i = 0; i < n - 1; i++)
        {
            for (Py_ssize_t j = i + 1; j < n; j++)
            {
                if ((crealf(data[i*n + j]) != 0.0f) || (cimagf(data[i*n + j]) != 0.0f))
                {
                    isSchur = 0;
                    break;
                }
            }
            if (!isSchur) { break; }
        }

        if (!isSchur)
        {
            cgees_("V", "N", NULL, &intn, data, &intn, &sdim, w, vs, &intn, work, &lwork, rwork, NULL, &info);
            if (info != 0)
            {
                free(buffer);
                *sq_info = -102;
                return;
            }
        } else {
            for (Py_ssize_t col = 0; col < n; col++)
            {
                w[col] = data[col*n + col];
            }
        }

        for (Py_ssize_t i = 0; i < n; i++)
        {
            if ((cimagf(w[i]) == 0.0f) && (crealf(w[i]) == 0.0f)) { *isSingular = 1; }
        }

        info = sqrtm_recursion_c(data, n, n);
        if (!isSchur)
        {
            cgemm_("N", "N", &intn, &intn, &intn, &cone, vs, &intn, data, &intn, &czero, &ret_data[idx*n*n], &intn);
            cgemm_("N", "C", &intn, &intn, &intn, &cone, &ret_data[idx*n*n], &intn, vs, &intn, &czero, data, &intn);
        }

        if (info != 0) { *isIllconditioned = 1; }
        swap_cf_c(data, &ret_data[idx*n*n], n, n, n);
    }
    free(buffer);
    return;
}


void
matrix_squareroot_z(const PyArrayObject* ap_Am, SCIPY_Z* restrict ret_data, int* isIllconditioned, int* isSingular, int* sq_info, int* unused)
{

    *isIllconditioned = 0;
    *isSingular = 0;

    SCIPY_Z* restrict Am_data = (SCIPY_Z*)PyArray_DATA(ap_Am);
    int ndim = PyArray_NDIM(ap_Am);
    npy_intp* shape = PyArray_SHAPE(ap_Am);
    npy_intp n = shape[ndim - 1];
    npy_intp* restrict strides = PyArray_STRIDES(ap_Am);

    npy_intp outer_size = 1;
    if (ndim > 2)
    {
        for (int i = 0; i < ndim - 2; i++) { outer_size *= shape[i];}
    }
    int info = 0, sdim = 0, lwork = -1, intn = (int)n;
    SCIPY_Z tmp_float = CPLX_Z(0.0f, 0.0f), cone = CPLX_Z(1.0f, 0.0f), czero = CPLX_Z(0.0f, 0.0f);
    zgees_("V", "N", NULL, &intn, NULL, &intn, &sdim, NULL, NULL, &intn, &tmp_float, &lwork, NULL, NULL, &info);
    if (info != 0) { *sq_info = -100; return; }

    lwork = (int)creal(tmp_float);
    size_t buffer_size = 2*n*n + 2*n + lwork;
    SCIPY_Z* buffer = malloc(buffer_size*sizeof(SCIPY_Z));
    if (buffer == NULL) { *sq_info = -101; return; }


    SCIPY_Z* restrict data = &buffer[0];
    SCIPY_Z* restrict vs = &buffer[n*n];
    SCIPY_Z* restrict w = &buffer[2*n*n];
    double* restrict rwork = &((double*)buffer)[2*n*n + n];
    SCIPY_Z* restrict work = &buffer[2*n*n + 2*n];

    for (npy_intp idx = 0; idx < outer_size; idx++) {
        npy_intp offset = 0;
        npy_intp temp_idx = idx;
        for (int i = ndim - 3; i >= 0; i--) {
            offset += (temp_idx % shape[i]) * strides[i];
            temp_idx /= shape[i];
        }
        SCIPY_Z* restrict slice_ptr = (SCIPY_Z*)(Am_data + (offset/sizeof(SCIPY_Z)));
        for (Py_ssize_t i = 0; i < n; i++) {
            for (Py_ssize_t j = 0; j < n; j++) {
                vs[i * n + j] = *(slice_ptr + (i*strides[ndim - 2]/sizeof(SCIPY_Z)) + (j*strides[ndim - 1]/sizeof(SCIPY_Z)));
            }
        }

        swap_cf_z(vs, data, n, n, n);

        // Check if array is upper triangular
        int isSchur = 1;
        for (Py_ssize_t i = 0; i < n - 1; i++)
        {
            for (Py_ssize_t j = i + 1; j < n; j++)
            {
                if ((creal(data[i*n + j]) != 0.0f) || (cimag(data[i*n + j]) != 0.0f))
                {
                    isSchur = 0;
                    break;
                }
            }
            if (!isSchur) { break; }
        }

        if (!isSchur)
        {
            zgees_("V", "N", NULL, &intn, data, &intn, &sdim, w, vs, &intn, work, &lwork, rwork, NULL, &info);
            if (info != 0)
            {
                free(buffer);
                *sq_info = -102;
                return;
            }
        } else {
            for (Py_ssize_t col = 0; col < n; col++)
            {
                w[col] = data[col*n + col];
            }
        }

        for (Py_ssize_t i = 0; i < n; i++)
        {
            if ((cimag(w[i]) == 0.0f) && (creal(w[i]) == 0.0f)) { *isSingular = 1; }
        }

        info = sqrtm_recursion_z(data, n, n);
        if (!isSchur)
        {
            zgemm_("N", "N", &intn, &intn, &intn, &cone, vs, &intn, data, &intn, &czero, &ret_data[idx*n*n], &intn);
            zgemm_("N", "C", &intn, &intn, &intn, &cone, &ret_data[idx*n*n], &intn, vs, &intn, &czero, data, &intn);
        }

        if (info != 0) { *isIllconditioned = 1; }
        swap_cf_z(data, &ret_data[idx*n*n], n, n, n);
    }
    free(buffer);
    return;
}


/*
 * sqrtm_recursion_(s|d|c|z) are functions that compute the square root of a
 * Fortran-ordered Schur matrix by recursion.
 * For triangular Schur form arrays, the method described in [1] is used.
 * For real valued quasi-upper triangular arrays, the method described in [2]
 * is integrated into the method of [1].
 *
 * [1] : https://doi.org/10.1007/978-3-642-36803-5_12
 * [2] : https://doi.org/10.1016/0024-3795(87)90118-2
 *
 * Input:
 * T : Pointer to the Schur array in Fortran layout
 * bign : Overall leading dimension of the input array
 * n : Current recursion size of the partial array
 *
 * Output:
 * info: 0 if successful, 1 if the Sylvester solver returned nonzero info.
 *
 */
int
sqrtm_recursion_s(float* T, npy_intp bign, npy_intp n)
{
    int i, i1 = 0, i2 = 0, info = 0, j, halfn, otherhalfn, int1 = 1, intbign = (int)bign;
    float scale = 0.0, a, b, c, d, alpha, theta, mu;
    if (n == 1)
    {
        // Scalar block
        T[0] = sqrtf(T[0]);
        return 0;
    } else if (n == 2) {
        a = T[0];
        c = T[1];
        b = T[bign];
        d = T[bign + 1];
        if (c == 0.0)
        {
            // Triangular 2x2 block
            if ((a == 0.0) && (d == 0.0) && (b == 0.0)) { return 0; }
            if ((a == 0.0) && (d == 0.0)) { T[bign] = INFINITY; return 0; }
            T[0] = sqrtf(a);
            T[bign + 1] = sqrtf(d);
            T[bign] = T[bign] / (T[0] + T[bign + 1]);
            return 0;
        } else {
            // Complex eigenvalues
            theta = a;            // since a == d
            mu = sqrtf(-b*c);     // since b*c < 0
            if ((theta == 0.0) && (mu == 0.0)) { return 0; }
            if (theta > 0.0)
            {
                alpha = sqrtf((theta + hypotf(theta, mu)) / 2);
            } else {
                alpha = mu / sqrtf(2 * (-theta + hypotf(theta, mu)));
            }
            T[0] = alpha;
            T[1] = T[1] / (2*alpha);
            T[bign] = T[bign] / (2*alpha);
            T[bign + 1] = alpha;
            return 0;
        }
    } else {
        // Recursion
        halfn = (int)(n / 2);
        // Don't split the 2x2 block, check if the entry under separation is zero
        if (T[(halfn-1)*bign + halfn] != 0.0) { halfn++; }
        otherhalfn = (int)(n - halfn);

        // Top left block
        i1 = sqrtm_recursion_s(T, bign, halfn);
        // Bottom right block
        i2 = sqrtm_recursion_s(&T[halfn*(bign + 1)], bign, otherhalfn);

        // Solve the Sylvester equation for the top right block
        strsyl_("N", "N", &int1, &halfn, &otherhalfn, T, &intbign, &T[halfn*(intbign + 1)], &intbign, &T[halfn*intbign], &intbign, &scale, &info);
        if (scale != 1.0)
        {
            for (i = 0; i < otherhalfn; i++)
            {
                for (j = 0; j < halfn; j++)
                {
                    T[(halfn + i)*bign + j] *= scale;
                }
            }
        }
    }
    // Propagate the Sylvester solver info up the recursion stack
    if (i1 || i2 || info)
    {
        return 1;
    } else {
        return 0;
    }
}


int
sqrtm_recursion_d(double* T, npy_intp bign, npy_intp n)
{
    int i, i1 = 0, i2 = 0, info = 0, j, halfn, otherhalfn, int1 = 1, intbign = (int)bign;
    double scale = 0.0, a, b, c, d, alpha, theta, mu;
    if (n == 1)
    {
        // Scalar block
        T[0] = sqrt(T[0]);
        return 0;
    } else if (n == 2) {
        a = T[0];
        c = T[1];
        b = T[bign];
        d = T[bign + 1];
        if (c == 0.0)
        {
            // Triangular 2x2 block
            if ((a == 0.0) && (d == 0.0) && (b == 0.0)) { return 0; }
            if ((a == 0.0) && (d == 0.0)) { T[bign] = INFINITY; return 0; }
            T[0] = sqrt(a);
            T[bign + 1] = sqrt(d);
            T[bign] = T[bign] / (T[0] + T[bign + 1]);
            return 0;
        } else {
            // Complex eigenvalues
            theta = a;            // since a == d
            mu = sqrt(-b*c);     // since b*c < 0
            if ((theta == 0.0) && (mu == 0.0)) { return 0; }
            if (theta > 0.0)
            {
                alpha = sqrt((theta + hypot(theta, mu)) / 2);
            } else {
                alpha = mu / sqrt(2 * (-theta + hypot(theta, mu)));
            }
            T[0] = alpha;
            T[1] = T[1] / (2*alpha);
            T[bign] = T[bign] / (2*alpha);
            T[bign + 1] = alpha;
            return 0;
        }
    } else {
        // Recursion
        halfn = (int)(n / 2);
        // Don't split the 2x2 block, check if the entry under separation is zero
        if (T[(halfn-1)*bign + halfn] != 0.0) { halfn++; }
        otherhalfn = (int)(n - halfn);

        // Top left block
        i1 = sqrtm_recursion_d(T, bign, halfn);
        // Bottom right block
        i2 = sqrtm_recursion_d(&T[halfn*(bign + 1)], bign, otherhalfn);

        // Solve the Sylvester equation for the top right block
        dtrsyl_("N", "N", &int1, &halfn, &otherhalfn, T, &intbign, &T[halfn*(intbign + 1)], &intbign, &T[halfn*intbign], &intbign, &scale, &info);

        if (scale != 1.0)
        {
            for (i = 0; i < otherhalfn; i++)
            {
                for (j = 0; j < halfn; j++)
                {
                    T[(halfn + i)*bign + j] *= scale;
                }
            }
        }
    }
    // Propagate the Sylvester solver info up the recursion stack
    if (i1 || i2 || info)
    {
        return 1;
    } else {
        return 0;
    }
}


int
sqrtm_recursion_c(SCIPY_C* T, npy_intp bign, npy_intp n)
{
    int i, i1 = 0, i2 = 0, info = 0, j, halfn, otherhalfn, int1 = 1, intbign = (int)bign;
    float scale = 0.0;
    if (n == 1)
    {
        // Scalar block
        T[0] = csqrtf(T[0]);
        return 0;
    } else if (n == 2) {
        // 2x2 upper triangular special case

        // Check if the 2x2 block is exactly 0.0 and if so do nothing.
        if ((cabsf(T[0]) == 0.0) && (cabsf(T[bign + 1]) == 0.0) && (cabsf(T[bign]) == 0.0)) { return 0; }

        T[0] = csqrtf(T[0]);
        // n + 1 is the next diagonal from T[0]
        T[bign + 1] = csqrtf(T[bign + 1]);
#if defined(_MSC_VER)
        // MSVC does not support complex division
        SCIPY_C tt = CPLX_C(crealf(T[0]) + crealf(T[bign + 1]), cimagf(T[0]) + cimagf(T[bign + 1]));
        float cc = cabsf(tt);
        T[bign] = _FCmulcr(_FCmulcc(T[bign], conjf(tt)), 1 / (cc*cc));
#else
        T[bign] = T[bign] / (T[0] + T[bign + 1]);
#endif
        return 0;
    } else {
        halfn = (int)(n / 2);
        otherhalfn = (int)(n - halfn);
        SCIPY_C* T11 = T;
        SCIPY_C* T22 = &T[halfn*(bign + 1)];
        SCIPY_C* T12 = &T[halfn*bign];
        i1 = sqrtm_recursion_c(T11, bign, halfn);
        i2 = sqrtm_recursion_c(T22, bign, otherhalfn);

        ctrsyl_("N", "N", &int1, &halfn, &otherhalfn, T11, &intbign, T22, &intbign, T12, &intbign, &scale, &info);
        if (scale != 1.0)
        {
            for (i = 0; i < otherhalfn; i++)
            {
                for (j = 0; j < halfn; j++)
                {
#if defined(_MSC_VER)
                    T[(halfn + i)*bign + j] = _FCmulcr(T[(halfn + i)*bign + j], scale);
#else
                    T[(halfn + i)*bign + j] *= scale;
#endif
                }
            }
        }
    }
    // Propagate the Sylvester solver info up the recursion stack
    if (i1 || i2 || info)
    {
        return 1;
    } else {
        return 0;
    }
}


int
sqrtm_recursion_z(SCIPY_Z* T, npy_intp bign, npy_intp n)
{
    int i, i1 = 0, i2 = 0, info = 0, j, halfn, otherhalfn, int1 = 1, intbign = (int)bign;
    double scale = 0.0;
    if (n == 1)
    {
        // Scalar block
        T[0] = csqrt(T[0]);
        return 0;
    } else if (n == 2) {
        // 2x2 upper triangular special case

        // Check if the 2x2 block is exactly 0.0 and if so do nothing.
        if ((cabs(T[0]) == 0.0) && (cabs(T[bign + 1]) == 0.0) && (cabs(T[bign]) == 0.0)) { return 0; }

        T[0] = csqrt(T[0]);
        // n + 1 is the next diagonal from T[0]
        T[bign + 1] = csqrt(T[bign + 1]);
#if defined(_MSC_VER)
        // MSVC does not support complex division
        SCIPY_Z tt = CPLX_Z(creal(T[0]) + creal(T[bign + 1]), cimag(T[0]) + cimag(T[bign + 1]));
        double cc = cabs(tt);
        T[bign] = _Cmulcr(_Cmulcc(T[bign], conj(tt)), 1 / (cc*cc));
#else
        T[bign] = T[bign] / (T[0] + T[bign + 1]);
#endif
        return 0;
    } else {
        halfn = (int)(n / 2);
        otherhalfn = (int)(n - halfn);
        SCIPY_Z* T11 = T;
        SCIPY_Z* T22 = &T[halfn*(bign + 1)];
        SCIPY_Z* T12 = &T[halfn*bign];
        i1 = sqrtm_recursion_z(T11, bign, halfn);
        i2 = sqrtm_recursion_z(T22, bign, otherhalfn);
        ztrsyl_("N", "N", &int1, &halfn, &otherhalfn, T11, &intbign, T22, &intbign, T12, &intbign, &scale, &info);
        if (scale != 1.0)
        {
            for (i = 0; i < otherhalfn; i++)
            {
                for (j = 0; j < halfn; j++)
                {
#if defined(_MSC_VER)
                    T[(halfn + i)*bign + j] = _Cmulcr(T[(halfn + i)*bign + j], scale);
#else
                    T[(halfn + i)*bign + j] *= scale;
#endif
                }
            }
        }
    }
    // Propagate the Sylvester solver info up the recursion stack
    if (i1 || i2 || info)
    {
        return 1;
    } else {
        return 0;
    }
}


#undef SCIPY_Z
#undef SCIPY_C
#undef CPLX_Z
#undef CPLX_C

