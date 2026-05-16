#include "gs.h"


void smgs(CBLAS_INT n, CBLAS_INT k, float* V, CBLAS_INT ldv, float* vnew, const CBLAS_INT* indices) {
    CBLAS_INT ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    CBLAS_INT idx = 0;
    CBLAS_INT start = indices[idx];
    CBLAS_INT end = indices[idx + 1];

    /**
     * PROPACK encodes sentinels in 1-index specific way.
     * Therefore, we have to guard for a few edge cases of 0-indexing.
     * TODO: Fix this by rewriting the indexing logic.
     */
    while ((start <= k) && ((start < end) || ((start == 0) && (end == 0) && (idx == 0)))) {
        // Orthogonalize against columns [start, end] (0-indexed)
        for (CBLAS_INT i = start; i <= end; i++) {
            // Compute projection coefficient: coef = V(:,i)' * vnew
            float coef = BLAS_FUNC(sdot)(&n, &V[i * ldv], &ione, vnew, &ione);

            // Orthogonalize: vnew = vnew - coef * V(:,i)
            float neg_coef = -coef;
            BLAS_FUNC(saxpy)(&n, &neg_coef, &V[i * ldv], &ione, vnew, &ione);
        }

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void dmgs(CBLAS_INT n, CBLAS_INT k, double* V, CBLAS_INT ldv, double* vnew, const CBLAS_INT* indices) {
    CBLAS_INT ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    CBLAS_INT idx = 0;
    CBLAS_INT start = indices[idx];
    CBLAS_INT end = indices[idx + 1];

    /**
     * PROPACK encodes sentinels in 1-index specific way.
     * Therefore, we have to guard for a few edge cases of 0-indexing.
     * TODO: Fix this by rewriting the indexing logic.
     */
    while ((start <= k) && ((start < end) || ((start == 0) && (end == 0) && (idx == 0)))) {
        // Orthogonalize against columns [start, end] (0-indexed)
        for (CBLAS_INT i = start; i <= end; i++) {
            // Compute projection coefficient: coef = V(:,i)' * vnew
            double coef = BLAS_FUNC(ddot)(&n, &V[i * ldv], &ione, vnew, &ione);

            // Orthogonalize: vnew = vnew - coef * V(:,i)
            double neg_coef = -coef;
            BLAS_FUNC(daxpy)(&n, &neg_coef, &V[i * ldv], &ione, vnew, &ione);
        }

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void cmgs(CBLAS_INT n, CBLAS_INT k, PROPACK_CPLXF_TYPE* V, CBLAS_INT ldv, PROPACK_CPLXF_TYPE* vnew, const CBLAS_INT* indices) {
    CBLAS_INT ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    CBLAS_INT idx = 0;
    CBLAS_INT start = indices[idx];
    CBLAS_INT end = indices[idx + 1];

    while ((start <= k) && ((start < end) || ((start == 0) && (end == 0) && (idx == 0)))) {
        // Orthogonalize against columns [start, end] (0-indexed)
        for (CBLAS_INT i = start; i <= end; i++) {
            // Compute projection coefficient: coef = V(:,i)^H * vnew (conjugate dot product)
            PROPACK_CPLXF_TYPE coef = cdotc_(&n, &V[i * ldv], &ione, vnew, &ione);

            // Orthogonalize: vnew = vnew - coef * V(:,i)
            PROPACK_CPLXF_TYPE neg_coef = PROPACK_cplxf(-crealf(coef), -cimagf(coef));
            BLAS_FUNC(caxpy)(&n, &neg_coef, &V[i * ldv], &ione, vnew, &ione);
        }

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void zmgs(CBLAS_INT n, CBLAS_INT k, PROPACK_CPLX_TYPE* V, CBLAS_INT ldv, PROPACK_CPLX_TYPE* vnew, const CBLAS_INT* indices) {
    CBLAS_INT ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    CBLAS_INT idx = 0;
    CBLAS_INT start = indices[idx];
    CBLAS_INT end = indices[idx + 1];

    while ((start <= k) && ((start < end) || ((start == 0) && (end == 0) && (idx == 0)))) {
        // Orthogonalize against columns [start, end] (0-indexed)
        for (CBLAS_INT i = start; i <= end; i++) {
            // Compute projection coefficient: coef = V(:,i)^H * vnew (conjugate dot product)
            PROPACK_CPLX_TYPE coef = zdotc_(&n, &V[i * ldv], &ione, vnew, &ione);

            // Orthogonalize: vnew = vnew - coef * V(:,i)
            PROPACK_CPLX_TYPE neg_coef = PROPACK_cplx(-creal(coef), -cimag(coef));
            BLAS_FUNC(zaxpy)(&n, &neg_coef, &V[i * ldv], &ione, vnew, &ione);
        }

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void scgs(CBLAS_INT n, CBLAS_INT k, float* V, CBLAS_INT ldv, float* vnew, const CBLAS_INT* indices, float* work) {
    CBLAS_INT ione = 1;
    float one = 1.0f;
    float zero = 0.0f;
    float neg_one = -1.0f;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    CBLAS_INT idx = 0;
    CBLAS_INT start = indices[idx];
    CBLAS_INT end = indices[idx + 1];

    while (start < k && start >= 0 && start <= end) {
        // Block size
        CBLAS_INT block_size = end - start + 1;

        // Compute all projection coefficients for this block: work = V_block^T * vnew
        BLAS_FUNC(sgemv)("T", &n, &block_size, &one, &V[start * ldv], &ldv, vnew, &ione, &zero, work, &ione);

        // Orthogonalize: vnew = vnew - V_block * work
        BLAS_FUNC(sgemv)("N", &n, &block_size, &neg_one, &V[start * ldv], &ldv, work, &ione, &one, vnew, &ione);

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void dcgs(CBLAS_INT n, CBLAS_INT k, double* V, CBLAS_INT ldv, double* vnew, const CBLAS_INT* indices, double* work) {
    CBLAS_INT ione = 1;
    double one = 1.0;
    double zero = 0.0;
    double neg_one = -1.0;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    CBLAS_INT idx = 0;
    CBLAS_INT start = indices[idx];
    CBLAS_INT end = indices[idx + 1];

    while (start < k && start >= 0 && start <= end) {
        // Block size
        CBLAS_INT block_size = end - start + 1;

        // Compute all projection coefficients for this block: work = V_block^T * vnew
        BLAS_FUNC(dgemv)("T", &n, &block_size, &one, &V[start * ldv], &ldv, vnew, &ione, &zero, work, &ione);

        // Orthogonalize: vnew = vnew - V_block * work
        BLAS_FUNC(dgemv)("N", &n, &block_size, &neg_one, &V[start * ldv], &ldv, work, &ione, &one, vnew, &ione);

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void ccgs(CBLAS_INT n, CBLAS_INT k, PROPACK_CPLXF_TYPE* V, CBLAS_INT ldv, PROPACK_CPLXF_TYPE* vnew, const CBLAS_INT* indices, PROPACK_CPLXF_TYPE* work) {
    CBLAS_INT ione = 1;
    PROPACK_CPLXF_TYPE one = PROPACK_cplxf(1.0f, 0.0f);
    PROPACK_CPLXF_TYPE zero = PROPACK_cplxf(0.0f, 0.0f);
    PROPACK_CPLXF_TYPE neg_one = PROPACK_cplxf(-1.0f, 0.0f);

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    CBLAS_INT idx = 0;
    CBLAS_INT start = indices[idx];
    CBLAS_INT end = indices[idx + 1];

    while (start < k && start >= 0 && start <= end) {
        // Block size
        CBLAS_INT block_size = end - start + 1;

        // Compute all projection coefficients for this block: work = V_block^H * vnew
        BLAS_FUNC(cgemv)("C", &n, &block_size, &one, &V[start * ldv], &ldv, vnew, &ione, &zero, work, &ione);

        // Orthogonalize: vnew = vnew - V_block * work
        BLAS_FUNC(cgemv)("N", &n, &block_size, &neg_one, &V[start * ldv], &ldv, work, &ione, &one, vnew, &ione);

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void zcgs(CBLAS_INT n, CBLAS_INT k, PROPACK_CPLX_TYPE* V, CBLAS_INT ldv, PROPACK_CPLX_TYPE* vnew, const CBLAS_INT* indices, PROPACK_CPLX_TYPE* work) {
    CBLAS_INT ione = 1;
    PROPACK_CPLX_TYPE one = PROPACK_cplx(1.0, 0.0);
    PROPACK_CPLX_TYPE zero = PROPACK_cplx(0.0, 0.0);
    PROPACK_CPLX_TYPE neg_one = PROPACK_cplx(-1.0, 0.0);

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    CBLAS_INT idx = 0;
    CBLAS_INT start = indices[idx];
    CBLAS_INT end = indices[idx + 1];

    while (start < k && start >= 0 && start <= end) {
        // Block size
        CBLAS_INT block_size = end - start + 1;

        // Compute all projection coefficients for this block: work = V_block^H * vnew
        BLAS_FUNC(zgemv)("C", &n, &block_size, &one, &V[start * ldv], &ldv, vnew, &ione, &zero, work, &ione);

        // Orthogonalize: vnew = vnew - V_block * work
        BLAS_FUNC(zgemv)("N", &n, &block_size, &neg_one, &V[start * ldv], &ldv, work, &ione, &one, vnew, &ione);

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void sreorth(CBLAS_INT n, CBLAS_INT k, float* V, CBLAS_INT ldv, float* vnew, float* normvnew, const CBLAS_INT* indices, float alpha, float* work, CBLAS_INT iflag) {
    const CBLAS_INT NTRY = 5;
    CBLAS_INT ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    for (CBLAS_INT itry = 0; itry < NTRY; itry++)
    {
        float normvnew_0 = *normvnew;

        if (iflag == 1) {
            scgs(n, k, V, ldv, vnew, indices, work);
        } else {
            smgs(n, k, V, ldv, vnew, indices);
        }

        *normvnew = BLAS_FUNC(snrm2)(&n, vnew, &ione);
        if (*normvnew > alpha * normvnew_0) { return; }
    }

    // vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    *normvnew = 0.0f;
    for (CBLAS_INT i = 0; i < n; i++) { vnew[i] = 0.0f; }
}


void dreorth(CBLAS_INT n, CBLAS_INT k, double* V, CBLAS_INT ldv, double* vnew, double* normvnew, const CBLAS_INT* indices, double alpha, double* work, CBLAS_INT iflag) {
    const CBLAS_INT NTRY = 5;
    CBLAS_INT ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    for (CBLAS_INT itry = 0; itry < NTRY; itry++)
    {
        double normvnew_0 = *normvnew;

        if (iflag == 1) {
            dcgs(n, k, V, ldv, vnew, indices, work);
        } else {
            dmgs(n, k, V, ldv, vnew, indices);
        }

        *normvnew = BLAS_FUNC(dnrm2)(&n, vnew, &ione);

        if (*normvnew > alpha * normvnew_0) { return; }
    }

    // vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    *normvnew = 0.0;
    for (CBLAS_INT i = 0; i < n; i++) { vnew[i] = 0.0; }
}


void creorth(CBLAS_INT n, CBLAS_INT k, PROPACK_CPLXF_TYPE* V, CBLAS_INT ldv, PROPACK_CPLXF_TYPE* vnew, float* normvnew, const CBLAS_INT* indices, float alpha, PROPACK_CPLXF_TYPE* work, CBLAS_INT iflag) {
    const CBLAS_INT NTRY = 5;
    CBLAS_INT ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    for (CBLAS_INT itry = 0; itry < NTRY; itry++)
    {
        float normvnew_0 = *normvnew;

        if (iflag == 1) {
            ccgs(n, k, V, ldv, vnew, indices, work);
        } else {
            cmgs(n, k, V, ldv, vnew, indices);
        }

        *normvnew = BLAS_FUNC(scnrm2)(&n, vnew, &ione);

        if (*normvnew > alpha * normvnew_0) { return; }
    }

    // vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    *normvnew = 0.0f;
    for (CBLAS_INT i = 0; i < n; i++) { vnew[i] = PROPACK_cplxf(0.0f, 0.0f); }
}


void zreorth(CBLAS_INT n, CBLAS_INT k, PROPACK_CPLX_TYPE* V, CBLAS_INT ldv, PROPACK_CPLX_TYPE* vnew, double* normvnew, const CBLAS_INT* indices, double alpha, PROPACK_CPLX_TYPE* work, CBLAS_INT iflag) {
    const CBLAS_INT NTRY = 5;
    CBLAS_INT ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    for (CBLAS_INT itry = 0; itry < NTRY; itry++)
    {
        double normvnew_0 = *normvnew;

        if (iflag == 1) {
            zcgs(n, k, V, ldv, vnew, indices, work);
        } else {
            zmgs(n, k, V, ldv, vnew, indices);
        }

        *normvnew = BLAS_FUNC(dznrm2)(&n, vnew, &ione);

        if (*normvnew > alpha * normvnew_0) { return; }
    }

    // vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    *normvnew = 0.0;
    for (CBLAS_INT i = 0; i < n; i++) { vnew[i] = PROPACK_cplx(0.0, 0.0); }
}
