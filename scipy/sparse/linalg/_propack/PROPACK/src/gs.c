#include "gs.h"


void smgs(int n, int k, float* V, int ldv, float* vnew, const int* indices) {
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    /**
     * PROPACK encodes sentinels in 1-index specific way.
     * Therefore, we have to guard for a few edge cases of 0-indexing.
     * TODO: Fix this by rewriting the indexing logic.
     */
    while ((start <= k) && ((start < end) || ((start == 0) && (end == 0) && (idx == 0)))) {
        // Orthogonalize against columns [start, end] (0-indexed)
        for (int i = start; i <= end; i++) {
            // Compute projection coefficient: coef = V(:,i)' * vnew
            float coef = sdot_(&n, &V[i * ldv], &ione, vnew, &ione);

            // Orthogonalize: vnew = vnew - coef * V(:,i)
            float neg_coef = -coef;
            saxpy_(&n, &neg_coef, &V[i * ldv], &ione, vnew, &ione);
        }

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void dmgs(int n, int k, double* V, int ldv, double* vnew, const int* indices) {
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    /**
     * PROPACK encodes sentinels in 1-index specific way.
     * Therefore, we have to guard for a few edge cases of 0-indexing.
     * TODO: Fix this by rewriting the indexing logic.
     */
    while ((start <= k) && ((start < end) || ((start == 0) && (end == 0) && (idx == 0)))) {
        // Orthogonalize against columns [start, end] (0-indexed)
        for (int i = start; i <= end; i++) {
            // Compute projection coefficient: coef = V(:,i)' * vnew
            double coef = ddot_(&n, &V[i * ldv], &ione, vnew, &ione);

            // Orthogonalize: vnew = vnew - coef * V(:,i)
            double neg_coef = -coef;
            daxpy_(&n, &neg_coef, &V[i * ldv], &ione, vnew, &ione);
        }

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void cmgs(int n, int k, PROPACK_CPLXF_TYPE* V, int ldv, PROPACK_CPLXF_TYPE* vnew, const int* indices) {
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    while ((start <= k) && ((start < end) || ((start == 0) && (end == 0) && (idx == 0)))) {
        // Orthogonalize against columns [start, end] (0-indexed)
        for (int i = start; i <= end; i++) {
            // Compute projection coefficient: coef = V(:,i)^H * vnew (conjugate dot product)
            PROPACK_CPLXF_TYPE coef = cdotc_(&n, &V[i * ldv], &ione, vnew, &ione);

            // Orthogonalize: vnew = vnew - coef * V(:,i)
            PROPACK_CPLXF_TYPE neg_coef = PROPACK_cplxf(-crealf(coef), -cimagf(coef));
            caxpy_(&n, &neg_coef, &V[i * ldv], &ione, vnew, &ione);
        }

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void zmgs(int n, int k, PROPACK_CPLX_TYPE* V, int ldv, PROPACK_CPLX_TYPE* vnew, const int* indices) {
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    while ((start <= k) && ((start < end) || ((start == 0) && (end == 0) && (idx == 0)))) {
        // Orthogonalize against columns [start, end] (0-indexed)
        for (int i = start; i <= end; i++) {
            // Compute projection coefficient: coef = V(:,i)^H * vnew (conjugate dot product)
            PROPACK_CPLX_TYPE coef = zdotc_(&n, &V[i * ldv], &ione, vnew, &ione);

            // Orthogonalize: vnew = vnew - coef * V(:,i)
            PROPACK_CPLX_TYPE neg_coef = PROPACK_cplx(-creal(coef), -cimag(coef));
            zaxpy_(&n, &neg_coef, &V[i * ldv], &ione, vnew, &ione);
        }

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void scgs(int n, int k, float* V, int ldv, float* vnew, const int* indices, float* work) {
    int ione = 1;
    float one = 1.0f;
    float zero = 0.0f;
    float neg_one = -1.0f;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    while (start < k && start >= 0 && start <= end) {
        // Block size
        int block_size = end - start + 1;

        // Compute all projection coefficients for this block: work = V_block^T * vnew
        sgemv_("T", &n, &block_size, &one, &V[start * ldv], &ldv, vnew, &ione, &zero, work, &ione);

        // Orthogonalize: vnew = vnew - V_block * work
        sgemv_("N", &n, &block_size, &neg_one, &V[start * ldv], &ldv, work, &ione, &one, vnew, &ione);

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void dcgs(int n, int k, double* V, int ldv, double* vnew, const int* indices, double* work) {
    int ione = 1;
    double one = 1.0;
    double zero = 0.0;
    double neg_one = -1.0;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    while (start < k && start >= 0 && start <= end) {
        // Block size
        int block_size = end - start + 1;

        // Compute all projection coefficients for this block: work = V_block^T * vnew
        dgemv_("T", &n, &block_size, &one, &V[start * ldv], &ldv, vnew, &ione, &zero, work, &ione);

        // Orthogonalize: vnew = vnew - V_block * work
        dgemv_("N", &n, &block_size, &neg_one, &V[start * ldv], &ldv, work, &ione, &one, vnew, &ione);

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void ccgs(int n, int k, PROPACK_CPLXF_TYPE* V, int ldv, PROPACK_CPLXF_TYPE* vnew, const int* indices, PROPACK_CPLXF_TYPE* work) {
    int ione = 1;
    PROPACK_CPLXF_TYPE one = PROPACK_cplxf(1.0f, 0.0f);
    PROPACK_CPLXF_TYPE zero = PROPACK_cplxf(0.0f, 0.0f);
    PROPACK_CPLXF_TYPE neg_one = PROPACK_cplxf(-1.0f, 0.0f);

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    while (start < k && start >= 0 && start <= end) {
        // Block size
        int block_size = end - start + 1;

        // Compute all projection coefficients for this block: work = V_block^H * vnew
        cgemv_("C", &n, &block_size, &one, &V[start * ldv], &ldv, vnew, &ione, &zero, work, &ione);

        // Orthogonalize: vnew = vnew - V_block * work
        cgemv_("N", &n, &block_size, &neg_one, &V[start * ldv], &ldv, work, &ione, &one, vnew, &ione);

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void zcgs(int n, int k, PROPACK_CPLX_TYPE* V, int ldv, PROPACK_CPLX_TYPE* vnew, const int* indices, PROPACK_CPLX_TYPE* work) {
    int ione = 1;
    PROPACK_CPLX_TYPE one = PROPACK_cplx(1.0, 0.0);
    PROPACK_CPLX_TYPE zero = PROPACK_cplx(0.0, 0.0);
    PROPACK_CPLX_TYPE neg_one = PROPACK_cplx(-1.0, 0.0);

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    // Process each block specified in indices array
    int idx = 0;
    int start = indices[idx];
    int end = indices[idx + 1];

    while (start < k && start >= 0 && start <= end) {
        // Block size
        int block_size = end - start + 1;

        // Compute all projection coefficients for this block: work = V_block^H * vnew
        zgemv_("C", &n, &block_size, &one, &V[start * ldv], &ldv, vnew, &ione, &zero, work, &ione);

        // Orthogonalize: vnew = vnew - V_block * work
        zgemv_("N", &n, &block_size, &neg_one, &V[start * ldv], &ldv, work, &ione, &one, vnew, &ione);

        idx += 2;  // Move to next block
        start = indices[idx];
        end = indices[idx + 1];
    }
}


void sreorth(int n, int k, float* V, int ldv, float* vnew, float* normvnew, const int* indices, float alpha, float* work, int iflag) {
    const int NTRY = 5;
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    for (int itry = 0; itry < NTRY; itry++)
    {
        float normvnew_0 = *normvnew;

        if (iflag == 1) {
            scgs(n, k, V, ldv, vnew, indices, work);
        } else {
            smgs(n, k, V, ldv, vnew, indices);
        }

        *normvnew = snrm2_(&n, vnew, &ione);
        if (*normvnew > alpha * normvnew_0) { return; }
    }

    // vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    *normvnew = 0.0f;
    for (int i = 0; i < n; i++) { vnew[i] = 0.0f; }
}


void dreorth(int n, int k, double* V, int ldv, double* vnew, double* normvnew, const int* indices, double alpha, double* work, int iflag) {
    const int NTRY = 5;
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    for (int itry = 0; itry < NTRY; itry++)
    {
        double normvnew_0 = *normvnew;

        if (iflag == 1) {
            dcgs(n, k, V, ldv, vnew, indices, work);
        } else {
            dmgs(n, k, V, ldv, vnew, indices);
        }

        *normvnew = dnrm2_(&n, vnew, &ione);

        if (*normvnew > alpha * normvnew_0) { return; }
    }

    // vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    *normvnew = 0.0;
    for (int i = 0; i < n; i++) { vnew[i] = 0.0; }
}


void creorth(int n, int k, PROPACK_CPLXF_TYPE* V, int ldv, PROPACK_CPLXF_TYPE* vnew, float* normvnew, const int* indices, float alpha, PROPACK_CPLXF_TYPE* work, int iflag) {
    const int NTRY = 5;
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    for (int itry = 0; itry < NTRY; itry++)
    {
        float normvnew_0 = *normvnew;

        if (iflag == 1) {
            ccgs(n, k, V, ldv, vnew, indices, work);
        } else {
            cmgs(n, k, V, ldv, vnew, indices);
        }

        *normvnew = scnrm2_(&n, vnew, &ione);

        if (*normvnew > alpha * normvnew_0) { return; }
    }

    // vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    *normvnew = 0.0f;
    for (int i = 0; i < n; i++) { vnew[i] = PROPACK_cplxf(0.0f, 0.0f); }
}


void zreorth(int n, int k, PROPACK_CPLX_TYPE* V, int ldv, PROPACK_CPLX_TYPE* vnew, double* normvnew, const int* indices, double alpha, PROPACK_CPLX_TYPE* work, int iflag) {
    const int NTRY = 5;
    int ione = 1;

    // Check for quick return
    if ((k < 0) || (n <= 0)) { return; }

    for (int itry = 0; itry < NTRY; itry++)
    {
        double normvnew_0 = *normvnew;

        if (iflag == 1) {
            zcgs(n, k, V, ldv, vnew, indices, work);
        } else {
            zmgs(n, k, V, ldv, vnew, indices);
        }

        *normvnew = dznrm2_(&n, vnew, &ione);

        if (*normvnew > alpha * normvnew_0) { return; }
    }

    // vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    *normvnew = 0.0;
    for (int i = 0; i < n; i++) { vnew[i] = PROPACK_cplx(0.0, 0.0); }
}
