#include "getu0.h"


void sgetu0(int transa, CBLAS_INT m, CBLAS_INT n, CBLAS_INT j, CBLAS_INT ntry, float* u0, float* u0norm, float* U, CBLAS_INT ldu,
           PROPACK_aprod_s aprod, float* dparm, CBLAS_INT* iparm, CBLAS_INT* ierr, CBLAS_INT icgs, float* anormest, float* work, uint64_t* rng_state) {

    const float kappa = sqrtf(2.0f) / 2.0f;
    CBLAS_INT ione = 1;
    CBLAS_INT rsize, usize;

    // Determine vector sizes based on transpose flag
    if (transa == 0) {
        rsize = n;  // Random vector size
        usize = m;  // Output vector size
    } else {
        rsize = m;  // Random vector size
        usize = n;  // Output vector size
    }

    *ierr = 0;

    for (CBLAS_INT itry = 0; itry < ntry; itry++) {
        // Generate random vector
        for (CBLAS_INT i = 0; i < rsize; i++) { work[i] = random_float(rng_state); }

        // Compute norm of random vector
        float nrm = BLAS_FUNC(snrm2)(&rsize, work, &ione);

        // Apply matrix operation: u0 = Op(A) * work
        aprod(transa, m, n, work, u0, dparm, iparm);

        // Compute norm of result and estimate operator norm
        *u0norm = BLAS_FUNC(snrm2)(&usize, u0, &ione);
        *anormest = *u0norm / nrm;

        // Orthogonalize against existing vectors if j >= 1
        if (j >= 1) {
            CBLAS_INT indices[4];
            indices[0] = 0;      // start index (0-based)
            indices[1] = j;  // end index (0-based, inclusive)
            indices[2] = j+1;      // terminator
            indices[3] = j+1;      // Extra termination value to prevent out-of-bounds access.

            sreorth(usize, j, U, ldu, u0, u0norm, indices, kappa, work, icgs);
        }

        // Check if we have a valid vector
        if (*u0norm > 0.0f) { return; }
    }

    // Failed to generate valid vector
    *ierr = -1;
}


void dgetu0(int transa, CBLAS_INT m, CBLAS_INT n, CBLAS_INT j, CBLAS_INT ntry, double* u0, double* u0norm, double* U, CBLAS_INT ldu,
           PROPACK_aprod_d aprod, double* dparm, CBLAS_INT* iparm, CBLAS_INT* ierr, CBLAS_INT icgs, double* anormest, double* work, uint64_t* rng_state) {

    const double kappa = sqrt(2.0) / 2.0;
    CBLAS_INT ione = 1;
    CBLAS_INT rsize, usize;

    // Determine vector sizes based on transpose flag
    if (transa == 0) {
        rsize = n;  // Random vector size
        usize = m;  // Output vector size
    } else {
        rsize = m;  // Random vector size
        usize = n;  // Output vector size
    }

    *ierr = 0;

    for (CBLAS_INT itry = 0; itry < ntry; itry++) {
        // Generate random vector
        for (CBLAS_INT i = 0; i < rsize; i++) { work[i] = random_double(rng_state); }

        // Compute norm of random vector
        double nrm = BLAS_FUNC(dnrm2)(&rsize, work, &ione);

        // Apply matrix operation: u0 = Op(A) * work
        aprod(transa, m, n, work, u0, dparm, iparm);

        // Compute norm of result and estimate operator norm
        *u0norm = BLAS_FUNC(dnrm2)(&usize, u0, &ione);
        *anormest = *u0norm / nrm;

        // Orthogonalize against existing vectors if j >= 1
        if (j >= 1) {
            CBLAS_INT indices[4];
            indices[0] = 0;      // start index (0-based)
            indices[1] = j;      // end index (0-based, inclusive)
            indices[2] = j + 1;  // terminator
            indices[3] = j + 1;  // Extra termination value to prevent out-of-bounds access.

            dreorth(usize, j, U, ldu, u0, u0norm, indices, kappa, work, icgs);
        }

        // Check if we have a valid vector
        if (*u0norm > 0.0) { return; }
    }

    // Failed to generate valid vector
    *ierr = -1;
}


void cgetu0(int transa, CBLAS_INT m, CBLAS_INT n, CBLAS_INT j, CBLAS_INT ntry, PROPACK_CPLXF_TYPE* u0, float* u0norm, PROPACK_CPLXF_TYPE* U, CBLAS_INT ldu,
           PROPACK_aprod_c aprod, PROPACK_CPLXF_TYPE* cparm, CBLAS_INT* iparm, CBLAS_INT* ierr, CBLAS_INT icgs, float* anormest, PROPACK_CPLXF_TYPE* work,
           uint64_t* rng_state)
{
    const float kappa = sqrtf(2.0f) / 2.0f;
    CBLAS_INT ione = 1;
    CBLAS_INT rsize, usize;

    // Determine vector sizes based on transpose flag
    if (transa == 0) {
        rsize = n;  // Random vector size
        usize = m;  // Output vector size
    } else {
        rsize = m;  // Random vector size
        usize = n;  // Output vector size
    }

    *ierr = 0;

    for (CBLAS_INT itry = 0; itry < ntry; itry++) {
        // Generate random complex vector
        for (CBLAS_INT i = 0; i < rsize; i++) {
            work[i] = PROPACK_cplxf(random_float(rng_state), random_float(rng_state));
        }

        // Compute norm of random vector
        float nrm = BLAS_FUNC(scnrm2)(&rsize, work, &ione);

        // Apply matrix operation: u0 = Op(A) * work
        aprod(transa, m, n, work, u0, cparm, iparm);

        // Compute norm of result and estimate operator norm
        *u0norm = BLAS_FUNC(scnrm2)(&usize, u0, &ione);
        *anormest = *u0norm / nrm;

        // Orthogonalize against existing vectors if j >= 1
        if (j >= 1) {
            CBLAS_INT indices[4];
            indices[0] = 0;      // start index (0-based)
            indices[1] = j;      // end index (0-based, inclusive)
            indices[2] = j + 1;  // terminator
            indices[3] = j + 1;  // Extra termination value to prevent out-of-bounds access.

            creorth(usize, j, U, ldu, u0, u0norm, indices, kappa, work, icgs);
        }

        // Check if we have a valid vector
        if (*u0norm > 0.0f) { return; }
    }

    // Failed to generate valid vector
    *ierr = -1;
}


void zgetu0(int transa, CBLAS_INT m, CBLAS_INT n, CBLAS_INT j, CBLAS_INT ntry, PROPACK_CPLX_TYPE* u0, double* u0norm, PROPACK_CPLX_TYPE* U, CBLAS_INT ldu,
           PROPACK_aprod_z aprod, PROPACK_CPLX_TYPE* zparm, CBLAS_INT* iparm, CBLAS_INT* ierr, CBLAS_INT icgs, double* anormest, PROPACK_CPLX_TYPE* work,
           uint64_t* rng_state)
{
    const double kappa = sqrt(2.0) / 2.0;
    CBLAS_INT ione = 1;
    CBLAS_INT rsize, usize;

    // Determine vector sizes based on transpose flag
    if (transa == 0) {
        rsize = n;  // Random vector size
        usize = m;  // Output vector size
    } else {
        rsize = m;  // Random vector size
        usize = n;  // Output vector size
    }

    *ierr = 0;

    for (CBLAS_INT itry = 0; itry < ntry; itry++) {
        // Generate random complex vector
        for (CBLAS_INT i = 0; i < rsize; i++) {
            work[i] = PROPACK_cplx(random_double(rng_state), random_double(rng_state));
        }

        // Compute norm of random vector
        double nrm = BLAS_FUNC(dznrm2)(&rsize, work, &ione);

        // Apply matrix operation: u0 = Op(A) * work
        aprod(transa, m, n, work, u0, zparm, iparm);

        // Compute norm of result and estimate operator norm
        *u0norm = BLAS_FUNC(dznrm2)(&usize, u0, &ione);
        *anormest = *u0norm / nrm;

        // Orthogonalize against existing vectors if j >= 1
        if (j >= 1) {
            CBLAS_INT indices[4];
            indices[0] = 0;      // start index (0-based)
            indices[1] = j;      // end index (0-based, inclusive)
            indices[2] = j + 1;  // terminator
            indices[3] = j + 1;  // Extra termination value to prevent out-of-bounds access.

            zreorth(usize, j, U, ldu, u0, u0norm, indices, kappa, work, icgs);
        }

        // Check if we have a valid vector
        if (*u0norm > 0.0) { return; }
    }

    // Failed to generate valid vector
    *ierr = -1;
}
