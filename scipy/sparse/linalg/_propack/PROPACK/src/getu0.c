#include "getu0.h"


void sgetu0(int transa, int m, int n, int j, int ntry, float* u0, float* u0norm, float* U, int ldu,
           PROPACK_aprod_s aprod, float* dparm, int* iparm, int* ierr, int icgs, float* anormest, float* work, uint64_t* rng_state) {

    const float kappa = sqrtf(2.0f) / 2.0f;
    int ione = 1;
    int rsize, usize;

    // Determine vector sizes based on transpose flag
    if (transa == 0) {
        rsize = n;  // Random vector size
        usize = m;  // Output vector size
    } else {
        rsize = m;  // Random vector size
        usize = n;  // Output vector size
    }

    *ierr = 0;

    for (int itry = 0; itry < ntry; itry++) {
        // Generate random vector
        for (int i = 0; i < rsize; i++) { work[i] = random_float(rng_state); }

        // Compute norm of random vector
        float nrm = snrm2_(&rsize, work, &ione);

        // Apply matrix operation: u0 = Op(A) * work
        aprod(transa, m, n, work, u0, dparm, iparm);

        // Compute norm of result and estimate operator norm
        *u0norm = snrm2_(&usize, u0, &ione);
        *anormest = *u0norm / nrm;

        // Orthogonalize against existing vectors if j >= 1
        if (j >= 1) {
            int indices[4];
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


void dgetu0(int transa, int m, int n, int j, int ntry, double* u0, double* u0norm, double* U, int ldu,
           PROPACK_aprod_d aprod, double* dparm, int* iparm, int* ierr, int icgs, double* anormest, double* work, uint64_t* rng_state) {

    const double kappa = sqrt(2.0) / 2.0;
    int ione = 1;
    int rsize, usize;

    // Determine vector sizes based on transpose flag
    if (transa == 0) {
        rsize = n;  // Random vector size
        usize = m;  // Output vector size
    } else {
        rsize = m;  // Random vector size
        usize = n;  // Output vector size
    }

    *ierr = 0;

    for (int itry = 0; itry < ntry; itry++) {
        // Generate random vector
        for (int i = 0; i < rsize; i++) { work[i] = random_double(rng_state); }

        // Compute norm of random vector
        double nrm = dnrm2_(&rsize, work, &ione);

        // Apply matrix operation: u0 = Op(A) * work
        aprod(transa, m, n, work, u0, dparm, iparm);

        // Compute norm of result and estimate operator norm
        *u0norm = dnrm2_(&usize, u0, &ione);
        *anormest = *u0norm / nrm;

        // Orthogonalize against existing vectors if j >= 1
        if (j >= 1) {
            int indices[4];
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


void cgetu0(int transa, int m, int n, int j, int ntry, PROPACK_CPLXF_TYPE* u0, float* u0norm, PROPACK_CPLXF_TYPE* U, int ldu,
           PROPACK_aprod_c aprod, PROPACK_CPLXF_TYPE* cparm, int* iparm, int* ierr, int icgs, float* anormest, PROPACK_CPLXF_TYPE* work,
           uint64_t* rng_state)
{
    const float kappa = sqrtf(2.0f) / 2.0f;
    int ione = 1;
    int rsize, usize;

    // Determine vector sizes based on transpose flag
    if (transa == 0) {
        rsize = n;  // Random vector size
        usize = m;  // Output vector size
    } else {
        rsize = m;  // Random vector size
        usize = n;  // Output vector size
    }

    *ierr = 0;

    for (int itry = 0; itry < ntry; itry++) {
        // Generate random complex vector
        for (int i = 0; i < rsize; i++) {
            work[i] = PROPACK_cplxf(random_float(rng_state), random_float(rng_state));
        }

        // Compute norm of random vector
        float nrm = scnrm2_(&rsize, work, &ione);

        // Apply matrix operation: u0 = Op(A) * work
        aprod(transa, m, n, work, u0, cparm, iparm);

        // Compute norm of result and estimate operator norm
        *u0norm = scnrm2_(&usize, u0, &ione);
        *anormest = *u0norm / nrm;

        // Orthogonalize against existing vectors if j >= 1
        if (j >= 1) {
            int indices[4];
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


void zgetu0(int transa, int m, int n, int j, int ntry, PROPACK_CPLX_TYPE* u0, double* u0norm, PROPACK_CPLX_TYPE* U, int ldu,
           PROPACK_aprod_z aprod, PROPACK_CPLX_TYPE* zparm, int* iparm, int* ierr, int icgs, double* anormest, PROPACK_CPLX_TYPE* work,
           uint64_t* rng_state)
{
    const double kappa = sqrt(2.0) / 2.0;
    int ione = 1;
    int rsize, usize;

    // Determine vector sizes based on transpose flag
    if (transa == 0) {
        rsize = n;  // Random vector size
        usize = m;  // Output vector size
    } else {
        rsize = m;  // Random vector size
        usize = n;  // Output vector size
    }

    *ierr = 0;

    for (int itry = 0; itry < ntry; itry++) {
        // Generate random complex vector
        for (int i = 0; i < rsize; i++) {
            work[i] = PROPACK_cplx(random_double(rng_state), random_double(rng_state));
        }

        // Compute norm of random vector
        double nrm = dznrm2_(&rsize, work, &ione);

        // Apply matrix operation: u0 = Op(A) * work
        aprod(transa, m, n, work, u0, zparm, iparm);

        // Compute norm of result and estimate operator norm
        *u0norm = dznrm2_(&usize, u0, &ione);
        *anormest = *u0norm / nrm;

        // Orthogonalize against existing vectors if j >= 1
        if (j >= 1) {
            int indices[4];
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
