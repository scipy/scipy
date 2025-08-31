#include <math.h>
#include <float.h>
#include "lanbpro.h"


static inline int int_min(const int a, const int b) { return a < b ? a : b; }
static inline int int_max(const int a, const int b) { return a > b ? a : b; }


/**
 * @brief Scale the vector x by 1/alpha avoiding unnecessary under- and overflow.
 *
 * @param n The number of elements in the vector x.
 * @param alpha The scaling factor.
 * @param x The vector to be scaled.
 */
static void ssafescal(int n, float alpha, float* x)
{
    // Scale the vector x by 1/alpha avoiding unnecessary under- and overflow.
    static float sfmin = FLT_MIN;
    int ione = 1;
    int info;

    if (fabsf(alpha) >= sfmin) {
        float inv_alpha = 1.0f / alpha;
        sscal_(&n, &inv_alpha, x, &ione);
    } else {
        // Use LAPACK's safe scaling for very small alpha values
        float one = 1.0f;
        slascl_("G", &ione, &ione, &alpha, &one, &n, &ione, x, &n, &info);
    }
}


/**
 * @brief Scale the vector x by 1/alpha avoiding unnecessary under- and overflow.
 *
 * @param n The number of elements in the vector x.
 * @param alpha The scaling factor.
 * @param x The vector to be scaled.
 */
void dsafescal(int n, double alpha, double* x)
{
    // Scale the vector x by 1/alpha avoiding unnecessary under- and overflow.
    static double sfmin = DBL_MIN;
    int ione = 1;
    int info;

    if (fabs(alpha) >= sfmin) {
        double inv_alpha = 1.0 / alpha;
        dscal_(&n, &inv_alpha, x, &ione);
    } else {
        // Use LAPACK's safe scaling for very small alpha values
        double one = 1.0;
        dlascl_("G", &ione, &ione, &alpha, &one, &n, &ione, x, &n, &info);
    }
}


/**
 * @brief Scale the complex vector x by 1/alpha avoiding unnecessary under- and overflow.
 *
 * @param n The number of elements in the vector x.
 * @param alpha The scaling factor.
 * @param x The complex vector to be scaled.
 */
void csafescal(int n, float alpha, PROPACK_CPLXF_TYPE* x)
{
    // Scale the complex vector x by 1/alpha avoiding unnecessary under- and overflow.
    static float sfmin = FLT_MIN;
    int ione = 1;
    int info;

    if (fabsf(alpha) >= sfmin) {
        float inv_alpha = 1.0f / alpha;
        csscal_(&n, &inv_alpha, x, &ione);
    } else {
        // Use LAPACK's safe scaling for very small alpha values
        float one = 1.0f;
        clascl_("G", &ione, &ione, &alpha, &one, &n, &ione, x, &n, &info);
    }
}


/**
 * @brief Scale the complex vector x by 1/alpha avoiding unnecessary under- and overflow.
 *
 * @param n The number of elements in the vector x.
 * @param alpha The scaling factor.
 * @param x The complex vector to be scaled.
 */
void zsafescal(int n, double alpha, PROPACK_CPLX_TYPE* x)
{
    // Scale the complex vector x by 1/alpha avoiding unnecessary under- and overflow.
    static double sfmin = DBL_MIN;
    int ione = 1;
    int info;

    if (fabs(alpha) >= sfmin) {
        double inv_alpha = 1.0 / alpha;
        zdscal_(&n, &inv_alpha, x, &ione);
    } else {
        // Use LAPACK's safe scaling for very small alpha values
        double one = 1.0;
        zlascl_("G", &ione, &ione, &alpha, &one, &n, &ione, x, &n, &info);
    }
}


/**
 * @brief Perform Lanczos bidiagonalization with partial reorthogonalization (single precision).
 *
 * This function computes a partial Lanczos bidiagonalization of a matrix A using
 * single precision arithmetic. The bidiagonalization is performed with partial
 * reorthogonalization to maintain numerical stability and efficiency.
 *
 * @param m         The number of rows in the matrix A.
 * @param n         The number of columns in the matrix A.
 * @param k0        The starting index for the bidiagonalization process.
 * @param k         Pointer to the number of Lanczos steps to perform.
 * @param aprod     Function pointer for matrix-vector multiplication with A or A^T.
 * @param U         Pointer to the left singular vectors (output).
 * @param ldu       The leading dimension of the array U.
 * @param V         Pointer to the right singular vectors (output).
 * @param ldv       The leading dimension of the array V.
 * @param B         Pointer to the bidiagonal matrix (output).
 * @param ldb       The leading dimension of the array B.
 * @param rnorm     Pointer to the norm of the residual vector (output).
 * @param doption   Array of double precision options for controlling the algorithm.
 * @param ioption   Array of integer options for controlling the algorithm.
 * @param work      Workspace array for intermediate computations.
 * @param iwork     Integer workspace array for intermediate computations.
 * @param dparm     Array of double precision parameters for the algorithm.
 * @param iparm     Array of integer parameters for the algorithm.
 * @param ierr      Pointer to the error flag (output).
 * @param rng_state Pointer to the random number generator state.
 */
void slanbpro(int m, int n, int k0, int* k, PROPACK_aprod_s aprod,
              float* U, int ldu, float* V, int ldv, float* B, int ldb,
              float* rnorm, float* doption, int* ioption, float* work, int* iwork,
              float* dparm, int* iparm, int* ierr, uint64_t* rng_state)
{

    // Constants
    const float FUDGE = 1.01f;
    const float kappa = sqrtf(2.0f) / 2.0f;

    // Machine constants
    const float eps = FLT_EPSILON*0.5f;
    const float eps34 = powf(eps, 0.75f);
    const float epsn = (float)int_max(m, n) * eps;
    const float epsn2 = sqrtf((float)int_max(m, n)) * eps;

    // Local variables
    int i, j, inu, imu, is, iidx, j0;
    float s, mumax, numax, alpha, beta, a1, b1, amax, anormest;
    int force_reorth, full_reorth;
    int ione = 1;

    // Set default parameters
    float delta, eta, anorm;
    int cgs, elr;

    if (doption[0] < 0.0f)
    {
        delta = sqrtf(eps / *k);
    } else {
        delta = doption[0];
    }

    if (doption[1] < 0.0f)
    {
        eta = eps34 / sqrtf((float)*k);
    } else {
        eta = doption[1];
    }

    if ((delta <= eta) || (delta == 0.0f))
    {
        full_reorth = 1;
    } else {
        full_reorth = 0;
    }

    if (doption[2] > 0.0f)
    {
        anorm = doption[2];
    } else if (k0 > 0) {
        anorm = hypotf(B[0], B[ldb]);
        if (anorm <= 0.0f) {
            *ierr = -1;
            doption[2] = anorm;
            return;
        }
    } else {
        anorm = 0.0f;
    }

    cgs = ioption[0];
    elr = ioption[1];
    *ierr = 0;

    // Get starting vector if needed
    if (*rnorm == 0.0f) {
        sgetu0(0, m, n, k0, 3, &U[k0 * ldu], rnorm, U, ldu, aprod, dparm, iparm, ierr, cgs, &anormest, work, rng_state);
        anorm = fmaxf(anorm, anormest);
    }

    // Set pointers into work array
    // work layout: [mu(k+1), nu(k+1), workspace...]
    imu = 0;
    inu = *k + 1;
    is = 2 * (*k + 1);
    iidx = 0;
    for (i = 0; i < int_max(m, n) + 2 * (*k) + 2; i++) { work[i] = 0.0f; }
    for (i = 0; i < 2 * (*k) + 1; i++) { iwork[i] = 0; }

    // Prepare Lanczos iteration
    if (k0 == 0) {
        amax = 0.0f;
        alpha = 0.0f;
        beta = *rnorm;
        force_reorth = 0;

        // Compute ||A x|| / ||x|| for a random vector x
        // to make it less likely that ||A|| is grossly
        // underestimated at the beginning of the iteration.
        if (n > m) {
            sgetu0(0, m, n, 0, 1, &work[is], &s, U, ldu, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is + m], rng_state);
        } else {
            sgetu0(1, m, n, 0, 1, &work[is], &s, V, ldv, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is + n], rng_state);
        }
        anorm = fmaxf(anorm, FUDGE * anormest);
        j0 = 0;

        if (beta != 0.0f) { ssafescal(m, beta, &U[0]); }
        work[imu] = 1.0f;
        work[inu] = 1.0f;
    } else {
        force_reorth = 1;
        alpha = B[k0 - 1];
        beta = *rnorm;

        if ((k0 < *k) && (beta * delta < anorm * eps)) {
            full_reorth = 1;
            *ierr = k0;
        }

        iwork[iidx] = 0;
        iwork[iidx + 1] = k0 - 1;
        iwork[iidx + 2] = k0;
        iwork[iidx + 3] = k0;

        sscal_(&m, rnorm, &U[k0 * ldu], &ione);
        sreorth(m, k0, U, ldu, &U[k0 * ldu], rnorm, &iwork[iidx], kappa, &work[is], cgs);
        ssafescal(m, *rnorm, &U[k0 * ldu]);
        sset_mu(k0, &work[imu], &iwork[iidx], epsn2);
        sset_mu(k0, &work[inu], &iwork[iidx], epsn2);
        beta = *rnorm;

        // Estimate ||B||_2^2 as ||B^T * B||_1
        B[k0 - 1 + ldb] = beta;
        amax = 0.0f;
        for (j = 0; j < k0; j++)
        {
            amax = fmaxf(amax, fmaxf(B[j], B[j + ldb]));
            if (j == 0)
            {
                anorm = fmaxf(anorm, FUDGE * alpha);
            } else if (j == 1) {
                a1 = B[ldb] / amax;
                a1 = FUDGE * amax * sqrtf(powf(B[0] / amax, 2) + a1 * a1 + B[1] / amax * a1);
                anorm = fmaxf(anorm, a1);
            } else {
                a1 = B[j - 1] / amax;
                b1 = B[j - 1 + ldb] / amax;
                a1 = FUDGE * amax * sqrtf(a1 * a1 + b1 * b1 + a1 * B[j - 2 + ldb] / amax + B[j] / amax * b1);
                anorm = fmaxf(anorm, a1);
            }
        }
        j0 = k0;
    }

    numax = 0.0f;
    mumax = 0.0f;

    // Start Lanczos bidiagonalization iteration
    for (j = j0; j < *k; j++)
    {
        // alpha_j * v_j = A^T * u_j - beta_j * v_{j-1}
        aprod(1, m, n, &U[j * ldu], &V[j * ldv], dparm, iparm);

        if (j == 0)
        {
            alpha = snrm2_(&n, &V[j * ldv], &ione);
            anorm = fmaxf(anorm, FUDGE * alpha);
        } else {
            float neg_beta = -beta;
            saxpy_(&n, &neg_beta, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
            alpha = snrm2_(&n, &V[j * ldv], &ione);

            // Extended local reorthogonalization
            if ((j > 0) && (elr > 0) && (alpha < kappa * beta))
            {
                for (i = 0; i < elr; i++)
                {
                    s = sdot_(&n, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
                    float neg_s = -s;
                    saxpy_(&n, &neg_s, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
                    if (beta != 0.0f)
                    {
                        beta = beta + s;
                        B[j - 1 + ldb] = beta;
                    }
                    s = snrm2_(&n, &V[j * ldv], &ione);
                    if (s >= kappa * alpha) { break; }
                    alpha = s;
                }
                work[inu + j - 1] = eps;
                alpha = s;
            }

            B[j] = alpha;
            amax = fmaxf(amax, alpha);

            // Update estimate of ||A||_2
            if (j == 1)
            {
                a1 = B[ldb] / amax;
                a1 = FUDGE * amax * sqrtf(powf(B[0] / amax, 2) + a1 * a1 + B[1] / amax * a1);
            } else {
                a1 = B[j - 1] / amax;
                b1 = B[j - 1 + ldb] / amax;
                a1 = FUDGE * amax * sqrtf(a1 * a1 + b1 * b1 + a1 * B[j - 2 + ldb] / amax + B[j] / amax * b1);
            }
            anorm = fmaxf(anorm, a1);
        }

        // Update the nu recurrence
        if ((!full_reorth) && (alpha != 0.0f))
        {
            supdate_nu(&numax, &work[imu], &work[inu], j, &B[0], &B[ldb], anorm, epsn2);
        }

        // Reorthogonalize if necessary
        if ((full_reorth || numax > delta || force_reorth) && alpha != 0.0f)
        {
            if ((full_reorth) || (eta == 0.0f))
            {
                iwork[iidx] = 0;
                iwork[iidx + 1] = j - 1;
                iwork[iidx + 2] = j;
                iwork[iidx + 3] = j;
            } else if (!force_reorth) {
                scompute_int(&work[inu], j - 1, delta, eta, &iwork[iidx]);
            }

            sreorth(n, j - 1, V, ldv, &V[j * ldv], &alpha, &iwork[iidx], kappa, &work[is], cgs);
            sset_mu(j - 1, &work[inu], &iwork[iidx], eps);
            numax = eta;

            force_reorth = (force_reorth ? 0 : 1);
        }

        // Check whether an invariant subspace was found
        if ((alpha < anorm * epsn) && (j < *k - 1))
        {
            *rnorm = alpha;
            alpha = 0.0f;

            // Try to build an orthogonal subspace starting with a random vector
            sgetu0(1, m, n, j - 1, 3, &V[j * ldv], &alpha, V, ldv, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is], rng_state);

            if (alpha == 0.0f)
            {
                // We failed to generate a new random vector
                // in span(A^T) orthogonal to span(V(:,1:j-1)).
                // Most likely span(V(:,1:j-1)) is an invariant
                // subspace.
                *k = j;
                *ierr = -j-1;
                doption[2] = anorm;
                return;
            } else {
                // We have managed to generate a random vector
                // in span(A^T) orthogonal to V(:,1:j-1), so we
                // can continue the LBD and "deflate" the subspace
                // by setting alpha_{j} = 0.

                ssafescal(n, alpha, &V[j * ldv]);
                alpha = 0.0f;
                force_reorth = 1;
                if (delta > 0.0f) { full_reorth = 0; }
            }

        } else if ((j > 0) && (!full_reorth) && (j < *k - 1) && (delta * alpha < anorm * eps)) {
            *ierr = j;
        }
        B[j] = alpha;

        if (alpha != 0.0f) {
            ssafescal(n, alpha, &V[j * ldv]);
        }

        // beta_{j+1} * u_{j+1} = A * v_j - alpha_j * u_j
        aprod(0, m, n, &V[j * ldv], &U[(j + 1) * ldu], dparm, iparm);

        float neg_alpha = -alpha;
        saxpy_(&m, &neg_alpha, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
        beta = snrm2_(&m, &U[(j + 1) * ldu], &ione);

        // Extended local reorthogonalization
        if ((elr > 0) && (beta < kappa * alpha))
        {
            for (i = 0; i < elr; i++)
            {
                s = sdot_(&m, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
                float neg_s = -s;
                saxpy_(&m, &neg_s, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
                if (alpha != 0.0f)
                {
                    alpha = alpha + s;
                    B[j] = alpha;
                }
                s = snrm2_(&m, &U[(j + 1) * ldu], &ione);
                if (s >= kappa * beta) { break; }
                beta = s;
            }
            work[imu + j] = eps;
            beta = s;
        }

        B[j + ldb] = beta;
        amax = fmaxf(amax, beta);

        // Update estimate of ||A||_2
        if (j <= 0)
        {
            a1 = hypotf(B[0], B[ldb]);
        } else {
            a1 = B[j] / amax;
            a1 = amax * sqrtf(a1 * a1 + powf(B[j + ldb] / amax, 2) + a1 * B[j - 1 + ldb] / amax);
        }
        anorm = fmaxf(anorm, a1);

        // Update the mu recurrence
        if ((!full_reorth) && (beta != 0.0f))
        {
            supdate_mu(&mumax, &work[imu], &work[inu], j, &B[0], &B[ldb], anorm, epsn2);
        }

        // Reorthogonalize u_{j+1} if necessary
        if ((full_reorth || (mumax > delta) || force_reorth) && beta != 0.0f)
        {
            if ((full_reorth) || (eta == 0.0f))
            {
                iwork[iidx] = 0;
                iwork[iidx + 1] = j;
                iwork[iidx + 2] = j + 1;
                iwork[iidx + 3] = j + 1;
            } else if (!force_reorth) {
                scompute_int(&work[imu], j, delta, eta, &iwork[iidx]);
            } else {
                for (i = 0; i < 2 * j + 1; i++)
                {
                    if (iwork[iidx + i] == j)
                    {
                        iwork[iidx + i] = j + 1;
                        iwork[iidx + i + 1] = j + 1;
                        break;
                    }
                }
            }

            sreorth(m, j, U, ldu, &U[(j + 1) * ldu], &beta, &iwork[iidx], kappa, &work[is], cgs);
            sset_mu(j, &work[imu], &iwork[iidx], eps);
            mumax = eta;
            force_reorth = (force_reorth ? 0 : 1);
        }

        // Check whether an invariant subspace was found
        if ((beta < anorm * epsn) && (j < *k - 1))
        {
            *rnorm = beta;
            beta = 0.0f;
            // Try to build an orthogonal subspace starting with a random vector
            sgetu0(0, m, n, j, 3, &U[(j + 1) * ldu], &beta, U, ldu, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is], rng_state);

            if (beta == 0.0f)
            {
                // We failed to generate a new random vector
                // in span(A) orthogonal to span(U(:,1:j)).
                // Most likely span(U(:,1:j)) is an invariant
                // subspace.
                *k = j+1;
                *ierr = -j-1;
                doption[2] = anorm;
                return;
            } else {
                // We have managed to generate a random vector
                // in span(A) orthogonal to U(:,1:j), so we can
                // continue the LBD and "deflate" the subspace by
                // setting beta_{j+1} = 0.
                ssafescal(m, beta, &U[(j + 1) * ldu]);
                beta = 0.0f;
                force_reorth = 1;
                if (delta > 0.0f) { full_reorth = 0; }
            }
        } else if (!full_reorth && (j < *k - 1) && (delta * beta < anorm * eps)) {
            *ierr = j;
        }

        B[j + ldb] = beta;

        if ((beta != 0.0f) && (beta != 1.0f))
        {
            ssafescal(m, beta, &U[(j + 1) * ldu]);
        }
        *rnorm = beta;
    }

    doption[2] = anorm;
}


/**
 * @brief Perform Lanczos bidiagonalization with partial reorthogonalization (double precision).
 *
 * This function computes a partial Lanczos bidiagonalization of a matrix A using
 * double precision arithmetic. The bidiagonalization is performed with partial
 * reorthogonalization to maintain numerical stability and efficiency.
 *
 * @param m         The number of rows in the matrix A.
 * @param n         The number of columns in the matrix A.
 * @param k0        The starting index for the bidiagonalization process.
 * @param k         Pointer to the number of Lanczos steps to perform.
 * @param aprod     Function pointer for matrix-vector multiplication with A or A^T.
 * @param U         Pointer to the left singular vectors (output).
 * @param ldu       The leading dimension of the array U.
 * @param V         Pointer to the right singular vectors (output).
 * @param ldv       The leading dimension of the array V.
 * @param B         Pointer to the bidiagonal matrix (output).
 * @param ldb       The leading dimension of the array B.
 * @param rnorm     Pointer to the norm of the residual vector (output).
 * @param doption   Array of double precision options for controlling the algorithm.
 * @param ioption   Array of integer options for controlling the algorithm.
 * @param work      Workspace array for intermediate computations.
 * @param iwork     Integer workspace array for intermediate computations.
 * @param dparm     Array of double precision parameters for the algorithm.
 * @param iparm     Array of integer parameters for the algorithm.
 * @param ierr      Pointer to the error flag (output).
 * @param rng_state Pointer to the random number generator state.
 */
void dlanbpro(int m, int n, int k0, int* k, PROPACK_aprod_d aprod,
              double* U, int ldu, double* V, int ldv, double* B, int ldb,
              double* rnorm, double* doption, int* ioption, double* work, int* iwork,
              double* dparm, int* iparm, int* ierr, uint64_t* rng_state) {

    // Constants
    const double FUDGE = 1.01;
    const double kappa = sqrt(2.0) / 2.0;

    // Machine constants
    const double eps = DBL_EPSILON*0.5;
    const double eps34 = pow(eps, 0.75);
    const double epsn = (double)int_max(m, n) * eps;
    const double epsn2 = sqrt((double)int_max(m, n)) * eps;

    // Local variables
    int i, j, inu, imu, is, iidx, j0;
    double s, mumax, numax, alpha, beta, a1, b1, amax, anormest;
    int force_reorth, full_reorth;
    int ione = 1;

    // Set default parameters
    double delta, eta, anorm;
    int cgs, elr;

    if (doption[0] < 0.0)
    {
        delta = sqrt(eps / *k);
    } else {
        delta = doption[0];
    }

    if (doption[1] < 0.0)
    {
        eta = eps34 / sqrt((double)*k);
    } else {
        eta = doption[1];
    }

    if ((delta <= eta) || (delta == 0.0))
    {
        full_reorth = 1;
    } else {
        full_reorth = 0;
    }

    if (doption[2] > 0.0)
    {
        anorm = doption[2];
    } else if (k0 > 0) {
        anorm = hypot(B[0], B[ldb]);
        if (anorm <= 0.0) {
            *ierr = -1;
            doption[2] = anorm;
            return;
        }
    } else {
        anorm = 0.0;
    }

    cgs = ioption[0];
    elr = ioption[1];
    *ierr = 0;

    // Get starting vector if needed
    if (*rnorm == 0.0) {
        dgetu0(0, m, n, k0, 3, &U[k0 * ldu], rnorm, U, ldu, aprod, dparm, iparm, ierr, cgs, &anormest, work, rng_state);
        anorm = fmax(anorm, anormest);
    }

    // Set pointers into work array
    // work layout: [mu(k+1), nu(k+1), workspace...]
    imu = 0;
    inu = *k + 1;
    is = 2 * (*k + 1);
    iidx = 0;
    for (i = 0; i < int_max(m, n) + 2 * (*k) + 2; i++) { work[i] = 0.0; }
    for (i = 0; i < 2 * (*k) + 1; i++) { iwork[i] = 0; }

    // Prepare Lanczos iteration
    if (k0 == 0) {
        amax = 0.0;
        alpha = 0.0;
        beta = *rnorm;
        force_reorth = 0;

        // Compute ||A x|| / ||x|| for a random vector x
        // to make it less likely that ||A|| is grossly
        // underestimated at the beginning of the iteration.
        if (n > m) {
            dgetu0(0, m, n, 0, 1, &work[is], &s, U, ldu, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is + m], rng_state);
        } else {
            dgetu0(1, m, n, 0, 1, &work[is], &s, V, ldv, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is + n], rng_state);
        }
        anorm = fmax(anorm, FUDGE * anormest);
        j0 = 0;

        if (beta != 0.0) { dsafescal(m, beta, &U[0]); }
        work[imu] = 1.0;
        work[inu] = 1.0;
    } else {
        force_reorth = 1;
        alpha = B[k0 - 1];
        beta = *rnorm;

        if ((k0 < *k) && (beta * delta < anorm * eps)) {
            full_reorth = 1;
            *ierr = k0;
        }

        iwork[iidx] = 0;
        iwork[iidx + 1] = k0 - 1;
        iwork[iidx + 2] = k0;
        iwork[iidx + 3] = k0;

        dscal_(&m, rnorm, &U[k0 * ldu], &ione);
        dreorth(m, k0, U, ldu, &U[k0 * ldu], rnorm, &iwork[iidx], kappa, &work[is], cgs);
        dsafescal(m, *rnorm, &U[k0 * ldu]);
        dset_mu(k0, &work[imu], &iwork[iidx], epsn2);
        dset_mu(k0, &work[inu], &iwork[iidx], epsn2);
        beta = *rnorm;

        // Estimate ||B||_2^2 as ||B^T * B||_1
        B[k0 - 1 + ldb] = beta;
        amax = 0.0;
        for (j = 0; j < k0; j++)
        {
            amax = fmax(amax, fmax(B[j], B[j + ldb]));
            if (j == 0)
            {
                anorm = fmax(anorm, FUDGE * alpha);
            } else if (j == 1) {
                a1 = B[ldb] / amax;
                a1 = FUDGE * amax * sqrt(pow(B[0] / amax, 2) + a1 * a1 + B[1] / amax * a1);
                anorm = fmax(anorm, a1);
            } else {
                a1 = B[j - 1] / amax;
                b1 = B[j - 1 + ldb] / amax;
                a1 = FUDGE * amax * sqrt(a1 * a1 + b1 * b1 + a1 * B[j - 2 + ldb] / amax + B[j] / amax * b1);
                anorm = fmax(anorm, a1);
            }
        }
        j0 = k0;
    }

    numax = 0.0;
    mumax = 0.0;

    // Start Lanczos bidiagonalization iteration
    for (j = j0; j < *k; j++)
    {
        // alpha_j * v_j = A^T * u_j - beta_j * v_{j-1}
        aprod(1, m, n, &U[j * ldu], &V[j * ldv], dparm, iparm);

        if (j == 0)
        {
            alpha = dnrm2_(&n, &V[j * ldv], &ione);
            anorm = fmax(anorm, FUDGE * alpha);
        } else {
            double neg_beta = -beta;
            daxpy_(&n, &neg_beta, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
            alpha = dnrm2_(&n, &V[j * ldv], &ione);

            // Extended local reorthogonalization
            if ((j > 0) && (elr > 0) && (alpha < kappa * beta))
            {
                for (i = 0; i < elr; i++)
                {
                    s = ddot_(&n, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
                    double neg_s = -s;
                    daxpy_(&n, &neg_s, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
                    if (beta != 0.0)
                    {
                        beta = beta + s;
                        B[j - 1 + ldb] = beta;
                    }
                    s = dnrm2_(&n, &V[j * ldv], &ione);
                    if (s >= kappa * alpha) { break; }
                    alpha = s;
                }
                work[inu + j - 1] = eps;
                alpha = s;
            }

            B[j] = alpha;
            amax = fmax(amax, alpha);

            // Update estimate of ||A||_2
            if (j == 1) {
                a1 = B[ldb] / amax;
                a1 = FUDGE * amax * sqrt(pow(B[0] / amax, 2) + a1 * a1 + B[1] / amax * a1);
            } else {
                a1 = B[j - 1] / amax;
                b1 = B[j - 1 + ldb] / amax;
                a1 = FUDGE * amax * sqrt(a1 * a1 + b1 * b1 + a1 * B[j - 2 + ldb] / amax + B[j] / amax * b1);
            }
            anorm = fmax(anorm, a1);
        }

        // Update the nu recurrence
        if ((!full_reorth) && (alpha != 0.0))
        {
            dupdate_nu(&numax, &work[imu], &work[inu], j, &B[0], &B[ldb], anorm, epsn2);
        }

        // Reorthogonalize if necessary
        if ((full_reorth || numax > delta || force_reorth) && alpha != 0.0)
        {
            if ((full_reorth) || (eta == 0.0))
            {
                iwork[iidx] = 0;
                iwork[iidx + 1] = j - 1;
                iwork[iidx + 2] = j;
                iwork[iidx + 3] = j;
            } else if (!force_reorth) {
                dcompute_int(&work[inu], j - 1, delta, eta, &iwork[iidx]);
            }
            dreorth(n, j - 1, V, ldv, &V[j * ldv], &alpha, &iwork[iidx], kappa, &work[is], cgs);
            dset_mu(j - 1, &work[inu], &iwork[iidx], eps);
            numax = eta;

            force_reorth = (force_reorth ? 0 : 1);
        }

        // Check whether an invariant subspace was found
        if ((alpha < anorm * epsn) && (j < *k - 1))
        {
            *rnorm = alpha;
            alpha = 0.0;

            // Try to build an orthogonal subspace starting with a random vector
            dgetu0(1, m, n, j - 1, 3, &V[j * ldv], &alpha, V, ldv, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is], rng_state);

            if (alpha == 0.0)
            {
                // We failed to generate a new random vector
                // in span(A^T) orthogonal to span(V(:,1:j-1)).
                // Most likely span(V(:,1:j-1)) is an invariant
                // subspace.
                *k = j;
                *ierr = -j-1;
                doption[2] = anorm;
                return;
            } else {
                // We have managed to generate a random vector
                // in span(A^T) orthogonal to V(:,1:j-1), so we
                // can continue the LBD and "deflate" the subspace
                // by setting alpha_{j} = 0.

                dsafescal(n, alpha, &V[j * ldv]);
                alpha = 0.0;
                force_reorth = 1;
                if (delta > 0.0) { full_reorth = 0; }
            }

        } else if ((j > 0) && (!full_reorth) && (j < *k - 1) && (delta * alpha < anorm * eps)) {
            *ierr = j;
        }
        B[j] = alpha;

        if (alpha != 0.0)
        {
            dsafescal(n, alpha, &V[j * ldv]);
        }

        // beta_{j+1} * u_{j+1} = A * v_j - alpha_j * u_j
        aprod(0, m, n, &V[j * ldv], &U[(j + 1) * ldu], dparm, iparm);

        double neg_alpha = -alpha;
        daxpy_(&m, &neg_alpha, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
        beta = dnrm2_(&m, &U[(j + 1) * ldu], &ione);

        // Extended local reorthogonalization
        if ((elr > 0) && (beta < kappa * alpha))
        {
            for (i = 0; i < elr; i++)
            {
                s = ddot_(&m, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
                double neg_s = -s;
                daxpy_(&m, &neg_s, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
                if (alpha != 0.0)
                {
                    alpha = alpha + s;
                    B[j] = alpha;
                }
                s = dnrm2_(&m, &U[(j + 1) * ldu], &ione);
                if (s >= kappa * beta) { break; }
                beta = s;
            }
            work[imu + j] = eps;
            beta = s;
        }

        B[j + ldb] = beta;
        amax = fmax(amax, beta);

        // Update estimate of ||A||_2
        if (j <= 0)
        {
            a1 = hypot(B[0], B[ldb]);
        } else {
            a1 = B[j] / amax;
            a1 = amax * sqrt(a1 * a1 + pow(B[j + ldb] / amax, 2) + a1 * B[j - 1 + ldb] / amax);
        }
        anorm = fmax(anorm, a1);
        // Update the mu recurrence
        if ((!full_reorth) && (beta != 0.0))
        {
            dupdate_mu(&mumax, &work[imu], &work[inu], j, &B[0], &B[ldb], anorm, epsn2);
        }

        // Reorthogonalize u_{j+1} if necessary
        if ((full_reorth || (mumax > delta) || force_reorth) && (beta != 0.0))
        {
            if ((full_reorth) || (eta == 0.0))
            {
                iwork[iidx] = 0;
                iwork[iidx + 1] = j;
                iwork[iidx + 2] = j + 1;
                iwork[iidx + 3] = j + 1;
            } else if (!force_reorth) {
                dcompute_int(&work[imu], j, delta, eta, &iwork[iidx]);
            } else {
                for (i = 0; i < 2 * j + 1; i++)
                {
                    if (iwork[iidx + i] == j)
                    {
                        iwork[iidx + i] = j + 1;
                        iwork[iidx + i + 1] = j + 1;
                        break;
                    }
                }
            }

            dreorth(m, j, U, ldu, &U[(j + 1) * ldu], &beta, &iwork[iidx], kappa, &work[is], cgs);
            dset_mu(j, &work[imu], &iwork[iidx], eps);
            mumax = eta;

            force_reorth = (force_reorth ? 0 : 1);
        }

        // Check whether an invariant subspace was found
        if ((beta < anorm * epsn) && (j < *k - 1))
        {
            *rnorm = beta;
            beta = 0.0;
            // Try to build an orthogonal subspace starting with a random vector
            dgetu0(0, m, n, j, 3, &U[(j + 1) * ldu], &beta, U, ldu, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is], rng_state);

            if (beta == 0.0)
            {
                // We failed to generate a new random vector
                // in span(A) orthogonal to span(U(:,1:j)).
                // Most likely span(U(:,1:j)) is an invariant
                // subspace.
                *k = j+1;
                *ierr = -j-1;
                doption[2] = anorm;
                return;
            } else {
                // We have managed to generate a random vector
                // in span(A) orthogonal to U(:,1:j), so we can
                // continue the LBD and "deflate" the subspace by
                // setting beta_{j+1} = 0.
                dsafescal(m, beta, &U[(j + 1) * ldu]);
                beta = 0.0;
                force_reorth = 1;
                if (delta > 0.0) { full_reorth = 0; }
            }
        } else if (!full_reorth && (j < *k - 1) && (delta * beta < anorm * eps)) {
            *ierr = j;
        }

        B[j + ldb] = beta;
        if ((beta != 0.0) && (beta != 1.0))
        {
            dsafescal(m, beta, &U[(j + 1) * ldu]);
        }
        *rnorm = beta;
    }

    doption[2] = anorm;
}


/**
 * @brief Perform Lanczos bidiagonalization with partial reorthogonalization (single precision complex).
 *
 * This function computes a partial Lanczos bidiagonalization of a matrix A using
 * single precision complex arithmetic. The bidiagonalization is performed with partial
 * reorthogonalization to maintain numerical stability and efficiency.
 *
 * @param m         The number of rows in the matrix A.
 * @param n         The number of columns in the matrix A.
 * @param k0        The starting index for the bidiagonalization process.
 * @param k         Pointer to the number of Lanczos steps to perform.
 * @param aprod     Function pointer for matrix-vector multiplication with A or A^T.
 * @param U         Pointer to the left singular vectors (output).
 * @param ldu       The leading dimension of the array U.
 * @param V         Pointer to the right singular vectors (output).
 * @param ldv       The leading dimension of the array V.
 * @param B         Pointer to the bidiagonal matrix (output).
 * @param ldb       The leading dimension of the array B.
 * @param rnorm     Pointer to the norm of the residual vector (output).
 * @param doption   Array of double precision options for controlling the algorithm.
 * @param ioption   Array of integer options for controlling the algorithm.
 * @param work      Workspace array for intermediate computations.
 * @param iwork     Integer workspace array for intermediate computations.
 * @param dparm     Array of double precision parameters for the algorithm.
 * @param iparm     Array of integer parameters for the algorithm.
 * @param ierr      Pointer to the error flag (output).
 * @param rng_state Pointer to the random number generator state.
 */
void clanbpro(
    int m, int n, int k0, int* k, PROPACK_aprod_c aprod, PROPACK_CPLXF_TYPE* U, int ldu,
    PROPACK_CPLXF_TYPE* V, int ldv, float* B, int ldb, float* rnorm, float* soption,
    int* ioption, float* swork, PROPACK_CPLXF_TYPE* cwork, int* iwork, PROPACK_CPLXF_TYPE* cparm,
    int* iparm, int* ierr, uint64_t* rng_state)
{

    // Constants
    const float FUDGE = 1.01f;
    const float kappa = sqrtf(2.0f) / 2.0f;

    // Machine constants
    const float eps = FLT_EPSILON*0.5;
    const float eps34 = powf(eps, 0.75f);
    const float epsn = (float)int_max(m, n) * eps;
    const float epsn2 = sqrtf((float)int_max(m, n)) * eps;

    // Local variables
    int i, j, inu, imu, is, iidx, j0;
    float mumax, numax, alpha, beta, a1, b1, amax, anormest, nrm;
    PROPACK_CPLXF_TYPE s;
    int force_reorth, full_reorth;
    int ione = 1;

    // Set default parameters
    float delta, eta, anorm;
    int cgs, elr;

    if (soption[0] < 0.0f)
    {
        delta = sqrtf(eps / *k);
    } else {
        delta = soption[0];
    }

    if (soption[1] < 0.0f)
    {
        eta = eps34 / sqrtf((float)*k);
    } else {
        eta = soption[1];
    }

    if ((delta <= eta) || (delta == 0.0f))
    {
        full_reorth = 1;
    } else {
        full_reorth = 0;
    }

    if (soption[2] > 0.0f)
    {
        anorm = soption[2];
    } else if (k0 > 0) {
        anorm = hypotf(B[0], B[ldb]);
        if (anorm <= 0.0f) {
            *ierr = -1;
            soption[2] = anorm;
            return;
        }
    } else {
        anorm = 0.0f;
    }

    cgs = ioption[0];
    elr = ioption[1];
    *ierr = 0;

    // Get starting vector if needed
    if (*rnorm == 0.0f)
    {
        cgetu0(0, m, n, k0, 3, &U[k0 * ldu], rnorm, U, ldu, aprod, cparm, iparm, ierr, cgs, &anormest, cwork, rng_state);
        anorm = fmaxf(anorm, anormest);
    }

    // Set pointers into work array
    // work layout: [mu(k+1), nu(k+1), workspace...]
    imu = 0;
    inu = *k + 1;
    is = 0;
    iidx = 0;
    for (i = 0; i < 2 * (*k) + 2; i++) { swork[i] = 0.0f; }
    for (i = 0; i < 2 * (*k) + 1; i++) { iwork[i] = 0; }
    for (i = 0; i < int_max(m, n); i++) { cwork[i] = PROPACK_cplxf(0.0f, 0.0f); }

    // Prepare Lanczos iteration
    if (k0 == 0) {
        amax = 0.0f;
        alpha = 0.0f;
        beta = *rnorm;
        force_reorth = 0;

        // Compute ||A x|| / ||x|| for a random vector x
        // to make it less likely that ||A|| is grossly
        // underestimated at the beginning of the iteration.
        if (n > m) {
            cgetu0(0, m, n, 0, 1, &cwork[is], &nrm, U, ldu, aprod, cparm, iparm, ierr, cgs, &anormest, &cwork[is + m], rng_state);
        } else {
            cgetu0(1, m, n, 0, 1, &cwork[is], &nrm, V, ldv, aprod, cparm, iparm, ierr, cgs, &anormest, &cwork[is + n], rng_state);
        }
        anorm = fmaxf(anorm, FUDGE * anormest);
        j0 = 0;

        if (beta != 0.0f) { csafescal(m, beta, &U[0]); }
        swork[imu] = 1.0f;
        swork[inu] = 1.0f;
    } else {
        force_reorth = 1;
        alpha = B[k0 - 1];
        beta = *rnorm;

        if ((k0 < *k) && (beta * delta < anorm * eps)) {
            full_reorth = 1;
            *ierr = k0;
        }

        iwork[iidx] = 0;
        iwork[iidx + 1] = k0 - 1;
        iwork[iidx + 2] = k0;
        iwork[iidx + 3] = k0;

        csscal_(&m, rnorm, &U[k0 * ldu], &ione);
        creorth(m, k0, U, ldu, &U[k0 * ldu], rnorm, &iwork[iidx], kappa, &cwork[is], cgs);
        csafescal(m, *rnorm, &U[k0 * ldu]);
        sset_mu(k0, &swork[imu], &iwork[iidx], epsn2);
        sset_mu(k0, &swork[inu], &iwork[iidx], epsn2);
        beta = *rnorm;

        // Estimate ||B||_2^2 as ||B^T * B||_1
        B[k0 - 1 + ldb] = beta;
        amax = 0.0f;
        for (j = 0; j < k0; j++)
        {
            amax = fmaxf(amax, fmaxf(B[j], B[j + ldb]));
            if (j == 0)
            {
                anorm = fmaxf(anorm, FUDGE * alpha);
            } else if (j == 1) {
                a1 = B[ldb] / amax;
                a1 = FUDGE * amax * sqrtf(powf(B[0] / amax, 2) + a1 * a1 + B[1] / amax * a1);
                anorm = fmaxf(anorm, a1);
            } else {
                a1 = B[j - 1] / amax;
                b1 = B[j - 1 + ldb] / amax;
                a1 = FUDGE * amax * sqrtf(a1 * a1 + b1 * b1 + a1 * B[j - 2 + ldb] / amax + B[j] / amax * b1);
                anorm = fmaxf(anorm, a1);
            }
        }
        j0 = k0;
    }

    numax = 0.0f;
    mumax = 0.0f;

    // Start Lanczos bidiagonalization iteration
    for (j = j0; j < *k; j++)
    {
        // alpha_j * v_j = A^T * u_j - beta_j * v_{j-1}
        aprod(1, m, n, &U[j * ldu], &V[j * ldv], cparm, iparm);

        if (j == 0)
        {
            alpha = scnrm2_(&n, &V[j * ldv], &ione);
            anorm = fmaxf(anorm, FUDGE * alpha);
        } else {
            PROPACK_CPLXF_TYPE neg_beta = PROPACK_cplxf(-beta, 0.0);
            caxpy_(&n, &neg_beta, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
            alpha = scnrm2_(&n, &V[j * ldv], &ione);

            // Extended local reorthogonalization
            if ((j > 0) && (elr > 0) && (alpha < kappa * beta))
            {
                for (i = 0; i < elr; i++)
                {
                    s = cdotc_(&n, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
                    PROPACK_CPLXF_TYPE neg_s = PROPACK_cplxf(-creal(s), -cimag(s));
                    caxpy_(&n, &neg_s, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
                    nrm = scnrm2_(&n, &V[j * ldv], &ione);
                    if (nrm >= kappa * alpha) { break; }
                    alpha = nrm;
                }
                swork[inu + j - 1] = eps;
                alpha = nrm;
            }

            B[j] = alpha;
            amax = fmaxf(amax, alpha);

            // Update estimate of ||A||_2
            if (j == 1)
            {
                a1 = B[ldb] / amax;
                a1 = FUDGE * amax * sqrtf(powf(B[0] / amax, 2) + a1 * a1 + B[1] / amax * a1);
            } else {
                a1 = B[j - 1] / amax;
                b1 = B[j - 1 + ldb] / amax;
                a1 = FUDGE * amax * sqrtf(a1 * a1 + b1 * b1 + a1 * B[j - 2 + ldb] / amax + B[j] / amax * b1);
            }
            anorm = fmaxf(anorm, a1);
        }

        // Update the nu recurrence
        if ((!full_reorth) && (alpha != 0.0f))
        {
            supdate_nu(&numax, &swork[imu], &swork[inu], j, &B[0], &B[ldb], anorm, epsn2);
        }

        // Reorthogonalize if necessary
        if ((full_reorth || (numax > delta) || force_reorth) && alpha != 0.0f)
        {
            if ((full_reorth) || (eta == 0.0f))
            {
                iwork[iidx] = 0;
                iwork[iidx + 1] = j - 1;
                iwork[iidx + 2] = j;
                iwork[iidx + 3] = j;
            } else if (!force_reorth) {
                scompute_int(&swork[inu], j - 1, delta, eta, &iwork[iidx]);
            }

            creorth(n, j - 1, V, ldv, &V[j * ldv], &alpha, &iwork[iidx], kappa, &cwork[is], cgs);
            sset_mu(j - 1, &swork[inu], &iwork[iidx], eps);
            numax = eta;

            force_reorth = (force_reorth ? 0 : 1);
        }

        // Check whether an invariant subspace was found
        if ((alpha < anorm * epsn) && (j < *k - 1))
        {
            *rnorm = alpha;
            alpha = 0.0f;

            // Try to build an orthogonal subspace starting with a random vector
            cgetu0(1, m, n, j - 1, 3, &V[j * ldv], &alpha, V, ldv, aprod, cparm, iparm, ierr, cgs, &anormest, &cwork[is], rng_state);

            if (alpha == 0.0f)
            {
                // We failed to generate a new random vector
                // in span(A^H) orthogonal to span(V(:,1:j-1)).
                // Most likely span(V(:,1:j-1)) is an invariant
                // subspace.
                *k = j;
                *ierr = -j-1;
                soption[2] = anorm;
                return;
            } else {
                // We have managed to generate a random vector
                // in span(A^H) orthogonal to V(:,1:j-1), so we
                // can continue the LBD and "deflate" the subspace
                // by setting alpha_{j} = 0.

                csafescal(n, alpha, &V[j * ldv]);
                alpha = 0.0f;
                force_reorth = 1;
                if (delta > 0.0f) { full_reorth = 0; }
            }

        } else if ((j > 0) && (!full_reorth) && (j < *k - 1) && (delta * alpha < anorm * eps)) {
            *ierr = j;
        }
        B[j] = alpha;

        if (alpha != 0.0f)
        {
            csafescal(n, alpha, &V[j * ldv]);
        }

        // beta_{j+1} * u_{j+1} = A * v_j - alpha_j * u_j
        aprod(0, m, n, &V[j * ldv], &U[(j + 1) * ldu], cparm, iparm);

        PROPACK_CPLXF_TYPE neg_alpha = PROPACK_cplxf(-alpha, 0.0);
        caxpy_(&m, &neg_alpha, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
        beta = scnrm2_(&m, &U[(j + 1) * ldu], &ione);

        // Extended local reorthogonalization
        if ((elr > 0) && (beta < kappa * alpha))
        {
            for (i = 0; i < elr; i++)
            {
                s = cdotc_(&m, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
                PROPACK_CPLXF_TYPE neg_s = PROPACK_cplxf(-creal(s), -cimag(s));
                caxpy_(&m, &neg_s, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
                nrm = scnrm2_(&m, &U[(j + 1) * ldu], &ione);
                if (nrm >= kappa * beta) { break; }
                beta = nrm;
            }
            swork[imu + j] = eps;
            beta = nrm;
        }

        B[j + ldb] = beta;
        amax = fmaxf(amax, beta);

        // Update estimate of ||A||_2
        if (j <= 0)
        {
            a1 = hypotf(B[0], B[ldb]);
        } else {
            a1 = B[j] / amax;
            a1 = amax * sqrtf(a1 * a1 + powf(B[j + ldb] / amax, 2) + a1 * B[j - 1 + ldb] / amax);
        }
        anorm = fmaxf(anorm, a1);

        // Update the mu recurrence
        if ((!full_reorth) && (beta != 0.0f))
        {
            supdate_mu(&mumax, &swork[imu], &swork[inu], j, &B[0], &B[ldb], anorm, epsn2);
        }

        // Reorthogonalize u_{j+1} if necessary
        if ((full_reorth || mumax > delta || force_reorth) && beta != 0.0f)
        {
            if ((full_reorth) || (eta == 0.0f))
            {
                iwork[iidx] = 0;
                iwork[iidx + 1] = j;
                iwork[iidx + 2] = j + 1;
                iwork[iidx + 3] = j + 1;
            } else if (!force_reorth) {
                scompute_int(&swork[imu], j, delta, eta, &iwork[iidx]);
            } else {
                for (i = 0; i < 2 * j + 1; i++)
                {
                    if (iwork[iidx + i] == j)
                    {
                        iwork[iidx + i] = j + 1;
                        iwork[iidx + i + 1] = j + 1;
                        break;
                    }
                }
            }

            creorth(m, j, U, ldu, &U[(j + 1) * ldu], &beta, &iwork[iidx], kappa, &cwork[is], cgs);
            sset_mu(j, &swork[imu], &iwork[iidx], eps);
            mumax = eta;

            force_reorth = (force_reorth ? 0 : 1);
        }

        // Check whether an invariant subspace was found
        if ((beta < anorm * epsn) && (j < *k - 1))
        {
            *rnorm = beta;
            beta = 0.0f;

            // Try to build an orthogonal subspace starting with a random vector
            cgetu0(0, m, n, j, 3, &U[(j + 1) * ldu], &beta, U, ldu, aprod, cparm, iparm, ierr, cgs, &anormest, &cwork[is], rng_state);

            if (beta == 0.0f)
            {
                // We failed to generate a new random vector
                // in span(A) orthogonal to span(U(:,1:j)).
                // Most likely span(U(:,1:j)) is an invariant
                // subspace.
                *k = j+1;
                *ierr = -j-1;
                soption[2] = anorm;
                return;
            } else {
                // We have managed to generate a random vector
                // in span(A) orthogonal to U(:,1:j), so we can
                // continue the LBD and "deflate" the subspace by
                // setting beta_{j+1} = 0.
                csafescal(m, beta, &U[(j + 1) * ldu]);
                beta = 0.0f;
                force_reorth = 1;
                if (delta > 0.0f) { full_reorth = 0; }
            }
        } else if (!full_reorth && (j < *k - 1) && (delta * beta < anorm * eps)) {
            *ierr = j;
        }

        B[j + ldb] = beta;
        if ((beta != 0.0f) && (beta != 1.0f))
        {
            csafescal(m, beta, &U[(j + 1) * ldu]);
        }
        *rnorm = beta;
    }

    soption[2] = anorm;
}


/**
 * @brief Perform Lanczos bidiagonalization with partial reorthogonalization (double precision complex).
 *
 * This function computes a partial Lanczos bidiagonalization of a matrix A using
 * double precision complex arithmetic. The bidiagonalization is performed with partial
 * reorthogonalization to maintain numerical stability and efficiency.
 *
 * @param m         The number of rows in the matrix A.
 * @param n         The number of columns in the matrix A.
 * @param k0        The starting index for the bidiagonalization process.
 * @param k         Pointer to the number of Lanczos steps to perform.
 * @param aprod     Function pointer for matrix-vector multiplication with A or A^T.
 * @param U         Pointer to the left singular vectors (output).
 * @param ldu       The leading dimension of the array U.
 * @param V         Pointer to the right singular vectors (output).
 * @param ldv       The leading dimension of the array V.
 * @param B         Pointer to the bidiagonal matrix (output).
 * @param ldb       The leading dimension of the array B.
 * @param rnorm     Pointer to the norm of the residual vector (output).
 * @param doption   Array of double precision options for controlling the algorithm.
 * @param ioption   Array of integer options for controlling the algorithm.
 * @param work      Workspace array for intermediate computations.
 * @param iwork     Integer workspace array for intermediate computations.
 * @param dparm     Array of double precision parameters for the algorithm.
 * @param iparm     Array of integer parameters for the algorithm.
 * @param ierr      Pointer to the error flag (output).
 * @param rng_state Pointer to the random number generator state.
 */
void zlanbpro(
    int m, int n, int k0, int* k, PROPACK_aprod_z aprod, PROPACK_CPLX_TYPE* U, int ldu,
    PROPACK_CPLX_TYPE* V, int ldv, double* B, int ldb, double* rnorm, double* doption,
    int* ioption, double* dwork, PROPACK_CPLX_TYPE* zwork, int* iwork, PROPACK_CPLX_TYPE* zparm,
    int* iparm, int* ierr, uint64_t* rng_state)
{

    // Constants
    const double FUDGE = 1.01;
    const double kappa = sqrt(2.0) / 2.0;

    // Machine constants
    const double eps = DBL_EPSILON*0.5;
    const double eps34 = pow(eps, 0.75);
    const double epsn = (double)int_max(m, n) * eps;
    const double epsn2 = sqrt((double)int_max(m, n)) * eps;

    // Local variables
    int i, j, inu, imu, is, iidx, j0;
    double mumax, numax, alpha, beta, a1, b1, amax, anormest, nrm;
    PROPACK_CPLX_TYPE s;
    int force_reorth, full_reorth;
    int ione = 1;

    // Set default parameters
    double delta, eta, anorm;
    int cgs, elr;

    if (doption[0] < 0.0)
    {
        delta = sqrt(eps / *k);
    } else {
        delta = doption[0];
    }

    if (doption[1] < 0.0)
    {
        eta = eps34 / sqrt((double)*k);
    } else {
        eta = doption[1];
    }

    if ((delta <= eta) || (delta == 0.0))
    {
        full_reorth = 1;
    } else {
        full_reorth = 0;
    }

    if (doption[2] > 0.0)
    {
        anorm = doption[2];
    } else if (k0 > 0) {
        anorm = hypot(B[0], B[ldb]);
        if (anorm <= 0.0) {
            *ierr = -1;
            doption[2] = anorm;
            return;
        }
    } else {
        anorm = 0.0;
    }

    cgs = ioption[0];
    elr = ioption[1];
    *ierr = 0;

    // Get starting vector if needed
    if (*rnorm == 0.0)
    {
        zgetu0(0, m, n, k0, 3, &U[k0 * ldu], rnorm, U, ldu, aprod, zparm, iparm, ierr, cgs, &anormest, zwork, rng_state);
        anorm = fmax(anorm, anormest);
    }

    // Set pointers into work array
    // work layout: [mu(k+1), nu(k+1), workspace...]
    imu = 0;
    inu = *k + 1;
    is = 0;
    iidx = 0;
    for (i = 0; i < 2 * (*k) + 2; i++) { dwork[i] = 0.0; }
    for (i = 0; i < 2 * (*k) + 1; i++) { iwork[i] = 0; }
    for (i = 0; i < int_max(m, n); i++) { zwork[i] = PROPACK_cplx(0.0, 0.0); }

    // Prepare Lanczos iteration
    if (k0 == 0) {
        amax = 0.0;
        alpha = 0.0;
        beta = *rnorm;
        force_reorth = 0;

        // Compute ||A x|| / ||x|| for a random vector x
        // to make it less likely that ||A|| is grossly
        // underestimated at the beginning of the iteration.
        if (n > m) {
            zgetu0(0, m, n, 0, 1, &zwork[is], &nrm, U, ldu, aprod, zparm, iparm, ierr, cgs, &anormest, &zwork[is + m], rng_state);
        } else {
            zgetu0(1, m, n, 0, 1, &zwork[is], &nrm, V, ldv, aprod, zparm, iparm, ierr, cgs, &anormest, &zwork[is + n], rng_state);
        }
        anorm = fmax(anorm, FUDGE * anormest);
        j0 = 0;

        if (beta != 0.0) { zsafescal(m, beta, &U[0]); }
        dwork[imu] = 1.0;
        dwork[inu] = 1.0;
    } else {
        force_reorth = 1;
        alpha = B[k0 - 1];
        beta = *rnorm;

        if ((k0 < *k) && (beta * delta < anorm * eps)) {
            full_reorth = 1;
            *ierr = k0;
        }

        iwork[iidx] = 0;
        iwork[iidx + 1] = k0 - 1;
        iwork[iidx + 2] = k0;
        iwork[iidx + 3] = k0;

        zdscal_(&m, rnorm, &U[k0 * ldu], &ione);
        zreorth(m, k0, U, ldu, &U[k0 * ldu], rnorm, &iwork[iidx], kappa, &zwork[is], cgs);
        zsafescal(m, *rnorm, &U[k0 * ldu]);
        dset_mu(k0, &dwork[imu], &iwork[iidx], epsn2);
        dset_mu(k0, &dwork[inu], &iwork[iidx], epsn2);
        beta = *rnorm;

        // Estimate ||B||_2^2 as ||B^T * B||_1
        B[k0 - 1 + ldb] = beta;
        amax = 0.0;
        for (j = 0; j < k0; j++)
        {
            amax = fmax(amax, fmax(B[j], B[j + ldb]));
            if (j == 0)
            {
                anorm = fmax(anorm, FUDGE * alpha);
            } else if (j == 1) {
                a1 = B[ldb] / amax;
                a1 = FUDGE * amax * sqrt(pow(B[0] / amax, 2) + a1 * a1 + B[1] / amax * a1);
                anorm = fmax(anorm, a1);
            } else {
                a1 = B[j - 1] / amax;
                b1 = B[j - 1 + ldb] / amax;
                a1 = FUDGE * amax * sqrt(a1 * a1 + b1 * b1 + a1 * B[j - 2 + ldb] / amax + B[j] / amax * b1);
                anorm = fmax(anorm, a1);
            }
        }
        j0 = k0;
    }

    numax = 0.0;
    mumax = 0.0;

    // Start Lanczos bidiagonalization iteration
    for (j = j0; j < *k; j++)
    {
        // alpha_j * v_j = A^T * u_j - beta_j * v_{j-1}
        aprod(1, m, n, &U[j * ldu], &V[j * ldv], zparm, iparm);

        if (j == 0)
        {
            alpha = dznrm2_(&n, &V[j * ldv], &ione);
            anorm = fmax(anorm, FUDGE * alpha);
        } else {
            PROPACK_CPLX_TYPE neg_beta = PROPACK_cplx(-beta, 0.0);
            zaxpy_(&n, &neg_beta, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
            alpha = dznrm2_(&n, &V[j * ldv], &ione);

            // Extended local reorthogonalization
            if ((j > 0) && (elr > 0) && (alpha < kappa * beta))
            {
                for (i = 0; i < elr; i++)
                {
                    s = zdotc_(&n, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
                    PROPACK_CPLX_TYPE neg_s = PROPACK_cplx(-creal(s), -cimag(s));
                    zaxpy_(&n, &neg_s, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
                    nrm = dznrm2_(&n, &V[j * ldv], &ione);
                    if (nrm >= kappa * alpha) { break; }
                    alpha = nrm;
                }
                dwork[inu + j - 1] = eps;
                alpha = nrm;
            }

            B[j] = alpha;
            amax = fmax(amax, alpha);

            // Update estimate of ||A||_2
            if (j == 1)
            {
                a1 = B[ldb] / amax;
                a1 = FUDGE * amax * sqrt(pow(B[0] / amax, 2) + a1 * a1 + B[1] / amax * a1);
            } else {
                a1 = B[j - 1] / amax;
                b1 = B[j - 1 + ldb] / amax;
                a1 = FUDGE * amax * sqrt(a1 * a1 + b1 * b1 + a1 * B[j - 2 + ldb] / amax + B[j] / amax * b1);
            }
            anorm = fmax(anorm, a1);
        }

        // Update the nu recurrence
        if ((!full_reorth) && (alpha != 0.0))
        {
            dupdate_nu(&numax, &dwork[imu], &dwork[inu], j, &B[0], &B[ldb], anorm, epsn2);
        }

        // Reorthogonalize if necessary
        if ((full_reorth || numax > delta || force_reorth) && alpha != 0.0)
        {
            if ((full_reorth) || (eta == 0.0))
            {
                iwork[iidx] = 0;
                iwork[iidx + 1] = j - 1;
                iwork[iidx + 2] = j;
                iwork[iidx + 3] = j;
            } else if (!force_reorth) {
                dcompute_int(&dwork[inu], j - 1, delta, eta, &iwork[iidx]);
            }

            zreorth(n, j - 1, V, ldv, &V[j * ldv], &alpha, &iwork[iidx], kappa, &zwork[is], cgs);
            dset_mu(j - 1, &dwork[inu], &iwork[iidx], eps);
            numax = eta;

            force_reorth = (force_reorth ? 0 : 1);
        }

        // Check whether an invariant subspace was found
        if ((alpha < anorm * epsn) && (j < *k - 1))
        {
            *rnorm = alpha;
            alpha = 0.0;

            // Try to build an orthogonal subspace starting with a random vector
            zgetu0(1, m, n, j - 1, 3, &V[j * ldv], &alpha, V, ldv, aprod, zparm, iparm, ierr, cgs, &anormest, &zwork[is], rng_state);

            if (alpha == 0.0)
            {
                // We failed to generate a new random vector
                // in span(A^H) orthogonal to span(V(:,1:j-1)).
                // Most likely span(V(:,1:j-1)) is an invariant
                // subspace.
                *k = j;
                *ierr = -j-1;
                doption[2] = anorm;
                return;
            } else {
                // We have managed to generate a random vector
                // in span(A^H) orthogonal to V(:,1:j-1), so we
                // can continue the LBD and "deflate" the subspace
                // by setting alpha_{j} = 0.

                zsafescal(n, alpha, &V[j * ldv]);
                alpha = 0.0;
                force_reorth = 1;
                if (delta > 0.0) { full_reorth = 0; }
            }

        } else if ((j > 0) && (!full_reorth) && (j < *k - 1) && (delta * alpha < anorm * eps)) {
            *ierr = j;
        }
        B[j] = alpha;

        if (alpha != 0.0)
        {
            zsafescal(n, alpha, &V[j * ldv]);
        }

        // beta_{j+1} * u_{j+1} = A * v_j - alpha_j * u_j
        aprod(0, m, n, &V[j * ldv], &U[(j + 1) * ldu], zparm, iparm);

        PROPACK_CPLX_TYPE neg_alpha = PROPACK_cplx(-alpha, 0.0);
        zaxpy_(&m, &neg_alpha, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
        beta = dznrm2_(&m, &U[(j + 1) * ldu], &ione);

        // Extended local reorthogonalization
        if ((elr > 0) && (beta < kappa * alpha))
        {
            for (i = 0; i < elr; i++)
            {
                s = zdotc_(&m, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
                PROPACK_CPLX_TYPE neg_s = PROPACK_cplx(-creal(s), -cimag(s));
                zaxpy_(&m, &neg_s, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
                nrm = dznrm2_(&m, &U[(j + 1) * ldu], &ione);
                if (nrm >= kappa * beta) { break; }
                beta = nrm;
            }
            dwork[imu + j] = eps;
            beta = nrm;
        }

        B[j + ldb] = beta;
        amax = fmax(amax, beta);

        // Update estimate of ||A||_2
        if (j <= 0)
        {
            a1 = hypot(B[0], B[ldb]);
        } else {
            a1 = B[j] / amax;
            a1 = amax * sqrt(a1 * a1 + pow(B[j + ldb] / amax, 2) + a1 * B[j - 1 + ldb] / amax);
        }
        anorm = fmax(anorm, a1);

        // Update the mu recurrence
        if ((!full_reorth) && (beta != 0.0))
        {
            dupdate_mu(&mumax, &dwork[imu], &dwork[inu], j, &B[0], &B[ldb], anorm, epsn2);
        }

        // Reorthogonalize u_{j+1} if necessary
        if ((full_reorth || (mumax > delta) || force_reorth) && (beta != 0.0))
        {
            if ((full_reorth) || (eta == 0.0))
            {
                iwork[iidx] = 0;
                iwork[iidx + 1] = j;
                iwork[iidx + 2] = j + 1;
                iwork[iidx + 3] = j + 1;
            } else if (!force_reorth) {
                dcompute_int(&dwork[imu], j, delta, eta, &iwork[iidx]);
            } else {
                for (i = 0; i < 2 * j + 1; i++)
                {
                    if (iwork[iidx + i] == j)
                    {
                        iwork[iidx + i] = j + 1;
                        iwork[iidx + i + 1] = j + 1;
                        break;
                    }
                }
            }

            zreorth(m, j, U, ldu, &U[(j + 1) * ldu], &beta, &iwork[iidx], kappa, &zwork[is], cgs);
            dset_mu(j, &dwork[imu], &iwork[iidx], eps);
            mumax = eta;

            force_reorth = (force_reorth ? 0 : 1);
        }

        // Check whether an invariant subspace was found
        if ((beta < anorm * epsn) && (j < *k - 1))
        {
            *rnorm = beta;
            beta = 0.0;

            // Try to build an orthogonal subspace starting with a random vector
            zgetu0(0, m, n, j, 3, &U[(j + 1) * ldu], &beta, U, ldu, aprod, zparm, iparm, ierr, cgs, &anormest, &zwork[is], rng_state);

            if (beta == 0.0)
            {
                // We failed to generate a new random vector
                // in span(A) orthogonal to span(U(:,1:j)).
                // Most likely span(U(:,1:j)) is an invariant
                // subspace.
                *k = j+1;
                *ierr = -j-1;
                doption[2] = anorm;
                return;
            } else {
                // We have managed to generate a random vector
                // in span(A) orthogonal to U(:,1:j), so we can
                // continue the LBD and "deflate" the subspace by
                // setting beta_{j+1} = 0.
                zsafescal(m, beta, &U[(j + 1) * ldu]);
                beta = 0.0;
                force_reorth = 1;
                if (delta > 0.0) { full_reorth = 0; }
            }
        } else if (!full_reorth && (j < *k - 1) && (delta * beta < anorm * eps)) {
            *ierr = j;
        }

        B[j + ldb] = beta;
        if ((beta != 0.0) && (beta != 1.0))
        {
            zsafescal(m, beta, &U[(j + 1) * ldu]);
        }
        *rnorm = beta;
    }

    doption[2] = anorm;
}


