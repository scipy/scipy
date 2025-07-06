#include "lanbpro.h"
#include "blaslapack_declarations.h"
#include "common.h"
#include "getu0.h"
#include "gs.h"
#include <math.h>
#include <float.h>


static void ssafescal(int n, float alpha, float* x) {
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


void dsafescal(int n, double alpha, double* x) {
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


void csafescal(int n, float alpha, PROPACK_CPLXF_TYPE* x) {
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


void zsafescal(int n, double alpha, PROPACK_CPLX_TYPE* x) {
    // Scale the complex vector x by 1/alpha avoiding unnecessary under- and overflow.
    static double sfmin = DBL_MIN;
    int ione = 1;
    int info;

    if (fabs(alpha) >= sfmin) {
        double inv_alpha = 1.0 / alpha;
        zscal_(&n, &inv_alpha, x, &ione);
    } else {
        // Use LAPACK's safe scaling for very small alpha values
        double one = 1.0;
        zlascl_("G", &ione, &ione, &alpha, &one, &n, &ione, x, &n, &info);
    }
}



void slanbpro(int m, int n, int k0, int* k, PROPACK_aprod_s aprod,
              float* U, int ldu, float* V, int ldv, float* B, int ldb,
              float* rnorm, float* doption, int* ioption, float* work, int* iwork,
              float* dparm, int* iparm, int* ierr, uint64_t* rng_state) {

    // Constants
    const float one = 1.0f;
    const float zero = 0.0f;
    const float FUDGE = 1.01f;
    const float kappa = sqrtf(2.0) / 2.0f;

    // Machine constants
    const float eps = FLT_EPSILON;
    const float eps34 = powf(eps, 0.75f);
    const float epsn = (float)fmaxf(m, n) * eps;
    const float epsn2 = sqrtf((float)fmaxf(m, n)) * eps;
    const float sfmin = FLT_MIN;

    // Local variables
    int i, j, inu, imu, is, iidx, j0;
    float s, mumax, numax, alpha, beta, a1, b1, amax, anormest;
    int force_reorth, full_reorth;
    int ione = 1;

    // Set default parameters
    float delta, eta, anorm;
    int cgs, elr;

    if (doption[0] < zero) {
        delta = sqrtf(eps / *k);
    } else {
        delta = doption[0];
    }
    if (doption[1] < zero) {
        eta = eps34 / sqrtf((float)*k);
    } else {
        eta = doption[1];
    }
    if (delta <= eta || delta == zero) {
        full_reorth = 1;
    } else {
        full_reorth = 0;
    }
    if (doption[2] > zero) {
        anorm = doption[2];
    } else if (k0 > 0) {
        anorm = hypotf(B[0], B[ldb]);
        if (anorm <= zero) {
            *ierr = -1;
            return;
        }
    } else {
        anorm = zero;
    }

    cgs = ioption[0];
    elr = ioption[1];
    *ierr = 0;

    // Get starting vector if needed
    if (*rnorm == zero) {
        sgetu0(0, m, n, k0, 3, &U[k0 * ldu], rnorm, U, ldu, aprod,
               dparm, iparm, ierr, cgs, &anormest, work, rng_state);
        anorm = fmaxf(anorm, anormest);
    }

    // Set pointers into work array
    // work layout: [mu(k+1), nu(k+1), workspace...]
    imu = 0;
    inu = *k + 1;
    is = 2 * (*k + 1);
    iidx = 0;

    // Initialize work arrays
    for (i = 0; i < fmaxf(m, n) + 2 * (*k); i++) { work[i] = zero; }
    for (i = 0; i < 2 * (*k) + 1; i++) { iwork[i] = 0; }

    // Prepare Lanczos iteration
    if (k0 == 0) {
        amax = zero;
        alpha = zero;
        beta = *rnorm;
        force_reorth = 0;

        // Compute ||A x|| / ||x|| for a random vector x to estimate ||A||
        if (n > m) {
            sgetu0(0, m, n, 0, 1, &work[is], &s, U, ldu, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is + m], rng_state);
        } else {
            sgetu0(1, m, n, 0, 1, &work[is], &s, V, ldv, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is + n], rng_state);
        }
        anorm = fmaxf(anorm, FUDGE * anormest);
        j0 = 0;

        if (beta != zero) {
            ssafescal(m, beta, &U[0]);
        }
        work[imu] = one;
        work[inu] = one;
    } else {
        force_reorth = 1;
        alpha = B[k0 - 1];
        beta = *rnorm;

        if (k0 < *k && beta * delta < anorm * eps) {
            full_reorth = 1;
            *ierr = k0;
        }

        iwork[iidx] = 0;
        iwork[iidx + 1] = k0 - 1;

        sscal_(&m, rnorm, &U[k0 * ldu], &ione);
        sreorth(m, k0, U, ldu, &U[k0 * ldu], rnorm, &iwork[iidx], kappa, &work[is], cgs);
        if (*rnorm != zero) {
            s = one / *rnorm;
            sscal_(&m, &s, &U[k0 * ldu], &ione);
        }
        sset_mu(k0, &work[imu], &iwork[iidx], epsn2);
        sset_mu(k0, &work[inu], &iwork[iidx], epsn2);
        beta = *rnorm;

        // Estimate ||B||_2^2 as ||B^T * B||_1
        B[k0 - 1 + ldb] = beta;
        amax = zero;
        for (j = 0; j < k0; j++) {
            amax = fmaxf(amax, fmaxf(B[j], B[j + ldb]));
            if (j == 0) {
                anorm = fmaxf(anorm, FUDGE * alpha);
            } else if (j == 1) {
                a1 = B[ldb] / amax;
                a1 = FUDGE * amax * sqrtf(powf(B[0] / amax, 2) + a1 * a1 +
                                         B[1] / amax * a1);
                anorm = fmaxf(anorm, a1);
            } else {
                a1 = B[j - 1] / amax;
                b1 = B[j - 1 + ldb] / amax;
                a1 = FUDGE * amax * sqrtf(a1 * a1 + b1 * b1 +
                                         a1 * B[j - 2 + ldb] / amax +
                                         B[j] / amax * b1);
                anorm = fmaxf(anorm, a1);
            }
        }
        j0 = k0;
    }

    numax = zero;
    mumax = zero;

    // Main Lanczos bidiagonalization iteration
    for (j = j0; j < *k; j++) {

        // alpha_j * v_j = A^T * u_j - beta_j * v_{j-1}
        aprod(1, m, n, &U[j * ldu], &V[j * ldv], dparm, iparm);

        if (j == 0) {
            alpha = snrm2_(&n, &V[j * ldv], &ione);
            anorm = fmaxf(anorm, FUDGE * alpha);
        } else {
            float neg_beta = -beta;
            saxpy_(&n, &neg_beta, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
            alpha = snrm2_(&n, &V[j * ldv], &ione);

            // Extended local reorthogonalization
            if (j > 0 && elr > 0 && alpha < kappa * beta) {
                for (i = 0; i < elr; i++) {
                    s = sdot_(&n, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
                    float neg_s = -s;
                    saxpy_(&n, &neg_s, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
                    if (beta != zero) {
                        beta = beta + s;
                        B[j - 1 + ldb] = beta;
                    }
                    s = snrm2_(&n, &V[j * ldv], &ione);
                    if (s >= kappa * alpha) break;
                    alpha = s;
                }
                work[inu + j - 1] = eps;
                alpha = s;
            }

            B[j] = alpha;
            amax = fmaxf(amax, alpha);

            // Update estimate of ||A||_2
            if (j == 1) {
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
        if (!full_reorth && alpha != zero) {
            supdate_nu(&numax, &work[imu], &work[inu], j, &B[0], &B[ldb], anorm, epsn2);
        }

        // Reorthogonalize if necessary
        if ((full_reorth || numax > delta || force_reorth) && alpha != zero) {
            if (full_reorth || eta == zero) {
                iwork[iidx] = 0;
                iwork[iidx + 1] = j - 1;
            } else if (!force_reorth) {
                scompute_mu(&work[inu], j - 1, delta, eta, &iwork[iidx]);
            }

            sreorth(n, j - 1, V, ldv, &V[j * ldv], &alpha, &iwork[iidx], kappa, &work[is], cgs);

            sset_mu(j - 1, &work[inu], &iwork[iidx], eps);
            numax = eta;

            if (force_reorth) {
                force_reorth = 0;
            } else {
                force_reorth = 1;
            }
        }

        // Check whether an invariant subspace was found
        if (alpha < anorm * epsn && j < *k - 1) {
            *rnorm = alpha;
            alpha = zero;

            // Try to build an orthogonal subspace starting with a random vector
            sgetu0(1, m, n, j - 1, 3, &V[j * ldv], &alpha, V, ldv, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is], rng_state);

            if (alpha == zero) {
                // Failed to generate new random vector - invariant subspace found
                *k = j - 1;
                *ierr = -j;
                return;
            } else {
                // Successfully generated random vector - continue with deflation
                ssafescal(n, alpha, &V[j * ldv]);
                alpha = zero;
                force_reorth = 1;
                if (delta > zero) {
                    full_reorth = 0;
                }
            }
        } else if ((j > 0) && (!full_reorth) && (j < *k - 1) && (delta * alpha < anorm * eps)) {
            *ierr = j;
        }

        B[j] = alpha;

        if (alpha != zero) {
            ssafescal(n, alpha, &V[j * ldv]);
        }

        // beta_{j+1} * u_{j+1} = A * v_j - alpha_j * u_j
        aprod(0, m, n, &V[j * ldv], &U[(j + 1) * ldu], dparm, iparm);

        float neg_alpha = -alpha;
        saxpy_(&m, &neg_alpha, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
        beta = snrm2_(&m, &U[(j + 1) * ldu], &ione);

        // Extended local reorthogonalization
        if ((elr > 0) && (beta < kappa * alpha)) {
            for (i = 0; i < elr; i++) {
                s = sdot_(&m, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
                float neg_s = -s;
                saxpy_(&m, &neg_s, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
                if (alpha != zero) {
                    alpha = alpha + s;
                    B[j] = alpha;
                }
                s = snrm2_(&m, &U[(j + 1) * ldu], &ione);
                if (s >= kappa * beta) break;
                beta = s;
            }
            work[imu + j] = eps;
            beta = s;
        }

        B[j + ldb] = beta;
        amax = fmaxf(amax, beta);

        // Update estimate of ||A||_2
        if (j <= 0) {
            a1 = hypotf(B[0], B[ldb]);
        } else {
            a1 = B[j] / amax;
            a1 = amax * sqrtf(a1 * a1 + powf(B[j + ldb] / amax, 2) +
                             a1 * B[j - 1 + ldb] / amax);
        }
        anorm = fmaxf(anorm, a1);

        // Update the mu recurrence
        if (!full_reorth && beta != zero) {
            supdate_mu(&mumax, &work[imu], &work[inu], j, &B[0], &B[ldb],
                       anorm, epsn2);
        }

        // Reorthogonalize u_{j+1} if necessary
        if ((full_reorth || mumax > delta || force_reorth) && beta != zero) {
            if (full_reorth || eta == zero) {
                iwork[iidx] = 0;
                iwork[iidx + 1] = j;
            } else if (!force_reorth) {
                scompute_mu(&work[imu], j, delta, eta, &iwork[iidx]);
            } else {
                for (i = 0; i < 2 * j + 1; i++) {
                    if (iwork[iidx + i] == j) {
                        iwork[iidx + i] = j + 1;
                        break;
                    }
                }
            }

            sreorth(m, j, U, ldu, &U[(j + 1) * ldu], &beta, &iwork[iidx],
                    kappa, &work[is], cgs);

            sset_mu(j, &work[imu], &iwork[iidx], eps);
            mumax = eta;

            if (force_reorth) {
                force_reorth = 0;
            } else {
                force_reorth = 1;
            }
        }

        // Check whether an invariant subspace was found
        if (beta < anorm * epsn && j < *k - 1) {
            *rnorm = beta;
            beta = zero;

            // Try to build an orthogonal subspace starting with a random vector
            sgetu0(0, m, n, j, 3, &U[(j + 1) * ldu], &beta, U, ldu, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is], rng_state);

            if (beta == zero) {
                // Failed to generate new random vector - invariant subspace found
                *k = j;
                *ierr = -j;
                return;
            } else {
                // Successfully generated random vector - continue with deflation
                ssafescal(m, beta, &U[(j + 1) * ldu]);
                beta = zero;
                force_reorth = 1;
                if (delta > zero) {
                    full_reorth = 0;
                }
            }
        } else if (!full_reorth && j < *k - 1 &&
                   (delta * beta < anorm * eps)) {
            *ierr = j;
        }

        B[j + ldb] = beta;
        if (beta != zero && beta != one) {
            ssafescal(m, beta, &U[(j + 1) * ldu]);
        }
        *rnorm = beta;
    }

    *ierr = 0;
}

void dlanbpro(int m, int n, int k0, int* k, PROPACK_aprod_d aprod,
              double* U, int ldu, double* V, int ldv, double* B, int ldb,
              double* rnorm, double* doption, int* ioption, double* work, int* iwork,
              double* dparm, int* iparm, int* ierr, uint64_t* rng_state) {

    // Constants
    const double one = 1.0;
    const double zero = 0.0;
    const double FUDGE = 1.01;
    const double kappa = sqrt(2.0) / 2.0;

    // Machine constants
    const double eps = DBL_EPSILON;
    const double eps34 = pow(eps, 0.75);
    const double epsn = (double)fmax(m, n) * eps;
    const double epsn2 = sqrt((double)fmax(m, n)) * eps;
    const double sfmin = DBL_MIN;

    // Local variables
    int i, j, inu, imu, is, iidx, j0;
    double s, mumax, numax, alpha, beta, a1, b1, amax, anormest;
    int force_reorth, full_reorth;
    int ione = 1;

    // Set default parameters
    double delta, eta, anorm;
    int cgs, elr;

    if (doption[0] < zero) {
        delta = sqrt(eps / *k);
    } else {
        delta = doption[0];
    }
    if (doption[1] < zero) {
        eta = eps34 / sqrt((double)*k);
    } else {
        eta = doption[1];
    }
    if (delta <= eta || delta == zero) {
        full_reorth = 1;
    } else {
        full_reorth = 0;
    }
    if (doption[2] > zero) {
        anorm = doption[2];
    } else if (k0 > 0) {
        anorm = hypot(B[0], B[ldb]);
        if (anorm <= zero) {
            *ierr = -1;
            return;
        }
    } else {
        anorm = zero;
    }

    cgs = ioption[0];
    elr = ioption[1];
    *ierr = 0;

    // Get starting vector if needed
    if (*rnorm == zero) {
        dgetu0(0, m, n, k0, 3, &U[k0 * ldu], rnorm, U, ldu, aprod,
               dparm, iparm, ierr, cgs, &anormest, work, rng_state);
        anorm = fmax(anorm, anormest);
    }

    // Set pointers into work array
    // work layout: [mu(k+1), nu(k+1), workspace...]
    imu = 0;
    inu = *k + 1;
    is = 2 * (*k + 1);
    iidx = 0;

    // Initialize work arrays
    for (i = 0; i < fmax(m, n) + 2 * (*k); i++) { work[i] = zero; }
    for (i = 0; i < 2 * (*k) + 1; i++) { iwork[i] = 0; }

    // Prepare Lanczos iteration
    if (k0 == 0) {
        amax = zero;
        alpha = zero;
        beta = *rnorm;
        force_reorth = 0;

        // Compute ||A x|| / ||x|| for a random vector x to estimate ||A||
        if (n > m) {
            dgetu0(0, m, n, 0, 1, &work[is], &s, U, ldu, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is + m], rng_state);
        } else {
            dgetu0(1, m, n, 0, 1, &work[is], &s, V, ldv, aprod, dparm, iparm, ierr, cgs, &anormest, &work[is + n], rng_state);
        }
        anorm = fmax(anorm, FUDGE * anormest);
        j0 = 0;

        if (beta != zero) {
            dsafescal(m, beta, &U[0]);
        }
        work[imu] = one;
        work[inu] = one;
    } else {
        force_reorth = 1;
        alpha = B[k0 - 1];
        beta = *rnorm;

        if (k0 < *k && beta * delta < anorm * eps) {
            full_reorth = 1;
            *ierr = k0;
        }

        iwork[iidx] = 0;
        iwork[iidx + 1] = k0 - 1;

        dscal_(&m, rnorm, &U[k0 * ldu], &ione);
        dreorth(m, k0, U, ldu, &U[k0 * ldu], rnorm, &iwork[iidx], kappa, &work[is], cgs);
        if (*rnorm != zero) {
            s = one / *rnorm;
            dscal_(&m, &s, &U[k0 * ldu], &ione);
        }
        dset_mu(k0, &work[imu], &iwork[iidx], epsn2);
        dset_mu(k0, &work[inu], &iwork[iidx], epsn2);
        beta = *rnorm;

        // Estimate ||B||_2^2 as ||B^T * B||_1
        B[k0 - 1 + ldb] = beta;
        amax = zero;
        for (j = 0; j < k0; j++) {
            amax = fmax(amax, fmax(B[j], B[j + ldb]));
            if (j == 0) {
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

    numax = zero;
    mumax = zero;

    // Main Lanczos bidiagonalization iteration
    for (j = j0; j < *k; j++) {

        // alpha_j * v_j = A^T * u_j - beta_j * v_{j-1}
        aprod(1, m, n, &U[j * ldu], &V[j * ldv], dparm, iparm);

        if (j == 0) {
            alpha = dnrm2_(&n, &V[j * ldv], &ione);
            anorm = fmax(anorm, FUDGE * alpha);
        } else {
            double neg_beta = -beta;
            daxpy_(&n, &neg_beta, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
            alpha = dnrm2_(&n, &V[j * ldv], &ione);

            // Extended local reorthogonalization
            if (j > 0 && elr > 0 && alpha < kappa * beta) {
                for (i = 0; i < elr; i++) {
                    s = ddot_(&n, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
                    double neg_s = -s;
                    daxpy_(&n, &neg_s, &V[(j - 1) * ldv], &ione, &V[j * ldv], &ione);
                    if (beta != zero) {
                        beta = beta + s;
                        B[j - 1 + ldb] = beta;
                    }
                    s = dnrm2_(&n, &V[j * ldv], &ione);
                    if (s >= kappa * alpha) break;
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
        if (!full_reorth && alpha != zero) {
            dupdate_nu(&numax, &work[imu], &work[inu], j, &B[0], &B[ldb], anorm, epsn2);
        }

        // Reorthogonalize if necessary
        if ((full_reorth || numax > delta || force_reorth) && alpha != zero) {
            if (full_reorth || eta == zero) {
                iwork[iidx] = 0;
                iwork[iidx + 1] = j - 1;
            } else if (!force_reorth) {
                dcompute_mu(&work[inu], j - 1, delta, eta, &iwork[iidx]);
            }

            dreorth(n, j - 1, V, ldv, &V[j * ldv], &alpha, &iwork[iidx], kappa, &work[is], cgs);

            dset_mu(j - 1, &work[inu], &iwork[iidx], eps);
            numax = eta;

            if (force_reorth) {
                force_reorth = 0;
            } else {
                force_reorth = 1;
            }
        }

        // Check whether an invariant subspace was found
        if (alpha < anorm * epsn && j < *k - 1) {
            *rnorm = alpha;
            alpha = zero;

            // Try to build an orthogonal subspace starting with a random vector
            dgetu0(1, m, n, j - 1, 3, &V[j * ldv], &alpha, V, ldv, aprod,
                   dparm, iparm, ierr, cgs, &anormest, &work[is], rng_state);

            if (alpha == zero) {
                // Failed to generate new random vector - invariant subspace found
                *k = j - 1;
                *ierr = -j;
                return;
            } else {
                // Successfully generated random vector - continue with deflation
                dsafescal(n, alpha, &V[j * ldv]);
                alpha = zero;
                force_reorth = 1;
                if (delta > zero) {
                    full_reorth = 0;
                }
            }
        } else if ((j > 0) && (!full_reorth) && (j < *k - 1) && (delta * alpha < anorm * eps)) {
            *ierr = j;
        }

        B[j] = alpha;

        if (alpha != zero) {
            dsafescal(n, alpha, &V[j * ldv]);
        }

        // beta_{j+1} * u_{j+1} = A * v_j - alpha_j * u_j
        aprod(0, m, n, &V[j * ldv], &U[(j + 1) * ldu], dparm, iparm);

        double neg_alpha = -alpha;
        daxpy_(&m, &neg_alpha, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
        beta = dnrm2_(&m, &U[(j + 1) * ldu], &ione);

        // Extended local reorthogonalization
        if ((elr > 0) && (beta < kappa * alpha)) {
            for (i = 0; i < elr; i++) {
                s = ddot_(&m, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
                double neg_s = -s;
                daxpy_(&m, &neg_s, &U[j * ldu], &ione, &U[(j + 1) * ldu], &ione);
                if (alpha != zero) {
                    alpha = alpha + s;
                    B[j] = alpha;
                }
                s = dnrm2_(&m, &U[(j + 1) * ldu], &ione);
                if (s >= kappa * beta) break;
                beta = s;
            }
            work[imu + j] = eps;
            beta = s;
        }

        B[j + ldb] = beta;
        amax = fmax(amax, beta);

        // Update estimate of ||A||_2
        if (j <= 0) {
            a1 = hypot(B[0], B[ldb]);
        } else {
            a1 = B[j] / amax;
            a1 = amax * sqrt(a1 * a1 + pow(B[j + ldb] / amax, 2) +
                             a1 * B[j - 1 + ldb] / amax);
        }
        anorm = fmax(anorm, a1);

        // Update the mu recurrence
        if (!full_reorth && beta != zero) {
            dupdate_mu(&mumax, &work[imu], &work[inu], j, &B[0], &B[ldb],
                       anorm, epsn2);
        }

        // Reorthogonalize u_{j+1} if necessary
        if ((full_reorth || mumax > delta || force_reorth) && beta != zero) {
            if (full_reorth || eta == zero) {
                iwork[iidx] = 1;
                iwork[iidx + 1] = j;
                iwork[iidx + 2] = j + 1;
            } else if (!force_reorth) {
                dcompute_mu(&work[imu], j, delta, eta, &iwork[iidx]);
            } else {
                for (i = 0; i < 2 * j + 1; i++) {
                    if (iwork[iidx + i] == j) {
                        iwork[iidx + i] = j + 1;
                        break;
                    }
                }
            }

            dreorth(m, j, U, ldu, &U[(j + 1) * ldu], &beta, &iwork[iidx],
                    kappa, &work[is], cgs);

            dset_mu(j, &work[imu], &iwork[iidx], eps);
            mumax = eta;

            if (force_reorth) {
                force_reorth = 0;
            } else {
                force_reorth = 1;
            }
        }

        // Check whether an invariant subspace was found
        if (beta < anorm * epsn && j < *k - 1) {
            *rnorm = beta;
            beta = zero;

            // Try to build an orthogonal subspace starting with a random vector
            dgetu0(0, m, n, j, 3, &U[(j + 1) * ldu], &beta, U, ldu, aprod,
                   dparm, iparm, ierr, cgs, &anormest, &work[is], rng_state);

            if (beta == zero) {
                // Failed to generate new random vector - invariant subspace found
                *k = j;
                *ierr = -j;
                return;
            } else {
                // Successfully generated random vector - continue with deflation
                dsafescal(m, beta, &U[(j + 1) * ldu]);
                beta = zero;
                force_reorth = 1;
                if (delta > zero) {
                    full_reorth = 0;
                }
            }
        } else if (!full_reorth && j < *k - 1 &&
                   (delta * beta < anorm * eps)) {
            *ierr = j;
        }

        B[j + ldb] = beta;
        if (beta != zero && beta != one) {
            dsafescal(m, beta, &U[(j + 1) * ldu]);
        }
        *rnorm = beta;
    }

    *ierr = 0;
}


