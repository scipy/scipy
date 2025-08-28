#include "../include/propack/propack.h"
#include <float.h>
#include <math.h>

#include "common.h"
#include "getu0.h"
#include "lanbpro.h"
#include "ritzvec.h"
#include "blaslapack_declarations.h"


static inline int int_min(const int a, const int b) { return a < b ? a : b; }
static inline int int_max(const int a, const int b) { return a > b ? a : b; }


void slansvd_irl(int which, int jobu, int jobv, int m, int n, int dim, int p, int *neig, int maxiter,
                 PROPACK_aprod_s aprod, float* U, int ldu, float* sigma, float* bnd, float* V, int ldv,
                 float tolin, float* work, int lwork, int* iwork, float* doption, int* ioption,
                 int* info, float* dparm, int* iparm, uint64_t* rng_state)
{
    // Parameters
    int int1 = 1, int0 = 0;

    // Local variables
    int k, i, ibnd, iwrk, ierr, ip, iq, nconv, lwrk, kold;
    int ialpha, ibeta, ialpha1, ibeta1, ishift, nshft, lapinfo;
    float eps, eps34, epsn2, epsn, sfmin, anorm, rnorm, tol;
    float shift, relgap;
    int iter;

    // Set machine dependent constants
    eps = FLT_EPSILON;
    eps34 = powf(eps, 3.0f/4.0f);
    epsn = (float)int_max(m, n) * eps / 2.0f;
    epsn2 = sqrtf((float)int_max(m, n)) * eps / 2.0f;
    sfmin = FLT_MIN;

    // Guard against absurd arguments
    dim = int_min(int_min(dim, n + 1), m + 1);
    k = dim - p;
    tol = fminf(1.0f, fmaxf(16.0f * eps, tolin));
    anorm = 0.0f;

    // Set pointers into work array
    ibnd = 0;
    ialpha = ibnd + dim + 1;
    ibeta = ialpha + dim;
    ialpha1 = ibeta + dim;
    ibeta1 = ialpha1 + dim;
    ishift = ibeta1 + dim;
    ip = ishift + dim;
    iq = ip + (dim + 1) * (dim + 1);
    iwrk = iq + dim * dim;
    lwrk = lwork - iwrk;

    // Zero out work array sections
    for (i = 0; i < 8 * dim + 3 + 2 * dim * dim; i++) { work[i] = 0.0f; }

    // Set up random starting vector if none is provided by the user
    rnorm = snrm2_(&m, &U[0], &int1);
    if (rnorm == 0.0f)
    {
        sgetu0(0, m, n, 0, 1, &U[0], &rnorm, U, ldu, aprod, dparm, iparm, &ierr, ioption[0], &anorm, &work[iwrk], rng_state);
    }

    iter = 0;
    *info = 0;
    nconv = 0;
    kold = 0;

    // Iterate until convergence...
    while ((nconv < *neig) && (iter < maxiter))
    {
        // Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j}
        slanbpro(m, n, kold, &dim, aprod, U, ldu, V, ldv, &work[ialpha], dim, &rnorm,
                 &doption[0], &ioption[0], &work[iwrk], iwork, dparm, iparm, &ierr, rng_state);
        kold = k;

        // Compute and analyze SVD(B) and error bounds
        scopy_(&dim, &work[ialpha], &int1, &work[ialpha1], &int1);
        scopy_(&dim, &work[ibeta], &int1, &work[ibeta1], &int1);

        // Zero out bounds array
        for (i = 0; i < dim + 1; i++) { work[ibnd + i] = 0.0f; }

        // QR factorization of bidiagonal matrix
        sbdqr((dim == int_min(m, n)), 0, dim, &work[ialpha1], &work[ibeta1],
              &work[ibnd + dim - 1], &work[ibnd + dim], &work[ip], dim + 1);

        // SVD of bidiagonal matrix
        sbdsqr_("U", &dim, &int0, &int1, &int0, &work[ialpha1], &work[ibeta1], work, &int1,
                &work[ibnd], &int1, work, &int1, &work[iwrk], &lapinfo);

        // Update anorm estimate
        if (dim > 5)
        {
            anorm = work[ialpha1];
        } else {
            anorm = fmaxf(anorm, work[ialpha1]);
        }

        // Scale error bounds by rnorm
        for (i = 0; i < dim; i++)
        {
            work[ibnd + i] = fabsf(rnorm * work[ibnd + i]);
        }

        // Refine error bounds using the "Gap theorem"
        if (which == 0)  // which == 'S' (smallest)
        {
            srefinebounds(int_min(m, n), dim, &work[ialpha1], &work[ibnd], epsn * anorm, eps34);
        } else {  // which == 'L' (largest)
            srefinebounds(int_min(m, n), int_min(dim, *neig), &work[ialpha1], &work[ibnd], epsn * anorm, eps34);
        }

        // Determine the number of converged singular values
        if (which == 0)  // smallest singular values
        {
            i = dim - *neig;
            nconv = 0;
            while (i < dim)
            {
                if (work[ibnd + i] <= tol * work[ialpha1])
                {
                    nconv++;
                    sigma[nconv - 1] = work[ialpha1 + i];
                    bnd[nconv - 1] = work[ibnd + i];
                }
                i++;
            }
        } else {  // largest singular values
            i = 0;
            nconv = 0;
            while (i < int_min(dim, *neig))
            {
                if (work[ibnd + i] <= tol * work[ialpha1 + i])
                {
                    sigma[nconv] = work[ialpha1 + i];
                    bnd[nconv] = work[ibnd + i];
                    nconv++;
                    i++;
                } else {
                    i = k + 1;
                }
            }
        }

        // Test if an invariant subspace has been found
        if (ierr < 0)
        {
            if (dim < k)
            {
                *info = dim;
            }
            break;
        }

        if (nconv < *neig)
        {
            // Implicit restart:
            // Apply shifts mu_1, mu_2,...,mu_p to "back up" the
            // bidiagonalization (dim-k) steps to
            //
            //      A V_{k}^{+} = U_{k+1}^{+} B_{k}^{+}
            // corresponding to starting vector
            //      u_1^{+} = \prod_{i=1}^{p} (A A^T - mu_i^2) u_1
            //
            // We use exact shifts mu_i for which the relative gap between mu_i
            // and the lower bound on the k'th Ritzvalue is larger than doption[3]

            // Zero out shift array
            for (i = 0; i < dim - k; i++) { work[ishift + i] = 0.0f; }

            nshft = 0;
            if (which == 0)  // smallest
            {
                for (i = 0; i < k; i++)
                {
                    relgap = (work[ialpha1 + i] - work[ibnd + i] - work[ialpha1 + dim - *neig - 1]);
                    if (relgap > doption[3] * work[ialpha1 + dim - *neig - 1])
                    {
                        work[ishift + nshft] = work[ialpha1 + i];
                    } else {
                        work[ishift + nshft] = work[ialpha1];
                    }
                    nshft++;
                }
            } else {  // largest
                for (i = dim - 1; i >= k; i--)
                {
                    relgap = work[ialpha1 + k - 1] - (work[ialpha1 + i] + work[ibnd + i]);
                    if (relgap > doption[3] * work[ialpha1 + k - 1])
                    {
                        work[ishift + nshft] = work[ialpha1 + i];
                    } else {
                        work[ishift + nshft] = 0.0f;
                    }
                    nshft++;
                }
            }

            // Apply shifts and accumulate rotations such that
            //   B_{dim}^{+} = P * B_{dim} * Q^T

            // Zero out rotation matrices and set to identity
            for (i = 0; i < (dim + 1) * (dim + 1); i++) { work[ip + i] = 0.0f; }
            for (i = 0; i < dim * dim; i++) { work[iq + i] = 0.0f; }

            // Set diagonal elements to 1 (identity matrices)
            for (i = 0; i < dim + 1; i++)
            {
                work[ip + i * (dim + 2)] = 1.0f;
            }
            for (i = 0; i < dim; i++)
            {
                work[iq + i * (dim + 1)] = 1.0f;
            }

            // Apply shifts
            for (i = dim - 1; i >= k; i--)
            {
                shift = work[ishift + dim - 1 - i];
                sbsvdstep(1, 1, dim + 1, dim, i + 1, shift, &work[ialpha], &work[ibeta],
                         &work[ip], dim + 1, &work[iq], dim);
            }

            // Compute first k+1 left and first k right updated Lanczos vectors
            //   U_{dim+1}^{+} = U_{dim+1} * P(:,1:k+1)
            //   V_{dim}^{+} = V_{dim} * Q(:,1:k)

            sgemm_ovwr_left(0, m, k + 1, dim + 1, 1.0f, U, ldu, &work[ip], dim + 1, &work[iwrk], lwrk / (k + 1));
            sgemm_ovwr_left(0, n, k, dim, 1.0f, V, ldv, &work[iq], dim, &work[iwrk], lwrk / k);

            rnorm = work[ibeta + k - 1];
        }
        iter++;
    }

    // Calculate singular vectors if requested
    if ((nconv >= *neig || *info > 0) && (jobu || jobv))
    {
        scopy_(&dim, &work[ialpha], &int1, &work[ialpha1], &int1);
        scopy_(&dim, &work[ibeta], &int1, &work[ibeta1], &int1);
        lwrk = lwrk + dim * dim + (dim + 1) * (dim + 1);
        sritzvec(which, jobu, jobv, m, n, nconv, dim, &work[ialpha1], &work[ibeta1], U, ldu, V, ldv, &work[ip], lwrk, iwork);
    }

    *neig = nconv;
}


void dlansvd_irl(int which, int jobu, int jobv, int m, int n, int dim, int p, int *neig, int maxiter,
                 PROPACK_aprod_d aprod, double* U, int ldu, double* sigma, double* bnd, double* V, int ldv,
                 double tolin, double* work, int lwork, int* iwork, double* doption, int* ioption,
                 int* info, double* dparm, int* iparm, uint64_t* rng_state)
{
    // Parameters
    int int1 = 1, int0 = 0;

    // Local variables
    int k, i, ibnd, iwrk, ierr, ip, iq, nconv, lwrk, kold;
    int ialpha, ibeta, ialpha1, ibeta1, ishift, nshft, lapinfo;
    double eps, eps34, epsn2, epsn, sfmin, anorm, rnorm, tol;
    double shift, relgap;
    int iter;

    // Set machine dependent constants
    eps = DBL_EPSILON;
    eps34 = pow(eps, 3.0/4.0);
    epsn = (double)fmax(m, n) * eps / 2.0;
    epsn2 = sqrt((double)fmax(m, n)) * eps / 2.0;
    sfmin = DBL_MIN;

    // Guard against absurd arguments
    dim = int_min(int_min(dim, n + 1), m + 1);
    k = dim - p;
    tol = fmin(1.0, fmax(16.0 * eps, tolin));
    anorm = 0.0;

    // Set pointers into work array
    ibnd = 0;
    ialpha = ibnd + dim + 1;
    ibeta = ialpha + dim;
    ialpha1 = ibeta + dim;
    ibeta1 = ialpha1 + dim;
    ishift = ibeta1 + dim;
    ip = ishift + dim;
    iq = ip + (dim + 1) * (dim + 1);
    iwrk = iq + dim * dim;
    lwrk = lwork - iwrk;

    // Zero out work array sections
    for (i = 0; i < 8 * dim + 3 + 2 * dim * dim; i++) { work[i] = 0.0; }

    // Set up random starting vector if none is provided by the user
    rnorm = dnrm2_(&m, &U[0], &int1);
    if (rnorm == 0.0)
    {
        dgetu0(0, m, n, 0, 1, &U[0], &rnorm, U, ldu, aprod, dparm, iparm, &ierr, ioption[0], &anorm, &work[iwrk], rng_state);
    }

    iter = 0;
    *info = 0;
    nconv = 0;
    kold = 0;

    // Iterate until convergence...
    while ((nconv < *neig) && (iter < maxiter))
    {
        // Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j}
        dlanbpro(m, n, kold, &dim, aprod, U, ldu, V, ldv, &work[ialpha], dim, &rnorm,
                 &doption[0], &ioption[0], &work[iwrk], iwork, dparm, iparm, &ierr, rng_state);
        kold = k;

        // Compute and analyze SVD(B) and error bounds
        dcopy_(&dim, &work[ialpha], &int1, &work[ialpha1], &int1);
        dcopy_(&dim, &work[ibeta], &int1, &work[ibeta1], &int1);

        // Zero out bounds array
        for (i = 0; i < dim + 1; i++) { work[ibnd + i] = 0.0; }

        // QR factorization of bidiagonal matrix
        dbdqr((dim == int_min(m, n)), 0, dim, &work[ialpha1], &work[ibeta1],
              &work[ibnd + dim - 1], &work[ibnd + dim], &work[ip], dim + 1);

        // SVD of bidiagonal matrix
        dbdsqr_("U", &dim, &int0, &int1, &int0, &work[ialpha1], &work[ibeta1], work, &int1,
                &work[ibnd], &int1, work, &int1, &work[iwrk], &lapinfo);

        // Update anorm estimate
        if (dim > 5)
        {
            anorm = work[ialpha1];
        } else {
            anorm = fmax(anorm, work[ialpha1]);
        }

        // Scale error bounds by rnorm
        for (i = 0; i < dim; i++)
        {
            work[ibnd + i] = fabs(rnorm * work[ibnd + i]);
        }

        // Refine error bounds using the "Gap theorem"
        if (which == 0)  // which == 'S' (smallest)
        {
            drefinebounds(int_min(m, n), dim, &work[ialpha1], &work[ibnd], epsn * anorm, eps34);
        } else {  // which == 'L' (largest)
            drefinebounds(int_min(m, n), int_min(dim, *neig), &work[ialpha1], &work[ibnd], epsn * anorm, eps34);
        }

        // Determine the number of converged singular values
        if (which == 0)  // smallest singular values
        {
            i = dim - *neig;
            nconv = 0;
            while (i < dim)
            {
                if (work[ibnd + i] <= tol * work[ialpha1])
                {
                    nconv++;
                    sigma[nconv - 1] = work[ialpha1 + i];
                    bnd[nconv - 1] = work[ibnd + i];
                }
                i++;
            }
        } else {  // largest singular values
            i = 0;
            nconv = 0;
            while (i < int_min(dim, *neig))
            {
                if (work[ibnd + i] <= tol * work[ialpha1 + i])
                {
                    sigma[nconv] = work[ialpha1 + i];
                    bnd[nconv] = work[ibnd + i];
                    nconv++;
                    i++;
                } else {
                    i = k + 1;
                }
            }
        }

        // Test if an invariant subspace has been found
        if (ierr < 0)
        {
            if (dim < k)
            {
                *info = dim;
            }
            break;
        }

        if (nconv < *neig)
        {
            // Implicit restart:
            // Apply shifts mu_1, mu_2,...,mu_p to "back up" the
            // bidiagonalization (dim-k) steps to
            //
            //      A V_{k}^{+} = U_{k+1}^{+} B_{k}^{+}
            // corresponding to starting vector
            //      u_1^{+} = \prod_{i=1}^{p} (A A^T - mu_i^2) u_1
            //
            // We use exact shifts mu_i for which the relative gap between mu_i
            // and the lower bound on the k'th Ritzvalue is larger than doption[3]

            // Zero out shift array
            for (i = 0; i < dim - k; i++) { work[ishift + i] = 0.0; }

            nshft = 0;
            if (which == 0)  // smallest
            {
                for (i = 0; i < k; i++)
                {
                    relgap = (work[ialpha1 + i] - work[ibnd + i] - work[ialpha1 + dim - *neig - 1]);
                    if (relgap > doption[3] * work[ialpha1 + dim - *neig - 1])
                    {
                        work[ishift + nshft] = work[ialpha1 + i];
                    } else {
                        work[ishift + nshft] = work[ialpha1];
                    }
                    nshft++;
                }
            } else {  // largest
                for (i = dim - 1; i >= k; i--)
                {
                    relgap = work[ialpha1 + k - 1] - (work[ialpha1 + i] + work[ibnd + i]);
                    if (relgap > doption[3] * work[ialpha1 + k - 1])
                    {
                        work[ishift + nshft] = work[ialpha1 + i];
                    } else {
                        work[ishift + nshft] = 0.0;
                    }
                    nshft++;
                }
            }

            // Apply shifts and accumulate rotations such that
            //   B_{dim}^{+} = P * B_{dim} * Q^T

            // Zero out rotation matrices and set to identity
            for (i = 0; i < (dim + 1) * (dim + 1); i++) { work[ip + i] = 0.0; }
            for (i = 0; i < dim * dim; i++) { work[iq + i] = 0.0; }

            // Set diagonal elements to 1 (identity matrices)
            for (i = 0; i < dim + 1; i++)
            {
                work[ip + i * (dim + 2)] = 1.0;
            }
            for (i = 0; i < dim; i++)
            {
                work[iq + i * (dim + 1)] = 1.0;
            }

            // Apply shifts
            for (i = dim - 1; i >= k; i--)
            {
                shift = work[ishift + dim - 1 - i];
                dbsvdstep(1, 1, dim + 1, dim, i + 1, shift, &work[ialpha], &work[ibeta],
                         &work[ip], dim + 1, &work[iq], dim);
            }

            // Compute first k+1 left and first k right updated Lanczos vectors
            //   U_{dim+1}^{+} = U_{dim+1} * P(:,1:k+1)
            //   V_{dim}^{+} = V_{dim} * Q(:,1:k)

            dgemm_ovwr_left(0, m, k + 1, dim + 1, 1.0, U, ldu, &work[ip], dim + 1, &work[iwrk], lwrk / (k + 1));
            dgemm_ovwr_left(0, n, k, dim, 1.0, V, ldv, &work[iq], dim, &work[iwrk], lwrk / k);

            rnorm = work[ibeta + k - 1];
        }
        iter++;
    }

    // Calculate singular vectors if requested
    if ((nconv >= *neig || *info > 0) && (jobu || jobv))
    {
        dcopy_(&dim, &work[ialpha], &int1, &work[ialpha1], &int1);
        dcopy_(&dim, &work[ibeta], &int1, &work[ibeta1], &int1);
        lwrk = lwrk + dim * dim + (dim + 1) * (dim + 1);
        dritzvec(which, jobu, jobv, m, n, nconv, dim, &work[ialpha1], &work[ibeta1], U, ldu, V, ldv, &work[ip], lwrk, iwork);
    }

    *neig = nconv;
}


void clansvd_irl(int which, int jobu, int jobv, int m, int n, int dim, int p, int *neig, int maxiter,
                 PROPACK_aprod_c aprod, PROPACK_CPLXF_TYPE* U, int ldu, float* sigma, float* bnd,
                 PROPACK_CPLXF_TYPE* V, int ldv, float tolin, float* work, int lwork,
                 PROPACK_CPLXF_TYPE* cwork, int lcwork, int* iwork, float* soption,
                 int* ioption, int* info, PROPACK_CPLXF_TYPE* cparm, int* iparm, uint64_t* rng_state)
{
    // Parameters
    int int1 = 1, int0 = 0;

    // Local variables
    int k, i, ibnd, iwrk, ierr, ip, iq, nconv, lwrk, kold;
    int ialpha, ibeta, ialpha1, ibeta1, ishift, nshft, lapinfo;
    float eps, eps34, epsn2, epsn, sfmin, anorm, rnorm, tol;
    float shift, relgap;
    int iter;

    // Set machine dependent constants
    eps = FLT_EPSILON;
    eps34 = powf(eps, 3.0f/4.0f);
    epsn = (float)fmaxf(m, n) * eps / 2.0f;
    epsn2 = sqrtf((float)fmaxf(m, n)) * eps / 2.0f;
    sfmin = FLT_MIN;

    // Guard against absurd arguments
    dim = fminf(fminf(dim, n + 1), m + 1);
    k = dim - p;
    tol = fminf(1.0f, fmaxf(16.0f * eps, tolin));
    anorm = 0.0f;

    // Set pointers into work array
    ibnd = 0;
    ialpha = ibnd + dim + 1;
    ibeta = ialpha + dim;
    ialpha1 = ibeta + dim;
    ibeta1 = ialpha1 + dim;
    ishift = ibeta1 + dim;
    ip = ishift + dim;
    iq = ip + (dim + 1) * (dim + 1);
    iwrk = iq + dim * dim;
    lwrk = lwork - iwrk + 1;

    // Zero out work array sections
    for (i = 0; i < 8 * dim + 3 + 2 * dim * dim; i++) { work[i] = 0.0f; }

    // Set up random starting vector if none is provided by the user
    rnorm = scnrm2_(&m, U, &int1);
    if (rnorm == 0.0f)
    {
        cgetu0(0, m, n, 0, 1, U, &rnorm, U, ldu, aprod, cparm, iparm, &ierr, ioption[0], &anorm, cwork, rng_state);
    }

    iter = 0;
    *info = 0;
    nconv = 0;
    kold = 0;

    // Iterate until convergence...
    while ((nconv < *neig) && (iter < maxiter))
    {
        // Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j}
        clanbpro(m, n, kold, &dim, aprod, U, ldu, V, ldv, &work[ialpha], dim, &rnorm,
                 &soption[0], &ioption[0], &work[iwrk], cwork, iwork, cparm, iparm, &ierr, rng_state);
        kold = k;

        // Compute and analyze SVD(B) and error bounds
        scopy_(&dim, &work[ialpha], &int1, &work[ialpha1], &int1);
        scopy_(&dim, &work[ibeta], &int1, &work[ibeta1], &int1);

        // Zero out bounds array
        for (i = 0; i < dim + 1; i++) { work[ibnd + i] = 0.0f; }

        // QR factorization of bidiagonal matrix
        sbdqr((dim == fminf(m, n)), 0, dim, &work[ialpha1], &work[ibeta1],
              &work[ibnd + dim - 1], &work[ibnd + dim], &work[ip], dim + 1);

        // SVD of bidiagonal matrix
        sbdsqr_("U", &dim, &int0, &int1, &int0, &work[ialpha1], &work[ibeta1], work, &int1,
                &work[ibnd], &int1, work, &int1, &work[iwrk], &lapinfo);

        // Update anorm estimate
        if (dim > 5)
        {
            anorm = work[ialpha1];
        } else {
            anorm = fmaxf(anorm, work[ialpha1]);
        }

        // Scale error bounds by rnorm
        for (i = 0; i < dim; i++)
        {
            work[ibnd + i] = fabsf(rnorm * work[ibnd + i]);
        }

        // Refine error bounds using the "Gap theorem"
        if (which == 0)  // which == 'S' (smallest)
        {
            srefinebounds(fminf(m, n), dim, &work[ialpha1], &work[ibnd], epsn * anorm, eps34);
        } else {  // which == 'L' (largest)
            srefinebounds(fminf(m, n), fminf(dim, *neig), &work[ialpha1], &work[ibnd], epsn * anorm, eps34);
        }

        // Determine the number of converged singular values
        if (which == 0)  // smallest singular values
        {
            i = dim - *neig;
            nconv = 0;
            while (i < dim)
            {
                if (work[ibnd + i] <= tol * work[ialpha1])
                {
                    nconv++;
                    sigma[nconv - 1] = work[ialpha1 + i];
                    bnd[nconv - 1] = work[ibnd + i];
                }
                i++;
            }
        } else {  // largest singular values
            i = 0;
            nconv = 0;
            while (i < fminf(dim, *neig))
            {
                if (work[ibnd + i] <= tol * work[ialpha1 + i])
                {
                    sigma[nconv] = work[ialpha1 + i];
                    bnd[nconv] = work[ibnd + i];
                    nconv++;
                    i++;
                } else {
                    i = k + 1;
                }
            }
        }

        // Test if an invariant subspace has been found
        if (ierr < 0)
        {
            if (dim < k)
            {
                *info = dim;
            }
            break;
        }

        if (nconv < *neig)
        {
            // Implicit restart:
            // Apply shifts mu_1, mu_2,...,mu_p to "back up" the
            // bidiagonalization (dim-k) steps to
            //
            //      A V_{k}^{+} = U_{k+1}^{+} B_{k}^{+}
            // corresponding to starting vector
            //      u_1^{+} = \prod_{i=1}^{p} (A A^T - mu_i^2) u_1
            //
            // We use exact shifts mu_i for which the relative gap between mu_i
            // and the lower bound on the k'th Ritzvalue is larger than soption[3]

            // Zero out shift array
            for (i = 0; i < dim - k; i++) { work[ishift + i] = 0.0f; }

            nshft = 0;
            if (which == 0)  // smallest
            {
                for (i = 0; i < k; i++)
                {
                    relgap = (work[ialpha1 + i] - work[ibnd + i] - work[ialpha1 + dim - *neig - 1]);
                    if (relgap > soption[3] * work[ialpha1 + dim - *neig - 1])
                    {
                        work[ishift + nshft] = work[ialpha1 + i];
                    } else {
                        work[ishift + nshft] = work[ialpha1];
                    }
                    nshft++;
                }
            } else {  // largest
                for (i = dim - 1; i >= k; i--)
                {
                    relgap = work[ialpha1 + k - 1] - (work[ialpha1 + i] + work[ibnd + i]);
                    if (relgap > soption[3] * work[ialpha1 + k - 1])
                    {
                        work[ishift + nshft] = work[ialpha1 + i];
                    } else {
                        work[ishift + nshft] = 0.0f;
                    }
                    nshft++;
                }
            }

            // Apply shifts and accumulate rotations such that
            //   B_{dim}^{+} = P * B_{dim} * Q^T

            // Zero out rotation matrices and set to identity
            for (i = 0; i < (dim + 1) * (dim + 1); i++) { work[ip + i] = 0.0f; }
            for (i = 0; i < dim * dim; i++) { work[iq + i] = 0.0f; }

            // Set diagonal elements to 1 (identity matrices)
            for (i = 0; i < dim + 1; i++)
            {
                work[ip + i * (dim + 2)] = 1.0f;
            }
            for (i = 0; i < dim; i++)
            {
                work[iq + i * (dim + 1)] = 1.0f;
            }

            // Apply shifts
            for (i = dim - 1; i >= k; i--)
            {
                shift = work[ishift + dim - 1 - i];
                sbsvdstep(1, 1, dim + 1, dim, i + 1, shift, &work[ialpha], &work[ibeta],
                         &work[ip], dim + 1, &work[iq], dim);
            }

            // Compute first k+1 left and first k right updated Lanczos vectors
            //   U_{dim+1}^{+} = U_{dim+1} * P(:,1:k+1)
            //   V_{dim}^{+} = V_{dim} * Q(:,1:k)

            csgemm_ovwr_left(0, m, k + 1, dim + 1, U, ldu, &work[ip], dim + 1, cwork, lcwork / (k+1));
            csgemm_ovwr_left(0, n, k, dim, V, ldv, &work[iq], dim, cwork, lcwork / k);

            rnorm = work[ibeta + k - 1];
        }
        iter++;
    }

    // Calculate singular vectors if requested
    if ((nconv >= *neig || *info > 0) && (jobu || jobv))
    {
        scopy_(&dim, &work[ialpha], &int1, &work[ialpha1], &int1);
        scopy_(&dim, &work[ibeta], &int1, &work[ibeta1], &int1);
        lwrk = lwrk + dim * dim + (dim + 1) * (dim + 1);
        critzvec(which, jobu, jobv, m, n, nconv, dim, &work[ialpha1], &work[ibeta1], U, ldu, V, ldv, &work[ip], lwrk, cwork, lcwork, iwork);
    }

    *neig = nconv;
}


void zlansvd_irl(int which, int jobu, int jobv, int m, int n, int dim, int p, int *neig, int maxiter,
                 PROPACK_aprod_z aprod, PROPACK_CPLX_TYPE* U, int ldu, double* sigma, double* bnd,
                 PROPACK_CPLX_TYPE* V, int ldv, double tolin, double* work, int lwork,
                 PROPACK_CPLX_TYPE* zwork, int lzwork, int* iwork, double* doption,
                 int* ioption, int* info, PROPACK_CPLX_TYPE* zparm, int* iparm, uint64_t* rng_state)
{
    // Parameters
    int int1 = 1, int0 = 0;

    // Local variables
    int k, i, ibnd, iwrk, ierr, ip, iq, nconv, lwrk, kold;
    int ialpha, ibeta, ialpha1, ibeta1, ishift, nshft, lapinfo;
    double eps, eps34, epsn2, epsn, sfmin, anorm, rnorm, tol;
    double shift, relgap;
    int iter;

    // Set machine dependent constants
    eps = DBL_EPSILON;
    eps34 = pow(eps, 3.0/4.0);
    epsn = (double)fmax(m, n) * eps / 2.0;
    epsn2 = sqrt((double)fmax(m, n)) * eps / 2.0;
    sfmin = DBL_MIN;

    // Guard against absurd arguments
    dim = fmin(fmin(dim, n + 1), m + 1);
    k = dim - p;
    tol = fmin(1.0, fmax(16.0 * eps, tolin));
    anorm = 0.0;

    // Set pointers into work array
    ibnd = 0;
    ialpha = ibnd + dim + 1;
    ibeta = ialpha + dim;
    ialpha1 = ibeta + dim;
    ibeta1 = ialpha1 + dim;
    ishift = ibeta1 + dim;
    ip = ishift + dim;
    iq = ip + (dim + 1) * (dim + 1);
    iwrk = iq + dim * dim;
    lwrk = lwork - iwrk + 1;

    // Zero out work array sections
    for (i = 0; i < 8 * dim + 3 + 2 * dim * dim; i++) { work[i] = 0.0; }

    // Set up random starting vector if none is provided by the user
    rnorm = dznrm2_(&m, U, &int1);
    if (rnorm == 0.0)
    {
        zgetu0(0, m, n, 0, 1, U, &rnorm, U, ldu, aprod, zparm, iparm, &ierr, ioption[0], &anorm, zwork, rng_state);
    }

    iter = 0;
    *info = 0;
    nconv = 0;
    kold = 0;

    // Iterate until convergence...
    while ((nconv < *neig) && (iter < maxiter))
    {
        // Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j}
        zlanbpro(m, n, kold, &dim, aprod, U, ldu, V, ldv, &work[ialpha], dim, &rnorm,
                 &doption[0], &ioption[0], &work[iwrk], zwork, iwork, zparm, iparm, &ierr, rng_state);
        kold = k;

        // Compute and analyze SVD(B) and error bounds
        dcopy_(&dim, &work[ialpha], &int1, &work[ialpha1], &int1);
        dcopy_(&dim, &work[ibeta], &int1, &work[ibeta1], &int1);

        // Zero out bounds array
        for (i = 0; i < dim + 1; i++) { work[ibnd + i] = 0.0; }

        // QR factorization of bidiagonal matrix
        dbdqr((dim == fmin(m, n)), 0, dim, &work[ialpha1], &work[ibeta1],
              &work[ibnd + dim - 1], &work[ibnd + dim], &work[ip], dim + 1);

        // SVD of bidiagonal matrix
        dbdsqr_("U", &dim, &int0, &int1, &int0, &work[ialpha1], &work[ibeta1], work, &int1,
                &work[ibnd], &int1, work, &int1, &work[iwrk], &lapinfo);

        // Update anorm estimate
        if (dim > 5)
        {
            anorm = work[ialpha1];
        } else {
            anorm = fmax(anorm, work[ialpha1]);
        }

        // Scale error bounds by rnorm
        for (i = 0; i < dim; i++)
        {
            work[ibnd + i] = fabs(rnorm * work[ibnd + i]);
        }

        // Refine error bounds using the "Gap theorem"
        if (which == 0)  // which == 'S' (smallest)
        {
            drefinebounds(fmin(m, n), dim, &work[ialpha1], &work[ibnd], epsn * anorm, eps34);
        } else {  // which == 'L' (largest)
            drefinebounds(fmin(m, n), fmin(dim, *neig), &work[ialpha1], &work[ibnd], epsn * anorm, eps34);
        }

        // Determine the number of converged singular values
        if (which == 0)  // smallest singular values
        {
            i = dim - *neig;
            nconv = 0;
            while (i < dim)
            {
                if (work[ibnd + i] <= tol * work[ialpha1])
                {
                    nconv++;
                    sigma[nconv - 1] = work[ialpha1 + i];
                    bnd[nconv - 1] = work[ibnd + i];
                }
                i++;
            }
        } else {  // largest singular values
            i = 0;
            nconv = 0;
            while (i < fmin(dim, *neig))
            {
                if (work[ibnd + i] <= tol * work[ialpha1 + i])
                {
                    sigma[nconv] = work[ialpha1 + i];
                    bnd[nconv] = work[ibnd + i];
                    nconv++;
                    i++;
                } else {
                    i = k + 1;
                }
            }
        }

        // Test if an invariant subspace has been found
        if (ierr < 0)
        {
            if (dim < k)
            {
                *info = dim;
            }
            break;
        }

        if (nconv < *neig)
        {
            // Implicit restart:
            // Apply shifts mu_1, mu_2,...,mu_p to "back up" the
            // bidiagonalization (dim-k) steps to
            //
            //      A V_{k}^{+} = U_{k+1}^{+} B_{k}^{+}
            // corresponding to starting vector
            //      u_1^{+} = \prod_{i=1}^{p} (A A^T - mu_i^2) u_1
            //
            // We use exact shifts mu_i for which the relative gap between mu_i
            // and the lower bound on the k'th Ritzvalue is larger than doption[3]

            // Zero out shift array
            for (i = 0; i < dim - k; i++) { work[ishift + i] = 0.0; }

            nshft = 0;
            if (which == 0)  // smallest
            {
                for (i = 0; i < k; i++)
                {
                    relgap = (work[ialpha1 + i] - work[ibnd + i] - work[ialpha1 + dim - *neig - 1]);
                    if (relgap > doption[3] * work[ialpha1 + dim - *neig - 1])
                    {
                        work[ishift + nshft] = work[ialpha1 + i];
                    } else {
                        work[ishift + nshft] = work[ialpha1];
                    }
                    nshft++;
                }
            } else {  // largest
                for (i = dim - 1; i >= k; i--)
                {
                    relgap = work[ialpha1 + k - 1] - (work[ialpha1 + i] + work[ibnd + i]);
                    if (relgap > doption[3] * work[ialpha1 + k - 1])
                    {
                        work[ishift + nshft] = work[ialpha1 + i];
                    } else {
                        work[ishift + nshft] = 0.0;
                    }
                    nshft++;
                }
            }

            // Apply shifts and accumulate rotations such that
            //   B_{dim}^{+} = P * B_{dim} * Q^T

            // Zero out rotation matrices and set to identity
            for (i = 0; i < (dim + 1) * (dim + 1); i++) { work[ip + i] = 0.0; }
            for (i = 0; i < dim * dim; i++) { work[iq + i] = 0.0; }

            // Set diagonal elements to 1 (identity matrices)
            for (i = 0; i < dim + 1; i++)
            {
                work[ip + i * (dim + 2)] = 1.0;
            }
            for (i = 0; i < dim; i++)
            {
                work[iq + i * (dim + 1)] = 1.0;
            }

            // Apply shifts
            for (i = dim - 1; i >= k; i--)
            {
                shift = work[ishift + dim - 1 - i];
                dbsvdstep(1, 1, dim + 1, dim, i + 1, shift, &work[ialpha], &work[ibeta],
                         &work[ip], dim + 1, &work[iq], dim);
            }

            // Compute first k+1 left and first k right updated Lanczos vectors
            //   U_{dim+1}^{+} = U_{dim+1} * P(:,1:k+1)
            //   V_{dim}^{+} = V_{dim} * Q(:,1:k)

            zdgemm_ovwr_left(0, m, k + 1, dim + 1, U, ldu, &work[ip], dim + 1, zwork, lzwork / (k+1));
            zdgemm_ovwr_left(0, n, k, dim, V, ldv, &work[iq], dim, zwork, lzwork / k);

            rnorm = work[ibeta + k - 1];
        }
        iter++;
    }

    // Calculate singular vectors if requested
    if ((nconv >= *neig || *info > 0) && (jobu || jobv))
    {
        dcopy_(&dim, &work[ialpha], &int1, &work[ialpha1], &int1);
        dcopy_(&dim, &work[ibeta], &int1, &work[ibeta1], &int1);
        lwrk = lwrk + dim * dim + (dim + 1) * (dim + 1);
        zritzvec(which, jobu, jobv, m, n, nconv, dim, &work[ialpha1], &work[ibeta1], U, ldu, V, ldv, &work[ip], lwrk, zwork, lzwork, iwork);
    }

    *neig = nconv;
}
