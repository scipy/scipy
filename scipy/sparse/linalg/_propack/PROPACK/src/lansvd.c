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


void slansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_s aprod,
             float* U, int ldu, float* sigma, float* bnd, float* V, int ldv,
             float tolin, float* work, int lwork, int* iwork,
             float* doption, int* ioption, int* info, float* dparm, int* iparm,
             uint64_t* rng_state)
{
    // Parameters
    int int1 = 1, int0 = 0;

    // Local variables
    int i, j, dj, jold, ibnd, ib, ib1, iwrk, ierr, ip, iq, neig, lwrk, lapinfo, lanmax, nlandim;
    float eps, eps34, epsn2, epsn, sfmin, anorm, rnorm, tol;

    // Set machine dependent constants
    eps = FLT_EPSILON*0.5f;
    eps34 = powf(eps, 3.0f/4.0f);
    epsn = (float)int_max(m, n) * eps / 2.0f;
    epsn2 = sqrtf((float)int_max(m, n)) * eps / 2.0f;
    sfmin = FLT_MIN;

    // Guard against absurd arguments
    lanmax = int_min(int_min(n + 1, m + 1), kmax);
    tol = fminf(1.0f, fmaxf(16.0f * eps, tolin));
    anorm = 0.0f;

    // Set pointers into work array
    ibnd = 0;
    ib = ibnd + lanmax + 1;
    ib1 = ib + 2 * lanmax;
    ip = ib1 + 2 * lanmax;
    iq = ip + (lanmax + 1) * (lanmax + 1);
    iwrk = iq + lanmax * lanmax;
    lwrk = lwork - iwrk;

    // Zero out work array sections
    for (i = 0; i < 7 * lanmax + 2 + 2 * lanmax * lanmax; i++) { work[i] = 0.0f; }

    // Set up random starting vector if none is provided by the user
    rnorm = snrm2_(&m, &U[0], &int1);  // U(:,0) in 0-based indexing
    if (rnorm == 0.0f)
    {
        sgetu0(0, m, n, 0, 1, &U[0], &rnorm, U, ldu, aprod, dparm, iparm, &ierr, ioption[0], &anorm, &work[iwrk], rng_state);
    }

    *info = 0;
    neig = 0;
    jold = 0;
    j = int_min(k + int_max(k, 8) + 1, lanmax);

    // Iterate until convergence...
    while (neig < k)
    {
        // Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j}
        slanbpro(m, n, jold, &j, aprod, U, ldu, V, ldv, &work[ib], lanmax, &rnorm, doption, ioption, &work[iwrk], iwork, dparm, iparm, &ierr, rng_state);
        jold = j;

        // Compute and analyze SVD(B) and error bounds
        int two_lanmax = 2 * lanmax;
        scopy_(&two_lanmax, &work[ib], &int1, &work[ib1], &int1);

        // Zero out bounds array
        for (i = 0; i < j + 1; i++) { work[ibnd + i] = 0.0f; }

        // QR factorization of bidiagonal matrix
        sbdqr((j == int_min(m, n)), 0, j, &work[ib1], &work[ib1 + lanmax], &work[ibnd + j - 1], &work[ibnd + j], &work[ip], lanmax + 1);

        // SVD of bidiagonal matrix
        sbdsqr_("U", &j, &int0, &int1, &int0, &work[ib1], &work[ib1 + lanmax], work, &int1, &work[ibnd], &int1, work, &int1, &work[iwrk], &lapinfo);

        // Update anorm estimate
        if (j > 5)
        {
            anorm = work[ib1];
        } else {
            anorm = fmaxf(anorm, work[ib1]);
        }

        // Scale error bounds by rnorm
        for (i = 0; i < j+1; i++) { work[ibnd + i] = fabsf(rnorm * work[ibnd + i]); }

        // Refine error bounds using the "Gap theorem"
        srefinebounds(int_min(m, n), j, &work[ib1], &work[ibnd], epsn * anorm, eps34);

        // Determine the number of converged singular values
        for (i = 0; i < int_min(j, k); i++) { bnd[i] = work[ibnd + i]; }

        i = 0;
        neig = 0;
        while (i < int_min(j, k))
        {
            if (work[ibnd + i] <= tol * work[ib1 + i])
            {
                sigma[neig] = work[ib1 + i];
                neig++;
                i++;
            } else {
                break;
            }
        }

        // Test if an invariant subspace has been found or
        // if the workspace has been exhausted.
        if (ierr < 0)
        {
            // Invariant subspace found
            if (j < k) { *info = j; }
            break;
        }

        if (j >= lanmax)
        {
            // Maximum dimension of Krylov subspace exceeded
            if (neig < k) { *info = -1; }
            break;
        }

        // Increase the dimension of the Krylov subspace.
        // If any Ritz values have converged then try to estimate the average
        // number of iterations per converged Ritz value.
        // Else increase the dimension by 50%.
        if (neig > 1)
        {
            dj = int_min(j / 2, ((k - neig) * (j - 6)) / (2 * neig + 1));
            dj = int_min(100, int_max(2, dj));
        } else {
            dj = j / 2;
            dj = int_min(100, int_max(10, dj));
        }
        j = int_min(j + dj, lanmax);
    }

    // Calculate singular vectors if requested
    if (((neig >= k) || (*info > 0)) && (jobu || jobv))
    {
        lwrk = lwrk + lanmax * lanmax + (lanmax + 1) * (lanmax + 1);
        sritzvec(1, jobu, jobv, m, n, neig, jold, &work[ib], &work[ib + lanmax], U, ldu, V, ldv, &work[ip], lwrk, iwork);
    }

    k = neig;
    nlandim = j;
}


void dlansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_d aprod,
             double* U, int ldu, double* sigma, double* bnd, double* V, int ldv,
             double tolin, double* work, int lwork, int* iwork,
             double* doption, int* ioption, int* info, double* dparm, int* iparm,
             uint64_t* rng_state)
{
    // Parameters
    int int1 = 1, int0 = 0;

    // Local variables
    int i, j, dj, jold, ibnd, ib, ib1, iwrk, ierr, ip, iq, neig, lwrk, lapinfo, lanmax, nlandim;
    double eps, eps34, epsn2, epsn, sfmin, anorm, rnorm, tol;

    // Set machine dependent constants
    eps = DBL_EPSILON*0.5;
    eps34 = pow(eps, 3.0/4.0);
    epsn = (double)int_max(m, n) * eps / 2.0;
    epsn2 = sqrt((double)int_max(m, n)) * eps / 2.0;
    sfmin = DBL_MIN;

    // Guard against absurd arguments
    lanmax = int_min(int_min(n + 1, m + 1), kmax);
    tol = fmin(1.0, fmax(16.0 * eps, tolin));
    anorm = 0.0;

    // Set pointers into work array
    ibnd = 0;
    ib = ibnd + lanmax + 1;
    ib1 = ib + 2 * lanmax;
    ip = ib1 + 2 * lanmax;
    iq = ip + (lanmax + 1) * (lanmax + 1);
    iwrk = iq + lanmax * lanmax;
    lwrk = lwork - iwrk;

    // Zero out work array sections
    for (i = 0; i < 7 * lanmax + 2 + 2 * lanmax * lanmax; i++) { work[i] = 0.0; }

    // Set up random starting vector if none is provided by the user
    rnorm = dnrm2_(&m, &U[0], &int1);
    if (rnorm == 0.0)
    {
        dgetu0(0, m, n, 0, 1, &U[0], &rnorm, U, ldu, aprod, dparm, iparm, &ierr, ioption[0], &anorm, &work[iwrk], rng_state);
    }

    *info = 0;
    neig = 0;
    jold = 0;
    j = int_min(k + int_max(k, 8) + 1, lanmax);

    // Iterate until convergence...
    while (neig < k)
    {
        // Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j}
        dlanbpro(m, n, jold, &j, aprod, U, ldu, V, ldv, &work[ib], lanmax, &rnorm, doption, ioption, &work[iwrk], iwork, dparm, iparm, &ierr, rng_state);
        jold = j;

        // Compute and analyze SVD(B) and error bounds
        int two_lanmax = 2 * lanmax;
        dcopy_(&two_lanmax, &work[ib], &int1, &work[ib1], &int1);

        // Zero out bounds array
        for (i = 0; i < j + 1; i++) { work[ibnd + i] = 0.0; }

        // QR factorization of bidiagonal matrix
        dbdqr((j == int_min(m, n)), 0, j, &work[ib1], &work[ib1 + lanmax], &work[ibnd + j - 1], &work[ibnd + j], &work[ip], lanmax + 1);

        // SVD of bidiagonal matrix
        dbdsqr_("U", &j, &int0, &int1, &int0, &work[ib1], &work[ib1 + lanmax], work, &int1, &work[ibnd], &int1, work, &int1, &work[iwrk], &lapinfo);

        // Update anorm estimate
        if (j > 5)
        {
            anorm = work[ib1];
        } else {
            anorm = fmax(anorm, work[ib1]);
        }

        // Scale error bounds by rnorm
        for (i = 0; i < j+1; i++) { work[ibnd + i] = fabs(rnorm * work[ibnd + i]); }

        // Refine error bounds using the "Gap theorem"
        drefinebounds(int_min(m, n), j, &work[ib1], &work[ibnd], epsn * anorm, eps34);

        // Determine the number of converged singular values
        for (i = 0; i < int_min(j, k); i++) { bnd[i] = work[ibnd + i]; }

        i = 0;
        neig = 0;
        while (i < int_min(j, k))
        {
            if (work[ibnd + i] <= tol * work[ib1 + i])
            {
                sigma[neig] = work[ib1 + i];
                neig++;
                i++;
            } else {
                break;
            }
        }

        // Test if an invariant subspace has been found or
        // if the workspace has been exhausted.
        if (ierr < 0)
        {
            // Invariant subspace found
            if (j < k) { *info = j; }
            break;
        }

        if (j >= lanmax)
        {
            // Maximum dimension of Krylov subspace exceeded
            if (neig < k) { *info = -1; }
            break;
        }

        // Increase the dimension of the Krylov subspace.
        // If any Ritz values have converged then try to estimate the average
        // number of iterations per converged Ritz value.
        // Else increase the dimension by 50%.
        if (neig > 1)
        {
            dj = int_min(j / 2, ((k - neig) * (j - 6)) / (2 * neig + 1));
            dj = int_min(100, int_max(2, dj));
        } else {
            dj = j / 2;
            dj = int_min(100, int_max(10, dj));
        }
        j = int_min(j + dj, lanmax);
    }

    // Calculate singular vectors if requested
    if (((neig >= k) || (*info > 0)) && (jobu || jobv))
    {
        lwrk = lwrk + lanmax * lanmax + (lanmax + 1) * (lanmax + 1);
        dritzvec(1, jobu, jobv, m, n, neig, jold, &work[ib], &work[ib + lanmax], U, ldu, V, ldv, &work[ip], lwrk, iwork);
    }

    k = neig;
    nlandim = j;
}


void clansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_c aprod,
             PROPACK_CPLXF_TYPE* U, int ldu, float* sigma, float* bnd, PROPACK_CPLXF_TYPE* V, int ldv,
             float tolin, float* work, int lwork, PROPACK_CPLXF_TYPE* cwork, int lcwork,
             int* iwork, float* soption, int* ioption, int* info,
             PROPACK_CPLXF_TYPE* cparm, int* iparm, uint64_t* rng_state)
{
    // Parameters
    int int1 = 1, int0 = 0;

    // Local variables
    int i, j, dj, jold, ibnd, ib, ib1, iwrk, ierr, ip, iq, neig, lwrk, lapinfo, lanmax, nlandim;
    float eps, eps34, epsn2, epsn, sfmin, anorm, rnorm, tol;

    // Set machine dependent constants
    eps = FLT_EPSILON*0.5f;
    eps34 = powf(eps, 3.0f/4.0f);
    epsn = (float)int_max(m, n) * eps / 2.0f;
    epsn2 = sqrtf((float)int_max(m, n)) * eps / 2.0f;
    sfmin = FLT_MIN;

    // Guard against absurd arguments
    lanmax = int_min(int_min(n + 1, m + 1), kmax);
    tol = fminf(1.0f, fmaxf(16.0f * eps, tolin));
    anorm = 0.0f;

    // Set pointers into work array
    ibnd = 0;
    ib = ibnd + lanmax + 1;
    ib1 = ib + 2 * lanmax;
    ip = ib1 + 2 * lanmax;
    iq = ip + (lanmax + 1) * (lanmax + 1);
    iwrk = iq + lanmax * lanmax;
    lwrk = lwork - iwrk;

    // Zero out work array sections
    for (i = 0; i < 7 * lanmax + 2 + 2 * lanmax * lanmax; i++) { work[i] = 0.0f; }

    // Zero out complex work array
    for (i = 0; i < lcwork; i++) { cwork[i] = PROPACK_cplxf(0.0f, 0.0f); }

    // Set up random starting vector if none is provided by the user
    rnorm = scnrm2_(&m, &U[0], &int1);
    if (rnorm == 0.0f)
    {
        cgetu0(0, m, n, 0, 1, &U[0], &rnorm, U, ldu, aprod, cparm, iparm, &ierr, ioption[0], &anorm, cwork, rng_state);
    }

    *info = 0;
    neig = 0;
    jold = 0;
    j = int_min(k + int_max(k, 8) + 1, lanmax);

    // Iterate until convergence...
    while (neig < k)
    {
        // Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j}
        clanbpro(m, n, jold, &j, aprod, U, ldu, V, ldv, &work[ib], lanmax, &rnorm, soption, ioption, &work[iwrk], cwork, iwork, cparm, iparm, &ierr, rng_state);
        jold = j;

        // Compute and analyze SVD(B) and error bounds
        int two_lanmax = 2 * lanmax;
        scopy_(&two_lanmax, &work[ib], &int1, &work[ib1], &int1);

        // Zero out bounds array
        for (i = 0; i < j + 1; i++) { work[ibnd + i] = 0.0f; }

        // QR factorization of bidiagonal matrix
        sbdqr((j == int_min(m, n)), 0, j, &work[ib1], &work[ib1 + lanmax], &work[ibnd + j - 1], &work[ibnd + j], &work[ip], lanmax + 1);

        // SVD of bidiagonal matrix
        sbdsqr_("U", &j, &int0, &int1, &int0, &work[ib1], &work[ib1 + lanmax], work, &int1, &work[ibnd], &int1, work, &int1, &work[iwrk], &lapinfo);

        // Update anorm estimate
        if (j > 5)
        {
            anorm = work[ib1];
        } else {
            anorm = fmaxf(anorm, work[ib1]);
        }

        // Scale error bounds by rnorm
        for (i = 0; i < j+1; i++) { work[ibnd + i] = fabsf(rnorm * work[ibnd + i]); }

        // Refine error bounds using the "Gap theorem"
        srefinebounds(int_min(m, n), j, &work[ib1], &work[ibnd], epsn * anorm, eps34);

        // Determine the number of converged singular values
        for (i = 0; i < int_min(j, k); i++) { bnd[i] = work[ibnd + i]; }

        i = 0;
        neig = 0;
        while (i < int_min(j, k))
        {
            if (work[ibnd + i] <= tol * work[ib1 + i])
            {
                sigma[neig] = work[ib1 + i];
                neig++;
                i++;
            } else {
                break;
            }
        }

        // Test if an invariant subspace has been found or
        // if the workspace has been exhausted.
        if (ierr < 0)
        {
            // Invariant subspace found
            if (j < k) { *info = j; }
            break;
        }

        if (j >= lanmax)
        {
            // Maximum dimension of Krylov subspace exceeded
            if (neig < k) { *info = -1; }
            break;
        }

        // Increase the dimension of the Krylov subspace.
        // If any Ritz values have converged then try to estimate the average
        // number of iterations per converged Ritz value.
        // Else increase the dimension by 50%.
        if (neig > 1)
        {
            dj = int_min(j / 2, ((k - neig) * (j - 6)) / (2 * neig + 1));
            dj = int_min(100, int_max(2, dj));
        } else {
            dj = j / 2;
            dj = int_min(100, int_max(10, dj));
        }
        j = int_min(j + dj, lanmax);
    }

    // Calculate singular vectors if requested
    if (((neig >= k) || (*info > 0)) && (jobu || jobv))
    {
        lwrk = lwrk + lanmax * lanmax + (lanmax + 1) * (lanmax + 1);
        critzvec(1, jobu, jobv, m, n, neig, jold, &work[ib], &work[ib + lanmax], U, ldu, V, ldv, &work[ip], lwrk, cwork, lcwork, iwork);
    }

    k = neig;
    nlandim = j;
}


void zlansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_z aprod,
             PROPACK_CPLX_TYPE* U, int ldu, double* sigma, double* bnd, PROPACK_CPLX_TYPE* V, int ldv,
             double tolin, double* work, int lwork, PROPACK_CPLX_TYPE* zwork, int lzwork,
             int* iwork, double* doption, int* ioption, int* info,
             PROPACK_CPLX_TYPE* zparm, int* iparm, uint64_t* rng_state)
{
    // Parameters
    int int1 = 1, int0 = 0;

    // Local variables
    int i, j, dj, jold, ibnd, ib, ib1, iwrk, ierr, ip, iq, neig, lwrk, lapinfo, lanmax, nlandim;
    double eps, eps34, epsn2, epsn, sfmin, anorm, rnorm, tol;

    // Set machine dependent constants
    eps = DBL_EPSILON*0.5;
    eps34 = pow(eps, 3.0/4.0);
    epsn = (double)int_max(m, n) * eps / 2.0;
    epsn2 = sqrt((double)int_max(m, n)) * eps / 2.0;
    sfmin = DBL_MIN;

    // Guard against absurd arguments
    lanmax = int_min(int_min(n + 1, m + 1), kmax);
    tol = fmin(1.0, fmax(16.0 * eps, tolin));
    anorm = 0.0;

    // Set pointers into work array
    ibnd = 0;
    ib = ibnd + lanmax + 1;
    ib1 = ib + 2 * lanmax;
    ip = ib1 + 2 * lanmax;
    iq = ip + (lanmax + 1) * (lanmax + 1);
    iwrk = iq + lanmax * lanmax;
    lwrk = lwork - iwrk;

    // Zero out work array sections
    for (i = 0; i < 7 * lanmax + 2 + 2 * lanmax * lanmax; i++) { work[i] = 0.0; }

    // Zero out complex work array
    for (i = 0; i < lzwork; i++) { zwork[i] = PROPACK_cplx(0.0, 0.0); }

    // Set up random starting vector if none is provided by the user
    rnorm = dznrm2_(&m, &U[0], &int1);
    if (rnorm == 0.0)
    {
        zgetu0(0, m, n, 0, 1, &U[0], &rnorm, U, ldu, aprod, zparm, iparm, &ierr, ioption[0], &anorm, zwork, rng_state);
    }

    *info = 0;
    neig = 0;
    jold = 0;
    j = int_min(k + int_max(k, 8) + 1, lanmax);

    // Iterate until convergence...
    while (neig < k)
    {
        // Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j}
        zlanbpro(m, n, jold, &j, aprod, U, ldu, V, ldv, &work[ib], lanmax, &rnorm, doption, ioption, &work[iwrk], zwork, iwork, zparm, iparm, &ierr, rng_state);
        jold = j;

        // Compute and analyze SVD(B) and error bounds
        int two_lanmax = 2 * lanmax;
        dcopy_(&two_lanmax, &work[ib], &int1, &work[ib1], &int1);

        // Zero out bounds array
        for (i = 0; i < j + 1; i++) { work[ibnd + i] = 0.0; }

        // QR factorization of bidiagonal matrix
        dbdqr((j == int_min(m, n)), 0, j, &work[ib1], &work[ib1 + lanmax], &work[ibnd + j - 1], &work[ibnd + j], &work[ip], lanmax + 1);

        // SVD of bidiagonal matrix
        dbdsqr_("U", &j, &int0, &int1, &int0, &work[ib1], &work[ib1 + lanmax], work, &int1, &work[ibnd], &int1, work, &int1, &work[iwrk], &lapinfo);

        // Update anorm estimate
        if (j > 5)
        {
            anorm = work[ib1];
        } else {
            anorm = fmax(anorm, work[ib1]);
        }

        // Scale error bounds by rnorm
        for (i = 0; i < j+1; i++) { work[ibnd + i] = fabs(rnorm * work[ibnd + i]); }

        // Refine error bounds using the "Gap theorem"
        drefinebounds(int_min(m, n), j, &work[ib1], &work[ibnd], epsn * anorm, eps34);

        // Determine the number of converged singular values
        for (i = 0; i < int_min(j, k); i++) { bnd[i] = work[ibnd + i]; }

        i = 0;
        neig = 0;
        while (i < int_min(j, k))
        {
            if (work[ibnd + i] <= tol * work[ib1 + i])
            {
                sigma[neig] = work[ib1 + i];
                neig++;
                i++;
            } else {
                break;
            }
        }

        // Test if an invariant subspace has been found or
        // if the workspace has been exhausted.
        if (ierr < 0)
        {
            // Invariant subspace found
            if (j < k) { *info = j; }
            break;
        }

        if (j >= lanmax)
        {
            // Maximum dimension of Krylov subspace exceeded
            if (neig < k) { *info = -1; }
            break;
        }

        // Increase the dimension of the Krylov subspace.
        // If any Ritz values have converged then try to estimate the average
        // number of iterations per converged Ritz value.
        // Else increase the dimension by 50%.
        if (neig > 1)
        {
            dj = int_min(j / 2, ((k - neig) * (j - 6)) / (2 * neig + 1));
            dj = int_min(100, int_max(2, dj));
        } else {
            dj = j / 2;
            dj = int_min(100, int_max(10, dj));
        }
        j = int_min(j + dj, lanmax);
    }

    // Calculate singular vectors if requested
    if (((neig >= k) || (*info > 0)) && (jobu || jobv))
    {
        lwrk = lwrk + lanmax * lanmax + (lanmax + 1) * (lanmax + 1);
        zritzvec(1, jobu, jobv, m, n, neig, jold, &work[ib], &work[ib + lanmax], U, ldu, V, ldv, &work[ip], lwrk, zwork, lzwork, iwork);
    }

    k = neig;
    nlandim = j;
}
