#include "__nnls.h"
#include <stdio.h>

/* Algorithm NNLS: NONNEGATIVE LEAST SQUARES
*
* Given an m by n matrix A, an m-vector B, and an n-vector X, compute an
* n-vector X which solves the least squares problem
*
*            a * x = b  subject to x >= 0
*
* This is a C translation of the original Fortran code, which was developed by
* Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
* 1973 JUN 15, and published in the book "SOLVING LEAST SQUARES PROBLEMS",
* Prentice-HalL, 1974. Revised FEB 1995 to accompany reprinting of the book
* (DOI: 10.1137/1.9781611971217) by SIAM.
*
*/
void
__nnls(const int m, const int n, double* restrict a, double* restrict b,
       double* restrict x, double* restrict w, double* restrict zz,
       int* restrict indices, const int maxiter, double* rnorm, int* info)
{
    int i = 0, ii = 0, ip = 0, indz = 0, iteration = 0, iz = 0, izmax = 0;
    int j = 0, jj = 0, k = 0, one = 1, tmpint = 0;
    double tau = 0.0, unorm = 0.0, ztest, alpha, cc, ss, wmax, T, tmp_work;
    double pivot = 1.0, pivot2 = 0.0, tmp = 0.0, spacing = 0.0;
    *info = 1;
    if (m <= 0 || n <= 0)
    {
        *info = 2;
        return;
    }

    // Initialize the indices and the solution vector x.
    for (i = 0; i < n; i++) { indices[i] = i; }
    for (i = 0; i < n; i++) { x[i] = 0.0; }

    // Outer loop
    while (indz < (m < n ? m : n))
    {
        // Compute the dual vector components in set Z.
        // Essentially a permuted gemv operation via BLAS ddot, in NumPy notation;
        // w[indices[indz:]] = A[indz:m, indices[indsz:]] @ b[indz:m]
        for (i = indz; i < n; i++)
        {
            j = indices[i];
            tmpint = m - indz;
            w[j] = ddot_(&tmpint, &a[indz + j*m], &one, &b[indz], &one);
        }

        // Find the next linearly independent column that corresponds to the
        // largest entry in the dual vector w.
        // ====================================================================
        while (1)
        {
            // Finding the largest w[j] and its index
            // izmax, wmax = argmax(w[indices[indz:]])
            wmax = 0.0;
            for (k = indz; k < n; k++)
            {
                j = indices[k];
                if (w[j] > wmax) { wmax = w[j]; izmax = k; }
            }
            // If wmax <= 0.0, terminate since this is a KKT certificate.
            if (wmax <= 0.0) { goto END; }
            iz = izmax;
            j = indices[iz];

            // The sign of wmax is OK for j to be moved to set p. Begin the
            // transformation and check new diagonal element to avoid near-linear
            // dependence.
            pivot = a[indz + j*m];
            tmpint = m - indz;
            dlarfgp_(&tmpint, &pivot, &a[indz + 1 + j*m], &one, &tau);

            // Compute the norm of a[0:indz, j] to check for linear dependence.
            unorm = (indz > 0 ? dnrm2_(&indz, &a[j*m], &one) : 0.0);
            // unorm is nonnegative
            spacing = (unorm > 0.0 ? nextafter(unorm, 2*unorm) - unorm : 0.0);

            // Test for independence by checking the pivot for zero.
            if (fabs(pivot) > 100.0*spacing)
            {
                // Column j is sufficiently independent. Copy b into zz and solve
                // for ztest which is the new prospective value for x[j].
                for (i = 0; i < m; i++) { zz[i] = b[i]; }
                tmpint = m - indz;
                pivot2 = a[indz + j*m];
                a[indz + j*m] = 1.0;
                dlarf_("L", &tmpint, &one, &a[indz + j*m], &one, &tau, &zz[indz], &tmpint, &tmp_work);
                // See if ztest is positive. This is from the original F77 code.
                // Probably better to use a sign test instead of a division.
                ztest = zz[indz] / pivot;
                if (ztest > 0.0)
                {
                    break;
                } else {
                    a[indz + j*m] = pivot2;
                }
            }
            // Reject j as a candidate to be moved from set z to set p.
            // a(indz,j) is restored, set w(j)=0., and loop back to test dual
            // coeffs again.
            w[j] = 0.0;
        }
        // ====================================================================

        // the index j=indices[iz]  has been selected to be moved from set z to
        // set p. Update b, update indices, apply householder transformations to
        // cols in new set z,  zero subdiagonal elements in col j,  set w(j)=0.
        for (i = 0; i < m; i++) { b[i] = zz[i]; }
        indices[iz] = indices[indz];
        indices[indz] = j;
        indz++;
        // Apply the householder transformation to the remaining columns.
        if (indz < n)
        {
            tmpint = m - indz + 1;
            for (k = indz; k < n; k++)
            {
                jj = indices[k];
                dlarf_("L", &tmpint, &one, &a[indz - 1 + j*m], &one, &tau, &a[indz - 1 + jj*m], &tmpint, &tmp_work);
            }
        }
        // Restore the pivot element into a, zero the subdiagonal elements in col j
        a[indz - 1 + j*m] = pivot;
        if (indz < m) { for (i = indz; i < m; i++) { a[j*m + i] = 0.0; } }
        // Zero the dual coefficient for the column.
        w[j] = 0.0;

        // Solve the permuted triangular system, store in zz.
        for (k = 0; k < indz; k++)
        {
            // ip traverses the indices of P set in reverse
            ip = indz - 1 - k;
            if (k != 0)
            {
                for (i = 0; i <= ip; i++)
                {
                    zz[i] = zz[i] - a[i + jj*m] * zz[ip + 1];
                }
            }
            jj = indices[ip];
            zz[ip] = zz[ip] / a[ip + jj*m];
        }

        // ****** Inner loop ******
        while (1)
        {
            iteration++;
            if (iteration >= maxiter) { *info = 3; goto END; }

            // See if all new constrained coefficients are feasible,
            // if not compute alpha
            alpha = 2.0;
            for (ip = 0; ip < indz; ip++)
            {
                k = indices[ip];
                if (zz[ip] <= 0.0)
                {
                    T = -x[k] / (zz[ip] - x[k]);
                    if (alpha > T)
                    {
                        alpha = T;
                        jj = ip;
                    }
                }
            }
            // If all new constrained coefficients are feasible, alpha is still
            // 2.0. If so exit from secondary loop to main loop.
            if (alpha == 2.0) { break; }  // Get back to outer loop

            // Otherwise interpolate between old x and zz.
            for (ip = 0; ip < indz; ip++)
            {
                k = indices[ip];
                x[k] = x[k] + alpha*(zz[ip] - x[k]);
            }

            // Modify a, b, and the indicies to move coefficient i from set p
            // to set z. While loop simulates a goto in the original F77 code.
            i = indices[jj];
            while (1)
            {
                x[i] = 0.0;

                if (jj != indz-1)
                {
                    jj++;
                    for (j = jj; j < indz; j++)
                    {
                        ii = indices[j];
                        indices[j-1] = ii;
                        dlartgp_(&a[j-1 + ii*m], &a[j + ii*m], &cc, &ss, &a[j-1 + ii*m]);
                        a[j + ii*m] = 0.0;
                        // Apply the Givens rotation to all columns except ii.
                        // Because the columns are not ordered we do it manually.
                        for (k = 0; k < n; k++)
                        {
                            if (k != ii)
                            {
                                tmp = a[j-1 + k*m];
                                a[j-1 + k*m] =  cc*tmp + ss*a[j + k*m];
                                a[j   + k*m] = -ss*tmp + cc*a[j + k*m];
                            }
                        }
                        tmp = b[j-1];
                        b[j-1] =  cc*tmp + ss*b[j];
                        b[j]   = -ss*tmp + cc*b[j];
                    }
                }
                indz--;
                indices[indz] = i;

                // See if remaining coefficients in set P are feasible
                // since determination of alpha guarantees it. If still
                // there are infeasible ones, they are due to numerical
                // noise. Any that are nonpositive will be set to zero
                // and moved from set p to set z.
                int nobreak = 0;
                for (jj = 0; jj < indz; jj++)
                {
                    i = indices[jj];
                    if (x[i] <= 0.0) { break; }
                    if (jj == indz - 1) { nobreak = 1; }
                }
                // If for loop completes without break, then leave the while loop
                if (nobreak) { break; }
            }

            for (i = 0; i < m; i++) { zz[i] = b[i]; }
            for (k = 0; k < indz; k++)
            {
                ip = indz - 1 - k;
                if (k != 0)
                {
                    for (i = 0; i <= ip; i++)
                    {
                        zz[i] = zz[i] - a[i + jj*m] * zz[ip + 1];

                    }
                }
                jj = indices[ip];
                zz[ip] = zz[ip] / a[ip + jj*m];
            }
            // ****** end of inner loop ******
        }

        // Back in the outer loop
        for (k = 0; k < indz; k++)
        {
            i = indices[k];
            x[i] = zz[k];
        }
        // ****** end of outer loop ******
    }
END:
    // Compute the residual vector and its norm.
    if (indz < m)
    {
        tmpint = m - indz;
        *rnorm = dnrm2_(&tmpint, &b[indz], &one);
    } else {
        for (i = 0; i < n; i++) { w[i] = 0.0; }
        *rnorm = 0.0;
    }
    return;
}
