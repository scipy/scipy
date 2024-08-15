# cython: boundscheck=False
# cython: initializedcheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: cpow=True


__all__ = ['_nnls']

from scipy.linalg.cython_lapack cimport dlarfgp, dlarf, dlartgp
from scipy.linalg.cython_blas cimport dnrm2
import numpy as np
cimport numpy as cnp
cnp.import_array()

def _nnls(cnp.ndarray[cnp.float64_t, ndim=2] A_in,
          cnp.ndarray[cnp.float64_t, ndim=1] b_in,
          int maxiter):
    # Make copies of the input to be mutated
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode='c'] A = A_in.copy(order='C')
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] b = b_in.copy()

    cdef int m = <int>A.shape[0], n = <int>A.shape[1]
    cdef int i = 0, ii = 0, ip = 0, iteration = 0, iz = 0, iz1 = 0, izmax = 0
    cdef int j = 0, jj = 0, k = 0
    cdef int col = 0, nrow = 0, nsetp = 0, one = 1, tmpint = 0
    cdef double tau = 0.0, unorm = 0.0, ztest, tmp, alpha, beta, cc, ss, wmax, T
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] w
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] x
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] work
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] zz
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode='c'] inds
    cdef bint skip = False

    inds = cnp.PyArray_Arange(0, n, 1, cnp.NPY_INT32)
    w = cnp.PyArray_EMPTY(1, [n], cnp.NPY_FLOAT64, 0)
    work = cnp.PyArray_EMPTY(1, [m], cnp.NPY_FLOAT64, 0)
    x = cnp.PyArray_ZEROS(1, [n], cnp.NPY_FLOAT64, 0)
    zz = cnp.PyArray_EMPTY(1, [m], cnp.NPY_FLOAT64, 0)

    # Quit if all coefficients are already in the solution or if m columns of A
    # have been triangularized.
    while (iz1 < n) and (nsetp < m):
        # simulating a goto from col independence check
        if skip:
            skip = False
        else:
            w[inds[iz1:]] = b[nrow:] @ A[nrow:, inds[iz1:]]

        # Find the largest w[j] and its index.
        wmax = 0.0
        for col in range(iz1, n):
            j = inds[col]
            if w[j] > wmax:
                wmax = w[j]
                izmax = col
        iz = izmax
        j = inds[iz]

        # If wmax <= 0.0, terminate since this is a KKT certificate.
        if wmax <= 0.0:
            break

        # The sign of wmax is OK for j to be moved to set p. Begin the transformation
        # and check new diagonal element to avoid near-linear dependence.
        work[nrow:] = A[nrow:, j]
        tmpint = m - nrow
        # DLARFGP( N, ALPHA, X, INCX, TAU )
        dlarfgp(&tmpint, &work[nrow], &work[nrow+1], &one, &tau)
        beta = work[nrow]
        work[nrow] = 1.
        unorm = 0.0
        if nsetp > 0:
            unorm = dnrm2(&nsetp, &A[0, j], &n)

        if ((unorm + abs(beta)*0.01) - unorm) > 0.0:
            # Column j is sufficiently independent. Copy b into zz and solve for
            # ztest which is the new prospective value for x[j].
            zz[:] = b[:]
            # dlarf(SIDE, M, N, V, INCV, TAU, C, LDC, WORK)
            dlarf(<char*>'L', &tmpint, &one, &work[nrow], &one, &tau,
                  &zz[nrow], &tmpint, &tmp)
            ztest = zz[nrow] / beta

            if ztest <= 0.0:
                # reject column j as a candidate to be moved from set z to set p.
                # Set w[j] to 0.0 and move to the next greatest entry in w.
                w[j] = 0.0
                skip = True
                continue
        else:
            # Column j is not numerically independent, reject column j
            w[j] = 0.0
            skip = True
            continue

        # column j accepted
        A[nrow, j] = beta
        b[:] = zz[:]
        inds[iz] = inds[iz1]
        inds[iz1] = j
        iz1 += 1
        nsetp += 1

        if iz1 < n:
            # Apply the householder trafo to remaining columns
            for col in inds[iz1:]:
                zz[nrow:] = A[nrow:, col]
                dlarf(<char*>'L', &tmpint, &one, &work[nrow], &one, &tau,
                      &zz[nrow], &tmpint, &tmp)
                A[nrow:, col] = zz[nrow:]

        nrow += 1

        if nsetp < m-1:
            A[nrow:, j] = 0.0

        w[j] = 0.0

        # Solve triangular system, store in zz
        zz[:] = b[:]
        for k in range(nsetp):
            ip = nsetp - k - 1
            if k != 0:
                for ii in range(ip+1):
                    zz[ii] -= A[ii, jj]*zz[ip+1]
            jj = inds[ip]
            zz[ip] /= A[ip, jj]

        # Inner loop
        while True:
            iteration += 1

            if iteration == maxiter:
                return x, 0.0, -1

            # See if all new constrained coefficients are feasible
            # otherwise compute alpha that should be in [0, 1]
            alpha = 2.0
            for ip in range(nsetp):
                k = inds[ip]
                if zz[ip] <= 0.0:
                    T = -x[k] / (zz[ip]-x[k])
                    if alpha > T:
                        alpha = T
                        jj = ip

            # If all new constrained coefficients are feasible
            # alpha is still 2 then exit otherwise interpolate
            # between old x and zz.
            if alpha == 2.0:
                break

            x[inds[:nsetp]] *= 1 - alpha
            x[inds[:nsetp]] += alpha*zz[:nsetp]

            # Modify A, B, and the indices to move coefficient
            # i from set p to set z. While loop simulates a goto
            i = inds[jj]
            while True:
                x[i] = 0.0

                if jj != nsetp:
                    jj += 1
                    for j in range(jj, nsetp):
                        ii = inds[j]
                        inds[j-1] = ii
                        dlartgp(&A[j-1, ii], &A[j, ii], &cc, &ss, &A[j-1, ii])
                        A[j, ii] = 0.0
                        # Apply Givens rotation to all cols except ii
                        for col in range(n):
                            if col != ii:
                                tmp = A[j-1, col]
                                A[j-1, col] = cc*tmp + ss*A[j, col]
                                A[j, col] = -ss*tmp + cc*A[j, col]

                        tmp = b[j-1]
                        b[j-1] = cc*tmp + ss*b[j]
                        b[j] = -ss*tmp + cc*b[j]

                nrow -= 1
                nsetp -= 1
                iz1 -= 1
                inds[iz1] = i

                # See if remaining coefficients in set P are feasible
                # since determination of alpha guarantees it. If still
                # there are infeasible ones, they are due to numerical
                # noise. Any that are nonpositive will be set to zero
                # and moved from set p to set z.
                for jj in range(nsetp):
                    i = inds[jj]
                    if x[i] <= 0.0:
                        # numerical noise; back to top of while loop
                        break
                else:
                    # No break; leave while loop
                    break

            zz[:] = b[:]
            for k in range(nsetp):
                ip = nsetp - k - 1
                if k != 0:
                    for ii in range(ip+1):
                        zz[ii] -= A[ii, jj]*zz[ip+1]
                jj = inds[ip]
                zz[ip] /= A[ip, jj]

            # Back to inner loop beginning

        # Back in outer loop
        x[inds[:nsetp]] = zz[:nsetp]

        # Back to the outer loop beginning

    return x, np.linalg.norm(b[nrow:]), 0
