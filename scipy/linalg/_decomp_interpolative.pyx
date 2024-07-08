# cython: boundscheck=False
# cython: initializedcheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: cpow=True

"""
This file is a Cython rewrite of the original Fortran code of "ID: A software package
for low-rank approximation of matrices via interpolative decompositions, Version 0.4",
written by Per-Gunnar Martinsson, Vladimir Rokhlin, Yoel Shkolnisky, and Mark Tygert.

The original Fortran code can be found at the last author's current website
http://tygert.com/software.html


References
----------

N. Halko, P.G. Martinsson, and J. A. Tropp, "Finding structure with randomness:
probabilistic algorithms for constructing approximate matrix decompositions",
SIAM Review, 53 (2011), pp. 217-288. DOI:10.1137/090771806

H. Cheng, Z. Gimbutas, P.G. Martinsson, V.Rokhlin, "On the Compression of Low
Rank Matrices", SIAM Journal of Scientific Computing, 2005, Vol.26(4),
DOI:10.1137/030602678



Copyright (C) 2024 SciPy developers

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

a. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
b. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
c. Names of the SciPy Developers may not be used to endorse or promote
   products derived from this software without specific prior written
   permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.


Notes
-----

The translated functions from the original Fortran77 code are as follows (with various
internal functions subsumed into respective functions):

    idd_diffsnorm
    idd_estrank
    idd_findrank
    idd_id2svd
    idd_ldiv
    idd_poweroftwo
    idd_reconid
    idd_snorm
    iddp_aid
    iddp_asvd
    iddp_id
    iddp_qrpiv
    iddp_rid
    iddp_rsvd
    iddp_svd
    iddr_aid
    iddr_asvd
    iddr_id
    iddr_qrpiv
    iddr_rid
    iddr_rsvd
    iddr_svd
    idz_diffsnorm
    idz_estrank
    idz_findrank
    idz_id2svd
    idz_reconid
    idz_snorm
    idzp_aid
    idzp_asvd
    idzp_id
    idzp_qrpiv
    idzp_rid
    idzp_rsvd
    idzp_svd
    idzr_aid
    idzr_asvd
    idzr_id
    idzr_rid
    idzr_rsvd
    idzr_qrpiv
    idzr_svd

"""

import numpy as np
from numpy.typing import NDArray
cimport numpy as cnp
cnp.import_array()

from cpython.mem cimport PyMem_Free, PyMem_Malloc, PyMem_Realloc
from libc.math cimport hypot

import scipy.linalg as la
from scipy.fft import rfft, fft
from scipy.sparse.linalg import LinearOperator

from scipy.linalg.cython_lapack cimport dlarfgp, dorm2r, zunm2r, zlarfgp
from scipy.linalg.cython_blas cimport dnrm2, dtrsm, dznrm2, ztrsm


__all__ = ['idd_estrank', 'idd_ldiv', 'idd_poweroftwo', 'idd_reconid', 'iddp_aid',
           'iddp_asvd', 'iddp_id', 'iddp_qrpiv', 'iddp_svd', 'iddr_aid', 'iddr_asvd',
           'iddr_id', 'iddr_qrpiv', 'iddr_svd', 'idz_estrank', 'idz_reconid',
           'idzp_aid', 'idzp_asvd', 'idzp_id', 'idzp_qrpiv', 'idzp_svd', 'idzr_aid',
           'idzr_asvd', 'idzr_id', 'idzr_qrpiv', 'idzr_svd', 'idd_id2svd', 'idz_id2svd'
           # LinearOperator funcs
           'idd_findrank', 'iddp_rid', 'iddp_rsvd', 'iddr_rid', 'iddr_rsvd',
           'idz_findrank', 'idzp_rid', 'idzp_rsvd', 'idzr_rid', 'idzr_rsvd',
           'idd_snorm', 'idz_snorm', 'idd_diffsnorm', 'idz_diffsnorm'
           ]


def idd_diffsnorm(A: LinearOperator, B: LinearOperator, int its=20, rng=None):
    cdef int n = A.shape[1], j = 0, intone = 1
    cdef cnp.float64_t snorm = 0.0
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] v1
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] v2
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] u1
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] u2

    if not rng:
        rng = np.random.default_rng()
    v1 = rng.uniform(low=-1., high=1., size=n)
    v1 /= dnrm2(&n, &v1[0], &intone)

    for j in range(its):
        u1 = A.matvec(v1)
        u2 = B.matvec(v1)
        u1 -= u2
        v1 = A.rmatvec(u1)
        v2 = B.rmatvec(u1)
        v1 -= v2

        snorm = dnrm2(&n, &v1[0], &intone)
        if snorm > 0.0:
            v1 /= snorm

        snorm = np.sqrt(snorm)

    return snorm


def idd_estrank(cnp.ndarray[cnp.float64_t, mode="c", ndim=2] a: NDArray, eps: float,
                rng=None):
    cdef int m = a.shape[0], n = a.shape[1]
    cdef int intone = 1, n2, nsteps = 3, row, r, nstep, cols, k, nulls
    cdef cnp.float64_t h, alpha, beta
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=3] albetas
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] tau_arr
    cdef cnp.ndarray[cnp.int64_t, mode='c', ndim=1] subselect
    cdef cnp.float64_t *aa
    cdef cnp.float64_t *ff
    cdef cnp.float64_t[:, ::1] Fmemview
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] giv2x2
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] rta
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] Fc
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] F

    if not rng:
        rng = np.random.default_rng()

    n2 = idd_poweroftwo(m)

    # This part is the initialization that is done via idd_frmi
    # for a Subsampled Randomized Fourier Transfmrom (SRFT).

    # Draw (nsteps x m x 2) arrays from [-1, 1) uniformly and scale
    # each 2-element row to unity norm
    albetas = rng.uniform(low=-1.0, high=1.0, size=[nsteps, m, 2])
    aa = <cnp.float64_t *>cnp.PyArray_DATA(albetas)
    # Walk over every 2D row and normalize
    for r in range(0, 2*nsteps*m, 2):
        h = 1/hypot(aa[r], aa[r+1])
        aa[r] *= h
        aa[r+1] *= h

    # idd_random_transf
    rta = a.copy()

    # Rotate and shuffle "a" nsteps-many times
    giv2x2 = cnp.PyArray_ZEROS(2, [2, 2], cnp.NPY_FLOAT64, 0)
    for nstep in range(nsteps):
        for row in range(m-1):
            alpha, beta = albetas[nstep, row, 0], albetas[nstep, row, 1]
            giv2x2[0, 0] = alpha
            giv2x2[0, 1] = beta
            giv2x2[1, 0] = -beta
            giv2x2[1, 1] = alpha
            np.matmul(giv2x2, rta[row:row+2, :], out=rta[row:row+2, :])

        rta = rta[rng.permutation(m), :]

    # idd_subselect pick randomly n2-many rows
    subselect = rng.choice(m, n2, replace=False)
    rta = rta[subselect, :]

    # Perform rfft on each column. Note that the first and the last
    # element of the result is real valued (n2 is power of 2).
    #
    # We view the complex valued entries as two consecutive doubles
    # (by also removing the 2nd and last all-0 rows -- see idd_frm).
    # Then after transpose we do a final row shuffle after transpose.
    Fc = rfft(rta.T, axis=1)
    # Move the first col to second col
    Fc[:, 0] *= 1.j
    # Perform the final permutation
    F = Fc.view(np.float64)[:, 1:-1].T[rng.permutation(n2), :]

    Fcopy = F.copy()
    cols = F.shape[1]
    row = F.shape[0]
    sssmax = 0.
    ff = <cnp.float64_t *>cnp.PyArray_DATA(F)
    for r in range(cols):
        h = dnrm2(&row, &ff[r], &cols)
        if h > sssmax:
            sssmax = h

    tau_arr = cnp.PyArray_ZEROS(1, [cols], cnp.NPY_FLOAT64, 0)
    k, nulls = 0, 0

    # In Fortran id_dist, F is transposed and works on the columns
    # Since we have a C-array we work directly on rows
    # The reflectors are overwritten on rows of F directly
    # Hence at any k'th step, we have
    #
    #            [ B  r  r  r  r  r  r  r ]
    #            [           ....         ]
    #            [           ....         ]
    #            [ x  x  x  B  r  r  r  r ]
    #            [ x  x  x  x  B  r  r  r ]
    #            [ x  x  x  x  x  B  r  r ]
    #            [ x  x  x  x  x  x  x  x ]
    #            [ x  x  x  x  x  x  x  x ]
    #

    # Loop until nulls = 7, or krank+nulls = n2, or krank+nulls = n.
    Fmemview = F
    while (nulls < 7) and (k+nulls < min(n, n2)):
        # Apply previous Householder reflectors
        if k > 0:
            for kk in range(k):
                F[k, kk:] -= tau_arr[kk]*(F[kk, kk:] @ F[k, kk:])*F[kk, kk:]

        # Get the next Householder reflector and store in F
        r = cols-k
        # n, alpha, x, incx, tau
        dlarfgp(&r, &Fmemview[k, k], &Fmemview[k, k+1], &intone, &tau_arr[k])
        beta = F[k, k]
        F[k, k] = 1

        if (beta <= eps*sssmax):
            nulls += 1
        k += 1

    if nulls < 7:
        k = 0

    return k, Fcopy


def idd_findrank(A: LinearOperator, cnp.float64_t eps, rng=None):
    # Estimate the rank of A by repeatedly using A.rmatvec(random vec)

    cdef int m = A.shape[0], n = A.shape[1], k = 0, kk = 0,r = n, krank
    cdef int no_of_cols = 4, intone = 1, info = 0
    cdef cnp.float64_t[::1] tau = cnp.PyArray_ZEROS(1, [min(m, n)], cnp.NPY_FLOAT64, 0)
    cdef cnp.float64_t[::1] y = cnp.PyArray_ZEROS(1, [n], cnp.NPY_FLOAT64, 0)
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] retarr

    # The size of the QR decomposition is rank dependent which is unknown
    # at runtime. Hence we don't want to allocate a dense version of the
    # linear operator which can be too big. Instead, a typical "realloc double
    # if run out of space" strategy is used here. Starts with 4*n
    # Also, we hold the A.T @ x results in a separate array to return
    # and do the same for that too.
    cdef cnp.float64_t *ra = <cnp.float64_t*>PyMem_Malloc(
        sizeof(cnp.float64_t)*no_of_cols*n
        )
    cdef cnp.float64_t *reallocated_ra
    cdef cnp.float64_t *ret = <cnp.float64_t*>PyMem_Malloc(
        sizeof(cnp.float64_t)*no_of_cols*n
        )
    cdef cnp.float64_t *reallocated_ret
    cdef cnp.float64_t enorm = 0.0

    if (not ra) or (not ret):
        raise MemoryError("Failed to allocate at least required memory "
                          f"{no_of_cols*n*8} bytes for"
                          "'scipy.linalg.interpolative.idd_findrank()' "
                          "function.")

    if not rng:
        rng = np.random.default_rng()

    krank = 0
    try:
        while True:

            # Generate random vector and rmatvec then save the result
            x = rng.uniform(size=m)
            y = A.rmatvec(x)
            for kk in range(n):
                ret[krank*n + kk] = y[kk]

            if krank == 0:
                enorm = dnrm2(&n, &y[0], &intone)
            else:  # krank > 0
                # Transpose-Apply previous Householder reflectors, if any
                # SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO
                dorm2r(<char*>'L',<char*>'T', &n, &intone, &krank, &ra[0], &n,
                       &tau[0], &y[0], &n, &ra[(no_of_cols-1)*n], &info)

            # Get the next Householder reflector
            r = n-krank
            # N, ALPHA, X, INCX, TAU
            dlarfgp(&r, &y[krank], &y[krank+1], &intone, &tau[krank])

            for kk in range(n):
                ra[krank*n + kk] = y[kk]

            # Running out of space; try to double the size of ra
            if krank == (no_of_cols-2):
                reallocated_ra = <cnp.float64_t *>PyMem_Realloc(
                    ra, sizeof(cnp.float64_t)*no_of_cols*n*2)
                reallocated_ret = <cnp.float64_t *>PyMem_Realloc(
                    ret, sizeof(cnp.float64_t)*no_of_cols*n*2)

                if reallocated_ra and reallocated_ret:
                    ra = reallocated_ra
                    ret = reallocated_ret
                    no_of_cols *= 2
                else:
                    raise MemoryError(
                        "'scipy.linalg.interpolative.idd_findrank()' failed to "
                        f"allocate the required memory,{no_of_cols*n*16} bytes "
                        "while trying to determine the rank (currently "
                        f"{krank}) of a LinearOperator with precision {eps}."
                    )
            krank += 1
            if (y[krank-1] < eps*enorm) or (krank >= min(m, n)):
                break
    finally:
        # Crashed or successfully ended up here
        # Discard Householder vectors
        PyMem_Free(ra)
        retarr = cnp.PyArray_EMPTY(2, [krank, n], cnp.NPY_FLOAT64, 0)
        for k in range(krank):
            for kk in range(n):
                retarr[k, kk] = ret[k*n+kk]
        PyMem_Free(ret)

    return krank, retarr


def idd_id2svd(
    cnp.ndarray[cnp.float64_t, mode='c', ndim=2] cols,
    cnp.ndarray[cnp.int64_t, mode='c', ndim=1] perms,
    cnp.ndarray[cnp.float64_t, ndim=2] proj,
    ):
    cdef int m = cols.shape[0], krank = cols.shape[1]
    cdef int n = proj.shape[1] + krank, info, ci
    cdef cnp.ndarray[cnp.float64_t, mode='fortran', ndim=2] C
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] tau1
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] tau2
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] UU
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] S
    cdef cnp.ndarray[cnp.float64_t, ndim=2] V
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] VV
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds1
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds2
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] p

    UU = cnp.PyArray_ZEROS(2, [m, krank], cnp.NPY_FLOAT64, 0)
    VV = cnp.PyArray_ZEROS(2, [n, krank], cnp.NPY_FLOAT64, 0)
    p = cnp.PyArray_ZEROS(2, [krank, n], cnp.NPY_FLOAT64, 0)

    # idd_reconint
    for ci in range(krank):
        p[ci, perms[ci]] = 1.0

    p[:, perms[krank:]] = proj[:, :]

    inds1, tau1 = iddr_qrpiv(cols, krank)
    # idd_rinqr and idd_rearr
    r = np.triu(cols[:krank, :])
    for ci in range(krank-1, -1, -1):
        r[:, [ci, inds1[ci]]] = r[:,  [inds1[ci], ci]]

    t = p.T.copy()
    inds2, tau2 = iddr_qrpiv(t, krank)
    r2 = np.triu(t[:krank, :])
    for ci in range(krank-1, -1, -1):
        r2[:, [ci, inds2[ci]]] = r2[:,  [inds2[ci], ci]]

    r3 = r @ r2.T
    UU[:krank, :krank], S, V = la.svd(r3,
                                      full_matrices=False,
                                      check_finite=False)

    # Apply Q of col to U from the left, use cols as scratch
    C = cols[:, :krank].copy(order='F')
    dorm2r(<char*>'R', <char*>'T',
           &krank, &m, &krank, &C[0, 0], &m, &tau1[0],
           &UU[0,0], &krank, &cols[0, 0], &info)

    VV[:krank, :krank] = V[:, :].T
    # Apply Q of t to V from the left
    C = t[:, :krank].copy(order='F')
    dorm2r(<char*>'R', <char*>'T',
           &krank, &n, &krank, &C[0, 0], &n, &tau2[0],
           &VV[0, 0], &krank, &cols[0, 0], &info)

    return UU, S, VV


cdef inline int idd_ldiv(int l, int n) noexcept nogil:
    cdef int m = l
    while (n % m != 0):
        m -= 1
    return m


cdef int idd_poweroftwo(int m) noexcept nogil:
    """
    Find the integer solution to l = floor(log2(m))
    """
    cdef int n = 1
    while (n < m):
        n <<= 1  # Times 2
    return n >> 1  # Divide by 2


def idd_reconid(B, idx, proj):
    cdef int m = B.shape[0], krank = B.shape[1]
    cdef int n = len(idx)
    approx = np.zeros([m, n], dtype=np.float64)

    approx[:, idx[:krank]] = B
    approx[:, idx[krank:]] = B @ proj

    return approx


def idd_snorm(A: LinearOperator, int its=20, rng=None):
    cdef int n = A.shape[1]
    cdef int j = 0, intone = 1
    cdef cnp.float64_t snorm = 0.0
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] v
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] u

    if not rng:
        rng = np.random.default_rng()
    v = rng.uniform(low=-1., high=1., size=n)
    v /= dnrm2(&n, &v[0], &intone)

    for j in range(its):
        u = A.matvec(v)
        v = A.rmatvec(u)
        snorm = dnrm2(&n, &v[0], &intone)
        if snorm > 0.0:
            v /= snorm

        snorm = np.sqrt(snorm)

    return snorm


def iddp_aid(cnp.ndarray[cnp.float64_t, ndim=2] a: NDArray, eps: float, rng=None):
    krank, proj = idd_estrank(a, eps, rng=rng)
    if krank != 0:
        proj = proj[:krank, :]
        return iddp_id(proj, eps=eps)

    return iddp_id(a, eps=eps)


def iddp_asvd(cnp.ndarray[cnp.float64_t, ndim=2] a: NDArray, eps: float, rng=None):
    cdef int m = a.shape[0], n = a.shape[1]
    cdef int krank, info, ci
    cdef cnp.ndarray[cnp.float64_t, mode='fortran', ndim=2] C
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] tau1
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] tau2
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] UU
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] S
    cdef cnp.ndarray[cnp.float64_t, ndim=2] V
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] VV
    cdef cnp.ndarray[cnp.float64_t, ndim=2] proj
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] perms
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds1
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds2
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] p
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] col

    krank, perms, proj = iddp_aid(a.copy(), eps, rng=rng)

    if krank > 0:
        UU = cnp.PyArray_ZEROS(2, [m, krank], cnp.NPY_FLOAT64, 0)
        VV = cnp.PyArray_ZEROS(2, [n, krank], cnp.NPY_FLOAT64, 0)

        p = cnp.PyArray_ZEROS(2, [krank, n], cnp.NPY_FLOAT64, 0)
        col = a[:, perms[:krank]].copy()

        # idd_reconint
        for ci in range(krank):
            p[ci, perms[ci]] = 1.0

        # p[np.arange(krank), perms[:krank]] = 1.
        p[:, perms[krank:]] = proj[:, :]

        inds1, tau1 = iddr_qrpiv(col, krank)
        # idd_rinqr and idd_rearr
        r = np.triu(col[:krank, :])
        for ci in range(krank-1, -1, -1):
            r[:, [ci, inds1[ci]]] = r[:,  [inds1[ci], ci]]

        t = p.T.copy()
        inds2, tau2 = iddr_qrpiv(t, krank)
        r2 = np.triu(t[:krank, :])
        for ci in range(krank-1, -1, -1):
            r2[:, [ci, inds2[ci]]] = r2[:,  [inds2[ci], ci]]

        r3 = r @ r2.T
        UU[:krank, :krank], S, V = la.svd(r3, full_matrices=False)

        # Apply Q of col to U from the left
        C = col[:, :krank].copy(order='F')
        dorm2r(<char*>'R', <char*>'T',
               &krank, &m, &krank, &C[0, 0], &m, &tau1[0],
               &UU[0,0], &krank, &a[0, 0], &info)

        VV[:krank, :krank] = V[:, :].T
        # Apply Q of t to V from the left
        C = t[:, :krank].copy(order='F')
        dorm2r(<char*>'R', <char*>'T',
               &krank, &n, &krank, &C[0, 0], &n, &tau2[0],
               &VV[0, 0], &krank, &a[0, 0], &info)

    return UU, S, VV


def iddp_id(cnp.ndarray[cnp.float64_t, ndim=2] a: NDArray, eps: float):
    cdef int n = a.shape[1], krank, tmp_int, p
    cdef cnp.float64_t one = 1
    krank, _, inds = iddp_qrpiv(a, eps)

    # Change pivots to permutation
    perms = cnp.PyArray_ZEROS(1, [n], cnp.NPY_INT64, 0)
    for p in range(n):
        perms[p] = p

    if krank > 0:
        for p in range(krank):
            # Apply pivots
            tmp_int = perms[p]
            perms[p] = perms[inds[p]]
            perms[inds[p]] = tmp_int
            # perms[[p, inds[p]]] = perms[[inds[p], p]]

    # Let A = [A1, A2] and A1 has krank cols and upper triangular.
    # Find X that satisfies A1 @ X = A2
    # In SciPy.linalg this amounts to;
    #
    # proj = la.solve_triangular(a[:krank, :krank], a[:krank, krank:],
    #                            lower=False, check_finite=False)
    #
    # Push into BLAS without transposes.
    # A1 = a[:krank, :krank]
    # A2 = a[:krank, krank:]
    # Instead solve X @ A1.T = A2.T
    # Fortran already sees A1 as A1.T and becomes lower tri, side = R

    tmp_int = n - krank
    # SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB
    dtrsm(<char*>'R', <char*>'L', <char*>'N', <char*>'N',
          &tmp_int, &krank, &one, &a[0, 0], &n, &a[0, krank], &n)

    return krank, np.array(perms), a[:krank, krank:]


def iddp_qrpiv(cnp.ndarray[cnp.float64_t, mode="c", ndim=2] a, cnp.float64_t eps):
    """
    This is a minimal version of ?GEQP3 from LAPACK with an
    additional early stopping criterion over given precision.

    This function overwrites entries of "a" !
    """

    cdef int m = a.shape[0], n = a.shape[1]
    cdef cnp.ndarray col_norms = cnp.PyArray_ZEROS(1, [n], cnp.NPY_FLOAT64, 0)
    cdef int k = 0, kpiv = 0, i = 0, tmp_int = 0, int_n = 0
    cdef cnp.float64_t tmp_sca = 0.
    cdef cnp.ndarray taus = cnp.PyArray_ZEROS(1, [m], cnp.NPY_FLOAT64, 0)
    cdef cnp.ndarray ind = cnp.PyArray_ZEROS(1, [n], cnp.NPY_INT64, 0)
    cdef cnp.float64_t[::1] taus_v = taus
    cdef cnp.float64_t feps = 0.1e-16  # np.finfo(np.float64).eps
    cdef cnp.float64_t ssmax, ssmaxin
    cdef int nupdate = 0

    for i in range(n):
        col_norms[i] = dnrm2(&m, &a[0, i], &n)**2

    kpiv = np.argmax(col_norms)
    ssmax = col_norms[kpiv]
    ssmaxin = ssmax

    for k in range(min(m, n)):

        # Pivoting
        ind[k] = kpiv
        # Swap columns a[:, k] and a[:, kpiv]
        a[:, [kpiv, k]] = a[:, [k, kpiv]]

        # Swap col_norms[krank] and col_norms[kpiv]
        col_norms[[kpiv, k]] = col_norms[[k, kpiv]]

        if k < m-1:
            # Compute the householder reflector for column k
            tmp_sca = a[k, k]
            # FIX: Convert these to F_INT
            tmp_int = <int>(m - k)
            int_n = <int>n
            dlarfgp(&tmp_int, &tmp_sca, &a[k+1, k], &int_n, &taus_v[k])

            # Overwrite with 1. for easy matmul
            a[k, k] = 1
            if k < n-1:
                # Apply the householder reflector to the rest on the right
                a[k:, k+1:] -= np.outer(taus[k]*a[k:, k], a[k:, k] @ a[k:, k+1:])

            # Put back the beta in place
            a[k, k] = tmp_sca

            # Update the norms
            col_norms[k] = 0
            col_norms[k+1:] -= a[k, k+1:]**2
            ssmax = 0
            kpiv = k+1
            if k < n-1:
                kpiv = np.argmax(col_norms[k+1:]) + (k + 1)
                ssmax = col_norms[kpiv]

            if (((ssmax < 1000*feps*ssmaxin) and (nupdate == 0)) or
                    ((ssmax < ((1000*feps)**2)*ssmaxin) and (nupdate == 1))):
                nupdate += 1
                ssmax = 0
                kpiv = k+1

                if k < n-1:
                    for i in range(k+1, n):
                        tmp_int = m-k-1
                        col_norms[i] = dnrm2(&tmp_int, &a[k+1, i], &n)**2
                    kpiv = np.argmax(col_norms[k+1:]) + (k + 1)
                    ssmax = col_norms[kpiv]
        if (ssmax <= (eps**2)*ssmaxin):
            break
    # a is overwritten; return numerical rank and pivots
    return k + 1, taus, ind


def iddp_rid(A: LinearOperator, cnp.float64_t eps, rng=None):
    _, ret = idd_findrank(A, eps, rng)
    return iddp_id(ret, eps)


def iddp_rsvd(A: LinearOperator, cnp.float64_t eps, rng=None):
    cdef int n = A.shape[1]
    cdef int krank, j
    cdef cnp.ndarray[cnp.int64_t, mode='c', ndim=1] perms
    cdef cnp.ndarray[cnp.float64_t, ndim=2] proj
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] col
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] x

    krank, perms, proj = iddp_rid(A, eps, rng)
    if krank > 0:
        # idd_getcols
        col = cnp.PyArray_EMPTY(2, [n, krank], cnp.NPY_FLOAT64, 0)
        x = cnp.PyArray_ZEROS(1, [n], cnp.NPY_FLOAT64, 0)

        for j in range(krank):
            x[perms[j]] = 1.
            col[:, j] = A.matvec(x)
            x[perms[j]] = 0.

        return idd_id2svd(cols=col, perms=perms, proj=proj)

    # TODO: figure out empty return
    return None


def iddp_svd(cnp.ndarray[cnp.float64_t, ndim=2] a: NDArray, eps: float):
    """a is overwritten"""
    cdef int m = a.shape[0], krank, info
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] taus
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] UU
    cdef cnp.ndarray[cnp.float64_t, mode='fortran', ndim=2] C

    # Get the pivoted QR
    krank, taus, inds = iddp_qrpiv(a, eps)

    if krank > 0:
        r = np.triu(a[:krank, :])
        # Apply pivots in reverse
        for p in range(krank-1, -1, -1):
            r[:, [p, inds[p]]] = r[:, [inds[p], p]]

        # JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO
        # dgesvd(<char*>'S', <char*>'O', &krank, &n)
        U, S, V = la.svd(r, full_matrices=False)

        # Apply Q to U via dorm2r
        # Possibly U is shorter than Q
        UU = np.zeros([m, krank], dtype=a.dtype)
        UU[:krank, :krank] = U
        # Do the transpose dance for C-layout, use a for scratch
        C = a[:, :krank].copy(order='F')
        dorm2r(<char*>'R', <char*>'T',
               &krank, &m, &krank, &C[0, 0], &m, &taus[0],
               &UU[0,0], &krank, &a[0, 0], &info)

    return UU, S, V


def iddr_aid(cnp.ndarray[cnp.float64_t, mode="c", ndim=2] a: NDArray, int krank,
             rng=None):
    cdef int m = a.shape[0], n = a.shape[1], n2, nsteps = 3, row, r, nstep, L
    cdef cnp.float64_t h, alpha, beta
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=3] albetas
    cdef cnp.ndarray[cnp.npy_int64, mode='c', ndim=1] subselect
    cdef cnp.float64_t *aa
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] giv2x2
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] rta
    cdef cnp.ndarray[cnp.npy_int64, mode='c', ndim=1] marker

    if not rng:
        rng = np.random.default_rng()

    # idd_aidi
    L = krank + 8
    n2 = 0
    if (L >= n2) or (L > m):
        inds, proj = iddr_id(a, krank)
        return inds, proj

    n2 = idd_poweroftwo(m)

    # idd_sfrmi
    # idd_pairsamps
    ind = rng.permutation(n2)
    ind2 = cnp.PyArray_ZEROS(1, [L], cnp.NPY_INT64, 0)

    marker = cnp.PyArray_ZEROS(1, [n2//2], cnp.NPY_INT64, 0)
    for k in range(L):
        marker[(ind[k]+1)//2] = marker[(ind[k]+1)//2]+1

    for r in range(n2//2):
        if marker[r] != 0:
            l2 += 1
            ind2[r] = r

    # Draw (nsteps x m x 2) arrays from [-1, 1) uniformly and scale
    # each 2-element row to unity norm
    albetas = rng.uniform(low=-1.0, high=1.0, size=[nsteps, m, 2])
    aa = <cnp.float64_t *>cnp.PyArray_DATA(albetas)
    # Walk over every 2D row and normalize
    for r in range(0, 2*nsteps*m, 2):
        # ignoring the improbable zero generation by rng.uniform
        h = 1.0/hypot(aa[r], aa[r+1])
        aa[r] *= h
        aa[r+1] *= h

    # idd_random_transf
    rta = a.copy()

    # Rotate and shuffle "a" nsteps-many times
    giv2x2 = cnp.PyArray_ZEROS(2, [2, 2], cnp.NPY_FLOAT64, 0)
    for nstep in range(nsteps):
        for row in range(m-1):
            alpha, beta = albetas[nstep, row, 0], albetas[nstep, row, 1]
            giv2x2[0, 0] = alpha
            giv2x2[0, 1] = beta
            giv2x2[1, 0] = -beta
            giv2x2[1, 1] = alpha
            np.matmul(giv2x2, rta[row:row+2, :], out=rta[row:row+2, :])

        rta = rta[rng.permutation(m), :]

    # idd_subselect pick randomly n2-many rows
    subselect = rng.choice(m, n2, replace=False)
    rta = rta[subselect, :]

    # idd_sffti
    twopi = 2*np.pi
    twopii = twopi*1.j
    nblock = idd_ldiv(l2, n2)
    fact = 1/np.sqrt(n2)

    if l2 == 1:
        wsave = np.exp(-twopii*k*ind2[0]/np.arange(1, n2+1))*fact
    else:
        m = n2//nblock

        wsave = np.empty(m*l2, dtype=complex)
        for j in range(l2):
            i = ind2[j]
            if (i+1) <= (n//2 - m//2):
                idivm = i // m
                imodm = i - m*idivm
                for k in range(m):
                    wsave[m*j+k] = (
                        np.exp(-twopii*(k)*imodm/m)*
                        np.exp(-twopii*(k)*(idivm+1)/n)*
                        fact
                        )
            else:
                idivm = (i+1)//(m//2)
                imodm = (i+1)-(m//2)*idivm
                for k in range(m):
                    wsave[m*j+k] = np.exp(-twopii*(k-1)*imodm/m)*fact

    # idd_sfft.f
    # There is some significant index olympics happening in the original Fortran code
    # however I could not reverse engineer it to understand what is happening and kept
    # as is with all its cryptic movements and their performance hits.
    # See DOI:10.1016/j.acha.2007.12.002 - Section 3.3

    # Perform partial FFT to each nblock
    F = rfft(rta.reshape(nblock, m, -1), order='F', axis=0)
    # Roll the first entry to the last in the first axis for
    # the real frequency components. (faster than np.roll)
    F = F[[x for x in range(1, F.shape[0])] + [0], :, :]
    # Convert back to 2D array
    F = F.reshape(F.shape[0]*F.shape[1], -1)

    csum = np.zeros_like(F[0, :])
    rsum = np.zeros_like(F[0, :])

    for j in range(l2):
        i = ind2[j]
        if (i+1) <= (n//2 - m//2):
            idivm = i // m
            imodm = i - m*idivm
            csum[:] = 0.0
            for k in range(m):
                csum += F[m*idivm+k, :] * wsave[m*j+k]
            rta[2*i, :] = csum.real
            rta[2*i+1, :] = csum.imag

        else:
            idivm = (i+1)//(m//2)
            imodm = (i+1)-(m//2)*idivm
            csum[:] = 0.0
            for k in range(m):
                csum += F[m*(nblock//2)+k, :] * wsave[m*j+k]
            rta[2*i, :] = csum.real
            rta[2*i+1, :] = csum.imag
            if i == (n//2) - 1:
                for k in range(m):
                    rsum += F[m*(nblock//2)+k, :]
                rta[n-2, :] = rsum
                rta[n-2, :] *= fact

                rsum[:] = 0.0
                for k in range(m//2):
                    rsum += F[m*(nblock//2)+2*k-1]
                    rsum -= F[m*(nblock//2)+2*k]
                rta[n-1, :] = rsum
                rta[n-1, :] *= fact

    # idd_subselect pick randomly l2-many rows
    subselect = rng.choice(n2, l2, replace=False)
    rta = rta[subselect, :]

    perms, proj = iddr_id(rta, krank)

    return perms, proj


def iddr_asvd(cnp.ndarray[cnp.float64_t, mode="c", ndim=2] a: NDArray, int krank,
              rng=None):
    cdef int m = a.shape[0], n = a.shape[1]
    cdef int info, ci
    cdef cnp.ndarray[cnp.float64_t, mode='fortran', ndim=2] C
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] tau1
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] tau2
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] UU
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] S
    cdef cnp.ndarray[cnp.float64_t, ndim=2] V
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] VV
    cdef cnp.ndarray[cnp.float64_t, ndim=2] proj
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] perms
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds1
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds2
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] p
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] col

    perms, proj = iddr_aid(a.copy(), krank=krank, rng=rng)

    UU = cnp.PyArray_ZEROS(2, [m, krank], cnp.NPY_FLOAT64, 0)
    VV = cnp.PyArray_ZEROS(2, [n, krank], cnp.NPY_FLOAT64, 0)

    p = cnp.PyArray_ZEROS(2, [krank, n], cnp.NPY_FLOAT64, 0)
    col = a[:, perms[:krank]].copy()

    # idd_reconint
    for ci in range(krank):
        p[ci, perms[ci]] = 1.0

    p[:, perms[krank:]] = proj[:, :]

    inds1, tau1 = iddr_qrpiv(col, krank)
    # idd_rinqr and idd_rearr
    r = np.triu(col[:krank, :])
    for ci in range(krank-1, -1, -1):
        r[:, [ci, inds1[ci]]] = r[:,  [inds1[ci], ci]]

    t = p.T.copy()
    inds2, tau2 = iddr_qrpiv(t, krank)
    r2 = np.triu(t[:krank, :])
    for ci in range(krank-1, -1, -1):
        r2[:, [ci, inds2[ci]]] = r2[:,  [inds2[ci], ci]]

    r3 = r @ r2.T
    UU[:krank, :krank], S, V = la.svd(r3, full_matrices=False)

    # Apply Q of col to U from the left
    C = col[:, :krank].copy(order='F')
    dorm2r(<char*>'R', <char*>'T',
           &krank, &m, &krank, &C[0, 0], &m, &tau1[0],
           &UU[0,0], &krank, &a[0, 0], &info)

    VV[:krank, :krank] = V[:, :].T
    # Apply Q of t to V from the left
    C = t[:, :krank].copy(order='F')
    dorm2r(<char*>'R', <char*>'T',
           &krank, &n, &krank, &C[0, 0], &n, &tau2[0],
           &VV[0, 0], &krank, &a[0, 0], &info)

    return UU, S, VV


def iddr_id(cnp.ndarray[cnp.float64_t, ndim=2] a, int krank):
    cdef int n = a.shape[1]
    cdef int tmp_int
    cdef cnp.float64_t one = 1.0
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] perms

    inds, _ = iddr_qrpiv(a, krank)
    perms = cnp.PyArray_Arange(0, n, 1, cnp.NPY_INT64)

    if krank > 0:
        for p in range(krank):
            # Apply pivots
            tmp_int = perms[p]
            perms[p] = perms[inds[p]]
            perms[inds[p]] = tmp_int

    # See iddp_id comments for below
    tmp_int = n - krank
    # SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB
    dtrsm(<char*>'R', <char*>'L', <char*>'N', <char*>'N',
          &tmp_int, &krank, &one, &a[0, 0], &n, &a[0, krank], &n)

    return perms, a[:krank, krank:]


def iddr_qrpiv(cnp.ndarray[cnp.float64_t, mode="c", ndim=2] a: NDArray, krank: int):
    cdef int m = a.shape[0], n = a.shape[1]
    cdef cnp.ndarray col_norms = cnp.PyArray_ZEROS(1, [n], cnp.NPY_FLOAT64, 0)
    cdef int loop = 0, loops, kpiv = 0, i = 0, tmp_int = 0, int_n = 0
    cdef cnp.float64_t tmp_sca = 0.
    cdef cnp.ndarray taus = cnp.PyArray_ZEROS(1, [m], cnp.NPY_FLOAT64, 0)
    cdef cnp.ndarray ind = cnp.PyArray_ZEROS(1, [n], cnp.NPY_INT64, 0)
    cdef cnp.float64_t[::1] taus_v = taus
    cdef cnp.float64_t feps = 0.1e-16  # np.finfo(np.float64).eps
    cdef cnp.float64_t ssmax, ssmaxin
    cdef int nupdate = 0

    loops = min(krank, min(m, n))
    for i in range(n):
        col_norms[i] = dnrm2(&m, &a[0, i], &n)**2

    kpiv = np.argmax(col_norms)
    ssmax = col_norms[kpiv]
    ssmaxin = ssmax

    for loop in range(loops):

        ind[loop] = kpiv
        # Swap columns a[:, k] and a[:, kpiv]
        a[:, [kpiv, loop]] = a[:, [loop, kpiv]]
        # Swap col_norms[krank] and col_norms[kpiv]
        col_norms[[kpiv, loop]] = col_norms[[loop, kpiv]]

        if loop < m-1:
            tmp_sca = a[loop, loop]
            # FIX: Convert these to F_INT
            tmp_int = <int>(m - loop)
            int_n = <int>n
            dlarfgp(&tmp_int, &tmp_sca, &a[loop+1, loop], &int_n, &taus_v[loop])

            # Overwrite with 1. for easy matmul
            a[loop, loop] = 1
            if loop < n-1:
                # Apply the householder reflector to the rest on the right
                a[loop:, loop+1:] -= np.outer(taus[loop]*a[loop:, loop],
                                              a[loop:, loop] @ a[loop:, loop+1:])

            # Put back the beta in place
            a[loop, loop] = tmp_sca

            # Update the norms
            col_norms[loop] = 0
            col_norms[loop+1:] -= a[loop, loop+1:]**2
            ssmax = 0
            kpiv = loop+1

            if loop < n-1:
                kpiv = np.argmax(col_norms[loop+1:]) + (loop + 1)
                ssmax = col_norms[kpiv]
            if (((ssmax < 1000*feps*ssmaxin) and (nupdate == 0)) or
                    ((ssmax < ((1000*feps)**2)*ssmaxin) and (nupdate == 1))):
                nupdate += 1
                ssmax = 0
                kpiv = loop+1

                if loop < n-1:
                    for i in range(loop+1, n):
                        tmp_int = m-loop-1
                        col_norms[i] = dnrm2(&tmp_int, &a[loop+1, i], &n)**2
                    kpiv = np.argmax(col_norms[loop+1:]) + (loop + 1)
                    ssmax = col_norms[kpiv]

    return ind, taus


def iddr_rid(A: LinearOperator, int krank, rng=None):
    cdef int m = A.shape[0], n = A.shape[1], k = 0
    cdef int L = min(krank+2, min(m, n))
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] r

    if not rng:
        rng = np.random.default_rng()

    r = cnp.PyArray_EMPTY(2, [L, n], cnp.NPY_FLOAT64, 0)
    for k in range(L):
        r[k, :] = A.rmatvec(rng.uniform(size=m))

    return iddr_id(a=r, krank=krank)


def iddr_rsvd(A: LinearOperator, int krank, rng=None):
    cdef int n = A.shape[1], j
    cdef cnp.ndarray[cnp.int64_t, mode='c', ndim=1] perms
    cdef cnp.ndarray[cnp.float64_t, ndim=2] proj
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] col

    perms, proj = iddr_rid(A, krank, rng)
    # idd_getcols
    col = cnp.PyArray_EMPTY(2, [n, krank], cnp.NPY_FLOAT64, 0)
    x = cnp.PyArray_ZEROS(1, [n], cnp.NPY_FLOAT64, 0)
    for j in range(krank):
        x[perms[j]] = 1.
        col[:, j] = A.matvec(x)
        x[perms[j]] = 0.

    return idd_id2svd(cols=col, perms=perms, proj=proj)


def iddr_svd(cnp.ndarray[cnp.float64_t, mode="c", ndim=2] a: NDArray, int krank):
    cdef int m = a.shape[0], info = 0
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] taus
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] UU
    cdef cnp.ndarray[cnp.float64_t, mode='fortran', ndim=2] C

    # Get the pivoted QR
    inds, taus = iddr_qrpiv(a, krank)

    r = np.triu(a[:krank, :])
    # Apply pivots in reverse
    for p in range(krank-1, -1, -1):
        r[:, [p, inds[p]]] = r[:, [inds[p], p]]

    # JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO
    # dgesvd(<char*>'S', <char*>'O', &krank, &n)
    U, S, V = la.svd(r, full_matrices=False)

    # Apply Q to U via dorm2r
    # Possibly U is shorter than Q
    UU = np.zeros([m, krank], dtype=a.dtype)
    UU[:krank, :krank] = U
    # Do the transpose dance for C-layout, use a for scratch
    C = a[:, :krank].copy(order='F')
    dorm2r(<char*>'R', <char*>'T',
           &krank, &m, &krank, &C[0, 0], &m, &taus[0],
           &UU[0,0], &krank, &a[0, 0], &info)

    return UU, S, V


def idz_diffsnorm(A: LinearOperator, B: LinearOperator, int its=20, rng=None):
    cdef int n = A.shape[1], j = 0, intone = 1
    cdef cnp.float64_t snorm = 0.0
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] v1
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] v2
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] u1
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] u2

    if not rng:
        rng = np.random.default_rng()
    v1 = rng.uniform(low=-1, high=1, size=(n, 2)).view(np.complex128).ravel()
    v1 /= dznrm2(&n, &v1[0], &intone)

    for j in range(its):
        u1 = A.matvec(v1)
        u2 = B.matvec(v1)
        u1 -= u2
        v1 = A.rmatvec(u1)
        v2 = B.rmatvec(u1)
        v1 -= v2

        snorm = dznrm2(&n, &v1[0], &intone)
        if snorm > 0.0:
            v1 /= snorm

        snorm = np.sqrt(snorm)

    return snorm


def idz_estrank(cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] a: NDArray, eps: float,
                rng=None):
    cdef int m = a.shape[0], n = a.shape[1], n2, nsteps = 3, row, r, nstep, cols, k
    cdef cnp.float64_t h, alpha, beta
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=3] albetas
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] tau_arr
    cdef cnp.ndarray[cnp.npy_int64, mode='c', ndim=1] subselect
    cdef double complex[:, ::1] ff
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=2] giv2x2
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] rta
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] F

    if not rng:
        rng = np.random.default_rng()

    n2 = idd_poweroftwo(m)
    # This part is the initialization that is done via idz_frmi
    # for a Subsampled Randomized Fourier Transfmrom (SRFT).

    # Draw (nsteps x m x 4) array from [0, 2)*pi uniformly for
    # random points on complex unit circle and unitary rotations
    albetas = np.empty([nsteps, m, 4])
    albetas[:, :, 2:] = rng.uniform(low=0.0, high=2.0, size=[nsteps, m, 2])
    albetas[:, :, 2:] *= np.pi
    np.cos(albetas[:, :, 2], out=albetas[:, :, 0])
    np.sin(albetas[:, :, 2], out=albetas[:, :, 1])
    np.cos(albetas[:, :, 3], out=albetas[:, :, 2])
    np.sin(albetas[:, :, 3], out=albetas[:, :, 3])

    # idd_random_transf
    rta = a.copy()

    # Rotate and shuffle "a" nsteps-many times
    giv2x2 = cnp.PyArray_ZEROS(2, [2, 2], cnp.NPY_FLOAT64, 0)
    for nstep in range(nsteps):
        # Multiply with a point on the unit circle
        rta *= albetas[nstep, :, 2:].view(np.complex128)
        # Rotate
        for row in range(m-1):
            alpha, beta = albetas[nstep, row, 0], albetas[nstep, row, 1]
            giv2x2[0, 0] = alpha
            giv2x2[0, 1] = beta
            giv2x2[1, 0] = -beta
            giv2x2[1, 1] = alpha
            np.matmul(giv2x2, rta[row:row+2, :], out=rta[row:row+2, :])

        rta = rta[rng.permutation(m), :]

    # idd_subselect pick randomly n2-many rows
    subselect = rng.choice(m, n2, replace=False)
    rta = rta[subselect, :]
    # Perform rfft on each column.
    F = fft(rta, axis=0)[rng.permutation(n2), :]

    Fcopy = F.copy()
    cols = F.shape[1]
    row = F.shape[0]
    sssmax = 0.

    for r in range(cols):
        h = dznrm2(&row, &F[0, r], &cols)
        if h > sssmax:
            sssmax = h

    tau_arr = cnp.PyArray_ZEROS(1, [cols], cnp.NPY_COMPLEX128, 0)
    k, nulls = 0, 0
    ff = F
    # Loop until nulls = 7, or krank+nulls = n2, or krank+nulls = n.
    while (nulls < 7) and (k+nulls < min(n, n2)):
        # Apply previous Householder reflectors
        if k > 0:
            for kk in range(k):
                F[k, kk:] -= (
                    np.conj(tau_arr[kk])*
                    (F[kk, kk:].conj() @ F[k, kk:])*
                    F[kk, kk:]
                    )

        # Get the next Householder reflector and store in F
        r = cols-k
        row = 1
        zlarfgp(&r, &ff[k, k], &ff[k, k+1], &row, &tau_arr[k])
        if (np.abs(F[k, k]) <= eps*sssmax):
            nulls += 1
        F[k, k] = 1
        k += 1

    if nulls < 7:
        k = 0

    return k, Fcopy


def idz_findrank(A: LinearOperator, cnp.float64_t eps, rng=None):
    # Estimate the rank of A by repeatedly using A.rmatvec(random vec)

    cdef int m = A.shape[0], n = A.shape[1], k = 0, kk = 0,r = n, krank
    cdef int no_of_cols = 4, intone = 1, info = 0
    cdef cnp.complex128_t[::1] tau = cnp.PyArray_ZEROS(1, [min(m, n)],
                                                       cnp.NPY_COMPLEX128, 0)
    cdef cnp.complex128_t[::1] y = cnp.PyArray_ZEROS(1, [n], cnp.NPY_COMPLEX128, 0)
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] retarr
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] x

    # The size of the QR decomposition is rank dependent which is unknown
    # at runtime. Hence we don't want to allocate a dense version of the
    # linear operator which can be too big. Instead, a typical "realloc double
    # if run out of space" strategy is used here. Starts with 4*n
    # Also, we hold the A.T @ x results in a separate array to return
    # and do the same for that too.
    cdef cnp.complex128_t *ra = <cnp.complex128_t*>PyMem_Malloc(
        sizeof(cnp.complex128_t)*no_of_cols*n
        )
    cdef cnp.complex128_t *reallocated_ra
    cdef cnp.complex128_t *ret = <cnp.complex128_t*>PyMem_Malloc(
        sizeof(cnp.complex128_t)*no_of_cols*n
        )
    cdef cnp.complex128_t *reallocated_ret
    cdef cnp.complex128_t enorm = 0.0

    if (not ra) or (not ret):
        raise MemoryError("Failed to allocate at least required memory "
                          f"{no_of_cols*n*8} bytes for"
                          "'scipy.linalg.interpolative.idz_findrank()' "
                          "function.")

    if not rng:
        rng = np.random.default_rng()

    krank = 0
    try:
        while True:

            # Generate random vector and rmatvec then save the result
            x = rng.uniform(size=(m,2)).view(np.complex128).ravel()
            y = A.rmatvec(x)

            for kk in range(n):
                ret[krank*n + kk] = y[kk]

            if krank == 0:
                enorm = dznrm2(&n, &y[0], &intone)
            else:  # krank > 0
                # Transpose-Apply previous Householder reflectors, if any
                # SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO
                zunm2r(<char*>'L',<char*>'C', &n, &intone, &krank, &ra[0], &n,
                       &tau[0], &y[0], &n, &ra[(no_of_cols-1)*n], &info)

            # Get the next Householder reflector
            r = n-krank
            # N, ALPHA, X, INCX, TAU
            zlarfgp(&r, &y[krank], &y[krank+1], &intone, &tau[krank])

            for kk in range(n):
                ra[krank*n + kk] = y[kk]

            # Running out of space; try to double the size of ra
            if krank == (no_of_cols-2):
                reallocated_ra = <cnp.complex128_t *>PyMem_Realloc(
                    ra, sizeof(cnp.complex128_t)*no_of_cols*n*2)
                reallocated_ret = <cnp.complex128_t *>PyMem_Realloc(
                    ret, sizeof(cnp.complex128_t)*no_of_cols*n*2)

                if reallocated_ra and reallocated_ret:
                    ra = reallocated_ra
                    ret = reallocated_ret
                    no_of_cols *= 2
                else:
                    raise MemoryError(
                        "'scipy.linalg.interpolative.idz_findrank()' failed to "
                        f"allocate the required memory,{no_of_cols*n*16} bytes "
                        "while trying to determine the rank (currently "
                        f"{krank}) of a LinearOperator with precision {eps}."
                    )
            krank += 1
            if (np.abs(y[krank-1]) < eps*enorm) or (krank >= min(m, n)):
                break
    finally:
        # Crashed or successfully ended up here
        # Discard Householder vectors
        PyMem_Free(ra)
        retarr = cnp.PyArray_EMPTY(2, [krank, n], cnp.NPY_COMPLEX128, 0)
        for k in range(krank):
            for kk in range(n):
                retarr[k, kk] = ret[k*n+kk]
        PyMem_Free(ret)

    return krank, retarr


def idz_id2svd(
    cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] cols,
    cnp.ndarray[cnp.int64_t, mode='c', ndim=1] perms,
    cnp.ndarray[cnp.complex128_t, ndim=2] proj,
    ):
    cdef int m = cols.shape[0], krank = cols.shape[1]
    cdef int n = proj.shape[1] + krank, info, ci
    cdef cnp.ndarray[cnp.complex128_t, mode='fortran', ndim=2] C
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] tau1
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] tau2
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] UU
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] S
    cdef cnp.ndarray[cnp.complex128_t, ndim=2] V
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] VV
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds1
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds2
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] p

    if krank > 0:
        UU = cnp.PyArray_ZEROS(2, [m, krank], cnp.NPY_COMPLEX128, 0)
        VV = cnp.PyArray_ZEROS(2, [n, krank], cnp.NPY_COMPLEX128, 0)
        p = cnp.PyArray_ZEROS(2, [krank, n], cnp.NPY_COMPLEX128, 0)

        # idd_reconint
        for ci in range(krank):
            p[ci, perms[ci]] = 1.0

        p[:, perms[krank:]] = proj[:, :]
        inds1, tau1 = idzr_qrpiv(cols, krank)
        # idz_rinqr and idz_rearr
        r = np.triu(cols[:krank, :])
        for ci in range(krank-1, -1, -1):
            r[:, [ci, inds1[ci]]] = r[:,  [inds1[ci], ci]]

        t = p.T.conj().copy()
        inds2, tau2 = idzr_qrpiv(t, krank)
        r2 = np.triu(t[:krank, :])
        for ci in range(krank-1, -1, -1):
            r2[:, [ci, inds2[ci]]] = r2[:,  [inds2[ci], ci]]

        r3 = r @ r2.T.conj()
        UU[:krank, :krank], S, V = la.svd(r3, full_matrices=False)

        # Apply Q of col to U from the left
        # But do the adjoint dance for LAPACK via U.H @ Q.H
        np.conjugate(tau1, out=tau1)
        C = cols[:, :krank].conj().copy(order='F')
        zunm2r(<char*>'R', <char*>'C',
            &krank, &m, &krank, &C[0, 0], &m, &tau1[0],
            &UU[0,0], &krank, &cols[0, 0], &info)

        VV[:krank, :krank] = V[:, :].conj().T

        # Apply Q of t to V from the left
        # But do the adjoint dance for LAPACK via V.H @ Q.H
        np.conjugate(tau2, out=tau2)
        C = t[:, :krank].conj().copy(order='F')
        zunm2r(<char*>'R', <char*>'C',
            &krank, &n, &krank, &C[0, 0], &n, &tau2[0],
            &VV[0, 0], &krank, &cols[0, 0], &info)

    return UU, S, VV


def idz_reconid(B, idx, proj):
    cdef int m = B.shape[0], krank = B.shape[1]
    cdef int n = len(idx)
    approx = np.zeros([m, n], dtype=np.complex128)

    approx[:, idx[:krank]] = B
    approx[:, idx[krank:]] = B @ proj

    return approx


def idz_snorm(A: LinearOperator, int its=20, rng=None):
    cdef int n = A.shape[1]
    cdef int j = 0, intone = 1
    cdef cnp.float64_t snorm = 0.0
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] v
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] u

    if not rng:
        rng = np.random.default_rng()

    v = rng.uniform(low=-1, high=1, size=(n, 2)).view(np.complex128).ravel()
    v /= dznrm2(&n, &v[0], &intone)

    for j in range(its):
        u = A.matvec(v)
        v = A.rmatvec(u)
        snorm = dznrm2(&n, &v[0], &intone)
        if snorm > 0.0:
            v /= snorm

        snorm = np.sqrt(snorm)

    return snorm


def idzp_aid(cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] a: NDArray, eps: float,
             rng=None):
    krank, proj = idz_estrank(a, eps=eps, rng=rng)
    if krank != 0:
        proj = proj[:krank, :]
        return idzp_id(proj, eps=eps)

    return idzp_id(a, eps=eps)


def idzp_asvd(cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] a, cnp.float64_t eps,
              rng=None):
    cdef int m = a.shape[0], n = a.shape[1]
    cdef int krank, info, ci
    cdef cnp.ndarray[cnp.complex128_t, mode='fortran', ndim=2] C
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] tau1
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] tau2
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] UU
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] S
    cdef cnp.ndarray[cnp.complex128_t, ndim=2] V
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] VV
    cdef cnp.ndarray[cnp.complex128_t, ndim=2] proj
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] perms
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds1
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds2
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] p
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] col

    krank, perms, proj = idzp_aid(a.copy(), eps, rng)

    if krank > 0:
        UU = cnp.PyArray_ZEROS(2, [m, krank], cnp.NPY_COMPLEX128, 0)
        VV = cnp.PyArray_ZEROS(2, [n, krank], cnp.NPY_COMPLEX128, 0)
        p = cnp.PyArray_ZEROS(2, [krank, n], cnp.NPY_COMPLEX128, 0)
        col = a[:, perms[:krank]].copy()

        # idd_reconint
        for ci in range(krank):
            p[ci, perms[ci]] = 1.0

        p[:, perms[krank:]] = proj[:, :]
        inds1, tau1 = idzr_qrpiv(col, krank)
        # idz_rinqr and idz_rearr
        r = np.triu(col[:krank, :])
        for ci in range(krank-1, -1, -1):
            r[:, [ci, inds1[ci]]] = r[:,  [inds1[ci], ci]]

        t = p.T.conj().copy()
        inds2, tau2 = idzr_qrpiv(t, krank)
        r2 = np.triu(t[:krank, :])
        for ci in range(krank-1, -1, -1):
            r2[:, [ci, inds2[ci]]] = r2[:,  [inds2[ci], ci]]

        r3 = r @ r2.T.conj()
        UU[:krank, :krank], S, V = la.svd(r3, full_matrices=False)

        # Apply Q of col to U from the left
        # But do the adjoint dance for LAPACK via U.H @ Q.H
        np.conjugate(tau1, out=tau1)
        C = col[:, :krank].conj().copy(order='F')
        zunm2r(<char*>'R', <char*>'C',
            &krank, &m, &krank, &C[0, 0], &m, &tau1[0],
            &UU[0,0], &krank, &a[0, 0], &info)

        VV[:krank, :krank] = V[:, :].conj().T

        # Apply Q of t to V from the left
        # But do the adjoint dance for LAPACK via V.H @ Q.H
        np.conjugate(tau2, out=tau2)
        C = t[:, :krank].conj().copy(order='F')
        zunm2r(<char*>'R', <char*>'C',
            &krank, &n, &krank, &C[0, 0], &n, &tau2[0],
            &VV[0, 0], &krank, &a[0, 0], &info)

    return UU, S, VV


def idzp_id(cnp.ndarray[cnp.complex128_t, mode="c", ndim=2] a, cnp.float64_t eps):
    cdef int n = a.shape[1], krank, tmp_int, p
    cdef double complex one = 1
    krank, _, inds = idzp_qrpiv(a, eps)

    # Change pivots to permutation
    perms = cnp.PyArray_Arange(0, n, 1, cnp.NPY_INT64)

    if krank > 0:
        for p in range(krank):
            # Apply pivots
            tmp_int = perms[p]
            perms[p] = perms[inds[p]]
            perms[inds[p]] = tmp_int

    tmp_int = n - krank
    # SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB
    ztrsm(<char*>'R', <char*>'L', <char*>'N', <char*>'N',
          &tmp_int, &krank, &one, &a[0, 0], &n, &a[0, krank], &n)

    return krank, perms, a[:krank, krank:]


def idzp_qrpiv(cnp.ndarray[cnp.complex128_t, mode="c", ndim=2] a, cnp.float64_t eps):
    cdef int m = a.shape[0], n = a.shape[1]
    cdef cnp.ndarray col_norms = cnp.PyArray_ZEROS(1, [n], cnp.NPY_FLOAT64, 0)
    cdef int k = 0, kpiv = 0, i = 0, tmp_int = 0, int_n = 0
    cdef double complex tmp_sca = 0.
    cdef cnp.ndarray taus = cnp.PyArray_ZEROS(1, [m], cnp.NPY_COMPLEX128, 0)
    cdef cnp.ndarray ind = cnp.PyArray_ZEROS(1, [n], cnp.NPY_INT64, 0)
    cdef double complex[::1] taus_v = taus
    cdef cnp.float64_t feps = 0.1e-16  # Smaller than np.finfo(np.float64).eps
    cdef cnp.float64_t ssmax, ssmaxin
    cdef int nupdate = 0

    for i in range(n):
        col_norms[i] = dznrm2(&m, &a[0, i], &n)**2

    kpiv = np.argmax(col_norms)
    ssmax = col_norms[kpiv]
    ssmaxin = ssmax

    for k in range(min(m, n)):

        # Pivoting
        ind[k] = kpiv
        # Swap columns a[:, k] and a[:, kpiv]
        a[:, [kpiv, k]] = a[:, [k, kpiv]]

        # Swap col_norms[krank] and col_norms[kpiv]
        col_norms[[kpiv, k]] = col_norms[[k, kpiv]]

        if k < m-1:
            # Compute the householder reflector for column k
            tmp_sca = a[k, k]
            # FIX: Convert these to F_INT
            tmp_int = <int>(m - k)
            int_n = <int>n
            zlarfgp(&tmp_int, &tmp_sca, &a[k+1, k], &int_n, &taus_v[k])

            # Overwrite with 1. for easy matmul
            a[k, k] = 1.0
            if k < n-1:
                # Apply the householder reflector to the rest on the right.
                # Note! Tau returned by zlarfgp is complex valued and thus,
                # reflector is not Hermitian, hence the conjugates. See the
                # documentation of zlarfgp.
                a[k:, k+1:] -= np.outer(taus[k].conj()*a[k:, k],
                                        a[k:, k].conj() @ a[k:, k+1:]
                                        )

            # Put back the beta in place
            a[k, k] = tmp_sca
            # Update the norms
            col_norms[k] = 0
            col_norms[k+1:] -= (a[k, k+1:] * a[k, k+1:].conj()).real
            ssmax = 0.0
            kpiv = k+1

            if k < n-1:
                kpiv = np.argmax(col_norms[k+1:]) + (k + 1)
                ssmax = col_norms[kpiv]

            if (((ssmax < 1000*feps*ssmaxin) and (nupdate == 0)) or
                    ((ssmax < ((1000*feps)**2)*ssmaxin) and (nupdate == 1))):
                nupdate += 1
                ssmax = 0
                kpiv = k+1
                if k < n-1:
                    for i in range(k+1, n):
                        tmp_int = m-k-1
                        col_norms[i] = dznrm2(&tmp_int, &a[k+1, i], &n)**2
                    kpiv = np.argmax(col_norms[k+1:]) + (k + 1)
                    ssmax = col_norms[kpiv]
        if (ssmax <= (eps**2)*ssmaxin):
            break
    # a is overwritten; return numerical rank and pivots

    return k+1, taus, ind


def idzp_rid(A: LinearOperator, cnp.float64_t eps, rng=None):
    _, ret = idz_findrank(A, eps, rng=rng)
    return idzp_id(ret, eps=eps)


def idzp_rsvd(A: LinearOperator, cnp.float64_t eps, rng=None):
    cdef int n = A.shape[1]
    cdef int krank, j
    cdef cnp.ndarray[cnp.int64_t, mode='c', ndim=1] perms
    cdef cnp.ndarray[cnp.complex128_t, ndim=2] proj
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] col
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] x

    krank, perms, proj = idzp_rid(A, eps, rng=rng)

    if krank > 0:
        # idd_getcols
        col = cnp.PyArray_EMPTY(2, [n, krank], cnp.NPY_COMPLEX128, 0)
        x = cnp.PyArray_ZEROS(1, [n], cnp.NPY_COMPLEX128, 0)

        for j in range(krank):
            x[perms[j]] = 1.
            col[:, j] = A.matvec(x)
            x[perms[j]] = 0.

        return idz_id2svd(cols=col, perms=perms, proj=proj)

    # TODO: figure out empty return
    return None


def idzp_svd(cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] a, cnp.float64_t eps):
    cdef int m = a.shape[0], krank, info
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] taus
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] UU
    cdef cnp.ndarray[cnp.complex128_t, ndim=2] V
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] r
    cdef cnp.ndarray[cnp.complex128_t, mode='fortran', ndim=2] C
    cdef cnp.ndarray[cnp.float64_t, ndim=1] S

    # Get the pivoted QR
    krank, taus, inds = idzp_qrpiv(a, eps)
    UU = cnp.PyArray_ZEROS(2, [m, krank], cnp.NPY_COMPLEX128, 0)

    if krank > 0:
        r = np.triu(a[:krank, :])

        for p in range(krank-1, -1, -1):
            r[:, [p, inds[p]]] = r[:, [inds[p], p]]

        UU[:krank, :krank], S, V = la.svd(r, full_matrices=False)
        # Apply Q to U via zunm2r
        np.conjugate(taus, out=taus)
        # But do the adjoint dance for LAPACK via U.H @ Q.H; use a for scratch
        C = a[:, :krank].conj().copy(order='F')
        zunm2r(<char*>'R', <char*>'C',
               &krank, &m, &krank, &C[0, 0], &m, &taus[0],
               &UU[0,0], &krank, &a[0, 0], &info)

    return UU, S, V


def idzr_aid(cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] a: NDArray, int krank,
             rng=None):
    cdef int m = a.shape[0], n2, L, nblock, nsteps = 3, mb
    cdef cnp.float64_t twopi = 2*np.pi, fact
    cdef double complex twopii = twopi*1.j
    cdef cnp.ndarray[cnp.npy_int64, mode='c', ndim=1] ind
    cdef cnp.ndarray[cnp.npy_int64, mode='c', ndim=1] subselect
    cdef cnp.ndarray[cnp.npy_float64, mode='c', ndim=1] dm1
    cdef cnp.ndarray[cnp.npy_float64, mode='c', ndim=1] dm2
    cdef cnp.ndarray[cnp.npy_float64, mode='c', ndim=3] albetas
    cdef cnp.ndarray[cnp.npy_float64, mode='c', ndim=2] rta
    cdef cnp.ndarray[cnp.npy_float64, mode='c', ndim=2] giv2x2

    if not rng:
        rng = np.random.default_rng()

    n2 = 0
    L = krank + 8
    if (L >= n2) or (L > m):
        inds, proj = idzr_id(a, krank)
        return inds, proj

    n2 = idd_poweroftwo(m)
    # This part is the initialization that is done via idz_frmi
    # for a Subsampled Randomized Fourier Transfmrom (SRFT).

    # Draw (nsteps x m x 4) array from [0, 2)*pi uniformly for
    # random points on complex unit circle and unitary rotations
    albetas = np.empty([nsteps, m, 4])
    albetas[:, :, 2:] = rng.uniform(low=0.0, high=2.0, size=[nsteps, m, 2])
    albetas[:, :, 2:] *= np.pi
    np.cos(albetas[:, :, 2], out=albetas[:, :, 0])
    np.sin(albetas[:, :, 2], out=albetas[:, :, 1])
    np.cos(albetas[:, :, 3], out=albetas[:, :, 2])
    np.sin(albetas[:, :, 3], out=albetas[:, :, 3])

    # idd_random_transf
    rta = a.copy()

    # Rotate and shuffle "a" nsteps-many times
    giv2x2 = np.array([[0., 0. ], [0., 0.]])
    for nstep in range(nsteps):
        # Multiply with a point on the unit circle
        rta *= albetas[nstep, :, 2:].view(np.complex128)
        # Rotate
        for row in range(m-1):
            alpha, beta = albetas[nstep, row, 0], albetas[nstep, row, 1]
            giv2x2[0, 0] = alpha
            giv2x2[0, 1] = beta
            giv2x2[1, 0] = -beta
            giv2x2[1, 1] = alpha
            np.matmul(giv2x2, rta[row:row+2, :], out=rta[row:row+2, :])

        rta = rta[rng.permutation(m), :]

    # idd_subselect pick randomly n2-many rows
    subselect = rng.choice(m, n2, replace=False)
    rta = rta[subselect, :]
    ind = rng.choice(n2, L, replace=False)

    nblock = idd_ldiv(L, n2)
    mb = n2 // nblock
    fact = 1.0 / np.sqrt(n2)

    # Create (L x mb) DFT matrix
    # wsave = np.empty([L, mb], dtype=np.complex128)
    dm1, dm2 = np.divmod(ind, mb, dtype=np.float64)
    dm1 /= n2
    dm1 += dm2 / mb
    wsave = np.outer(dm1, -twopii*np.arange(mb))
    np.exp(wsave, out=wsave)
    wsave *= fact

    # Perform partial FFT to each nblock then swap first two axes for transposition
    # and subsample by ind // mb. This is basically a few options combined into one
    # First we view each column as (nblock x mb) then take fft of each mb-long chunk.
    # Then we transpose and multiply with DFT matrix and subselect.
    # See DOI:10.1016/j.acha.2007.12.002 - Section 3.3

    # Original fortran code does this single column at a time. We do a bit of array
    # manipulation to do it in one go for all columns at once.
    F = np.swapaxes(
          fft(rta.reshape(nblock, mb, -1, order='F'), axis=0), 0, 1
          )[:, ind // mb, :]
    # Perform direct calculation with DFT matrix
    V = np.einsum('ij,jim->im', wsave, F)

    return idzr_id(V, krank)


def idzr_asvd(cnp.ndarray[cnp.complex128_t, mode="c", ndim=2] a, int krank, rng=None):
    cdef int m = a.shape[0], n = a.shape[1]
    cdef int info, ci
    cdef cnp.ndarray[cnp.complex128_t, mode='fortran', ndim=2] C
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] tau1
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] tau2
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] UU
    cdef cnp.ndarray[cnp.float64_t, mode='c', ndim=1] S
    cdef cnp.ndarray[cnp.complex128_t, ndim=2] V
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] VV
    cdef cnp.ndarray[cnp.complex128_t, ndim=2] proj
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] perms
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds1
    cdef cnp.ndarray[cnp.npy_int64, ndim=1] inds2
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] p
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] col
    UU = cnp.PyArray_ZEROS(2, [m, krank], cnp.NPY_COMPLEX128, 0)
    VV = cnp.PyArray_ZEROS(2, [n, krank], cnp.NPY_COMPLEX128, 0)
    p = cnp.PyArray_ZEROS(2, [krank, n], cnp.NPY_COMPLEX128, 0)

    perms, proj = idzr_aid(a.copy(), krank=krank, rng=rng)
    col = a[:, perms[:krank]].copy()

    # idd_reconint
    for ci in range(krank):
        p[ci, perms[ci]] = 1.0

    p[:, perms[krank:]] = proj[:, :]
    inds1, tau1 = idzr_qrpiv(col, krank)
    # idz_rinqr and idz_rearr
    r = np.triu(col[:krank, :])
    for ci in range(krank-1, -1, -1):
        r[:, [ci, inds1[ci]]] = r[:,  [inds1[ci], ci]]

    t = p.T.conj().copy()
    inds2, tau2 = idzr_qrpiv(t, krank)
    r2 = np.triu(t[:krank, :])
    for ci in range(krank-1, -1, -1):
        r2[:, [ci, inds2[ci]]] = r2[:,  [inds2[ci], ci]]

    r3 = r @ r2.T.conj()
    UU[:krank, :krank], S, V = la.svd(r3, full_matrices=False)

    # Apply Q of col to U from the left
    # But do the adjoint dance for LAPACK via U.H @ Q.H
    np.conjugate(tau1, out=tau1)
    C = col[:, :krank].conj().copy(order='F')
    zunm2r(<char*>'R', <char*>'C',
           &krank, &m, &krank, &C[0, 0], &m, &tau1[0],
           &UU[0,0], &krank, &a[0, 0], &info)

    VV[:krank, :krank] = V[:, :].conj().T

    # Apply Q of t to V from the left
    # But do the adjoint dance for LAPACK via V.H @ Q.H
    np.conjugate(tau2, out=tau2)
    C = t[:, :krank].conj().copy(order='F')
    zunm2r(<char*>'R', <char*>'C',
           &krank, &n, &krank, &C[0, 0], &n, &tau2[0],
           &VV[0, 0], &krank, &a[0, 0], &info)

    return UU, S, VV


def idzr_id(cnp.ndarray[cnp.complex128_t, ndim=2] a, int krank):
    cdef int n = a.shape[1], tmp_int, p
    cdef double complex one = 1.0
    cdef cnp.ndarray[cnp.int64_t, ndim=1] inds
    cdef cnp.ndarray[cnp.int64_t, ndim=1] perms

    inds, _ = idzr_qrpiv(a, krank)
    perms = cnp.PyArray_Arange(0, n, 1, cnp.NPY_INT64)

    if krank > 0:
        for p in range(krank):
            # Apply pivots
            tmp_int = perms[p]
            perms[p] = perms[inds[p]]
            perms[inds[p]] = tmp_int
    tmp_int = n - krank
    # SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB
    ztrsm(<char*>'R', <char*>'L', <char*>'N', <char*>'N',
          &tmp_int, &krank, &one, &a[0, 0], &n, &a[0, krank], &n)

    return perms, a[:krank, krank:]


def idzr_qrpiv(cnp.ndarray[cnp.complex128_t, mode="c", ndim=2] a, int krank):
    cdef int m = a.shape[0], n = a.shape[1]
    cdef int loop = 0, loops, kpiv = 0, i = 0, tmp_int = 0
    cdef cnp.ndarray col_norms = cnp.PyArray_ZEROS(1, [n], cnp.NPY_FLOAT64, 0)
    cdef double complex tmp_sca = 0.
    cdef cnp.ndarray taus = cnp.PyArray_ZEROS(1, [m], cnp.NPY_COMPLEX128, 0)
    cdef cnp.ndarray ind = cnp.PyArray_ZEROS(1, [n], cnp.NPY_INT64, 0)
    cdef double complex[::1] taus_v = taus
    cdef cnp.float64_t feps = 0.1e-16  # Smaller than np.finfo(np.float64).eps
    cdef cnp.float64_t ssmax, ssmaxin
    cdef int nupdate = 0

    loops = min(krank, min(m, n))
    for i in range(n):
        col_norms[i] = dznrm2(&m, &a[0, i], &n)**2

    kpiv = np.argmax(col_norms)
    ssmax = col_norms[kpiv]
    ssmaxin = ssmax

    for loop in range(loops):

        ind[loop] = kpiv
        # Swap columns a[:, k] and a[:, kpiv]
        a[:, [kpiv, loop]] = a[:, [loop, kpiv]]
        # Swap col_norms[krank] and col_norms[kpiv]
        col_norms[[kpiv, loop]] = col_norms[[loop, kpiv]]

        if loop < m-1:
            tmp_sca = a[loop, loop]
            # FIX: Convert these to F_INT
            tmp_int = (m - loop)
            zlarfgp(&tmp_int, &tmp_sca, &a[loop+1, loop], &n, &taus_v[loop])

            # Overwrite with 1. for easy matmul
            a[loop, loop] = 1
            if loop < n-1:
                # Apply the householder reflector to the rest on the right
                a[loop:, loop+1:] -= np.outer(
                    np.conj(taus[loop])*a[loop:, loop],
                    a[loop:, loop].conj() @ a[loop:, loop+1:]
                    )
            # Put back the beta in place
            a[loop, loop] = tmp_sca

            # Update the norms
            col_norms[loop] = 0
            col_norms[loop+1:] -= (a[loop, loop+1:]*a[loop, loop+1:].conj()).real
            ssmax = 0
            kpiv = loop+1

            if loop < n-1:
                kpiv = np.argmax(col_norms[loop+1:]) + (loop + 1)
                ssmax = col_norms[kpiv]
            if (((ssmax < 1000*feps*ssmaxin) and (nupdate == 0)) or
                    ((ssmax < ((1000*feps)**2)*ssmaxin) and (nupdate == 1))):
                nupdate += 1
                ssmax = 0
                kpiv = loop+1

                if loop < n-1:
                    for i in range(loop+1, n):
                        tmp_int = m-loop-1
                        col_norms[i] = dznrm2(&tmp_int, &a[loop+1, i], &n)**2
                    kpiv = np.argmax(col_norms[loop+1:]) + (loop + 1)
                    ssmax = col_norms[kpiv]

    return ind, taus


def idzr_rid(A: LinearOperator, int krank, rng=None):
    cdef int m = A.shape[0], n = A.shape[1], k = 0
    cdef int L = min(krank+2, min(m, n))
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] r

    if not rng:
        rng = np.random.default_rng()

    r = cnp.PyArray_EMPTY(2, [L, n], cnp.NPY_COMPLEX128, 0)
    for k in range(L):
        r[k, :] = A.rmatvec(rng.uniform(size=(m,2)).view(np.complex128).ravel())

    return idzr_id(a=r.conj(), krank=krank)


def idzr_rsvd(A: LinearOperator, int krank, rng=None):
    cdef int n = A.shape[1], j
    cdef cnp.ndarray[cnp.int64_t, mode='c', ndim=1] perms
    cdef cnp.ndarray[cnp.complex128_t, ndim=2] proj
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] col

    perms, proj = idzr_rid(A, krank, rng)
    # idd_getcols
    col = cnp.PyArray_EMPTY(2, [n, krank], cnp.NPY_COMPLEX128, 0)
    x = cnp.PyArray_ZEROS(1, [n], cnp.NPY_COMPLEX128, 0)
    for j in range(krank):
        x[perms[j]] = 1.
        col[:, j] = A.matvec(x)
        x[perms[j]] = 0.

    return idz_id2svd(cols=col, perms=perms, proj=proj)


def idzr_svd(cnp.ndarray[cnp.complex128_t, mode="c", ndim=2] a, int krank):
    cdef int m = a.shape[0], n = a.shape[1], info = 0
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=1] taus
    cdef cnp.ndarray[cnp.int64_t, mode='c', ndim=1] inds
    cdef cnp.ndarray[cnp.complex128_t, mode='c', ndim=2] UU
    cdef cnp.ndarray[cnp.complex128_t, mode='fortran', ndim=2] C
    UU = cnp.PyArray_ZEROS(2, [m, krank], cnp.NPY_COMPLEX128, 0)

    krank = min(krank, min(m, n))
    # Get the pivoted QR
    inds, taus = idzr_qrpiv(a, krank)
    r = np.triu(a[:krank, :])
    # Apply pivots in reverse
    for p in range(krank-1, -1, -1):
        r[:, [p, inds[p]]] = r[:, [inds[p], p]]

    # JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO
    # zgesvd()
    UU[:krank, :krank], S, V = la.svd(r, full_matrices=False)

    # Apply Q to U via zunm2r
    np.conjugate(taus, out=taus)
    # But do the adjoint dance for LAPACK via U.H @ Q.H; use a for scratch
    C = a[:, :krank].conj().copy(order='F')
    zunm2r(<char*>'R', <char*>'C',
           &krank, &m, &krank, &C[0, 0], &m, &taus[0],
           &UU[0,0], &krank, &a[0, 0], &info)

    return UU, S, V
