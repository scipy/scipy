"""
Cythonized kernels of _qmnvt.py
"""

from libc.math cimport INFINITY, modf
from numpy cimport npy_intp

import numpy as np

from scipy.special.cython_special cimport ndtri as phinv
from scipy.special.cython_special cimport ndtr as phi

import cython


#pythran export _qmvn_inner(float[:], float[:, :], int64, int64, float[:, :], float[:], float[:])  # noqa: E501
@cython.boundscheck(False)
@cython.cdivision(True)
def _qmvn_inner(double[::1] q, double[:, ::1] rndm, int n_qmc_samples, int n_batches, double[:, ::1] cho, double[::1] lo, double[::1] hi):

    cdef:
        npy_intp  n = cho.shape[0]
        npy_intp i, j, k
        double prob, error_var, r, qq, z_k, res_k, scrtch, lo_i, hi_i, ct, c_k, d_k, s_k, ci

        double[::1] c, d, dc, pv, n_qmc_zeros
        double[:, ::1] y

    ct = cho[0, 0]
    ci = phi[double](lo[0] / ct)
    dci = phi(hi[0] / ct) - phi(lo[0] / ct)
    prob = 0.0
    error_var = 0.0

    y = np.zeros((n - 1, n_qmc_samples))
    n_qmc_zeros = np.zeros(n_qmc_samples)

    for j in range(n_batches):
        c = np.full(n_qmc_samples, ci)
        dc = np.full(n_qmc_samples, dci)
        pv = dc.copy()
        d = n_qmc_zeros.copy()

        for i in range(1, n):
            r = rndm[j, i]
            qq = q[i - 1]
            ct = cho[i, i]
            lo_i, hi_i = lo[i], hi[i]
            for k in range(n_qmc_samples):
                # Pseudorandomly-shifted lattice coordinate.
                z_k = qq * (k + 1) + r
                z_k = modf(z_k, &scrtch)   # remainder(z, 1.0)

                # Tent periodization transform.
                x_k = abs(2 * z_k - 1)

                res_k = phinv(c[k] + x_k*dc[k])

                # work around inf - inf = nan
                if res_k == INFINITY:
                    res_k = 1e100
                elif res_k == -INFINITY:
                    res_k = -1e100

                y[i-1, k] = res_k

                s_k = 0.0
                for ii in range(i):
                    s_k += cho[i, ii] * y[ii, k]

                c[k] = phi((lo_i - s_k) / ct)
                d[k] = phi((hi_i - s_k) / ct)
                dc[k] = d[k] - c[k]
                pv[k] = pv[k] * dc[k]

        # Accumulate the mean and error variances with online formulations.
        dd = (np.mean(pv) - prob) / (j + 1)
        prob += dd
        error_var = (j - 1) * error_var / (j + 1) + dd * dd

    # Error bounds are 3 times the standard error of the estimates.
    est_error = 3 * np.sqrt(error_var)
    n_samples = n_qmc_samples * n_batches

    return prob, est_error, n_samples


