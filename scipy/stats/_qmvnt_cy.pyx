"""
Cythonized kernels of _qmnvt.py
"""
import numpy as np
import cython

from scipy.special import gammaincinv

from libc.math cimport INFINITY, modf
from numpy cimport npy_intp
from scipy.special.cython_special cimport ndtri as phinv
from scipy.special.cython_special cimport ndtr as phi


#pythran export _qmvn_inner(float[:], float[:, :], int64, int64, float[:, :], float[:], float[:])  # noqa: E501
@cython.boundscheck(False)
@cython.cdivision(True)
def _qmvn_inner(double[::1] q,
                double[:, ::1] rndm,
                int n_qmc_samples,
                int n_batches,
                double[:, ::1] cho,
                double[::1] lo,
                double[::1] hi):

    cdef:
        npy_intp n = cho.shape[0]
        npy_intp i, j, k
        double prob, error_var, r, qq, z_k, res_k, scrtch, lo_i, hi_i, ct, s_k, ci

        double[::1] c, d, dc, pv, n_qmc_zeros
        double[:, ::1] y

    ct = cho[0, 0]
    ci = phi(lo[0] / ct)
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


#pythran export _qmvt_inner(float[], float[:, :], int, int, float[:, :], float[], float[], float)  # noqa: E501
def _qmvt_inner(double[::1] q,
                double[:, ::1] rndm,
                int n_qmc_samples,
                int n_batches,
                double[:, ::1] cho,
                double[::1] lo,
                double [::1] hi,
                double nu):

    cdef:
        double prob = 0.0
        double error_var = 0.0
        double scrtch, qq, rr, z_k, c_k, d_k, lois_k, hiis_k
        int n = cho.shape[0]
        int i, j, k

        double[::1] y = np.zeros(n_qmc_samples)
        double[::1] c = np.ones(n_qmc_samples)
        double[::1] d = np.ones(n_qmc_samples)
        double[::1] dc = np.zeros(n_qmc_samples)
        double[::1] pv = np.zeros(n_qmc_samples)
        double[:, ::1] s = np.zeros((n, n_qmc_samples))
        double[::1] r

    i_samples = np.arange(n_qmc_samples) + 1

    for j in range(n_batches):
        pv = np.ones(n_qmc_samples)
        s = np.zeros((n, n_qmc_samples))

        # i == 0 special actions: We'll use one of the QR variates to pull out the
        # t-distribution scaling.
        z = q[0]*i_samples + rndm[j, 0]
        z -= z.astype(int)
        x = abs(2*z - 1)
        r = np.sqrt(2 * gammaincinv(nu / 2, x)) if nu > 0 else np.ones_like(x)

        for i in range(0, n):
            qq = q[i]
            rr = rndm[j, i]
            for k in range(n_qmc_samples):
                # Pseudorandomly-shifted lattice coordinate.
                z_k = qq * (k + 1) + rr
                z_k = modf(z_k, &scrtch)   # remainder(z, 1.0)
                # Tent periodization transform.
                x_k = abs(2 * z_k - 1)

                if i > 0:
                    y[k] = phinv(c[k] + x_k*dc[k])

                    # s[i:, :] += cho[i:, i - 1][:, np.newaxis] * y
                    for ip in range(i, n):
                        s[ip, k] += cho[ip, i - 1] * y[k]

            for k in range(n_qmc_samples):
                lois_k = lo[i] * r[k] - s[i, k]
                if lois_k < -9:
                    c_k = 0.0
                elif lois_k < 9:
                    c_k = phi(lois_k)
                else:
                    c_k = 1.0

                hiis_k = hi[i] * r[k] - s[i, k]
                if hiis_k < -9:
                    d_k = 0.0
                elif hiis_k < 9:
                    d_k = phi(hiis_k)
                else:
                    d_k = 1.0

                c[k] = c_k
                d[k] = d_k

                dc[k] = (d_k - c_k)
                pv[k] *= dc[k]

        # Accumulate the mean and error variances with online formulations.
        dd = (np.mean(pv) - prob) / (j + 1)
        prob += dd
        error_var = (j - 1) * error_var / (j + 1) + dd * dd

    # Error bounds are 3 times the standard error of the estimates.
    est_error = 3.0 * np.sqrt(error_var)
    n_samples = n_qmc_samples * n_batches

    return prob, est_error, n_samples
