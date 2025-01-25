"""
Pythranized kernels of _qmnvt.py
"""
import numpy as np

from scipy.special import gammaincinv, ndtr as phi, ndtri as phinv


#pythran export _qmvt_inner(float[], float[:, :], int, int, float[:, :], float[], float[], float)  # noqa: E501
def _qmvt_inner(q, rndm, n_qmc_samples, n_batches, cho, lo, hi, nu):

    prob = 0.0
    error_var = 0.0

    i_samples = np.arange(n_qmc_samples) + 1
    n = cho.shape[0]

    dc = np.zeros(n_qmc_samples)
    c = np.ones(n_qmc_samples)
    d = np.ones(n_qmc_samples)

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
            # Pseudorandomly-shifted lattice coordinate.
            z = q[i] * i_samples + rndm[j, i]
            # Fast remainder(z, 1.0)
            z -= z.astype(int)
            # Tent periodization transform.
            x = abs(2 * z - 1)
            if i > 0:
                y = phinv(c + x * dc)
                s[i:, :] += cho[i:, i - 1][:, np.newaxis] * y
            si = s[i, :]

            c = np.ones(n_qmc_samples)
            d = np.ones(n_qmc_samples)
            lois = lo[i] * r - si
            hiis = hi[i] * r - si
            c[lois < -9] = 0.0
            d[hiis < -9] = 0.0
            lo_mask = abs(lois) < 9
            hi_mask = abs(hiis) < 9
            c[lo_mask] = phi(lois[lo_mask])
            d[hi_mask] = phi(hiis[hi_mask])

            dc = d - c
            pv *= dc

        # Accumulate the mean and error variances with online formulations.
        d = (np.mean(pv) - prob) / (j + 1)
        prob += d
        error_var = (j - 1) * error_var / (j + 1) + d * d

    # Error bounds are 3 times the standard error of the estimates.
    est_error = 3.0 * np.sqrt(error_var)
    n_samples = n_qmc_samples * n_batches

    return prob, est_error, n_samples


#pythran export _qmvn_inner(float[:], float[:, :], int, int, float[:, :], float[:], float[:])  # noqa: E501
def _qmvn_inner(q, rndm, n_qmc_samples, n_batches, cho, lo, hi):

    n = cho.shape[0]
    ct = cho[0, 0]
    c = phi(lo[0] / ct)
    d = phi(hi[0] / ct)
    ci = c
    dci = d - ci
    prob = 0.0
    error_var = 0.0

    n = cho.shape[0]
    y = np.zeros((n - 1, n_qmc_samples))
    i_samples = np.arange(n_qmc_samples) + 1

    for j in range(n_batches):
        c = np.full(n_qmc_samples, ci)
        dc = np.full(n_qmc_samples, dci)
        pv = dc.copy()
        for i in range(1, n):
            # Pseudorandomly-shifted lattice coordinate.
            z = q[i - 1] * i_samples + rndm[j, i]  #rng.random()
            # Fast remainder(z, 1.0)
            z -= z.astype(int)
            # Tent periodization transform.
            x = abs(2 * z - 1)

            # work around inf - inf = nan in the matmul below
            res =  phinv(c + x * dc)
            res[res == -np.inf] = -1e100
            res[res == np.inf] = 1e100
            y[i - 1, :] = res

            # s = cho[i, :i] @ y[:i, :] : XXX unroll for now (BLAS)
            s = np.zeros(n_qmc_samples)
            for pp in range(n_qmc_samples):
                for qq in range(i):
                    s[pp] += cho[i, qq] * y[qq, pp]

            ct = cho[i, i]
            c = phi((lo[i] - s) / ct)
            d = phi((hi[i] - s) / ct)
            dc = d - c
            pv = pv * dc

        # Accumulate the mean and error variances with online formulations.
        d = (pv.mean() - prob) / (j + 1)
        prob += d
        error_var = (j - 1) * error_var / (j + 1) + d * d

    # Error bounds are 3 times the standard error of the estimates.
    est_error = 3 * np.sqrt(error_var)
    n_samples = n_qmc_samples * n_batches

    return prob, est_error, n_samples


