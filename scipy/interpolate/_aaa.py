# Copyright (c) 2017, The Chancellor, Masters and Scholars of the University
# of Oxford, and the Chebfun Developers. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the University of Oxford nor the names of its
#       contributors may be used to endorse or promote products derived from
#       this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from functools import cached_property
import warnings

import numpy as np
import scipy


__all__ = ["AAA"]


class AAA:
    r"""
    AAA real or complex rational approximation.

    As described in [1]_, the AAA algorithm is a greedy algorithm for approximation by
    rational functions on a real or complex set of points. The rational approximation is
    represented in a barycentric form from which the roots (zeros), poles, and residues
    can be computed.

    Parameters
    ----------
    f : array_like or callable
        Function values at `z` or a function that will be called ``f(z)``.
    z : array_like
        Values at which `f` is provided.
    rtol : float, optional
        Relative tolerance, defaults to 1e-13.
    max_terms : int, optional
        Maximum number of terms in the barycentric representation, defaults to 100.

    Attributes
    ----------
    support_points : array
        Support points of the approximation.
    values : array
        Value of the approximation at the `support_points`.
    weights : array
        Weights of the barycentric approximation.
    errors : array
        Error in the successive iterations of AAA.
    poles : array
        Poles of the AAA approximation.
    residues : array
        Residues associated with the poles of the approximation.
    roots : array
        Roots (zeros) of the AAA approximation.

    Warns
    -----
    RuntimeWarning
        If `rtol` is not achieved in `max_terms` iterations.

    See Also
    --------
    pade : Padé approximation

    Notes
    -----
    At the :math:`m` th iteration, the rational approximation in the AAA algorithm takes
    the barycentric form

    .. math::

        r(z) = \frac{n(z)}{d(z)} =
        \frac{\sum_{j=1}^m\frac{w_jf_j}{z-z_j}}{\sum_{j=1}^m\frac{w_j}{z-z_j}},

    where :math:`z_1,\dots,z_m` are real or complex support points selected from `z`,
    :math:`f_1,\dots,f_m` are a set of real or complex data values, and
    :math:`w_1,\dots,w_m` are real or complex weights. The algorithm then proceeds to
    select the next support point :math:`z_{m+1}` is selected from the remaining
    unselected points in `z` such that the nonlinear residual :math:`|f(z) - n(z)/d(z)|`
    is maximised. The weights are selected to solve the least-squares problem

    .. math::

        \text{minimise}\lVert fd - n \rVert \quad \text{subject to} \quad
        \sum_{i=1}^{m+1}w_i = 1,

    over the unselected points in `z`.

    References
    ----------
    .. [1] Y. Nakatsukasa, O. Sete, and L. N. Trefethen, "The AAA algorithm for
            rational approximation", SIAM J. Sci. Comp. 40 (2018), A1494-A1522.

    Examples
    --------

    Here we reproduce a number of the numerical examples from [1]_ as demonstration
    of the functionality offered by this method.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.interpolate import AAA

    For the first example we approximate the gamma function on ``[-3.5, 4.5]`` by
    extrapolating from 100 samples in ``[-1.5, 1.5]``.

    >>> from scipy.special import gamma
    >>> sample_points = np.linspace(-1.5, 1.5, num=100)
    >>> r = AAA(gamma, sample_points)
    >>> z = np.linspace(-3.5, 4.5, num=1000)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(z, gamma(z), label="Gamma")
    >>> ax.plot(sample_points, gamma(sample_points), label="Sample points")
    >>> ax.plot(z, r(z).real, '--', label="AAA approximation")
    >>> ax.set(xlabel="z", ylabel="r(z)", ylim=[-8, 8], xlim=[-3.5, 4.5])
    >>> ax.legend()
    >>> plt.show()

    We can also view the poles of the rational approximation and their residue:

    >>> for i in np.argsort(r.poles):
    ...     print(f"{r.poles[i]=:.3e}, {r.residues[i]=:.3e}")
    r.poles[i]=-3.816e+00+0.000e+00j, r.residues[i]=3.658e-02+0.000e+00j
    r.poles[i]=-3.003e+00+0.000e+00j, r.residues[i]=-1.692e-01-0.000e+00j
    r.poles[i]=-2.000e+00+0.000e+00j, r.residues[i]=5.000e-01+0.000e+00j
    r.poles[i]=-1.000e+00+0.000e+00j, r.residues[i]=-1.000e+00+0.000e+00j
    r.poles[i]=5.858e-17+0.000e+00j, r.residues[i]=1.000e+00+0.000e+00j  # may vary
    r.poles[i]=4.775e+00-3.069e+00j, r.residues[i]=-8.113e-01-2.302e+00j
    r.poles[i]=4.775e+00+3.069e+00j, r.residues[i]=-8.113e-01+2.302e+00j
    r.poles[i]=5.291e+00-9.737e-01j, r.residues[i]=8.733e-01+1.070e+01j
    r.poles[i]=5.291e+00+9.737e-01j, r.residues[i]=8.733e-01-1.070e+01j

    For the second example, we call `AAA` with a spiral of 1000 points wind 7.5 times
    around the origin in the complex plane.

    >>> z = np.exp(np.linspace(-0.5, 0.5 + 15j*np.pi, 1000))
    >>> r = AAA(lambda z: np.tan(np.pi*z/2), z)

    We see that AAA takes 12 steps to converge with the following errors:

    >>> r.errors.size
    12
    >>> r.errors
    array([2.49261500e+01+0.j, 4.28045609e+01+0.j, 1.71346935e+01+0.j,
        8.65055336e-02+0.j, 1.27106444e-02+0.j, 9.90889874e-04+0.j,
        5.86910543e-05+0.j, 1.28735561e-06+0.j, 3.57007532e-08+0.j,
        6.37019486e-10+0.j, 1.67119305e-11+0.j, 1.17128023e-13+0.j])

    We can also plot the computed poles:

    >>> fig, ax = plt.subplots()
    >>> ax.plot(z.real, z.imag, '.', markersize=2, label="Sample points")
    >>> ax.plot(r.poles.real, r.poles.imag, '.', markersize=5, label="Computed poles")
    >>> ax.set(xlim=[-3.5, 3.5], ylim=[-3.5, 3.5], aspect="equal")
    >>> ax.legend()
    >>> plt.show()
    """
    def __init__(self, f, z, *, rtol=1e-13, max_terms=100):
        # input validation
        z = np.asarray(z).ravel()
        if callable(f):
            f = f(z)
        else:
            f = np.asarray(f).ravel()

            if f.size != z.size:
                raise ValueError("`f` and `z` must be the same size.")

        # Remove infinite or NaN function values and repeated entries
        to_keep = (np.isfinite(f)) & (~np.isnan(f))
        f = f[to_keep]
        z = z[to_keep]
        z, uni = np.unique(z, return_index=True)
        f = f[uni]

        # Initialization for AAA iteration
        M = np.size(z)
        atol = rtol * np.linalg.norm(f, ord=np.inf)
        mask = np.ones(M, dtype=np.bool_)
        zj = np.empty(max_terms, dtype=np.complex128)
        fj = np.empty(max_terms, dtype=np.complex128)
        # Cauchy matrix
        C = np.empty((M, max_terms), dtype=np.complex128)
        # Loewner matrix
        A = np.empty((M, max_terms), dtype=np.complex128)
        errors = np.empty(max_terms, dtype=np.complex128)
        R = np.repeat(np.mean(f), M)

        # AAA iteration
        for m in range(max_terms):
            # Introduce next support point
            # Select next support point
            jj = np.argmax(np.abs(f[mask] - R[mask]))
            # Update support points
            zj[m] = z[mask][jj]
            # Update data values
            fj[m] = f[mask][jj]
            # Next column of Cauchy matrix
            with np.errstate(divide="ignore", invalid="ignore"):
                C[:, m] = 1 / (z - z[mask][jj])
            # Update mask
            mask[np.flatnonzero(mask)[jj]] = False
            # Update Loewner matrix
            with np.errstate(invalid="ignore"):
                A[:, m] = (f - fj[m]) * C[:, m]

            # Compute weights
            rows = mask.sum()
            if rows >= m + 1:
                # The usual tall-skinny case
                # Reduced SVD, use ?gesvd rather than ?gesdd to avoid convergence issues
                _, s, V = scipy.linalg.svd(
                    A[mask, : m + 1], full_matrices=False, check_finite=False,
                    lapack_driver="gesvd"
                )
                # Treat case of multiple min singular values
                mm = s == np.min(s)
                # Aim for non-sparse weight vector
                wj = V.conj()[mm, :].sum(axis=0) / np.sqrt(mm.sum())
            else:
                # Fewer rows than columns
                V = scipy.linalg.null_space(A[mask, : m + 1])
                nm = V.shape[-1]
                # Aim for non-sparse wt vector
                wj = V.sum(axis=-1) / np.sqrt(nm)

            # Compute rational approximant
            # Omit columns with `wj == 0`
            i0 = wj != 0
            with np.errstate(invalid="ignore"):
                # Numerator
                N = C[:, : m + 1][:, i0] @ (wj[i0] * fj[: m + 1][i0])
                # Denominator
                D = C[:, : m + 1][:, i0] @ wj[i0]
            # Interpolate at support points with `wj !=0`
            D_inf = np.isinf(D) | np.isnan(D)
            D[D_inf] = 1
            N[D_inf] = f[D_inf]
            R = N / D

            # Check if converged
            max_error = np.linalg.norm(f - R, ord=np.inf)
            errors[m] = max_error
            if max_error <= atol:
                break

        if m == max_terms - 1:
            warnings.warn(f"AAA failed to converge within {max_terms} iterations.",
                          RuntimeWarning, stacklevel=2)

        # Trim off unused array allocation
        zj = zj[: m + 1]
        fj = fj[: m + 1]

        # Remove support points with zero weight
        i_non_zero = wj != 0
        self.support_points = zj[i_non_zero]
        self.values = fj[i_non_zero]
        self.weights = wj[i_non_zero]
        self.errors = errors[: m + 1]

    def __call__(self, z):
        """Evaluate the rational approximation at given values.

        Parameters
        ----------
        z : array_like
            Input values.
        """
        # evaluate rational function in barycentric form.
        z = np.asarray(z)
        zv = np.ravel(z)

        # Cauchy matrix
        with np.errstate(invalid="ignore", divide="ignore"):
            CC = 1 / np.subtract.outer(zv, self.support_points)
        # Vector of values
        with np.errstate(invalid="ignore"):
            r = CC @ (self.weights * self.values) / (CC @ self.weights)

        # Deal with input inf: `r(inf) = lim r(z) = sum(w*f) / sum(w)`
        with np.errstate(divide="ignore"):
            r[np.isinf(zv)] = np.sum(self.weights * self.values) / np.sum(self.weights)

        # Deal with NaN
        ii = np.flatnonzero(np.isnan(r))
        for jj in ii:
            if np.isnan(zv[jj]) or not np.any(zv[jj] == self.support_points):
                # r(NaN) = NaN is fine.
                # The second case may happen if `r(zv[ii]) = 0/0` at some point.
                pass
            else:
                # Clean up values `NaN = inf/inf` at support points.
                # Find the corresponding node and set entry to correct value:
                r[jj] = self.values[zv[jj] == self.support_points].squeeze()

        return np.reshape(r, z.shape)

    @cached_property
    def poles(self):
        # Compute poles via generalized eigenvalue problem
        m = self.weights.size
        B = np.eye(m + 1)
        B[0, 0] = 0

        E = np.zeros_like(B, dtype=np.complex128)
        E[0, 1:] = self.weights
        E[1:, 0] = 1
        np.fill_diagonal(E[1:, 1:], self.support_points)

        pol = scipy.linalg.eigvals(E, B)
        return pol[np.isfinite(pol)]

    @cached_property
    def residues(self):
        # Compute residues via formula for res of quotient of analytic functions
        with np.errstate(invalid="ignore", divide="ignore"):
            N = (1 / (self.poles[:, np.newaxis] - self.support_points)) @ (
                self.values * self.weights
            )
            Ddiff = (
                -((1 / np.subtract.outer(self.poles, self.support_points)) ** 2)
                @ self.weights
            )
            return N / Ddiff

    @cached_property
    def roots(self):
        # Compute zeros via generalized eigenvalue problem
        m = self.weights.size
        B = np.eye(m + 1)
        B[0, 0] = 0
        E = np.zeros_like(B, dtype=np.complex128)
        E[0, 1:] = self.weights * self.values
        E[1:, 0] = 1
        np.fill_diagonal(E[1:, 1:], self.support_points)

        zer = scipy.linalg.eigvals(E, B)
        return zer[np.isfinite(zer)]
