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

import warnings
import operator

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
    points : 1D array_like, shape (n,)
        1-D array containing values of the independent variable. Values may be real or
        complex but must be finite.
    values : 1D array_like, shape (n,)
        Function values ``f(z)`` at `points`. Infinite and NaN values of `values` and
        corresponding values of `points` will be discarded.
    rtol : float, optional
        Relative tolerance, defaults to ``eps**0.75``. If a small subset of the entries
        in `values` are much larger than the rest the default tolerance may be too
        loose.
    max_terms : int, optional
        Maximum number of terms in the barycentric representation, defaults to ``100``.
        Must be greater than or equal to one.

    Attributes
    ----------
    support_points : array
        Support points of the approximation. These are a subset of the provided
        `points` at which the approximation strictly interpolates the provided `values`.
        See notes for more details.
    support_values : array
        Value of the approximation at the `support_points`.
    weights : array
        Weights of the barycentric approximation.
    errors : array
        Error :math:`|f(z) - r(z)|_\infty` over `points` in the successive iterations
        of AAA.

    Warns
    -----
    RuntimeWarning
        If `rtol` is not achieved in `max_terms` iterations.

    See Also
    --------
    pade : Padé approximation

    Notes
    -----
    At iteration :math:`m` (at which point there are :math:`m` terms in the both the
    numerator and denominator of the approximation), the
    rational approximation in the AAA algorithm takes the barycentric form

    .. math::

        r(z) = n(z)/d(z) =
        \frac{\sum_{j=1}^m\ w_j f_j / (z - z_j)}{\sum_{j=1}^m w_j / (z - z_j)},

    where :math:`z_1,\dots,z_m` are real or complex support points selected from
    `points`, :math:`f_1,\dots,f_m` are the corresponding real or complex data values
    from `values`, and :math:`w_1,\dots,w_m` are real or complex weights.

    Each iteration of the algorithm has two parts: the greedy selection the next support
    point and the computation of the weights. The first part of each iteration is to
    select the next support point to be added :math:`z_{m+1}` from the remaining
    unselected `points`, such that the nonlinear residual
    :math:`|f(z_{m+1}) - n(z_{m+1})/d(z_{m+1})|` is maximised. The algorithm terminates
    when this maximum is less than ``rtol * np.linalg.norm(f, ord=np.inf)``. This means
    the interpolation property is only satisfied up to a tolerance, except at the
    support points where approximation exactly interpolates the supplied data.

    In the second part of each iteration, the weights :math:`w_j` are selected to solve
    the least-squares problem

    .. math::

        \text{minimise}_{w_j}|fd - n| \quad \text{subject to} \quad
        \sum_{j=1}^{m+1} w_j = 1,

    over the unselected elements of `points`.

    References
    ----------
    .. [1] Y. Nakatsukasa, O. Sete, and L. N. Trefethen, "The AAA algorithm for
            rational approximation", SIAM J. Sci. Comp. 40 (2018), A1494-A1522.
            :doi:`10.1137/16M110612`

    Examples
    --------

    Here we reproduce a number of the numerical examples from [1]_ as a demonstration
    of the functionality offered by this method.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.interpolate import AAA

    For the first example we approximate the gamma function on ``[-3.5, 4.5]`` by
    extrapolating from 100 samples in ``[-1.5, 1.5]``.

    >>> from scipy.special import gamma
    >>> sample_points = np.linspace(-1.5, 1.5, num=100)
    >>> r = AAA(sample_points, gamma(sample_points))
    >>> z = np.linspace(-3.5, 4.5, num=1000)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(z, gamma(z), label="Gamma")
    >>> ax.plot(sample_points, gamma(sample_points), label="Sample points")
    >>> ax.plot(z, r(z).real, '--', label="AAA approximation")
    >>> ax.set(xlabel="z", ylabel="r(z)", ylim=[-8, 8], xlim=[-3.5, 4.5])
    >>> ax.legend()
    >>> plt.show()

    We can also view the poles of the rational approximation and their residues:

    >>> order = np.argsort(r.poles())
    >>> r.poles()[order]
    array([-3.81591039e+00+0.j        , -3.00269049e+00+0.j        ,
           -1.99999988e+00+0.j        , -1.00000000e+00+0.j        ,
            5.85842812e-17+0.j        ,  4.77485458e+00-3.06919376j,
            4.77485458e+00+3.06919376j,  5.29095868e+00-0.97373072j,
            5.29095868e+00+0.97373072j])
    >>> r.residues()[order]
    array([ 0.03658074 +0.j        , -0.16915426 -0.j        ,
            0.49999915 +0.j        , -1.         +0.j        ,
            1.         +0.j        , -0.81132013 -2.30193429j,
           -0.81132013 +2.30193429j,  0.87326839+10.70148546j,
            0.87326839-10.70148546j])

    For the second example, we call `AAA` with a spiral of 1000 points that wind 7.5
    times around the origin in the complex plane.

    >>> z = np.exp(np.linspace(-0.5, 0.5 + 15j*np.pi, 1000))
    >>> r = AAA(z, np.tan(np.pi*z/2), rtol=1e-13)

    We see that AAA takes 12 steps to converge with the following errors:

    >>> r.errors.size
    12
    >>> r.errors
    array([2.49261500e+01, 4.28045609e+01, 1.71346935e+01, 8.65055336e-02,
           1.27106444e-02, 9.90889874e-04, 5.86910543e-05, 1.28735561e-06,
           3.57007424e-08, 6.37007837e-10, 1.67103357e-11, 1.17112299e-13])

    We can also plot the computed poles:

    >>> fig, ax = plt.subplots()
    >>> ax.plot(z.real, z.imag, '.', markersize=2, label="Sample points")
    >>> ax.plot(r.poles().real, r.poles().imag, '.', markersize=5,
    ...         label="Computed poles")
    >>> ax.set(xlim=[-3.5, 3.5], ylim=[-3.5, 3.5], aspect="equal")
    >>> ax.legend()
    >>> plt.show()
    """
    def __init__(self, points, values, *, rtol=None, max_terms=100):
        # input validation
        z = np.asarray(points)
        f = np.asarray(values)

        if z.ndim != 1 or f.ndim != 1:
            raise ValueError("`points` and `values` must be 1-D.")

        if f.size != z.size:
            raise ValueError("`points` and `values` must be the same size.")

        if not np.all(np.isfinite(z)):
            raise ValueError("`points` must be finite.")

        max_terms = operator.index(max_terms)
        if max_terms < 1:
            raise ValueError("`max_terms` must be an integer value greater than or "
                             "equal to one.")

        # Remove infinite or NaN function values and repeated entries
        to_keep = (np.isfinite(f)) & (~np.isnan(f))
        f = f[to_keep]
        z = z[to_keep]
        z, uni = np.unique(z, return_index=True)
        f = f[uni]

        self._compute_weights(f, z, rtol, max_terms)

        # only compute once
        self._poles = None
        self._residues = None
        self._roots = None

    def _compute_weights(self, f, z, rtol, max_terms):
        # Initialization for AAA iteration
        M = np.size(z)
        mask = np.ones(M, dtype=np.bool_)
        dtype = np.result_type(z, f, 1.0)
        rtol = np.finfo(dtype).eps**0.75 if rtol is None else rtol
        atol = rtol * np.linalg.norm(f, ord=np.inf)
        zj = np.empty(max_terms, dtype=dtype)
        fj = np.empty(max_terms, dtype=dtype)
        # Cauchy matrix
        C = np.empty((M, max_terms), dtype=dtype)
        # Loewner matrix
        A = np.empty((M, max_terms), dtype=dtype)
        errors = np.empty(max_terms, dtype=A.real.dtype)
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
            # Ignore errors as we manually interpolate at support points
            with np.errstate(divide="ignore", invalid="ignore"):
                C[:, m] = 1 / (z - z[mask][jj])
            # Update mask
            mask[np.flatnonzero(mask)[jj]] = False
            # Update Loewner matrix
            # Ignore errors as inf values will be masked out in SVD call
            with np.errstate(invalid="ignore"):
                A[:, m] = (f - fj[m]) * C[:, m]

            # Compute weights
            rows = mask.sum()
            if rows >= m + 1:
                # The usual tall-skinny case
                _, s, V = scipy.linalg.svd(
                    A[mask, : m + 1], full_matrices=False, check_finite=False,
                )
                # Treat case of multiple min singular values
                mm = s == np.min(s)
                # Aim for non-sparse weight vector
                wj = (V.conj()[mm, :].sum(axis=0) / np.sqrt(mm.sum())).astype(dtype)
            else:
                # Fewer rows than columns
                V = scipy.linalg.null_space(A[mask, : m + 1], check_finite=False)
                nm = V.shape[-1]
                # Aim for non-sparse wt vector
                wj = V.sum(axis=-1) / np.sqrt(nm)

            # Compute rational approximant
            # Omit columns with `wj == 0`
            i0 = wj != 0
            # Ignore errors as we manually interpolate at support points
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
        self.support_values = fj[i_non_zero]
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
        # Ignore errors due to inf/inf at support points, these will be fixed later
        with np.errstate(invalid="ignore", divide="ignore"):
            CC = 1 / np.subtract.outer(zv, self.support_points)
            # Vector of values
            r = CC @ (self.weights * self.support_values) / (CC @ self.weights)

        # Deal with input inf: `r(inf) = lim r(z) = sum(w*f) / sum(w)`
        if np.any(np.isinf(zv)):
            r[np.isinf(zv)] = (np.sum(self.weights * self.support_values)
                               / np.sum(self.weights))

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
                r[jj] = self.support_values[zv[jj] == self.support_points].squeeze()

        return np.reshape(r, z.shape)

    def poles(self):
        """Compute the poles of the rational approximation.

        Returns
        -------
        poles : array
            Poles of the AAA approximation, repeated according to their multiplicity
            but not in any specific order.
        """
        if self._poles is None:
            # Compute poles via generalized eigenvalue problem
            m = self.weights.size
            B = np.eye(m + 1, dtype=self.weights.dtype)
            B[0, 0] = 0

            E = np.zeros_like(B, dtype=np.result_type(self.weights,
                                                      self.support_points))
            E[0, 1:] = self.weights
            E[1:, 0] = 1
            np.fill_diagonal(E[1:, 1:], self.support_points)

            pol = scipy.linalg.eigvals(E, B)
            self._poles = pol[np.isfinite(pol)]
        return self._poles

    def residues(self):
        """Compute the residues of the poles of the approximation.

        Returns
        -------
        residues : array
            Residues associated with the `poles` of the approximation
        """
        if self._residues is None:
            # Compute residues via formula for res of quotient of analytic functions
            N = (1/(np.subtract.outer(self.poles(), self.support_points))) @ (
                self.support_values * self.weights
            )
            Ddiff = (
                -((1/np.subtract.outer(self.poles(), self.support_points))**2)
                @ self.weights
            )
            self._residues = N / Ddiff
        return self._residues

    def roots(self):
        """Compute the zeros of the rational approximation.

        Returns
        -------
        zeros : array
            Zeros of the AAA approximation, repeated according to their multiplicity
            but not in any specific order.
        """
        if self._roots is None:
            # Compute zeros via generalized eigenvalue problem
            m = self.weights.size
            B = np.eye(m + 1, dtype=self.weights.dtype)
            B[0, 0] = 0
            E = np.zeros_like(B, dtype=np.result_type(self.weights, self.support_values,
                                                      self.support_points))
            E[0, 1:] = self.weights * self.support_values
            E[1:, 0] = 1
            np.fill_diagonal(E[1:, 1:], self.support_points)

            zer = scipy.linalg.eigvals(E, B)
            self._roots = zer[np.isfinite(zer)]
        return self._roots
