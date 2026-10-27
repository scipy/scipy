import os
import warnings

import numpy as np
from numpy.linalg import LinAlgError
from scipy._lib._ccallback_c import get_capsule_signature

from scipy.linalg import LinAlgWarning
from scipy.linalg.lapack import get_lapack_funcs
from ._rbfinterp_common import _monomial_powers_impl

from . import _rbfinterp_pythran as _pythran_mod
from ._rbfinterp_pythran import (
    _build_system_with_kernel as _pythran_build_system_with_kernel,
    _build_evaluation_coefficients_with_kernel as
        _pythran_build_evaluation_coefficients_with_kernel,
    _polynomial_matrix as _pythran_polynomial_matrix
)


dsysv, dsytrs = get_lapack_funcs(
    ('sysv', 'sytrs'), dtype=np.float64, ilp64="preferred"
)


def _get_kernel_capsule(kernel):
    """Return a float64(float64) capsule for given name.

    Parameters
    ----------
    kernel : str or LowLevelCallable
        Either the name of a built-in RBF (str) or a ``LowLevelCallable``
        wrapping a user-compiled ``double (*)(double)`` kernel.

    Returns
    -------
    capsule : PyCapsule
        Opaque capsule suitable for passing to
        ``_build_system_with_kernel`` /
        ``_build_evaluation_coefficients_with_kernel``.

    Notes
    -----
    For the ``str`` path, the compiled capsule is retrieved as a module
    attribute of ``_bfinterp_pythran``

    For the ``LowLevelCallable`` path, element 0 of the LLC tuple is the
    normalised PyCapsule.  Pythran ignores the capsule name and trusts the
    ``float64(float64)`` annotation in the export line, so no ``signature=`` coercion
    is needed on the LLC side.
    """
    if isinstance(kernel, str):
        # Built in kernel
        return getattr(_pythran_mod, kernel)
    else:
        # LowLevelCallable is a tuple subclass; element 0 is the capsule
        capsule =  tuple.__getitem__(kernel, 0)
        sig = get_capsule_signature(capsule)
        if sig != "double (double)":
            raise ValueError(
                f"LowLevelCallable kernel must have signature \"double (double)\", got"
                f" \"{sig}\"."
                f"Construct with: "
                f"LowLevelCallable(fn, signature=\"double (double)\")"
            )
        return capsule


# trampolines for pythran-compiled functions to drop the `xp` argument
def _build_evaluation_coefficients(
    x, y, kernel, epsilon, powers, shift, scale, xp
):
    capsule = _get_kernel_capsule(kernel)
    return _pythran_build_evaluation_coefficients_with_kernel(
        x, y, capsule, epsilon, powers, shift, scale
    )


def polynomial_matrix(x, powers, xp):
    return _pythran_polynomial_matrix(x, powers)


def _monomial_powers(ndim, degree, xp, device=None):
    out = _monomial_powers_impl(ndim, degree)
    out = np.asarray(out, dtype=np.int64)
    if len(out) == 0:
        out = out.reshape(0, ndim)
    return out


def _build_system(y, d, smoothing, kernel, epsilon, powers, xp):
    capsule = _get_kernel_capsule(kernel)
    return _pythran_build_system_with_kernel(
        y, d, smoothing, capsule, epsilon, powers
    )


def _scaled_norm(lhs, nobs):
    """Return the 1-norm of ``D @ lhs @ D`` and the scale ``kmax`` used in ``D``.

    The condition number is measured with the polynomial blocks scaled by the
    largest kernel entry, so that it does not depend on how large the kernel
    values are compared to the monomials (gh-14150). The system itself is solved
    unscaled.
    """
    kabs = np.abs(lhs[:nobs, :nobs])
    kmax = kabs.max()
    if kmax == 0.0:
        kmax = 1.0
    colsum = np.empty(lhs.shape[1])
    colsum[:nobs] = kabs.sum(axis=0) + kmax*np.abs(lhs[nobs:, :nobs]).sum(axis=0)
    colsum[nobs:] = kmax*np.abs(lhs[:nobs, nobs:]).sum(axis=0)
    return colsum.max(), kmax


def _scaled_rcond(lu, ipiv, anorm, nobs, kmax):
    """Estimate the reciprocal 1-norm condition number of ``D @ lhs @ D``.

    ``D`` is diagonal with ones for the RBF coefficients and `kmax` for the
    monomial coefficients, and `anorm` is the 1-norm of ``D @ lhs @ D``. The
    norm of the inverse is estimated from the Bunch-Kaufman factorization of
    ``lhs`` with a simplified version of the Hager/Higham method used by
    LAPACK's ``?lacn2``.
    """
    n = lu.shape[0]
    dinv = np.ones(n)
    dinv[nobs:] = 1.0/kmax

    def solve(v):
        # The matrix is symmetric, so this is also the transpose solve.
        x, _ = dsytrs(lu, ipiv, dinv*v)
        return dinv*x

    x = np.full(n, 1.0/n)
    est = 0.0
    for _ in range(5):
        y = solve(x)
        new_est = np.abs(y).sum()
        if new_est <= est:
            break
        est = new_est
        z = solve(np.where(y >= 0, 1.0, -1.0))
        j = np.argmax(np.abs(z))
        if abs(z[j]) <= z @ x:
            break
        x = np.zeros(n)
        x[j] = 1.0

    # Higham's alternative estimate, which catches cases where the above fails
    alt = (-1.0)**np.arange(n) * (1.0 + np.arange(n)/max(n - 1, 1))
    est = max(est, 2*np.abs(solve(alt)).sum()/(3*n))
    return 1.0/(anorm*est)


def _build_and_solve_system(y, d, smoothing, kernel, epsilon, powers, xp,
                            check_cond=True):
    """Build and solve the RBF interpolation system of equations.

    Parameters
    ----------
    y : (P, N) float ndarray
        Data point coordinates.
    d : (P, S) float ndarray
        Data values at `y`.
    smoothing : (P,) float ndarray
        Smoothing parameter for each data point.
    kernel : str
        Name of the RBF.
    epsilon : float
        Shape parameter.
    powers : (R, N) int ndarray
        The exponents for each monomial in the polynomial.
    check_cond : bool, optional
        Whether to warn if the system is ill-conditioned. The condition number
        estimate costs about as much as a few solves, which is negligible for
        one large system but not for the many small ones solved when
        ``neighbors`` is given.

    Returns
    -------
    coeffs : (P + R, S) float ndarray
        Coefficients for each RBF and monomial.
    shift : (N,) float ndarray
        Domain shift used to create the polynomial matrix.
    scale : (N,) float ndarray
        Domain scaling used to create the polynomial matrix.

    """
    lhs, rhs, shift, scale = _build_system(
        y, d, smoothing, kernel, epsilon, powers, xp
        )
    nobs = y.shape[0]
    if check_cond:
        # `dsysv` overwrites `lhs`, so take its norm first.
        anorm, kmax = _scaled_norm(lhs, nobs)

    lu, ipiv, coeffs, info = dsysv(lhs, rhs, overwrite_a=True, overwrite_b=True)
    if info < 0:
        raise ValueError(f"The {-info}-th argument had an illegal value.")
    elif info > 0:
        msg = "Singular matrix."
        nmonos = powers.shape[0]
        if nmonos > 0:
            pmat = polynomial_matrix((y - shift)/scale, powers, xp)
            rank = np.linalg.matrix_rank(pmat)
            if rank < nmonos:
                msg = (
                    "Singular matrix. The matrix of monomials evaluated at "
                    "the data point coordinates does not have full column "
                    f"rank ({rank}/{nmonos})."
                    )

        raise LinAlgError(msg)

    if check_cond:
        rcond = _scaled_rcond(lu, ipiv, anorm, nobs, kmax)
        if rcond < np.finfo(np.float64).eps:
            # show the warning at the user's call site
            _warn_skips = (os.path.dirname(os.path.dirname(__file__)),)
            warnings.warn(
                "The RBF interpolation matrix is ill-conditioned "
                f"(rcond={rcond:.2e}), the result may not be accurate.",
                LinAlgWarning, skip_file_prefixes=_warn_skips,
            )

    return shift, scale, coeffs

def compute_interpolation(x, y, kernel, epsilon, powers, shift, scale, coeffs, xp):
    vec = _build_evaluation_coefficients(
        x, y, kernel, epsilon, powers, shift, scale, xp
    )
    return vec @ coeffs
