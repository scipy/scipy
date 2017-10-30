"""quasi-Newton Hessian approximations"""

from __future__ import division, print_function, absolute_import
import numpy as np
from numpy.linalg import norm
from scipy.linalg import get_blas_funcs

class QuasiNewtonApprox(object):
    """Quasi-Newton approximation.

    Parameters
    ----------
    approx_type : {'hess', 'inv_hess'}, optional
        Define what is being approximated. Set it
        to 'hess' to store and update the Hessian
        matrix and to 'inv_hess' to store and update
        the inverse Hessian matrix.
    """
    syr = get_blas_funcs('syr', dtype='d')  # Symetric rank 1 update
    syr2 = get_blas_funcs('syr2', dtype='d')  # Symetric rank 2 update
    symv = get_blas_funcs('symv', dtype='d')  # Symetric matrix-vector product

    def __init__(self, approx_type='hess'):
        if approx_type not in ('hess', 'inv_hess'):
            ValueError("Unknown approx_type.")
        self.approx_type = approx_type

    def instanciate_matrix(self, s, y):
        """Instanciate unscaled matrix."""
        n = len(y)
        # Create matrix
        if self.approx_type == 'hess':
            self.B = np.eye(n, dtype=float)
        else:
            self.H = np.eye(n, dtype=float)

    def scale_matrix(self, s, y):
        """Scale matrix.

        Heuristic to scale matrix at first iteration described
        in Nocedal and Wright "Numerical Optimization"
        p.143 formula (6.20)
        """
        y_norm2 = np.dot(y, y)
        ys = np.dot(y, s)
        if ys == 0.0:
            scale = 1.0
        else:
            if self.approx_type == 'hess':
                scale = y_norm2 / ys
            else:
                scale = ys / y_norm2

        if self.approx_type == 'hess':
            self.B *= scale
        else:
            self.H *= scale

    def dot(self, x):
        """Matrix-vector multiplication.

        Parameters
        ----------
        x : array_like
            1-d array representing a vector.

        Returns
        -------
        Hx : array
            1-d  represents the result of multiplying the approximation matrix
            by vector x.
        """
        if self.approx_type == 'hess':
            return self.symv(1, self.B, x)
        else:
            return self.symv(1, self.H, x)

    def get_matrix(self):
        """Return current approximation matrix."""
        if self.approx_type == 'hess':
            M = self.B
        else:
            M = self.H
        M_triu = np.triu(M)
        M_diag = np.diag(M)
        return M_triu + M_triu.T - np.diag(M_diag)


class BFGS(QuasiNewtonApprox):
    """Broyden-Fletcher-Goldfarb-Shanno (BFGS) Hessian matrix approximation.
    Parameters
    ----------
    approx_type : {'hess', 'inv_hess'}, optional
        Define what is being approximated. Set it
        to 'hess' to store and update the Hessian
        matrix and to 'inv_hess' to store and update
        the inverse Hessian matrix.
    exception_strategy : {'skip_update', 'damped_bfgs'}, optional
        Define how to proceed when the curvature condition is violated.
        Set it to 'skip_update' to just skip the update. Or, alternatively,
        set it to 'damped_bfgs' to interpolate between the actual BFGS
        result and the unmodified matrix. Both methods are explained
        in [1]_, p.536-537.
    min_curvature : float
        Define the minimum allowed value of ``delta_grad.T delta_x``
        allowed to go unnafected by the exception strategy. By
        default is equal to ``1e-2`` when ``exception_strategy = 'skip_update'``
        and equal to ``0.2`` when ``exception_strategy = 'damped_bfgs'``.

    References
    ----------
    .. [1] Nocedal, Jorge, and Stephen J. Wright. "Numerical optimization"
           Second Edition (2006).
    """

    def __init__(self, approx_type='hess', exception_strategy='skip_update',
                 min_curvature=None):
        if exception_strategy == 'skip_update':
            self.min_curvature = min_curvature if min_curvature is not None else 1e-2
        elif exception_strategy == 'damped_bfgs':
            self.min_curvature = min_curvature if min_curvature is not None else 0.2
        else:
            ValueError("Unknown approx_type.")
        super(BFGS, self).__init__(approx_type)
        self.exception_strategy = exception_strategy

    def _update_inverse_hessian(self, ys, Hy, yHy, s):
        """Update inverse Hessian matrix.

        BFGS update using the formula:

            H <- H + ((H*y).T*y + s.T*y)/(s.T*y)^2 * (s*s.T) - 1/(s.T*y) * ((H*y)*s.T + s*(H*y).T)

        where ``s = delta_x`` and ``y = delta_grad``. This formula is
        equivalent to (6.17) in [1]_ written in a more efficient way
        for implementation.

        References
        ----------
        .. [1] Nocedal, Jorge, and Stephen J. Wright. "Numerical optimization"
               Second Edition (2006).
        """
        self.H = self.syr2(-1.0 / ys, s, Hy, a=self.H)
        self.H = self.syr((ys+yHy)/ys**2, s, a=self.H)
        return

    def _update_hessian(self, ys, Bs, sBs, y):
        """Update Hessian matrix.

        BFGS update using the formula:

            B <- B - (B*s)*(B*s).T/s.T*(B*s) + y*y^T/s.T*y

        where ``s`` is short for ``delta_x`` and ``y`` is short
        for ``delta_grad``. Formula (6.19) in [1]_.

        References
        ----------
        .. [1] Nocedal, Jorge, and Stephen J. Wright. "Numerical optimization"
               Second Edition (2006).
        """
        self.B = self.syr(1.0 / ys, y, a=self.B)
        self.B = self.syr(-1.0 / sBs, Bs, a=self.B)
        return

    def update(self, delta_x, delta_grad):
        """Update approximation matrix.

        Parameters
        ----------
        delta_x : ndarray
            The difference between two points the gradient
            function have been evaluated: ``delta_x = x2 - x1``.
        delta_grad : ndarray
            The difference between the gradient evaluated
            in two points: ``delta_grad = grad(x2) - grad(x1)``.
        """
        # Auxiliar variables w and z
        if self.approx_type == 'hess':
            w = delta_x
            z = delta_grad
        else:
            w = delta_grad
            z = delta_x
        # Do some common operations
        wz = np.dot(w, z)
        Mw = self.dot(w)
        wMw = Mw.dot(w)
        if wMw == 0.0:
            return
        # Check if curvature condition is violated
        if wz < self.min_curvature * wMw:
            # If the option 'skip_update' is set
            # we just skip the update when the condion
            # is violated.
            if self.exception_strategy == 'skip_update':
                return
            # If the option 'damped_bfgs' is set we
            # interpolate between the actual BFGS
            # result and the unmodified matrix.
            elif self.exception_strategy == 'damped_bfgs':
                update_factor = (1-self.min_curvature) / (1 - wz/wMw)
                z = update_factor*z + (1-update_factor)*Mw
                wz = np.dot(w, z)
                Mw = self.dot(w)
                wMw = Mw.dot(w)
        # Update matrix
        if self.approx_type == 'hess':
            return self._update_hessian(wz, Mw, wMw, z)
        else:
            return self._update_inverse_hessian(wz, Mw, wMw, z)


class SR1(QuasiNewtonApprox):
    """Symmetric-rank-1 Hessian matrix approximation.

    Parameters
    ----------
    approx_type : {'hess', 'inv_hess'}, optional
        Define what is being approximated. Set it
        to 'hess' to store and update the Hessian
        matrix and to 'inv_hess' to store and update
        the inverse Hessian matrix.
    min_curvature : float
        Define the minimum allowed value of the denominator
        of the update. When the condition is violated we skip
        the update. By default uses ``1e-8``.

    References
    ----------
    .. [1] Nocedal, Jorge, and Stephen J. Wright. "Numerical optimization"
           Second Edition (2006).
    """

    def __init__(self, approx_type='hess', min_curvature=1e-8):
        self.min_curvature = min_curvature
        super(SR1, self).__init__(approx_type)

    def update(self, delta_x, delta_grad):
        """Update approximation matrix.

        Parameters
        ----------
        delta_x : ndarray
            The difference between two points the gradient
            function have been evaluated: ``delta_x = x2 - x1``.
        delta_grad : ndarray
            The difference between the gradient evaluated
            in two points: ``delta_grad = grad(x2) - grad(x1)``.
        """
        # Auxiliar variables w and z
        if self.approx_type == 'hess':
            w = delta_x
            z = delta_grad
        else:
            w = delta_grad
            z = delta_x
        # Do some common operations
        Mw = self.dot(w)
        z_minus_Mw = z - Mw
        denominator = np.dot(w, z_minus_Mw)
        # If the denominator is too small
        # we just skip the update.
        if np.abs(denominator) < self.min_curvature*norm(w)*norm(z_minus_Mw) or denominator == 0.0:
            return
        # Update matrix
        if self.approx_type == 'hess':
            self.B = self.syr(1/denominator, z_minus_Mw, a=self.B)
        else:
            self.H = self.syr(1/denominator, z_minus_Mw, a=self.H)
        return
    
