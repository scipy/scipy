"""Dog-leg trust-region optimization."""
from __future__ import division, print_function, absolute_import

import numpy as np
import scipy.linalg
from ._trustregion import _help_solve_subproblem, _minimize_trust_region

__all__ = []

def _minimize_dogleg(
        fun,
        x0,
        args=(),
        jac=None,
        hess=None,
        **trust_region_options):
    """
    Minimization of scalar function of one or more variables using
    the dog-leg trust-region algorithm.

    This function is called by the `minimize` function.
    It is not supposed to be called directly.
    """
    if jac is None:
        raise ValueError('Jacobian is required for dogleg minimization')
    if hess is None:
        raise ValueError('Hessian is required for dogleg minimization')
    return scipy.optimize._trustregion._minimize_trust_region(
            fun, x0,
            args=args,
            jac=jac,
            hess=hess,
            solve_subproblem=_solve_subproblem_dogleg,
            **trust_region_options)


def _solve_subproblem_dogleg(m, trust_radius):
    """
    Minimize a function using the dog-leg trust-region algorithm.

    This algorithm requires function values and first and second derivatives.
    It also performs a costly Hessian decomposition for most iterations,
    and the Hessian is required to be positive definite.

    Parameters
    ----------
    m : LazyLocalQuadraticModel
        The quadratic model of the objective function.
    trust_radius : float
        We are allowed to wander only this far away from the origin.

    Returns
    -------
    p : ndarray
        The proposed step.
    hits_boundary : bool
        True if the proposed step is on the boundary of the trust region.

    Notes
    -----
    The Hessian is required to be positive definite.

    This function is called by the `_minimize_trust_region` function.
    It is not supposed to be called directly.

    References
    ----------
    .. [1] Jorge Nocedal and Stephen Wright,
           Numerical Optimization, second edition,
           Springer-Verlag, 2006, page 73.

    """

    # Compute the Newton point.
    # This is the optimum for the quadratic model function.
    # If it is inside the trust radius then return this point.
    p_best = m.newton_point()
    if scipy.linalg.norm(p_best) < trust_radius:
        hits_boundary = False
        return p_best, hits_boundary

    # Compute the Cauchy point.
    # This is the predicted optimum along the direction of steepest descent.
    p_u = m.cauchy_point()

    # If the Cauchy point is outside the trust region,
    # then return the point where the path intersects the boundary.
    p_u_norm = scipy.linalg.norm(p_u)
    if p_u_norm >= trust_radius:
        p_boundary = p_u * (trust_radius / p_u_norm)
        hits_boundary = True
        return p_boundary, hits_boundary

    # Compute the intersection of the trust region boundary
    # and the line segment connecting the Cauchy and Newton points.
    # This requires solving a quadratic equation.
    # ||p_u + t*(p_best - p_u)||**2 == trust_radius**2
    # Solve this for positive time t using the quadratic formula.
    ta, tb = _help_solve_subproblem(p_u, p_best - p_u, trust_radius)
    p_boundary = p_u + tb * (p_best - p_u)
    hits_boundary = True
    return p_boundary, hits_boundary

