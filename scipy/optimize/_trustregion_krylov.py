from ._trustregion import (_minimize_trust_region)
from ._trlib import (TRLIBQuadraticSubproblem)

__all__ = ['_minimize_trust_krylov']

def _minimize_trust_krylov(fun, x0, args=(), jac=None, hess=None,
        hessp=None, inexact=True, **trust_region_options):
    """
    Minimization of a scalar function of one or more variables using
    a nearly exact trust-region algorithm that only requires matrix
    vector products with the hessian matrix.

    Options
    -------
    inexact : bool, optional
        Accuracy to solve subproblems. If True requires less nonlinear
        iterations, but more vector products.

    .. versionadded:: 1.0.0
    """

    if jac is None:
        raise ValueError('Jacobian is required for trust region ',
                         'exact minimization.')
    if hess is None and hessp is None:
        raise ValueError('Either the Hessian or the Hessian-vector product '
                         'is required for Krylov trust-region minimization')
    if inexact:
        return _minimize_trust_region(fun, x0, args=args, jac=jac,
                hess=hess, hessp=hessp,
                subproblem=TRLIBQuadraticSubproblem(
                    tol_rel_i=-2.0, tol_rel_b=-3.0),
                **trust_region_options)
    else:
        return _minimize_trust_region(fun, x0, args=args, jac=jac,
                hess=hess, hessp=hessp,
                subproblem=TRLIBQuadraticSubproblem(
                    tol_rel_i=1e-8, tol_rel_b=1e-6),
                **trust_region_options)
