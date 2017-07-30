from ._trlib import TRLIBQuadraticSubproblem

__all__ = ['TRLIBQuadraticSubproblem', 'get_trlib_quadratic_subproblem']


def get_trlib_quadratic_subproblem(tol_rel_i=-2.0, tol_rel_b=-3.0):
    def subproblem_factory(x, fun, jac, hess, hessp):
        return TRLIBQuadraticSubproblem(x, fun, jac, hess, hessp,
                                        tol_rel_i=tol_rel_i,
                                        tol_rel_b=tol_rel_b)
    return subproblem_factory
