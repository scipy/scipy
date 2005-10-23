"""This module defines encapsulation classes for all iterative solvers
implemented in the itsolvers module.

All classes provide a method "solve(b, x)" for approximatively compute
x = inv(A)*b.

The classes are intended to replace superlu.superlu_context in cases
where appropriate."""

import itsolvers

class ItSolver:
    "abstract base class for iteravtive solver classes"
    def __init__(self, A, tol, maxit, K, debug):
        self.A = A
        self.tol = tol
        self.maxit = maxit
        self.K = K
        self.nofCalled = 0
        self.totalIterations = 0
        self.lastIterations = 0
        self.lastInfo = 0
        self.debug = debug
        
    def solve(self, b, x):
        "solve A*x = b iteratively with zero initial guess"
        x[:] = 0
        info, iter, relres = self.itsolver(self.A, b, x, self.tol, self.maxit, self.K)
        self.nofCalled += 1
        self.totalIterations += iter
        self.lastIterations = iter
        self.lastInfo = info
        if self.debug:
            print 'iterative solver returned:', info, iter, relres
        if info < 0:
            raise 'error', ('iterative solver %s returned error code %d' % (self.__class__.__name__, info), info, iter, relres)

    def __str__(self):
        s = '<%s.%s instance>\n\n' % (self.__class__.__module__, self.__class__.__name__)
        for name in ['nofCalled', 'totalIterations', 'lastIterations', 'lastInfo']:
            s += '   %s: %s\n' % (name, getattr(self, name))
        s += '\n'
        return s

    def __repr__(self):
        return self.__str__()

class Pcg(ItSolver):
    def __init__(self, A, tol, maxit, K=None, debug=0):
        ItSolver.__init__(self, A, tol, maxit, K, debug)
        self.itsolver = itsolvers.pcg

class Minres(ItSolver):
    def __init__(self, A, tol, maxit, K=None, debug=0):
        ItSolver.__init__(self, A, tol, maxit, K, debug)
        self.itsolver = itsolvers.minres
        
class Qmrs(ItSolver):
    def __init__(self, A, tol, maxit, K=None, debug=0):
        ItSolver.__init__(self, A, tol, maxit, K, debug)
        self.itsolver = itsolvers.qmrs

class Cgs(ItSolver):
    """wrapper class for the itsolvers.cgs iterative solver

Cgs(A, tol, maxit, K=None) constructs the Cgs object.

methods:

solve(b, x) solves the linear system A*x = b with a zero initial guess
"""
    def __init__(self, A, tol, maxit, K=None, debug=0):
        ItSolver.__init__(self, A, tol, maxit, K, debug)
        self.itsolver = itsolvers.cgs
                        
if __name__ == '__main__':
    import math
    import Numeric
    import precon, poisson

    A = poisson.poisson2d_sym(100)
    n = A.shape[0];
    b = Numeric.ones(n, 'd'); b = b / math.sqrt(Numeric.dot(b, b))
    x = Numeric.zeros(n, 'd')

    def resid(A, b, x):
        r = x.copy()
        A.matvec(x, r)
        r = b - r
        return math.sqrt(Numeric.dot(r, r))
    
    solver = Pcg(A, 1e-10, 300)
    solver.solve(b, x)
    print resid(A, b, x), solver.nofCalled, solver.totalIterations

    solver.K = precon.ssor(A.to_sss())
    solver.solve(b, x)
    print resid(A, b, x), solver.nofCalled, solver.totalIterations
    
