from __future__ import division, print_function, absolute_import
"""
Unit test for constraint conversion
"""

import numpy as np
from numpy.testing import (assert_, assert_array_almost_equal,
                           assert_allclose, assert_equal, TestCase)
from scipy._lib._numpy_compat import suppress_warnings
from scipy.optimize import (NonlinearConstraint, LinearConstraint, Bounds,
                            OptimizeWarning, minimize)
from .test_minimize_constrained import (Maratos, HyperbolicIneq, Rosenbrock,
                                        IneqRosenbrock, EqIneqRosenbrock,
                                        BoundedRosenbrock, Elec)
from scipy._lib._numpy_compat import _assert_warns, suppress_warnings


class TestOldToNew(object):
    x0 = (2, 0)
    bnds = ((0, None), (0, None))
    method = "trust-constr"

    def test_constraint_dictionary_1(self):
        fun = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2
        cons = ({'type': 'ineq', 'fun': lambda x: x[0] - 2 * x[1] + 2},
                {'type': 'ineq', 'fun': lambda x: -x[0] - 2 * x[1] + 6},
                {'type': 'ineq', 'fun': lambda x: -x[0] + 2 * x[1] + 2})

        with suppress_warnings() as sup:
            sup.filter(UserWarning, "delta_grad == 0.0")
            res = minimize(fun, self.x0, method=self.method,
                           bounds=self.bnds, constraints=cons)
        assert_allclose(res.x, [1.4, 1.7], rtol=1e-4)
        assert_allclose(res.fun, 0.8, rtol=1e-4)

    def test_constraint_dictionary_2(self):
        fun = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2
        cons = {'type': 'eq',
                'fun': lambda x, p1, p2: p1*x[0] - p2*x[1],
                'args': (1, 1.1),
                'jac': lambda x, p1, p2: np.array([[p1, -p2]])}
        with suppress_warnings() as sup:
            sup.filter(UserWarning, "delta_grad == 0.0")
            res = minimize(fun, self.x0, method=self.method,
                           bounds=self.bnds, constraints=cons)
        assert_allclose(res.x, [1.7918552, 1.62895927])
        assert_allclose(res.fun, 1.3857466063348418)

    def test_constraint_dictionary_3(self):
        fun = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2
        cons = [{'type': 'ineq', 'fun': lambda x: x[0] - 2 * x[1] + 2},
                NonlinearConstraint(lambda x: x[0] - x[1], 0, 0)]

        with suppress_warnings() as sup:
            sup.filter(UserWarning, "delta_grad == 0.0")
            res = minimize(fun, self.x0, method=self.method,
                           bounds=self.bnds, constraints=cons)
        assert_allclose(res.x, [1.75, 1.75], rtol=1e-4)
        assert_allclose(res.fun, 1.125, rtol=1e-4)


class TestNewToOld(object):
    method = 'slsqp'
    elec = Elec(n_electrons=2)
    elec.x_opt = np.array([-0.58438468,  0.58438466,  0.73597047,
                           -0.73597044,  0.34180668, -0.34180667])
    brock = BoundedRosenbrock()
    brock.x_opt = [0, 0]
    list_of_problems = [Maratos(),
                        HyperbolicIneq(),
                        Rosenbrock(),
                        IneqRosenbrock(),
                        EqIneqRosenbrock(),
                        elec,
                        brock
                        ]

    def test_list_of_problems(self):

        for prob in self.list_of_problems:

            with suppress_warnings() as sup:
                sup.filter(UserWarning)
                result = minimize(prob.fun, prob.x0,
                                  method=self.method,
                                  bounds=prob.bounds,
                                  constraints=prob.constr)

            assert_array_almost_equal(result.x, prob.x_opt, decimal=3)

    def test_warn_mixed_constraints(self):
        # warns about inefficiency of mixed equality/inequality constraints
        fun = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2 + (x[2] - 0.75)**2
        cons = NonlinearConstraint(lambda x: [x[0] - x[1], x[1] - x[2]],
                                   [1.1, .8], [1.1, 1.4])
        bnds = ((0, None), (0, None), (0, None))
        _assert_warns(OptimizeWarning, minimize, fun, (2, 0, 1),
                      method=self.method, bounds=bnds, constraints=cons)

    def test_warn_ignored_options(self):
        # warns about constraint options being ignored
        fun = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2 + (x[2] - 0.75)**2
        x0 = (2, 0, 1)

        if self.method == "slsqp":
            bnds = ((0, None), (0, None), (0, None))
        else:
            bnds = None

        cons = NonlinearConstraint(lambda x: x[0], 2, np.inf)
        res = minimize(fun, x0, method=self.method,
                       bounds=bnds, constraints=cons)
        # no warnings without constraint options
        assert_allclose(res.fun, 1)

        cons = LinearConstraint([1, 0, 0], 2, np.inf)
        res = minimize(fun, x0, method=self.method,
                       bounds=bnds, constraints=cons)
        # no warnings without constraint options
        assert_allclose(res.fun, 1)

        cons = []
        cons.append(NonlinearConstraint(lambda x: x[0], 2, np.inf,
                                        keep_feasible=True))
        cons.append(NonlinearConstraint(lambda x: x[0], 2, np.inf,
                                        hess=42))
        cons.append(NonlinearConstraint(lambda x: x[0], 2, np.inf,
                                        finite_diff_jac_sparsity=42))
        cons.append(NonlinearConstraint(lambda x: x[0], 2, np.inf,
                                        finite_diff_rel_step=42))
        cons.append(LinearConstraint([1, 0, 0], 2, np.inf,
                                     keep_feasible=True))
        for con in cons:
            _assert_warns(OptimizeWarning, minimize, fun, x0,
                          method=self.method, bounds=bnds, constraints=cons)
