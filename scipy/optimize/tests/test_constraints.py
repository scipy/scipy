from __future__ import division, print_function, absolute_import
import numpy as np
import scipy.sparse as spc
from scipy.optimize._constraints import (BoxConstraint,
                                         LinearConstraint,
                                         NonlinearConstraint)
from numpy.testing import (TestCase, assert_array_almost_equal,
                           assert_array_equal, assert_array_less,
                           assert_raises, assert_equal, assert_,
                           run_module_suite, assert_allclose, assert_warns,
                           dec)


class TestConversions(TestCase):

    def test_box_to_linear_conversion(self):
        box = BoxConstraint(("interval", [10, 20, 30], [50, np.inf, 70]))
        linear = box.to_linear(sparse_jacobian=True)
        assert_array_equal(linear.A.todense(), np.eye(3))

    def test_linear_to_nonlinear_conversion(self):
        A = np.array([[1, 2, 3, 4], [5, 0, 0, 6], [7, 0, 8, 0]])
        linear = LinearConstraint(A, ("interval", [10, 20, 30], [50, np.inf, 70]))
        nonlinear = linear.to_nonlinear()
        x = [1, 2, 3, 4]
        assert_array_equal(nonlinear.fun(x), A.dot(x))
        assert_array_equal(nonlinear.jac(x), A)



class TestReinforceBoxConstraints(TestCase):

    def test_reinforce_box_constraints(self):
        lb = np.array([0, 20, 30])
        ub = np.array([0.5, np.inf, 70])
        enforce_feasibility = np.array([True, False, True],
                                   dtype=bool)
        box = BoxConstraint(("interval", lb, ub),
                            enforce_feasibility)
        x0 = [1, 2, 3]

        x0 = box.reinforce_constraint(x0)
        assert_array_less(lb[enforce_feasibility], x0[enforce_feasibility])
        assert_array_less(x0[enforce_feasibility], ub[enforce_feasibility])
