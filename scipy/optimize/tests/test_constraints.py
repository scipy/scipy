from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.optimize import (parse_constraint,
                            BoxConstraint,
                            LinearConstraint,
                            NonlinearConstraint,
                            CanonicalConstraint)
from numpy.testing import (TestCase, assert_array_almost_equal,
                           assert_array_equal, assert_array_less,
                           assert_raises, assert_equal, assert_,
                           run_module_suite, assert_allclose, assert_warns,
                           dec)


class TestParseConstraint(TestCase):

    def test_equality_constraint(self):
        kind = ("equals", [10, 20, 30])
        eq, ineq, val_eq, val_ineq, sign, fun_len = parse_constraint(kind)
        assert_array_equal(eq, [0, 1, 2])
        assert_array_equal(val_eq, [10, 20, 30])
        assert_array_equal(ineq, [])
        assert_array_equal(val_ineq, [])
        assert_array_equal(sign, [])

    def test_greater_constraint(self):
        kind = ("greater", [10, 20, 30])
        eq, ineq, val_eq, val_ineq, sign, fun_len = parse_constraint(kind)
        assert_array_equal(eq, [])
        assert_array_equal(val_eq, [])
        assert_array_equal(ineq, [0, 1, 2])
        assert_array_equal(val_ineq, [10, 20, 30])
        assert_array_equal(sign, [-1, -1, -1])

        kind = ("greater", [10, np.inf, 30])
        eq, ineq, val_eq, val_ineq, sign, fun_len = parse_constraint(kind)
        assert_array_equal(eq, [])
        assert_array_equal(val_eq, [])
        assert_array_equal(ineq, [0, 2])
        assert_array_equal(val_ineq, [10, 30])
        assert_array_equal(sign, [-1, -1])

    def test_less_constraint(self):
        kind = ("less", [10, 20, 30])
        eq, ineq, val_eq, val_ineq, sign, fun_len = parse_constraint(kind)
        assert_array_equal(eq, [])
        assert_array_equal(val_eq, [])
        assert_array_equal(ineq, [0, 1, 2])
        assert_array_equal(val_ineq, [10, 20, 30])
        assert_array_equal(sign, [1, 1, 1])

        kind = ("less", [10, np.inf, 30])
        eq, ineq, val_eq, val_ineq, sign, fun_len = parse_constraint(kind)
        assert_array_equal(eq, [])
        assert_array_equal(val_eq, [])
        assert_array_equal(ineq, [0, 2])
        assert_array_equal(val_ineq, [10, 30])
        assert_array_equal(sign, [1, 1])

    def test_interval_constraint(self):
        kind = ("interval", [10, 20, 30], [50, 60, 70])
        eq, ineq, val_eq, val_ineq, sign, fun_len = parse_constraint(kind)
        assert_array_equal(eq, [])
        assert_array_equal(val_eq, [])
        assert_array_equal(ineq, [0, 1, 2, 0, 1, 2])
        assert_array_equal(val_ineq, [10, 20, 30, 50, 60, 70])
        assert_array_equal(sign, [-1, -1, -1, 1, 1, 1])

        kind = ("interval", [10, 20, 30], [50, 20, 70])
        eq, ineq, val_eq, val_ineq, sign, fun_len = parse_constraint(kind)
        assert_array_equal(eq, [1])
        assert_array_equal(val_eq, [20])
        assert_array_equal(ineq, [0, 2, 0, 2])
        assert_array_equal(val_ineq, [10, 30, 50, 70])
        assert_array_equal(sign, [-1, -1, 1, 1])

        kind = ("interval", [10, 20, 30], [50, 20, np.inf])
        eq, ineq, val_eq, val_ineq, sign, fun_len = parse_constraint(kind)
        assert_array_equal(eq, [1])
        assert_array_equal(val_eq, [20])
        assert_array_equal(ineq, [0, 2, 0])
        assert_array_equal(val_ineq, [10, 30, 50])
        assert_array_equal(sign, [-1, -1, 1])

        kind = ("interval", [-np.inf, 20, 30], [50, 20, np.inf])
        eq, ineq, val_eq, val_ineq, sign, fun_len = parse_constraint(kind)
        assert_array_equal(eq, [1])
        assert_array_equal(val_eq, [20])
        assert_array_equal(ineq, [2, 0])
        assert_array_equal(val_ineq, [30, 50])
        assert_array_equal(sign, [-1, 1])

    def test_exceptions(self):
        assert_raises(ValueError, parse_constraint, ("blbalbda",))
        assert_raises(ValueError, parse_constraint, ("interval", [1, 2, 3], [1, 2]))
        assert_raises(ValueError, parse_constraint, ("interval", [1, 2, 3], [1, 2, 1]))


class TestConversions(TestCase):

    def test_box_to_linear_conversion(self):
        box = BoxConstraint(("interval", [10, 20, 30], [50, np.inf, 70]))
        linear = box.to_linear(sparse=True)
        assert_array_equal(linear.A.todense(), np.eye(3))

    def test_linear_to_nonlinear_conversion(self):
        A = np.array([[1, 2, 3, 4], [5, 0, 0, 6], [7, 0, 8, 0]])
        linear = LinearConstraint(A, ("interval", [10, 20, 30], [50, np.inf, 70]))
        nonlinear = linear.to_nonlinear()
        x = [1, 2, 3, 4]
        assert_array_equal(nonlinear.fun(x), A.dot(x))
        assert_array_equal(nonlinear.jac(x), A)

    def test_box_to_canonical_conversion(self):
        box = BoxConstraint(("interval", [10, 20, 30], [50, np.inf, 70]))
        canonical = box.to_canonical(sparse=True)

        x = [1, 2, 3]
        assert_array_equal(canonical.n_eq, 0)
        assert_array_equal(canonical.constr_eq(x), [])
        assert_array_equal(canonical.jac_eq(x), np.empty((0, 3)))
        assert_array_equal(canonical.n_ineq, 5)
        assert_array_equal(canonical.constr_ineq(x), [10-1,
                                                      20-2,
                                                      30-3,
                                                      1-50,
                                                      3-70])
        assert_array_equal(canonical.jac_ineq(x).todense(), [[-1, 0, 0],
                                                             [0, -1, 0],
                                                             [0, 0, -1],
                                                             [1, 0, 0],
                                                             [0, 0, 1]])
        assert_array_equal(canonical.hess, None)

    def test_linear_to_canonical_conversion(self):
        A = np.array([[1, 2, 3, 4], [5, 0, 0, 6], [7, 0, 8, 0]])
        linear = LinearConstraint(A, ("interval", [10, 20, 30], [10, np.inf, 70]))
        canonical = linear.to_canonical()

        x = [1, 2, 3, 4]
        assert_array_equal(canonical.n_eq, 1)
        assert_array_equal(canonical.constr_eq(x), [1+4+9+16-10])
        assert_array_equal(canonical.jac_eq(x), [[1, 2, 3, 4]])
        assert_array_equal(canonical.n_ineq, 3)
        assert_array_equal(canonical.constr_ineq(x), [20-5*1-6*4,
                                                      30-7*1-8*3,
                                                      7*1+8*3-70])
        assert_array_equal(canonical.jac_ineq(x), [[-5, 0, 0, -6],
                                                   [-7, 0, -8, 0],
                                                   [7, 0, 8, 0]])
        assert_array_equal(canonical.hess, None)

    def test_nonlinear_to_canonical_conversion(self):
        f1 = 10
        g1 = np.array([1, 2, 3, 4])
        H1 = np.eye(4)

        f2 = 1
        g2 = np.array([1, 1, 1, 1])
        H2 = np.zeros((4, 4))

        f3 = 12
        g3 = np.array([1, 0, 0, 1])
        H3 = np.diag([1, 2, 3, 4])

        def fun(x):
            return np.array([f1 + g1.dot(x) + 1/2*H1.dot(x).dot(x),
                             f2 + g2.dot(x) + 1/2*H2.dot(x).dot(x),
                             f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x)])

        def jac(x):
            return np.vstack([g1 + H1.dot(x),
                              g2 + H2.dot(x),
                              g3 + H3.dot(x)])

        def hess(x, v):
            return v[0]*H1 + v[1]*H2 + v[2]*H3

        nonlinear = NonlinearConstraint(fun, jac, hess, ("interval", [10, 20, 30],
                                                                     [10, np.inf, 70]))
        canonical = nonlinear.to_canonical()
        x = [1, 2, 3, 4]
        assert_array_equal(canonical.n_eq, 1)
        assert_array_equal(canonical.constr_eq(x), [f1 + g1.dot(x) + 1/2*H1.dot(x).dot(x) - 10])
        assert_array_equal(canonical.jac_eq(x), np.atleast_2d(g1 + H1.dot(x)))
        assert_array_equal(canonical.n_ineq, 3)
        assert_array_equal(canonical.constr_ineq(x), [20-(f2 + g2.dot(x) + 1/2*H2.dot(x).dot(x)),
                                                      30-(f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x)),
                                                      f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x) - 70])
        assert_array_equal(canonical.jac_ineq(x), np.vstack([-(g2 + H2.dot(x)),
                                                             -(g3 + H3.dot(x)),
                                                             g3 + H3.dot(x)]))
        v_eq = np.array([10])
        v_ineq = np.array([5, 6, 3])
        assert_array_equal(canonical.hess(x, v_eq, v_ineq), 10*H1 + (-5)*H2 + (-6+3)*H3)
        v_eq = np.array([50])
        v_ineq = np.array([4, -2, 30])
        assert_array_equal(canonical.hess(x, v_eq, v_ineq), 50*H1 + (-4)*H2 + (2+30)*H3)

        nonlinear = NonlinearConstraint(fun, jac, hess, ("interval", [10, 20, 30],
                                                                     [20, 20, 70]))
        canonical = nonlinear.to_canonical()
        x = [1, 2, 3, 4]
        assert_array_equal(canonical.n_eq, 1)
        assert_array_equal(canonical.constr_eq(x), [f2+ g2.dot(x) + 1/2*H2.dot(x).dot(x) - 20])
        assert_array_equal(canonical.jac_eq(x), np.atleast_2d(g2 + H2.dot(x)))
        assert_array_equal(canonical.n_ineq, 4)
        assert_array_equal(canonical.constr_ineq(x), [10-(f1 + g1.dot(x) + 1/2*H1.dot(x).dot(x)),
                                                      30-(f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x)),
                                                      f1 + g1.dot(x) + 1/2*H1.dot(x).dot(x) - 20,
                                                      f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x) - 70])
        assert_array_equal(canonical.jac_ineq(x), np.vstack([-(g1 + H1.dot(x)),
                                                             -(g3 + H3.dot(x)),
                                                             g1 + H1.dot(x),
                                                             g3 + H3.dot(x)]))
        v_eq = np.array([10])
        v_ineq = np.array([5, 6, 3, 12])
        assert_array_equal(canonical.hess(x, v_eq, v_ineq), (-5+3)*H1 + 10*H2 + (-6+12)*H3)
        v_eq = np.array([50])
        v_ineq = np.array([4, -2, 30, 2])
        assert_array_equal(canonical.hess(x, v_eq, v_ineq), (-4+30)*H1 + 50*H2 + (2+2)*H3)
