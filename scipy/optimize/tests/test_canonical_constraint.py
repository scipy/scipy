from __future__ import division, print_function, absolute_import
import numpy as np
import scipy.sparse as spc
from scipy.optimize._constraints import (NonlinearConstraint,
                                         LinearConstraint,
                                         BoxConstraint)
from scipy.optimize._canonical_constraint import (parse_constraint,
                                                  CanonicalConstraint,
                                                  box_to_canonical,
                                                  linear_to_canonical,
                                                  nonlinear_to_canonical,
                                                  concatenate_canonical_constraints)
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


class TestConversions(TestCase):

    def test_box_to_canonical_conversion(self):
        box = BoxConstraint(("interval", [10, 20, 30], [50, np.inf, 70]),
                                  [True, False, True])
        canonical = box_to_canonical(box, sparse_jacobian=True)

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
        assert_array_equal(canonical.enforce_feasibility,
                           [True, False, True, True, True])

    def test_linear_to_canonical_conversion(self):
        A = np.array([[1, 2, 3, 4], [5, 0, 0, 6], [7, 0, 8, 0]])
        linear = LinearConstraint(A, ("interval",
                                      [10, 20, 30],
                                      [10, np.inf, 70]),
                                  [True, False, True])
        canonical = linear_to_canonical(linear)

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
        assert_array_equal(canonical.enforce_feasibility, [False, True, True])

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

        nonlinear = NonlinearConstraint(fun,
                                        ("interval",
                                         [10, 20, 30],
                                         [10, np.inf, 70]),
                                        jac, hess,
                                        [True, False, True])
        canonical = nonlinear_to_canonical(nonlinear)
        x = [1, 2, 3, 4]
        assert_array_equal(canonical.n_eq, 1)
        assert_array_equal(canonical.constr_eq(x),
                           [f1 + g1.dot(x) + 1/2*H1.dot(x).dot(x) - 10])
        assert_array_equal(canonical.jac_eq(x), np.atleast_2d(g1 + H1.dot(x)))
        assert_array_equal(canonical.n_ineq, 3)
        assert_array_equal(canonical.constr_ineq(x),
                           [20-(f2 + g2.dot(x) + 1/2*H2.dot(x).dot(x)),
                            30-(f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x)),
                            f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x) - 70])
        assert_array_equal(canonical.jac_ineq(x), np.vstack([-(g2 + H2.dot(x)),
                                                             -(g3 + H3.dot(x)),
                                                             g3 + H3.dot(x)]))
        v_eq = np.array([10])
        v_ineq = np.array([5, 6, 3])
        assert_array_equal(canonical.hess(x, v_eq, v_ineq),
                           10*H1 + (-5)*H2 + (-6+3)*H3)
        v_eq = np.array([50])
        v_ineq = np.array([4, -2, 30])
        assert_array_equal(canonical.hess(x, v_eq, v_ineq),
                           50*H1 + (-4)*H2 + (2+30)*H3)
        assert_array_equal(canonical.enforce_feasibility, [False, True, True])

        nonlinear = NonlinearConstraint(fun,
                                        ("interval",
                                         [10, 20, 30],
                                         [20, 20, 70]),
                                        jac, hess,
                                        [True, False, True])
        canonical = nonlinear_to_canonical(nonlinear)
        x = [1, 2, 3, 4]
        assert_array_equal(canonical.n_eq, 1)
        assert_array_equal(canonical.constr_eq(x),
                           [f2 + g2.dot(x) + 1/2*H2.dot(x).dot(x) - 20])
        assert_array_equal(canonical.jac_eq(x), np.atleast_2d(g2 + H2.dot(x)))
        assert_array_equal(canonical.n_ineq, 4)
        assert_array_equal(canonical.constr_ineq(x),
                           [10-(f1 + g1.dot(x) + 1/2*H1.dot(x).dot(x)),
                            30-(f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x)),
                            f1 + g1.dot(x) + 1/2*H1.dot(x).dot(x) - 20,
                            f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x) - 70])
        assert_array_equal(canonical.jac_ineq(x), np.vstack([-(g1 + H1.dot(x)),
                                                             -(g3 + H3.dot(x)),
                                                             g1 + H1.dot(x),
                                                             g3 + H3.dot(x)]))
        v_eq = np.array([10])
        v_ineq = np.array([5, 6, 3, 12])
        assert_array_equal(canonical.hess(x, v_eq, v_ineq),
                           (-5+3)*H1 + 10*H2 + (-6+12)*H3)
        v_eq = np.array([50])
        v_ineq = np.array([4, -2, 30, 2])
        assert_array_equal(canonical.hess(x, v_eq, v_ineq),
                           (-4+30)*H1 + 50*H2 + (2+2)*H3)
        assert_array_equal(canonical.enforce_feasibility, [True, True, True, True])


class TestConcatenateConstraints(TestCase):

    def test_concatenate_constraints_sparse(self):
        # Define first constraint
        A = np.array([[1, 2, 3, 4],
                      [5, 0, 0, 6],
                      [7, 0, 8, 0]])
        linear = LinearConstraint(A, ("interval",
                                      [10, 20, 30],
                                      [10, np.inf, 70]),
                                  [True, False, True])

        # Define second constraint
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

        nonlinear = NonlinearConstraint(fun,
                                        ("interval",
                                         [10, 20, 30],
                                         [10, np.inf, 70]),
                                        jac, hess,
                                        [True, False, True])

        # Define third constraint
        box = BoxConstraint(("interval",
                             [10, 20, 30, -np.inf],
                             [50, np.inf, 70, np.inf]), True)

        list_configurations = [(None, None, None, None, None),
                               (True, True, True, True, True),
                               (True, True, True, True, None),
                               (False, False, False, False, False),
                               (False, False, False, False, None),
                               (False, False, True, False, None),
                               (True, False, True, False, False)]

        for conf in list_configurations:
            sparse_jacobian = conf[4]
            # Concatenate constraint
            canonical = concatenate_canonical_constraints(
                [linear_to_canonical(linear, conf[0]),
                 nonlinear_to_canonical(nonlinear, conf[1]),
                 box_to_canonical(box, conf[2]),
                 nonlinear_to_canonical(nonlinear, conf[3])],
                sparse_jacobian=conf[4])

            # Test number of constraints
            assert_equal(canonical.n_eq, 1 + 1 + 0 + 1)
            assert_equal(canonical.n_ineq, 3 + 3 + 5 + 3)

            # Test constraint evaluation
            x = [1, 2, 3, 4]
            assert_array_almost_equal(canonical.constr_eq(x),
                                      [1+4+9+16-10,
                                       f1 + g1.dot(x) + 1/2*H1.dot(x).dot(x) - 10,
                                       f1 + g1.dot(x) + 1/2*H1.dot(x).dot(x) - 10])
            assert_array_almost_equal(canonical.constr_ineq(x),
                                      [20-5*1-6*4,
                                       30-7*1-8*3,
                                       7*1+8*3-70,
                                       20-(f2 + g2.dot(x) + 1/2*H2.dot(x).dot(x)),
                                       30-(f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x)),
                                       f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x) - 70,
                                       10-1,
                                       20-2,
                                       30-3,
                                       1-50,
                                       3-70,
                                       20-(f2 + g2.dot(x) + 1/2*H2.dot(x).dot(x)),
                                       30-(f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x)),
                                       f3 + g3.dot(x) + 1/2*H3.dot(x).dot(x) - 70])

            # Test Jacobian Evaluation
            j_eq = canonical.jac_eq(x)
            if sparse_jacobian is not None:
                assert_(spc.issparse(j_eq) == sparse_jacobian)
            if spc.issparse(j_eq):
                j_eq = j_eq.todense()
            assert_array_almost_equal(j_eq,
                                      np.vstack([np.array([1, 2, 3, 4]),
                                                 g1 + H1.dot(x),
                                                 g1 + H1.dot(x)]))
            j_ineq = canonical.jac_ineq(x)
            if sparse_jacobian is not None:
                assert_(spc.issparse(j_ineq) == sparse_jacobian)
            if spc.issparse(j_ineq):
                j_ineq = j_ineq.todense()
            assert_array_almost_equal(j_ineq,
                                      np.vstack([np.array([[-5, 0, 0, -6],
                                                           [-7, 0, -8, 0],
                                                           [7, 0, 8, 0]]),
                                                 -(g2 + H2.dot(x)),
                                                 -(g3 + H3.dot(x)),
                                                 g3 + H3.dot(x),
                                                 np.array([[-1, 0, 0, 0],
                                                           [0, -1, 0, 0],
                                                           [0, 0, -1, 0],
                                                           [1, 0, 0, 0],
                                                           [0, 0, 1, 0]]),
                                                 -(g2 + H2.dot(x)),
                                                 -(g3 + H3.dot(x)),
                                                 g3 + H3.dot(x)]))

            # Test Hessian Evaluation
            v_eq = np.array([1, 2, 3])
            v_ineq = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
            H = canonical.hess(x, v_eq, v_ineq)
            H_expected = 2*H1 + 3*H1 - 4*H2 - 5*H3 + \
                         6*H3 - 12*H2 - 13*H3 + 14*H3
            np.random.seed(1)
            for i in range(10):
                p = np.random.uniform(size=4)
                assert_array_almost_equal(H.dot(p), H_expected.dot(p))

            # Test feasible constraint list evaluation
            assert_array_equal(canonical.enforce_feasibility,
                               [False, True, True,
                                False, True, True,
                                True, True, True, True, True,
                                False, True, True])
