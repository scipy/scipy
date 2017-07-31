from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.optimize._constraints import (BoxConstraint,
                                         LinearConstraint,
                                         NonlinearConstraint,
                                         _check_kind,
                                         _check_enforce_feasibility,
                                         _reinforce_box_constraint)
from numpy.testing import (TestCase, assert_array_almost_equal,
                           assert_array_equal, assert_array_less,
                           assert_equal, assert_,
                           run_module_suite, assert_allclose, assert_warns,
                           dec)
import warnings
import pytest


class TestCheckKind(TestCase):

    def test_kind_wrong_type(self):
        with pytest.raises(ValueError):
            _check_kind(1, "bla")

    def test_kind_empty(self):
        with pytest.raises(ValueError):
            _check_kind(1, [])

    def test_kind_invalid_format(self):
        with pytest.raises(ValueError):
            _check_kind(3, ["interval", [1, 2, 3]])

    def test_kind_mismatching_ub_lb(self):
        with pytest.raises(ValueError):
            _check_kind(3, ["interval", [1, 2, 3], [1, 2]])

    def test_kind_ub_smaller_than_lb(self):
        with pytest.raises(ValueError):
            _check_kind(3, ["interval", [1, 2, 3], [1, 2, 1]])

    def test_string(self):
        keyword, lb = _check_kind("greater", 3)
        assert_equal(keyword, "greater")
        assert_equal(lb, [0, 0, 0])

    def test_broadcast(self):
        keyword, lb = _check_kind(("greater", 1), 3)
        assert_equal(keyword, "greater")
        assert_equal(lb, [1, 1, 1])


class TestCheckEnforceFeasibility(TestCase):

    def test_wrong_size(self):
        with pytest.raises(ValueError):
            _check_enforce_feasibility([True, True], 3)

    def test_single_value(self):
        f = _check_enforce_feasibility(True, 3)
        assert_array_equal(f, [True, True, True])


class TestReinforceBoxConstraints(TestCase):

    def test_reinforce_box_constraints(self):
        lb = np.array([0, 20, 30])
        ub = np.array([0.5, np.inf, 70])
        enforce_feasibility = np.array([True, False, True],
                                       dtype=bool)
        kind = ("interval", lb, ub)
        x0 = [1, 2, 3]
        x0 = _reinforce_box_constraint(kind, enforce_feasibility, x0)
        assert_array_less(lb[enforce_feasibility], x0[enforce_feasibility])
        assert_array_less(x0[enforce_feasibility], ub[enforce_feasibility])


class TestBoxConstraint(TestCase):

    def test_unfeasible_initial_point(self):
        lb = np.array([0, 20, 30])
        ub = np.array([0.5, np.inf, 70])
        x0 = np.array([1, 2, 3])
        enforce_feasibility = np.array([False, True, True],
                                       dtype=bool)
        kind = ("interval", lb, ub)
        box = BoxConstraint(kind, enforce_feasibility)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            x0_new = box.evaluate_and_initialize(x0)
        assert_((lb[enforce_feasibility] <= x0_new[enforce_feasibility]).all())
        assert_((x0_new[enforce_feasibility] <= ub[enforce_feasibility]).all())

    def test_box_to_linear_conversion(self):
        box = BoxConstraint(("interval", [10, 20, 30], [50, np.inf, 70]))
        x0 = np.array([1, 2, 3])
        x0 = box.evaluate_and_initialize(x0)
        linear = box.to_linear()
        assert_array_equal(linear.A.todense(), np.eye(3))

    def test_box_to_nonlinear_conversion(self):
        box = BoxConstraint(("interval", [10, 20, 30], [50, np.inf, 70]))
        x0 = np.array([1, 2, 3])
        x0 = box.evaluate_and_initialize(x0)
        nonlinear = box.to_nonlinear()
        assert_array_equal(nonlinear.fun(x0), x0)
        assert_array_equal(nonlinear.jac(x0).todense(), np.eye(3))


class TestLinearConstraint(TestCase):

    def test_unfeasible_initial_point(self):
        x0 = np.array([1, 2, 3, 4])
        A = np.array([[1, 2, 3, 4], [5, 0, 0, 6], [7, 0, 8, 0]])
        enforce_feasibility = np.array([True, True, True],
                                       dtype=bool)
        kind = ("less",)
        box = LinearConstraint(A, kind, enforce_feasibility)
        with pytest.raises(ValueError):
            box.evaluate_and_initialize(x0)

    def test_linear_to_nonlinear_conversion(self):
        x0 = np.array([1, 2, 3, 4])
        A = np.array([[1, 2, 3, 4], [5, 0, 0, 6], [7, 0, 8, 0]])
        enforce_feasibility = np.array([False, False, False],
                                       dtype=bool)
        kind = ("less",)
        linear = LinearConstraint(A, kind, enforce_feasibility)
        x0 = linear.evaluate_and_initialize(x0)
        nonlinear = linear.to_nonlinear()
        assert_array_equal(nonlinear.fun(x0), A.dot(x0))
        assert_array_equal(nonlinear.jac(x0), A)


class TestNonlinearConstraint(TestCase):

    def test_unfeasible_initial_point(self):
        x0 = np.array([1, 2, 3, 4])
        A = np.array([[1, 2, 3, 4], [5, 0, 0, 6], [7, 0, 8, 0]])

        def fun(x):
            return A.dot(x)

        def jac(x):
            return A

        enforce_feasibility = np.array([True, True, True],
                                       dtype=bool)
        kind = ("less",)
        box = NonlinearConstraint(fun, kind, jac, None, enforce_feasibility)
        with pytest.raises(ValueError):
            box.evaluate_and_initialize(x0)
