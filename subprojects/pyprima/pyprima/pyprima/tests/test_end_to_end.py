from pyprima import minimize, Bounds, LinearConstraint, NonlinearConstraint
import numpy as np

def obj(x):
    return (x[0] - 1)**2 + (x[1] - 2.5)**2
obj.x0 = np.array([5, 5])
obj.optimal = np.array([1, 2.5])


def test_end_to_end_no_constraints():
    result = minimize(obj, obj.x0, method='cobyla')
    assert np.allclose(result.x, obj.optimal, atol=1e-3)


def test_end_to_end_bounds():
    bounds = Bounds([-5, 10], [5, 10])
    result = minimize(obj, obj.x0, method='cobyla', bounds=bounds)
    assert np.allclose(result.x, np.array([1, 10]), atol=1e-3)


def test_end_to_end_linear_constraints(minimize_with_debugging):
    # x1 + x2 = 5
    A = np.array([[1, 1]])
    b = np.array([5])
    result = minimize_with_debugging(obj, obj.x0, method='cobyla', constraints=[LinearConstraint(A, b, b)])
    assert np.allclose(result.x, np.array([1.75, 3.25]), atol=1e-3)
    assert np.allclose(A @ result.x, b, atol=1e-3)


def test_end_to_end_nonlinear_constraint():
    def cons(x):
        # x1**2 - 0.25 <= 0, i.e. x1 <= 0.5
        return x[0]**2 - 0.25
    result = minimize(obj, obj.x0, method='cobyla', constraints=[NonlinearConstraint(cons, 0, 0)])
    assert np.allclose(result.x, np.array([0.5, 2.5]), atol=1e-3)
    assert np.isclose(cons(result.x), 0, atol=1e-8)
