from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (assert_, assert_allclose, run_module_suite,
                           assert_equal, assert_raises)
from scipy.integrate import solve_ivp


def fun_rational(x, y):
    return np.array([y[1] / x,
                     y[1] * (y[0] + 2 * y[1] - 1) / (x * (y[0] - 1))])


def jac_rational(x, y):
    return np.array([
        [0, 1 / x],
        [-2 * y[1] ** 2 / (x * (y[0] - 1) ** 2),
         (y[0] + 4 * y[1] - 1) / (x * (y[0] - 1))]
    ])


def sol_rational(x):
    return np.vstack((x / (x + 10), 10 * x / (x + 10)**2))


# Swapped y_1 and y_2.
def fun_rational_swapped(x, y):
    return np.array([y[0] / x,
                     y[0] * (y[1] + 2 * y[0] - 1) / (x * (y[1] - 1))])


def jac_rational_swapped(x, y):
    return np.array([
        [0, 1 / x],
        [-2 * y[0] ** 2 / (x * (y[1] - 1) ** 2),
         (y[1] + 4 * y[0] - 1) / (x * (y[1] - 1))]
    ])


def event_rational_1(x, y):
    return y[0] - y[1] ** 0.7


def event_rational_2(x, y):
    return y[1] ** 0.6 - y[0]


def event_rational_3(x, y):
    return x - 7.4


def compute_error(y, y_true, rtol, atol):
    e = (y - y_true) / (atol + rtol * y_true)
    return np.sqrt(np.sum(e**2, axis=0) / e.shape[0])


def test_rk():
    rtol = 1e-3
    atol = 1e-6
    for method in ['RK23', 'RK45']:
        for x_span in ([5, 9], [5, 1]):
            res = solve_ivp(fun_rational, x_span, [1/3, 2/9], rtol=rtol,
                            atol=atol, method=method)
            assert_equal(res.x[0], x_span[0])
            assert_equal(res.x[-1], x_span[-1])
            assert_(res.x_events is None)
            assert_(res.success)
            assert_equal(res.status, 0)

            y_true = sol_rational(res.x)
            e = compute_error(res.y, y_true, rtol, atol)
            assert_(np.all(e < 0.2))

            xc = np.linspace(*x_span)
            yc_true = sol_rational(xc)
            yc = res.sol(xc)

            e = compute_error(yc, yc_true, rtol, atol)
            assert_(np.all(e < 0.2))

            assert_allclose(res.sol(res.x), res.y, rtol=1e-15, atol=1e-15)
            assert_allclose(res.sol(res.x, 1), res.yp, rtol=1e-15, atol=1e-13)


def test_radau():
    rtol = 1e-3
    atol = 1e-6

    for x_span in ([5, 9], [5, 1]):
        for jac in [None, jac_rational]:
            res = solve_ivp(fun_rational, x_span, [1/3, 2/9], rtol=rtol,
                            atol=atol, method='Radau', jac=jac)
            assert_equal(res.x[0], x_span[0])
            assert_equal(res.x[-1], x_span[-1])
            assert_(res.x_events is None)
            assert_(res.success)
            assert_equal(res.status, 0)

            y_true = sol_rational(res.x)
            e = compute_error(res.y, y_true, rtol, atol)
            assert_(np.all(e < 0.2))

            xc = np.linspace(*x_span)
            yc_true = sol_rational(xc)
            yc = res.sol(xc)

            e = compute_error(yc, yc_true, rtol, atol)
            assert_(np.all(e < 0.2))

            assert_allclose(res.sol(res.x), res.y, rtol=1e-15, atol=1e-15)

    # This matrix swaps the variables.
    M = np.array([[0, 1], [1, 0]])
    for x_span in ([5, 9], [5, 1]):
        for jac in [None, jac_rational_swapped]:
            res = solve_ivp(fun_rational_swapped, x_span, [2/9, 1/3],
                            rtol=rtol, atol=atol, method='Radau', M=M, jac=jac)
            assert_equal(res.x[0], x_span[0])
            assert_equal(res.x[-1], x_span[-1])
            assert_(res.x_events is None)
            assert_(res.success)
            assert_equal(res.status, 0)

            y_true = sol_rational(res.x)
            y_true[[0, 1]] = y_true[[1, 0]]
            e = compute_error(res.y, y_true, rtol, atol)
            assert_(np.all(e < 0.2))

            xc = np.linspace(*x_span)
            yc_true = sol_rational(xc)
            yc_true[[0, 1]] = yc_true[[1, 0]]
            yc = res.sol(xc)

            e = compute_error(yc, yc_true, rtol, atol)
            assert_(np.all(e < 0.2))

            assert_allclose(res.sol(res.x), res.y, rtol=1e-15, atol=1e-15)


def test_events():
    event_rational_3.terminate = True

    for method in ['RK23', 'RK45']:
        res = solve_ivp(fun_rational, [5, 8], [1/3, 2/9], method=method,
                        events=(event_rational_1, event_rational_2))
        assert_equal(res.status, 0)
        assert_equal(res.x_events[0].size, 1)
        assert_equal(res.x_events[1].size, 1)
        assert_(5.3 < res.x_events[0][0] < 5.7)
        assert_(7.3 < res.x_events[1][0] < 7.7)

        event_rational_1.direction = 1
        event_rational_2.direction = 1
        res = solve_ivp(fun_rational, [5, 8], [1 / 3, 2 / 9], method=method,
                        events=(event_rational_1, event_rational_2))
        assert_equal(res.status, 0)
        assert_equal(res.x_events[0].size, 1)
        assert_equal(res.x_events[1].size, 0)
        assert_(5.3 < res.x_events[0][0] < 5.7)

        event_rational_1.direction = -1
        event_rational_2.direction = -1
        res = solve_ivp(fun_rational, [5, 8], [1 / 3, 2 / 9], method=method,
                        events=(event_rational_1, event_rational_2))
        assert_equal(res.status, 0)
        assert_equal(res.x_events[0].size, 0)
        assert_equal(res.x_events[1].size, 1)
        assert_(7.3 < res.x_events[1][0] < 7.7)

        event_rational_1.direction = 0
        event_rational_2.direction = 0

        res = solve_ivp(fun_rational, [5, 8], [1 / 3, 2 / 9], method=method,
                        events=(event_rational_1, event_rational_2,
                                event_rational_3))
        assert_equal(res.status, 1)
        assert_equal(res.x_events[0].size, 1)
        assert_equal(res.x_events[1].size, 0)
        assert_equal(res.x_events[2].size, 1)
        assert_(5.3 < res.x_events[0][0] < 5.7)
        assert_(7.3 < res.x_events[2][0] < 7.5)

        # Also test that termination by event doesn't break interpolants.
        xc = np.linspace(res.x[0], res.x[-1])
        yc_true = sol_rational(xc)
        yc = res.sol(xc)
        e = compute_error(yc, yc_true, 1e-3, 1e-6)
        assert_(np.all(e < 0.2))

    # Test in backward direction.
    event_rational_1.direction = 0
    event_rational_2.direction = 0
    for method in ['RK23', 'RK45', 'Radau']:
        res = solve_ivp(fun_rational, [8, 5], [4/9, 20/81], method=method,
                        events=(event_rational_1, event_rational_2))
        assert_equal(res.status, 0)
        assert_equal(res.x_events[0].size, 1)
        assert_equal(res.x_events[1].size, 1)
        assert_(5.3 < res.x_events[0][0] < 5.7)
        assert_(7.3 < res.x_events[1][0] < 7.7)

        event_rational_1.direction = -1
        event_rational_2.direction = -1
        res = solve_ivp(fun_rational, [8, 5], [4/9, 20/81], method=method,
                        events=(event_rational_1, event_rational_2))
        assert_equal(res.status, 0)
        assert_equal(res.x_events[0].size, 1)
        assert_equal(res.x_events[1].size, 0)
        assert_(5.3 < res.x_events[0][0] < 5.7)

        event_rational_1.direction = 1
        event_rational_2.direction = 1
        res = solve_ivp(fun_rational, [8, 5], [4/9, 20/81], method=method,
                        events=(event_rational_1, event_rational_2))
        assert_equal(res.status, 0)
        assert_equal(res.x_events[0].size, 0)
        assert_equal(res.x_events[1].size, 1)
        assert_(7.3 < res.x_events[1][0] < 7.7)

        event_rational_1.direction = 0
        event_rational_2.direction = 0

        res = solve_ivp(fun_rational, [8, 5], [4/9, 20/81], method=method,
                        events=(event_rational_1, event_rational_2,
                                event_rational_3))
        assert_equal(res.status, 1)
        assert_equal(res.x_events[0].size, 0)
        assert_equal(res.x_events[1].size, 1)
        assert_equal(res.x_events[2].size, 1)
        assert_(7.3 < res.x_events[1][0] < 7.7)
        assert_(7.3 < res.x_events[2][0] < 7.5)

        # Also test that termination by event doesn't break interpolants.
        xc = np.linspace(res.x[-1], res.x[0])
        yc_true = sol_rational(xc)
        yc = res.sol(xc)
        e = compute_error(yc, yc_true, 1e-3, 1e-6)
        assert_(np.all(e < 0.2))


def test_parameters_validation():
    assert_raises(ValueError, solve_ivp, fun_rational, [1, 2], [0, 0],
                  method='M')
    assert_raises(ValueError, solve_ivp, fun_rational, [1, 1], [0, 0])
    assert_raises(ValueError, solve_ivp, fun_rational, [1, 2], [[0, 0]])
    assert_raises(ValueError, solve_ivp, lambda x, y: np.zeros(3), [1, 2],
                  [0, 0])
    assert_raises(ValueError, solve_ivp, fun_rational, [1, 2], [0, 0],
                  M=np.identity(2), method='RK23')
    assert_raises(ValueError, solve_ivp, fun_rational, [1, 2], [0, 0],
                  method='Radau', M=np.identity(3))
    assert_raises(ValueError, solve_ivp, fun_rational, [1, 2], [0, 0],
                  method='Radau', jac=lambda x, y: np.identity(3))
    assert_raises(ValueError, solve_ivp, fun_rational, [1, 2], [0, 0],
                  method='Radau', jac=np.identity(3))


if __name__ == '__main__':
    run_module_suite()
