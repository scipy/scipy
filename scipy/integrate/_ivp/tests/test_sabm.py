from itertools import product
from numpy.testing import (assert_, assert_allclose,
                           assert_equal, assert_no_warnings, suppress_warnings)
import pytest
from pytest import raises as assert_raises
import numpy as np
from scipy.integrate import solve_ivp, SABM
from scipy.integrate import OdeSolution
from scipy.integrate._ivp.base import ConstantDenseOutput


def fun_zero(t, y):
    return np.zeros_like(y)


def fun_linear(t, y):
    return np.array([-y[0] - 5 * y[1], y[0] + y[1]])


def sol_linear(t):
    return np.vstack((-5 * np.sin(2 * t),
                      2 * np.cos(2 * t) + np.sin(2 * t)))


def fun_rational(t, y):
    return np.array([y[1] / t,
                     y[1] * (y[0] + 2 * y[1] - 1) / (t * (y[0] - 1))])


def sol_rational(t):
    return np.asarray((t / (t + 10), 10 * t / (t + 10) ** 2))


def fun_rational_vectorized(t, y):
    return np.vstack((y[1] / t,
                      y[1] * (y[0] + 2 * y[1] - 1) / (t * (y[0] - 1))))


def fun_dae(t, y, x):
    return np.array([-2*y[0] + 10])


def gun_dae(t, y, x):
    return np.array([x[0]+y[0] - 5 - np.exp(-t) - np.exp(-2*t)])


def sol_dae(t):
    return np.array([5+np.exp(-2*t),
                     np.exp(-t)])


def compute_error(y, y_true, rtol, atol):
    e = (y - y_true) / (atol + rtol * np.abs(y_true))
    return np.linalg.norm(e, axis=0) / np.sqrt(e.shape[0])


def test_integration_ODE():
    rtol = 1e-3
    atol = 1e-6
    y0 = [1/3, 2/9]
    method = 'SABM'

    for vectorized, t_span, mode in product(
            [True, False],
            [[5, 9], [5, 1]],
            ["Explicit", "Implicit"]):

        if vectorized:
            fun = fun_rational_vectorized
        else:
            fun = fun_rational

        with suppress_warnings() as sup:
            sup.filter(UserWarning,
                       "The following arguments have no effect for a chosen "
                       "solver: `jac`")
            res = solve_ivp(fun, t_span, y0, rtol=rtol, mode=mode,
                            atol=atol, method=method, dense_output=True,
                            vectorized=vectorized)

        assert_equal(res.t[0], t_span[0])
        assert_(res.t_events is None)
        assert_(res.y_events is None)
        assert_(res.success)
        assert_equal(res.status, 0)

        assert_equal(res.nlu, 0)

        y_true = sol_rational(res.t)
        e = compute_error(res.y, y_true, rtol, atol)
        assert_(np.all(e < 5))

        tc = np.linspace(*t_span)
        yc_true = sol_rational(tc)
        yc = res.sol(tc)

        e = compute_error(yc, yc_true, rtol, atol)
        assert_(np.all(e < 5))

        tc = (t_span[0] + t_span[-1]) / 2
        yc_true = sol_rational(tc)
        yc = res.sol(tc)

        e = compute_error(yc, yc_true, rtol, atol)
        assert_(np.all(e < 5))


def test_integration_DAE():
    rtol = 1e-3
    atol = 1e-6
    method = 'SABM'

    for vectorized, t_span, mode in product(
            [True, False],
            [[5, 7], [5, 3]],
            ["Explicit", "Implicit"]):

        y0 = [5+np.exp(-2*t_span[0]), np.exp(-t_span[0])]
        fun = fun_dae
        gun = gun_dae
        with suppress_warnings() as sup:
            sup.filter(UserWarning,
                       "The following arguments have no effect for a chosen "
                       "solver: `jac`")
            res = solve_ivp(fun, t_span, y0, method=method, gun=gun,
                            num_diff=1, rtol=rtol, atol=atol, mode=mode,
                            dense_output=True, vectorized=vectorized)

        assert_equal(res.t[0], t_span[0])
        assert_(res.t_events is None)
        assert_(res.y_events is None)
        assert_(res.success)
        assert_equal(res.status, 0)

        assert_equal(res.nlu, 0)

        y_true = sol_dae(res.t)
        e = compute_error(res.y, y_true, rtol, atol)
        assert_(np.all(e < 5))

        tc = np.linspace(*t_span)
        yc_true = sol_dae(tc)
        yc = res.sol(tc)

        e = compute_error(yc, yc_true, rtol, atol)
        assert_(np.all(e < 5))

        tc = (t_span[0] + t_span[-1]) / 2
        yc_true = sol_dae(tc)
        yc = res.sol(tc)

        e = compute_error(yc, yc_true, rtol, atol)
        assert_(np.all(e < 5))


def test_events():
    def event_rational_1(t, y):
        return y[0] - y[1] ** 0.7

    def event_rational_2(t, y):
        return y[1] ** 0.6 - y[0]

    def event_rational_3(t, y):
        return t - 7.4

    event_rational_3.terminal = True
    method = 'SABM'

    res = solve_ivp(fun_rational, [5, 8], [1/3, 2/9], method=method,
                    events=(event_rational_1, event_rational_2))
    assert_equal(res.status, 0)
    assert_equal(res.t_events[0].size, 1)
    assert_equal(res.t_events[1].size, 1)
    assert_(5.3 < res.t_events[0][0] < 5.7)
    assert_(7.3 < res.t_events[1][0] < 7.7)

    assert_equal(res.y_events[0].shape, (1, 2))
    assert_equal(res.y_events[1].shape, (1, 2))
    assert np.isclose(
            event_rational_1(res.t_events[0][0], res.y_events[0][0]), 0)
    assert np.isclose(
            event_rational_2(res.t_events[1][0], res.y_events[1][0]), 0)

    event_rational_1.direction = 1
    event_rational_2.direction = 1
    res = solve_ivp(fun_rational, [5, 8], [1 / 3, 2 / 9], method=method,
                    events=(event_rational_1, event_rational_2))
    assert_equal(res.status, 0)
    assert_equal(res.t_events[0].size, 1)
    assert_equal(res.t_events[1].size, 0)
    assert_(5.3 < res.t_events[0][0] < 5.7)
    assert_equal(res.y_events[0].shape, (1, 2))
    assert_equal(res.y_events[1].shape, (0,))
    assert np.isclose(
            event_rational_1(res.t_events[0][0], res.y_events[0][0]), 0)

    event_rational_1.direction = -1
    event_rational_2.direction = -1
    res = solve_ivp(fun_rational, [5, 8], [1 / 3, 2 / 9], method=method,
                    events=(event_rational_1, event_rational_2))
    assert_equal(res.status, 0)
    assert_equal(res.t_events[0].size, 0)
    assert_equal(res.t_events[1].size, 1)
    assert_(7.3 < res.t_events[1][0] < 7.7)
    assert_equal(res.y_events[0].shape, (0,))
    assert_equal(res.y_events[1].shape, (1, 2))
    assert np.isclose(
            event_rational_2(res.t_events[1][0], res.y_events[1][0]), 0)

    event_rational_1.direction = 0
    event_rational_2.direction = 0

    res = solve_ivp(fun_rational, [5, 8], [1 / 3, 2 / 9], method=method,
                    events=(event_rational_1, event_rational_2,
                            event_rational_3), dense_output=True)
    assert_equal(res.status, 1)
    assert_equal(res.t_events[0].size, 1)
    assert_equal(res.t_events[1].size, 0)
    assert_equal(res.t_events[2].size, 1)
    assert_(5.3 < res.t_events[0][0] < 5.7)
    assert_(7.3 < res.t_events[2][0] < 7.5)
    assert_equal(res.y_events[0].shape, (1, 2))
    assert_equal(res.y_events[1].shape, (0,))
    assert_equal(res.y_events[2].shape, (1, 2))
    assert np.isclose(
            event_rational_1(res.t_events[0][0], res.y_events[0][0]), 0)
    assert np.isclose(
            event_rational_3(res.t_events[2][0], res.y_events[2][0]), 0)

    res = solve_ivp(fun_rational, [5, 8], [1 / 3, 2 / 9], method=method,
                    events=event_rational_1, dense_output=True)
    assert_equal(res.status, 0)
    assert_equal(res.t_events[0].size, 1)
    assert_(5.3 < res.t_events[0][0] < 5.7)

    assert_equal(res.y_events[0].shape, (1, 2))
    assert np.isclose(
            event_rational_1(res.t_events[0][0], res.y_events[0][0]), 0)

    # Also test that termination by event doesn't break interpolants.
    tc = np.linspace(res.t[0], res.t[-1])
    yc_true = sol_rational(tc)
    yc = res.sol(tc)
    e = compute_error(yc, yc_true, 1e-3, 1e-6)
    assert_(np.all(e < 5))

    # Test that the y_event matches solution
    assert np.allclose(
            sol_rational(res.t_events[0][0]), res.y_events[0][0], rtol=1e-3,
            atol=1e-6)

    # Test in backward direction.
    event_rational_1.direction = 0
    event_rational_2.direction = 0

    res = solve_ivp(fun_rational, [8, 5], [4/9, 20/81], method=method,
                    events=(event_rational_1, event_rational_2))
    assert_equal(res.status, 0)
    assert_equal(res.t_events[0].size, 1)
    assert_equal(res.t_events[1].size, 1)
    assert_(5.3 < res.t_events[0][0] < 5.7)
    assert_(7.3 < res.t_events[1][0] < 7.7)

    assert_equal(res.y_events[0].shape, (1, 2))
    assert_equal(res.y_events[1].shape, (1, 2))
    assert np.isclose(
            event_rational_1(res.t_events[0][0], res.y_events[0][0]), 0)
    assert np.isclose(
            event_rational_2(res.t_events[1][0], res.y_events[1][0]), 0)

    event_rational_1.direction = -1
    event_rational_2.direction = -1
    res = solve_ivp(fun_rational, [8, 5], [4/9, 20/81], method=method,
                    events=(event_rational_1, event_rational_2))
    assert_equal(res.status, 0)
    assert_equal(res.t_events[0].size, 1)
    assert_equal(res.t_events[1].size, 0)
    assert_(5.3 < res.t_events[0][0] < 5.7)

    assert_equal(res.y_events[0].shape, (1, 2))
    assert_equal(res.y_events[1].shape, (0,))
    assert np.isclose(
            event_rational_1(res.t_events[0][0], res.y_events[0][0]), 0)

    event_rational_1.direction = 1
    event_rational_2.direction = 1
    res = solve_ivp(fun_rational, [8, 5], [4/9, 20/81], method=method,
                    events=(event_rational_1, event_rational_2))
    assert_equal(res.status, 0)
    assert_equal(res.t_events[0].size, 0)
    assert_equal(res.t_events[1].size, 1)
    assert_(7.3 < res.t_events[1][0] < 7.7)

    assert_equal(res.y_events[0].shape, (0,))
    assert_equal(res.y_events[1].shape, (1, 2))
    assert np.isclose(
            event_rational_2(res.t_events[1][0], res.y_events[1][0]), 0)

    event_rational_1.direction = 0
    event_rational_2.direction = 0

    res = solve_ivp(fun_rational, [8, 5], [4/9, 20/81], method=method,
                    events=(event_rational_1, event_rational_2,
                            event_rational_3), dense_output=True)
    assert_equal(res.status, 1)
    assert_equal(res.t_events[0].size, 0)
    assert_equal(res.t_events[1].size, 1)
    assert_equal(res.t_events[2].size, 1)
    assert_(7.3 < res.t_events[1][0] < 7.7)
    assert_(7.3 < res.t_events[2][0] < 7.5)

    assert_equal(res.y_events[0].shape, (0,))
    assert_equal(res.y_events[1].shape, (1, 2))
    assert_equal(res.y_events[2].shape, (1, 2))
    assert np.isclose(
            event_rational_2(res.t_events[1][0], res.y_events[1][0]), 0)
    assert np.isclose(
            event_rational_3(res.t_events[2][0], res.y_events[2][0]), 0)

    # Also test that termination by event doesn't break interpolants.
    tc = np.linspace(res.t[-1], res.t[0])
    yc_true = sol_rational(tc)
    yc = res.sol(tc)
    e = compute_error(yc, yc_true, 1e-3, 1e-6)
    assert_(np.all(e < 5))


def test_max_step():
    rtol = 1e-3
    atol = 1e-6
    y0 = [1/3, 2/9]
    method = SABM
    for t_span in ([5, 9], [5, 1]):
        res = solve_ivp(fun_rational, t_span, y0, rtol=rtol,
                        max_step=0.5, atol=atol, method=method,
                        dense_output=True)
        assert_equal(res.t[0], t_span[0])
        assert_equal(res.t[-1], t_span[-1])
        assert_(np.all(np.abs(np.diff(res.t)) <= 0.5 + 1e-15))
        assert_(res.t_events is None)
        assert_(res.success)
        assert_equal(res.status, 0)

        y_true = sol_rational(res.t)
        e = compute_error(res.y, y_true, rtol, atol)
        assert_(np.all(e < 5))

        tc = np.linspace(*t_span)
        yc_true = sol_rational(tc)
        yc = res.sol(tc)

        e = compute_error(yc, yc_true, rtol, atol)
        assert_(np.all(e < 5))

        # See comment in test_integration.
        assert_allclose(res.sol(res.t), res.y, rtol=rtol, atol=atol)

        assert_raises(ValueError, method, fun_rational, t_span[0], y0,
                      t_span[1], max_step=-1)

        solver = method(fun_rational, t_span[0], y0, t_span[1],
                        rtol=rtol, atol=atol, max_step=1e-20)
        message = solver.step()

        assert_equal(solver.status, 'failed')
        assert_("step size is less" in message)
        assert_raises(RuntimeError, solver.step)


def test_first_step():
    rtol = 1e-3
    atol = 1e-4
    y0 = [1/3, 2/9]
    first_step = 0.1
    method = SABM
    for t_span in ([5, 9], [5, 1]):
        res = solve_ivp(fun_rational, t_span, y0, rtol=rtol,
                        max_step=0.5, atol=atol, method=method,
                        dense_output=True, first_step=first_step)

        assert_equal(res.t[0], t_span[0])
        assert_equal(res.t[-1], t_span[-1])
        assert_allclose(first_step, np.abs(res.t[1] - 5))
        assert_(res.t_events is None)
        assert_(res.success)
        assert_equal(res.status, 0)

        y_true = sol_rational(res.t)
        e = compute_error(res.y, y_true, rtol, atol)
        assert_(np.all(e < 5))

        tc = np.linspace(*t_span)
        yc_true = sol_rational(tc)
        yc = res.sol(tc)

        e = compute_error(yc, yc_true, rtol, atol)
        assert_(np.all(e < 5))

        assert_allclose(res.sol(res.t), res.y, rtol=rtol, atol=atol)

        assert_raises(ValueError, method, fun_rational, t_span[0], y0,
                      t_span[1], first_step=-1)
        assert_raises(ValueError, method, fun_rational, t_span[0], y0,
                      t_span[1], first_step=5)


def test_t_eval():
    rtol = 1e-3
    atol = 1e-6
    y0 = [1/3, 2/9]
    method = 'SABM'
    for t_span in ([5, 9], [5, 1]):
        t_eval = np.linspace(t_span[0], t_span[1], 10)
        res = solve_ivp(fun_rational, t_span, y0, rtol=rtol, atol=atol,
                        method=method, t_eval=t_eval)
        assert_equal(res.t, t_eval)
        assert_(res.t_events is None)
        assert_(res.success)
        assert_equal(res.status, 0)

        y_true = sol_rational(res.t)
        e = compute_error(res.y, y_true, rtol, atol)
        assert_(np.all(e < 5))

    t_eval = [5, 5.01, 7, 8, 8.01, 9]
    res = solve_ivp(fun_rational, [5, 9], y0, rtol=rtol, atol=atol,
                    method=method, t_eval=t_eval)
    assert_equal(res.t, t_eval)
    assert_(res.t_events is None)
    assert_(res.success)
    assert_equal(res.status, 0)

    y_true = sol_rational(res.t)
    e = compute_error(res.y, y_true, rtol, atol)
    assert_(np.all(e < 5))

    t_eval = [5, 4.99, 3, 1.5, 1.1, 1.01, 1]
    res = solve_ivp(fun_rational, [5, 1], y0, rtol=rtol, atol=atol,
                    method=method, t_eval=t_eval)
    assert_equal(res.t, t_eval)
    assert_(res.t_events is None)
    assert_(res.success)
    assert_equal(res.status, 0)

    t_eval = [5.01, 7, 8, 8.01]
    res = solve_ivp(fun_rational, [5, 9], y0, rtol=rtol, atol=atol,
                    method=method, t_eval=t_eval)
    assert_equal(res.t, t_eval)
    assert_(res.t_events is None)
    assert_(res.success)
    assert_equal(res.status, 0)

    y_true = sol_rational(res.t)
    e = compute_error(res.y, y_true, rtol, atol)
    assert_(np.all(e < 5))

    t_eval = [4.99, 3, 1.5, 1.1, 1.01]
    res = solve_ivp(fun_rational, [5, 1], y0, rtol=rtol, atol=atol,
                    method=method, t_eval=t_eval)
    assert_equal(res.t, t_eval)
    assert_(res.t_events is None)
    assert_(res.success)
    assert_equal(res.status, 0)

    t_eval = [4, 6]
    assert_raises(ValueError, solve_ivp, fun_rational, [5, 9], y0,
                  rtol=rtol, atol=atol, t_eval=t_eval)


def test_t_eval_dense_output():
    rtol = 1e-3
    atol = 1e-6
    y0 = [1/3, 2/9]
    method = 'SABM'
    t_span = [5, 9]
    t_eval = np.linspace(t_span[0], t_span[1], 10)
    res = solve_ivp(fun_rational, t_span, y0, rtol=rtol, atol=atol,
                    method=method, t_eval=t_eval)
    res_d = solve_ivp(fun_rational, t_span, y0, rtol=rtol, atol=atol,
                      method=method, t_eval=t_eval, dense_output=True)
    assert_equal(res.t, t_eval)
    assert_(res.t_events is None)
    assert_(res.success)
    assert_equal(res.status, 0)

    assert_equal(res.t, res_d.t)
    assert_equal(res.y, res_d.y)
    assert_(res_d.t_events is None)
    assert_(res_d.success)
    assert_equal(res_d.status, 0)

    # if t and y are equal only test values for one case
    y_true = sol_rational(res.t)
    e = compute_error(res.y, y_true, rtol, atol)
    assert_(np.all(e < 5))


def test_no_integration():
    method = 'SABM'
    sol = solve_ivp(lambda t, y: -y, [4, 4], [2, 3],
                    method=method, dense_output=True)
    assert_equal(sol.sol(4), [2, 3])
    assert_equal(sol.sol([4, 5, 6]), [[2, 2, 2], [3, 3, 3]])


def test_no_integration_class():
    method = SABM
    solver = method(lambda t, y: -y, 0.0, [10.0, 0.0], 0.0)
    solver.step()
    assert_equal(solver.status, 'finished')
    sol = solver.dense_output()
    assert_equal(sol(0.0), [10.0, 0.0])
    assert_equal(sol([0, 1, 2]), [[10, 10, 10], [0, 0, 0]])

    solver = method(lambda t, y: -y, 0.0, [], np.inf)
    solver.step()
    assert_equal(solver.status, 'finished')
    sol = solver.dense_output()
    assert_equal(sol(100.0), [])
    assert_equal(sol([0, 1, 2]), np.empty((0, 3)))


def test_empty():
    def fun(t, y):
        return np.zeros((0,))

    y0 = np.zeros((0,))
    method = 'SABM'

    sol = assert_no_warnings(solve_ivp, fun, [0, 10], y0,
                             method=method, dense_output=True)
    assert_equal(sol.sol(10), np.zeros((0,)))
    assert_equal(sol.sol([1, 2, 3]), np.zeros((0, 3)))

    sol = assert_no_warnings(solve_ivp, fun, [0, np.inf], y0,
                             method=method, dense_output=True)
    assert_equal(sol.sol(10), np.zeros((0,)))
    assert_equal(sol.sol([1, 2, 3]), np.zeros((0, 3)))


def test_ConstantDenseOutput():
    sol = ConstantDenseOutput(0, 1, np.array([1, 2]))
    assert_allclose(sol(1.5), [1, 2])
    assert_allclose(sol([1, 1.5, 2]), [[1, 1, 1], [2, 2, 2]])

    sol = ConstantDenseOutput(0, 1, np.array([]))
    assert_allclose(sol(1.5), np.empty(0))
    assert_allclose(sol([1, 1.5, 2]), np.empty((0, 3)))


def test_classes():
    y0 = [1 / 3, 2 / 9]
    cls = SABM
    solver = cls(fun_rational, 5, y0, np.inf)
    assert_equal(solver.n, 2)
    assert_equal(solver.status, 'running')
    assert_equal(solver.t_bound, np.inf)
    assert_equal(solver.direction, 1)
    assert_equal(solver.t, 5)
    assert_equal(solver.y, y0)
    assert_(solver.step_size is None)

    assert_(solver.nfev > 0)
    assert_(solver.njev >= 0)
    assert_equal(solver.nlu, 0)

    assert_raises(RuntimeError, solver.dense_output)

    message = solver.step()
    assert_equal(solver.status, 'running')
    assert_equal(message, None)
    assert_equal(solver.n, 2)
    assert_equal(solver.t_bound, np.inf)
    assert_equal(solver.direction, 1)
    assert_(solver.t > 5)
    assert_(not np.all(np.equal(solver.y, y0)))
    assert_(solver.step_size > 0)
    assert_(solver.nfev > 0)
    assert_(solver.njev >= 0)
    assert_(solver.nlu >= 0)
    sol = solver.dense_output()
    assert_allclose(sol(5), y0, rtol=1e-8, atol=0)


def test_OdeSolution():
    ts = np.array([0, 2, 5], dtype=float)
    s1 = ConstantDenseOutput(ts[0], ts[1], np.array([-1]))
    s2 = ConstantDenseOutput(ts[1], ts[2], np.array([1]))

    sol = OdeSolution(ts, [s1, s2])

    assert_equal(sol(-1), [-1])
    assert_equal(sol(1), [-1])
    assert_equal(sol(2), [-1])
    assert_equal(sol(3), [1])
    assert_equal(sol(5), [1])
    assert_equal(sol(6), [1])

    assert_equal(sol([0, 6, -2, 1.5, 4.5, 2.5, 5, 5.5, 2]),
                 np.array([[-1, 1, -1, -1, 1, 1, 1, 1, -1]]))

    ts = np.array([10, 4, -3])
    s1 = ConstantDenseOutput(ts[0], ts[1], np.array([-1]))
    s2 = ConstantDenseOutput(ts[1], ts[2], np.array([1]))

    sol = OdeSolution(ts, [s1, s2])
    assert_equal(sol(11), [-1])
    assert_equal(sol(10), [-1])
    assert_equal(sol(5), [-1])
    assert_equal(sol(4), [-1])
    assert_equal(sol(0), [1])
    assert_equal(sol(-3), [1])
    assert_equal(sol(-4), [1])

    assert_equal(sol([12, -5, 10, -3, 6, 1, 4]),
                 np.array([[-1, 1, -1, 1, -1, 1, -1]]))

    ts = np.array([1, 1])
    s = ConstantDenseOutput(1, 1, np.array([10]))
    sol = OdeSolution(ts, [s])
    assert_equal(sol(0), [10])
    assert_equal(sol(1), [10])
    assert_equal(sol(2), [10])

    assert_equal(sol([2, 1, 0]), np.array([[10, 10, 10]]))


@pytest.mark.parametrize('method', ['SABM'])
def test_integration_zero_rhs(method):
    result = solve_ivp(fun_zero, [0, 10], np.ones(3), method=method)
    assert_(result.success)
    assert_equal(result.status, 0)
    assert_allclose(result.y, 1.0, rtol=1e-15)
