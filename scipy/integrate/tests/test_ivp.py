from __future__ import division, print_function, absolute_import
from numpy.testing import (assert_, assert_allclose, run_module_suite,
                           assert_equal, assert_raises, assert_no_warnings)
import numpy as np
from scipy.integrate import solve_ivp, RK23, RK45, Radau, BDF
from scipy.integrate import DenseOutput, OdeSolution
from scipy.integrate._py.common import num_jac


def fun_rational(t, y):
    return np.array([y[1] / t,
                     y[1] * (y[0] + 2 * y[1] - 1) / (t * (y[0] - 1))])


def jac_rational(t, y):
    return np.array([
        [0, 1 / t],
        [-2 * y[1] ** 2 / (t * (y[0] - 1) ** 2),
         (y[0] + 4 * y[1] - 1) / (t * (y[0] - 1))]
    ])


def sol_rational(t):
    return np.vstack((t / (t + 10), 10 * t / (t + 10) ** 2))


def event_rational_1(t, y):
    return y[0] - y[1] ** 0.7


def event_rational_2(t, y):
    return y[1] ** 0.6 - y[0]


def event_rational_3(t, y):
    return t - 7.4


def compute_error(y, y_true, rtol, atol):
    e = (y - y_true) / (atol + rtol * np.abs(y_true))
    return np.sqrt(np.sum(e**2, axis=0) / e.shape[0])


def test_integration():
    rtol = 1e-3
    atol = 1e-6
    y0 = [1/3, 2/9]
    for method in ['RK23', 'RK45', 'Radau', 'BDF']:
        for t_span in ([5, 9], [5, 1]):
            res = solve_ivp(fun_rational, t_span, y0, rtol=rtol,
                            atol=atol, method=method, dense_output=True)
            assert_equal(res.t[0], t_span[0])
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

            assert_allclose(res.sol(res.t), res.y, rtol=1e-15, atol=1e-15)


def test_events():
    event_rational_3.terminal = True

    for method in ['RK23', 'RK45', 'Radau', 'BDF']:
        res = solve_ivp(fun_rational, [5, 8], [1/3, 2/9], method=method,
                        events=(event_rational_1, event_rational_2))
        assert_equal(res.status, 0)
        assert_equal(res.t_events[0].size, 1)
        assert_equal(res.t_events[1].size, 1)
        assert_(5.3 < res.t_events[0][0] < 5.7)
        assert_(7.3 < res.t_events[1][0] < 7.7)

        event_rational_1.direction = 1
        event_rational_2.direction = 1
        res = solve_ivp(fun_rational, [5, 8], [1 / 3, 2 / 9], method=method,
                        events=(event_rational_1, event_rational_2))
        assert_equal(res.status, 0)
        assert_equal(res.t_events[0].size, 1)
        assert_equal(res.t_events[1].size, 0)
        assert_(5.3 < res.t_events[0][0] < 5.7)

        event_rational_1.direction = -1
        event_rational_2.direction = -1
        res = solve_ivp(fun_rational, [5, 8], [1 / 3, 2 / 9], method=method,
                        events=(event_rational_1, event_rational_2))
        assert_equal(res.status, 0)
        assert_equal(res.t_events[0].size, 0)
        assert_equal(res.t_events[1].size, 1)
        assert_(7.3 < res.t_events[1][0] < 7.7)

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

        # Also test that termination by event doesn't break interpolants.
        tc = np.linspace(res.t[0], res.t[-1])
        yc_true = sol_rational(tc)
        yc = res.sol(tc)
        e = compute_error(yc, yc_true, 1e-3, 1e-6)
        assert_(np.all(e < 5))

    # Test in backward direction.
    event_rational_1.direction = 0
    event_rational_2.direction = 0
    for method in ['RK23', 'RK45', 'Radau', 'BDF']:
        res = solve_ivp(fun_rational, [8, 5], [4/9, 20/81], method=method,
                        events=(event_rational_1, event_rational_2))
        assert_equal(res.status, 0)
        assert_equal(res.t_events[0].size, 1)
        assert_equal(res.t_events[1].size, 1)
        assert_(5.3 < res.t_events[0][0] < 5.7)
        assert_(7.3 < res.t_events[1][0] < 7.7)

        event_rational_1.direction = -1
        event_rational_2.direction = -1
        res = solve_ivp(fun_rational, [8, 5], [4/9, 20/81], method=method,
                        events=(event_rational_1, event_rational_2))
        assert_equal(res.status, 0)
        assert_equal(res.t_events[0].size, 1)
        assert_equal(res.t_events[1].size, 0)
        assert_(5.3 < res.t_events[0][0] < 5.7)

        event_rational_1.direction = 1
        event_rational_2.direction = 1
        res = solve_ivp(fun_rational, [8, 5], [4/9, 20/81], method=method,
                        events=(event_rational_1, event_rational_2))
        assert_equal(res.status, 0)
        assert_equal(res.t_events[0].size, 0)
        assert_equal(res.t_events[1].size, 1)
        assert_(7.3 < res.t_events[1][0] < 7.7)

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
    for method in ['RK23', 'RK45', 'Radau', 'BDF']:
        for t_span in ([5, 9], [5, 1]):
            res = solve_ivp(fun_rational, t_span, y0, rtol=rtol,
                            max_step=0.5, atol=atol, method=method,
                            dense_output=True)
            assert_equal(res.t[0], t_span[0])
            assert_equal(res.t[-1], t_span[-1])
            assert_(np.all(np.abs(np.diff(res.t)) <= 0.5))
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

            assert_allclose(res.sol(res.t), res.y, rtol=1e-15, atol=1e-15)


def test_t_eval():
    rtol = 1e-3
    atol = 1e-6
    y0 = [1/3, 2/9]
    for t_span in ([5, 9], [5, 1]):
        t_eval = np.linspace(*t_span, 10)
        res = solve_ivp(fun_rational, t_span, y0, rtol=rtol, atol=atol,
                        t_eval=t_eval)
        assert_equal(res.t, t_eval)
        assert_(res.t_events is None)
        assert_(res.success)
        assert_equal(res.status, 0)

        y_true = sol_rational(res.t)
        e = compute_error(res.y, y_true, rtol, atol)
        assert_(np.all(e < 5))

    t_eval = [5, 5.01, 7, 8, 8.01, 9]
    res = solve_ivp(fun_rational, [5, 9], y0, rtol=rtol, atol=atol,
                    t_eval=t_eval)
    assert_equal(res.t, t_eval)
    assert_(res.t_events is None)
    assert_(res.success)
    assert_equal(res.status, 0)

    y_true = sol_rational(res.t)
    e = compute_error(res.y, y_true, rtol, atol)
    assert_(np.all(e < 5))

    t_eval = [5, 4.99, 3, 1.5, 1.1, 1.01, 1]
    res = solve_ivp(fun_rational, [5, 1], y0, rtol=rtol, atol=atol,
                    t_eval=t_eval)
    assert_equal(res.t, t_eval)
    assert_(res.t_events is None)
    assert_(res.success)
    assert_equal(res.status, 0)

    t_eval = [5.01, 7, 8, 8.01]
    res = solve_ivp(fun_rational, [5, 9], y0, rtol=rtol, atol=atol,
                    t_eval=t_eval)
    assert_equal(res.t, t_eval)
    assert_(res.t_events is None)
    assert_(res.success)
    assert_equal(res.status, 0)

    y_true = sol_rational(res.t)
    e = compute_error(res.y, y_true, rtol, atol)
    assert_(np.all(e < 5))

    t_eval = [4.99, 3, 1.5, 1.1, 1.01]
    res = solve_ivp(fun_rational, [5, 1], y0, rtol=rtol, atol=atol,
                    t_eval=t_eval)
    assert_equal(res.t, t_eval)
    assert_(res.t_events is None)
    assert_(res.success)
    assert_equal(res.status, 0)

    t_eval = [4, 6]
    assert_raises(ValueError, solve_ivp, fun_rational, [5, 9], y0,
                  rtol=rtol, atol=atol, t_eval=t_eval)


def test_no_integration():
    for method in ['RK23', 'RK45', 'Radau', 'BDF']:
        sol = solve_ivp(lambda t, y: -y, [4, 4], [2, 3], method=method, dense_output=True)
        assert_equal(sol.sol(4), [2, 3])


def test_empty():
    def fun(t, y):
        return np.zeros((0,))

    ic = np.zeros((0,))

    for method in ['RK23', 'RK45', 'Radau', 'BDF']:
        sol = assert_no_warnings(solve_ivp, fun, [0, 10], ic, method=method, dense_output=True)
        assert_equal(sol.sol(10), np.zeros((0,)))

    for method in ['RK23', 'RK45', 'Radau', 'BDF']:
        sol = assert_no_warnings(solve_ivp, fun, [0, np.inf], ic, method=method, dense_output=True)
        assert_equal(sol.sol(10), np.zeros((0,)))


def test_classes():
    y0 = [1 / 3, 2 / 9]
    for cls in [RK23, RK45, Radau, BDF]:
        solver = cls(fun_rational, 5, y0, np.inf)
        assert_equal(solver.n, 2)
        assert_equal(solver.status, 'running')
        assert_equal(solver.t_crit, np.inf)
        assert_equal(solver.direction, 1)
        assert_equal(solver.t, 5)
        assert_equal(solver.y, y0)
        assert_(solver.step_size is None)
        assert_raises(RuntimeError, solver.dense_output)

        message = solver.step()
        assert_equal(solver.status, 'running')
        assert_equal(message, None)
        assert_equal(solver.n, 2)
        assert_equal(solver.t_crit, np.inf)
        assert_equal(solver.direction, 1)
        assert_(solver.t > 5)
        assert_(not np.all(np.equal(solver.y, y0)))
        assert_(solver.step_size > 0)

        assert_raises(ValueError, solver.step, max_step=-1)

        message = solver.step(max_step=1e-20)
        assert_equal(solver.status, 'failed')
        assert_("step size is less" in message)
        assert_raises(RuntimeError, solver.step)


class ConstDenseOutput(DenseOutput):
    def __init__(self, t_old, t, c):
        super(ConstDenseOutput, self).__init__(t_old, t)
        self.c = c

    def _call_impl(self, t):
        if t.ndim == 0:
            return self.c
        else:
            y = np.empty_like(t)
            y.fill(self.c)
            y = np.atleast_2d(y)

        return y


def test_OdeSolution():
    ts = np.array([0, 2, 5], dtype=float)
    s1 = ConstDenseOutput(ts[0], ts[1], -1)
    s2 = ConstDenseOutput(ts[1], ts[2], 1)

    sol = OdeSolution(ts, [s1, s2])

    assert_equal(sol(-1), -1)
    assert_equal(sol(1), -1)
    assert_equal(sol(3), 1)
    assert_equal(sol(6), 1)

    assert_equal(sol([6, -2, 1.5, 4.5, 2.5, 5.5]),
                 np.array([[1, -1, -1, 1, 1, 1]]))


def test_numjac():
    def fun(t, y):
        return np.array([
            -0.04 * y[0] + 1e4 * y[1] * y[2],
            0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] ** 2,
            3e7 * y[1] ** 2
        ])

    def jac(t, y):
        return np.array([
            [-0.04, 1e4 * y[2], 1e4 * y[1]],
            [0.04, -1e4 * y[2] - 6e7 * y[1], -1e4 * y[1]],
            [0, 6e7 * y[1], 0]
        ])

    t = 1
    y = np.array([1, 0, 0])
    J_true = jac(t, y)
    threshold = 1e-5

    J_num, factor = num_jac(fun, t, y, fun(t, y), threshold, None)
    assert_allclose(J_num, J_true, rtol=1e-5, atol=1e-5)

    J_num, factor = num_jac(fun, t, y, fun(t, y), threshold, factor)
    assert_allclose(J_num, J_true, rtol=1e-5, atol=1e-5)


if __name__ == '__main__':
    run_module_suite()
