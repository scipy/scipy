import numpy as np
from numpy.testing import (assert_, assert_array_equal, assert_allclose,
                           run_module_suite, assert_raises, assert_equal)
from scipy.sparse import coo_matrix
from scipy.integrate._bvp import (modify_mesh, compute_fun_jac,
                                  compute_bc_jac, compute_jac_indices,
                                  construct_global_jac, solve_bvp)


def measles_fun(x, y):
    beta = 1575 * (1 + np.cos(2 * np.pi * x))
    return np.vstack((
        0.02 - beta * y[0] * y[2],
        beta * y[0] * y[2] - y[1] / 0.0279,
        y[1] / 0.0279 - y[2] / 0.01
    ))


def measles_df_dy(x, y):
    n, m = y.shape
    beta = 1575 * (1 + np.cos(2 * np.pi * x))
    J = np.empty((n, n, m))
    J[0, 0] = -beta * y[2]
    J[0, 1] = 0
    J[0, 2] = -beta * y[0]
    J[1, 0] = beta * y[2]
    J[1, 1] = -1 / 0.0279
    J[1, 2] = beta * y[0]
    J[2, 0] = 0
    J[2, 1] = 1 / 0.0279
    J[2, 2] = -1 / 0.01
    return J


def exp_fun(x, y):
    return np.vstack((y[1], y[0]))


def exp_bc(ya, yb):
    return np.hstack((ya[0] - 1, yb[0]))


def exp_sol(x):
    return (np.exp(-x) - np.exp(x - 2)) / (1 - np.exp(-2))


def sl_fun(x, y, p):
    return np.vstack((y[1], -p[0]**2 * y[0]))


def sl_fun_jac(x, y, p):
    n, m = y.shape
    Jy = np.empty((n, 2, m))
    Jy[0, 0] = 0
    Jy[0, 1] = 1
    Jy[1, 0] = -p[0]**2
    Jy[1, 1] = 0

    Jp = np.empty((n, 1, m))
    Jp[0, 0] = 0
    Jp[1, 0] = -2 * p[0] * y[0]

    return Jy, Jp


def sl_bc(ya, yb, p):
    return np.hstack((ya[0], yb[0], ya[1] - p[0]))


def sl_bc_jac(ya, yb, p):
    Jya = np.zeros((3, 2))
    Jyb = np.zeros((3, 2))
    Jp = np.zeros((3, 1))

    Jya[0, 0] = 1
    Jya[2, 1] = 1

    Jyb[1, 0] = 1

    Jp[2, 0] = -1

    return Jya, Jyb, Jp


def sl_sol(x, p):
    return np.sin(p[0] * x)


def test_modify_mesh():
    x = np.array([0, 1, 3, 9], dtype=float)
    x_new = modify_mesh(x, np.array([0]), np.array([2]))
    assert_array_equal(x_new, np.array([0, 0.5, 1, 3, 5, 7, 9]))


def test_compute_fun_jac():
    x = np.linspace(0, 1, 5)
    y = np.empty((3, x.shape[0]))
    y[0] = 0.01
    y[1] = 0.02
    y[2] = 0.03
    p = np.array([])
    Jy, Jp = compute_fun_jac(lambda x, y, p: measles_fun(x, y), x, y, p)
    Jy_an = measles_df_dy(x, y)
    assert_allclose(Jy, Jy_an)
    assert_(Jp is None)

    x = np.linspace(0, np.pi, 5)
    y = np.empty((2, x.shape[0]))
    y[0] = np.sin(x)
    y[1] = np.cos(x)
    p = np.array([1.0])
    Jy, Jp = compute_fun_jac(sl_fun, x, y, p)
    Jy_an, Jp_an = sl_fun_jac(x, y, p)
    assert_allclose(Jy, Jy_an)
    assert_allclose(Jp, Jp_an)


def test_compute_jac_indices():
    n = 2
    m = 4
    k = 2
    i, j = compute_jac_indices(n, m, k)
    s = coo_matrix((np.ones_like(i), (i, j))).toarray()
    s_true = np.array([
        [1, 1, 1, 1, 0, 0, 0, 0, 1, 1],
        [1, 1, 1, 1, 0, 0, 0, 0, 1, 1],
        [0, 0, 1, 1, 1, 1, 0, 0, 1, 1],
        [0, 0, 1, 1, 1, 1, 0, 0, 1, 1],
        [0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
        [0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
        [1, 1, 0, 0, 0, 0, 1, 1, 1, 1],
        [1, 1, 0, 0, 0, 0, 1, 1, 1, 1],
        [1, 1, 0, 0, 0, 0, 1, 1, 1, 1],
        [1, 1, 0, 0, 0, 0, 1, 1, 1, 1],
    ])
    assert_array_equal(s, s_true)


def test_compute_bc_jac():
    ya = np.array([0.0, 1])
    yb = np.array([0.0, -1])
    p = np.array([0.5])
    Jya, Jyb, Jp = compute_bc_jac(sl_bc, ya, yb, p)
    Jya_an, Jyb_an, Jp_an = sl_bc_jac(ya, yb, p)
    assert_allclose(Jya, Jya_an)
    assert_allclose(Jyb, Jyb_an)
    assert_allclose(Jp, Jp_an)


def test_compute_global_jac():
    n = 2
    m = 5
    k = 1
    i_jac, j_jac = compute_jac_indices(2, 5, 1)
    x = np.linspace(0, 1, 5)
    h = np.diff(x)
    y = np.vstack((np.sin(np.pi * x), np.pi * np.cos(np.pi * x)))
    p = np.array([3.0])

    f = sl_fun(x, y, p)

    x_middle = x[:-1] + 0.5 * h
    y_middle = 0.5 * (y[:, :-1] + y[:, 1:]) - h/8 * (f[:, 1:] - f[:, :-1])

    df_dy, df_dp = sl_fun_jac(x, y, p)
    df_dy_middle, df_dp_middle = sl_fun_jac(x_middle, y_middle, p)
    dbc_dya, dbc_dyb, dbc_dp = sl_bc_jac(y[:, 0], y[:, -1], p)

    J = construct_global_jac(n, m, k, i_jac, j_jac, h, df_dy, df_dy_middle,
                             df_dp, df_dp_middle, dbc_dya, dbc_dyb, dbc_dp)
    J = J.toarray()

    def J_block(h, p):
        return np.array([
            [h**2*p**2/12 - 1, -0.5*h, -h**2*p**2/12 + 1, -0.5*h],
            [0.5*h*p**2, h**2*p**2/12 - 1, 0.5*h*p**2, 1 - h**2*p**2/12]
        ])

    J_true = np.zeros((m * n + k, m * n + k))
    for i in range(m - 1):
        J_true[i * n: (i + 1) * n, i * n: (i + 2) * n] = J_block(h[i], p)

    J_true[:(m - 1) * n:2, -1] = p * h**2/6 * (y[0, :-1] - y[0, 1:])
    J_true[1:(m - 1) * n:2, -1] = p * (h * (y[0, :-1] + y[0, 1:]) +
                                       h**2/6 * (y[1, :-1] - y[1, 1:]))

    J_true[8, 0] = 1
    J_true[9, 8] = 1
    J_true[10, 1] = 1
    J_true[10, 10] = -1

    assert_allclose(J, J_true, rtol=1e-10)

    df_dy, df_dp = compute_fun_jac(sl_fun, x, y, p)
    df_dy_middle, df_dp_middle = compute_fun_jac(sl_fun, x_middle, y_middle, p)
    dbc_dya, dbc_dyb, dbc_dp = compute_bc_jac(sl_bc, y[:, 0], y[:, -1], p)
    J = construct_global_jac(n, m, k, i_jac, j_jac, h, df_dy, df_dy_middle,
                             df_dp, df_dp_middle, dbc_dya, dbc_dyb, dbc_dp)
    J = J.toarray()
    assert_allclose(J, J_true, rtol=1e-8, atol=1e-9)


def test_parameter_validation():
    x = [0, 1, 0.5]
    y = np.zeros((2, 3))
    assert_raises(ValueError, solve_bvp, exp_fun, exp_bc, x, y)

    x = np.linspace(0, 1, 5)
    y = np.zeros((2, 4))
    assert_raises(ValueError, solve_bvp, exp_fun, exp_bc, x, y)

    fun = lambda x, y, p: exp_fun(x, y)
    bc = lambda ya, yb, p: exp_bc(ya, yb)

    y = np.zeros((2, x.shape[0]))
    assert_raises(ValueError, solve_bvp, fun, bc, x, y, p=[1])

    def wrong_shape_fun(x, y):
        return np.zeros(3)

    assert_raises(ValueError, solve_bvp, wrong_shape_fun, bc, x, y)


def test_no_params():
    x = np.linspace(0, 1, 5)
    x_test = np.linspace(0, 1, 100)
    y = np.zeros((2, x.shape[0]))
    sol = solve_bvp(exp_fun, exp_bc, x, y)

    assert_equal(sol.status, 0)
    assert_(sol.success)

    sol_test = sol.sol(x_test)

    assert_allclose(sol_test[0], exp_sol(x_test), atol=1e-5)

    f_test = exp_fun(x_test, sol_test)
    res = sol.sol(x_test, 1) - f_test
    rel_res = res / (1 + np.abs(f_test))
    norm_res = np.sum(rel_res**2, axis=0)**0.5
    assert_(np.all(norm_res < 1e-3))

    assert_(np.all(sol.res < 1e-3))
    assert_allclose(sol.sol(sol.x), sol.y, rtol=1e-10, atol=1e-10)
    assert_allclose(sol.sol(sol.x, 1), sol.f, rtol=1e-10, atol=1e-10)


def test_with_params():
    x = np.linspace(0, np.pi, 5)
    x_test = np.linspace(0, np.pi, 100)
    y = np.ones((2, x.shape[0]))
    sol = solve_bvp(sl_fun, sl_bc, x, y, p=[0.5])

    assert_equal(sol.status, 0)
    assert_(sol.success)

    assert_allclose(sol.p, [1.0], rtol=1e-4)

    sol_test = sol.sol(x_test)

    assert_allclose(sol_test[0], sl_sol(x_test, [1]), rtol=1e-4, atol=1e-4)

    f_test = sl_fun(x_test, sol_test, [1])
    res = sol.sol(x_test, 1) - f_test
    rel_res = res / (1 + np.abs(f_test))
    norm_res = np.sum(rel_res ** 2, axis=0) ** 0.5
    assert_(np.all(norm_res < 1e-3))

    assert_(np.all(sol.res < 1e-3))
    assert_allclose(sol.sol(sol.x), sol.y, rtol=1e-10, atol=1e-10)
    assert_allclose(sol.sol(sol.x, 1), sol.f, rtol=1e-10, atol=1e-10)


if __name__ == '__main__':
    run_module_suite()
