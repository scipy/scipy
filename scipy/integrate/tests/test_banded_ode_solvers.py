
from __future__ import division, print_function, absolute_import

import itertools
import numpy as np
from numpy.testing import run_module_suite, assert_allclose
from scipy.integrate import ode


def _band_count(a):
    """Returns ml and mu, the lower and upper band sizes of a."""
    nrows, ncols = a.shape
    ml = 0
    for k in range(-nrows+1, 0):
        if np.diag(a, k).any():
            ml = -k
            break
    mu = 0
    for k in range(nrows-1, 0, -1):
        if np.diag(a, k).any():
            mu = k
            break
    return ml, mu


def _linear_func(t, y, a):
    """Linear system dy/dt = a * y"""
    return a.dot(y)


def _linear_jac(t, y, a):
    """Jacobian of a * y is a."""
    return a


def _linear_banded_jac(t, y, a):
    """Banded Jacobian."""
    ml, mu = _band_count(a)
    bjac = []
    for k in range(mu, 0, -1):
        bjac.append(np.r_[[0] * k, np.diag(a, k)])
    bjac.append(np.diag(a))
    for k in range(-1, -ml-1, -1):
        bjac.append(np.r_[np.diag(a, k), [0] * (-k)])
    return bjac


def _solve_linear_sys(a, y0, tend=1, dt=0.1,
                      solver=None, method='bdf', use_jac=True,
                      with_jacobian=False, banded=False):
    """Use scipy.integrate.ode to solve a linear system of ODEs."""
    if banded:
        lband, uband = _band_count(a)
    else:
        lband = None
        uband = None

    if use_jac:
        if banded:
            r = ode(_linear_func, _linear_banded_jac)
        else:
            r = ode(_linear_func, _linear_jac)
    else:
        r = ode(_linear_func)

    if solver is None:
        if np.iscomplexobj(a):
            solver = "zvode"
        else:
            solver = "vode"

    r.set_integrator(solver,
                     with_jacobian=with_jacobian,  # Is this redundant?
                     method=method,
                     lband=lband, uband=uband,
                     rtol=1e-9, atol=1e-10,
                     )
    t0 = 0
    r.set_initial_value(y0, t0)
    r.set_f_params(a)
    r.set_jac_params(a)

    t = [t0]
    y = [y0]
    while r.successful() and r.t < tend:
        r.integrate(r.t + dt)
        t.append(r.t)
        y.append(r.y)

    t = np.array(t)
    y = np.array(y)
    return t, y


def _analytical_solution(a, y0, t):
    """
    Analytical solution to the linear differential equations dy/dt = a*y.

    The solution is only valid if `a` is diagonalizable.

    Returns a 2-d array with shape (len(t), len(y0)).
    """
    lam, v = np.linalg.eig(a)
    c = np.linalg.solve(v, y0)
    e = c * np.exp(lam * t.reshape(-1, 1))
    sol = e.dot(v.T)
    return sol


def test_banded_ode_solvers():
    # Test the "lsoda", "vode" and "zvode" solvers of the `ode` class
    # with a system that has a banded Jacobian matrix.

    # lband = 2, uband = 1:
    a_real = np.array([[-0.6,  0.1,    0,    0,    0],
                       [ 0.2, -0.5,  0.9,    0,    0],
                       [ 0.1,  0.1, -0.4,  0.1,    0],
                       [   0,  0.3, -0.1, -0.9, -0.3],
                       [   0,    0,  0.1,  0.1, -0.7]])
    # lband = 0, uband = 1:
    a_real_upper = np.triu(a_real)
    # lband = 2, uband = 0:
    a_real_lower = np.tril(a_real)
    # lband = 0, uband = 0
    a_real_diag = np.triu(a_real_lower)

    # complex, lband = 2, uband = 1:
    a_complex = a_real - 0.5j * a_real
    # complex, lband = 0, uband = 0
    a_complex_diag = np.diag(np.diag(a_complex))

    tend = 1.0
    dt = 0.25
    t_exact = np.arange(0, tend + 0.5*dt, dt)

    # Real systems
    for a in [a_real, a_real_upper, a_real_lower, a_real_diag]:
        y0 = np.arange(1, a.shape[0]+1)
        y_exact = _analytical_solution(a, y0, t_exact)
        p = [['vode', 'lsoda'],  # solver
             ['bdf', 'adams'],   # method
             [False, True],      # use_jac
             [False, True],      # with_jacobian
             [False, True],      # banded
            ]
        for solver, method, use_jac, with_jacobian, banded in \
                                                 itertools.product(*p):
            t, y = _solve_linear_sys(a, y0, tend=tend, dt=dt,
                                     solver=solver,
                                     method=method,
                                     use_jac=use_jac,
                                     with_jacobian=with_jacobian,
                                     banded=banded)
            yield assert_allclose, t, t_exact
            yield assert_allclose, y, y_exact

    # Complex systems (solver is "zvode")
    for a in [a_complex, a_complex_diag]:
        y0 = np.arange(1, a.shape[0]+1) + 1j
        y_exact = _analytical_solution(a, y0, t_exact)
        p = [['bdf', 'adams'],  # method
             [False, True],     # use_jac
             [False, True],     # with_jacobian
             [False, True],     # banded
            ]
        for method, use_jac, with_jacobian, banded in itertools.product(*p):
            t, y = _solve_linear_sys(a, y0, tend=tend, dt=dt,
                                     solver="zvode",
                                     method=method,
                                     use_jac=use_jac,
                                     with_jacobian=with_jacobian,
                                     banded=banded)
            yield assert_allclose, t, t_exact
            yield assert_allclose, y, y_exact


if __name__ == "__main__":
    run_module_suite()
