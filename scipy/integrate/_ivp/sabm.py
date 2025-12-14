from collections import deque
import numpy as np
from scipy.optimize import fsolve
from .base import OdeSolver, DenseOutput
from .common import (validate_max_step, validate_tol,
                     select_initial_step, validate_first_step, warn_extraneous)


MAX_ORDER = 5
MIN_H = 1e-12

# bashforth coefficients
B = [[1, 0, 0, 0, 0, 0],
     [3/2, -1/2, 0, 0, 0, 0],
     [23/12, -16/12, 5/12, 0, 0, 0],
     [55/24, -59/24, 37/24, -9/24, 0, 0],
     [1901/720, -2774/720, 2616/720, -1274/720, 251/720, 0],
     [4277/1440, -7923/1440, 9982/1440, -7298/1440, 2877/1440, -475/1440]]

# moulton coefficients
M = [[1, 0, 0, 0, 0, 0],
     [1/2, 1/2, 0, 0, 0, 0],
     [5/12, 2/3, -1/12, 0, 0, 0],
     [3/8, 19/24, -5/24, 1/24, 0, 0],
     [251/720, 646/720, -264/720, 106/720, -19/720, 0],
     [95/288, 1427/1440, -133/240, 241/720, 173/1440, 3/160]]


class SABM(OdeSolver):
    """Semi-explicit and semi-implicit method based on
    Adams-Bashforth-Moulton formulas [1]_.

    This method solve differential-algebraic equations (DAE) [2]_ and ODE.
    This is fixed order method with the order that can be choosen between
    1 to 5. The general framework of the SABM algorithm is described in [1]_.
    This class implements a varying step size which depends on the error
    estimation. The strategy used for the error is the difference between
    the predictor and the corrector.

    Can not be applied in the complex domain.

    Note that the nfev is over evaluated. Most of the evaluations are
    done by the constraints (function gun). This may be not efficient if the
    constraints are difficult to solve.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system for the differential variables. There is
        two option for the signature of fun : It can be ``fun(t, y, x)`` to
        solve DAE. Here ``t`` is a scalar, ``y`` is a ndarray with shape (n,)
        and ``x`` is ndarray with shape (m,; then ``fun`` must return
        array_like with shape (n,). Alternatively it can be ``fun(t,y)``
        to solve an ODE.
    t0 : float
        Initial time.
    y0 : array_like, shape (n,) if ODE and shape(n+m,) if DAE
        Initial state.
    t_bound : float
        Boundary time - the integration won't continue beyond it. It also
        determines the direction of the integration.
    gun : callable, optional
        Algebraic equations. Default is ``None`` which means that is an ODE.
        Otherwise, it is the right-hand side of the system for the algebraic
        variables. The calling signature is ``gun(t, y, x)``. Here ``t``
        is a scalar, ``y`` is a ndarray with shape (n,), and ``x`` is a
        ndarray with shape (m,); then ``gun`` must return array_like with
        shape (m,).
    mode : string, optional
        Choose between ``Explicit`` and ``Implicit`` form of the method.
        Default is ``Explicit``.
    order : int, optional
        Choose the order of the method. Default is 2. The order is bounded by
        ``MAX_ORDER`` which is egal to 5.
    num_diff : int, optional
        The number of differential variables. Must be set at 'n' for a DAE.
    max_step : float, optional
        Maximum allowed step size. Default is np.inf, i.e., the step size is
        not bounded and determined solely by the solver.
    rtol : float, optional
        Relative tolerance. This tolerance is use to solve the algebraic
        equations with fsolve. Default is the default relative tolerance
        use by fsolve.
    atol : float, array_like shape (n,), optional
        Absolute tolerance. The step will terminate if the absolute error
        between the predictor and the corrector for each variable is at most
        atol. Default value is 1e-6.
    vectorized : bool, optional
        Whether `fun` is implemented in a vectorized fashion. Default is False.
    first_step : float or None, optional
        Initial step size. Default is ``None`` which means that the algorithm
        should choose.

    Attributes
    ----------
    n : int
        Number of equations.
    status : string
        Current status of the solver: 'running', 'finished' or 'failed'.
    t_bound : float
        Boundary time.
    direction : float
        Integration direction: +1 or -1.
    t : float
        Current time.
    y : ndarray
        Current state.
    t_old : float
        Previous time. None if no steps were made yet.
    step_size : float
        Size of the last successful step. None if no steps were made yet.
    nfev : int
        Number of evaluations of the right-hand side.
    njev : int
        Number of evaluations of the Jacobian.
    nlu : int
        Number of LU decompositions.

    References
    ----------
    .. [1] Tutueva, A; Karimov, T.; Butusov, D. Semi-Implicit and Semi-Explicit
           Adams-Bashforth-Moulton Methods. Mathematics 2020, 8, 780.
    .. [2] E. Hairer, G. Wanner, "Solving Ordinary Differential Equations II:
           Stiff and Differential-Algebraic Problems".
    """
    def __init__(self, fun, t0, y0, t_bound, gun=None, mode="Explicit",
                 order=2, num_diff=None, max_step=np.inf, rtol=1.49012e-08,
                 atol=1e-6, vectorized=False, first_step=None, **extraneous):

        warn_extraneous(extraneous)
        self.ode = gun is None
        y0 = np.array(y0)
        if self.ode:
            f = fun
            num_diff = len(y0)
            f0 = fun(t0, y0)
        else:
            if num_diff is None:
                raise ValueError("num_diff must be set")

            def f(t, x):
                return np.concatenate((fun(t, x[:num_diff], x[num_diff:]),
                                       gun(t, x[:num_diff], x[num_diff:])),
                                      axis=None)
            f0 = fun(t0, y0[:num_diff], y0[num_diff:])

        super().__init__(f, t0, y0, t_bound, vectorized,
                         support_complex=False)
        self.max_step = validate_max_step(max_step)

        if np.asarray(atol).ndim == 0:
            atol = np.ones(num_diff)*atol
        self.rtol, self.atol = validate_tol(rtol, atol, num_diff)

        if first_step is None:
            if self.ode:
                self.h_abs = select_initial_step(fun, self.t, self.y, f0,
                                                 self.direction, 1, self.rtol,
                                                 self.atol)
            else:
                def fun_0(t, y):
                    return fun(t, y, y0[num_diff:])
                self.h_abs = select_initial_step(fun_0, self.t,
                                                 self.y[:num_diff], f0,
                                                 self.direction, 1,
                                                 self.rtol, self.atol)
        else:
            self.h_abs = validate_first_step(first_step, t0, t_bound)

        if self.ode:
            if vectorized:
                self.f = lambda t, y, x: np.hstack(fun(t, y))
            else:
                self.f = lambda t, y, x: np.array(fun(t, y))
        else:
            if vectorized:
                self.f = lambda t, y, x: np.hstack(fun(t, y, x))
            else:
                self.f = lambda t, y, x: np.array(fun(t, y, x))

        if vectorized:
            self.g = lambda t, y, x: np.hstack(gun(t, y, x))
        else:
            self.g = lambda t, y, x: np.array(gun(t, y, x))

        # number of differential equation
        self.num_diff = num_diff

        # order method
        self.p = 1
        self.order = min(order, MAX_ORDER)

        # implicit / explicit
        self.mode = mode

        self.prev_f_y = deque(maxlen=order)
        self.prev_f_y.append(self.f(t0, y0[:num_diff], y0[num_diff:]))
        self.nfev += 1

    def _step_impl(self):
        if self.mode == "Implicit":
            (y_new, x_new, err) = self.Semi_Implicit_Adams_Bashforth_Moulton()
        elif self.mode == "Explicit":
            (y_new, x_new, err) = self.Semi_Explicit_Adams_Bashforth_Moulton()
        else:
            return False, "Error mode. Need 'Implicit' or 'Explicit' for mode"

        if y_new is False:
            return False, x_new

        h = self.h_abs
        i = np.argmax(err/y_new)
        if err[i] < 1e-14:
            h_new = min([3/2 * h, self.max_step])

        else:
            h_hat = h*(self.atol[i]/err[i])**(1/(self.p+1))
            if h_hat > h:
                h_new = min([h_hat, 3/2 * h, self.max_step])
            else:
                h_new = max([h_hat, 1/2 * h, MIN_H])

        if h_new < MIN_H:
            return False, self.TOO_SMALL_STEP

        h_new = min(h_new, abs(self.t_bound - (self.t + h*self.direction)))
        self.t += self.h_abs*self.direction
        self.h_abs = h_new
        self.prev_f_y.append(self.f(self.t, y_new, x_new))
        self.nfev += 1
        if x_new is None:
            self.y = y_new
        else:
            self.y = np.concatenate((y_new, x_new), axis=None)
        if self.p < self.order and self.p < MAX_ORDER:
            self.p += 1

        return True, None

    def _dense_output_impl(self):
        prev_f = np.array([self.prev_f_y[-2], self.prev_f_y[-1]])
        return SABMDenseOutput(self.t_old, self.t, self.y, prev_f,
                               self.num_diff, self.g, self.ode)

    def solve_constraints(self, y, x0):
        def equations(z):
            return self.g(self.t + self.h_abs * self.direction, y, z)

        (x, dic, flag, message) = fsolve(equations,
                                         x0=x0, xtol=self.rtol, full_output=1)
        if flag != 1:
            return (False, "Solving the constraints fail to converge.")
        self.nfev += dic["nfev"]

        return (True, x)

    def Explicit_Adams_Bashforth(self):
        p = self.p
        coeff = B[p-1]
        h = self.h_abs

        y_p = np.copy(self.y[:self.num_diff])
        for i in range(0, p):
            y_p += h*self.direction*coeff[i]*self.prev_f_y[p-i-1]

        if self.ode:
            return (y_p, None)

        (flag, x_p) = self.solve_constraints(y_p, self.y[self.num_diff:])
        if not flag:
            return (flag, x_p)
        return (y_p, x_p)

    def Semi_Explicit_Adams_Bashforth_Moulton(self):
        y_p, x_p = self.Explicit_Adams_Bashforth()
        if y_p is False:
            return (y_p, x_p, False)
        h = self.h_abs
        p = self.p
        coeff = M[p-1]
        y_new = np.copy(self.y[:self.num_diff])
        direction = self.direction
        for i in range(self.num_diff):
            P = np.concatenate((y_new[0:i], y_p[i:]), axis=None)
            y_new[i] += h*direction*coeff[0] * \
                self.f(self.t+h*direction, P, x_p)[i]
            self.nfev += 1
            for j in range(1, p):
                y_new[i] += h*direction*coeff[j]*self.prev_f_y[p-j][i]

            error = abs(y_new[i] - y_p[i])
            if error > self.atol[i] and h > MIN_H:
                self.h_abs = h/2
                return self.Semi_Explicit_Adams_Bashforth_Moulton()
            elif error > self.atol[i] and h < MIN_H:
                return (False, """The minimal step size is reached. The method
                        doesn't converge.""", None)

        error = np.abs(y_new - y_p)
        if self.ode:
            return (y_new, None, error)

        (flag, x_new) = self.solve_constraints(y_new, x_p)
        if not flag:
            return (flag, x_new, None)

        return (y_new, x_new, error)

    def Semi_Implicit_Adams_Bashforth_Moulton(self):
        y_p, x_p = self.Explicit_Adams_Bashforth()
        if y_p is False:
            return (y_p, x_p, False)
        h = self.h_abs
        p = self.p
        coeff = M[p-1]
        y_new = np.copy(self.y[:self.num_diff])
        direction = self.direction
        for i in range(self.num_diff):
            def equations(z):
                P = np.concatenate((y_new[0:i], z, y_p[i+1:]), axis=None)
                res = y_new[i] + h*direction*coeff[0] * \
                    self.f(self.t + h*direction, P, x_p)[i]
                for j in range(1, p):
                    res += h*direction*coeff[j]*self.prev_f_y[p-j][i]
                return z - res
            (y_new[i], dic, flag, message) = fsolve(
                equations, x0=y_p[i], xtol=self.rtol, full_output=1)
            self.nfev += dic["nfev"]
            if flag != 1:
                error = self.atol[i] + 1
            else:
                error = abs(y_new[i] - y_p[i])
            if error > self.atol[i] and h > MIN_H:
                self.h_abs = h/2
                return self.Semi_Implicit_Adams_Bashforth_Moulton()
            elif error > self.atol[i] and h < MIN_H:
                return (False, """The minimal step size is reached. The method
                        doesn't converge.""", None)

        error = np.abs(y_new - y_p)
        if self.ode:
            return (y_new, None, error)

        (flag, x_new) = self.solve_constraints(y_new, x_p)
        if not flag:
            return (flag, x_new, None)

        return (y_new, x_new, error)


class SABMDenseOutput(DenseOutput):
    def __init__(self, t_old, t, y, fy, num_diff, g, ode):
        super().__init__(t_old, t)
        self.y = y[:num_diff]
        self.x = y[num_diff:]
        self.fy = fy
        self.t0 = t_old
        self.t1 = t
        self.g = g
        self.ode = ode

    def _call_impl(self, t):

        a = (self.fy[1] - self.fy[0])/(2*(self.t1 - self.t0))
        b = self.fy[0] - 2*a*self.t0
        c = self.y - a*self.t1**2 - b*self.t1

        num_diff = len(self.y)
        if t.ndim == 0:
            y_h = a*t**2 + b*t + c
            if self.ode:
                return y_h

            def equations(z):
                return self.g(t, y_h, z)
            x_h = fsolve(equations, x0=self.x)
            return np.concatenate((y_h, x_h))

        else:
            sol = np.zeros((num_diff+len(self.x), len(t)))
            sol[:num_diff, :] = np.outer(
                a, t**2) + np.outer(b, t) + np.outer(c, np.ones_like(t))

            if self.ode:
                return sol

            for i in range(len(t)):
                def equations(z):
                    return self.g(t[i], sol[:num_diff, i], z)
                sol[num_diff:, i] = fsolve(equations, x0=self.x)

        return sol
