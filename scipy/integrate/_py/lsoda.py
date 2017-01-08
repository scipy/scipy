import numpy as np
from scipy.integrate import ode
from .common import validate_tol, warn_extraneous
from .base import OdeSolver, DenseOutput


class LSODA(OdeSolver):
    """Lawrence Livermore solver for ODEs with automatic detection of nonstiff
    and stiff problems

    This solver switches automatically between the nonstiff Adams method and
    the stiff BDF method. The method was originally detailed in [1]_. This
    wrapper calls the compiled Fortran implementation included in ODEPACK [2]_.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system. The calling signature is ``fun(t, y)``.
        Here ``t`` is a scalar and there are two options for ndarray ``y``.
        It can either have shape (n,), then ``fun`` must return array_like with
        shape (n,). Or alternatively it can have shape (n, k), then ``fun``
        must return array_like with shape (n, k), i.e. each column
        corresponds to a single column in ``y``. The choice between the two
        options is determined by `vectorized` argument (see below). The
        vectorized implementation allows faster approximation of the Jacobian
        by finite differences.
    t0 : float
        Initial time.
    y0 : array_like, shape (n,)
        Initial state.
    t_crit : float
        Boundary time --- the integration won't continue beyond it. It also
        determines the direction of the integration.
    first_step : float or None, optional
        Initial step size. Default is ``None`` which means that the algorithm
        should choose.
    min_step : float, optional
        Minimum allowed step size. Default is 0.0, i.e. the step is not
        bounded and determined solely by the solver.
    max_step : float, optional
        Maximum allowed step size. Default is ``np.inf``, i.e. the step is not
        bounded and determined solely by the solver.
    rtol, atol : float and array_like, optional
        Relative and absolute tolerances. The solver keeps the local error
        estimates less than ``atol + rtol * abs(y)``. Here `rtol` controls a
        relative accuracy (number of correct digits). But if a component of `y`
        is approximately below `atol` then the error only needs to fall within
        the same `atol` threshold, and the number of correct digits is not
        guaranteed. If components of y have different scales, it might be
        beneficial to set different `atol` values for different components by
        passing array_like with shape (n,) for `atol`. Default values are
        1e-3 for `rtol` and 1e-6 for `atol`.
    jac : None or callable, optional
        Jacobian matrix of the right-hand side of the system with respect to
        ``y``. The Jacobian matrix has shape (n, n) and its element (i, j) is
        equal to ``d f_i / d y_j``. The function will be called as
        ``jac(t, y)``. If ``None`` (default), then the Jacobian will be
        approximated by finite differences. It is generally recommended to
        provide the Jacobian rather than relying on a finite difference
        approximation.
    lband, uband : int or None, optional
        Jacobian band width, ``jac[i,j] != 0 for i-lband <= j <= i+uband``.
        Setting these requires your jac routine to return the jacobian in packed
        format, ``jac_packed[i-j+uband, j] = jac[i,j]``.
    vectorized : bool, optional
        Whether `fun` is implemented in a vectorized fashion. Default is False.

    Attributes
    ----------
    n : int
        Number of equations.
    status : string
        Current status of the solver: 'running', 'finished' or 'failed'.
    t_crit : float
        Boundary time.
    direction : float
        Integration direction: +1 or -1.
    t : float
        Current time.
    y : ndarray
        Current state.
    t_old : float
        Previous time. None if no steps were made yet.
    nfev : int
        Number of the system's rhs evaluations.
    njev : int
        Number of the Jacobian evaluations.

    References
    ----------
    .. [1] L. Petzold, "Automatic selection of methods for solving stiff and
    nonstiff systems of ordinary differential equations", SIAM Journal on
    Scientific and Statistical Computing, Vol. 4, No. 1, pp. 136-148, 1983.
    .. [2] A. C. Hindmarsh, "ODEPACK, A Systematized Collection of ODE Solvers,"
    IMACS Transactions on Scientific Computation, Vol 1., pp. 55-64, 1983.
    """

    def __init__(self, fun, t0, y0, t_crit, first_step=None, min_step=0.0,
                 max_step=np.inf, rtol=1e-3, atol=1e-6, jac=None, lband=None,
                 uband=None, vectorized=False, **extraneous):
        warn_extraneous(extraneous)
        super(LSODA, self).__init__(fun, t0, y0, t_crit, vectorized)

        if first_step is None:
            # 0.0 is LSODA value for default
            first_step = 0.0
        elif first_step <= 0.0:
            raise ValueError("`first_step` must be positive or None.")

        if max_step == np.inf:
            # 0.0 is LSODA value for default
            max_step = 0.0
        elif max_step <= 0.0:
            raise ValueError("`max_step` must be positive.")

        if min_step < 0.0:
            raise ValueError("`min_step` must be nonnegative.")

        rtol, atol = validate_tol(rtol, atol, self.n)

        self.first_step = first_step
        self.max_step = max_step
        self.min_step = min_step
        self.rtol = rtol
        self.atol = atol
        self.jac = self._validate_jac(jac)
        self.lband = lband
        self.uband = uband

        solver = ode(self.fun, self.jac)
        solver.set_integrator('lsoda', rtol=rtol, atol=atol, max_step=max_step,
                              min_step=min_step, first_step=first_step,
                              lband=lband, uband=uband)
        solver.set_initial_value(y0, t0)

        # Inject t_crit into rwork array as needed for itask=5
        solver._integrator.rwork[0] = self.t_crit
        solver._integrator.call_args[4] = solver._integrator.rwork

        self._lsoda_solver = solver

    def _validate_jac(self, jac):
        if jac is None:
            jac_wrapped = None
        elif callable(jac):
            def jac_wrapped(t, y):
                self.njev += 1
                return np.asarray(jac(t, y), dtype=float)
        else:
            raise ValueError("`jac` is expected to be `None` or `callable`, "
                             "but is actually type {}".format(type(jac)))

        return jac_wrapped

    def _step_impl(self):
        solver = self._lsoda_solver
        integrator = solver._integrator

        # From lsoda.step and lsoda.integrate
        # itask=5 means take a single step and do not go past t_crit
        itask = integrator.call_args[2]
        integrator.call_args[2] = 5
        solver._y, solver.t = integrator.run(solver.f,
                                             solver.jac or (lambda: None),
                                             solver._y, solver.t, self.t_crit,
                                             solver.f_params, solver.jac_params)
        integrator.call_args[2] = itask

        if solver.successful():
            self.t = solver.t
            self.y = solver._y
            self.njev = integrator.iwork[12]
            return True, None
        else:
            return False, 'Unexpected istate in LSODA.'

    def _dense_output_impl(self):
        iwork = self._lsoda_solver._integrator.iwork
        rwork = self._lsoda_solver._integrator.rwork

        order = iwork[14]
        h = rwork[11]
        yh = np.reshape(rwork[20:20 + (order + 1) * self.n],
                        (self.n, order + 1), order='F').copy()

        return LsodaDenseOutput(self.t_old, self.t, h, order, yh)


class LsodaDenseOutput(DenseOutput):
    def __init__(self, t_old, t, h, order, yh):
        super(LsodaDenseOutput, self).__init__(t_old, t)
        self.h = h
        self.yh = yh
        self.p = np.arange(order + 1)

    def _call_impl(self, t):
        if t.ndim == 0:
            x = ((t - self.t) / self.h) ** self.p
        else:
            x = ((t - self.t) / self.h) ** self.p[:, None]

        return np.dot(self.yh, x)
