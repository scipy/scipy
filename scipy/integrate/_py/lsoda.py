import numpy as np
from scipy.integrate import ode
from .common import validate_tol, warn_extraneous
from .base import OdeSolver, DenseOutput


class LSODA(OdeSolver):
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
        self.jac = jac
        self.lband = lband
        self.uband = uband

        solver = ode(self.fun, jac)
        solver.set_integrator('lsoda', rtol=rtol, atol=atol, max_step=max_step,
                              min_step=min_step, first_step=first_step,
                              lband=lband, uband=uband)
        solver.set_initial_value(y0, t0)

        # Inject t_crit into rwork array as needed for itask=5
        solver._integrator.rwork[0] = self.t_crit
        solver._integrator.call_args[4] = solver._integrator.rwork

        self._lsoda_solver = solver
        self.y_old = None

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
            self.y_old = self.y
            self.y = solver._y
            self.njev = integrator.iwork[12]
            return True, None
        else:
            return False, 'Unexpected istate in LSODA.'

    def _dense_output_impl(self):
        return LsodaDenseOutput(self.t_old, self.t, self.y_old, self.y)


class LsodaDenseOutput(DenseOutput):
    def __init__(self, t_old, t, y_old, y):
        super(LsodaDenseOutput, self).__init__(t_old, t)
        self.y_old = y_old
        self.y = y

    def _call_impl(self, t):
        if np.ndim(t) > 0:
            y_old = np.expand_dims(self.y_old, 1)
            y = np.expand_dims(self.y, 1)
        else:
            y_old = self.y_old
            y = self.y

        return y_old + (t - self.t_old) * (y - y_old) / (self.t - self.t_old)
