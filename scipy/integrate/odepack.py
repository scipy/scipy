# Author: Travis Oliphant
#
# Aug 2013: Juan Luis Cano
#   Rewritten odeint loop in Fortran.

from __future__ import division, print_function, absolute_import

__all__ = ['odeint']

import warnings
from collections import namedtuple
from . import _pyodepack

_msgs = {
    2: "Integration successful.",
    -1: "Excess work done on this call (perhaps wrong dfunc type).",
    -2: "Excess accuracy requested (tolerances too small).",
    -3: "Illegal input detected (internal error).",
    -4: "Repeated error test failures (internal error).",
    -5: "Repeated convergence failures (perhaps bad Jacobian or tolerances).",
    -6: "Error weight became zero during problem.",
    -7: "Internal workspace insufficient to finish (internal error)."
}

Result = namedtuple('Result', 'y success infodict')


def odeint(func, y0, t, dfunc=None, col_deriv=0,
           ml=None, mu=None, rtol=None, atol=None, tcrit=None, first_step=0.0,
           max_step=0.0, min_step=0.0, ixpr=0, max_nosteps=0, max_msgs=0,
           max_order_ns=12, max_order_s=5):
    """Integrate a system of ordinary differential equations.

    Solve a system of ordinary differential equations using lsoda from the
    FORTRAN library odepack.

    Solves the initial value problem for stiff or non-stiff systems
    of first order ode-s::

        dy/dt = func(t0, y, ...)

    where y can be a vector.

    Parameters
    ----------
    func : callable(t0, y, ...)
        Computes the derivative of y at t0.
    y0 : array
        Initial condition on y (can be a vector).
    t : array
        A sequence of time points for which to solve for y.  The initial
        value point should be the first element of this sequence.
    dfunc : callable(t0, y, ...)
        Gradient (Jacobian) of `func`.
    col_deriv : bool, optional
        True if `Dfun` defines derivatives down columns (faster),
        otherwise `Dfun` should define derivatives across rows.

    Returns
    -------
    Result
        The integration result represented as a `Result` object.
        Its attributes are: `y`, the solution array, `success`, a boolean
        flag indicating if the solver exited successfully and `infodict`
        with additional output information:

        =======  ============================================================
        key      meaning
        =======  ============================================================
        'hu'     vector of step sizes successfully used for each time step.
        'tcur'   vector with the value of t reached for each time step.
                 (will always be at least as large as the input times).
        'tolsf'  vector of tolerance scale factors, greater than 1.0,
                 computed when a request for too much accuracy was detected.
        'tsw'    value of t at the time of the last method switch
                 (given for each time step)
        'nst'    cumulative number of time steps
        'nfe'    cumulative number of function evaluations for each time step
        'nje'    cumulative number of jacobian evaluations for each time step
        'nqu'    a vector of method orders for each successful step.
        'imxer'  index of the component of largest magnitude in the
                 weighted local error vector (e / ewt) on an error return, -1
                 otherwise.
        'lenrw'  the length of the double work array required.
        'leniw'  the length of integer work array required.
        'mused'  a vector of method indicators for each successful time step:
                 1: adams (nonstiff), 2: bdf (stiff)
        =======  ============================================================

    Other Parameters
    ----------------
    ml, mu : int, optional
        If either of these are not None or non-negative, then the
        Jacobian is assumed to be banded.  These give the number of
        lower and upper non-zero diagonals in this banded matrix.
        For the banded case, `Dfun` should return a matrix whose
        columns contain the non-zero bands (starting with the
        lowest diagonal).  Thus, the return matrix from `Dfun` should
        have shape ``len(y0) * (ml + mu + 1)`` when ``ml >=0`` or ``mu >=0``.
    rtol, atol : float, optional
        The input parameters `rtol` and `atol` determine the error
        control performed by the solver.  The solver will control the
        vector, e, of estimated local errors in y, according to an
        inequality of the form ``max-norm of (e / ewt) <= 1``,
        where ewt is a vector of positive error weights computed as
        ``ewt = rtol * abs(y) + atol``.
        rtol and atol can be either vectors the same length as y or scalars.
        Defaults to 1.49012e-8.
    tcrit : ndarray, optional
        Vector of critical points (e.g. singularities) where integration
        care should be taken.
    first_step : float, (0: solver-determined), optional
        The step size to be attempted on the first step.
    max_step : float, (0: solver-determined), optional
        The maximum absolute step size allowed.
    min_step : float, (0: solver-determined), optional
        The minimum absolute step size allowed.
    ixpr : bool, optional
        Whether to generate extra printing at method switches.
    max_nosteps : int, (0: solver-determined), optional
        Maximum number of (internally defined) steps allowed for each
        integration point in t.
    max_msgs : int, (0: solver-determined), optional
        Maximum number of messages printed.
    max_order_ns : int, (0: solver-determined), optional
        Maximum order to be allowed for the non-stiff (Adams) method.
    max_order_s : int, (0: solver-determined), optional
        Maximum order to be allowed for the stiff (BDF) method.

    See Also
    --------
    ode : a more object-oriented integrator based on VODE.
    quad : for finding the area under a curve.

    """

    # TODO: Jacobian
    # TODO: Critical points
    tol = 1.49012e-8
    if rtol is None:
        rtol = tol
    if atol is None:
        atol = tol
    if ml is None:
        ml = -1
    if mu is None:
        mu = -1

    jt = 1
    if dfunc is None:
        jt = 2
        dfunc = lambda t, y: None

    odeint = _pyodepack.pyodepack.odeint
    y, iostate, rout, iout = odeint(func, y0, t, rtol, atol, tcrit,
                                    first_step, max_step, min_step,
                                    dfunc, jt, ml, mu,
                                    ixpr, max_nosteps, max_msgs,
                                    max_order_ns, max_order_s)
    if iostate == 0:
        raise MemoryError('Could not allocate work arrays')
    elif iostate < 0:
        with warnings.catch_warnings():
            warnings.simplefilter("always")
            warnings.warn(RuntimeWarning(_msgs[iostate]))

    infodict = {
        'hu': rout[:, 0],
        'tcur': rout[:, 2],
        'tolsf': rout[:, 3],
        'tsw': rout[:, 4],
        'nst': iout[:, 0],
        'nfe': iout[:, 1],
        'nje': iout[:, 2],
        'nqu': iout[:, 3],
        'imxer': iout[:, 5],
        'lenrw': iout[:, 6],
        'leniw': iout[:, 7],
        'mused': iout[:, 8],
    }
    return Result(y, iostate == 2, infodict)
