## Automatically adapted for scipy Oct 21, 2005 by

# Author: Travis Oliphant

__all__ = ['odeint']

import _odepack
from copy import copy

_msgs = {2: "Integration successful.",
         -1: "Excess work done on this call (perhaps wrong Dfun type).",
         -2: "Excess accuracy requested (tolerances too small).",
         -3: "Illegal input detected (internal error).",
         -4: "Repeated error test failures (internal error).",
         -5: "Repeated convergence failures (perhaps bad Jacobian or tolerances).",
         -6: "Error weight became zero during problem.",
         -7: "Internal workspace insufficient to finish (internal error)."
         }

def odeint(func, y0, t, args=(), Dfun=None, col_deriv=0, full_output=0,
           ml=None, mu=None, rtol=None, atol=None, tcrit=None, h0=0.0,
           hmax=0.0, hmin=0.0, ixpr=0, mxstep=0, mxhnil=0, mxordn=12,
           mxords=5, printmessg=0):
    """Integrate a system of ordinary differential equations.

    Solve a system of ordinary differential equations using lsoda from the
    FORTRAN library odepack.

    Solves the initial value problem for stiff or non-stiff systems
    of first order ode-s::

        dy/dt = func(y,t0,...)

    where y can be a vector.

    Parameters
    ----------
    func : callable(y, t0, ...)
        Computes the derivative of y at t0.
    y0 : array
        Initial condition on y (can be a vector).
    t : array
        A sequence of time points for which to solve for y.  The initial
        value point should be the first element of this sequence.
    args : tuple
        Extra arguments to pass to function.
    Dfun : callable(y, t0, ...)
        Gradient (Jacobian) of func.
    col_deriv : boolean
        True if Dfun defines derivatives down columns (faster),
        otherwise Dfun should define derivatives across rows.
    full_output : boolean
        True if to return a dictionary of optional outputs as the second output
    printmessg : boolean
        Whether to print the convergence message

    Returns
    -------
    y : array, shape (len(y0), len(t))
        Array containing the value of y for each desired time in t,
        with the initial value y0 in the first row.

    infodict : dict, only returned if full_output == True
        Dictionary containing additional output information

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
    ml, mu : integer
        If either of these are not-None or non-negative, then the
        Jacobian is assumed to be banded.  These give the number of
        lower and upper non-zero diagonals in this banded matrix.
        For the banded case, Dfun should return a matrix whose
        columns contain the non-zero bands (starting with the
        lowest diagonal).  Thus, the return matrix from Dfun should
        have shape len(y0) * (ml + mu + 1) when ml >=0 or mu >=0
    rtol, atol : float
        The input parameters rtol and atol determine the error
        control performed by the solver.  The solver will control the
        vector, e, of estimated local errors in y, according to an
        inequality of the form::
            max-norm of (e / ewt) <= 1
        where ewt is a vector of positive error weights computed as::
            ewt = rtol * abs(y) + atol
        rtol and atol can be either vectors the same length as y or scalars.
    tcrit : array
        Vector of critical points (e.g. singularities) where integration
        care should be taken.
    h0 : float, (0: solver-determined)
        The step size to be attempted on the first step.
    hmax : float, (0: solver-determined)
        The maximum absolute step size allowed.
    hmin : float, (0: solver-determined)
        The minimum absolute step size allowed.
    ixpr : boolean
        Whether to generate extra printing at method switches.
    mxstep : integer, (0: solver-determined)
        Maximum number of (internally defined) steps allowed for each
        integration point in t.
    mxhnil : integer, (0: solver-determined)
        Maximum number of messages printed.
    mxordn : integer, (0: solver-determined)
        Maximum order to be allowed for the nonstiff (Adams) method.
    mxords : integer, (0: solver-determined)
        Maximum order to be allowed for the stiff (BDF) method.

    See Also
    --------
    ode : a more object-oriented integrator based on VODE
    quad : for finding the area under a curve

    """

    if ml is None:
        ml = -1 # changed to zero inside function call
    if mu is None:
        mu = -1 # changed to zero inside function call
    t = copy(t)
    y0 = copy(y0)
    output = _odepack.odeint(func, y0, t, args, Dfun, col_deriv, ml, mu,
                             full_output, rtol, atol, tcrit, h0, hmax, hmin,
                             ixpr, mxstep, mxhnil, mxordn, mxords)
    if output[-1] < 0:
        print _msgs[output[-1]]
        print "Run with full_output = 1 to get quantitative information."
    else:
        if printmessg:
            print _msgs[output[-1]]

    if full_output:
        output[1]['message'] = _msgs[output[-1]]

    output = output[:-1]
    if len(output) == 1:
        return output[0]
    else:
        return output
