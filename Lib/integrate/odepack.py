# Author: Travis Oliphant

__all__ = ['odeint']

import _odepack

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

    Description:

      Solve a system of ordinary differential equations Using lsoda from the
      FORTRAN library odepack.
      
      Solves the initial value problem for stiff or non-stiff systems
      of first order ode-s:
           dy/dt = func(y,t0,...) where y can be a vector.

    Inputs:

      func -- func(y,t0,...) computes the derivative of y at t0.
      y0   -- initial condition on y (can be a vector).
      t    -- a sequence of time points for which to solve for y.  The intial 
              value point should be the first element of this sequence.
      args -- extra arguments to pass to function.
      Dfun -- the gradient (Jacobian) of func (same input signature as func).
      col_deriv -- non-zero implies that Dfun defines derivatives down
                   columns (faster), otherwise Dfun should define derivatives
                   across rows.
      full_output -- non-zero to return a dictionary of optional outputs as
                     the second output.
      printmessg -- print the convergence message.

    Outputs: (y, {infodict,})

      y -- a rank-2 array containing the value of y in each row for each
           desired time in t (with the initial value y0 in the first row).

      infodict -- a dictionary of optional outputs:
        'hu'    : a vector of step sizes successfully used for each time step.
        'tcur'  : a vector with the value of t reached for each time step.
                  (will always be at least as large as the input times).
        'tolsf' : a vector of tolerance scale factors, greater than 1.0,
                  computed when a request for too much accuracy was detected.
        'tsw'   : the value of t at the time of the last method switch
                  (given for each time step).
        'nst'   : the cumulative number of time steps.
        'nfe'   : the cumulative number of function evaluations for eadh
                  time step.
        'nje'   : the cumulative number of jacobian evaluations for each
                  time step.
        'nqu'   : a vector of method orders for each successful step.
        'imxer' : index of the component of largest magnitude in the
                   weighted local error vector (e / ewt) on an error return.
        'lenrw' : the length of the double work array required.
        'leniw' : the length of integer work array required.
        'mused' : a vector of method indicators for each successful time step:
                  1 -- adams (nonstiff)
                  2 -- bdf (stiff)

    Additional Inputs:

      ml, mu -- If either of these are not-None or non-negative, then the
                Jacobian is assumed to be banded.  These give the number of
                lower and upper non-zero diagonals in this banded matrix.
                For the banded case, Dfun should return a matrix whose
                columns contain the non-zero bands (starting with the
                lowest diagonal).  Thus, the return matrix from Dfun should
                have shape len(y0) x (ml + mu + 1) when ml >=0 or mu >=0
      rtol -- The input parameters rtol and atol determine the error
      atol    control performed by the solver.  The solver will control the
              vector, e, of estimated local errors in y, according to an
              inequality of the form
                   max-norm of (e / ewt) <= 1
              where ewt is a vector of positive error weights computed as
                   ewt = rtol * abs(y) + atol
              rtol and atol can be either vectors the same length as y or
              scalars.
      tcrit -- a vector of critical points (e.g. singularities) where
               integration care should be taken.

       (For the next inputs a zero default means the solver determines it).

      h0 -- the step size to be attempted on the first step.
      hmax -- the maximum absolute step size allowed.
      hmin -- the minimum absolute step size allowed.
      ixpr -- non-zero to generate extra printing at method switches.
      mxstep -- maximum number of (internally defined) steps allowed
                for each integration point in t.
      mxhnil -- maximum number of messages printed.
      mxordn -- maximum order to be allowed for the nonstiff (Adams) method.
      mxords -- maximum order to be allowed for the stiff (BDF) method.

    """

    if ml is None:
        ml = -1 # changed to zero inside function call
    if mu is None:
        mu = -1 # changed to zero inside function call
    t = t.copy()
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



