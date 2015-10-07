"""
Functions
---------
.. autosummary::
   :toctree: generated/

    fmin_bfgs_h

"""
# ******NOTICE***************
# bfgs_h.py module by Henry C. Herbol
#
# You may copy and use this module as you see fit with no
# guarantee implied provided you keep this notice in all copies.
#
# Notice header adapted from optimize.py by Travis E. Oliphant
# Code adapted from scipy.optimize._minmize_bfgs by Travis E. Oliphant
# *****END NOTICE************
#
# An adapted bfgs optimization algorithm to optimize independent of
# target functions.

import numpy
from numpy import (Inf, sqrt, isinf)
from .optimize import (vecnorm, wrap_function, _check_unknown_options, OptimizeResult, approx_fprime)

__all__ = ['fmin_bfgs_h', '_minimize_bfgs_h']

_epsilon = sqrt(numpy.finfo(float).eps)
# standard status messages of optimizers
_status_message = {'success': 'Optimization terminated successfully.',
                   'maxfev': 'Maximum number of function evaluations has '
                              'been exceeded.',
                   'maxiter': 'Maximum number of iterations has been '
                              'exceeded.',
                   'pr_loss': 'Desired error not necessarily achieved due '
                              'to precision loss.'}

def fmin_bfgs_h(f, x0, fprime=None, args=(), gtol=1e-5, alpha=0.5, beta=0.7,
              H_reset=True, norm=Inf, epsilon=_epsilon, maxiter=None,
              full_output=0, disp=True, retall=False, callback=None):
    """
    Gradient Minimization independent of target function using the BFGS
    algorithm while applying the inverse hessian directly to steps.  This
    uses the backtracking line search method to adjust step size alpha if
    function (fun) is defined.  If no objective function is defined, passing
    None for function will run the minimization for a constant alpha for 
    maxiter iterations.

    Parameters
    ----------
    f : callable f(x,*args)
        Objective function to be minimized. If none available, pass None.
    x0 : ndarray
        Initial guess.
    fprime : callable f'(x,*args), optional
        Gradient of f.
    args : tuple, optional
        Extra arguments passed to f and fprime.
    gtol : float, optional
        Gradient norm must be less than gtol before successful termination.
    alpha : float, optional
        BFGS step size.
    beta : float, optional
        Adjustment to step size alpha whenever f(x0) indicates alpha is too
        large.
    H_reset : bool, optional
        Reset inverse hessian to the identity if True.
    norm : float, optional
        Order of norm (Inf is max, -Inf is min)
    epsilon : int or ndarray, optional
        If fprime is approximated, use this value for the step size.
    maxiter : int, optional
        Maximum number of iterations to perform.
    full_output : bool, optional
        If True,return fopt, func_calls, grad_calls, and warnflag
        in addition to xopt.
    disp : bool, optional
        Print convergence message if True.
    retall : bool, optional
        Return a list of results at each iteration if True.
    callback : callable, optional
        An optional user-supplied function to call after each
        iteration.  Called as callback(xk), where xk is the
        current parameter vector.

    Returns
    -------
    xopt : ndarray
        Parameters which minimize f, i.e. f(xopt) == fopt.
    fopt : float
        Minimum value. If function f was None, then fopt == None.
    gopt : ndarray
        Value of gradient at minimum, f'(xopt), which should be near 0.
    Bopt : ndarray
        Value of 1/f''(xopt), i.e. the inverse hessian matrix.
    func_calls : int
        Number of function_calls made. If function f was None, then 0.
    grad_calls : int
        Number of gradient calls made.
    warnflag : integer
        1 : Maximum number of iterations exceeded.
        2 : Gradient and/or function calls not changing.
    allvecs  :  list
        `OptimizeResult` at each iteration.  Only returned if retall is True.

    See also
    --------
    minimize: Interface to minimization algorithms for multivariate
        functions. See the 'BFGS-H' `method` in particular.

    Notes
    -----
    Optimize the function, f, whose gradient is given by fprime
    using the quasi-Newton method of Broyden, Fletcher, Goldfarb,
    and Shanno (BFGS).  This method is target function independent
    so it adjusts step size accordingly.

    Optimization of the function using the Broyden-Fletcher-Goldfarb-
    Shanno (BFGS) algorithm independent of target functions.  If an objective
    function is passed, the line search is calculated by the backtracking
    line search method (if None, then a constant alpha is used for maxiter
    loops):

    * alpha is a constant that is changed by beta whenever f(x0)
      indicates the maximum has increased.
    * rho, a scalar used in scipy.optimize._minimize_bfgs, was removed as
      alpha is now manually adjusted.

    Method is an application of the BFGS(Hess) function described by Sheppard
    et al.

    References
    ----------
    * Wright, and Nocedal 'Numerical Optimization', 1999, pg. 198.
    * Fletcher, R., and Powell, M. J. D. (1963) ( A rapidly convergent descent
      method for minimization.) The Computer Journal, Volume 6, pp. 163-168
    * Broyden, C. G. (1966) Quasi-Newton Methods and their Application to
      Functional Minimisation Math. Comp., Volume 21, pp. 368-381
    * Broyden, C. G. (1970) The Convergence of Single-Rank Quasi-Newton
      Methods Math. Comp., Volume 24, pp. 365-382
    * Sheppard, D., Terrell, R., and Henkelman, G. (2007) Optimization Methods
      for Finding Minimum Energy Paths. doi: 10.1063/1.2841941

    .. versionadded:: 0.17.0
    """
    opts = {'gtol': gtol,
            'alpha': alpha,
            'beta': beta,
            'H_reset': H_reset,
            'norm': norm,
            'eps': epsilon,
            'disp': disp,
            'maxiter': maxiter,
            'return_all': retall}

    res = _minimize_bfgs_h(f, x0, args, fprime, callback=callback, **opts)

    if full_output:
        retlist = (res['x'], res['fun'], res['jac'], res['hess_inv'],
                   res['nfev'], res['njev'], res['status'])
        if retall:
            retlist += (res['allvecs'], )
        return retlist
    else:
        if retall:
            return res['x'], res['allvecs']
        else:
            return res['x']

def _minimize_bfgs_h(fun, x0, args=(), jac=None, callback=None, 
    gtol=1e-5, alpha=0.5, beta=0.7, H_reset=True, norm=Inf,
    eps=_epsilon, maxiter=None, disp=False, return_all=False,
    **unknown_options):
    """
    Gradient Minimization independent of target function using the BFGS
    algorithm while applying the inverse hessian directly to steps.  This
    uses the backtracking line search method to adjust step size alpha if
    function (fun) is defined.  If no objective function is defined, passing
    None for function will run the minimization for a constant alpha for 
    maxiter iterations.

    Options
    -------
    gtol : float, optional
        Gradient norm must be less than gtol for successful convergence.
    alpha : float, optional
        BFGS step size.
    beta : float, optional
        Adjustment to step size alpha whenever f(x0) indicates alpha is too
        large.
    H_reset : bool, optional
        Reset inverse hessian to the identity if True.
    norm : float, optional
        Order of norm (Inf is max, -Inf is min)
    eps : float or ndarray
        If `jac` is approximated, use this value for the step size.
    disp : bool
        Set to True to print convergence messages.
    maxiter : int
        Maximum number of iterations to perform.
    return_all : bool, optional
        Returns all display output in list if True
   
    Returns
    -------
    xopt : ndarray
        Parameters which minimize f, i.e. f(xopt) == fopt.
    fopt : float
        Minimum value. If function f was None, then fopt == None.
    gopt : ndarray
        Value of gradient at minimum, f'(xopt), which should be near 0.
    Bopt : ndarray
        Value of 1/f''(xopt), i.e. the inverse hessian matrix.
    func_calls : int
        Number of function_calls made. If function f was None, then 0.
    grad_calls : int
        Number of gradient calls made.
    warnflag : integer
        1 : Maximum number of iterations exceeded.
        2 : Gradient and/or function calls not changing.
    allvecs  :  list
        `OptimizeResult` at each iteration.  Only returned if retall is True.

    Notes
    -----
    Optimization of the function using the Broyden-Fletcher-Goldfarb-
    Shanno (BFGS) algorithm independent of target functions.  Instead
    of using a line search for the step size alpha, the following 
    adjustments have been made:

    * alpha is a constant that is changed by beta whenever f(x0)
      indicates the maximum has increased.
    * rho, a scalar used in scipy.optimize._minimize_bfgs, was removed as
      alpha is now manually adjusted.

    Method is an application of the BFGS(Hess) function described by Sheppard
    et al.

    References
    ----------
    * Fletcher, R., and Powell, M. J. D. (1963) ( A rapidly convergent descent
      method for minimization.) The Computer Journal, Volume 6, pp. 163-168
    * Broyden, C. G. (1966) Quasi-Newton Methods and their Application to
      Functional Minimisation Math. Comp., Volume 21, pp. 368-381
    * Broyden, C. G. (1970) The Convergence of Single-Rank Quasi-Newton
      Methods Math. Comp., Volume 24, pp. 365-382
    * Sheppard, D., Terrell, R., and Henkelman, G. (2007) Optimization Methods
      for Finding Minimum Energy Paths. doi: 10.1063/1.2841941

    .. versionadded:: 0.17.0
    """

    # Interface code with scipy output
    ##################################
    _check_unknown_options(unknown_options)
    f = fun
    fprime = jac
    epsilon = eps
    retall = return_all

    if maxiter is None:
        maxiter = len(x0) * 200

    func_calls, f = wrap_function(f, args)
    if fprime is None:
        grad_calls, myfprime = wrap_function(approx_fprime, (f, epsilon))
    else:
        grad_calls, myfprime = wrap_function(fprime, args)
    ##################################

    if beta > 1:
        if disp:
            print("Warning - Unreasonable Beta (must be less than or\
                        equal to 1). Setting to 1.")
        beta = 1.0

    # Get x0 as a flat array
    x0 = numpy.asarray(x0).flatten()
    if x0.ndim == 0:
        x0.shape = (1,)

    gfk = fprime(x0)

    k,N = 0,len(x0)
    I = numpy.eye(N, dtype=int)

    # Initialize inv Hess as Identity matrix
    Hk = I

    # Store your old func_max
    if f is not None:
        old_fval = f(x0, *args)

    xk = x0

    # Criteria on gradient for continuing simulation
    gnorm = vecnorm(gfk, ord=norm)
    
    warnflag = 0
    if retall:
        allvecs = [x0]
    # Main Loop:
    if disp:
        print("Alpha, Beta, H_reset = %lg, %lg, %s"
                   % (alpha,beta,str(H_reset)))
    while (gnorm > gtol) and (k < maxiter):
        if disp:
            print("Step %d, " % k),
        # Get your step direction
        pk = -numpy.dot(Hk, gfk)

        # If we are doing unreasonably small step sizes, quit
        if numpy.linalg.norm(pk*alpha) < 1E-7:
            if disp:
                print("Error - Step size unreasonable (%lg)" 
                            % numpy.linalg.norm(pk*alpha))
            warnflag = 2
            break

        # Hold new parameters
        xkp1 = xk + alpha * pk

        # Recalculate sk to maintain the secant condition
        sk = xkp1 - xk

        # Get the new gradient
        gfkp1 = fprime(xkp1)

        # Check if max has increased
        if f is not None:
            fval = f(xkp1, *args)
        if f is not None and fval > old_fval:
            # Step taken overstepped the minimum.  Lowering step size
            if disp:
                print("\tResetting System as %lg > %lg!"
                        % (fval, old_fval))
                print("\talpha: %lg" % alpha),
            alpha *= float(beta)
            if disp:
                print("-> %lg\n" % alpha)

            # Reset the Inverse Hessian if desired - This is recommended!
            if H_reset:
                Hk = I
            continue
        
        # Store new parameters, as it has passed the check
        # (fval < old_fval is True)
        xk = xkp1

        # Store for output if desired
        if retall:
            allvecs.append(xk)
        
        # Store new max value in old_max for future comparison
        if f is not None:
            old_fval = fval

        # Get difference in gradients for further calculations
        yk = gfkp1 - gfk
        # Store new gradient in old gradient
        gfk = gfkp1

        # Update the conditional check
        gnorm = vecnorm(gfk, ord=norm)

        if disp:
            print("gnorm %lg" % gnorm)

        # If callback is desired
        if callback is not None:
            callback(xk)

        # Increment the loop counter
        k += 1
        gnorm = vecnorm(gfk, ord=norm)
        if (gnorm <= gtol):
            break

        try:  # this was handled in numeric, let it remaines for more safety
            rhok = 1.0 / (numpy.dot(yk, sk))
        except ZeroDivisionError:
            rhok = 1000.0
            if disp:
                print("Divide-by-zero encountered: rhok assumed large")
        if isinf(rhok):  # this is patch for numpy
            rhok = 1000.0
            if disp:
                print("Divide-by-zero encountered: rhok assumed large")
        # Run BFGS Update for the Inverse Hessian
        A1 = I - sk[:, numpy.newaxis] * yk[numpy.newaxis, :] * rhok
        A2 = I - yk[:, numpy.newaxis] * sk[numpy.newaxis, :] * rhok
        Hk = numpy.dot(A1, numpy.dot(Hk, A2)) + \
             (rhok * sk[:, numpy.newaxis] * sk[numpy.newaxis, :])

    if f is not None:
        fval = old_fval
    else:
        fval = float('NaN')

    if numpy.isnan(fval):
        # This can happen if the first call to f returned NaN;
        # the loop is then never entered.
        warnflag = 2

    if warnflag == 2:
        msg = _status_message['pr_loss']
        if disp:
            print("Warning: " + msg)
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % k)
            print("         Function evaluations: %d" % func_calls[0])
            print("         Gradient evaluations: %d" % grad_calls[0])

    elif k >= maxiter:
        warnflag = 1
        msg = _status_message['maxiter']
        if disp:
            print("Warning: " + msg)
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % k)
            print("         Function evaluations: %d" % func_calls[0])
            print("         Gradient evaluations: %d" % grad_calls[0])
    else:
        msg = _status_message['success']
        if disp:
            print(msg)
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % k)
            print("         Function evaluations: %d" % func_calls[0])
            print("         Gradient evaluations: %d" % grad_calls[0])

    result = OptimizeResult(fun=fval, jac=gfk, hess_inv=Hk, nfev=func_calls[0],
                            njev=grad_calls[0], status=warnflag,
                            success=(warnflag == 0), message=msg, x=xk,
                            nit=k)
    if retall:
        result['allvecs'] = allvecs
    return result
