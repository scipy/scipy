# ******NOTICE***************
# optimize.py module by Travis E. Oliphant
#
# You may copy and use this module as you see fit with no
# guarantee implied provided you keep this notice in all copies.
# *****END NOTICE************

# A collection of optimization algorithms.  Version 0.4.1
# CHANGES
#  Added fminbound (July 2001)

# Minimization routines
"""optimize.py

A collection of general-purpose optimization routines using Numeric

fmin        ---      Nelder-Mead Simplex algorithm (uses only function calls).
fminBFGS    ---      Quasi-Newton method (uses function and gradient).
fminNCG     ---      Line-search Newton Conjugate Gradient (uses function, 
                     gradient and hessian (if it's provided)).
fminbound   ---      Bounded minimization for scalar functions.

"""

import Numeric
import MLab
try:
    from fastumath import absolute, sqrt
except ImporError:
    from scipy.fastumath import absolute, sqrt
Num = Numeric
max = MLab.max
min = MLab.min
abs = absolute
__version__="0.4.1"

def rosen(x):  # The Rosenbrock function
    return MLab.sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)

def rosen_der(x):
    xm = x[1:-1]
    xm_m1 = x[:-2]
    xm_p1 = x[2:]
    der = MLab.zeros(x.shape,x.typecode())
    der[1:-1] = 200*(xm-xm_m1**2) - 400*(xm_p1 - xm**2)*xm - 2*(1-xm)
    der[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])
    der[-1] = 200*(x[-1]-x[-2]**2)
    return der

def rosen_hess(x):
        x = Num.asarray(x)
        H = MLab.diag(-400*x[:-1],1) - MLab.diag(400*x[:-1],-1)
        diagonal = Num.zeros(len(x),x.typecode())
        diagonal[0] = 1200*x[0]-400*x[1]+2
        diagonal[-1] = 200
        diagonal[1:-1] = 202 + 1200*x[1:-1]**2 - 400*x[2:]
        H = H + MLab.diag(diagonal)
        return H

def rosen_hess_p(x,p):
        x = Num.asarray(x)
        Hp = Num.zeros(len(x),x.typecode())
        Hp[0] = (1200*x[0]**2 - 400*x[1] + 2)*p[0] - 400*x[0]*p[1]
        Hp[1:-1] = -400*x[:-2]*p[:-2]+(202+1200*x[1:-1]**2-400*x[2:])*p[1:-1] \
                   -400*x[1:-1]*p[2:]
        Hp[-1] = -400*x[-2]*p[-2] + 200*p[-1]
        return Hp

        
def fmin(func, x0, args=(), xtol=1e-4, ftol=1e-4, maxiter=None, maxfun=None, 
         full_output=0, printmessg=1):
    """Minimize a function using the simplex algorithm.

    Description:
    
      Uses a Nelder-Mead simplex algorithm to find the minimum of function
      of one or more variables.

    Inputs:

      func -- the Python function or method to be minimized.
      x0 -- the initial guess.
      args -- extra arguments for func.

    Outputs: (xopt, {fopt, warnflag})

      xopt -- minimizer of function

      fopt -- value of function at minimum: fopt = func(xopt)
      warnflag -- Integer warning flag:
                  1 : 'Maximum number of function evaluations.'
                  2 : 'Maximum number of iterations.'

    Additional Inputs:

      xtol -- acceptable relative error in xopt for convergence.
      ftol -- acceptable relative error in func(xopt) for convergence.
      maxiter -- the maximum number of iterations to perform.
      maxfun -- the maximum number of function evaluations.
      full_output -- non-zero if fval and warnflag outputs are desired.
      printmessg -- non-zero to print convergence messages.
      
      """
    x0 = Num.asarray(x0)
    assert (len(x0.shape)==1)
    N = len(x0)
    if maxiter is None:
        maxiter = N * 200
    if maxfun is None:
        maxfun = N * 200

    rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
    one2np1 = range(1,N+1)

    sim = Num.zeros((N+1,N),x0.typecode())
    fsim = Num.zeros((N+1,),'d')
    sim[0] = x0
    fsim[0] = apply(func,(x0,)+args)
    nonzdelt = 0.05
    zdelt = 0.00025
    for k in range(0,N):
        y = Num.array(x0,copy=1)
        if y[k] != 0:
            y[k] = (1+nonzdelt)*y[k]
        else:
            y[k] = zdelt

        sim[k+1] = y
        f = apply(func,(y,)+args)
        fsim[k+1] = f

    ind = Num.argsort(fsim)
    fsim = Num.take(fsim,ind)  # sort so sim[0,:] has the lowest function value
    sim = Num.take(sim,ind,0)
    
    iterations = 1
    funcalls = N+1
    
    while (funcalls < maxfun and iterations < maxiter):
        if (max(Num.ravel(abs(sim[1:]-sim[0]))) <= xtol \
            and max(abs(fsim[0]-fsim[1:])) <= ftol):
            break

        xbar = Num.add.reduce(sim[:-1],0) / N
        xr = (1+rho)*xbar - rho*sim[-1]
        fxr = apply(func,(xr,)+args)
        funcalls = funcalls + 1
        doshrink = 0

        if fxr < fsim[0]:
            xe = (1+rho*chi)*xbar - rho*chi*sim[-1]
            fxe = apply(func,(xe,)+args)
            funcalls = funcalls + 1

            if fxe < fxr:
                sim[-1] = xe
                fsim[-1] = fxe
            else:
                sim[-1] = xr
                fsim[-1] = fxr
        else: # fsim[0] <= fxr
            if fxr < fsim[-2]:
                sim[-1] = xr
                fsim[-1] = fxr
            else: # fxr >= fsim[-2]
                # Perform contraction
                if fxr < fsim[-1]:
                    xc = (1+psi*rho)*xbar - psi*rho*sim[-1]
                    fxc = apply(func,(xc,)+args)
                    funcalls = funcalls + 1

                    if fxc <= fxr:
                        sim[-1] = xc
                        fsim[-1] = fxc
                    else:
                        doshrink=1
                else:
                    # Perform an inside contraction
                    xcc = (1-psi)*xbar + psi*sim[-1]
                    fxcc = apply(func,(xcc,)+args)
                    funcalls = funcalls + 1

                    if fxcc < fsim[-1]:
                        sim[-1] = xcc
                        fsim[-1] = fxcc
                    else:
                        doshrink = 1

                if doshrink:
                    for j in one2np1:
                        sim[j] = sim[0] + sigma*(sim[j] - sim[0])
                        fsim[j] = apply(func,(sim[j],)+args)
                    funcalls = funcalls + N

        ind = Num.argsort(fsim)
        sim = Num.take(sim,ind,0)
        fsim = Num.take(fsim,ind)
        iterations = iterations + 1

    x = sim[0]
    fval = min(fsim)
    warnflag = 0

    if funcalls >= maxfun:
        warnflag = 1
        if printmessg:
            print "Warning: Maximum number of function evaluations has "\
                  "been exceeded."
    elif iterations >= maxiter:
        warnflag = 2
        if printmessg:
            print "Warning: Maximum number of iterations has been exceeded"
    else:
        if printmessg:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % iterations
            print "         Function evaluations: %d" % funcalls

    if full_output:
        return x, fval, warnflag
    else:        
        return x


def zoom(a_lo, a_hi):
    pass

    

def line_search(f, fprime, xk, pk, gfk, args=(), c1=1e-4, c2=0.9, amax=50):
    """Minimize the function f(xk+alpha pk)

    Uses the line search algorithm of
    Wright and Nocedal in 'Numerical Optimization', 1999, pg. 59-60

    Outputs: (alpha0, gc, fc)
    """

    fc = 0
    gc = 0
    alpha0 = 1.0
    phi0  = apply(f,(xk,)+args)
    phi_a0 = apply(f,(xk+alpha0*pk,)+args)
    fc = fc + 2
    derphi0 = Num.dot(gfk,pk)
    derphi_a0 = Num.dot(apply(fprime,(xk+alpha0*pk,)+args),pk)
    gc = gc + 1

    # check to see if alpha0 = 1 satisfies Strong Wolfe conditions.
    if (phi_a0 <= phi0 + c1*alpha0*derphi0) \
       and (abs(derphi_a0) <= c2*abs(derphi0)):
        return alpha0, fc, gc

    alpha0 = 0
    alpha1 = 1
    phi_a1 = phi_a0
    phi_a0 = phi0

    i = 1
    while 1:
        if (phi_a1 > phi0 + c1*alpha1*derphi0) or \
           ((phi_a1 >= phi_a0) and (i > 1)):
            return zoom(alpha0, alpha1)

        derphi_a1 = Num.dot(apply(fprime,(xk+alpha1*pk,)+args),pk)
        gc = gc + 1
        if (abs(derphi_a1) <= -c2*derphi0):
            return alpha1

        if (derphi_a1 >= 0):
            return zoom(alpha1, alpha0)

        alpha2 = (amax-alpha1)*0.25 + alpha1
        i = i + 1
        alpha0 = alpha1
        alpha1 = alpha2
        phi_a0 = phi_a1
        phi_a1 = apply(f,(xk+alpha1*pk,)+args)

    

def line_search_BFGS(f, xk, pk, gfk, args=(), c1=1e-4, alpha0=1):
    """Minimize over alpha, the function f(xk+alpha pk)

    Uses the interpolation algorithm (Armiijo backtracking) as suggested by
    Wright and Nocedal in 'Numerical Optimization', 1999, pg. 56-57
    
    Outputs: (alpha, fc, gc)
    """

    fc = 0
    phi0 = apply(f,(xk,)+args)               # compute f(xk)
    phi_a0 = apply(f,(xk+alpha0*pk,)+args)     # compute f
    fc = fc + 2
    derphi0 = Num.dot(gfk,pk)

    if (phi_a0 <= phi0 + c1*alpha0*derphi0):
        return alpha0, fc, 0

    # Otherwise compute the minimizer of a quadratic interpolant:

    alpha1 = -(derphi0) * alpha0**2 / 2.0 / (phi_a0 - phi0 - derphi0 * alpha0)
    phi_a1 = apply(f,(xk+alpha1*pk,)+args)
    fc = fc + 1

    if (phi_a1 <= phi0 + c1*alpha1*derphi0):
        return alpha1, fc, 0

    # Otherwise loop with cubic interpolation until we find an alpha which 
    # satifies the first Wolfe condition (since we are backtracking, we will
    # assume that the value of alpha is not too small and satisfies the second
    # condition.

    while 1:       # we are assuming pk is a descent direction
        factor = alpha0**2 * alpha1**2 * (alpha1-alpha0)
        a = alpha0**2 * (phi_a1 - phi0 - derphi0*alpha1) - \
            alpha1**2 * (phi_a0 - phi0 - derphi0*alpha0)
        a = a / factor
        b = -alpha0**3 * (phi_a1 - phi0 - derphi0*alpha1) + \
            alpha1**3 * (phi_a0 - phi0 - derphi0*alpha0)
        b = b / factor

        alpha2 = (-b + Num.sqrt(abs(b**2 - 3 * a * derphi0))) / (3.0*a)
        phi_a2 = apply(f,(xk+alpha2*pk,)+args)
        fc = fc + 1

        if (phi_a2 <= phi0 + c1*alpha2*derphi0):
            return alpha2, fc, 0

        if (alpha1 - alpha2) > alpha1 / 2.0 or (1 - alpha2/alpha1) < 0.96:
            alpha2 = alpha1 / 2.0

        alpha0 = alpha1
        alpha1 = alpha2
        phi_a0 = phi_a1
        phi_a1 = phi_a2

epsilon = 1e-8

def approx_fprime(xk,f,*args):
    f0 = apply(f,(xk,)+args)
    grad = Num.zeros((len(xk),),'d')
    ei = Num.zeros((len(xk),),'d')
    for k in range(len(xk)):
        ei[k] = 1.0
        grad[k] = (apply(f,(xk+epsilon*ei,)+args) - f0)/epsilon
        ei[k] = 0.0
    return grad

def approx_fhess_p(x0,p,fprime,*args):
    f2 = apply(fprime,(x0+epsilon*p,)+args)
    f1 = apply(fprime,(x0,)+args)
    return (f2 - f1)/epsilon


def fminBFGS(f, x0, fprime=None, args=(), avegtol=1e-5, maxiter=None, 
             full_output=0, printmessg=1):
    """Minimize a function using the BFGS algorithm.

    Description:

      Optimize the function, f, whose gradient is given by fprime using the
      quasi-Newton method of Broyden, Fletcher, Goldfarb, and Shanno (BFGS)
      See Wright, and Nocedal 'Numerical Optimization', 1999, pg. 198.

    Inputs:

      f -- the Python function or method to be minimized.
      x0 -- the initial guess for the minimizer.
      fprime -- a function to compute the gradient of f.
      args -- extra arguments to f and fprime.

    Outputs: (xopt, {fopt, func_calls, grad_calls, warnflag})

      xopt -- the minimizer of f.

      fopt -- the value of f(xopt).
      func_calls -- the number of function_calls.
      grad_calls -- the number of gradient calls.
      warnflag -- an integer warning flag:
                  1 : 'Maximum number of iterations exceeded.'

    Additional Inputs:

      avegtol -- the minimum occurs when fprime(xopt)==0.  This specifies how
                 close to zero the average magnitude of fprime(xopt) needs
                 to be.
      maxiter -- the maximum number of iterations.
      full_output -- if non-zero then return fopt, func_calls, grad_calls,
                     and warnflag in addition to xopt.
      printmessg -- print convergence message if non-zero.                
      """
    app_fprime = 0
    if fprime is None:
        app_fprime = 1

    x0 = Num.asarray(x0)
    if maxiter is None:
        maxiter = len(x0)*200
    func_calls = 0
    grad_calls = 0
    k = 0
    N = len(x0)
    gtol = N*avegtol
    I = MLab.eye(N)
    Hk = I

    if app_fprime:
        gfk = apply(approx_fprime,(x0,f)+args)
        func_calls = func_calls + len(x0) + 1
    else:
        gfk = apply(fprime,(x0,)+args)
        grad_calls = grad_calls + 1
    xk = x0
    sk = [2*gtol]
    while (Num.add.reduce(abs(gfk)) > gtol) and (k < maxiter):
        pk = -Num.dot(Hk,gfk)
        alpha_k, fc, gc = line_search_BFGS(f,xk,pk,gfk,args)
        func_calls = func_calls + fc
        xkp1 = xk + alpha_k * pk
        sk = xkp1 - xk
        xk = xkp1
        if app_fprime:
            gfkp1 = apply(approx_fprime,(xkp1,f)+args)
            func_calls = func_calls + gc + len(x0) + 1
        else:
            gfkp1 = apply(fprime,(xkp1,)+args)
            grad_calls = grad_calls + gc + 1

        yk = gfkp1 - gfk
        k = k + 1

        rhok = 1 / Num.dot(yk,sk)
        A1 = I - sk[:,Num.NewAxis] * yk[Num.NewAxis,:] * rhok
        A2 = I - yk[:,Num.NewAxis] * sk[Num.NewAxis,:] * rhok
        Hk = Num.dot(A1,Num.dot(Hk,A2)) + rhok * sk[:,Num.NewAxis] \
                                               * sk[Num.NewAxis,:]
        gfk = gfkp1


    if printmessg or full_output:
        fval = apply(f,(xk,)+args)
    if k >= maxiter:
        warnflag = 1
        if printmessg:
            print "Warning: Maximum number of iterations has been exceeded"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls
    else:
        warnflag = 0
        if printmessg:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls

    if full_output:
        return xk, fval, func_calls, grad_calls, warnflag
    else:        
        return xk


def fminNCG(f, x0, fprime, fhess_p=None, fhess=None, args=(), avextol=1e-5,
            maxiter=None, full_output=0, printmessg=1):
    """Minimze 
  Description:

    Optimize the function, f, whose gradient is given by fprime using the
    Newton-CG method.  fhess_p must compute the hessian times an arbitrary
    vector. If it is not given, finite-differences on fprime are used to
    compute it. See Wright, and Nocedal 'Numerical Optimization', 1999,
    pg. 140.

  Inputs:

    f -- the Python function or method to be minimized.
    x0 -- the initial guess for the minimizer.
    fprime -- a function to compute the gradient of f: fprime(x, *args)
    fhess_p -- a function to compute the Hessian of f times an
               arbitrary vector: fhess_p (x, p, *args)
    fhess -- a function to compute the Hessian matrix of f.
    args -- extra arguments for f, fprime, fhess_p, and fhess (the same
            set of extra arguments is supplied to all of these functions).

  Outputs: (xopt, {fopt, fcalls, gcalls, hcalls, warnflag})

    xopt -- the minimizer of f
    
    fopt -- the value of the function at xopt: fopt = f(xopt)
    fcalls -- the number of function calls.
    gcalls -- the number of gradient calls.
    hcalls -- the number of hessian calls.
    warnflag -- algorithm warnings:
                1 : 'Maximum number of iterations exceeded.'

  Additional Inputs:

    avextol -- Convergence is assumed when the average relative error in
               the minimizer falls below this amount.  
    maxiter -- Maximum number of iterations to allow.
    full_output -- If non-zero return the optional outputs.
    printmessg -- If non-zero print convergence message.                
    
  Remarks:

    Only one of fhess_p or fhess need be given.  If fhess is provided,
    then fhess_p will be ignored.  If neither fhess nor fhess_p is
    provided, then the hessian product will be approximated using finite
    differences on fprime.

    """

    x0 = Num.asarray(x0)
    fcalls = 0
    gcalls = 0
    hcalls = 0
    if maxiter is None:
        maxiter = len(x0)*200
    
    xtol = len(x0)*avextol
    update = [2*xtol]
    xk = x0
    k = 0
    while (Num.add.reduce(abs(update)) > xtol) and (k < maxiter):
        # Compute a search direction pk by applying the CG method to
        #  del2 f(xk) p = - grad f(xk) starting from 0.
        b = -apply(fprime,(xk,)+args)
        gcalls = gcalls + 1
        maggrad = Num.add.reduce(abs(b))
        eta = min([0.5,Num.sqrt(maggrad)])
        termcond = eta * maggrad
        xsupi = 0
        ri = -b
        psupi = -ri
        i = 0
        dri0 = Num.dot(ri,ri)

        if fhess is not None:             # you want to compute hessian once.
            A = apply(fhess,(xk,)+args)
            hcalls = hcalls + 1

        while Num.add.reduce(abs(ri)) > termcond:
            if fhess is None:
                if fhess_p is None:
                    Ap = apply(approx_fhess_p,(xk,psupi,fprime)+args)
                    gcalls = gcalls + 2
                else:
                    Ap = apply(fhess_p,(xk,psupi)+args)
                    hcalls = hcalls + 1
            else:
                Ap = Num.dot(A,psupi)
            # check curvature
            curv = Num.dot(psupi,Ap)
            if (curv <= 0):
                if (i > 0):
                    break
                else:
                    xsupi = xsupi + dri0/curv * psupi
                    break
            alphai = dri0 / curv
            xsupi = xsupi + alphai * psupi
            ri = ri + alphai * Ap
            dri1 = Num.dot(ri,ri)
            betai = dri1 / dri0
            psupi = -ri + betai * psupi
            i = i + 1
            dri0 = dri1          # update Num.dot(ri,ri) for next time.
    
        pk = xsupi  # search direction is solution to system.
        gfk = -b    # gradient at xk
        alphak, fc, gc = line_search_BFGS(f,xk,pk,gfk,args)
        fcalls = fcalls + fc
        gcalls = gcalls + gc

        update = alphak * pk
        xk = xk + update
        k = k + 1

    if printmessg or full_output:
        fval = apply(f,(xk,)+args)
    if k >= maxiter:
        warnflag = 1
        if printmessg:
            print "Warning: Maximum number of iterations has been exceeded"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % fcalls
            print "         Gradient evaluations: %d" % gcalls
            print "         Hessian evaluations: %d" % hcalls
    else:
        warnflag = 0
        if printmessg:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % fcalls
            print "         Gradient evaluations: %d" % gcalls
            print "         Hessian evaluations: %d" % hcalls
            
    if full_output:
        return xk, fval, fcalls, gcalls, hcalls, warnflag
    else:        
        return xk

def fminbound(func, x1, x2, args=(), xtol=1e-5, maxfun=500, 
              full_output=0, disp=1):
    """Bounded minimization for scalar functions.

    Description:

      Finds a local minimizer of the scalar function func in the interval
      x1 < xopt < x2.

    Inputs:

      func -- the function to be minimized (must accept scalar input and return
              scalar output).
      x1, x2 -- the optimization bounds.
      args -- extra arguments to pass to function.
      xtol -- the convergence tolerance.
      maxfun -- maximum function evaluations.
      full_output -- Non-zero to return optional outputs. 
      disp -- Non-zero to print messages.
              0 : no message printing.
              1 : non-convergence notification messages only.
              2 : print a message on convergence too.
              3 : print iteration results. 


    Outputs: (xopt, {fval, ierr, numfunc})

      xopt -- The minimizer of the function over the interval.
      fval -- The function value at the minimum point.
      ierr -- An error flag (1 if converged, 0 if maximum number of
              function calls reached).
      numfunc -- The number of function calls.
    """

    if x1 > x2:
        raise ValueError, "The lower bound exceeds the upper bound."

    flag = 1
    header = ' Func-count     x          f(x)          Procedure'
    step='       initial'

    sqrt_eps = sqrt(2.2e-16)
    c0 = 0.5*(3.0-sqrt(5.0))
    a, b = x1, x2
    fulc = a + c0*(b-a)
    nfc, xf = fulc, fulc
    rat = e = 0.0
    x = xf
    fx = func(x,*args)
    num = 1
    fmin_data = (1, xf, fx)

    ffulc = fnfc = fx
    xm = 0.5*(a+b)
    tol1 = sqrt_eps*abs(xf) + xtol / 3.0
    tol2 = 2.0*tol1

    if disp > 2:
        print (" ")
        print (header)
        print "%5.0f   %12.6g %12.6g %s" % (fmin_data + (step,))

    # Main loop

    while ( abs(xf-xm) > (tol2 - 0.5*(b-a)) ):
        golden = 1
        # Check for parabolic fit
        if abs(e) > tol1:
            golden = 0
            r = (xf-nfc)*(fx-ffulc)
            q = (xf-fulc)*(fx-fnfc)
            p = (xf-fulc)*q - (xf-nfc)*r
            q = 2.0*(q-r)
            if q > 0.0: p = -p
            q = abs(q)
            r = e
            e = rat

            # Check for acceptability of parabola
            if ( (abs(p) < abs(0.5*q*r)) and (p > q*(a-xf)) and (p < q*(b-xf))):
                rat = (p+0.0) / q;
                x = xf + rat
                step = '       parabolic'

                if ((x-a) < tol2) or ((b-x) < tol2):
                    si = Numeric.sign(xm-xf) + ((xm-xf)==0)
                    rat = tol1*si
            else:      # do a golden section step
                golden = 1

        if golden:  # Do a golden-section step
            if xf >= xm:
                e=a-xf
            else:
                e=b-xf
            rat = c0*e
            step = '       golden'

        si = Numeric.sign(rat) + (rat == 0)
        x = xf + si*max([abs(rat), tol1])
        fu = func(x,*args)
        num += 1
        fmin_data = (num, x, fu)
        if disp > 2:
            print "%5.0f   %12.6g %12.6g %s" % (fmin_data + (step,))
                
        if fu <= fx:
            if x >= xf:
                a = xf
            else:
                b = xf
            fulc, ffulc = nfc, fnfc
            nfc, fnfc = xf, fx
            xf, fx = x, fu
        else:
            if x < xf:
                a = x
            else:
                b = x
            if (fu <= fnfc) or (nfc == xf):
                fulc, ffulc = nfc, fnfc
                nfc, fnfc = x, fu
            elif (fu <= ffulc) or (fulc == xf) or (fulc == nfc):
                fulc, ffulc = x, fu

        xm = 0.5*(a+b)
        tol1 = sqrt_eps*abs(xf) + xtol/3.0
        tol2 = 2.0*tol1

        if num >= maxfun:
            flag = 0
            fval = fx
            if disp > 0:
                _endprint(x, flag, fval, maxfun, tol, disp)
            if full_output:
                return xf, fval, flag, num
            else:
                return xf

    fval = fx
    if disp > 0:
        _endprint(x, flag, fval, maxfun, xtol, disp)
    
    if full_output:
        return xf, fval, flag, num
    else:
        return xf


def _endprint(x, flag, fval, maxfun, xtol, disp):
    if flag == 1:
        if disp > 1:
            print "\nOptimization terminated successfully;\n the returned value" + \
                  " satisfies the termination criteria\n (using xtol = ", xtol, ")"
    if flag == 0:
        print "\nMaximum number of function evaluations exceeded --- increase maxfun argument.\n"
    return
            
if __name__ == "__main__":
    import string
    import time

    
    times = []
    algor = []
    x0 = [0.8,1.2,0.7]
    start = time.time()
    x = fmin(rosen,x0)
    print x
    times.append(time.time() - start)
    algor.append('Nelder-Mead Simplex\t')

    start = time.time()
    x = fminBFGS(rosen, x0, fprime=rosen_der, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('BFGS Quasi-Newton\t')

    start = time.time()
    x = fminBFGS(rosen, x0, avegtol=1e-4, maxiter=100)
    print x
    times.append(time.time() - start)
    algor.append('BFGS without gradient\t')


    start = time.time()
    x = fminNCG(rosen, x0, rosen_der, fhess_p=rosen_hess_p, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('Newton-CG with hessian product')
    

    start = time.time()
    x = fminNCG(rosen, x0, rosen_der, fhess=rosen_hess, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('Newton-CG with full hessian')

    print "\nMinimizing the Rosenbrock function of order 3\n"
    print " Algorithm \t\t\t       Seconds"
    print "===========\t\t\t      ========="
    for k in range(len(algor)):
        print algor[k], "\t -- ", times[k]
        















