__version__ = "$Revision$"[10:-1]
import _minpack
from common_routines import *

def fsolve(func,x0,args=(),Dfun=None,full_output=0,col_deriv=0,xtol=1.49012e-8,maxfev=0,band=None,epsfcn=0.0,factor=100,diag=None):
    """[x,infodict,ier,msg] = fsolve(func, x0, args=(), Dfun=None, full_output=0,
                      col_deriv=0,xtol=1.49012e-8, maxfev=0, band=None, epsfcn=0.0,
                      factor=100, diag=None)

    return the roots of the N (non-linear) equations defined by func(x)=0
    given a starting estimate x0.

    Optional Inputs:
      args         additional arguments to the function (single value or TUPLE)
      Dfun         the jacobian of the function (derivatives across rows)
      full_output  non-zero to return the optional
                     outputs (otherwise only x is returned)
        
    other optional inputs available (see documentation).
    """
    x0 = myasarray(x0)
    n = len(x0)
    if type(args) != type(()): args = (args,)
    check_func(func,x0,args,n,(n,))
    if Dfun == None:
        if band == None:
            ml,mu = -10,-10
        else:
            ml,mu = band[:2]
        if (maxfev == 0):
            maxfev = 200*(n+1)
        retval = _minpack._hybrd(func,x0,args,full_output,xtol,maxfev,ml,mu,epsfcn,factor,diag)
    else:
        check_func(Dfun,x0,args,n,(n,n))
        if (maxfev == 0):
            maxfev = 100*(n+1)
        retval = _minpack._hybrj(func,Dfun,x0,args,full_output,col_deriv,xtol,maxfev,factor,diag)

    errors = {0:["Improper input parameters were entered.",TypeError],
              1:["The solution converged.",None],
              2:["The number of calls to function has reached maxfev = %d." % maxfev, ValueError],
              3:["xtol=%f is too small, no further improvement in the approximate\n  solution is possible." % xtol, ValueError],
              4:["The iteration is not making good progress, as measured by the \n  improvement from the last five Jacobian evaluations.", ValueError],
              5:["The iteration is not making good progress, as measured by the \n  improvement from the last ten iterations.", ValueError],
              'unknown': ["An error occurred.", TypeError]}

    info = retval[-1]    # The FORTRAN return value
    if (info != 1 and not full_output):
        if info in [2,3,4,5]:
            print "Warning: " + errors[info][0]
        else:
            try:
                raise errors[info][1], errors[info][0]
            except KeyError:
                raise errors['unknown'][1], errors['unknown'][0]

    if n == 1:
        retval = (retval[0][0],) + retval[1:]

    if full_output:
        try:
            return retval + (errors[info][0],)  # Return all + the message
        except KeyError:
            return retval + (errors['unknown'][0],)
    else:
        return retval[0]


def leastsq(func,x0,args=(),Dfun=None,full_output=0,col_deriv=0,ftol=1.49012e-8,xtol=1.49012e-8,gtol=0.0,maxfev=0,epsfcn=0.0,factor=100,diag=None):
    """[x,infodict] = leastsq(func, x0, args=(), Dfun=None, full_output=0,
                             ftol=1.49012e-8, xtol=1.49012e-8, gtol=0.0,
                             maxfev=0, epsfcn=0.0, factor=100, diag=None)

    return the point which minimizes the sum of squares of M (non-linear)
    equations in N unknowns given a starting estimate, x0.

                    x = arg min(sum(func(y)**2))
                             y
             
    Optional Inputs:
      args         additional arguments to the function (single or TUPLE)
      Dfun         the jacobian of the function (derivatives across rows)
      full_output  non-zero to return the dictionary infodict with optional
                     outputs (otherwise only x is returned)
        
    other optional inputs available (see documentation).
    """
    x0 = myasarray(x0)
    n = len(x0)
    if type(args) != type(()): args = (args,)
    m = check_func(func,x0,args,n)[0]
    if Dfun == None:
        if (maxfev == 0):
            maxfev = 200*(n+1)
        retval = _minpack._lmdif(func,x0,args,full_output,ftol,xtol,gtol,maxfev,epsfcn,factor,diag)
    else:
        if col_deriv:
            check_func(func,x0,args,n,(n,m))
        else:
            check_func(func,x0,args,n,(m,n))
        if (maxfev == 0):
            maxfev = 100*(n+1)
        retval = _minpack._lmder(func,Dfun,x0,args,full_output,col_deriv,ftol,xtol,gtol,maxfev,factor,diag)

    errors = {0:["Improper input parameters.", TypeError],
              1:["Both actual and predicted relative reductions in the sum of squares\n  are at most %f" % ftol, None],
              2:["The relative error between two consecutive iterates is at most %f" % xtol, None],
              3:["Both actual and predicted relative reductions in the sum of squares\n  are at most %f and the relative error between two consecutive iterates is at \n  most %f" % (ftol,xtol), None],
              4:["The cosine of the angle between func(x) and any column of the\n  Jacobian is at most %f in absolute value", gtol, None],
              5:["Number of calls to function has reached maxfev = %d." % maxfev, ValueError],
              6:["ftol=%f is too small, no further reduction in the sum of squares\n  is possible.""" % ftol, ValueError],
              7:["xtol=%f is too small, no further improvement in the approximate\n  solution is possible." % xtol, ValueError],
              8:["gtol=%f is too small, func(x) is orthogonal to the columns of\n  the Jacobian to machine precision." % gtol, ValueError],
              'unknown':["Unknown error.", TypeError]}

    info = retval[-1]    # The FORTRAN return value

    if (info not in [1,2,3,4] and not full_output):
        if info in [5,6,7,8]:
            print "Warning: " + errors[info][0]
        else:
            try:
                raise errors[info][1], errors[info][0]
            except KeyError:
                raise errors['unknown'][1], errors['unknown'][0]

    if n == 1:
        retval = (retval[0][0],) + retval[1:]

    if full_output:
        return retval + (errors[info][0],)
    else:
        return retval[:-1] + (errors[info][0],)


def check_gradient(fcn,Dfcn,x0,col_deriv=0):
    """good,err = check_gradient(fun,Dfun,x0,col_deriv=0)"""

    x = myasarray(x0)
    n = len(x)
    x.shape = (n,)
    fvec = myasarray(fcn(x))
    if 1 not in fvec.shape:
        raise ValueError, "Function does not return a 1-D array."
    m = len(fvec)
    fvec.shape = (m,)
    ldfjac = m
    fjac = myasarray(Dfcn(x))
    fjac.shape = (m,n)
    if col_deriv == 0:
        fjac = transpose(fjac)

    xp = zeros((n,),Float64)
    err = zeros((m,),Float64)
    fvecp = None
    _minpack._chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,1,err)
    
    fvecp = myasarray(fcn(xp))
    fvecp.shape = (m,)
    _minpack._chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,2,err)
    
    good = (product(greater(err,0.5)))

    return (good,err)













