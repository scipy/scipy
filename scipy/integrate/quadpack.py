# Author: Travis Oliphant 2001

import _quadpack
import sys
import numpy
from numpy import inf, Inf

__all__ = ['quad', 'dblquad', 'tplquad', 'quad_explain']


error = _quadpack.error

def quad_explain(output=sys.stdout):
    """
    Print extra information about integrate.quad() parameters and returns.

    Parameters
    ----------
    output : instance with "write" method
        Information about `quad` is passed to ``output.write()``.
        Default is ``sys.stdout``.

    Returns
    -------
    None

    """
    output.write("""
Extra information for quad() inputs and outputs:

  If full_output is non-zero, then the third output argument (infodict)
  is a dictionary with entries as tabulated below.  For infinite limits, the
  range is transformed to (0,1) and the optional outputs are given with
  respect to this transformed range.  Let M be the input argument limit and
  let K be infodict['last'].  The entries are:

  'neval' : The number of function evaluations.
  'last'  : The number, K, of subintervals produced in the subdivision process.
  'alist' : A rank-1 array of length M, the first K elements of which are the
            left end points of the subintervals in the partition of the
            integration range.
  'blist' : A rank-1 array of length M, the first K elements of which are the
            right end points of the subintervals.
  'rlist' : A rank-1 array of length M, the first K elements of which are the
            integral approximations on the subintervals.
  'elist' : A rank-1 array of length M, the first K elements of which are the
            moduli of the absolute error estimates on the subintervals.
  'iord'  : A rank-1 integer array of length M, the first L elements of
            which are pointers to the error estimates over the subintervals
            with L=K if K<=M/2+2 or L=M+1-K otherwise. Let I be the sequence
            infodict['iord'] and let E be the sequence infodict['elist'].
            Then E[I[1]], ..., E[I[L]] forms a decreasing sequence.

  If the input argument points is provided (i.e. it is not None), the
  following additional outputs are placed in the output dictionary.  Assume the
  points sequence is of length P.

  'pts'   : A rank-1 array of length P+2 containing the integration limits
            and the break points of the intervals in ascending order.
            This is an array giving the subintervals over which integration
            will occur.
  'level' : A rank-1 integer array of length M (=limit), containing the
            subdivision levels of the subintervals, i.e., if (aa,bb) is a
            subinterval of (pts[1], pts[2]) where pts[0] and pts[2] are
            adjacent elements of infodict['pts'], then (aa,bb) has level l if
            |bb-aa|=|pts[2]-pts[1]| * 2**(-l).
  'ndin'  : A rank-1 integer array of length P+2.  After the first integration
            over the intervals (pts[1], pts[2]), the error estimates over some
            of the intervals may have been increased artificially in order to
            put their subdivision forward.  This array has ones in slots
            corresponding to the subintervals for which this happens.

Weighting the integrand:

  The input variables, weight and wvar, are used to weight the integrand by
  a select list of functions.  Different integration methods are used
  to compute the integral with these weighting functions.  The possible values
  of weight and the corresponding weighting functions are.

  'cos'     : cos(w*x)                            : wvar = w
  'sin'     : sin(w*x)                            : wvar = w
  'alg'     : g(x) = ((x-a)**alpha)*((b-x)**beta) : wvar = (alpha, beta)
  'alg-loga': g(x)*log(x-a)                       : wvar = (alpha, beta)
  'alg-logb': g(x)*log(b-x)                       : wvar = (alpha, beta)
  'alg-log' : g(x)*log(x-a)*log(b-x)              : wvar = (alpha, beta)
  'cauchy'  : 1/(x-c)                             : wvar = c

  wvar holds the parameter w, (alpha, beta), or c depending on the weight
  selected.  In these expressions, a and b are the integration limits.

  For the 'cos' and 'sin' weighting, additional inputs and outputs are
  available.

  For finite integration limits, the integration is performed using a
  Clenshaw-Curtis method which uses Chebyshev moments.  For repeated
  calculations, these moments are saved in the output dictionary:

  'momcom' : The maximum level of Chebyshev moments that have been computed,
             i.e., if M_c is infodict['momcom'] then the moments have been
             computed for intervals of length |b-a|* 2**(-l), l=0,1,...,M_c.
  'nnlog'  : A rank-1 integer array of length M(=limit), containing the
             subdivision levels of the subintervals, i.e., an element of this
             array is equal to l if the corresponding subinterval is
             |b-a|* 2**(-l).
  'chebmo' : A rank-2 array of shape (25, maxp1) containing the computed
             Chebyshev moments.  These can be passed on to an integration
             over the same interval by passing this array as the second
             element of the sequence wopts and passing infodict['momcom'] as
             the first element.

  If one of the integration limits is infinite, then a Fourier integral is
  computed (assuming w neq 0).  If full_output is 1 and a numerical error
  is encountered, besides the error message attached to the output tuple,
  a dictionary is also appended to the output tuple which translates the
  error codes in the array info['ierlst'] to English messages.  The output
  information dictionary contains the following entries instead of 'last',
  'alist', 'blist', 'rlist', and 'elist':

  'lst'    : The number of subintervals needed for the integration (call it K_f).
  'rslst'  : A rank-1 array of length M_f=limlst, whose first K_f elements
             contain the integral contribution over the interval (a+(k-1)c,
             a+kc) where c = (2*floor(|w|) + 1) * pi / |w| and k=1,2,...,K_f.
  'erlst'  : A rank-1 array of length M_f containing the error estimate
             corresponding to the interval in the same position in
             infodict['rslist'].
  'ierlst' : A rank-1 integer array of length M_f containing an error flag
             corresponding to the interval in the same position in
             infodict['rslist'].  See the explanation dictionary (last entry
             in the output tuple) for the meaning of the codes.
""")
    return


def quad(func, a, b, args=(), full_output=0, epsabs=1.49e-8, epsrel=1.49e-8,
         limit=50, points=None, weight=None, wvar=None, wopts=None, maxp1=50,
         limlst=50):
    """
    Compute a definite integral.

    Integrate func from a to b (possibly infinite interval) using a technique
    from the Fortran library QUADPACK.

    If func takes many arguments, it is integrated along the axis corresponding
    to the first argument. Use the keyword argument `args` to pass the other
    arguments.

    Run scipy.integrate.quad_explain() for more information on the
    more esoteric inputs and outputs.

    Parameters
    ----------

    func : function
        A Python function or method to integrate.
    a : float
        Lower limit of integration (use -scipy.integrate.Inf for -infinity).
    b : float
        Upper limit of integration (use scipy.integrate.Inf for +infinity).
    args : tuple, optional
        extra arguments to pass to func
    full_output : int
        Non-zero to return a dictionary of integration information.
        If non-zero, warning messages are also suppressed and the
        message is appended to the output tuple.

    Returns
    -------

    y : float
        The integral of func from a to b.
    abserr : float
        an estimate of the absolute error in the result.

    infodict : dict
        a dictionary containing additional information.
        Run scipy.integrate.quad_explain() for more information.
    message :
        a convergence message.
    explain :
        appended only with 'cos' or 'sin' weighting and infinite
        integration limits, it contains an explanation of the codes in
        infodict['ierlst']

    Other Parameters
    ----------------
    epsabs :
        absolute error tolerance.
    epsrel :
        relative error tolerance.
    limit :
        an upper bound on the number of subintervals used in the adaptive
        algorithm.
    points :
        a sequence of break points in the bounded integration interval
        where local difficulties of the integrand may occur (e.g.,
        singularities, discontinuities). The sequence does not have
        to be sorted.
    weight :
        string indicating weighting function.
    wvar :
        variables for use with weighting functions.
    limlst :
        Upper bound on the number of cylces (>=3) for use with a sinusoidal
        weighting and an infinite end-point.
    wopts :
        Optional input for reusing Chebyshev moments.
    maxp1 :
        An upper bound on the number of Chebyshev moments.

    See Also
    --------
    dblquad, tplquad : double and triple integrals
    fixed_quad : fixed-order Gaussian quadrature
    quadrature : adaptive Gaussian quadrature
    odeint, ode : ODE integrators
    simps, trapz, romb : integrators for sampled data
    scipy.special : for coefficients and roots of orthogonal polynomials

    Examples
    --------
    Calculate :math:`\\int^4_0 x^2 dx` and compare with an analytic result

    >>> from scipy import integrate
    >>> x2 = lambda x: x**2
    >>> integrate.quad(x2,0.,4.)
    (21.333333333333332, 2.3684757858670003e-13)
    >> print 4.**3/3
    21.3333333333

    Calculate :math:`\\int^\\infty_0 e^{-x} dx`

    >>> invexp = lambda x: exp(-x)
    >>> integrate.quad(invexp,0,inf)
    (0.99999999999999989, 5.8426061711142159e-11)

    >>> f = lambda x,a : a*x
    >>> y, err = integrate.quad(f, 0, 1, args=(1,))
    >>> y
    0.5
    >>> y, err = integrate.quad(f, 0, 1, args=(3,))
    >>> y
    1.5

    """
    if type(args) != type(()): args = (args,)
    if (weight is None):
        retval = _quad(func,a,b,args,full_output,epsabs,epsrel,limit,points)
    else:
        retval = _quad_weight(func,a,b,args,full_output,epsabs,epsrel,limlst,limit,maxp1,weight,wvar,wopts)

    ier = retval[-1]
    if ier == 0:
        return retval[:-1]

    msgs = {80: "A Python error occurred possibly while calling the function.",
             1: "The maximum number of subdivisions (%d) has been achieved.\n  If increasing the limit yields no improvement it is advised to analyze \n  the integrand in order to determine the difficulties.  If the position of a \n  local difficulty can be determined (singularity, discontinuity) one will \n  probably gain from splitting up the interval and calling the integrator \n  on the subranges.  Perhaps a special-purpose integrator should be used." % limit,
             2: "The occurrence of roundoff error is detected, which prevents \n  the requested tolerance from being achieved.  The error may be \n  underestimated.",
             3: "Extremely bad integrand behavior occurs at some points of the\n  integration interval.",
             4: "The algorithm does not converge.  Roundoff error is detected\n  in the extrapolation table.  It is assumed that the requested tolerance\n  cannot be achieved, and that the returned result (if full_output = 1) is \n  the best which can be obtained.",
             5: "The integral is probably divergent, or slowly convergent.",
             6: "The input is invalid.",
             7: "Abnormal termination of the routine.  The estimates for result\n  and error are less reliable.  It is assumed that the requested accuracy\n  has not been achieved.",
            'unknown': "Unknown error."}

    if weight in ['cos','sin'] and (b == Inf or a == -Inf):
        msgs[1] = "The maximum number of cycles allowed has been achieved., e.e.\n  of subintervals (a+(k-1)c, a+kc) where c = (2*int(abs(omega)+1))\n  *pi/abs(omega), for k = 1, 2, ..., lst.  One can allow more cycles by increasing the value of limlst.  Look at info['ierlst'] with full_output=1."
        msgs[4] = "The extrapolation table constructed for convergence acceleration\n  of the series formed by the integral contributions over the cycles, \n  does not converge to within the requested accuracy.  Look at \n  info['ierlst'] with full_output=1."
        msgs[7] = "Bad integrand behavior occurs within one or more of the cycles.\n  Location and type of the difficulty involved can be determined from \n  the vector info['ierlist'] obtained with full_output=1."
        explain = {1: "The maximum number of subdivisions (= limit) has been \n  achieved on this cycle.",
                   2: "The occurrence of roundoff error is detected and prevents\n  the tolerance imposed on this cycle from being achieved.",
                   3: "Extremely bad integrand behavior occurs at some points of\n  this cycle.",
                   4: "The integral over this cycle does not converge (to within the required accuracy) due ot roundoff in the extrapolation procedure invoked on this cycle.  It is assumed that the result on this interval is the best which can be obtained.",
                   5: "The integral over this cycle is probably divergent or slowly convergent."}

    try:
        msg = msgs[ier]
    except KeyError:
        msg = msgs['unknown']

    if ier in [1,2,3,4,5,7]:
        if full_output:
            if weight in ['cos','sin'] and (b == Inf or a == Inf):
                return retval[:-1] + (msg, explain)
            else:
                return retval[:-1] + (msg,)
        else:
            import warnings
            warnings.warn(msg)
            return retval[:-1]
    else:
        raise ValueError(msg)


def _quad(func,a,b,args,full_output,epsabs,epsrel,limit,points):
    infbounds = 0
    if (b != Inf and a != -Inf):
        pass   # standard integration
    elif (b == Inf and a != -Inf):
        infbounds = 1
        bound = a
    elif (b == Inf and a == -Inf):
        infbounds = 2
        bound = 0     # ignored
    elif (b != Inf and a == -Inf):
        infbounds = -1
        bound = b
    else:
        raise RuntimeError("Infinity comparisons don't work for you.")

    if points is None:
        if infbounds == 0:
            return _quadpack._qagse(func,a,b,args,full_output,epsabs,epsrel,limit)
        else:
            return _quadpack._qagie(func,bound,infbounds,args,full_output,epsabs,epsrel,limit)
    else:
        if infbounds !=0:
            raise ValueError("Infinity inputs cannot be used with break points.")
        else:
            nl = len(points)
            the_points = numpy.zeros((nl+2,), float)
            the_points[:nl] = points
            return _quadpack._qagpe(func,a,b,the_points,args,full_output,epsabs,epsrel,limit)


def _quad_weight(func,a,b,args,full_output,epsabs,epsrel,limlst,limit,maxp1,weight,wvar,wopts):

    if weight not in ['cos','sin','alg','alg-loga','alg-logb','alg-log','cauchy']:
        raise ValueError("%s not a recognized weighting function." % weight)

    strdict = {'cos':1,'sin':2,'alg':1,'alg-loga':2,'alg-logb':3,'alg-log':4}

    if weight in ['cos','sin']:
        integr = strdict[weight]
        if (b != Inf and a != -Inf):  # finite limits
            if wopts is None:         # no precomputed chebyshev moments
                return _quadpack._qawoe(func,a,b,wvar,integr,args,full_output,epsabs,epsrel,limit,maxp1,1)
            else:                     # precomputed chebyshev moments
                momcom = wopts[0]
                chebcom = wopts[1]
                return _quadpack._qawoe(func,a,b,wvar,integr,args,full_output,epsabs,epsrel,limit,maxp1,2,momcom,chebcom)

        elif (b == Inf and a != -Inf):
            return _quadpack._qawfe(func,a,wvar,integr,args,full_output,epsabs,limlst,limit,maxp1)
        elif (b != Inf and a == -Inf):  # remap function and interval
            if weight == 'cos':
                def thefunc(x,*myargs):
                    y = -x
                    func = myargs[0]
                    myargs = (y,) + myargs[1:]
                    return apply(func,myargs)
            else:
                def thefunc(x,*myargs):
                    y = -x
                    func = myargs[0]
                    myargs = (y,) + myargs[1:]
                    return -apply(func,myargs)
            args = (func,) + args
            return _quadpack._qawfe(thefunc,-b,wvar,integr,args,full_output,epsabs,limlst,limit,maxp1)
        else:
            raise ValueError("Cannot integrate with this weight from -Inf to +Inf.")
    else:
        if a in [-Inf,Inf] or b in [-Inf,Inf]:
            raise ValueError("Cannot integrate with this weight over an infinite interval.")

        if weight[:3] == 'alg':
            integr = strdict[weight]
            return _quadpack._qawse(func,a,b,wvar,integr,args,full_output,epsabs,epsrel,limit)
        else:  # weight == 'cauchy'
            return _quadpack._qawce(func,a,b,wvar,args,full_output,epsabs,epsrel,limit)

def _infunc(x,func,gfun,hfun,more_args):
    a = gfun(x)
    b = hfun(x)
    myargs = (x,) + more_args
    return quad(func,a,b,args=myargs)[0]

def dblquad(func, a, b, gfun, hfun, args=(), epsabs=1.49e-8, epsrel=1.49e-8):
    """
    Compute a double integral.

    Return the double (definite) integral of func(y,x) from x=a..b and
    y=gfun(x)..hfun(x).

    Parameters
    -----------
    func : callable
        A Python function or method of at least two variables: y must be the
        first argument and x the second argument.
    (a,b) : tuple
        The limits of integration in x: a < b
    gfun : callable
        The lower boundary curve in y which is a function taking a single
        floating point argument (x) and returning a floating point result: a
        lambda function can be useful here.
    hfun : callable
        The upper boundary curve in y (same requirements as `gfun`).
    args : sequence, optional
        Extra arguments to pass to `func2d`.
    epsabs : float, optional
        Absolute tolerance passed directly to the inner 1-D quadrature
        integration. Default is 1.49e-8.
    epsrel : float
        Relative tolerance of the inner 1-D integrals. Default is 1.49e-8.

    Returns
    -------
    y : float
        The resultant integral.
    abserr : float
        An estimate of the error.

    See also
    --------
    quad : single integral
    tplquad : triple integral
    fixed_quad : fixed-order Gaussian quadrature
    quadrature : adaptive Gaussian quadrature
    odeint, ode : ODE integrators
    simps, trapz, romb : integrators for sampled data
    scipy.special : for coefficients and roots of orthogonal polynomials

    """
    return quad(_infunc,a,b,(func,gfun,hfun,args),epsabs=epsabs,epsrel=epsrel)

def _infunc2(y,x,func,qfun,rfun,more_args):
    a2 = qfun(x,y)
    b2 = rfun(x,y)
    myargs = (y,x) + more_args
    return quad(func,a2,b2,args=myargs)[0]

def tplquad(func, a, b, gfun, hfun, qfun, rfun, args=(), epsabs=1.49e-8,
            epsrel=1.49e-8):
    """
    Compute a triple (definite) integral.

    Return the triple integral of func(z, y, x) from
    x=a..b, y=gfun(x)..hfun(x), and z=qfun(x,y)..rfun(x,y)

    Parameters
    ----------
    func : function
        A Python function or method of at least three variables in the
        order (z, y, x).
    (a,b) : tuple
        The limits of integration in x: a < b
    gfun : function
        The lower boundary curve in y which is a function taking a single
        floating point argument (x) and returning a floating point result:
        a lambda function can be useful here.
    hfun : function
        The upper boundary curve in y (same requirements as gfun).
    qfun : function
        The lower boundary surface in z.  It must be a function that takes
        two floats in the order (x, y) and returns a float.
    rfun : function
        The upper boundary surface in z. (Same requirements as qfun.)
    args : Arguments
        Extra arguments to pass to func3d.
    epsabs : float, optional
        Absolute tolerance passed directly to the innermost 1-D quadrature
        integration. Default is 1.49e-8.
    epsrel : float, optional
        Relative tolerance of the innermost 1-D integrals. Default is 1.49e-8.

    Returns
    -------
    y : float
        The resultant integral.
    abserr : float
        An estimate of the error.

    See Also
    --------
    quad: Adaptive quadrature using QUADPACK
    quadrature: Adaptive Gaussian quadrature
    fixed_quad: Fixed-order Gaussian quadrature
    dblquad: Double integrals
    romb: Integrators for sampled data
    trapz: Integrators for sampled data
    simps: Integrators for sampled data
    ode: ODE integrators
    odeint: ODE integrators
    scipy.special: For coefficients and roots of orthogonal polynomials

    """
    return dblquad(_infunc2,a,b,gfun,hfun,(func,qfun,rfun,args),epsabs=epsabs,epsrel=epsrel)
