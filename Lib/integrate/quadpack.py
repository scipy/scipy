# Author: Travis Oliphant 2001

__all__ = ['quad', 'dblquad', 'tplquad', 'quad_explain', 'Inf','inf']
import _quadpack
from common_routines import *
import sys
import MLab

error = _quadpack.error


def quad_explain(output=sys.stdout):
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

  'lst' -- The number of subintervals needed for the integration (call it K_f).
  'rslst' -- A rank-1 array of length M_f=limlst, whose first K_f elements
             contain the integral contribution over the interval (a+(k-1)c,
             a+kc) where c = (2*floor(|w|) + 1) * pi / |w| and k=1,2,...,K_f.
  'erlst' -- A rank-1 array of length M_f containing the error estimate
             corresponding to the interval in the same position in
             infodict['rslist'].
  'ierlst' -- A rank-1 integer array of length M_f containing an error flag
              corresponding to the interval in the same position in
              infodict['rslist'].  See the explanation dictionary (last entry
              in the output tuple) for the meaning of the codes.
  """)
    return
    

Inf = MLab.array(1e308)**10
inf = Inf
	
def quad(func, a, b, args=(), full_output=0, epsabs=1.49e-8, epsrel=1.49e-8,
         limit=50, points=None, weight=None, wvar=None, wopts=None, maxp1=50,
         limlst=50):
    """Compute a definite integral.

  Description:
  
    Integrate func from a to b (possibly infinite interval) using a technique
    from the Fortran library QUADPACK.  Run scipy.integrate.quad_explain()
    for more information on the more esoteric inputs and outputs.

  Inputs:

    func -- a Python function or method to integrate.
    a -- lower limit of integration (use -scipy.integrate.Inf for -infinity).
    b -- upper limit of integration (use scipy.integrate.Inf for +infinity).
    args -- extra arguments to pass to func.
    full_output -- non-zero to return a dictionary of integration information.
                   If non-zero, warning messages are also suppressed and the
                   message is appended to the output tuple.

  Outputs: (y, abserr, {infodict, message, explain})

    y -- the integral of func from a to b.
    abserr -- an estimate of the absolute error in the result.
    
    infodict -- a dictionary containing additional information.
                Run scipy.integrate.quad_explain() for more information.
    message -- a convergence message.    
    explain -- appended only with 'cos' or 'sin' weighting and infinite
               integration limits, it contains an explanation of the codes in
               infodict['ierlst']    

  Additional Inputs:

    epsabs -- absolute error tolerance.
    epsrel -- relative error tolerance.
    limit -- an upper bound on the number of subintervals used in the adaptive
             algorithm.
    points -- a sequence of break points in the bounded integration interval
              where local difficulties of the integrand may occur (e.g.,
              singularities, discontinuities). The sequence does not have
              to be sorted.

         **
         ** Run scipy.integrate.quad_explain() for more information
         ** on the following inputs
         **
    weight -- string indicating weighting function.
    wvar -- variables for use with weighting functions.
    limlst -- Upper bound on the number of cylces (>=3) for use with a sinusoidal
              weighting and an infinite end-point.
    wopts -- Optional input for reusing Chebyshev moments.
    maxp1 -- An upper bound on the number of Chebyshev moments.

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
             2: "The ocurrence of roundoff error is detected, which prevents \n  the requested tolerance from being achieved.  The error may be \n  underestimated.",
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
            print "Warning: " + msg
            return retval[:-1]
    else:
        raise ValueError, msg


def _quad(func,a,b,args,full_output,epsabs,epsrel,limit,points):
    inf = 0
    if (b != Inf and a != -Inf):
        pass   # standard integration
    elif (b == Inf and a != -Inf):
        inf = 1
        bound = a
    elif (b == Inf and a == -Inf):
        inf = 2
        bound = 0     # ignored
    elif (b != Inf and a == -Inf):
        inf = -1
        bound = b
    else:
        raise RunTimeError, "Infinity comparisons don't work for you."

    if points is None:
        if inf == 0:
            return _quadpack._qagse(func,a,b,args,full_output,epsabs,epsrel,limit)
        else:
            return _quadpack._qagie(func,bound,inf,args,full_output,epsabs,epsrel,limit)
    else:
        if inf !=0:
            raise ValueError, "Infinity inputs cannot be used with break points."
        else:
            nl = len(myasarray(points))
            the_points = zeros((nl+2,),'d')
            the_points[:nl] = points
            return _quadpack._qagpe(func,a,b,the_points,args,full_output,epsabs,epsrel,limit)


def _quad_weight(func,a,b,args,full_output,epsabs,epsrel,limlst,limit,maxp1,weight,wvar,wopts):

    if weight not in ['cos','sin','alg','alg-loga','alg-logb','alg-log','cauchy']:
        raise ValueError, "%s not a recognized weighting function." % weight

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
            raise ValueError, "Cannot integrate with this weight from -Inf to +Inf."
    else:
        if a in [-Inf,Inf] or b in [-Inf,Inf]:
            raise ValueError, "Cannot integrate with this weight over an infinite interval."

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
    """Compute a double (definite) integral.

  Description:

    Return the double integral of func2d(y,x) from x=a..b and y=gfun(x)..hfun(x).

  Inputs:

    func2d -- a Python function or method of at least two variables: y must be
              the first argument and x the second argument.
    (a,b) -- the limits of integration in x: a < b
    gfun -- the lower boundary curve in y which is a function taking a single
              floating point argument (x) and returning a floating point result:
              a lambda function can be useful here.  
    hfun -- the upper boundary curve in y (same requirements as gfun).
    args -- extra arguments to pass to func2d.
    epsabs -- absolute tolerance passed directly to the inner 1-D quadrature
              integration.
    epsrel -- relative tolerance of the inner 1-D integrals.

  Outputs: (y, abserr)

    y -- the resultant integral.
    abserr -- an estimate of the error.

    """
    return quad(_infunc,a,b,(func,gfun,hfun,args),epsabs=epsabs,epsrel=epsrel)

def _infunc2(y,x,func,qfun,rfun,more_args):
    a2 = qfun(x,y)
    b2 = rfun(x,y)
    myargs = (y,x) + more_args
    return quad(func,a2,b2,args=myargs)[0]
             
def tplquad(func, a, b, gfun, hfun, qfun, rfun, args=(), epsabs=1.49e-8,
            epsrel=1.49e-8):
    """Compute a triple (definite) integral.
    
  Description:

    Return the triple integral of func3d(z, y,x) from x=a..b, y=gfun(x)..hfun(x),
    and z=qfun(x,y)..rfun(x,y)

  Inputs:

    func3d -- a Python function or method of at least three variables in the
              order (z, y, x).
    (a,b) -- the limits of integration in x: a < b
    gfun -- the lower boundary curve in y which is a function taking a single
              floating point argument (x) and returning a floating point result:
              a lambda function can be useful here.  
    hfun -- the upper boundary curve in y (same requirements as gfun).
    qfun -- the lower boundary surface in z.  It must be a function that takes
            two floats in the order (x, y) and returns a float.
    rfun -- the upper boundary surface in z. (Same requirements as qfun.)
    args -- extra arguments to pass to func3d.
    epsabs -- absolute tolerance passed directly to the innermost 1-D quadrature
              integration.
    epsrel -- relative tolerance of the innermost 1-D integrals.

  Outputs: (y, abserr)

    y -- the resultant integral.
    abserr -- an estimate of the error.

    """
    return dblquad(_infunc2,a,b,gfun,hfun,(func,qfun,rfun,args),epsabs=epsabs,epsrel=epsrel)

if __name__ == '__main__':
    # Some test cases:  Note that the general purpose integrator performs
    # rather well, you don't often have to resort to the special cases.
    # I've done it here only for testing.

    test_com, comp_res, tabl_res = [],[],[]

    # 1) Typical function with two extra arguments:
    def myfunc(x,n,z):       # Bessel function integrand
        return cos(n*x-z*sin(x))/pi

    test_com.append('Typical function with two extra arguments:')
    comp_res.append(quad(myfunc,0,pi,(2,1.8)))
    tabl_res.append(0.306143535325)

    # 2) Infinite integration limits --- Euler's constant
    def myfunc(x):           # Euler's constant integrand
        return -exp(-x)*log(x)

    test_com.append('Infinite integration limits:')
    comp_res.append(quad(myfunc,0,Inf))
    tabl_res.append(0.577215664901532860606512)

    # 3) Singular points in region of integration.
    def myfunc(x):
        if x > 0 and x < 2.5:
            return sin(x)
        elif x>= 2.5 and x <= 5.0:
            return exp(-x)
        else:
            return 0.0
        
    test_com.append('Singular points in region of integration:')
    comp_res.append(quad(myfunc,0,10,points=[2.5,5.0]))
    tabl_res.append(1 - cos(2.5) + exp(-2.5) - exp(-5.0)  )

    # 4) Sine weighted integral (finite limits)
    def myfunc(x,a):
        return exp(a*(x-1))

    ome = 2.0**3.4
    comp_res.append(quad(myfunc,0,1,args=20,weight='sin',wvar=ome))
    tabl_res.append((20*sin(ome)-ome*cos(ome)+ome*exp(-20))/(20**2 + ome**2))
    test_com.append('Sine weighted integral (finite limits)')
    
    # 5) Sine weighted integral (infinite limits)
    def myfunc(x,a):           
        return exp(-x*a)

    a = 4.0
    ome = 3.0
    comp_res.append(quad(myfunc,0,Inf,args=a,weight='sin',wvar=ome))
    tabl_res.append(ome/(a**2 + ome**2))
    test_com.append('Sine weighted integral (infinite limits):')

    # 6) Cosine weighted integral (negative infinite limits)
    def myfunc(x,a):
        return exp(x*a)

    a = 2.5
    ome = 2.3
    comp_res.append(quad(myfunc,-Inf,0,args=a,weight='cos',wvar=ome))
    tabl_res.append(a/(a**2 + ome**2))
    test_com.append('Sine weighted integral (negative infinite limits):')
    
    # 6) Algebraic-logarithmic weight.
    def myfunc(x,a):
        return 1/(1+x+2**(-a))

    a = 1.5
    comp_res.append(quad(myfunc,-1,1,args=a,weight='alg',wvar=(-0.5,-0.5)))
    tabl_res.append(pi/sqrt((1+2**(-a))**2 - 1))
    test_com.append('Algebraic-logarithmic weight:')

    # 7) Cauchy prinicpal value weighting w(x) = 1/(x-c)
    def myfunc(x,a):
        return 2.0**(-a)/((x-1)**2+4.0**(-a))

    comp_res.append(quad(myfunc,0,5,args=0.4,weight='cauchy',wvar=2.0) )
    a= 0.4
    tabl_res.append((2.0**(-0.4)*log(1.5)-2.0**(-1.4)*log((4.0**(-a)+16)/(4.0**(-a)+1)) - arctan(2.0**(a+2)) - arctan(2.0**a))/(4.0**(-a) + 1))
    test_com.append('Cauchy prinicpal value with weighting w(x) = 1/(x-c):')

    # 8) Double Integral test
    def simpfunc(y,x):       # Note order of arguments.
        return x+y

    a, b = 1.0, 2.0
    comp_res.append(dblquad(simpfunc,a,b,lambda x: x, lambda x: 2*x))
    tabl_res.append(5/6.0 * (b**3.0-a**3.0))
    test_com.append('Double integral:')

    # 9) Triple Integral test
    def simpfunc(z,y,x):      # Note order of arguments.
        return x+y+z

    a, b = 1.0, 2.0
    comp_res.append(tplquad(simpfunc,a,b,lambda x: x, lambda x: 2*x, lambda x,y: x-y, lambda x,y: x+y))
    tabl_res.append(8/3.0 * (b**4.0 - a**4.0))
    test_com.append('Triple integral:')

    for k in range(len(comp_res)):
        print test_com[k]
        print "   Computed:", comp_res[k][0], "\t\t Actual:", tabl_res[k]
        print "   Error Estimate:", comp_res[k][1], "\t Actual Error:", abs(tabl_res[k]-comp_res[k][0])
        print




        





