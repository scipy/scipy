# Author: Travis Oliphant


from orthogonal import P_roots
from Numeric import sum

def fixed_quad(func,a,b,args=(),n=5):
    """Compute a definite integral using fixed-order Gaussian quadrature.

  Description:

    Integrate func from a to b using Gaussian quadrature of order n.

  Inputs:

    func -- a Python function or method to integrate.
    a -- lower limit of integration
    b -- upper limit of integration
    args -- extra arguments to pass to function.
    n -- order of quadrature integration.

  Outputs: (val, None)

    val -- Gaussian quadrature approximation to the integral.
    
    """
    [x,w] = P_roots(n)
    y = (b-a)*(x+1)/2.0 + a
    return (b-a)/2.0*sum(w*func(y,*args)), None

def quadrature(func,a,b,args=(),tol=1e-7,maxiter=50):
    """Compute a definite integral using fixed-tolerance Gaussian quadrature.

  Description:

    Integrate func from a to b using Gaussian quadrature 
    with absolute tolerance tol.

  Inputs:

    func -- a Python function or method to integrate.
    a -- lower limit of integration.
    b -- upper limit of integration.
    args -- extra arguments to pass to function.
    tol -- iteration stops when error between last two iterates is less than
           tolerance.
    maxiter -- maximum number of iterations.

  Outputs: (val, err)

    val -- Gaussian quadrature approximation (within tolerance) to integral.
    err -- Difference between last two estimates of the integral.

    """
    err = 100.0
    val = err
    n = 1
    while (err > tol) and (n < maxiter):
        newval = fixed_quad(func,a,b,args,n)[0]
        err = abs(newval-val)
        val = newval
        n = n + 1
    if (n==maxiter):
        print "maxiter (%d) exceeded. Latest difference = %e" % (n,err)
    else:
        print "Took %d points." % n
    return val, err


# Romberg quadratures for numeric integration.
#
# Written by Scott M. Ransom <ransom@cfa.harvard.edu>
# last revision: 14 Nov 98
#
# Cosmetic changes by Konrad Hinsen <hinsen@cnrs-orleans.fr>
# last revision: 1999-7-21
#
# Adapted to scipy by Travis Oliphant <oliphant.travis@ieee.org>
# last revision: Dec 2001

def _difftrap(function, interval, numtraps):
    # Perform part of the trapezoidal rule to integrate a function.
    # Assume that we had called difftrap with all lower powers-of-2
    # starting with 1.  Calling difftrap only returns the summation
    # of the new ordinates.  It does _not_ multiply by the width
    # of the trapezoids.  This must be performed by the caller.
    #     'function' is the function to evaluate.
    #     'interval' is a sequence with lower and upper limits
    #                of integration.
    #     'numtraps' is the number of trapezoids to use (must be a
    #                power-of-2).
    if numtraps<=0:
        print "numtraps must be > 0 in difftrap()."
        return
    elif numtraps==1:
        return 0.5*(function(interval[0])+function(interval[1]))
    else:
        numtosum = numtraps/2
        h = float(interval[1]-interval[0])/numtosum
        lox = interval[0] + 0.5 * h;
        sum = 0.0
        for i in range(0, numtosum):
            sum = sum + function(lox + i*h)
        return sum

def _romberg_diff(b, c, k):
    # Compute the differences for the Romberg quadrature corrections.
    # See Forman Acton's "Real Computing Made Real," p 143.
    tmp = 4.0**k
    return (tmp * c - b)/(tmp - 1.0)

def _printresmat(function, interval, resmat):
    # Print the Romberg result matrix.
    i = j = 0
    print 'Romberg integration of', `function`,
    print 'from', interval
    print ''
    print '%6s %9s %9s' % ('Steps', 'StepSize', 'Results')
    for i in range(len(resmat)):
        print '%6d %9f' % (2**i, (interval[1]-interval[0])/(i+1.0)),
        for j in range(i+1):
            print '%9f' % (resmat[i][j]),
        print ''
    print ''
    print 'The final result is', resmat[i][j],
    print 'after', 2**(len(resmat)-1)+1, 'function evaluations.'

def romberg(function, a, b, tol=1.0E-7, show=0, divmax=10):
    """Romberg integration of a callable function or method.

    Returns the integral of |function| (a function of one variable)
    over |interval| (a sequence of length two containing the lower and
    upper limit of the integration interval), calculated using
    Romberg integration up to the specified |accuracy|. If |show| is 1,
    the triangular array of the intermediate results will be printed.
    """
    i = n = 1
    interval = [a,b]
    intrange = b-a
    ordsum = _difftrap(function, interval, n)
    result = intrange * ordsum
    resmat = [[result]]
    lastresult = result + tol * 2.0
    while (abs(result - lastresult) > tol) and (i <= divmax):
        n = n * 2
        ordsum = ordsum + _difftrap(function, interval, n)
        resmat.append([])
        resmat[i].append(intrange * ordsum / n)
        for k in range(i):
            resmat[i].append(_romberg_diff(resmat[i-1][k],
                                          resmat[i][k], k+1))
        result = resmat[i][i]
        lastresult = resmat[i-1][i-1]
        i = i + 1
    if show: _printresmat(function, interval, resmat)
    return result
    
    

