from orthogonal import P_roots
from Numeric import sum

def gauss_quad(func,a,b,args=(),n=5):
    """Compute a definite integral using fixed-order Gaussian quadrature.

  Description:

    Integrate func from a to b using Gaussian quadrature of order n.

  Inputs:

    func -- a Python function or method to integrate.
    a -- lower limit of integration
    b -- upper limit of integration
    args -- extra arguments to pass to function.
    n -- order of quadrature integration.

  Outputs: (val,)

    val -- Gaussian quadrature approximation to the integral.
    
    """
    [x,w] = P_roots(n)
    y = (b-a)*(x+1)/2.0 + a
    return (b-a)/2.0*sum(w*func(y,*args))

def gauss_quadtol(func,a,b,args=(),tol=1e-7,maxiter=50):
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

  Outputs: (val, )

    val -- Gaussian quadrature approximation (within tolerance) to integral.

    """
    err = 100.0
    val = err
    n = 1
    while (err > tol) and (n < maxiter):
        newval = gauss_quad(func,a,b,args,n)
        err = abs(newval-val)
        val = newval
        n = n + 1
    if (n==maxiter):
        print "maxiter (%d) exceeded. Latest difference = %e" % (n,err)
    else:
        print "Took %d points." % n
    return val
