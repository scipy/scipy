from orthogonal import P_roots
from Numeric import sum

def gauss_quad(func,a,b,n=5):
    """val = gauss_quad(func,a,b,n=5)

    Integrate func(x) from a to b using Gaussian Quadrature of order n.
    """
    [x,w] = P_roots(n)
    y = (b-a)*(x+1)/2.0 + a
    return (b-a)/2.0*sum(w*func(y))

def gauss_quadtol(func,a,b,tol=1e-7,NMAX=50):
    """val = quadtol(func,a,b,tol=1e-7,NMAX=30)

    Integrate func(x) from a to b using Gaussian Quadrature 
    with absolute tolerance tol.
    """
    err = 100.0
    val = err
    n = 1
    while (err > tol) and (n < NMAX):
        newval = quad(func,a,b,n)
        err = abs(newval-val)
        val = newval
        n = n + 1
    if (n==NMAX):
        print "NMAX (%d) exceeded. Latest difference = %e" % (n,err)
    else:
        print "Took %d points." % n
    return val
