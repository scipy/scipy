import _zeros

_iter = 100
_xtol = 1e-12
_args = ()

__all__ = ['bisect','ridder','brentq','brenth']
    
def bisect(f,a,b,xtol=_xtol,maxiter=_iter,args=_args, full_output=0, disp=1) :
    """ Find root of f in [a,b]

        Basic bisection routine to find a zero of the function
        f between the arguments a and b. f(a) and f(b) can not
        have the same signs. Slow but sure.

        f : Python function returning a number.

        a : Number, one end of the bracketing interval.

        b : Number, the other end of the bracketing interval.

        xtol : Number, the routine converges when a root is known
            to lie within xtol of the value return. Should be >= 0.
            The routine modifies this to take into account the relative
            precision of doubles.
            
        maxiter : Number, if convergence is not achieved in
            maxiter iterations, and error is raised. Must be
            >= 0.

        args : tuple containing extra arguments for the function f.
            f is called by apply(f,(x)+args).
    """
    if type(args) != type(()) :
	    args = (args,)
    return _zeros._bisect(f,a,b,xtol,maxiter,args,full_output,disp)

def ridder(f,a,b,xtol=_xtol,maxiter=_iter,args=_args,full_output=0,disp=1) :
    """ Find root of f in [a,b]

        Ridder routine to find a zero of the function
        f between the arguments a and b. f(a) and f(b) can not
        have the same signs. Faster than bisection, but not
        generaly as fast as the brent rountines. A description
        may be found in a recent edition of Numerical Recipes.
        The routine here is modified a bit to be more careful
        of tolerance.

        f : Python function returning a number.

        a : Number, one end of the bracketing interval.

        b : Number, the other end of the bracketing interval.

        xtol : Number, the routine converges when a root is known
            to lie within xtol of the value return. Should be >= 0.
            The routine modifies this to take into account the relative
            precision of doubles.
            
        maxiter : Number, if convergence is not achieved in
            maxiter iterations, and error is raised. Must be
            >= 0.

        args : tuple containing extra arguments for the function f.
            f is called by apply(f,(x)+args).
    """
    if type(args) != type(()) :
	    args = (args,)
    return _zeros._ridder(f,a,b,xtol,maxiter,args,full_output,disp)

def brentq(f,a,b,xtol=_xtol,maxiter=_iter,args=_args,full_output=0,disp=1) :
    """ Find root of f in [a,b]

        The classic Brent routine to find a zero of the function
        f between the arguments a and b. f(a) and f(b) can not
        have the same signs. Generally the best of the routines here.
        It is a safe version of the secant method that uses inverse
        quadratic extrapolation. The version here is a slight
        modification that uses a different formula in the extrapolation
        step. A description may be found in Numerical Recipes, but the
        code here is probably easier to understand.

        f : Python function returning a number.

        a : Number, one end of the bracketing interval.

        b : Number, the other end of the bracketing interval.

        xtol : Number, the routine converges when a root is known
            to lie within xtol of the value return. Should be >= 0.
            The routine modifies this to take into account the relative
            precision of doubles.
            
        maxiter : Number, if convergence is not achieved in
            maxiter iterations, and error is raised. Must be
            >= 0.

        args : tuple containing extra arguments for the function f.
            f is called by apply(f,(x)+args).
    """
    if type(args) != type(()) :
	    args = (args,)
    return _zeros._brentq(f,a,b,xtol,maxiter,args,full_output,disp)

def brenth(f,a,b,xtol=_xtol,maxiter=_iter,args=_args,full_output=0,disp=1) :
    """ Find root of f in [a,b]

        A variation on the classic Brent routine to find a zero
        of the function f between the arguments a and b that uses
        hyperbolic extrapolation instead of inverse quadratic
        extrapolation. There was a paper back in the 1980's ...
        f(a) and f(b) can not have the same signs. Generally on a
        par with the brent routine, but not as heavily tested.
        It is a safe version of the secant method that uses hyperbolic
        extrapolation. The version here is by Chuck Harris.

        f : Python function returning a number.

        a : Number, one end of the bracketing interval.

        b : Number, the other end of the bracketing interval.

        xtol : Number, the routine converges when a root is known
            to lie within xtol of the value return. Should be >= 0.
            The routine modifies this to take into account the relative
            precision of doubles.
            
        maxiter : Number, if convergence is not achieved in
            maxiter iterations, and error is raised. Must be
            >= 0.

        args : tuple containing extra arguments for the function f.
            f is called by apply(f,(x)+args).
    """
    if type(args) != type(()) :
	    args = (args,)
    return _zeros._brenth(f,a,b,xtol,maxiter,args,full_output,disp)

