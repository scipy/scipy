import numpy as np
from scipy.misc import factorial

__all__ = ["KroghInterpolator", "krogh_interpolate", "BarycentricInterpolator", "barycentric_interpolate", "PiecewisePolynomial", "piecewise_polynomial_interpolate","approximate_taylor_polynomial", "pchip"]

class KroghInterpolator(object):
    """
    The interpolating polynomial for a set of points

    Constructs a polynomial that passes through a given set of points,
    optionally with specified derivatives at those points.
    Allows evaluation of the polynomial and all its derivatives.
    For reasons of numerical stability, this function does not compute
    the coefficients of the polynomial, although they can be obtained
    by evaluating all the derivatives.

    Be aware that the algorithms implemented here are not necessarily
    the most numerically stable known. Moreover, even in a world of
    exact computation, unless the x coordinates are chosen very
    carefully - Chebyshev zeros (e.g. cos(i*pi/n)) are a good choice -
    polynomial interpolation itself is a very ill-conditioned process
    due to the Runge phenomenon. In general, even with well-chosen
    x values, degrees higher than about thirty cause problems with
    numerical instability in this code.

    Based on [1]_.

    Parameters
    ----------
    xi : array-like, length N
        Known x-coordinates
    yi : array-like, N by R
        Known y-coordinates, interpreted as vectors of length R,
        or scalars if R=1. When an xi occurs two or more times in
        a row, the corresponding yi's represent derivative values.

    References
    ----------
    .. [1] Krogh, "Efficient Algorithms for Polynomial Interpolation
        and Numerical Differentiation", 1970.

    """
    def __init__(self, xi, yi):
        """Construct an interpolator passing through the specified points

        The polynomial passes through all the pairs (xi,yi). One may additionally
        specify a number of derivatives at each point xi; this is done by
        repeating the value xi and specifying the derivatives as successive
        yi values.

        Parameters
        ----------
        xi : array-like, length N
            known x-coordinates
        yi : array-like, N by R
            known y-coordinates, interpreted as vectors of length R,
            or scalars if R=1. When an xi occurs two or more times in
            a row, the corresponding yi's represent derivative values.

        Examples
        --------
        To produce a polynomial that is zero at 0 and 1 and has
        derivative 2 at 0, call

        >>> KroghInterpolator([0,0,1],[0,2,0])

        This constructs the quadratic 2*X**2-2*X. The derivative condition
        is indicated by the repeated zero in the xi array; the corresponding
        yi values are 0, the function value, and 2, the derivative value.

        For another example, given xi, yi, and a derivative ypi for each 
        point, appropriate arrays can be constructed as:

        >>> xi_k, yi_k = np.repeat(xi, 2), np.ravel(np.dstack((yi,ypi)))
        >>> KroghInterpolator(xi_k, yi_k)

        To produce a vector-valued polynomial, supply a higher-dimensional
        array for yi:

        >>> KroghInterpolator([0,1],[[2,3],[4,5]])

        This constructs a linear polynomial giving (2,3) at 0 and (4,5) at 1.

        """
        self.xi = np.asarray(xi)
        self.yi = np.asarray(yi)
        if len(self.yi.shape)==1:
            self.vector_valued = False
            self.yi = self.yi[:,np.newaxis]
        elif len(self.yi.shape)>2:
            raise ValueError("y coordinates must be either scalars or vectors")
        else:
            self.vector_valued = True

        n = len(xi)
        self.n = n
        nn, r = self.yi.shape
        if nn!=n:
            raise ValueError("%d x values provided and %d y values; must be equal" % (n, nn))
        self.r = r

        c = np.zeros((n+1,r))
        c[0] = yi[0]
        Vk = np.zeros((n,r))
        for k in xrange(1,n):
            s = 0
            while s<=k and xi[k-s]==xi[k]:
                s += 1
            s -= 1
            Vk[0] = yi[k]/float(factorial(s))
            for i in xrange(k-s):
                assert xi[i]!=xi[k]
                if s==0:
                    Vk[i+1] = (c[i]-Vk[i])/(xi[i]-xi[k])
                else:
                    Vk[i+1] = (Vk[i+1]-Vk[i])/(xi[i]-xi[k])
            c[k] = Vk[k-s]
        self.c = c

    def __call__(self,x):
        """Evaluate the polynomial at the point x

        Parameters
        ----------
        x : scalar or array-like of length N

        Returns
        -------
        y : scalar, array of length R, array of length N, or array of length N by R
            If x is a scalar, returns either a vector or a scalar depending on
            whether the interpolator is vector-valued or scalar-valued.
            If x is a vector, returns a vector of values.
        """
        if _isscalar(x):
            scalar = True
            m = 1
        else:
            scalar = False
            m = len(x)
        x = np.asarray(x)

        n = self.n
        pi = 1
        p = np.zeros((m,self.r))
        p += self.c[0,np.newaxis,:]
        for k in xrange(1,n):
            w = x - self.xi[k-1]
            pi = w*pi
            p = p + np.multiply.outer(pi,self.c[k])
        if not self.vector_valued:
            if scalar:
                return p[0,0]
            else:
                return p[:,0]
        else:
            if scalar:
                return p[0]
            else:
                return p

    def derivatives(self,x,der=None):
        """Evaluate many derivatives of the polynomial at the point x

        Produce an array of all derivative values at the point x.

        Parameters
        ----------
        x : scalar or array-like of length N
            Point or points at which to evaluate the derivatives
        der : None or integer
            How many derivatives to extract; None for all potentially
            nonzero derivatives (that is a number equal to the number
            of points). This number includes the function value as 0th
            derivative.
        Returns
        -------
        d : array
            If the interpolator's values are R-dimensional then the
            returned array will be der by N by R. If x is a scalar,
            the middle dimension will be dropped; if R is 1 then the
            last dimension will be dropped.

        Example
        -------
        >>> KroghInterpolator([0,0,0],[1,2,3]).derivatives(0)
        array([1.0,2.0,3.0])
        >>> KroghInterpolator([0,0,0],[1,2,3]).derivatives([0,0])
        array([[1.0,1.0],
               [2.0,2.0],
               [3.0,3.0]])
        """
        if _isscalar(x):
            scalar = True
            m = 1
        else:
            scalar = False
            m = len(x)
        x = np.asarray(x)

        n = self.n
        r = self.r

        if der is None:
            der = self.n
        dern = min(self.n,der)
        pi = np.zeros((n,m))
        w = np.zeros((n,m))
        pi[0] = 1
        p = np.zeros((m,self.r))
        p += self.c[0,np.newaxis,:]

        for k in xrange(1,n):
            w[k-1] = x - self.xi[k-1]
            pi[k] = w[k-1]*pi[k-1]
            p += np.multiply.outer(pi[k],self.c[k])

        cn = np.zeros((max(der,n+1),m,r))
        cn[:n+1,...] += self.c[:n+1,np.newaxis,:]
        cn[0] = p
        for k in xrange(1,n):
            for i in xrange(1,n-k+1):
                pi[i] = w[k+i-1]*pi[i-1]+pi[i]
                cn[k] = cn[k]+pi[i,:,np.newaxis]*cn[k+i]
            cn[k]*=factorial(k)

        cn[n,...] = 0
        if not self.vector_valued:
            if scalar:
                return cn[:der,0,0]
            else:
                return cn[:der,:,0]
        else:
            if scalar:
                return cn[:der,0]
            else:
                return cn[:der]
    def derivative(self,x,der):
        """Evaluate one derivative of the polynomial at the point x

        Parameters
        ----------
        x : scalar or array-like of length N
            Point or points at which to evaluate the derivatives
        der : None or integer
            Which derivative to extract. This number includes the
            function value as 0th derivative.
        Returns
        -------
        d : array
            If the interpolator's values are R-dimensional then the
            returned array will be N by R. If x is a scalar,
            the middle dimension will be dropped; if R is 1 then the
            last dimension will be dropped.

        Notes
        -----
        This is computed by evaluating all derivatives up to the desired
        one (using self.derivatives()) and then discarding the rest.
        """
        return self.derivatives(x,der=der+1)[der]

def krogh_interpolate(xi,yi,x,der=0):
    """Convenience function for polynomial interpolation.

    Constructs a polynomial that passes through a given set of points,
    optionally with specified derivatives at those points.
    Evaluates the polynomial or some of its derivatives.
    For reasons of numerical stability, this function does not compute
    the coefficients of the polynomial, although they can be obtained
    by evaluating all the derivatives.

    Be aware that the algorithms implemented here are not necessarily
    the most numerically stable known. Moreover, even in a world of
    exact computation, unless the x coordinates are chosen very
    carefully - Chebyshev zeros (e.g. cos(i*pi/n)) are a good choice -
    polynomial interpolation itself is a very ill-conditioned process
    due to the Runge phenomenon. In general, even with well-chosen
    x values, degrees higher than about thirty cause problems with
    numerical instability in this code.

    Based on Krogh 1970, "Efficient Algorithms for Polynomial Interpolation
    and Numerical Differentiation"

    The polynomial passes through all the pairs (xi,yi). One may additionally
    specify a number of derivatives at each point xi; this is done by
    repeating the value xi and specifying the derivatives as successive
    yi values.

    Parameters
    ----------
    xi : array-like, length N
        known x-coordinates
    yi : array-like, N by R
        known y-coordinates, interpreted as vectors of length R,
        or scalars if R=1
    x : scalar or array-like of length N
        Point or points at which to evaluate the derivatives
    der : integer or list
        How many derivatives to extract; None for all potentially
        nonzero derivatives (that is a number equal to the number
        of points), or a list of derivatives to extract. This number
        includes the function value as 0th derivative.
    Returns
    -------
    d : array
        If the interpolator's values are R-dimensional then the
        returned array will be the number of derivatives by N by R.
        If x is a scalar, the middle dimension will be dropped; if
        the yi are scalars then the last dimension will be dropped.

    Notes
    -----
    Construction of the interpolating polynomial is a relatively expensive
    process. If you want to evaluate it repeatedly consider using the class
    KroghInterpolator (which is what this function uses).
    """
    P = KroghInterpolator(xi, yi)
    if der==0:
        return P(x)
    elif _isscalar(der):
        return P.derivative(x,der=der)
    else:
        return P.derivatives(x,der=np.amax(der)+1)[der]




def approximate_taylor_polynomial(f,x,degree,scale,order=None):
    """
    Estimate the Taylor polynomial of f at x by polynomial fitting.

    Parameters
    ----------
    f : callable
        The function whose Taylor polynomial is sought. Should accept
        a vector of x values.
    x : scalar
        The point at which the polynomial is to be evaluated.
    degree : int
        The degree of the Taylor polynomial
    scale : scalar
        The width of the interval to use to evaluate the Taylor polynomial.
        Function values spread over a range this wide are used to fit the
        polynomial. Must be chosen carefully.
    order : int or None
        The order of the polynomial to be used in the fitting; f will be
        evaluated ``order+1`` times. If None, use `degree`.

    Returns
    -------
    p : poly1d instance
        The Taylor polynomial (translated to the origin, so that
        for example p(0)=f(x)).

    Notes
    -----
    The appropriate choice of "scale" is a trade-off; too large and the
    function differs from its Taylor polynomial too much to get a good
    answer, too small and round-off errors overwhelm the higher-order terms.
    The algorithm used becomes numerically unstable around order 30 even
    under ideal circumstances.

    Choosing order somewhat larger than degree may improve the higher-order
    terms.

    """
    if order is None:
        order=degree

    n = order+1
    # Choose n points that cluster near the endpoints of the interval in
    # a way that avoids the Runge phenomenon. Ensure, by including the
    # endpoint or not as appropriate, that one point always falls at x
    # exactly.
    xs = scale*np.cos(np.linspace(0,np.pi,n,endpoint=n%1)) + x

    P = KroghInterpolator(xs, f(xs))
    d = P.derivatives(x,der=degree+1)

    return np.poly1d((d/factorial(np.arange(degree+1)))[::-1])



class BarycentricInterpolator(object):
    """The interpolating polynomial for a set of points

    Constructs a polynomial that passes through a given set of points.
    Allows evaluation of the polynomial, efficient changing of the y
    values to be interpolated, and updating by adding more x values.
    For reasons of numerical stability, this function does not compute
    the coefficients of the polynomial.

    This class uses a "barycentric interpolation" method that treats
    the problem as a special case of rational function interpolation.
    This algorithm is quite stable, numerically, but even in a world of
    exact computation, unless the x coordinates are chosen very
    carefully - Chebyshev zeros (e.g. cos(i*pi/n)) are a good choice -
    polynomial interpolation itself is a very ill-conditioned process
    due to the Runge phenomenon.

    Based on Berrut and Trefethen 2004, "Barycentric Lagrange Interpolation".
    """
    def __init__(self, xi, yi=None):
        """Construct an object capable of interpolating functions sampled at xi

        The values yi need to be provided before the function is evaluated,
        but none of the preprocessing depends on them, so rapid updates
        are possible.

        Parameters
        ----------
        xi : array-like of length N
            The x coordinates of the points the polynomial should pass through
        yi : array-like N by R or None
            The y coordinates of the points the polynomial should pass through;
            if R>1 the polynomial is vector-valued. If None the y values
            will be supplied later.
        """
        self.n = len(xi)
        self.xi = np.asarray(xi)
        if yi is not None and len(yi)!=len(self.xi):
            raise ValueError("yi dimensions do not match xi dimensions")
        self.set_yi(yi)
        self.wi = np.zeros(self.n)
        self.wi[0] = 1
        for j in xrange(1,self.n):
            self.wi[:j]*=(self.xi[j]-self.xi[:j])
            self.wi[j] = np.multiply.reduce(self.xi[:j]-self.xi[j])
        self.wi**=-1

    def set_yi(self, yi):
        """Update the y values to be interpolated

        The barycentric interpolation algorithm requires the calculation
        of weights, but these depend only on the xi. The yi can be changed
        at any time.

        Parameters
        ----------
        yi : array-like N by R
            The y coordinates of the points the polynomial should pass through;
            if R>1 the polynomial is vector-valued. If None the y values
            will be supplied later.
        """
        if yi is None:
            self.yi = None
            return
        yi = np.asarray(yi)
        if len(yi.shape)==1:
            self.vector_valued = False
            yi = yi[:,np.newaxis]
        elif len(yi.shape)>2:
            raise ValueError("y coordinates must be either scalars or vectors")
        else:
            self.vector_valued = True

        n, r = yi.shape
        if n!=len(self.xi):
            raise ValueError("yi dimensions do not match xi dimensions")
        self.yi = yi
        self.r = r


    def add_xi(self, xi, yi=None):
        """Add more x values to the set to be interpolated

        The barycentric interpolation algorithm allows easy updating by
        adding more points for the polynomial to pass through.

        Parameters
        ----------
        xi : array-like of length N1
            The x coordinates of the points the polynomial should pass through
        yi : array-like N1 by R or None
            The y coordinates of the points the polynomial should pass through;
            if R>1 the polynomial is vector-valued. If None the y values
            will be supplied later. The yi should be specified if and only if
            the interpolator has y values specified.
        """
        if yi is not None:
            if self.yi is None:
                raise ValueError("No previous yi value to update!")
            yi = np.asarray(yi)
            if len(yi.shape)==1:
                if self.vector_valued:
                    raise ValueError("Cannot extend dimension %d y vectors with scalars" % self.r)
                yi = yi[:,np.newaxis]
            elif len(yi.shape)>2:
                raise ValueError("y coordinates must be either scalars or vectors")
            else:
                n, r = yi.shape
                if r!=self.r:
                    raise ValueError("Cannot extend dimension %d y vectors with dimension %d y vectors" % (self.r, r))

            self.yi = np.vstack((self.yi,yi))
        else:
            if self.yi is not None:
                raise ValueError("No update to yi provided!")
        old_n = self.n
        self.xi = np.concatenate((self.xi,xi))
        self.n = len(self.xi)
        self.wi**=-1
        old_wi = self.wi
        self.wi = np.zeros(self.n)
        self.wi[:old_n] = old_wi
        for j in xrange(old_n,self.n):
            self.wi[:j]*=(self.xi[j]-self.xi[:j])
            self.wi[j] = np.multiply.reduce(self.xi[:j]-self.xi[j])
        self.wi**=-1

    def __call__(self, x):
        """Evaluate the interpolating polynomial at the points x

        Parameters
        ----------
        x : scalar or array-like of length M

        Returns
        -------
        y : scalar or array-like of length R or length M or M by R
            The shape of y depends on the shape of x and whether the
            interpolator is vector-valued or scalar-valued.

        Notes
        -----
        Currently the code computes an outer product between x and the
        weights, that is, it constructs an intermediate array of size
        N by M, where N is the degree of the polynomial.
        """
        scalar = _isscalar(x)
        x = np.atleast_1d(x)
        c = np.subtract.outer(x,self.xi)
        z = c==0
        c[z] = 1
        c = self.wi/c
        p = np.dot(c,self.yi)/np.sum(c,axis=-1)[:,np.newaxis]
        i, j = np.nonzero(z)
        p[i] = self.yi[j]
        if not self.vector_valued:
            if scalar:
                return p[0,0]
            else:
                return p[:,0]
        else:
            if scalar:
                return p[0]
            else:
                return p
def barycentric_interpolate(xi, yi, x):
    """Convenience function for polynomial interpolation

    Constructs a polynomial that passes through a given set of points,
    then evaluates the polynomial. For reasons of numerical stability,
    this function does not compute the coefficients of the polynomial.

    This function uses a "barycentric interpolation" method that treats
    the problem as a special case of rational function interpolation.
    This algorithm is quite stable, numerically, but even in a world of
    exact computation, unless the x coordinates are chosen very
    carefully - Chebyshev zeros (e.g. cos(i*pi/n)) are a good choice -
    polynomial interpolation itself is a very ill-conditioned process
    due to the Runge phenomenon.

    Based on Berrut and Trefethen 2004, "Barycentric Lagrange Interpolation".

    Parameters
    ----------
    xi : array-like of length N
        The x coordinates of the points the polynomial should pass through
    yi : array-like N by R
        The y coordinates of the points the polynomial should pass through;
        if R>1 the polynomial is vector-valued.
    x : scalar or array-like of length M

    Returns
    -------
    y : scalar or array-like of length R or length M or M by R
        The shape of y depends on the shape of x and whether the
        interpolator is vector-valued or scalar-valued.

    Notes
    -----

    Construction of the interpolation weights is a relatively slow process.
    If you want to call this many times with the same xi (but possibly
    varying yi or x) you should use the class BarycentricInterpolator.
    This is what this function uses internally.
    """
    return BarycentricInterpolator(xi, yi)(x)


class PiecewisePolynomial(object):
    """Piecewise polynomial curve specified by points and derivatives

    This class represents a curve that is a piecewise polynomial. It
    passes through a list of points and has specified derivatives at
    each point. The degree of the polynomial may very from segment to
    segment, as may the number of derivatives available. The degree
    should not exceed about thirty.

    Appending points to the end of the curve is efficient.
    """
    def __init__(self, xi, yi, orders=None, direction=None):
        """Construct a piecewise polynomial

        Parameters
        ----------
        xi : array-like of length N
            a sorted list of x-coordinates
        yi : list of lists of length N
            yi[i] is the list of derivatives known at xi[i]
        orders : list of integers, or integer
            a list of polynomial orders, or a single universal order
        direction : {None, 1, -1}
            indicates whether the xi are increasing or decreasing
            +1 indicates increasing
            -1 indicates decreasing
            None indicates that it should be deduced from the first two xi

        Notes
        -----
        If orders is None, or orders[i] is None, then the degree of the
        polynomial segment is exactly the degree required to match all i
        available derivatives at both endpoints. If orders[i] is not None,
        then some derivatives will be ignored. The code will try to use an
        equal number of derivatives from each end; if the total number of
        derivatives needed is odd, it will prefer the rightmost endpoint. If
        not enough derivatives are available, an exception is raised.
        """
        yi0 = np.asarray(yi[0])
        if len(yi0.shape)==2:
            self.vector_valued = True
            self.r = yi0.shape[1]
        elif len(yi0.shape)==1:
            self.vector_valued = False
            self.r = 1
        else:
            raise ValueError("Each derivative must be a vector, not a higher-rank array")

        self.xi = [xi[0]]
        self.yi = [yi0]
        self.n = 1

        self.direction = direction
        self.orders = []
        self.polynomials = []
        self.extend(xi[1:],yi[1:],orders)

    def _make_polynomial(self,x1,y1,x2,y2,order,direction):
        """Construct the interpolating polynomial object

        Deduces the number of derivatives to match at each end
        from order and the number of derivatives available. If
        possible it uses the same number of derivatives from
        each end; if the number is odd it tries to take the
        extra one from y2. In any case if not enough derivatives
        are available at one end or another it draws enough to
        make up the total from the other end.
        """
        n = order+1
        n1 = min(n//2,len(y1))
        n2 = min(n-n1,len(y2))
        n1 = min(n-n2,len(y1))
        if n1+n2!=n:
            raise ValueError("Point %g has %d derivatives, point %g has %d derivatives, but order %d requested" % (x1, len(y1), x2, len(y2), order))
        assert n1<=len(y1)
        assert n2<=len(y2)

        xi = np.zeros(n)
        if self.vector_valued:
            yi = np.zeros((n,self.r))
        else:
            yi = np.zeros((n,))

        xi[:n1] = x1
        yi[:n1] = y1[:n1]
        xi[n1:] = x2
        yi[n1:] = y2[:n2]

        return KroghInterpolator(xi,yi)

    def append(self, xi, yi, order=None):
        """Append a single point with derivatives to the PiecewisePolynomial

        Parameters
        ----------
        xi : float
        yi : array-like
            yi is the list of derivatives known at xi
        order : integer or None
            a polynomial order, or instructions to use the highest
            possible order
        """

        yi = np.asarray(yi)
        if self.vector_valued:
            if (len(yi.shape)!=2 or yi.shape[1]!=self.r):
                raise ValueError("Each derivative must be a vector of length %d" % self.r)
        else:
            if len(yi.shape)!=1:
                raise ValueError("Each derivative must be a scalar")

        if self.direction is None:
            self.direction = np.sign(xi-self.xi[-1])
        elif (xi-self.xi[-1])*self.direction < 0:
            raise ValueError("x coordinates must be in the %d direction: %s" % (self.direction, self.xi))

        self.xi.append(xi)
        self.yi.append(yi)


        if order is None:
            n1 = len(self.yi[-2])
            n2 = len(self.yi[-1])
            n = n1+n2
            order = n-1

        self.orders.append(order)
        self.polynomials.append(self._make_polynomial(
            self.xi[-2], self.yi[-2],
            self.xi[-1], self.yi[-1],
            order, self.direction))
        self.n += 1


    def extend(self, xi, yi, orders=None):
        """Extend the PiecewisePolynomial by a list of points

        Parameters
        ----------
        xi : array-like of length N1
            a sorted list of x-coordinates
        yi : list of lists of length N1
            yi[i] is the list of derivatives known at xi[i]
        orders : list of integers, or integer
            a list of polynomial orders, or a single universal order
        direction : {None, 1, -1}
            indicates whether the xi are increasing or decreasing
            +1 indicates increasing
            -1 indicates decreasing
            None indicates that it should be deduced from the first two xi
        """

        for i in xrange(len(xi)):
            if orders is None or _isscalar(orders):
                self.append(xi[i],yi[i],orders)
            else:
                self.append(xi[i],yi[i],orders[i])

    def __call__(self, x):
        """Evaluate the piecewise polynomial

        Parameters
        ----------
        x : scalar or array-like of length N

        Returns
        -------
        y : scalar or array-like of length R or length N or N by R
        """
        if _isscalar(x):
            pos = np.clip(np.searchsorted(self.xi, x) - 1, 0, self.n-2)
            y = self.polynomials[pos](x)
        else:
            x = np.asarray(x)
            m = len(x)
            pos = np.clip(np.searchsorted(self.xi, x) - 1, 0, self.n-2)
            if self.vector_valued:
                y = np.zeros((m,self.r))
            else:
                y = np.zeros(m)
            for i in xrange(self.n-1):
                c = pos==i
                y[c] = self.polynomials[i](x[c])
        return y

    def derivative(self, x, der):
        """Evaluate a derivative of the piecewise polynomial

        Parameters
        ----------
        x : scalar or array-like of length N
        der : integer
            which single derivative to extract

        Returns
        -------
        y : scalar or array-like of length R or length N or N by R

        Notes
        -----
        This currently computes (using self.derivatives()) all derivatives
        of the curve segment containing each x but returns only one.
        """
        return self.derivatives(x,der=der+1)[der]

    def derivatives(self, x, der):
        """Evaluate a derivative of the piecewise polynomial
        Parameters
        ----------
        x : scalar or array-like of length N
        der : integer
            how many derivatives (including the function value as
            0th derivative) to extract

        Returns
        -------
        y : array-like of shape der by R or der by N or der by N by R

        """
        if _isscalar(x):
            pos = np.clip(np.searchsorted(self.xi, x) - 1, 0, self.n-2)
            y = self.polynomials[pos].derivatives(x,der=der)
        else:
            x = np.asarray(x)
            m = len(x)
            pos = np.clip(np.searchsorted(self.xi, x) - 1, 0, self.n-2)
            if self.vector_valued:
                y = np.zeros((der,m,self.r))
            else:
                y = np.zeros((der,m))
            for i in xrange(self.n-1):
                c = pos==i
                y[:,c] = self.polynomials[i].derivatives(x[c],der=der)
        return y


def piecewise_polynomial_interpolate(xi,yi,x,orders=None,der=0):
    """Convenience function for piecewise polynomial interpolation

    Parameters
    ----------
    xi : array-like of length N
        a sorted list of x-coordinates
    yi : list of lists of length N
        yi[i] is the list of derivatives known at xi[i]
    x : scalar or array-like of length M
    orders : list of integers, or integer
        a list of polynomial orders, or a single universal order
    der : integer
        which single derivative to extract

    Returns
    -------
    y : scalar or array-like of length R or length M or M by R

    Notes
    -----
    If orders is None, or orders[i] is None, then the degree of the
    polynomial segment is exactly the degree required to match all i
    available derivatives at both endpoints. If orders[i] is not None,
    then some derivatives will be ignored. The code will try to use an
    equal number of derivatives from each end; if the total number of
    derivatives needed is odd, it will prefer the rightmost endpoint. If
    not enough derivatives are available, an exception is raised.

    Construction of these piecewise polynomials can be an expensive process;
    if you repeatedly evaluate the same polynomial, consider using the class
    PiecewisePolynomial (which is what this function does).
    """

    P = PiecewisePolynomial(xi, yi, orders)
    if der==0:
        return P(x)
    elif _isscalar(der):
        return P.derivative(x,der=der)
    else:
        return P.derivatives(x,der=np.amax(der)+1)[der]

def _isscalar(x):
    """Check whether x is if a scalar type, or 0-dim"""
    return np.isscalar(x) or hasattr(x, 'shape') and x.shape == ()

def _edge_case(m0, d1):
    return np.where((d1==0) | (m0==0), 0.0, 1.0/(1.0/m0+1.0/d1))

def _find_derivatives(x, y):
    # Determine the derivatives at the points y_k, d_k, by using
    #  PCHIP algorithm is:
    # We choose the derivatives at the point x_k by
    # Let m_k be the slope of the kth segment (between k and k+1)
    # If m_k=0 or m_{k-1}=0 or sgn(m_k) != sgn(m_{k-1}) then d_k == 0
    # else use weighted harmonic mean:
    #   w_1 = 2h_k + h_{k-1}, w_2 = h_k + 2h_{k-1}
    #   1/d_k = 1/(w_1 + w_2)*(w_1 / m_k + w_2 / m_{k-1})
    #   where h_k is the spacing between x_k and x_{k+1}

    hk = x[1:] - x[:-1]
    mk = (y[1:] - y[:-1]) / hk
    smk = np.sign(mk)
    condition = ((smk[1:] != smk[:-1]) | (mk[1:]==0) | (mk[:-1]==0))

    w1 = 2*hk[1:] + hk[:-1]
    w2 = hk[1:] + 2*hk[:-1]
    whmean = 1.0/(w1+w2)*(w1/mk[1:] + w2/mk[:-1])
    
    dk = np.zeros_like(y)
    dk[1:-1][condition] = 0.0
    dk[1:-1][~condition] = 1.0/whmean[~condition]

    # For end-points choose d_0 so that 1/d_0 = 1/m_0 + 1/d_1 unless
    #  one of d_1 or m_0 is 0, then choose d_0 = 0

    dk[0] = _edge_case(mk[0],dk[1])
    dk[-1] = _edge_case(mk[-1],dk[-2])
    return dk
    

def pchip(x, y):
    """PCHIP 1-d monotonic cubic interpolation 
    
    Description
    -----------
    x and y are arrays of values used to approximate some function f:
       y = f(x)
    This class factory function returns a callable class whose __call__ method 
    uses monotonic cubic, interpolation to find the value of new points.

    Parameters
    ----------
    x : array
        A 1D array of monotonically increasing real values.  x cannot
        include duplicate values (otherwise f is overspecified)
    y : array
        A 1-D array of real values.  y's length along the interpolation
        axis must be equal to the length of x.
 
    Assumes x is sorted in monotonic order (e.g. x[1] > x[0])
    """
    derivs = _find_derivatives(x,y)
    return PiecewisePolynomial(x, zip(y, derivs), orders=3, direction=None)

    
