
import numpy as N
import numpy.linalg as L

from scipy.optimize import golden
from scipy.sandbox.models import _bspline
from scipy.linalg import solveh_banded

def _upper2lower(ub):
    """
    Convert upper triangular banded matrix to lower banded form.
    """

    lb = N.zeros(ub.shape, ub.dtype)
    nrow, ncol = ub.shape
    for i in range(ub.shape[0]):
        lb[i,0:(ncol-i)] = ub[nrow-1-i,i:ncol]
        lb[i,(ncol-i):] = ub[nrow-1-i,0:i]
    return lb

def _lower2upper(lb):
    """
    Convert upper triangular banded matrix to lower banded form.
    """

    ub = N.zeros(lb.shape, lb.dtype)
    nrow, ncol = lb.shape
    for i in range(lb.shape[0]):
        ub[nrow-1-i,i:ncol] = lb[i,0:(ncol-i)]
        ub[nrow-1-i,0:i] = lb[i,(ncol-i):]
    return ub

def _triangle2unit(tb, lower=0):
    """
    Take a banded triangular matrix and return its diagonal and the unit matrix:
    the banded triangular matrix with 1's on the diagonal.
    """
    
    if lower: d = tb[0].copy()
    else: d = tb[-1].copy()

    if lower: return d, (tb / d)
    else:
        l = _upper2lower(tb)
        return d, _lower2upper(l / d)
    
def _trace_symbanded(a,b, lower=0):
    """
    Compute the trace(a*b) for two upper or lower banded real symmetric matrices.
    """

    if lower:
	t = _zero_triband(a * b, lower=1)
	return t[0].sum() + 2 * t[1:].sum()
    else:
	t = _zero_triband(a * b, lower=0)
	return t[-1].sum() + 2 * t[:-1].sum()


def _zero_triband(a, lower=0):
    """
    Zero out unnecessary elements of a real symmetric banded matrix.
    """

    nrow, ncol = a.shape
    if lower: 
	for i in range(nrow): a[i,(ncol-i):] = 0.
    else:
	for i in range(nrow): a[i,0:i] = 0.
    return a

def _zerofunc(x):
    return N.zeros(x.shape, N.float)


class BSpline:

    """
    knots should be sorted, knots[0] is lower boundary, knots[1] is upper boundary
    knots[1:-1] are internal knots
    """

    def __init__(self, knots, order=4, coef=None, M=None, eps=0.0):
        knots = N.squeeze(N.unique(N.asarray(knots)))

        if knots.ndim != 1:
            raise ValueError, 'expecting 1d array for knots'

        self.m = order 
        if M is None:
            M = self.m
        self.M = M
#         if self.M < self.m:
#             raise 'multiplicity of knots, M, must be at least equal to order, m'

        self.tau = N.hstack([[knots[0]-eps]*(self.M-1), knots, [knots[-1]+eps]*(self.M-1)])
        self.K = knots.shape[0] - 2
        if coef is None:
            self.coef = N.zeros((self.K + 2 * self.M - self.m), N.float64)
        else:
            self.coef = N.squeeze(coef)
            if self.coef.shape != (self.K + 2 * self.M - self.m):
                raise ValueError, 'coefficients of Bspline have incorrect shape'

    def __call__(self, x):
        b = N.asarray(self.basis(x)).T
        return N.squeeze(N.dot(b, self.coef)) 
    
    def basis_element(self, x, i, d=0):
        x = N.asarray(x, N.float64)
        _shape = x.shape
        if _shape == ():
            x.shape = (1,)
        x.shape = (N.product(_shape,axis=0),)
        if i < self.tau.shape[0] - 1:
            ## TODO: OWNDATA flags...
            v = _bspline.evaluate(x, self.tau, self.m, d, i, i+1)
        else:
            return N.zeros(x.shape, N.float64)

        if (i == self.tau.shape[0] - self.m):
            v = N.where(N.equal(x, self.tau[-1]), 1, v)
        v.shape = _shape
        return v

    def basis(self, x, d=0, upper=None, lower=None):
        x = N.asarray(x)
        _shape = x.shape
        if _shape == ():
            x.shape = (1,)
        x.shape = (N.product(_shape,axis=0),)

        if upper is None:
            upper = self.tau.shape[0] - self.m 
        if lower is None:
            lower = 0
        upper = min(upper, self.tau.shape[0] - self.m)
        lower = max(0, lower)
        
	d = N.asarray(d)
	if d.shape == ():
	    v = _bspline.evaluate(x, self.tau, self.m, int(d), lower, upper)
	else:
	    if d.shape[0] != 2:
		raise ValueError, "if d is not an integer, expecting a jx2 array with first row indicating order \
		of derivative, second row coefficient in front."

	    v = 0
	    for i in range(d.shape[1]):
		v += d[1,i] * _bspline.evaluate(x, self.tau, self.m, d[0,i], lower, upper)

        v.shape = (upper-lower,) + _shape
        if upper == self.tau.shape[0] - self.m:
            v[-1] = N.where(N.equal(x, self.tau[-1]), 1, v[-1])
        return v

    def gram(self, d=0, full=False):
        """
        Compute Gram inner product matrix.
        """
        
	d = N.squeeze(d)
	if N.asarray(d).shape == ():
	    self.g = _bspline.gram(self.tau, self.m, int(d), int(d))
	else:
	    d = N.asarray(d)
	    if d.shape[0] != 2:
		raise ValueError, "if d is not an integer, expecting a jx2 array with first row indicating order \
		of derivative, second row coefficient in front."
	    if d.shape == (2,):
		d.shape = (2,1)
	    self.g = 0
	    for i in range(d.shape[1]):
                for j in range(d.shape[1]):
		    self.g += d[1,i]* d[1,j] * _bspline.gram(self.tau, self.m, int(d[0,i]), int(d[0,j]))
	self.g = self.g.T
	self.d = d
	return N.nan_to_num(self.g)

class SmoothingSpline(BSpline):

    penmax = 30.
    method = "target_df"
    target_df = 5
    default_pen = 1.0e-03

    def smooth(self, y, x=None, weights=None):
	if self.method == "target_df":
	    self.fit_target_df(y, x=x, weights=weights, df=self.target_df)
	elif self.method == "optimize_gcv":
	    self.fit_optimize_gcv(y, x=x, weights=weights)

    def fit(self, y, x=None, weights=None, pen=0.):
        banded = True

        if x is None:
            x = self.tau[(self.M-1):-(self.M-1)] # internal knots

        if pen == 0.: # can't use cholesky for singular matrices
            banded = False
            
        if x.shape != y.shape:
            raise ValueError, 'x and y shape do not agree, by default x are the Bspline\'s internal knots'
        
	bt = self.basis(x)
	if pen >= self.penmax:
	    pen = self.penmax


        if weights is not None:
	    self.weights = weights 
	else:
	    self.weights = 1. 

	_w = N.sqrt(self.weights)
	bt *= _w

        # throw out rows with zeros (this happens at boundary points!)

        mask = N.flatnonzero(1 - N.alltrue(N.equal(bt, 0), axis=0))

        bt = bt[:,mask]
        y = y[mask]

	self.df_total = y.shape[0]
        bty = N.dot(bt, _w * y)
	self.N = y.shape[0]
        if not banded:
            self.btb = N.dot(bt, bt.T)
	    _g = band2array(self.g, lower=1, symmetric=True)
            self.coef, _, self.rank = L.lstsq(self.btb + pen*_g, bty)[0:3]
	    self.rank = min(self.rank, self.btb.shape[0])
        else:
            self.btb = N.zeros(self.g.shape, N.float64)
            nband, nbasis = self.g.shape
            for i in range(nbasis):
                for k in range(min(nband, nbasis-i)):
		    self.btb[k,i] = (bt[i] * bt[i+k]).sum()

	    bty.shape = (1,bty.shape[0])
            self.chol, self.coef = solveh_banded(self.btb + 
						 pen*self.g,
						 bty, lower=1)

	self.coef = N.squeeze(self.coef)
	self.resid = y * self.weights - N.dot(self.coef, bt)	    
	self.pen = pen

    def gcv(self):
	"""
	Generalized cross-validation score of current fit.
	"""

	norm_resid = (self.resid**2).sum()
	return norm_resid / (self.df_total - self.trace())

    def df_resid(self):
	"""
	self.N - self.trace()

	where self.N is the number of observations of last fit.
	"""
	
	return self.N - self.trace()

    def df_fit(self):
	"""
	= self.trace()

	How many degrees of freedom used in the fit? 
	"""
	return self.trace()

    def trace(self):
	"""
	Trace of the smoothing matrix S(pen)
	"""

	if self.pen > 0:
	    _invband = _bspline.invband(self.chol.copy())
	    tr = _trace_symbanded(_invband, self.btb, lower=1)
	    return tr
	else:
	    return self.rank

    def fit_target_df(self, y, x=None, df=None, weights=None, tol=1.0e-03):
	"""
	Fit smoothing spline with approximately df degrees of freedom
        used in the fit, i.e. so that self.trace() is approximately df.

	In general, df must be greater than the dimension of the null space
	of the Gram inner product. For cubic smoothing splines, this means
	that df > 2.

        """

	df = df or self.target_df

	apen, bpen = 0, 1.0e-03
	olddf = y.shape[0] - self.m

	if hasattr(self, "pen"):
	    self.fit(y, x=x, weights=weights, pen=self.pen)
	    curdf = self.trace()
	    if N.fabs(curdf - df) / df < tol: 
		return
	    if curdf > df:
		apen, bpen = self.pen, 2 * self.pen
	    else:
		apen, bpen = 0., self.pen

	while True:
	    curpen = 0.5 * (apen + bpen)
	    self.fit(y, x=x, weights=weights, pen=curpen)
	    curdf = self.trace()
	    if curdf > df:
		apen, bpen = curpen, 2 * curpen
	    else:
		apen, bpen = apen, curpen
	    if apen >= self.penmax:
		raise ValueError, "penalty too large, try setting penmax higher or decreasing df"
	    if N.fabs(curdf - df) / df < tol: 
		break

    def fit_optimize_gcv(self, y, x=None, weights=None, tol=1.0e-03,
		     bracket=(0,1.0e-03)):
	"""
	Fit smoothing spline trying to optimize GCV.

	Try to find a bracketing interval for scipy.optimize.golden
	based on bracket.

	It is probably best to use target_df instead, as it is 
	sometimes difficult to find a bracketing interval.

        """
    
	def _gcv(pen, y, x):
	    self.fit(y, x=x, pen=N.exp(pen))
	    a = self.gcv()
	    return a

	a = golden(_gcv, args=(y,x), brack=(-100,20), tol=tol)


def band2array(a, lower=0, symmetric=False, hermitian=False):
    """
    Take an upper or lower triangular banded matrix and return a matrix using
    LAPACK storage convention. For testing banded Cholesky decomposition, etc.
    """

    n = a.shape[1]
    r = a.shape[0]
    _a = 0

    if not lower:
        for j in range(r):
            _b = N.diag(a[r-1-j],k=j)[j:(n+j),j:(n+j)]
            _a += _b 
            if symmetric and j > 0: _a += _b.T
            elif hermitian and j > 0: _a += _b.conjugate().T
    else:
        for j in range(r):
            _b = N.diag(a[j],k=j)[0:n,0:n]
            _a += _b 
            if symmetric and j > 0: _a += _b.T
            elif hermitian and j > 0: _a += _b.conjugate().T
	_a = _a.T

    return _a

