import numpy as N
import numpy.linalg as L
import scipy.integrate
import _bspline

# note to self: check out eig_banded! in linalg.decomp?

def _zerofunc(x):
    return N.zeros(x.shape, N.float)

class BSpline:

    """
    knots should be sorted, knots[0] is lower boundary, knots[1] is upper boundary
    knots[1:-1] are internal knots
    """

    def __init__(self, knots, order=4, coef=None, M=None, eps=0.0):
        self.knots = N.squeeze(N.asarray(knots))

        if self.knots.ndim != 1:
            raise ValueError, 'expecting 1d array for knots'

        self.m = order 
        if M is None:
            M = self.m
        self.M = M
        if self.M < self.m:
            raise 'multiplicity of knots, M, must be at least equal to order, m'

        self.tau = N.hstack([[knots[0]-eps]*(self.M-1), knots, [knots[-1]+eps]*(self.M-1)])
        self.K = self.knots.shape[0] - 2
        if coef is None:
            coef = N.zeros((self.K + 2 * self.M - self.m), N.float64)
        else:
            self.coef = N.squeeze(coef)
            if self.coef.shape != (self.K + 2 * self.M - self.m):
                raise ValueError, 'coefficients of Bspline have incorrect shape'
    def __call__(self, x):
        v = 0
        b = self.basis(x)
        for i in range(self.coef.shape[0]):
            v += b[i] * self.coef[i]
        return v
    
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
        
        v = _bspline.evaluate(x, self.tau, self.m, d, lower, upper)
        v.shape = (upper-lower,) + _shape
        if upper == self.tau.shape[0] - self.m:
            v[-1] = N.where(N.equal(x, self.tau[-1]), 1, v[-1])
        return v

##         x = N.asarray(x)
##         if upper == None:
##             upper = self.tau.shape[0]-1
##         if lower == None:
##             lower = 0
##         lower = max(0, lower); upper = min(self.tau.shape[0]-1, upper)
##         which = [lower, upper]
##         which.sort()
##         lower, upper = which

##         if m is None:
##             m = self.m

##         if m == 1:
##             nbasis = upper - lower
##             v = N.zeros((nbasis,) + x.shape, N.float64)
##             for i in range(nbasis):
##                 if self.tau[i+lower] == self.tau[i+lower+1]:
##                     v[i] = N.zeros(x.shape, N.float64)
##                 else:
##                     if d <= 0:
##                         v[i] = (N.greater_equal(x, self.tau[i+lower]) *
##                                 N.less(x, self.tau[i+lower+1]))
##             return v
##         else:
##             b = self.basis(x, d=d-1, m=m-1, lower=lower,
##                            upper=upper+1)
##             nbasis = b.shape[0] - 1

##             v = N.zeros((nbasis,) + x.shape, N.float64)

##             for i in range(nbasis):

##                 if self.tau[i+lower+m-1] != self.tau[i+lower]:
##                     if d <= 0:
##                         f1 = (x - self.tau[i+lower]) / (self.tau[i+lower+m-1] - self.tau[i+lower])
##                     else:
##                         f1 = (m-1) / (self.tau[i+lower+m-1] - self.tau[i+lower])
##                 else:
##                     f1 = 0

##                 if self.tau[i+lower+m] != self.tau[i+lower+1]:
##                     if d <= 0:
##                         f2 = (self.tau[i+lower+m] - x) / (self.tau[i+lower+m] - self.tau[i+lower+1])
##                     else:
##                         f2 = -(m-1) / (self.tau[i+lower+m] - self.tau[i+lower+1])
##                 else:
##                     f2 = 0
                    
##                 v[i] = f1*b[i] + f2*b[i+1]

    def gram(self, dl=0, dr=0, full=False):
        """
        Approximate Gram inner product matrix using n
        equally spaced sample points.

        """
        
        self.g = N.nan_to_num(N.transpose(_bspline.gram(self.tau, self.m, dl, dr)))
        return self.g

class SmoothingSpline(BSpline):

    def fit(self, y, x=None, weights=None, pen=0., compute_gram=True):
        banded = True
        if x is None:
            x = self.knots # internal knots
            
        if pen == 0.: # can't do pinv for singular matrices
            banded = False
            
        if x.shape != y.shape:
            raise ValueError, 'x and y shape do not agree, by default x are the Bspline\'s internal knots'
        
        bt = self.basis(x)

        if weights is not None:
            bt *= N.sqrt(weights)

        # throw out rows with zeros (this happens at boundary points!)

        mask = N.flatnonzero(1 - N.alltrue(N.equal(bt, 0), axis=0))

        bt = bt[:,mask]
        y = y[mask]

        if compute_gram:
            self.g = self.gram(dr=2, dl=2)

        if not banded:
            btb = N.dot(bt, N.transpose(bt))
        else:
            btb = N.zeros(self.g.shape, N.float64)
            nband, nbasis = self.g.shape
            for i in range(nbasis):
                for k in range(nband):
                    j = i + self.m - 1 - k
                    if j >= 0 and j < nbasis:
                        btb[k,i] = (bt[i] * bt[j]).sum()
                        btb[nband-1-k,i] = btb[k,i]
        
        bty = N.dot(bt, y)
        if not banded:
            self.coef, r, self.rank = L.lstsq(btb + pen*self.g, bty)[0:3]
        else:
            self.coef = scipy.linalg.solve_banded((self.m-1,)*2,
                                                  btb + pen*self.g,
                                                  bty)


## s = BSpline(N.linspace(0,1,11), order=4)
## x = N.linspace(0,1,5001)
## y = s.basis(x, d=2)
## import pylab
## print s.tau
## def f(x):
##     return s.basis_element(x,6) * s.basis_element(x, 8)
## print scipy.integrate.romberg(f, 0, 1), 'integral'

## pylab.plot(x, y[6]*y[8])
## pylab.show()


import pylab
import time, gc
toc = time.time()
for i in range(1000):
    s = SmoothingSpline(N.linspace(0,1,51), order=4, M=4)
    f = s.gram(dr=2, dl=2)

gc.collect()
tic = time.time()
print (tic-toc) / 1000

## reader = csv.reader(file('/home/jtaylo/Desktop/bspline.csv'))
## v = []
## for row in reader:
##     v.append([float(x) for x in row])
## v = N.array(v)
            
import numpy.random as R
## import pylab

y = N.arange(51) + 10 * R.standard_normal((51,))
x = N.linspace(0,1,51) #s.knots[1:-1]
toc = time.time()
#s.fit(y, x=x, pen=1000.0, compute_gram=True)

y[-1] = 150.
for i in range(10):
    s.fit(y, x=x, pen=1.0e-20, compute_gram=False)
tic = time.time()
print (tic-toc) / 10

pylab.plot(x, y, 'bo')
x = N.linspace(-0.1,1.1,501)
pylab.plot(x, s(x), 'r')
pylab.show()

## print N.allclose(y, s(x)), N.add.reduce((y-s(x))**2)
## m = 0
## x = N.linspace(0,1,101)
## b = s.basis(x)
## X = []; Y=[]

## def itergrad(x, delta=1, n=1):
##     if n > 1:
##         z = N.gradient(x, delta)
##         return itergrad(z, delta=delta,n=n-1)
##     else:
##         return N.gradient(x, delta)
    
## d = 3
## for i in [0]:
##     #x = N.linspace(s.tau[i] + 0.01, s.tau[i+3] - 0.01, 1253)
##     x = N.linspace(0, 1, 1001)
##     db = s.basis(x,d=d)
##     b = s.basis(x)
##     z = itergrad(b[i], n=d, delta=x[1]-x[0])
##     X += [db[i]]
##     Y += [z]
##     pylab.plot(x, db[i])
##     #pylab.plot(x, z, 'ro')
## #    pylab.show()

## X = N.hstack(X); Y=N.hstack(Y)

## g = s.gram(dleft=2,dright=2)
## ## x = N.linspace(0,1,1000)
## ## ss = s.basis(x)
## ## G = N.zeros((g.shape[0],)*2, N.float64)
## ## for i in range(g.shape[0]):
## ##     print G.shape
## ##     for j in range(g.shape[0]):
## ##         G[i,j] = scipy.trapz(ss[i]*ss[j],x=x)
        
## ## print g.shape
## ## print N.corrcoef(X,Y)

