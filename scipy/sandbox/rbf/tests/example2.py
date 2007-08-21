from numpy import sin, asarray, exp, random, mgrid, pi, cos, sqrt, ones
from scipy.sandbox.rbf import Rbf
import pylab as pl

def truth_2d(x,y,w=2*pi):
    "moguls"
    return sin(w*x)*cos(w*y)

def truth_nd(*args):
    "a gausian sphere"
    x = asarray(list(args), 'float64')
    return exp( -sqrt((x**2).sum(axis=0)) )

# 2D example
N = 300
xi = random.rand(N)
yi = random.rand(N)
di = truth_2d(xi, yi)
xa, ya = mgrid[0:1:50j, 0:1:50j]
s = Rbf(xi, yi, di)
da = s(xa, ya)
pl.figure()
n = pl.normalize(-1., 1.)
pl.pcolor(xa, ya, da, norm=n, cmap=pl.cm.jet)
pl.scatter(xi, yi, 100, di, norm=n, cmap=pl.cm.jet)
pl.axis([0., 1., 0., 1.])
pl.colorbar()
pl.draw()
# 3d example
N = 300
xi = 2.*random.randn(N)
yi = 2.*random.randn(N)
zi = 2.*random.randn(N)
di = truth_nd(xi, yi, zi)
zas = [-0.25, 0.0, 0.25, 0.75]
xa, ya = mgrid[-1:1:50j, -1:1:50j]
s = Rbf(xi, yi, zi, di)
fig = pl.figure(figsize=(12, 3))
for idx, za in enumerate(zas):
    da = s(xa, ya, za*ones(xa.shape, 'f'))
    ax = fig.add_subplot(1,4,idx+1)
    ax.pcolor(xa, ya, da, norm=pl.normalize(0, 1), \
	  shading='flat', cmap=pl.cm.jet)
    ax.set_aspect('equal')

pl.show()

