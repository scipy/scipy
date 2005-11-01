"""Some test functions for bivariate interpolation.

Most of these have been yoinked from ACM TOMS 792.
http://netlib.org/toms/792
"""

import scipy as sp
from triangulate import Triangulation

def exponential(x, y):
    x = x*9
    y = y*9
    x1 = x+1.0
    x2 = x-2.0
    x4 = x-4.0
    x7 = x-7.0
    y1 = x+1.0
    y2 = y-2.0
    y3 = y-3.0
    y7 = y-7.0
    f = (0.75 * sp.exp(-(x2*x2+y2*y2)/4.0) +
         0.75 * sp.exp(-x1*x1/49.0 - y1/10.0) +
         0.5 * sp.exp(-(x7*x7 + y3*y3)/4.0) - 
         0.2 * sp.exp(-x4*x4 -y7*y7))
    return f
exponential.title = 'Exponential and Some Gaussians'

def cliff(x, y):
    f = sp.tanh(9.0*(y-x) + 1.0)/9.0
    return f
cliff.title = 'Cliff'

def saddle(x, y):
    f = (1.25 + sp.cos(5.4*y))/(6.0 + 6.0*(3*x-1.0)**2)
    return f
saddle.title = 'Saddle'

def gentle(x, y):
    f = sp.exp(-5.0625*((x-0.5)**2+(y-0.5)**2))/3.0
    return f
gentle.title = 'Gentle Peak'

def steep(x, y):
    f = sp.exp(-20.25*((x-0.5)**2+(y-0.5)**2))/3.0
    return f
steep.title = 'Steep Peak'

def sphere(x, y):
    circle = 64-81*((x-0.5)**2 + (y-0.5)**2)
    f = sp.where(circle >= 0, sp.sqrt(sp.clip(circle,0,100)) - 0.5, 0.0)
    return f
sphere.title = 'Sphere'

def trig(x, y):
    f = 2.0*sp.cos(10.0*x)*sp.sin(10.0*y) + sp.sin(10.0*x*y)
    return f
trig.title = 'Cosines and Sines'

def gauss(x, y):
    x = 5.0-10.0*x
    y = 5.0-10.0*y
    g1 = sp.exp(-x*x/2)
    g2 = sp.exp(-y*y/2)
    f = g1 + 0.75*g2*(1 + g1)
    return f
gauss.title = 'Gaussian Peak and Gaussian Ridges'

def cloverleaf(x, y):
    ex = sp.exp((10.0-20.0*x)/3.0)
    ey = sp.exp((10.0-20.0*y)/3.0)
    logitx = 1.0/(1.0+ex)
    logity = 1.0/(1.0+ey)
    f = (((20.0/3.0)**3 * ex*ey)**2 * (logitx*logity)**5 * 
        (ex-2.0*logitx)*(ey-2.0*logity))
    return f
cloverleaf.title = 'Cloverleaf'

def cosine_peak(x, y):
    circle = sp.hypot(80*x-40.0, 90*y-45.)
    f = sp.exp(-0.4*circle) * sp.cos(0.15*circle)
    return f
cosine_peak.title = 'Cosine Peak'

allfuncs = [exponential, cliff, saddle, gentle, steep, sphere, trig, gauss, cloverleaf, cosine_peak]


class LinearTester(object):
    name = 'Linear'
    def __init__(self, xrange=(-0.5, 1.5), yrange=(-0.5, 1.5), nrange=101, npoints=250):
        self.xrange = xrange
        self.yrange = yrange
        self.nrange = nrange
        self.npoints = npoints

        rng = sp.random.RandomState(1234567890)
        self.x = rng.uniform(xrange[0], xrange[1], size=npoints)
        self.y = rng.uniform(yrange[0], yrange[1], size=npoints)
        self.tri = Triangulation(self.x, self.y)

    def interpolator(self, func):
        z = func(self.x, self.y)
        return self.tri.linear_interpolator(z)

    def plot(self, func, interp=True, plotter='imshow', margin=0.4):
        import matplotlib as mpl
        from matplotlib import pylab as pl
        if interp:
            lpi = self.interpolator(func)
            z = lpi[self.yrange[0]+margin:self.yrange[1]-margin:complex(0,self.nrange),
                    self.xrange[0]+margin:self.xrange[1]-margin:complex(0,self.nrange)]
        else:
            y, x = sp.mgrid[self.yrange[0]+margin:self.yrange[1]-margin:complex(0,self.nrange),
                            self.xrange[0]+margin:self.xrange[1]-margin:complex(0,self.nrange)]
            z = func(x, y)

        extent = (self.xrange[0]+margin, self.xrange[1]-margin,
            self.yrange[0]+margin, self.yrange[1]-margin)
        pl.ioff()
        pl.clf()
        pl.hot() # Some like it hot
        if plotter == 'imshow':
            pl.imshow(z, interpolation='nearest', extent=extent, origin='lower', aspect='preserve')
        elif plotter == 'contour':
            Y, X = sp.ogrid[self.yrange[0]+margin:self.yrange[1]-margin:complex(0,self.nrange),
                self.xrange[0]+margin:self.xrange[1]-margin:complex(0,self.nrange)]
            pl.contour(sp.ravel(X), sp.ravel(Y), z, 20)
        x = self.x
        y = self.y
        lc = mpl.collections.LineCollection(sp.array([((x[i], y[i]), (x[j], y[j]))
            for i, j in self.tri.edge_db]), colors=[(0,0,0,0.2)])
        pl.gca().add_collection(lc)
        if interp:
            title = '%s Interpolant' % self.name
        else:
            title = 'Reference'
        if hasattr(func, 'title'):
            pl.title('%s: %s' % (func.title, title))
        else:
            pl.title(title)

        pl.show()
        pl.ion()

class NNTester(LinearTester):
    name = 'Natural Neighbors'
    def interpolator(self, func):
        z = func(self.x, self.y)
        return self.tri.nn_interpolator(z)

def plotallfuncs(allfuncs=allfuncs):
    from matplotlib import pylab as pl
    pl.ioff()
    nnt = NNTester(npoints=1000)
    lpt = LinearTester(npoints=1000)
    for func in allfuncs:
        print func.title
        nnt.plot(func, interp=False, plotter='imshow')
        pl.savefig('%s-ref-img.png' % func.func_name)
        nnt.plot(func, interp=True, plotter='imshow')
        pl.savefig('%s-nn-img.png' % func.func_name)
        lpt.plot(func, interp=True, plotter='imshow')
        pl.savefig('%s-lin-img.png' % func.func_name)
        nnt.plot(func, interp=False, plotter='contour')
        pl.savefig('%s-ref-con.png' % func.func_name)
        nnt.plot(func, interp=True, plotter='contour')
        pl.savefig('%s-nn-con.png' % func.func_name)
        lpt.plot(func, interp=True, plotter='contour')
        pl.savefig('%s-lin-con.png' % func.func_name)
    pl.ion()
