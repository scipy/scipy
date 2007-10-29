import scipy as s
import scipy.interpolate

from scipy.sandbox.rbf import Rbf

import matplotlib
matplotlib.use('Agg')
import pylab as p

# 1d tests - setup data
x = s.linspace(0,10,9)
y = s.sin(x)
xi = s.linspace(0,10,101)

# use interpolate methods
ius = s.interpolate.InterpolatedUnivariateSpline(x,y)
yi = ius(xi)
p.subplot(2,1,1)
p.plot(x,y,'o',xi,yi, xi, s.sin(xi),'r')
p.title('Interpolation using current scipy fitpack2')

# use RBF method
rbf = Rbf(x, y)
fi = rbf(xi)
p.subplot(2,1,2)
p.plot(x,y,'bo',xi,fi,'g',xi, s.sin(xi),'r')
p.title('RBF interpolation - multiquadrics')
p.savefig('rbf1d.png')
p.close()

# 2-d tests - setup scattered data
x = s.rand(100)*4.0-2.0
y = s.rand(100)*4.0-2.0
z = x*s.exp(-x**2-y**2)
ti = s.linspace(-2.0,2.0,100)
(XI,YI) = s.meshgrid(ti,ti)

# use RBF
rbf = Rbf(x,y,z,epsilon=2)
ZI = rbf(XI, YI)

# plot the result
n = p.normalize(-2., 2.)
p.subplot(1,1,1)
p.pcolor(XI,YI,ZI,cmap=p.cm.jet)
p.scatter(x,y,100,z,cmap=p.cm.jet)
p.title('RBF interpolation - multiquadrics')
p.xlim(-2,2)
p.ylim(-2,2)
p.colorbar()
p.savefig('rbf2d.png')
p.close()
