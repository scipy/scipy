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
p.show()

# 2-d tests - setup scattered data
x = s.rand(50,1)*4-2
y = s.rand(50,1)*4-2
z = x*s.exp(-x**2-y**2)
ti = s.linspace(-2.0,2.0,81)
(XI,YI) = s.meshgrid(ti,ti)

# use RBF
rbf = Rbf(x.flatten(),y.flatten(),z.flatten(),eps=2)
ZI = rbf(XI.flatten(), YI.flatten())
ZI.shape = XI.shape

# plot the result
from enthought.tvtk.tools import mlab
f=mlab.figure(browser=False)
su=mlab.Surf(XI,YI,ZI,ZI,scalar_visibility=True)
f.add(su)
su.lut_type='blue-red'
f.objects[0].axis.z_label='value'
pp = mlab.Spheres(s.c_[x.flatten(), y.flatten(), z.flatten()],radius=0.03)
f.add(pp)
