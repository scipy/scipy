#!/usr/bin/env python
# Created by Pearu Peterson, Aug 2003
""" Test xplt based demos for interpolate.fitpack2 module
"""
__usage__ = """
Build interpolate:
  python setup_interpolate.py build
Run demos (assumes that scipy is installed):
  python -i tests/demos_xplt.py
"""

import sys
from scipy_test.testing import set_package_path
set_package_path()
from interpolate.fitpack2 import UnivariateSpline,LSQUnivariateSpline,\
     InterpolatedUnivariateSpline
from interpolate.fitpack2 import LSQBivariateSpline, SmoothBivariateSpline
del sys.path[0]

from scipy import *

def demo1():
    x = arange(0,2*pi+pi/4,2*pi/8)
    xnew = arange(-pi/10,2*pi+pi/4+pi/10,pi/50)
    y = sin(x)


    def make_plot():
        xplt.plot(x,y,'x',xnew,spline(xnew),x,y,'b',xnew,sin(xnew),
                  spline.get_knots(),spline(spline.get_knots()),'o')        

    spline = UnivariateSpline(x,y,k=1)
    assert isinstance(spline,LSQUnivariateSpline)
    print 'Linear LSQ approximation of sin(x):',spline.__class__.__name__
    make_plot()
    print 'Residual=',spline.get_residual()
    raw_input('Press any key to continue..')

    spline.set_smoothing_factor(0)
    assert isinstance(spline,InterpolatedUnivariateSpline)
    print 'Linear interpolation of sin(x):',spline.__class__.__name__
    make_plot()
    print 'Residual=',spline.get_residual()
    raw_input('Press any key to continue..')

    spline = UnivariateSpline(x,y,k=1,s=0.1)
    print 'Linear smooth approximation of sin(x):',spline.__class__.__name__
    assert isinstance(spline,UnivariateSpline)
    make_plot()
    print 'Residual=',spline.get_residual()
    raw_input('Press any key to continue..')
    
if __name__ == "__main__":
    demo1()
