""" sample operation script
    Note that in the plot, quadratic, cubic and
    quintic lines blur together.  You can comment
    two out to see one clearly.
"""
import numpy as np
import interpolate1d as I
import matplotlib.pyplot as P
import interpolate_wrapper
import fitpack_wrapper
import time


## Interpolating in-range data.  Basic operation
if True:
    
    N = 10.0
    x = np.arange(N)
    x[1] = 1.2 # make it grid non-regular
    y = np.sin(x)    
    newx = np.arange(.05, N, .05)

    # block interpolation
    interp = I.Interpolate1d(x, y, 'block')
    y_block = interp(newx)

    # linear interpolation
    interp = I.Interpolate1d(x, y, 'linear')
    y_linear = interp(newx)

    # 2nd order spline
    interp = I.Interpolate1d(x, y, 'quad')
    y_quad = interp(newx)

    # 3rd order spline
    interp = I.Interpolate1d(x, y, 'cubic')
    y_cubic = interp(newx)
    
    # 4th order spline
    interp = I.Interpolate1d(x, y, 'quintic')
    y_quintic = interp(newx)

    # plot result
    print "plotting results"
    P.hold(True)
    P.plot(newx, y_block, 'g')
    P.plot(newx, y_linear, 'b')
    P.plot(newx, y_quad, 'r')
    P.plot(newx, y_cubic, 'm')
    P.plot(newx, y_quintic, 'y')
    P.title( "interpolating in-range data with Interpolate1d class" )
    P.show()
    print "plotted results"
    
    time.sleep(3)

## demoing some of the other interfac features
if True:
    N = 10.0
    x = np.arange(N)
    x[1] = 1.2 # make it grid non-regular
    y = np.sin(x)    
    newx = np.arange(.05, N, .05)

    # block interpolation with interp1d class
    y_block2 = I.interp1d(x, y, newx, 'block')

    # linear interpolation with a function passed in
    interp = I.Interpolate1d(x, y, interpolate_wrapper.linear)
    y_linear2 = interp(newx)

    # 2nd order spline, with partially-instantiated class
    interp = I.Interpolate1d(x, y, fitpack_wrapper.Spline(k=2))
    y_quad2 = interp(newx)

    # 3rd order spline with additional keyword arguments
    interp = I.Interpolate1d(x, y, fitpack_wrapper.Spline, kindkw = {'k':3})
    y_cubic2 = interp(newx)

    # 4th order spline
    interp = I.Interpolate1d(x, y, 'quintic')
    y_quintic2 = interp(newx)

    # plot result
    print "plotting results"
    P.hold(True)
    P.plot(newx, y_block2, 'g')
    P.plot(newx, y_linear2, 'b')
    P.plot(newx, y_quad2, 'r')
    P.plot(newx, y_cubic2, 'm')
    P.plot(newx, y_quintic2, 'y')
    P.title( "same data through different interface" )
    P.show()
    print "plotted results"
    
    #time.sleep(20)
    
