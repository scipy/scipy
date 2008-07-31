""" sample operation script

    Creates a sample dataset, performs several
    basic interpolations on it, and plots results
    for comparison.  Pauses, then does the same
    thing using different user-input options.

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
    interp = I.Interpolate1d(x, y, 'quartic')
    y_quartic = interp(newx)

    # plot result
    print "plotting results"
    P.hold(True)
    P.plot(newx, y_block, 'g')
    P.plot(newx, y_linear, 'b')
    P.plot(newx, y_quad, 'r')
    P.plot(newx, y_cubic, 'm')
    P.plot(newx, y_quartic, 'y')
    P.title( "interpolating in-range data with Interpolate1d class" )
    P.show()
    print "plotted results"
    
    time.sleep(5)


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

    # 4th order spline
    interp = I.Interpolate1d(x, y, 'Quartic')
    y_quartic2 = interp(newx)
    
    # 5th order spline
    interp = I.Interpolate1d(x, y, 'Quintic')
    y_quintic2 = interp(newx)

    # plot result
    print "plotting results"
    P.hold(True)
    P.plot(newx, y_block2, 'g')
    P.plot(newx, y_linear2, 'b')
    P.plot(newx, y_quad2, 'r')
    P.plot(newx, y_quartic2, 'y')
    P.plot(newx, y_quintic2, 'y')
    P.title( "same data through different interface" )
    P.show()
    print "plotted results"
        

# demoing block_average_above and logarithmic
if False:
    N = 10.0
    x = np.arange(N)
    x[1] = 1.2 # make it grid non-regular
    y = np.sin(x)    
    newx = np.arange(.05, N, .05)

    # block interpolation
    # FIXME : I'm not really familiar with logarithmic
    #   interpolation and thus can't comment on what these
    #   data should look like, but this looks weird.
    interp = I.Interpolate1d(x, y, 'logarithmic')
    y_logarithmic = interp(newx)

    # linear interpolation
    interp = I.Interpolate1d(x, y, 'block_average_above')
    y_blockavgabove = interp(newx)
    
    # plotting the results
    P.hold(true)
    P.plot(newx, y_logarithmic, 'g')
    P.plot(newx, y_blockavgabove, 'r')
    P.show()