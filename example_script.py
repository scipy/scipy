# sample operation script
import numpy as np
import interpolate1d as I
import matplotlib.pyplot as P


N = 10.0
x = np.arange(N)
y = np.sin(x)


## Interpolating in-range data
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


# plot result
print "plotting results"
P.hold(True)
P.plot(newx, y_block)
P.plot(newx, y_linear)
P.plot(newx, y_quad)
P.plot(newx, y_cubic)
P.show()
print "plotted results"

