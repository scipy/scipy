
# information about included functions
from info import __doc__

# primary usage for 1, 2, and N-dimensional interpolation
from interpolate1d import Interpolate1d, interp1d
from interpolate2d import Interpolate2d, interp2d
from interpolateNd import InterpolateNd, interpNd

# basic interpolation routines
# wrapped by interpolate*.py files
from interpolate_wrapper import linear, logarithmic, block, block_average_above

# support for spline interpolation
# wrapped by interpolate*.py files
from fitpack_wrapper import Spline

# wrapper around Fortran implementing Tom's algorithm 526
from algorithm526_wrapper import algorithm526

