
# information about included functions
from info import __doc__

# basic interpolation routines
from interpolate_wrapper import linear, logarithmic, block, block_average_above

# support for spline interpolation
from fitpack_wrapper import Spline

# wrapping for all supported interpolation types.
from interpolate1d import Interpolate1d, interp1d