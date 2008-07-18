#FIXME : better docstring
"""
Interpolation of 1D data

This module provides several functions and classes for interpolation
and extrapolation of 1D data (1D in both input and output).  The
primary function provided is:

    interp1d(x, y, new_x) : from data points x and y, interpolates
                                    values for points in new_x and
                                    returns them as an array.

Classes provided include:

    Interpolate1D  :   an object for interpolation of
                            various kinds.  interp1d is a wrapper
                            around this class.
                            
    Spline : an object for spline interpolation
    
Functions provided include:

    linear : linear interpolation
    logarithmic :  logarithmic interpolation
    block : block interpolation
    block_average_above : block average above interpolation

"""

from interpolate_wrapper import linear, logarithmic, block, block_average_above
from fitpack_wrapper import Spline
from Interpolate1D import Interpolate1D, interp1d