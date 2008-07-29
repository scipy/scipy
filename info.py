# FIXME : better docstring.  This needs updating as features change,
#       and it also discusses technical points as well as user-interface.

__doc__ = \
"""
    Interpolation of 1D data

    This module provides several functions and classes for interpolation
    and extrapolation of 1D (in both input and output) real-valued.  The
    primary function provided is:

        interp1d(x, y, new_x) : from data points (x[i], y[i]), interpolates
                                        values for points in new_x and
                                        returns them as an array.  x and new_x
                                        must both be in sorted order.

    Classes provided include:

        Interpolate1d  :   an object for interpolation of
                                various kinds.  interp1d is a wrapper
                                around this class.
                                
        Spline : an object for spline interpolation.  Interpolate1d
                                wraps this class if spline interpolation
                                is used.  However, not all functionality
                                of Spline is available through Interpolate1d.
        
    Functions provided include:

        linear : linear interpolation
        logarithmic :  logarithmic interpolation
        block : block interpolation
        block_average_above : block average above interpolation
        
    The dependency/interface architecture is as follows:
        interpolate1d.py is viewed by the user (through interp1d and Interpolate1d)
            It depends on fitpack_wrapper.py and interpolate_wrapper.py
        fitpack_wrapper is viewed by the user (through Spline)
            It depends on dfitpack.pyd, a Fortran extension module
        interpolate_wrapper is viewed by he user (through functions linear, etc)
            It depends on _interpolate.pyd, a C extension module.

"""
        
postpone_import = 1