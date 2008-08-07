# FIXME : better docstring.  This needs updating as features change,
#       and it also discusses technical points as well as user-interface.

__doc__ = \
"""
    This module provides several functions and classes for interpolation
    and extrapolation of real-valued data.  The primary interface is through
    the functions:

        interp1d(x, y, newx) :  from data points (x[i], y[i]), interpolates
                                        values for points in new_x and
                                        returns them as an array.  x and newx
                                        must both be in sorted order.
        
        interp2d(x, y, z, newx, newy): from data points (x[i], y[i], z[i]), interpolates
                                        values for points (newx[i], newy[i]) and
                                        returns them as an array.
                                        
        interpNd(data, coordinates): data contains known values (which are
                                        assume to be uniformly spaced), and coordinates
                                        indicates where to interpolate new values.
                                        
    Each function defaults to linear interpolation for in-range points, and
    returning a NaN when points are out of range.  However, by specifying keywords
    a variety of techniques can be chosen or, for the first wo functions, even explicitly
    written by the user.
    
    The following callable classes are also provided:

        Interpolate1d  :   an object for interpolation of
                                various kinds.  interp1d is a wrapper
                                around this class.
                                
        Spline : an object for spline interpolation.  Interpolate1d
                                wraps this class if spline interpolation
                                is used.  However, not all functionality
                                of Spline is available through Interpolate1d.
                                
        Interpolate2d  :
        
        Spline 2d  :  
        
        InterpolateNd  :
        
    These functions and classes constitute the primary api.  However, several
    additional functions are also provided (the primary api in many cases calls
    these functions):

        linear : linear interpolation
        logarithmic :  logarithmic interpolation
        block : block interpolation
        block_average_above : block average above interpolation

"""
        
postpone_import = 1