README


This file contains information about the architecture (not the api) of the module
and license information.  It is designed mostly for developers who need to
understand the organization of the files, extensions, etc and the logic (or lack thereof)
behind them so that they can continue development of the module.

The key the understanding the module is to understand where it came from.
Interpolation code in a variety of different forms for a variety of different purposes
had been spread throughout the various parts of Enthought code.  The goal was to
collect these various pieces of code and combine them into one module with a
reasonable, intuitive api.  But under the hood, there are several extension
modules whose functionalities overlap or do not naturally dovetail.

Main Files:
interpolate1d.py : 
    Contains the Interpolate1d class and the interp1d functional wrapper
    around it.  This mostly provides a user-interface; the user passes in
    keywords, and according to those keywords, Interpolate1d will call functions
    and classes in the wrapper files.
interpolate2d.py : 
    Completely analogous to interpolate1d.py.  A user interface that calls machinery
    in the wrapper files.  Its organization is also almost completely analogous to
    that of interpolate1d.py.
interpolateNd.py : 
    This file doubles as a user interface and wrapper file around the _nd_image
    extension.  Interpolate1d and 2d are both operated by inputting lists of
    points, which is pretty generic and lends itself to removing bad data etc.
    But the _nd_image extension requires a uniform grid, so this file only performs
    interpolation with the _nd_image extension.  _nd_image interpolates array entries,
    so this file handles 1) formatting, and 2) translations between spatial coordinates and
    array indices.
    
    The spline filtering step is carried over from ndimage, and is necessary to make
    the extension module work properly.  In the future, it could be good to modify
    the _nd_image module to include this annoying step itself and hide it from the
    developer.

Wrapper files:
fitpack_wrapper.py :
    This file provides the Spline and Spline2d classes which provide a variety of 1d and
    2d spline functions.  It is these classes which are accessed by Interpolate1d and
    Interpolate2d, but note that only part of their functionality is accessed.  Things like
    smoothing and seeing spline coefficients are supported by Spline and Spline2d but to
    exported.  Internally, these classes call the _dfitpack extension.
    
    This is based on code that was in scipy.interpolate.  Much of the functionality that was
    in scipy.interpolate is not reproduced here.
    
interpolate_wrapper.py :
    A variety of 1d interpolation functions.  Most of them are just wrappers around _interpolate,
    but others are stand-alone.
    
algorithm526_wrapper.py :
    very simple interface to the _interp_526 module.  The main reason for the new file is so
    that the imported function is recognized by Interpolate2d as a function, rather than a
    Fortran object so that _init_xyz works.

Extensions:
_dfitpack :
    Fortran extension module.  This wraps part of the functionality of the fitpack library.
_interpolate :
    C extension module with basic functions.
_nd_image : 
    C extensions module.  It uses spline interpolation of various orders to interpolate entries
    in an array.
_interp_526 :
    Fortran extension module implementing TOMS Algorithm 526.  Taken from the Pore Pressure
    project.

