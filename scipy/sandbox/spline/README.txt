This module is meant to replace the spline functionality of scipy.interpolate.

Changes so far:
    1. removed c wrappers - only usig f2py interface now
    2. several bugs fixed (interface to percur, size of c array in tck)
    3. added 11 unit tests - now fully covered
    3. basic cleanup of the module files and docs

Planned changes are:
    1. possibly rationalise the interface to spline
    2. add further functions from fitpack fortran routines
    3. broaden documentation
    4. more tests

Hopefully this module will end up in scipy, with scipy.interpolate only then
containing the interpolation functions interp1d, interp2d etc. which may depend
on the spline module and the new delaunay module (for scattered data).

John Travers