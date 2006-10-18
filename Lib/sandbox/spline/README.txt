This module is meant to replace the spline functionality of scipy.interpolate.
At the moment the code base is identical with a few minor cosmetic changes.
It does not compile yet!

Planned changes are:

1. remove dependence of f2c generated interface (use only f2py)
2. cleanup!
3. rationalise the interface (particularly the new spline/fitpack2 interface)
4. build comprehensive unit test suite

Hopefully this module will end up in scipy, with scipy.interpolate only then
containing the interpolation functions interp1d, interp2d etc. which may depend
on the spline module and the new delaunay module (for scattered data).

John Travers
