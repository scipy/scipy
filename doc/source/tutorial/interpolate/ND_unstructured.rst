.. _tutorial-interpolate_NDunstructured:

.. currentmodule:: scipy.interpolate

===============================================
Scattered data interpolation (:func:`griddata`)
===============================================

Suppose you have multidimensional data, for instance, for an underlying
function :math:`f(x, y)` you only know the values at points ``(x[i], y[i])``
that do not form a regular grid.

.. plot::
    :alt: " "

    Suppose we want to interpolate the 2-D function

    >>> def func(x, y):
    ...     return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

    on a grid in [0, 1]x[0, 1]

    >>> grid_x, grid_y = np.meshgrid(np.linspace(0, 1, 100),
    ...                              np.linspace(0, 1, 200), indexing='ij')

    but we only know its values at 1000 data points:

    >>> rng = np.random.default_rng()
    >>> points = rng.random((1000, 2))
    >>> values = func(points[:,0], points[:,1])

    This can be done with `griddata` -- below, we try out all of the
    interpolation methods:

    >>> from scipy.interpolate import griddata
    >>> grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
    >>> grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
    >>> grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')

    One can see that the exact result is reproduced by all of the
    methods to some degree, but for this smooth function the piecewise
    cubic interpolant gives the best results (black dots show the data being
    interpolated):

    >>> import matplotlib.pyplot as plt
    >>> plt.subplot(221)
    >>> plt.imshow(func(grid_x, grid_y).T, extent=(0, 1, 0, 1), origin='lower')
    >>> plt.plot(points[:, 0], points[:, 1], 'k.', ms=1)   # data
    >>> plt.title('Original')
    >>> plt.subplot(222)
    >>> plt.imshow(grid_z0.T, extent=(0, 1, 0, 1), origin='lower')
    >>> plt.title('Nearest')
    >>> plt.subplot(223)
    >>> plt.imshow(grid_z1.T, extent=(0, 1, 0, 1), origin='lower')
    >>> plt.title('Linear')
    >>> plt.subplot(224)
    >>> plt.imshow(grid_z2.T, extent=(0, 1, 0, 1), origin='lower')
    >>> plt.title('Cubic')
    >>> plt.gcf().set_size_inches(6, 6)
    >>> plt.show()

For each interpolation method, this function delegates to a corresponding
class object --- these classes can be used directly as well ---
`NearestNDInterpolator`, `LinearNDInterpolator` and `CloughTocher2DInterpolator`                                                     
for piecewise cubic interpolation in 2D.

All these interpolation methods rely on triangulation of the data using the
``QHull`` library wrapped in `scipy.spatial`.

.. note::

    `griddata` is based on triangulation, hence is appropriate for unstructured,
    scattered data. If your data is on a full grid,  the `griddata` function ---
    despite its name --- is not the right tool. Use `RegularGridInterpolator`
    instead.



.. _tutorial-interpolate_RBF:

========================================================
Using radial basis functions for smoothing/interpolation
========================================================

Radial basis functions can be used for smoothing/interpolating scattered
data in N dimensions, but should be used with caution for extrapolation
outside of the observed data range.

1-D Example
-----------

This example compares the usage of the `Rbf` and `UnivariateSpline` classes
from the scipy.interpolate module.

.. plot::
    :alt: " "

    >>> import numpy as np
    >>> from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
    >>> import matplotlib.pyplot as plt

    >>> # setup data
    >>> x = np.linspace(0, 10, 9)
    >>> y = np.sin(x)
    >>> xi = np.linspace(0, 10, 101)

    >>> # use fitpack2 method
    >>> ius = InterpolatedUnivariateSpline(x, y)
    >>> yi = ius(xi)

    >>> plt.subplot(2, 1, 1)
    >>> plt.plot(x, y, 'bo')
    >>> plt.plot(xi, yi, 'g')
    >>> plt.plot(xi, np.sin(xi), 'r')
    >>> plt.title('Interpolation using univariate spline')

    >>> # use RBF method
    >>> rbf = Rbf(x, y)
    >>> fi = rbf(xi)

    >>> plt.subplot(2, 1, 2)
    >>> plt.plot(x, y, 'bo')
    >>> plt.plot(xi, fi, 'g')
    >>> plt.plot(xi, np.sin(xi), 'r')
    >>> plt.title('Interpolation using RBF - multiquadrics')
    >>> plt.show()

..   :caption: Example of a 1-D RBF interpolation.

2-D Example
-----------

This example shows how to interpolate scattered 2-D data:

.. plot::
    :alt: " "

    >>> import numpy as np
    >>> from scipy.interpolate import Rbf
    >>> import matplotlib.pyplot as plt

    >>> # 2-d tests - setup scattered data
    >>> rng = np.random.default_rng()
    >>> x = rng.random(100)*4.0-2.0
    >>> y = rng.random(100)*4.0-2.0
    >>> z = x*np.exp(-x**2-y**2)
    >>> edges = np.linspace(-2.0, 2.0, 101)
    >>> centers = edges[:-1] + np.diff(edges[:2])[0] / 2.
    >>> XI, YI = np.meshgrid(centers, centers)

    >>> # use RBF
    >>> rbf = Rbf(x, y, z, epsilon=2)
    >>> ZI = rbf(XI, YI)

    >>> # plot the result
    >>> plt.subplot(1, 1, 1)
    >>> X_edges, Y_edges = np.meshgrid(edges, edges)
    >>> lims = dict(cmap='RdBu_r', vmin=-0.4, vmax=0.4)
    >>> plt.pcolormesh(X_edges, Y_edges, ZI, shading='flat', **lims)
    >>> plt.scatter(x, y, 100, z, edgecolor='w', lw=0.1, **lims)
    >>> plt.title('RBF interpolation - multiquadrics')
    >>> plt.xlim(-2, 2)
    >>> plt.ylim(-2, 2)
    >>> plt.colorbar()
