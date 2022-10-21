========================================
Interpolation (:mod:`scipy.interpolate`)
========================================

.. sectionauthor:: Travis E. Oliphant

.. sectionauthor:: Pauli Virtanen

.. currentmodule:: scipy.interpolate

.. contents::

There are several general facilities available in SciPy for interpolation and
smoothing for data in 1, 2, and higher dimensions. The choice of a specific
interpolation routine depends on whether your data is

- :ref:`One-dimensional <tutorial-interpolate_1Dsection>` (1-D);

- Is given on a :ref:`structured grid <tutorial-interpolate_regular_grid_interpolator>`
  in arbitrary (N) dimensions;

- Is :ref:`unstructured <tutorial-interpolate_NDunstructured>`;

For data smoothing, :ref:`functions are provided <tutorial-interpolate_fitpack>`
for 1- and 2-D (smoothed)
cubic-spline interpolation, based on the FORTRAN library FITPACK. There are
two interfaces for the FITPACK library, a procedural interface and an
object-oriented interface.

Additionally, routines are provided for interpolation / smoothing using
:ref:`radial basis functions <tutorial-interpolate_RBF>` with several kernels.

.. toctree::
   :maxdepth: 3

   interpolate/splines_and_polynomials
   interpolate/smoothing_splines


Missing data
============

Before describing various interpolators in detail, we note the (lack of) support
for interpolation with missing data. Two popular ways of representing missing
data are using masked arrays of the `numpy.ma` library, and encoding missing
values as not-a-number, ``NaN``. 

Neither of these two approaches are directly suppored in `scipy.interpolate`.
Individual routines may offer partial support, and/or workarounds, but in
general the library firmly adheres to the IEEE 754 semantics where a ``NaN``
means *not-a-number*, i.e. a result of an illegal mathematical operation
(think a division by zero), not *missing*.


.. _tutorial-interpolate_1Dsection:

1-D interpolation
=================

Piecewise linear interpolation
------------------------------

If all you need is a linear (a.k.a. broken line) interpolation, you can use
the `numpy.interp` routine. It takes two arrays of data to interpolate, ``x``,
and ``y``, and a third array, ``xnew``, of points to evaluate the interpolation on:


.. plot::

   >>> import numpy as np

   >>> x = np.linspace(0, 10, num=11)
   >>> y = np.cos(-x**2 / 9.0)

   Construct the interpolation

   >>> xnew = np.linspace(0, 10, num=1001)
   >>> ynew = np.interp(xnew, x, y)

   And plot it

   >>> import matplotlib.pyplot as plt
   >>> plt.plot(xnew, ynew, '-', label='linear interp')
   >>> plt.plot(x, y, 'o', label='data')
   >>> plt.legend(loc='best')
   >>> plt.show()

..   :caption: One-dimensional interpolation using `numpy.interp`


Cubic splines
-------------

Of course, piecewise linear interpolation produces corners at data points,
where linear pieces join. To produce a smoother curve, you can use cubic
splines, where the interpolating curve is made of cubic pieces with matching
first and second derivatives. In code, these objects are represented via the
``CubicSpline`` class instances. An instance is constructed with the ``x`` and
``y`` arrays of data, and then it can be evaluated using the target ``xnew``
values:

    >>> from scipy.interpolate import CubicSpline
    >>> spl = CubicSpline([1, 2, 3, 4, 5, 6], [1, 4, 8, 16, 25, 36])
    >>> spl(2.5)
    5.57

A `CubicSpline` object's ``__call__`` method accepts both scalar values and
arrays. It also accepts a second argument, ``nu``, to evaluate the 
derivative of order ``nu``. As an example, we plot the derivatives of a spline:

.. plot::

    >>> from scipy.interpolate import CubicSpline
    >>> x = np.linspace(0, 10, num=11)
    >>> y = np.cos(-x**2 / 9.)
    >>> spl = CubicSpline(x, y)

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(4, 1, figsize=(5, 7))
    >>> xnew = np.linspace(0, 10, num=1001)
    >>> ax[0].plot(xnew, spl(xnew))
    >>> ax[0].plot(x, y, 'o', label='data')
    >>> ax[1].plot(xnew, spl(xnew, nu=1), '--', label='1st derivative')
    >>> ax[2].plot(xnew, spl(xnew, nu=2), '--', label='2nd derivative')
    >>> ax[3].plot(xnew, spl(xnew, nu=3), '--', label='3rd derivative')
    >>> for j in range(4):
    ...     ax[j].legend(loc='best')
    >>> plt.tight_layout()
    >>> plt.show()

Note that the first and second derivatives are continuous by construction, and
the third derivative jumps at data points. 


Monotone interpolants
---------------------

Cubic splines are by construction twice continuously differentiable. This may
lead to the spline function oscillating and ''overshooting'' in between the
data points. In these situations, an alternative is to use the so-called
*monotone* cubic interpolants: these are constructed to be only once
continuously differentiable, and attempt to preserve the local shape implied
by the data. There are two objects of this class in `scipy.interpolate` :
`PchipInterpolator` and `Akima1DInterpolator` . To illustrate, let's consider
data with an outlier:

.. plot::

    >>> from scipy.interpolate import CubicSpline, PchipInterpolator, Akima1DInterpolator
    >>> x = np.array([1., 2., 3., 4., 4.5, 5., 6., 7., 8])
    >>> y = x**2
    >>> y[4] += 101

    >>> import matplotlib.pyplot as plt
    >>> xx = np.linspace(1, 8, 51)
    >>> plt.plot(xx, CubicSpline(x, y)(xx), '--', label='spline')
    >>> plt.plot(xx, Akima1DInterpolator(x, y)(xx), '-', label='Akima1D')
    >>> plt.plot(xx, PchipInterpolator(x, y)(xx), '-', label='pchip')
    >>> plt.plot(x, y, 'o')
    >>> plt.legend()
    >>> plt.show()



.. _tutorial-interpolate_interp1d:

Legacy interface for 1-D interpolation (:class:`interp1d`)
----------------------------------------------------------

.. note::
    `interp1d` is considered legacy API and is not recommended for use in new
    code. Consider using more specific interpolators instead. 

The `interp1d` class in `scipy.interpolate` is a convenient method to
create a function based on fixed data points, which can be evaluated
anywhere within the domain defined by the given data using linear
interpolation. An instance of this class is created by passing the 1-D
vectors comprising the data. The instance of this class defines a
``__call__`` method and can therefore by treated like a function which
interpolates between known data values to obtain unknown values.
Behavior at the boundary can be
specified at instantiation time. The following example demonstrates
its use, for linear and cubic spline interpolation:

.. plot::
   :alt: "This code generates an X-Y plot of a time-series with amplitude on the Y axis and time on the X axis. The original time-series is shown as a series of blue markers roughly defining some kind of oscillation. An orange trace showing the linear interpolation is drawn atop the data forming a jagged representation of the original signal. A dotted green cubic interpolation is also drawn that appears to smoothly represent the source data."

   >>> from scipy.interpolate import interp1d

   >>> x = np.linspace(0, 10, num=11, endpoint=True)
   >>> y = np.cos(-x**2/9.0)
   >>> f = interp1d(x, y)
   >>> f2 = interp1d(x, y, kind='cubic')

   >>> xnew = np.linspace(0, 10, num=41, endpoint=True)
   >>> import matplotlib.pyplot as plt
   >>> plt.plot(x, y, 'o', xnew, f(xnew), '-', xnew, f2(xnew), '--')
   >>> plt.legend(['data', 'linear', 'cubic'], loc='best')
   >>> plt.show()

..   :caption: One-dimensional interpolation using the
..             class :obj:`interpolate.interp1d` with
..             kind equals `linear` and `cubic`.


Another set of interpolations in `interp1d` is `nearest`, `previous`, and
`next`, where they return the nearest, previous, or next point along the
x-axis. Nearest and next can be thought of as a special case of a causal
interpolating filter. The following example demonstrates their use, using the
same data as in the previous example:

.. plot::
   :alt: "This code generates an X-Y plot of a time-series with amplitude on the Y axis and time on the X axis. The original time-series is shown as a series of blue markers roughly defining some kind of oscillation. An orange trace showing the nearest neighbor interpolation is drawn atop the original with a stair-like appearance where the original data is right in the middle of each stair step. A green trace showing the previous neighbor interpolation looks similar to the orange trace but the original data is at the back of each stair step. Similarly a dotted red trace showing the next neighbor interpolation goes through each of the previous points, but it is centered at the front edge of each stair."

   >>> from scipy.interpolate import interp1d

   >>> x = np.linspace(0, 10, num=11, endpoint=True)
   >>> y = np.cos(-x**2/9.0)
   >>> f1 = interp1d(x, y, kind='nearest')
   >>> f2 = interp1d(x, y, kind='previous')
   >>> f3 = interp1d(x, y, kind='next')

   >>> xnew = np.linspace(0, 10, num=1001, endpoint=True)
   >>> import matplotlib.pyplot as plt
   >>> plt.plot(x, y, 'o')
   >>> plt.plot(xnew, f1(xnew), '-', xnew, f2(xnew), '--', xnew, f3(xnew), ':')
   >>> plt.legend(['data', 'nearest', 'previous', 'next'], loc='best')
   >>> plt.show()

..   :caption: One-dimensional interpolation using the
..             class :obj:`interpolate.interp1d` with
..             kind equals `nearest`, `previous`, and
..             `next`.


.. _tutorial-interpolate_regular_grid_interpolator:

Multivariate data interpolation on a regular grid  (:class:`RegularGridInterpolator`)
======================================================================================

Suppose you have n-dimensional data on a regular grid, and you want to interpolate it.
In such a case, :class:`RegularGridInterpolator` can be useful.

Strictly speaking, this class efficiently handles data given on *rectilinear*
grids: hypercubic lattices with possibly unequal spacing between points.

The following example demonstrates its use, and compares the interpolation results
using each method.

.. plot::

   >>> import matplotlib.pyplot as plt
   >>> from scipy.interpolate import RegularGridInterpolator

   Suppose we want to interpolate this 2-D function.

   >>> def F(u, v):
   ...     return u * np.cos(u * v) + v * np.sin(u * v)

   Suppose we only know some data on a regular grid.

   >>> fit_points = [np.linspace(0, 3, 8), np.linspace(0, 3, 8)]
   >>> values = F(*np.meshgrid(*fit_points, indexing='ij'))

   Creating test points and true values for evaluations.

   >>> ut, vt = np.meshgrid(np.linspace(0, 3, 80), np.linspace(0, 3, 80), indexing='ij')
   >>> true_values = F(ut, vt)
   >>> test_points = np.array([ut.ravel(), vt.ravel()]).T

   We can creat interpolator and interpolate test points using each method.

   >>> interp = RegularGridInterpolator(fit_points, values)
   >>> fig, axes = plt.subplots(2, 3, figsize=(10, 6))
   >>> axes = axes.ravel()
   >>> fig_index = 0
   >>> for method in ['linear', 'nearest', 'slinear', 'cubic', 'quintic']:
   ...     im = interp(test_points, method=method).reshape(80, 80)
   ...     axes[fig_index].imshow(im)
   ...     axes[fig_index].set_title(method)
   ...     axes[fig_index].axis("off")
   ...     fig_index += 1
   >>> axes[fig_index].imshow(true_values)
   >>> axes[fig_index].set_title("True values")
   >>> fig.tight_layout()
   >>> fig.show()

   As expected, the higher degree spline interpolations are closest to the
   true values, though are more expensive to compute than with `linear`
   or `nearest`. The `slinear` interpolation also matches the `linear`
   interpolation.


If you prefer a funcional interface to explicitly creating a class instance,
the `interpn` convenience function offers the equivalent functionality.

Specifically, these two forms give identical results:

    >>> from scipy.interpolate import interpn
    >>> rgi = RegularGridInterpolator(fit_points, values)
    >>> result_rgi = rgi(test_points)

and

    >>> result_interpn = interpn(fit_points, values, test_points)
    >>> np.allclose(result_rgi, result_interpn, atol=1e-15)
    True


Finally, we note that if you are dealing with data on Cartesian grids with
integer coordinates, e.g. resampling image data, these routines may not be the
optimal choice. Consider using `scipy.ndimage.map_coordinates` instead.


.. _tutorial-interpolate_NDunstructured:

Scattered data interpolation (:func:`griddata`)
===============================================

Suppose you have multidimensional data, for instance, for an underlying
function *f(x, y)* you only know the values at points ``(x[i], y[i])``
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

    `griddata` is based on triangulation, hence is appropriate for unstructures,
    scattered data. If your data is on a full grid,  the `griddata` function ---
    despite its name --- is not the right tool. Use `RegularGridInterpolator`
    instead.



.. _tutorial-interpolate_RBF:

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
    >>> from matplotlib import cm

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
