.. _tutorial-interpolate_regular_grid_interpolator:

.. currentmodule:: scipy.interpolate

=====================================================================================
Multivariate data interpolation on a regular grid  (:class:`RegularGridInterpolator`)
=====================================================================================

Suppose you have N-dimensional data on a regular grid, and you
want to interpolate it. In such a case, :class:`RegularGridInterpolator` can be
useful. Several interpolation strategies are supported: nearest-neighbor,
linear, and tensor product splines of odd degree.

Strictly speaking, this class efficiently handles data given on *rectilinear*
grids: hypercubic lattices with possibly unequal spacing between points.
The number of points per dimension can be different for different dimensions.

The following example demonstrates its use, and compares the interpolation results
using each method.

.. plot::

   >>> import matplotlib.pyplot as plt
   >>> from scipy.interpolate import RegularGridInterpolator

   Suppose we want to interpolate this 2-D function.

   >>> def F(u, v):
   ...     return u * np.cos(u * v) + v * np.sin(u * v)

   Suppose we only know some data on a regular grid.

   >>> fit_points = [np.linspace(0, 3, 8), np.linspace(0, 3, 11)]
   >>> values = F(*np.meshgrid(*fit_points, indexing='ij'))

   Creating test points and true values for evaluations.

   >>> ut, vt = np.meshgrid(np.linspace(0, 3, 80), np.linspace(0, 3, 80), indexing='ij')
   >>> true_values = F(ut, vt)
   >>> test_points = np.array([ut.ravel(), vt.ravel()]).T

   We can create the interpolator and interpolate test points using each method.

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

If your data is such that spline methods produce ringing, you may consider
using `method="pchip"`, which uses the tensor product of PCHIP interpolators,
a `PchipInterpolator` per dimension.

If you prefer a functional interface opposed to explicitly creating a class instance,
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

