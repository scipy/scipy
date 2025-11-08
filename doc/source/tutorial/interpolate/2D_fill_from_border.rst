.. _tutorial-intepolate_2dfillfromborder:

.. currentmodule:: scipy.interpolate

==============================
Filling an image from a border
==============================

Suppose we have values for a function only on the border of a region,
either to estimate the background of an image from the parts not
covered by the sample, or as part of a boundary-value problem.

In particular, consider a square with sides of length :math:`L`, with
function values on the boundary given by

.. math::

   \sin(6\pi x/L), y = 0 \text{or} y = L
   \sin(6\pi y/L), x = 0 \text{or} x = L

That is, a sine wave with three full periods on each side, identical
on opposite sides of the square.

------------------
Analytic solutions
------------------

With this simple boundary condition, we can find a few analytic
solutions as a baseline.

On noting that each boundary condition is zero when the respective
coordinate is zero or L, we find one solution by adding the two
boundary conditions together:

.. math:: \sin(6\pi x/L) + \sin(6\pi y/L)

Since the sine is two-pi periodic and the boundary conditions reverse
sign at two of the corners, we may also find a solution by summing the
coordinates before taking the sine:

.. math:: \sin(6\pi (x+y)/L)

Another option is the solution to the Laplace equation:

.. math:: (\sin(6\pi x/L) \cosh(6\pi (y-L/2)/L) + \sin(6\pi y/L)\cosh(6\pi (x-L/2)/L))/\cosh(3\pi)

.. plot::
   :alt: "A visual comparison of the functions given above"
   :context: close-figs

   First, a bit of setup.  I choose :math:`L = 600` because it has
   many divisors.

   >>> import numpy as np
   >>> L = 600
   >>> grid_x, grid_y = np.mgrid[:L + 1, :L + 1]

   Next, we define a boolean mask to check values later, and define
   the function values on the boundary.

   >>> mask = np.ones_like(grid_x, dtype=bool)
   >>> mask[1:-1, 1:-1] = False
   >>> edge_x = grid_x[mask]
   >>> edge_y = grid_y[mask]
   >>> border_values = np.sin(6 * np.pi * edge_x / L) + np.sin(6 * np.pi * edge_y / L)

   Next the analytic solutions and plots, with a function to ensure
   the different functions are plotted the same way:

   >>> import matplotlib.pyplot as plt
   >>> def plot_and_check_filled_image(values, name):
   ...     plt.imshow(values, vmin=-2, vmax=2, cmap="RdBu_r")
   ...     plt.title(name)
   ...     error = values[mask] - border_values
   ...     print(f"{name:14s}: Mean error: {np.mean(error):.3g} Error std: {np.std(error):.3g}")
   ...
   >>> sum_of_sines = np.sin(6 * np.pi * grid_x / L) + np.sin(6 * np.pi * grid_y / L)
   >>> plt.subplot(131)
   >>> plot_and_check_filled_image(sum_of_sines, "Sum of sines")
   Sum of sines  : Mean error: 0 Error std: 0
   >>> sine_of_sum = np.sin(6 * np.pi * (grid_x + grid_y) / L)
   >>> plt.subplot(132)
   >>> plot_and_check_filled_image(sine_of_sum, "Sine of sum")
   Sine of sum   : Mean error: 4e-16 Error std: 1.28e-15
   >>> laplace_solution = (
   ...     np.sin(6 * np.pi * grid_x / L) * np.cosh(6 * np.pi * (grid_y - L / 2) / L)
   ...     + np.sin(6 * np.pi * grid_y / L) * np.cosh(6 * np.pi * (grid_x - L / 2) / L)
   ... ) / np.cosh(3 * np.pi)
   >>> plt.subplot(133)
   >>> plot_and_check_filled_image(laplace_solution, "Laplace solution")
   Laplace solution: Mean error: 3.4e-16 Error std: 3.61e-16

Not every boundary so readily admits analytic solutions, so we broaden
our search to other options.

---------------------------------
Simple unstructured interpolators
---------------------------------

One option for filling the border is to interpolate between the border
points.  The simplest SciPy interpolation routines are the zero-order
nearest-neighbor interpolator, the first-order linear interpolator,
and the third-order Clough-Tocher interpolator.  The latter two don't
handle the situation with dense borders and empty interior terribly
well and produce NaNs in the output.  Using only every other point
along the border reduces but does not eliminate the problem, but using
only every third point allows all interpolators to provide a value for
every point.

.. plot::
   :alt: "Use SciPy interpolators to fill from the boundary"
   :context: close-figs

   >>> from scipy.interpolate import (
   ...     NearestNDInterpolator, LinearNDInterpolator, CloughTocher2DInterpolator
   ... )
   >>> for i in range(3):
   ...     sparse_mask = np.zeros_like(grid_x, dtype=bool)
   ...     sparse_mask[::i+1, ::i+1] = True
   ...     sparse_mask[1:-1, 1:-1] = False
   ...     sparse_x = grid_x[sparse_mask]
   ...     sparse_y = grid_y[sparse_mask]
   ...     sparse_values = sum_of_sines[sparse_mask]
   ...     for j, interpolator_class in enumerate(
   ...         [NearestNDInterpolator, LinearNDInterpolator, CloughTocher2DInterpolator]
   ...     ):
   ...         interpolator = interpolator_class((sparse_x, sparse_y), sparse_values)
   ...         values = interpolator(grid_x, grid_y)
   ...         plt.subplot(3, 3, 3 * i + j + 1)
   ...         interpolator_name = interpolator_class.__name__.split("D")[0][:-1]
   ...         plot_and_check_filled_image(values, f"{interpolator_name:s}\n1 pt. in {i+1:d}")
   ...         print("Number of NaNs:", np.count_nonzero(np.isnan(values)), end="\n\n")
   ...
    Nearest
    1 pt. in 1: Mean error: 0 Error std: 0
    Number of NaNs: 0
    <BLANKLINE>
    Linear
    1 pt. in 1: Mean error: nan Error std: nan
    Number of NaNs: 104
    <BLANKLINE>
    CloughTocher
    1 pt. in 1: Mean error: nan Error std: nan
    Number of NaNs: 104
    <BLANKLINE>
    Nearest
    1 pt. in 2: Mean error: 1.85e-05 Error std: 0.0157
    Number of NaNs: 0
    <BLANKLINE>
    Linear
    1 pt. in 2: Mean error: -1.51e-17 Error std: 0.000247
    Number of NaNs: 2
    <BLANKLINE>
    CloughTocher
    1 pt. in 2: Mean error: -2.79e-10 Error std: 2.18e-06
    Number of NaNs: 2
    <BLANKLINE>
    Nearest
    1 pt. in 3: Mean error: -5.66e-18 Error std: 0.0181
    Number of NaNs: 0
    <BLANKLINE>
    Linear
    1 pt. in 3: Mean error: -3.37e-17 Error std: 0.00057
    Number of NaNs: 0
    <BLANKLINE>
    CloughTocher
    1 pt. in 3: Mean error: -2.98e-09 Error std: 6.26e-06
    Number of NaNs: 0

In other words, the linear interpolator has problems filling
rectangles with more than two hundred points on a side.

Simplex tolerance
-----------------

The difficulties had with the linear and cubic interpolators arise
from a floating-point tolerance value: increasing this tolerance will
also reduce the problem:

.. plot::
   :context: close-figs

   >>> for i, interpolator_class in enumerate(
   ...     [LinearNDInterpolator, CloughTocher2DInterpolator]
   ... ):
   ...     interpolator = interpolator_class((edge_x, edge_y), border_values)
   ...     multiplier = 3
   ...     values = interpolator(grid_x, grid_y, simplex_tolerance=multiplier)
   ...     plt.subplot(1, 2, i+1)
   ...     interpolator_name = interpolator_class.__name__.split("D")[0][:-1]
   ...     plot_and_check_filled_image(values, f"{interpolator_name:s}\nmul = {multiplier}")
   ...     print("Number of NaNs:", np.count_nonzero(np.isnan(values)), end="\n\n")
   ...

----------------------
Radial basis functions
----------------------

The radial-basis-function interpolators use all input points to
determine interpolated values, weighting nearer points more and
farther points less.  The kernel determines exactly how much more or
less weight a point will have based on distance: we will use the
linear, thin-plate spline, cubic, and quintic kernels here because
they don't require additional configuration for interpolation.

Since this weighting must be done for each target point, this
interpolation can get a bit slow, so we thin the target grid by a
factor of 8 (noting that 8 divides L, so that all borders will be
included in the output) in each direction, then use a regular-grid
interpolator to get back to the desired density.

.. plot::
   :context: close-figs

   >>> from scipy.interpolate import RBFInterpolator, RegularGridInterpolator
   >>> for i, kernel in enumerate(["linear", "thin_plate_spline", "cubic", "quintic"]):
   ...     rbf_interpolator = RBFInterpolator(
   ...         np.column_stack([edge_x, edge_y]), border_values, kernel=kernel
   ...     )
   ...     sparse_factor = 8
   ...     sparse_x = grid_x[::sparse_factor, ::sparse_factor].reshape(-1)
   ...     sparse_y = grid_y[::sparse_factor, ::sparse_factor].reshape(-1)
   ...     values_vector = rbf_interpolator(np.column_stack([sparse_x, sparse_y]))
   ...     sparse_side = L // sparse_factor + 1
   ...     sparse_values = values_vector.reshape(sparse_side, sparse_side)
   ...     refiner = RegularGridInterpolator(
   ...         (grid_x[::sparse_factor, 0], grid_y[0, ::sparse_factor]),
   ...         sparse_values
   ...     )
   ...     values = refiner(np.stack([grid_x, grid_y], -1))
   ...     plot_and_check_filled_image(values, kernel)
   ...
   linear        : Mean error: -1.57e-14 Error std: 0.00407
   thin_plate_spline: Mean error: -2.57e-13 Error std: 0.00407
   cubic         : Mean error: -2.86e-13 Error std: 0.00407
   quintic       : Mean error: -3.92e-08 Error std: 0.00407

-----------------------
Multi-linear regression
-----------------------

We are trying to predict the values of a function based on observed
values, which is what linear regression is trying to do.  In
particular, it is trying to find a linear combination of various
functions that fits the desired function, then uses that same linear
combination of the same provided functions at new points to predict
the desired function there.

If we collect the values of the desired function into a column vector
with one point per line and call that :math:`y`, and similarly collect
the provided functions into column vectors and collect those columns
into a matrix :math:`X`, we seek a vector :math:`\beta` such that
:math:`y = X \beta` is as close to satisfied as possible.  Statsmodels
gives the least-squares value of :math:`\beta` as
:math:`\beta = (X'X)^{-1} X' y`.

For the columns of :math:`X`, we start with polynomials from the NumPy
polynomials module: specifically the Legendre polynomials, which are
orthogonal on the interval :math:`[-1, 1]`, which should help the
numerics of the inverse, and the Chebyshev polynomials which are
frequently used for interpolation.

A quadratic polynomial can produce one extrema, a cubic polynomial
two, and a quartic three, so a seventh-degree polynomial should be
able to produce the six extrema of the boundary conditions.

.. plot::
   :context: close-figs

   >>> from numpy.polynomial.legendre import legval2d, leggrid2d
   >>> from numpy.polynomial.chebyshev import chebval2d, chebgrid2d
   >>> import scipy.linalg as la

   It will be helpful to map the x and y coordinates to the interval
   :math:`[0, 1]`.

   >>> scaled_grid_x = grid_x * 2 / L - 1
   >>> scaled_grid_y = grid_y * 2 / L - 1
   >>> scaled_edge_x = edge_x * 2 / L - 1
   >>> scaled_edge_y = edge_y * 2 / L - 1

   A seventh-degree polynomial will require eight coefficients,
   including the zero-order term.  We try eleventh- and
   fifteenth-degree polynomials as well, to see if that improves the
   fit.

   Since we provide function values only along the border of the
   region, the function values may be similar: we use the
   pseudo-inverse to avoid problems from the duplication.

   >>> for i, num_coeffs in enumerate([8, 12, 16]):
   ...     for j, (polyval2d, polygrid2d) in enumerate(
   ...         [(legval2d, leggrid2d), (chebval2d, chebgrid2d)]
   ...     ):
   ...         poly_family = polyval2d.__module__.rsplit(".", 1)[-1]
   ...         poly_edge_values = polyval2d(
   ...             scaled_edge_x, scaled_edge_y,
   ...             np.eye(num_coeffs**2).reshape(num_coeffs, num_coeffs, num_coeffs**2)
   ...         )
   ...         coeffs = la.pinvh(
   ...             poly_edge_values @ poly_edge_values.T,
   ...             rtol=np.finfo(poly_edge_values.dtype).eps * L * num_coeffs
   ...         ) @ (poly_edge_values @ border_values)
   ...         values = polygrid2d(
   ...             scaled_grid_x[:, 0], scaled_grid_y[0, :], coeffs.reshape(num_coeffs, num_coeffs)
   ...         )
   ...         plot_and_check_filled_image(values, f"{poly_family:s} degree {num_coeffs - 1:d}")
   legendre degree 7: Mean error: 1.21e-16 Error std: 0.408
   chebyshev degree 7: Mean error: -2.66e-17 Error std: 0.408
   legendre degree 11: Mean error: -9.44e-18 Error std: 0.0241
   chebyshev degree 11: Mean error: -6.65e-17 Error std: 0.0241
   legendre degree 15: Mean error: -1.01e-16 Error std: 0.000299
   chebyshev degree 15: Mean error: -5.93e-17 Error std: 0.000299

The higher-order polynomials produce interesting shapes in the
interior of the domain, which are pretty, but slightly concerning.  We
can further restrict the set of polynomials to those which introduce
no new extrema in the interior of the domain, by insisting all
polynomials must have one degree less than two.  We notice that this
description applies to the "sum of sines" described and plotted above,
as it is the sum of one function constant in x and one function
constant in y.  It is therefore likely that these regressions will
produce results similar to that function.


.. plot::
   :context: close-figs

   >>> for i, num_coeffs in enumerate([8, 12, 16]):
   ...     for j, (polyval2d, polygrid2d) in enumerate(
   ...         [(legval2d, leggrid2d), (chebval2d, chebgrid2d)]
   ...     ):
   ...         poly_family = polyval2d.__module__.rsplit(".", 1)[-1]
   ...         coeff_multiplier = np.ones((num_coeffs, num_coeffs))
   ...         coeff_multiplier[2:, 2:] = 0
   ...         poly_edge_values = polyval2d(
   ...             scaled_edge_x, scaled_edge_y,
   ...             np.eye(num_coeffs**2).reshape(num_coeffs, num_coeffs, num_coeffs**2)
   ...             * coeff_multiplier[:, :, np.newaxis]
   ...         )
   ...         coeffs = la.pinvh(
   ...             poly_edge_values @ poly_edge_values.T,
   ...             rtol=np.finfo(poly_edge_values.dtype).eps * L * num_coeffs
   ...         ) @ (poly_edge_values @ border_values)
   ...         values = polygrid2d(
   ...             scaled_grid_x[:, 0], scaled_grid_y[0, :],
   ...             coeffs.reshape(num_coeffs, num_coeffs) * coeff_multiplier[:, :]
   ...         )
   ...         plot_and_check_filled_image(values, f"{poly_family:s} degree {num_coeffs - 1:d}")
   legendre degree 7: Mean error: -3.23e-16 Error std: 0.408
   chebyshev degree 7: Mean error: 6.81e-17 Error std: 0.408
   legendre degree 11: Mean error: 3.42e-16 Error std: 0.0241
   chebyshev degree 11: Mean error: -2.64e-16 Error std: 0.0241
   legendre degree 15: Mean error: 1.41e-16 Error std: 0.000299
   chebyshev degree 15: Mean error: 7.43e-17 Error std: 0.000299


A different way to restrict the polynomials to avoid interior extrema
is to insist that the individual polynomials each satisfy the Laplace
equation, as do, for example,
:math:`1,x,y,x^2-y^2,xy,x^3-3xy^2,y^3-3yx^2,x^4-6x^2y^2+y^4,x^3y-xy^3,\dots`.

.. plot::
   :context: close-figs

   >>> def laplace_polynomials(x, y):
   ...     return np.stack([
   ...         np.ones_like(x),
   ...         x,
   ...         y,
   ...         x**2 - y**2,
   ...         x * y,
   ...         x**3 - 3 * x * y**2,
   ...         y**3 - 3 * y * x**2,
   ...         x**4 - 6 * x**2 * y**2 + y**4,
   ...         x**3 * y - x * y**3,
   ...         x**5 - 10 * x**3 * y**2 + 5 * x * y**4,
   ...         y**5 - 10 * y**3 * x**2 + 5 * y * x**4,
   ...         x**6 - 15 * x**4 * y**2 + 15 * x**2 * y**4 - y**6,
   ...         3 * x**5 * y - 10 * x**3 * y**3 + 3 * x * y**5,
   ...         x**7 - 21 * x**5 * y**2 + 35 * x**3 * y**4 - 7 * x * y**6,
   ...         y**7 - 21 * y**5 * x**2 + 35 * y**3 * x**4 - 7 * y * x**6
   ...     ])
   ...
   >>> poly_edge_values = laplace_polynomials(scaled_edge_x, scaled_edge_y)
   >>> coeffs = la.pinvh(
   ...     poly_edge_values @ poly_edge_values.T,
   ...     rtol=np.finfo(poly_edge_values.dtype).eps * L * num_coeffs
   ... ) @ (poly_edge_values @ border_values)
   >>> values = np.tensordot(
   ...     coeffs, laplace_polynomials(scaled_grid_x, scaled_grid_y), 1
   ... )
   >>> plot_and_check_filled_image(values, "Laplace polynomial")
   Laplace polynomial: Mean error: -3.55e-17 Error std: 0.663

--------------
Other packages
--------------

There are packages that will solve Laplace equation with given boundary conditions.

Scikit-Image has an ``inpaint_biharmonic`` function that will solve
also fill an image from the boundary conditions, in a manner related
to the Laplace Equation.


Statsmodels offers a regression interface that hides the linear
algebra and offers additional diagnostics.
