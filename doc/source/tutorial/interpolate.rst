.. _tutorial-interpolate:

========================================
Interpolation (:mod:`scipy.interpolate`)
========================================

.. sectionauthor:: Travis E. Oliphant

.. sectionauthor:: Pauli Virtanen

.. sectionauthor:: Evgeni Burovski

.. currentmodule:: scipy.interpolate


There are several general facilities available in SciPy for interpolation and
smoothing for data in 1, 2, and higher dimensions. The choice of a specific
interpolation routine depends on the data: whether it is one-dimensional,
is given on a structured grid, or is unstructured. One other factor is the
desired smoothness of the interpolator. In short, routines recommended *for
interpolation* can be summarized as follows:

+------------------+-------------------------+------------------------------+------------------------+---------------------------------------+
|                  | **kind**                | **routine**                  | **continuity**         | **comment**                           |
+==================+=========================+==============================+========================+=======================================+
|                  | linear                  | `numpy.interp`               | piecewise continuous   | Alternatively,                        |
|                  |                         |                              |                        | ``make_interp_spline(..., k=1)``      |
+                  +-------------------------+------------------------------+------------------------+---------------------------------------+
|                  | cubic spline            | `CubicSpline`                | 2nd derivative         |                                       |
+                  +-------------------------+------------------------------+------------------------+---------------------------------------+
|       1D         | monotone cubic spline   | `PchipInterpolator`          | 1st derivative         | non-overshooting                      |
+                  +-------------------------+------------------------------+------------------------+---------------------------------------+
|                  | non-cubic spline        | `make_interp_spline`         | (k-1)th derivative     | ``k=3`` is equivalent to `CubicSpline`|
+                  +-------------------------+------------------------------+------------------------+---------------------------------------+
|                  | nearest                 | `interp1d`                   |                        | kind='nearest', 'previous', 'next'    |
+------------------+-------------------------+------------------------------+------------------------+---------------------------------------+
| N-D curve        | nearest, linear, spline | `make_interp_spline`         | (k-1)th derivative     | use N-dim `y` array                   |
+------------------+-------------------------+------------------------------+------------------------+---------------------------------------+
|                  | nearest                 |                              |                        |  method='nearest'                     |
+                  +-------------------------+                              +                        +---------------------------------------+
| N-D regular      | linear                  | `RegularGridInterpolator`    |                        |  method='linear'                      |
+ (*rectilinear*)  +-------------------------+                              +------------------------+---------------------------------------+
| grid             | splines                 |                              | 2nd derivatives        |  method='cubic', 'quintic'            |
+                  +-------------------------+                              +------------------------+---------------------------------------+
|                  | monotone splines        |                              | 1st derivatives        |  method='pchip'                       |
+------------------+-------------------------+------------------------------+------------------------+---------------------------------------+
|                  | nearest                 | `NearestNDInterpolator`      |                        |                                       |
+                  +-------------------------+------------------------------+                        +                                       +
|  N-D scattered   | linear                  | `LinearNDInterpolator`       |                        | alias: `griddata`                     |
+                  +-------------------------+------------------------------+------------------------+                                       +
|                  | cubic (2D only)         | `CloughTocher2DInterpolator` | 1st derivatives        |                                       |
+                  +-------------------------+------------------------------+------------------------+---------------------------------------+
|                  | radial basis function   | `RBFInterpolator`            |                        |                                       |
+------------------+-------------------------+------------------------------+------------------------+---------------------------------------+


Smoothing and approximation of data
===================================

+-------------------------------+-------------------------+-------------------------------+
|                               | `make_smoothing_spline` |  classic smoothing splines,   |
|                               |                         |  GCV penalty                  |
| 1D spline functions           +-------------------------+-------------------------------+
|                               | `make_splrep`           | automated/semi-automated knot |
|                               |                         | selection                     |
+-------------------------------+-------------------------+-------------------------------+
| spine curves in N-D           | `make_splprep`          |                               |
+-------------------------------+-------------------------+-------------------------------+
| unconstrained least squares   | `make_lsq_spline`       |                               |
| spline fit                    |                         |                               |
+-------------------------------+-------------------------+-------------------------------+
| 2D smoothing surfaces         | `bisplrep`              | scattered data                |
+-------------------------------+-------------------------+-------------------------------+
|                               | `RectBivariateSpline`   | gridded data                  |
+-------------------------------+-------------------------+-------------------------------+
| Radial basis functions in N-D | `RBFInterpolator`       |                               |
+-------------------------------+-------------------------+-------------------------------+


**Further details are given in the links below**


.. toctree::
   :maxdepth: 3

   interpolate/1D
   interpolate/splines_and_polynomials
   interpolate/smoothing_splines
   interpolate/ND_regular_grid
   interpolate/ND_unstructured
   interpolate/extrapolation_examples
   interpolate/interp_transition_guide.md
