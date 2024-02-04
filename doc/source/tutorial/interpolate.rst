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
|                  | linear                  | `numpy.interp`               | piecewise continuous   | comes from numpy                      |
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


For data smoothing, :ref:`functions are provided <tutorial-interpolate_fitpack>`
for 1- and 2-D data using cubic splines, based on the FORTRAN library FITPACK. 

Additionally, routines are provided for interpolation / smoothing using
:ref:`radial basis functions <tutorial-interpolate_RBF>` with several kernels.

Further details are given in the links below. 

.. toctree::
   :maxdepth: 3

   interpolate/1D
   interpolate/splines_and_polynomials
   interpolate/smoothing_splines
   interpolate/ND_regular_grid
   interpolate/ND_unstructured
   interpolate/extrapolation_examples

