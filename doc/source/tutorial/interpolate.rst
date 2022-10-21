========================================
Interpolation (:mod:`scipy.interpolate`)
========================================

.. sectionauthor:: Travis E. Oliphant

.. sectionauthor:: Pauli Virtanen

.. sectionauthor:: Evgeni Burovski

.. currentmodule:: scipy.interpolate


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

+------------------+-------------------------+------------------------------+------------------------+---------------------------------------+
|                  | **kind**                | **routine**                   | **continuity**        | **comment**                           |
+==================+=========================+==============================+========================+=======================================+
|                  | linear                  | `numpy.interp`               | broken line            | comes from numpy                      |
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
+------------------+-------------------------+------------------------------+------------------------+---------------------------------------+




.. toctree::
   :maxdepth: 3

   interpolate/1D
   interpolate/splines_and_polynomials
   interpolate/smoothing_splines
   interpolate/ND_regular_grid
   interpolate/ND_unstructured


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







