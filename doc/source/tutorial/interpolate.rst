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

For data smoothing, functions are provided for 1- and 2-D (smoothed)
cubic-spline interpolation, based on the FORTRAN library FITPACK. There are
two interfaces for the FITPACK library, a procedural interface and an
object-oriented interface.  **TODO**


Missing data
============

Before describing various interpolators in detail, we note the (lack of) support
for interpolation with missing data. Two popular ways of representing missing
data are using masked arrays of the `numpy.ma` library, and encoding missing
values as not-a-number, ``NaN``. 

None of these two approaches is directly suppored in `scipy.interpolate`.
Individual routines may offer partial support, and/or workarounds, but in
general the library firmly adheres to the IEEE 754 semantics where a ``NaN``
means *not-a-number*, i.e. a result of an illegal mathematical operation
(think a division by zero), not *missing*.


1-D interpolation
=================

.. _tutorial-interpolate_1Dsection:

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

Of course, piecewise linear interpolation produces cusps at data points,
where linear pieces join. To produce a smoother curve, you can use cubic
splines, where the interpolating curve is made of cubic pieces with matching
first and second derivatives. In code, these objects are represented via the
`CubicSpline`` class instances. An instance is constructed with the ``x`` and
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
a data with an outlier:

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


Piecewise polynomials
---------------------

Internally, `CubicSpline` and monotone interpolants are represented as instances
of a `PPoly` class, which represents pieciwise polynomials in terms of
breakpoints and coefficients (`PPoly` objects can represent polynomials of 
arbitrary orders, not only cubics). For the data array ``x``, breakpoints are at
the data points, and the array of coefficients, ``c`` , define cubic polynomials
such that ``c[k, j]`` is a coefficient for ``(x - x[j])**(3-k)`` on the segment
between ``x[j]`` and ``x[j+1]`` .


Manipulating ``PPoly`` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _tutorial-interpolate_ppoly:

`PPoly` objects have convenient methods for constructing derivatives
and antiderivatives, computing integrals and root-finding. For example, we
tabulate the sine function and find the roots of its derivative.

    >>> from scipy.interpolate import CubicSpline
    >>> x = np.linspace(0, 10, 71)
    >>> y = np.sin(x)
    >>> spl = CubicSpline(x, y)

Now, differentiate the spline:

    >>> dspl = spl.derivative()

Here ``dspl`` is a `PPoly` instance which represents a polynomial approximation
to the derivative of the original object, ``spl`` . Evaluating ``dspl`` at a
fixed argument is equivalent to evaluating the original spline with the ``nu=1``
argument:

    >>> dspl(1.1), spl(1.1, nu=1)
    (0.45361436, 0.45361436)

However with the ``dspl`` object, we find the zeros of the derivative of ``spl``:

    >>> dspl.root() / np.pi
    array([-0.45480801,  0.50000034,  1.50000099,  2.5000016 ,  3.46249993])

This agrees well with roots :math:`\pi/2 + \pi\,n` of
:math:`\cos(x) = \sin'(x)`.
Note that by default it computed the roots *extrapolated* to the outside of
the interpolation interval :math:`0 \leqslant x \leqslant 10`, and that
the extrapolated results (the first and last values) are much less accurate.
We can switch off the extrapolation and limit the root-finding to the
interpolation interval:

    >>> dspl.root(extrapolate=False) / np.pi
    array([0.50000034,  1.50000099,  2.5000016])

In fact, the ``root`` method is a special case of a more general ``solve``
method which finds for a given value of :math:`y` the solutions of the
equation :math:`f(x) = y` , where :math:`f(x)` is a piecewise polynomial:

    >>> dspl.solve(0.5, extrapolate=False)
    array([0.33332755, 1.66667195, 2.3333271])

which agrees well with the expected values of  :math:`\pm\arccos(1/2) + 2\pi\,n`.

Integrals of piecewise polynomials can be computed using the ``.integrate``
method which accepts the lower and the upper limits of integration. As an
example, we compute an approximation to the complete elliptic integral
:math:`K(m) = \int_0^{\pi/2} [1 - m\sin^2 x]^{-1/2} dx`:

    >>> from scipy.special import ellipk
    >>> ellipk(0.5)
    1.8540746773013719

To this end, we tabulate the integrand, interpolate using the monotone
interpolant (we could as well used a `CubicSpline`):

    >>> from scipy.interpolate import PchipInterpolator
    >>> x = np.linspace(0, np.pi/2, 70)
    >>> y = 1 / np.sqrt(1 - 0.5*np.sin(x)**2)
    >>> spl = PchipInterpolator(x, y)

and integrate

    >>> spl.integrate(0, np.pi/2)
    1.854074674965991

which is indeed close to the value computed by `scipy.special.ellipk`.

All piecewise polynomials can be constructed with N-dimensional ``y`` values.
If ``y.ndim > 1``, it is understood as a stack of 1D ``y`` values, which are
arranged along the interpolation axis (with the default value if 0).
The latter is specified via the ``axis`` argument, and the invariant is that
``len(x) == y.shape[axis]``. As an example, we extend the elliptic integral
example above to compute the approximation for a range of ``m`` values, using
the NumPy broadcasting:

.. plot::

    >>> from scipy.interpolate import PchipInterpolator
    >>> m = np.linspace(0, 0.9, 11)
    >>> x = np.linspace(0, np.pi/2, 70)
    >>> y = 1 / np.sqrt(1 - m[:, None]*np.sin(x)**2)

    Now the ``y`` array has the shape ``(11, 70)``, so that the values of ``y``
    for fixed value of ``m`` are along the second axis of the ``y`` array.

    >>> spl = PchipInterpolator(x, y, axis=1)  # the default is axis=0
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(m, spl.integrate(0, np.pi/2), '--')

    >>> from scipy.special import ellipk
    >>> plt.plot(m, ellipk(m), 'o')
    >>> plt.show()


Interpolation with B-splines
----------------------------

A polynomial of degree :math:`k` can be thought of as a linear combination of
:math:`k+1` monomial basis elements, :math:`1, x, x^2, \cdots, x^k`. B-splines
form an alternative (if equivalent) :ref:`basis <tutorial-interpolate_bspl_basis>`
of degree-:math:`k` piecewise polynomials. 

As an illustration, we construct
the interpolation of a :ref:`derivative of a sine function<tutorial-interpolate_ppoly>`
with b-splines:

.. plot::

    >>> x = np.linspace(0, 3/2, 7)
    >>> y = np.sin(np.pi*x)

    To construct the interpolating objects given data arrays, ``x`` and ``y``,
    we use the `make_interp_spline` function:

    >>> from scipy.interpolate import make_interp_spline
    >>> bspl = make_interp_spline(x, y)

    Now ``bspl`` is a `BSpline` object which has an interface similar to `PPoly`.
    In particular, it can be evaluated at a data point and differentiated:

    >>> der = bspl.derivative()      # a BSpline representing the derivative
    >>> import matplotlib.pyplot as plt
    >>> xx = np.linspace(0, 3/2, 51)
    >>> plt.plot(xx, bspl(xx), '--', label=r'$\sin(\pi x)$ approx')
    >>> plt.plot(x, y, 'o', label='data')
    >>> plt.plot(xx, der(xx)/np.pi, '--', label='$d \sin(\pi x)/dx / \pi$ approx')
    >>> plt.legend()
    >>> plt.show()

Note that by default `make_interp_spline` constructs a cubic spline:

    >>> der.k
    3

This way, the default result of ``make_interp_spline(x, y)`` is equivalent to
``CubicSpline(x, y)``. The difference is that the former allows several optional
capabilities: it can construct splines of various degrees (via the optional
argument ``k``) and predefined knots (via the optional argument ``t``). 

A b-spline of degree ``k`` is defined by its knots and coefficients. The knots
are availble as the ``t`` attribute of a `BSpline` instance:

    >>> print(bspl.t)
    [0.  0.  0.  0.        0.5  0.75  1.        1.5  1.5  1.5  1.5 ]
    >>> print(x)
    [            0.  0.25  0.5  0.75  1.  1.25  1.5 ]

We see that the knot vector by default is constructed from the input
array ``x``: first, it is made :math:`(k+1)` -regular (it has ``k+1``
repeated knots appended and prepended); then, the second and
second-to-last points of the input array are removed---this is the so-called
*not-a-knot* boundary condition. 

In general, an interpolating spline of degree ``k`` needs
``len(t) - len(x) - k - 1`` boundary conditions. For cubic splines with
``(k+1)``-regular knot arrays this means two boundary conditions---or
removing two values from the ``x`` array. Various boundary conditions can be
requested using the optional ``bc_type`` argument of `make_interp_spline`.

The b-spline coefficients are accessed via the ``c`` attribute of a `BSpline`
object:

    >>> len(bspl.c)
    7

The convention is that for ``len(t)`` knots there are ``len(t) - k - 1``
coefficients. Note that the coeffients are given in the
:ref:`b-spline basis <tutorial-interpolate_bspl_basis>`, not the power basis
of :math:`1, x, \cdots, x^k`.


B-spline basis elements
-----------------------

.. _tutorial-interpolate_bspl_basis:

B-splines are piecewise polynomials, represented as linear combinations of
*b-spline basis elements* --- which themselves are certain linear combinations
of usual monomials, :math:`x^m` with :math:`m=0, 1, \dots, k`.

The b-spline basis generally more computationally stable then the power basis
and is useful for variety of applications which include interpolation, regression
and curve representation. The main feature is that these basis elements are
*localized* and equal zero outside of an interval defined by the *knot array*.

Specifically, a b-spline basis element of degree ``k`` (e.g. ``k=3`` for cubics)
is defined by :math:`k+2` knots and is zero outside of these knots.
To illustrate, plot a collection of non-zero basis elements on a certain
interval:

.. plot ::

    >>> k = 3      # cubic splines
    >>> t = [0., 1.4, 2., 3.1, 5.]   # internal knots
    >>> t = np.r_[[0]*k, t, [5]*k]   # add boundary knots

    >>> from scipy.interpolate import BSpline
    >>> import matplotlib.pyplot as plt
    >>> for j in [-2, -1, 0, 1, 2]:
    ...     a, b = t[k+j], t[-k+j-1]
    ...     xx = np.linspace(a, b, 101)
    ...     bspl = BSpline.basis_element(t[k+j:-k+j])
    ...     plt.plot(xx, bspl(xx), label=f'j = {j}')
    >>> plt.legend(loc='best')
    >>> plt.show()

Here `BSpline.basis_element` is essentially a shorthand for constructing a spline
with only a single non-zero coefficient. For instance, the ``j=2`` element in
the above example is equivalent to

    >>> c = np.zeros(t.size - k - 1)
    >>> c[-2] = 1
    >>> b = BSpline(t, c, k)
    >>> np.allclose(b(xx), bspl(xx))
    True

If desired, a b-spline can be converted into a `PPoly` object using
`PPoly.from_spline` method which accepts a `BSpline` instance and returns a
`PPoly` instance. The reverse conversion is performed by the
`BSpline.from_power_basis` method. However, conversions between bases is best
avoided because it accumulates rounding errors.


Parametric spline curves
------------------------

So far we considered spline *functions*, where the data, ``y``, is expected to
depend explicitly on the independent variable ``x``---so that the interpolating
function satisfies :math:`f(x_j) = y_j`. Spline *curves* treat
the ``x`` and ``y`` arrays as coordinates of points, :math:`\mathbf{p}_j` on a
plane, and an interpolating curve which passes through these points is
parameterized by some additional parameter (typically called ``u``). Note that
this construction readily generalizes to higher dimensions where
:math:`\mathbf{p}_j` are points in an N-dimensional space.

The choice of parametrization is problem-dependent and different parametrizations
may produce vastly different curves. As an example, we consider three
parametrizations of (a somewhat difficult) dataset, which we take from
Chapter 6 of Ref [1] listed in the `BSpline` docstring:

.. plot ::

    >>> x = [0, 1, 2, 3, 4, 5, 6]
    >>> y = [0, 0, 0, 9, 0, 0, 0]
    >>> p = np.stack((x, y))
    >>> p
    array([[0, 1, 2, 3, 4, 5, 6],
           [0, 0, 0, 9, 0, 0, 0]])

    We take elements of the ``p`` array as coordinates of seven point on the
    plane, where ``p[:, j]`` gives the coordinates of the point
    :math:`\mathbf{p}_j`.

    First, consider the *uniform* parametrization, :math:`u_j = j`:

    >>> u_unif = x

    Second, we consider the so-called *cord length* parametrization, which is
    nothing but a cumulative length of straight line segments connecting the
    data points:

    .. math::
        
        u_j = u_{j-1} + |\mathbf{p}_j - \mathbf{p}_{j-1}|

    for :math:`j=1, 2, \dots` and :math:`u_0 = 0`. Here :math:`| \cdots |` is the
    length between the consecutive points :math:`p_j` on the plane.

    >>> dp = p[:, 1:] - p[:, :-1]      # 2-vector distances between points
    >>> l = (dp**2).sum(axis=0)        # squares of lengths of 2-vectors between points
    >>> u_cord = np.sqrt(l).cumsum()   # cumulative sums of 2-norms
    >>> u_cord = np.r_[0, u_cord]      # first point is parametrized at zero

    Finally, we consider what is sometimes called the *centripetal*
    parametrization: :math:`u_j = u_{j-1} + |\mathbf{p}_j - \mathbf{p}_{j-1}|^{1/2}`.
    Due to the extra square root, the difference between consecutive values
    :math:`u_j - u_{j-1}` will be smaller than for the cord length parametrization: 

    >>> u_c = np.r_[0, np.cumsum((dp**2).sum(axis=0)**0.25)]

    Now plot the resulting curves:

    >>> from scipy.interpolate import make_interp_spline
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots(1, 3, figsize=(8, 3))
    >>> parametrizations = ['uniform', 'cord length', 'centripetal']
    >>>
    >>> for j, u in enumerate([u_unif, u_cord, u_c]):
    ...    spl = make_interp_spline(u, p, axis=1)    # note p is a 2D array
    ...    
    ...    uu = np.linspace(u[0], u[-1], 51)
    ...    xx, yy = spl(uu)
    ...    
    ...    ax[j].plot(xx, yy, '--')
    ...    ax[j].plot(p[0, :], p[1, :], 'o')
    ...    ax[j].set_title(parametrizations[j])
    >>> plt.show()


Legacy interface for 1-D interpolation (:class:`interp1d`)
----------------------------------------------------------

.. _tutorial-interpolate_interp1d:

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



Multivariate data interpolation on a regular grid  (:class:`RegularGridInterpolator`)
======================================================================================

.. _tutorial-interpolate_regular_grid_interpolator:

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
optimal choice. Consider using `ndimage.map_coordinates` instead.



- A class representing an interpolant (:class:`interp1d`) in 1-D,
  offering several interpolation methods.

- Convenience function :func:`griddata` offering a simple interface to
  interpolation in N dimensions (N = 1, 2, 3, 4, ...).
  Object-oriented interface for the underlying routines is also
  available.

- :class:`RegularGridInterpolator` provides several interpolation methods
  on a regular grid in arbitrary (N) dimensions,

- Functions for 1- and 2-D (smoothed) cubic-spline
  interpolation, based on the FORTRAN library FITPACK. They are both
  procedural and object-oriented interfaces for the FITPACK library.

- Interpolation using radial basis functions.



Multivariate data interpolation (:func:`griddata`)
==================================================

.. _tutorial-interpolate_NDunstructured:

Suppose you have multidimensional data, for instance, for an underlying
function *f(x, y)* you only know the values at points *(x[i], y[i])*
that do not form a regular grid.

.. plot::
    :alt: " "

    Suppose we want to interpolate the 2-D function

    >>> def func(x, y):
    ...     return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

    on a grid in [0, 1]x[0, 1]

    >>> grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]

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
    cubic interpolant gives the best results:

    >>> import matplotlib.pyplot as plt
    >>> plt.subplot(221)
    >>> plt.imshow(func(grid_x, grid_y).T, extent=(0,1,0,1), origin='lower')
    >>> plt.plot(points[:,0], points[:,1], 'k.', ms=1)
    >>> plt.title('Original')
    >>> plt.subplot(222)
    >>> plt.imshow(grid_z0.T, extent=(0,1,0,1), origin='lower')
    >>> plt.title('Nearest')
    >>> plt.subplot(223)
    >>> plt.imshow(grid_z1.T, extent=(0,1,0,1), origin='lower')
    >>> plt.title('Linear')
    >>> plt.subplot(224)
    >>> plt.imshow(grid_z2.T, extent=(0,1,0,1), origin='lower')
    >>> plt.title('Cubic')
    >>> plt.gcf().set_size_inches(6, 6)
    >>> plt.show()


.. _tutorial-interpolate_regular_grid_interpolator:

Multivariate data interpolation on a regular grid  (:class:`RegularGridInterpolator`)
======================================================================================

Suppose you have n-dimensional data on a regular grid, and you want to interpolate it.
In such a case, :class:`RegularGridInterpolator` can be useful.
The following example demonstrates its use, and compares the interpolation results
using each method.

.. plot::
   :alt: " "

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


Spline interpolation
====================

.. _tutorial-interpolate_splXXX:

Spline interpolation in 1-D: Procedural (interpolate.splXXX)
------------------------------------------------------------

Spline interpolation requires two essential steps: (1) a spline
representation of the curve is computed, and (2) the spline is
evaluated at the desired points. In order to find the spline
representation, there are two different ways to represent a curve and
obtain (smoothing) spline coefficients: directly and parametrically.
The direct method finds the spline representation of a curve in a 2-D
plane using the function :obj:`splrep`. The
first two arguments are the only ones required, and these provide the
:math:`x` and :math:`y` components of the curve. The normal output is
a 3-tuple, :math:`\left(t,c,k\right)` , containing the knot-points,
:math:`t` , the coefficients :math:`c` and the order :math:`k` of the
spline. The default spline order is cubic, but this can be changed
with the input keyword, *k.*

For curves in N-D space the function
:obj:`splprep` allows defining the curve
parametrically. For this function only 1 input argument is
required. This input is a list of :math:`N`-arrays representing the
curve in N-D space. The length of each array is the
number of curve points, and each array provides one component of the
N-D data point. The parameter variable is given
with the keyword argument, *u,*, which defaults to an equally-spaced
monotonic sequence between :math:`0` and :math:`1` . The default
output consists of two objects: a 3-tuple, :math:`\left(t,c,k\right)`
, containing the spline representation and the parameter variable
:math:`u.`

The keyword argument, *s* , is used to specify the amount of smoothing
to perform during the spline fit. The default value of :math:`s` is
:math:`s=m-\sqrt{2m}` where :math:`m` is the number of data-points
being fit. Therefore, **if no smoothing is desired a value of**
:math:`\mathbf{s}=0` **should be passed to the routines.**

Once the spline representation of the data has been determined,
functions are available for evaluating the spline
(:func:`splev`) and its derivatives
(:func:`splev`, :func:`spalde`) at any point
and the integral of the spline between any two points (
:func:`splint`). In addition, for cubic splines ( :math:`k=3`
) with 8 or more knots, the roots of the spline can be estimated (
:func:`sproot`). These functions are demonstrated in the
example that follows.

.. plot::
   :alt: " "

   >>> import numpy as np
   >>> import matplotlib.pyplot as plt
   >>> from scipy import interpolate

   Cubic-spline

   >>> x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/8)
   >>> y = np.sin(x)
   >>> tck = interpolate.splrep(x, y, s=0)
   >>> xnew = np.arange(0, 2*np.pi, np.pi/50)
   >>> ynew = interpolate.splev(xnew, tck, der=0)

   >>> plt.figure()
   >>> plt.plot(x, y, 'x', xnew, ynew, xnew, np.sin(xnew), x, y, 'b')
   >>> plt.legend(['Linear', 'Cubic Spline', 'True'])
   >>> plt.axis([-0.05, 6.33, -1.05, 1.05])
   >>> plt.title('Cubic-spline interpolation')
   >>> plt.show()

   Derivative of spline

   >>> yder = interpolate.splev(xnew, tck, der=1)
   >>> plt.figure()
   >>> plt.plot(xnew, yder, xnew, np.cos(xnew),'--')
   >>> plt.legend(['Cubic Spline', 'True'])
   >>> plt.axis([-0.05, 6.33, -1.05, 1.05])
   >>> plt.title('Derivative estimation from spline')
   >>> plt.show()

   All derivatives of spline

   >>> yders = interpolate.spalde(xnew, tck)
   >>> plt.figure()
   >>> for i in range(len(yders[0])):
   ...    plt.plot(xnew, [d[i] for d in yders], '--', label=f"{i} derivative")
   >>> plt.legend()
   >>> plt.axis([-0.05, 6.33, -1.05, 1.05])
   >>> plt.title('All derivatives of a B-spline')
   >>> plt.show()

   Integral of spline

   >>> def integ(x, tck, constant=-1):
   ...     x = np.atleast_1d(x)
   ...     out = np.zeros(x.shape, dtype=x.dtype)
   ...     for n in range(len(out)):
   ...         out[n] = interpolate.splint(0, x[n], tck)
   ...     out += constant
   ...     return out

   >>> yint = integ(xnew, tck)
   >>> plt.figure()
   >>> plt.plot(xnew, yint, xnew, -np.cos(xnew), '--')
   >>> plt.legend(['Cubic Spline', 'True'])
   >>> plt.axis([-0.05, 6.33, -1.05, 1.05])
   >>> plt.title('Integral estimation from spline')
   >>> plt.show()

   Roots of spline

   >>> interpolate.sproot(tck)
   array([3.1416])  # may vary

   Notice that `sproot` may fail to find an obvious solution at the edge of the
   approximation interval, :math:`x = 0`. If we define the spline on a slightly
   larger interval, we recover both roots :math:`x = 0` and :math:`x = 2\pi`:

   >>> x = np.linspace(-np.pi/4, 2.*np.pi + np.pi/4, 21)
   >>> y = np.sin(x)
   >>> tck = interpolate.splrep(x, y, s=0)
   >>> interpolate.sproot(tck)
   array([0., 3.1416])

   Parametric spline

   >>> t = np.arange(0, 1.1, .1)
   >>> x = np.sin(2*np.pi*t)
   >>> y = np.cos(2*np.pi*t)
   >>> tck, u = interpolate.splprep([x, y], s=0)
   >>> unew = np.arange(0, 1.01, 0.01)
   >>> out = interpolate.splev(unew, tck)
   >>> plt.figure()
   >>> plt.plot(x, y, 'x', out[0], out[1], np.sin(2*np.pi*unew), np.cos(2*np.pi*unew), x, y, 'b')
   >>> plt.legend(['Linear', 'Cubic Spline', 'True'])
   >>> plt.axis([-1.05, 1.05, -1.05, 1.05])
   >>> plt.title('Spline of parametrically-defined curve')
   >>> plt.show()

Spline interpolation in 1-d: Object-oriented (:class:`UnivariateSpline`)
------------------------------------------------------------------------

The spline-fitting capabilities described above are also available via
an objected-oriented interface.  The 1-D splines are
objects of the `UnivariateSpline` class, and are created with the
:math:`x` and :math:`y` components of the curve provided as arguments
to the constructor.  The class defines :meth:`__call__ <UnivariateSpline.__call__>`, allowing the object
to be called with the x-axis values, at which the spline should be
evaluated, returning the interpolated y-values.  This is shown in
the example below for the subclass `InterpolatedUnivariateSpline`.
The :meth:`integral <UnivariateSpline.integral>`,
:meth:`derivatives <UnivariateSpline.derivatives>`, and
:meth:`roots <UnivariateSpline.roots>` methods are also available
on `UnivariateSpline` objects, allowing definite integrals,
derivatives, and roots to be computed for the spline.

The UnivariateSpline class can also be used to smooth data by
providing a non-zero value of the smoothing parameter `s`, with the
same meaning as the `s` keyword of the :obj:`splrep` function
described above.  This results in a spline that has fewer knots
than the number of data points, and hence is no longer strictly
an interpolating spline, but rather a smoothing spline.  If this
is not desired, the `InterpolatedUnivariateSpline` class is available.
It is a subclass of `UnivariateSpline` that always passes through all
points (equivalent to forcing the smoothing parameter to 0). This
class is demonstrated in the example below.

The `LSQUnivariateSpline` class is the other subclass of `UnivariateSpline`.
It allows the user to specify the number and location of internal
knots explicitly with the parameter `t`.  This allows for the creation
of customized splines with non-linear spacing, to interpolate in
some domains and smooth in others, or change the character of the
spline.


.. plot::
   :alt: " "

   >>> import numpy as np
   >>> import matplotlib.pyplot as plt
   >>> from scipy import interpolate

   InterpolatedUnivariateSpline

   >>> x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/8)
   >>> y = np.sin(x)
   >>> s = interpolate.InterpolatedUnivariateSpline(x, y)
   >>> xnew = np.arange(0, 2*np.pi, np.pi/50)
   >>> ynew = s(xnew)

   >>> plt.figure()
   >>> plt.plot(x, y, 'x', xnew, ynew, xnew, np.sin(xnew), x, y, 'b')
   >>> plt.legend(['Linear', 'InterpolatedUnivariateSpline', 'True'])
   >>> plt.axis([-0.05, 6.33, -1.05, 1.05])
   >>> plt.title('InterpolatedUnivariateSpline')
   >>> plt.show()

   LSQUnivarateSpline with non-uniform knots

   >>> t = [np.pi/2-.1, np.pi/2+.1, 3*np.pi/2-.1, 3*np.pi/2+.1]
   >>> s = interpolate.LSQUnivariateSpline(x, y, t, k=2)
   >>> ynew = s(xnew)

   >>> plt.figure()
   >>> plt.plot(x, y, 'x', xnew, ynew, xnew, np.sin(xnew), x, y, 'b')
   >>> plt.legend(['Linear', 'LSQUnivariateSpline', 'True'])
   >>> plt.axis([-0.05, 6.33, -1.05, 1.05])
   >>> plt.title('Spline with Specified Interior Knots')
   >>> plt.show()

.. _tutorial-interpolate_2d_spline:

2-D spline representation: Procedural (:func:`bisplrep`)
--------------------------------------------------------------------

For (smooth) spline-fitting to a 2-D surface, the function
:func:`bisplrep` is available. This function takes as required inputs
the **1-D** arrays *x*, *y*, and *z*, which represent points on the
surface :math:`z=f\left(x,y\right).` The default output is a list
:math:`\left[tx,ty,c,kx,ky\right]` whose entries represent
respectively, the components of the knot positions, the coefficients
of the spline, and the order of the spline in each coordinate. It is
convenient to hold this list in a single object, *tck,* so that it can
be passed easily to the function :obj:`bisplev`. The
keyword, *s* , can be used to change the amount of smoothing performed
on the data while determining the appropriate spline. The default
value is :math:`s=m-\sqrt{2m}`, where :math:`m` is the number of data
points in the *x, y,* and *z* vectors. As a result, if no smoothing is
desired, then :math:`s=0` should be passed to
:obj:`bisplrep`.

To evaluate the 2-D spline and its partial derivatives
(up to the order of the spline), the function
:obj:`bisplev` is required. This function takes as the
first two arguments **two 1-D arrays** whose cross-product specifies
the domain over which to evaluate the spline. The third argument is
the *tck* list returned from :obj:`bisplrep`. If desired,
the fourth and fifth arguments provide the orders of the partial
derivative in the :math:`x` and :math:`y` direction, respectively.

It is important to note that 2-D interpolation should not
be used to find the spline representation of images. The algorithm
used is not amenable to large numbers of input points. The signal-processing
toolbox contains more appropriate algorithms for finding
the spline representation of an image. The 2-D
interpolation commands are intended for use when interpolating a 2-D
function as shown in the example that follows. This
example uses the :obj:`mgrid <numpy.mgrid>` command in NumPy which is
useful for defining a "mesh-grid" in many dimensions. (See also the
:obj:`ogrid <numpy.ogrid>` command if the full-mesh is not
needed). The number of output arguments and the number of dimensions
of each argument is determined by the number of indexing objects
passed in :obj:`mgrid <numpy.mgrid>`.

.. plot::
   :alt: " "

   >>> import numpy as np
   >>> from scipy import interpolate
   >>> import matplotlib.pyplot as plt

   Define function over a sparse 20x20 grid

   >>> x_edges, y_edges = np.mgrid[-1:1:21j, -1:1:21j]
   >>> x = x_edges[:-1, :-1] + np.diff(x_edges[:2, 0])[0] / 2.
   >>> y = y_edges[:-1, :-1] + np.diff(y_edges[0, :2])[0] / 2.
   >>> z = (x+y) * np.exp(-6.0*(x*x+y*y))

   >>> plt.figure()
   >>> lims = dict(cmap='RdBu_r', vmin=-0.25, vmax=0.25)
   >>> plt.pcolormesh(x_edges, y_edges, z, shading='flat', **lims)
   >>> plt.colorbar()
   >>> plt.title("Sparsely sampled function.")
   >>> plt.show()

   Interpolate function over a new 70x70 grid

   >>> xnew_edges, ynew_edges = np.mgrid[-1:1:71j, -1:1:71j]
   >>> xnew = xnew_edges[:-1, :-1] + np.diff(xnew_edges[:2, 0])[0] / 2.
   >>> ynew = ynew_edges[:-1, :-1] + np.diff(ynew_edges[0, :2])[0] / 2.
   >>> tck = interpolate.bisplrep(x, y, z, s=0)
   >>> znew = interpolate.bisplev(xnew[:,0], ynew[0,:], tck)

   >>> plt.figure()
   >>> plt.pcolormesh(xnew_edges, ynew_edges, znew, shading='flat', **lims)
   >>> plt.colorbar()
   >>> plt.title("Interpolated function.")
   >>> plt.show()

..   :caption: Example of a 2-D spline interpolation.


2-D spline representation: Object-oriented (:class:`BivariateSpline`)
---------------------------------------------------------------------------------

The :class:`BivariateSpline` class is the 2-D analog of the
:class:`UnivariateSpline` class.  It and its subclasses implement
the FITPACK functions described above in an object-oriented fashion,
allowing objects to be instantiated that can be called to compute
the spline value by passing in the two coordinates as the two
arguments.


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
