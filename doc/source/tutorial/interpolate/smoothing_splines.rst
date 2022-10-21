.. _tutorial-interpolate_fitpack:

.. currentmodule:: scipy.interpolate

=================
Smoothing splines
=================

For the interpolation problem, the task is to construct a curve which passes
through a given set of data points. This may be not appropriate if the data is
noisy: we then want to construct a smooth curve, ``g(x)``, which *approximates*
the input data without passing through each point exactly.
To this end, `scipy.interpolate` allows constructing *smoothing splines*, based
on the Fortran library FITPACK by P. Dierckx.

Specifically, given the data arrays ``x`` and ``y`` and the array of
non-negative *weights*, ``w``, we look for a spline function ``g(x)`` which
satisfies

.. math::

    \sum_j \left[ w_j (g(x_j) - y_j)\right]^2 \leqslant s

where ``s`` is the input parameter which controls the interplay between the
smoothness of the resulting curve ``g(x)`` and the quality of the approximation
of the data (i.e., the differences between :math:`g(x_j)` and :math:`y_j`).

Note that the limit ``s = 0`` corresponds to the interpolation problem where
:math:`g(x_j) = y_j`. Increasing ``s`` leads to smoother fits, and in the limit
of a very large ``s``, :math:`g(x)` degenerates into a single best-fit polynomial.

Finding a good value of the ``s`` parameter is a trial-and-error process. If
the weights correspond to the inverse of standard deviations of the input data,
the "good" value of ``s`` is expected to be somewhere between :math:`m - \sqrt{2m}`
and :math:`m + \sqrt{2m}`, where :math:`m` is the number of the data points.
If all weights equal unity, a reasonable choice might be around :math:`s \sim m\,\sigma^2`,
where :math:`\sigma` is an estimate for the standard deviation of the data. 

Internally, the FITPACK library works by adding internal knots to the spline
fit ``g(x)``, so that **the resulting knots do not necessarily coincide with the input data**.


Spline smoothing in 1-D
-----------------------

`scipy.interpolate` provides two interfaces for the FITPACK library, a functional
interface and an object-oriented interface. While equivalent, these interfaces
have different defaults. Below we discuss them in turn, starting from the
functional interface --- which we recommend for use in new code.


.. _tutorial-interpolate_splXXX:

Procedural (interpolate.splXXX)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

The knot array defines the interpolation interval to be ``t[k:-k]``, so that 
first :math:`k+1` and last :math:`k+1` entries of the ``t`` array define
*boundary knots*. The coefficients are a 1D array of length at least
``len(t) - k - 1``. Some routines pad this array to have ``len(c) == len(t)``---
these additional coefficients are ignored for the spline evaluation. 

The ``tck``-tuple format is compatible with
:ref:`interpolating b-splines <tutorial-interpolate_bspl_basis>`: the output of
`splrep` can be wrapped into a `BSpline` object, e.g. ``BSpline(*tck)``, and
the evaluation/integration/root-finding routines described below
can use ``tck``-tuples and `BSpline` objects interchangeably.

For curves in N-D space the function
:obj:`splprep` allows defining the curve
parametrically. For this function only 1 input argument is
required. This input is a list of :math:`N` arrays representing the
curve in N-D space. The length of each array is the
number of curve points, and each array provides one component of the
N-D data point. The parameter variable is given
with the keyword argument, ``u``, which defaults to an equally-spaced
monotonic sequence between :math:`0` and :math:`1` (i.e., the :ref:`uniform
parametrization <tutorial-interpolate_parametric>`).

The output consists of two objects: a 3-tuple, :math:`\left(t,c,k\right)`
, containing the spline representation and the parameter variable
:math:`u.`

The coeffients are a list of :math:`N` arrays, where each array corresponds to
a dimension of the input data. Note that the knots, ``t`` correspond to the
parametrization of the curve ``u``. 

The keyword argument, ``s`` , is used to specify the amount of smoothing
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

   Cubic spline

   >>> x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/8)
   >>> y = np.sin(x)
   >>> tck = interpolate.splrep(x, y, s=0)
   >>> xnew = np.arange(0, 2*np.pi, np.pi/50)
   >>> ynew = interpolate.splev(xnew, tck, der=0)

   Note that the last line is equivalent to ``BSpline(*tck)(xnew)``.

   >>> plt.figure()
   >>> plt.plot(x, y, 'x', xnew, ynew, xnew, np.sin(xnew), x, y, 'b')
   >>> plt.legend(['Linear', 'Cubic Spline', 'True'])
   >>> plt.axis([-0.05, 6.33, -1.05, 1.05])
   >>> plt.title('Cubic-spline interpolation')
   >>> plt.show()

   Derivative of spline

   >>> yder = interpolate.splev(xnew, tck, der=1)   # or BSpline(*tck)(xnew, 1)
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

Note that in the last example, `splprep` returns the spline coefficients as a
list of arrays, where each array corresponds to a dimension of the input data.
Thus to wrap its output to a `BSpline`, we need to transpose the coefficients
(or use ``BSpline(..., axis=1)``):

  >>> tt, cc, k = tck
  >>> cc = np.array(cc)
  >>> bspl = BSpline(tt, cc.T, k)    # note the transpose
  >>> xy = bspl(u)
  >>> xx, yy = xy.T   # transpose to unpack into a pair of arrays
  >>> np.allclose(x, xx)
  True
  >>> np.allclose(y, yy)
  True

So far all examples constructed interpolating splines with ``s=0``. To illustrate
the effect of the value of ``s`` for noisy data:

.. plot::

   >>> import numpy as np
   >>> from scipy.interpolate import splrep, splev
   >>> x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/16)
   >>> rng = np.random.default_rng()
   >>> y =  np.sin(x) + 0.4*rng.standard_normal(size=len(x))
   >>> tck = splrep(x, y, s=0)
   >>> tck_s = splrep(x, y, s=len(x))
   >>> import matplotlib.pyplot as plt
   >>> xnew = np.arange(0, 2*np.pi, np.pi/50)
   >>> plt.plot(xnew, np.sin(xnew), '-.', label='sin(x)')
   >>> plt.plot(xnew, splev(xnew, tck), '-', label='s=0')
   >>> plt.plot(xnew, splev(xnew, tck_s), '-', label=f's={len(x)}')
   >>> plt.plot(x, y, 'o')
   >>> plt.legend()
   >>> plt.show()

We see that the ``s=0`` curve follows the (random) fluctuations of the data points,
while the ``s > 0`` curve is close to the underlying sine function.

The default value of ``s`` depends on whether the weights are supplied or not,
and also differs for `splrep` and `splprep`. Therefore, we recommend to always
supply the value of ``s`` explicitly.


Object-oriented (:class:`UnivariateSpline`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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



2-D smoothing splines
---------------------

In addition to smoothing 1-D splines, the FITPACK library provides the means of
fitting 2-D *surfaces* to two-dimensional data, represented as tensor products
of 1-D splines. There are also two interfaces: a procedural interface and an
object-oriented interface.


.. _tutorial-interpolate_2d_spline:

Procedural (:func:`bisplrep`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For (smooth) spline-fitting to a 2-D surface, the function
:func:`bisplrep` is available. This function takes as required inputs
the **1-D** arrays ``x``, ``y``, and ``z``, which represent points on the
surface :math:`z=f(x, y).` The default output is a list
``[tx ,ty, c, kx, ky]`` whose entries represent
respectively, the components of the knot positions, the coefficients
of the spline, and the order of the spline in each coordinate. It is
convenient to hold this list in a single object, ``tck``, so that it can
be passed easily to the function :obj:`bisplev`. The
keyword, ``s`` , can be used to change the amount of smoothing performed
on the data while determining the appropriate spline. The default
value is :math:`s=m-\sqrt{2m}`, where :math:`m` is the number of data
points in the ``x``, ``y``, and ``z`` vectors. As a result, if no smoothing is
desired, then ``s=0`` should be passed to :obj:`bisplrep`.

To evaluate the 2-D spline and its partial derivatives
(up to the order of the spline), the function
:obj:`bisplev` is required. This function takes as the
first two arguments **two 1-D arrays** whose cross-product specifies
the domain over which to evaluate the spline. The third argument is
the ``tck`` list returned from :obj:`bisplrep`. If desired,
the fourth and fifth arguments provide the orders of the partial
derivative in the :math:`x` and :math:`y` direction, respectively.

It is important to note that 2-D interpolation should not
be used to find the spline representation of images. The algorithm
used is not amenable to large numbers of input points. `scipy.signal`
and `scipy.ndimage` contain more appropriate algorithms for finding
the spline representation of an image. 

The 2-D interpolation commands are intended for use when interpolating a 2-D
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


Object-oriented (:class:`BivariateSpline`)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`BivariateSpline` class is the 2-D analog of the
:class:`UnivariateSpline` class.  It and its subclasses implement
the FITPACK functions described above in an object-oriented fashion,
allowing objects to be instantiated that can be called to compute
the spline value by passing in the two coordinates as the two
arguments.

