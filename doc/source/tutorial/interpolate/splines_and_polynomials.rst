.. _tutorial-interpolate_splines_and_poly:

.. currentmodule:: scipy.interpolate

=====================
Piecewise polynomials
=====================

Internally, `CubicSpline` and monotone interpolants are represented as instances
of a `PPoly` class, which represents piecewise polynomials in terms of
breakpoints and coefficients (`PPoly` objects can represent polynomials of 
arbitrary orders, not only cubics). For the data array ``x``, breakpoints are at
the data points, and the array of coefficients, ``c`` , define cubic polynomials
such that ``c[k, j]`` is a coefficient for ``(x - x[j])**(3-k)`` on the segment
between ``x[j]`` and ``x[j+1]`` .


.. _tutorial-interpolate_ppoly:

Manipulating ``PPoly`` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    >>> dspl.roots() / np.pi
    array([-0.45480801,  0.50000034,  1.50000099,  2.5000016 ,  3.46249993])

This agrees well with roots :math:`\pi/2 + \pi\,n` of
:math:`\cos(x) = \sin'(x)`.
Note that by default it computed the roots *extrapolated* to the outside of
the interpolation interval :math:`0 \leqslant x \leqslant 10`, and that
the extrapolated results (the first and last values) are much less accurate.
We can switch off the extrapolation and limit the root-finding to the
interpolation interval:

    >>> dspl.roots(extrapolate=False) / np.pi
    array([0.50000034,  1.50000099,  2.5000016])

In fact, the ``root`` method is a special case of a more general ``solve``
method which finds for a given constant :math:`y` the solutions of the
equation :math:`f(x) = y` , where :math:`f(x)` is the piecewise polynomial:

    >>> dspl.solve(0.5, extrapolate=False) / np.pi
    array([0.33332755, 1.66667195, 2.3333271])

which agrees well with the expected values of  :math:`\pm\arccos(1/2) + 2\pi\,n`.

Integrals of piecewise polynomials can be computed using the ``.integrate``
method which accepts the lower and the upper limits of integration. As an
example, we compute an approximation to the complete elliptic integral
:math:`K(m) = \int_0^{\pi/2} [1 - m\sin^2 x]^{-1/2} dx`:

    >>> from scipy.special import ellipk
    >>> m = 0.5
    >>> ellipk(m)
    1.8540746773013719

To this end, we tabulate the integrand, interpolate using the monotone
interpolant (we could as well used a `CubicSpline`):

    >>> from scipy.interpolate import PchipInterpolator
    >>> x = np.linspace(0, np.pi/2, 70)
    >>> y = (1 - m*np.sin(x)**2))**(-1/2)
    >>> spl = PchipInterpolator(x, y)

and integrate

    >>> spl.integrate(0, np.pi/2)
    1.854074674965991

which is indeed close to the value computed by `scipy.special.ellipk`.

All piecewise polynomials can be constructed with N-dimensional ``y`` values.
If ``y.ndim > 1``, it is understood as a stack of 1D ``y`` values, which are
arranged along the interpolation axis (with the default value of 0).
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
    >>> plt.legend(['`ellipk`', 'integrated piecewise polynomial'])
    >>> plt.show()


Interpolation with B-splines
----------------------------

A polynomial of degree :math:`k` can be thought of as a linear combination of
:math:`k+1` monomial basis elements, :math:`1, x, x^2, \cdots, x^k`. B-splines
form an alternative (if equivalent) :ref:`basis <tutorial-interpolate_bspl_basis>`
of degree-:math:`k` piecewise polynomials. 

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
coefficients. Some routines (see the :ref:`Smoothing splines section
<tutorial-interpolate_fitpack>`) zero-pad the ``c`` arrays so that
``len(c) == len(t)``. These additional coefficients are ignored for evaluation.

Also note that the coeffients are given in the
:ref:`b-spline basis <tutorial-interpolate_bspl_basis>`, not the power basis
of :math:`1, x, \cdots, x^k`.


.. _tutorial-interpolate_bspl_basis:

B-spline basis elements
-----------------------

B-splines are piecewise polynomials, represented as linear combinations of
*b-spline basis elements* --- which themselves are certain linear combinations
of usual monomials, :math:`x^m` with :math:`m=0, 1, \dots, k`.

The b-spline basis are generally more computationally stable then the power basis
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


.. _tutorial-interpolate_parametric:

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

