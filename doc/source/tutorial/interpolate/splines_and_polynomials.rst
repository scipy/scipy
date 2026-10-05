.. _tutorial-interpolate_splines_and_poly:

.. currentmodule:: scipy.interpolate

=================================
Piecewise polynomials and splines
=================================

1D interpolation routines :ref:`discussed in the previous section
<tutorial-interpolate_1Dsection>`, work by constructing certain *piecewise
polynomials*: the interpolation range is split into intervals by the so-called
*breakpoints*, and there is a certain polynomial on each interval. These
polynomial pieces then match at the breakpoints with a predefined smoothness:
the second derivatives for cubic splines, the first derivatives for monotone
interpolants and so on.

A polynomial of degree :math:`k` can be thought of as a linear combination of
:math:`k+1` monomial basis elements, :math:`1, x, x^2, \cdots, x^k`. 
In some applications, it is useful to consider alternative (if formally
equivalent) bases. Two popular bases, implemented in `scipy.interpolate` are
B-splines (`BSpline`) and Bernstein polynomials (`BPoly`).
B-splines are often used, for example, in non-parametric regression problems,
and Bernstein polynomials are used for constructing Bezier curves.

`PPoly` objects represent piecewise polynomials in the 'usual' power basis.
This is the case for `CubicSpline` instances and monotone interpolants.
In general, `PPoly` objects can represent polynomials of 
arbitrary orders, not only cubics. For the data array ``x``, breakpoints are at
the data points, and the array of coefficients, ``c`` , define polynomials of
degree :math:`k`, such that ``c[i, j]`` is a coefficient for
``(x - x[j])**(k-i)`` on the segment between ``x[j]`` and ``x[j+1]`` .

`BSpline` objects represent B-spline functions --- linear combinations of
:ref:`b-spline basis elements <tutorial-interpolate_bspl_basis>`. 
These objects can be instantiated directly or constructed from data with the
`make_interp_spline` factory function.

Finally, Bernstein polynomials are represented as instances of the `BPoly` class.

All these classes implement a (mostly) similar interface, `PPoly` being the most
feature-complete. We next consider the main features of this interface and
discuss some details of the alternative bases for piecewise polynomials.


.. _tutorial-interpolate_ppoly:

Manipulating `PPoly` objects
============================

`PPoly` objects have convenient methods for constructing derivatives
and antiderivatives, computing integrals and root-finding. For example, we
tabulate the sine function and find the roots of its derivative.

    >>> import numpy as np
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

Note that the second form above evaluates the derivative in place, while with
the ``dspl`` object, we can find the zeros of the derivative of ``spl``:

    >>> dspl.roots() / np.pi
    array([-0.45480801,  0.50000034,  1.50000099,  2.5000016 ,  3.46249993])

This agrees well with the roots :math:`\pi/2 + \pi\,n` of
:math:`\cos(x) = \sin'(x)`.
Note that by default it computes the roots *extrapolated* to the outside of
the interpolation interval :math:`0 \leqslant x \leqslant 10`, and that
the extrapolated results (the first and last values) are much less accurate.
We can switch off the extrapolation and limit the root-finding to the
interpolation interval:

    >>> dspl.roots(extrapolate=False) / np.pi
    array([0.50000034,  1.50000099,  2.5000016])

In fact, the ``roots`` method is a special case of a more general ``solve``
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

To this end, we tabulate the integrand and interpolate it using the monotone
PCHIP interpolant (we could as well have used a `CubicSpline`):

    >>> from scipy.interpolate import PchipInterpolator
    >>> x = np.linspace(0, np.pi/2, 70)
    >>> y = (1 - m*np.sin(x)**2)**(-1/2)
    >>> spl = PchipInterpolator(x, y)

and integrate

    >>> spl.integrate(0, np.pi/2)
    1.854074674965991

which is indeed close to the value computed by `scipy.special.ellipk`.

All piecewise polynomials can be constructed with N-dimensional ``y`` values.
If ``y.ndim > 1``, it is interpreted as a stack of 1D ``y`` values, which are
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


B-splines: knots and coefficients
=================================

A b-spline function --- for instance, constructed from data via a
`make_interp_spline` call --- is defined by the so-called *knots* and coefficients.

As an illustration, let us again construct the interpolation of a sine function. 
The knots are available as the ``t`` attribute of a `BSpline` instance:

    >>> x = np.linspace(0, 3/2, 7)
    >>> y = np.sin(np.pi*x)
    >>> from scipy.interpolate import make_interp_spline
    >>> bspl = make_interp_spline(x, y, k=3)
    >>> print(bspl.t)
    [0.  0.  0.  0.        0.5  0.75  1.        1.5  1.5  1.5  1.5 ]
    >>> print(x)
    [            0.  0.25  0.5  0.75  1.  1.25  1.5 ]

We see that the knot vector by default is constructed from the input
array ``x``: first, it is made :math:`(k+1)` -regular (it has ``k``
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

We stress that the coefficients are given in the
:ref:`b-spline basis <tutorial-interpolate_bspl_basis>`, not the power basis
of :math:`1, x, \cdots, x^k`.


.. _tutorial-interpolate_bspl_control_points:

Coefficients and control points
-------------------------------

The b-spline coefficients have a simple geometric interpretation: they are
the *control points* of the spline. The polyline connecting the control
points --- the *control polygon* --- follows the shape of the spline, and the
spline lies within the convex hull of its control points. How exactly the
coefficients map onto the control points depends on whether the spline
represents a function or a curve.

**Spline functions.** For a function, :math:`y = f(x)`, the coefficients are
scalars, and they only give the :math:`y`-coordinates of the control points.
The matching :math:`x`-coordinates are the so-called *Greville abscissae*,
which are averages of ``k`` consecutive knots,

.. math::

    \xi_j = \frac{t_{j+1} + t_{j+2} + \dots + t_{j+k}}{k} ,
    \qquad j = 0, 1, \dots, n-1 ,

where ``n = len(t) - k - 1`` is the number of coefficients.

.. plot::

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.interpolate import make_interp_spline
    >>> x = np.linspace(0, 2*np.pi, 8)
    >>> spl = make_interp_spline(x, np.sin(x), k=3)
    >>> n = len(spl.t) - spl.k - 1
    >>> xg = np.array([spl.t[j+1:j+spl.k+1].mean() for j in range(n)])

    The control points are ``(xg[j], spl.c[j])``. Plot them together with the
    spline:

    >>> xx = np.linspace(x[0], x[-1], 200)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(xx, spl(xx), label='spline')
    >>> ax.plot(xg, spl.c, 'o--', label='control polygon')
    >>> ax.plot(x, np.sin(x), 'kx', label='data')
    >>> ax.set_title('y = f(x): coefficients are heights')
    >>> ax.legend()
    >>> plt.show()

**Parametric curves.** For a curve in :math:`d` dimensions, e.g. constructed
by `make_splprep`, each coefficient is itself a point: the coefficient array
has shape ``(n, d)``, and its rows are the control points.

.. plot::

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.interpolate import BSpline, make_splprep
    >>> theta = np.linspace(0, 1.5*np.pi, 10)
    >>> data = [theta*np.cos(theta), theta*np.sin(theta)]
    >>> spl, u = make_splprep(data, s=0)
    >>> spl.c.shape
    (10, 2)

    Since a basis element is non-zero on at most ``k+1`` knot intervals,
    moving a single control point only changes the curve locally:

    >>> c_moved = spl.c.copy()
    >>> c_moved[5] += [1.5, 1.5]
    >>> spl_moved = BSpline(spl.t, c_moved, spl.k)

    Plot the curve together with its control polygon (left), and the effect of
    moving a control point (right):

    >>> uu = np.linspace(0, 1, 300)
    >>> fig, axs = plt.subplots(1, 2, figsize=(7, 5), layout='constrained')
    >>> ax = axs[0]
    >>> ax.plot(*spl(uu), label='spline')
    >>> ax.plot(*spl.c.T, 'o--', label='control polygon')
    >>> ax.plot(*data, 'kx', label='data')
    >>> ax.set_title('parametric curve: rows of c are points')
    >>> ax.set_aspect('equal')
    >>> ax.legend(loc='center')
    >>> ax = axs[1]
    >>> ax.plot(*spl(uu), color='C0', alpha=0.4, label='original')
    >>> ax.plot(*spl_moved(uu).T, color='C0', label='moved')
    >>> ax.plot(*c_moved.T, 'o--', color='C1', label='control polygon')
    >>> ax.plot(*spl.c[5], 's', color='C1', alpha=0.4)
    >>> ax.set_title('moving control point 5')
    >>> ax.set_aspect('equal')
    >>> ax.legend(loc='center')
    >>> plt.show()

Note the transposes: ``spl.c`` has shape ``(n, d)``, so ``spl.c.T`` unpacks
into the ``x`` and ``y`` coordinates of the control points. ``spl(uu)``
already has shape ``(d, len(uu))``, because `make_splprep` constructs the
spline with ``axis=1``, while ``spl_moved``, constructed directly from the
``(n, d)`` coefficient array, evaluates to shape ``(len(uu), d)``.

**Spline surfaces.** A tensor product spline surface, :math:`z = f(x, y)`, e.g.
constructed by `RectBivariateSpline`, is a linear combination of products of
b-spline basis elements in the two directions,

.. math::

    f(x, y) = \sum_{i=0}^{n_x-1} \sum_{j=0}^{n_y-1} c_{ij} B_i(x) B_j(y) ,

where ``tx`` and ``ty`` are the knot vectors, and ``kx`` and ``ky`` are the
degrees in the :math:`x` and :math:`y` directions. Just like in the 1D case,
each direction has its own b-spline basis, with ``nx = len(tx) - kx - 1``
elements :math:`B_i(x)` and ``ny = len(ty) - ky - 1`` elements :math:`B_j(y)`.
The coefficients therefore form a 2D array, :math:`c_{ij}`, of shape
``(nx, ny)``: there is one coefficient for each pair of basis elements.

The surface is the graph of a function of two variables, so the geometric
picture is a direct generalization of the 1D spline functions above. Each
coefficient :math:`c_{ij}` is a scalar, and it only gives the height, i.e. the
:math:`z`-coordinate, of a control point. The :math:`x`- and
:math:`y`-coordinates are the Greville abscissae of the two knot vectors,
computed separately in each direction,

.. math::

    \xi_i = \frac{t^x_{i+1} + \dots + t^x_{i+k_x}}{k_x} ,
    \qquad
    \eta_j = \frac{t^y_{j+1} + \dots + t^y_{j+k_y}}{k_y} ,

for :math:`i = 0, \dots, n_x - 1` and :math:`j = 0, \dots, n_y - 1`. The
control points are thus

.. math::

    P_{ij} = (\xi_i, \eta_j, c_{ij}) ,

and they sit over a rectangular grid in the :math:`(x, y)` plane, the tensor
product of the 1D Greville grids. Connecting neighboring control points in each
direction gives the *control net*, which is the 2D analog of the control
polygon: it follows the shape of the surface, the surface lies in the convex
hull of the control points, and changing :math:`c_{ij}` only changes the
surface locally, on the rectangle where :math:`B_i(x) B_j(y)` is non-zero.

.. plot::

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.interpolate import RectBivariateSpline, BSpline
    >>> x = np.linspace(0, 4, 7)
    >>> y = np.linspace(0, 3, 6)
    >>> X, Y = np.meshgrid(x, y, indexing="ij")
    >>> rbs = RectBivariateSpline(x, y, np.sin(X) * np.cos(Y))

    The ``get_coeffs`` method returns a flat array of length ``nx * ny``, with
    the :math:`x` index varying slowest. Reshape it into an ``(nx, ny)`` grid:

    >>> tx, ty = rbs.get_knots()
    >>> kx, ky = rbs.degrees
    >>> nx, ny = len(tx) - kx - 1, len(ty) - ky - 1
    >>> rbs.get_coeffs().shape
    (42,)
    >>> c = rbs.get_coeffs().reshape(nx, ny)
    >>> c.shape
    (7, 6)

    To check that ``c[i, j]`` is the coefficient of :math:`B_i(x) B_j(y)`,
    evaluate the sum above using the design matrices of the 1D b-spline bases,
    ``Bx[p, i] = B_i(xx[p])`` and ``By[q, j] = B_j(yy[q])``:

    >>> xx = np.linspace(0, 4, 41)
    >>> yy = np.linspace(0, 3, 31)
    >>> Bx = BSpline.design_matrix(xx, tx, kx).toarray()
    >>> By = BSpline.design_matrix(yy, ty, ky).toarray()
    >>> np.allclose(Bx @ c @ By.T, rbs(xx, yy))
    True

    The control points sit at the Greville abscissae in each direction:

    >>> xg = np.array([tx[i+1:i+kx+1].mean() for i in range(nx)])
    >>> yg = np.array([ty[j+1:j+ky+1].mean() for j in range(ny)])
    >>> XG, YG = np.meshgrid(xg, yg, indexing="ij")

    As a check, for data sampled from a plane, the coefficients equal the
    plane evaluated at the Greville abscissae, so that the control net lies
    exactly on the surface:

    >>> plane = RectBivariateSpline(x, y, 1 + 2*X - 3*Y)
    >>> c_plane = plane.get_coeffs().reshape(nx, ny)
    >>> np.allclose(c_plane, 1 + 2*XG - 3*YG)
    True

    In general, the control net only approximates the surface. Since the
    boundary knots are repeated ``k+1`` times, the first and last Greville
    abscissae are the end points of the interval, and only a single basis
    element is non-zero there. Hence the control net touches the surface at
    the four corners:

    >>> xg[[0, -1]], yg[[0, -1]]
    (array([0., 4.]), array([0., 3.]))
    >>> ix, jy = [0, 0, -1, -1], [0, -1, 0, -1]
    >>> np.allclose(c[ix, jy], rbs(xg[ix], yg[jy], grid=False))
    True

    Finally, plot the surface together with its control net. The control
    points, ``(XG[i, j], YG[i, j], c[i, j])``, sit over the Greville grid,
    ``(xg, yg)``, rather than over the grid of data points, ``(x, y)``. The two
    grids are related, but they are not the same. Since the spline interpolates
    the data, there is one control point per data point, so both grids have the
    same size, and both start and end at the boundaries of the data. The
    interior Greville abscissae are, however, averages of the knots, and are
    shifted with respect to the data points:

    >>> x
    array([0.        , 0.66666667, 1.33333333, 2.        , 2.66666667,
           3.33333333, 4.        ])
    >>> xg
    array([0.        , 0.44444444, 1.11111111, 2.        , 2.88888889,
           3.55555556, 4.        ])

    Lines of the net connect neighboring control points in the :math:`x` and
    :math:`y` directions. The plot shows the properties discussed above:

    - the control net follows the overall shape of the surface, i.e. its
      peaks, valleys and saddle;
    - the net is not on the surface: it exaggerates the shape, so that the
      control points rise above the peaks and dip below the valleys. Indeed,
      ``c.max()`` and ``c.min()`` exceed the extremes of the surface itself;
    - the surface is smoother than the net, and it lies within the convex hull
      of the control points;
    - the control points at the four corners lie exactly on the surface.

    >>> XX, YY = np.meshgrid(xx, yy, indexing="ij")
    >>> fig = plt.figure(figsize=(7, 5))
    >>> ax = fig.add_subplot(projection="3d")
    >>> ax.plot_surface(XX, YY, rbs(xx, yy), cmap="viridis", alpha=0.6)
    >>> ax.plot_wireframe(XG, YG, c, color="C1")
    >>> ax.scatter(XG, YG, c, color="C1", label="control points")
    >>> ax.set_title("surface: c[i, j] are heights over the Greville grid")
    >>> ax.legend()
    >>> plt.show()


.. _tutorial-interpolate_bspl_basis:

B-spline basis elements
-----------------------

The b-spline basis is used in a variety of applications which include interpolation,
regression and curve representation.
B-splines are piecewise polynomials, represented as linear combinations of
*b-spline basis elements* --- which themselves are certain linear combinations
of usual monomials, :math:`x^m` with :math:`m=0, 1, \dots, k`.

The properties of b-splines are well described in the literature (see, for example,
references listed in the `BSpline` docstring). For our purposes, it is enough to know
that a b-spline function is uniquely defined by an array of coefficients and
an array of the so-called *knots*, which may or may not coincide with the data points,
``x``.

Specifically, a b-spline basis element of degree ``k`` (e.g. ``k=3`` for cubics)
is defined by :math:`k+2` knots and is zero outside of these knots.
To illustrate, plot a collection of non-zero basis elements on a certain
interval:

.. plot::

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


.. _tutorial-interpolate_bspl_design_matrix:

Design matrices in the B-spline basis
-------------------------------------

One common application of b-splines is in non-parametric regression. The reason
is that the localized nature of the b-spline basis elements makes linear
algebra banded. This is because at most :math:`k+1` basis elements are non-zero
at a given evaluation point, thus a design matrix built on b-splines has at most
:math:`k+1` diagonals.

As an illustration, we consider a toy example. Suppose our data are
one-dimensional and are confined to an interval :math:`[0, 6]`.
We construct a 4-regular knot vector which corresponds to 7 data points and
cubic, ``k=3``, splines:

>>> t = [0., 0., 0., 0., 2., 3., 4., 6., 6., 6., 6.]

Next, take 'observations' to be

>>> xnew = [1, 2, 3]

and construct the design matrix in the sparse CSR format

>>> from scipy.interpolate import BSpline
>>> mat = BSpline.design_matrix(xnew, t, k=3)
>>> mat
<Compressed Sparse Row sparse array of dtype 'float64'
	with 12 stored elements and shape (3, 7)>

Here each row of the design matrix corresponds to a value in the ``xnew`` array,
and a row has no more than ``k+1 = 4`` non-zero elements; row ``j``
contains basis elements evaluated at ``xnew[j]``:

>>> with np.printoptions(precision=3):
...     print(mat.toarray())
[[0.125 0.514 0.319 0.042 0.    0.    0.   ]
 [0.    0.111 0.556 0.333 0.    0.    0.   ]
 [0.    0.    0.125 0.75  0.125 0.    0.   ]]


Bernstein polynomials, ``BPoly``
================================

For :math:`t \in [0, 1]`, Bernstein basis polynomials of degree :math:`k` are defined via

.. math::

    b(t; k, a) = C_k^a t^a (1-t)^{k - a}

where :math:`C_k^a` is the binomial coefficient, and :math:`a=0, 1, \dots, k`, so that
there are :math:`k+1` basis polynomials of degree :math:`k`.

A ``BPoly`` object represents a *piecewise* Bernstein polynomial in terms of
breakpoints, ``x``, and coefficients, ``c``: ``c[a, j]`` gives the coefficient for
:math:`b(t; k, a)` for ``t`` on the interval between ``x[j]`` and ``x[j+1]``.

The user interface of `BPoly` objects is very similar to that of `PPoly` objects:
both can be evaluated, differentiated and integrated.

One additional feature of `BPoly` objects is the alternative constructor,
`BPoly.from_derivatives`, which constructs a `BPoly` object from data values and derivatives.
Specifically, ``b = BPoly.from_derivatives(x, y)`` returns a callable that interpolates
the provided values, ``b(x[i]) == y[i])``, and has the provided derivatives,
``b(x[i], nu=j) == y[i][j]``.

This operation is similar to `CubicHermiteSpline`, but it is more flexible in that
it can handle varying numbers of derivatives at different data points; i.e., the ``y``
argument can be a list of arrays of different lengths. See `BPoly.from_derivatives`
for further discussion and examples.


Conversion between bases
========================

In principle, all three bases for piecewise polynomials (the power basis, the Bernstein
basis, and b-splines) are equivalent, and a polynomial in one basis can be converted
into a different basis. One reason for converting between bases is that not all bases
implement all operations. For instance, root-finding is only implemented for `PPoly`,
and therefore to find roots of a `BSpline` object, you need to convert to `PPoly` first.
See methods `PPoly.from_bernstein_basis`, `PPoly.from_spline`,
`BPoly.from_power_basis`, and `BSpline.from_power_basis` for details about conversion.

In floating-point arithmetic, though, conversions always incur some precision loss.
Whether this is significant is problem-dependent, so it is therefore recommended to
exercise caution when converting between bases.
