.. _tutorial-interpolate_lsq-ndbspline:

Least-squares tensor-product splines
====================================

`make_lsq_ndbspline` fits a tensor-product B-spline to scattered data in
one or more dimensions. It returns an `NdBSpline`, which can then be evaluated
at new coordinate points. The fit and evaluation are separate operations:

.. code-block:: python

    spl = make_lsq_ndbspline(x, y, t, k=1)
    y_new = spl(x_new)

Here, the rows of ``x`` are the coordinates of the observations, and the first
axis of ``y`` identifies the corresponding observed values. For example, a
two-dimensional fit has ``x.shape == (n_samples, 2)``.

Knots and boundary convention
-----------------------------

The user supplies one complete knot vector ``t[d]`` for each dimension. The
function does not add boundary knots or otherwise pad these vectors. For
degree ``k[d]``, ``t[d]`` must be a finite, one-dimensional, non-decreasing
array containing at least ``2*k[d] + 2`` values, with
``t[d][k[d]] < t[d][-k[d] - 1]``. The latter two values delimit the base
interval in that dimension. A vector of length ``nt`` defines
``nt - k[d] - 1`` one-dimensional basis functions.

A common :math:`(k+1)`-regular construction repeats each endpoint
``k + 1`` times and places the interior knots between them. For example, a
cubic (``k=3``) knot vector on ``[0, 1]`` with three interior knots is:

.. code-block:: python

    >>> import numpy as np
    >>> k_cubic = 3
    >>> interior_knots = np.array([0.25, 0.5, 0.75])
    >>> t_cubic = np.concatenate((
    ...     np.repeat(0.0, k_cubic + 1),
    ...     interior_knots,
    ...     np.repeat(1.0, k_cubic + 1),
    ... ))
    >>> t_cubic.size - k_cubic - 1
    7

This convention differs from interfaces that add the repeated boundary knots
internally. For example, R's `splines::bs
<https://stat.ethz.ch/R-manual/R-devel/library/splines/html/bs.html>`_
accepts interior knots through ``knots`` and the two boundaries separately
through ``Boundary.knots``. Here, all of these values must already be combined
in ``t[d]``.

The knot vectors determine the locations and boundary behavior of the basis
functions. The following example uses linear splines and therefore repeats
each boundary twice:

.. code-block:: python

    >>> from scipy.interpolate import make_lsq_ndbspline
    >>> rng = np.random.default_rng(1234)
    >>> x = rng.uniform(0.0, 1.0, size=(100, 2))
    >>> y = (np.sin(np.pi*x[:, 0]) + x[:, 1]
    ...      + rng.normal(0.0, 0.05, x.shape[0]))
    >>> t = (
    ...     np.array([0.0, 0.0, 0.5, 1.0, 1.0]),
    ...     np.array([0.0, 0.0, 0.5, 1.0, 1.0]),
    ... )
    >>> spl = make_lsq_ndbspline(x, y, t, k=1)

Coordinates for evaluation follow the same row-wise convention:

.. code-block:: python

    >>> x_new = np.array([[0.25, 0.75], [0.50, 0.25]])
    >>> spl(x_new).shape
    (2,)

Multiple outputs
----------------

Trailing dimensions of ``y`` represent multiple outputs observed at the same
coordinates. They are fitted with the same design matrix. For example, two
outputs can be stacked as columns, giving ``y_multi.shape == (100, 2)``:

.. code-block:: python

    >>> y_multi = np.column_stack((
    ...     np.sin(np.pi*x[:, 0]) + x[:, 1],
    ...     x[:, 0] - x[:, 1]**2,
    ... ))
    >>> spl_multi = make_lsq_ndbspline(x, y_multi, t, k=1)
    >>> spl_multi(x_new).shape
    (2, 2)

The first result dimension corresponds to the two evaluation points and the
trailing dimension corresponds to the two fitted outputs. More generally,
``y`` may have shape ``(n_samples, ...)``, and evaluating the resulting spline
at ``n_eval`` points produces an array with shape ``(n_eval, ...)``.

Weighted fitting
----------------

For design matrix :math:`A`, coefficient vector :math:`c`, observations
:math:`y`, and :math:`W = \operatorname{diag}(w)`, the fitted coefficients
minimize

.. math::

   \lVert W(Ac-y) \rVert_2.

Consequently, ``w`` multiplies the residual directly. A point with zero weight
does not contribute to the fit. Each basis function must still be supported by
at least one point with positive weight.

Data on a rectilinear grid
--------------------------

For values sampled on a complete rectilinear grid,
`make_lsq_ndbspline_from_grid` accepts the one-dimensional grid axes directly.
There is no need to construct the coordinate rows manually:

.. code-block:: python

    >>> from scipy.interpolate import make_lsq_ndbspline_from_grid
    >>> points = (np.linspace(0.0, 1.0, 8), np.linspace(-1.0, 1.0, 9))
    >>> x0, x1 = np.meshgrid(*points, indexing="ij")
    >>> values = np.sin(np.pi*x0) + x1
    >>> t_grid = (
    ...     np.array([0.0, 0.0, 0.5, 1.0, 1.0]),
    ...     np.array([-1.0, -1.0, 0.0, 1.0, 1.0]),
    ... )
    >>> spl_grid = make_lsq_ndbspline_from_grid(
    ...     points, values, t_grid, k=1
    ... )
    >>> spl_grid([[0.25, 0.5]]).shape
    (1,)

In ``N`` dimensions, ``points`` contains ``N`` one-dimensional arrays and
``values.shape[:N]`` must match their lengths. Both constructors also accept
trailing dimensions in the data values, allowing several outputs to be fitted
with the same coordinates and basis. For example, values with shape
``(m1, ..., mN, 2)`` fit two outputs on the same grid.
