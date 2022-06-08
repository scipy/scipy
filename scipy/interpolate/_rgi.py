__all__ = ['RegularGridInterpolator', 'interpn']

import itertools

import numpy as np

from .interpnd import _ndim_coords_from_arrays
from ._cubic import PchipInterpolator
from ._bsplines import make_interp_spline
from ._fitpack2 import RectBivariateSpline


def _make_points_and_values_ascending(points, values):
    # create ascending points
    sorted_indexes = tuple(np.argsort(point) for point in points)
    points_asc = tuple(
        np.asarray(point)[sort_index] for (point, sort_index) in
        zip(points, sorted_indexes))

    # create ascending values
    ordered_indexes = tuple([*range(len(x))] for x in sorted_indexes)
    ordered_indexes_array = np.array(
        [i.flatten() for i in np.meshgrid(*ordered_indexes)]).transpose()
    sorted_indexes_array = np.array(
        [i.flatten() for i in np.meshgrid(*sorted_indexes)]).transpose()

    values_asc = np.zeros_like(np.asarray(values))
    for o, s in zip(ordered_indexes_array, sorted_indexes_array):
        values_asc[tuple(o)] = values[tuple(s)]

    return points_asc, values_asc


class RegularGridInterpolator:
    """
    Interpolation on a regular or rectilinear grid in arbitrary dimensions.

    The data must be defined on a rectilinear grid; that is, a rectangular
    grid with even or uneven spacing. Linear, nearest-neighbor, spline
    interpolations are supported. After setting up the interpolator object,
    the interpolation method may be chosen at each evaluation.

    Parameters
    ----------
    points : tuple of ndarray of float, with shapes (m1, ), ..., (mn, )
        The points defining the regular grid in n dimensions. The points in
        each dimension (i.e. every elements of the points tuple) must be
        strictly ascending or descending.

    values : array_like, shape (m1, ..., mn, ...)
        The data on the regular grid in n dimensions. Complex data can be
        acceptable.

    method : str, optional
        The method of interpolation to perform. Supported are "linear",
        "nearest", "slinear", "cubic", and "quintic". This parameter will
        become the default for the object's ``__call__`` method.
        Default is "linear".

    bounds_error : bool, optional
        If True, when interpolated values are requested outside of the
        domain of the input data, a ValueError is raised.
        If False, then `fill_value` is used.
        Default is True.

    fill_value : float or None, optional
        The value to use for points outside of the interpolation domain.
        If None, values outside the domain are extrapolated.
        Default is ``np.nan``.

    Methods
    -------
    __call__

    Attributes
    ----------
    grid : tuple of ndarrays
        The points defining the regular grid in n dimensions.
        This tuple defines the full grid via
        ``np.meshgrid(*grid, indexing='ij')``
    values : ndarray
        Data values at the grid.
    method : str
        Interpolation method.
    fill_value : float or ``None``
        Use this value for out-of-bounds arguments to `__call__`.
    bounds_error : bool
        If ``True``, out-of-bounds argument raise a ``ValueError``.

    Notes
    -----
    Contrary to `LinearNDInterpolator` and `NearestNDInterpolator`, this class
    avoids expensive triangulation of the input data by taking advantage of the
    regular grid structure.

    In other words, this class assumes that the data is defined on a
    *rectilinear* grid.

    .. versionadded:: 0.14

    The 'slinear'(k=1), 'cubic'(k=3), and 'quintic'(k=5) methods are
    tensor-product spline interpolators, where `k` is the spline degree,
    If any dimension has fewer points than `k` + 1, an error will be raised.

    .. versionadded:: 1.9

    Examples
    --------
    **Evaluate a function on the points of a 3-D grid**

    As a first example, we evaluate a simple example function on the points of
    a 3-D grid:

    >>> from scipy.interpolate import RegularGridInterpolator
    >>> def f(x, y, z):
    ...     return 2 * x**3 + 3 * y**2 - z
    >>> x = np.linspace(1, 4, 11)
    >>> y = np.linspace(4, 7, 22)
    >>> z = np.linspace(7, 9, 33)
    >>> xg, yg ,zg = np.meshgrid(x, y, z, indexing='ij', sparse=True)
    >>> data = f(xg, yg, zg)

    ``data`` is now a 3-D array with ``data[i, j, k] = f(x[i], y[j], z[k])``.
    Next, define an interpolating function from this data:

    >>> interp = RegularGridInterpolator((x, y, z), data)

    Evaluate the interpolating function at the two points
    ``(x,y,z) = (2.1, 6.2, 8.3)`` and ``(3.3, 5.2, 7.1)``:

    >>> pts = np.array([[2.1, 6.2, 8.3],
    ...                 [3.3, 5.2, 7.1]])
    >>> interp(pts)
    array([ 125.80469388,  146.30069388])

    which is indeed a close approximation to

    >>> f(2.1, 6.2, 8.3), f(3.3, 5.2, 7.1)
    (125.54200000000002, 145.894)

    **Interpolate and extrapolate a 2D dataset**

    As a second example, we interpolate and extrapolate a 2D data set:

    >>> x, y = np.array([-2, 0, 4]), np.array([-2, 0, 2, 5])
    >>> def ff(x, y):
    ...     return x**2 + y**2

    >>> xg, yg = np.meshgrid(x, y, indexing='ij')
    >>> data = ff(xg, yg)
    >>> interp = RegularGridInterpolator((x, y), data,
    ...                                  bounds_error=False, fill_value=None)

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(projection='3d')
    >>> ax.scatter(xg.ravel(), yg.ravel(), data.ravel(),
    ...            s=60, c='k', label='data')

    Evaluate and plot the interpolator on a finer grid

    >>> xx = np.linspace(-4, 9, 31)
    >>> yy = np.linspace(-4, 9, 31)
    >>> X, Y = np.meshgrid(xx, yy, indexing='ij')

    >>> # interpolator
    >>> ax.plot_wireframe(X, Y, interp((X, Y)), rstride=3, cstride=3,
    ...                   alpha=0.4, color='m', label='linear interp')

    >>> # ground truth
    >>> ax.plot_wireframe(X, Y, ff(X, Y), rstride=3, cstride=3,
    ...                   alpha=0.4, label='ground truth')
    >>> plt.legend()
    >>> plt.show()

    Other examples are given
    :ref:`in the tutorial <tutorial-interpolate_regular_grid_interpolator>`.

    See Also
    --------
    NearestNDInterpolator : Nearest neighbor interpolation on *unstructured*
                            data in N dimensions

    LinearNDInterpolator : Piecewise linear interpolant on *unstructured* data
                           in N dimensions

    interpn : a convenience function which wraps `RegularGridInterpolator`

    scipy.ndimage.map_coordinates : interpolation on grids with equal spacing
                                    (suitable for e.g., N-D image resampling)

    References
    ----------
    .. [1] Python package *regulargrid* by Johannes Buchner, see
           https://pypi.python.org/pypi/regulargrid/
    .. [2] Wikipedia, "Trilinear interpolation",
           https://en.wikipedia.org/wiki/Trilinear_interpolation
    .. [3] Weiser, Alan, and Sergio E. Zarantonello. "A note on piecewise linear
           and multilinear table interpolation in many dimensions." MATH.
           COMPUT. 50.181 (1988): 189-196.
           https://www.ams.org/journals/mcom/1988-50-181/S0025-5718-1988-0917826-0/S0025-5718-1988-0917826-0.pdf
           :doi:`10.1090/S0025-5718-1988-0917826-0`

    """
    # this class is based on code originally programmed by Johannes Buchner,
    # see https://github.com/JohannesBuchner/regulargrid

    _SPLINE_DEGREE_MAP = {"slinear": 1, "cubic": 3, "quintic": 5, }
    _SPLINE_METHODS = list(_SPLINE_DEGREE_MAP.keys())
    _ALL_METHODS = ["linear", "nearest"] + _SPLINE_METHODS

    def __init__(self, points, values, method="linear", bounds_error=True,
                 fill_value=np.nan):
        if method not in self._ALL_METHODS:
            raise ValueError("Method '%s' is not defined" % method)
        elif method in self._SPLINE_METHODS:
            self._validate_grid_dimensions(points, method)
        self.method = method
        self.bounds_error = bounds_error

        if not hasattr(values, 'ndim'):
            # allow reasonable duck-typed values
            values = np.asarray(values)

        if len(points) > values.ndim:
            raise ValueError("There are %d point arrays, but values has %d "
                             "dimensions" % (len(points), values.ndim))

        if hasattr(values, 'dtype') and hasattr(values, 'astype'):
            if not np.issubdtype(values.dtype, np.inexact):
                values = values.astype(float)

        self.fill_value = fill_value
        if fill_value is not None:
            fill_value_dtype = np.asarray(fill_value).dtype
            if (hasattr(values, 'dtype') and not
                    np.can_cast(fill_value_dtype, values.dtype,
                                casting='same_kind')):
                raise ValueError("fill_value must be either 'None' or "
                                 "of a type compatible with values")

        for i, p in enumerate(points):
            diff_p = np.diff(p)
            if not np.all(diff_p > 0.):
                if np.all(diff_p < 0.):
                    # input is descending, so make it ascending
                    points, values = _make_points_and_values_ascending(
                        points, values)
                else:
                    raise ValueError(
                        "The points in dimension %d must be strictly "
                        "ascending or descending" % i)
            if not np.asarray(p).ndim == 1:
                raise ValueError("The points in dimension %d must be "
                                 "1-dimensional" % i)
            if not values.shape[i] == len(p):
                raise ValueError("There are %d points and %d values in "
                                 "dimension %d" % (len(p), values.shape[i], i))
        self.grid = tuple([np.asarray(p) for p in points])
        self.values = values

    def __call__(self, xi, method=None):
        """
        Interpolation at coordinates.

        Parameters
        ----------
        xi : ndarray of shape (..., ndim)
            The coordinates to evaluate the interpolator at.

        method : str
            The method of interpolation to perform. Supported are "linear" and
            "nearest".

        Examples
        --------
        Here we define a nearest-neighbor interpolator of a simple function

        >>> x, y = np.array([0, 1, 2]), np.array([1, 3, 7])
        >>> def f(x, y):
        ...     return x**2 + y**2
        >>> data = f(*np.meshgrid(x, y, indexing='ij', sparse=True))
        >>> from scipy.interpolate import RegularGridInterpolator
        >>> interp = RegularGridInterpolator((x, y), data, method='nearest')

        By construction, the interpolator uses the nearest-neighbor
        interpolation

        >>> interp([[1.5, 1.3], [0.3, 4.5]])
        array([2., 9.])

        We can however evaluate the linear interpolant by overriding the
        `method` parameter

        >>> interp([[1.5, 1.3], [0.3, 4.5]], method='linear')
        array([ 4.7, 24.3])
        """
        is_method_changed = self.method != method
        method = self.method if method is None else method
        if method not in self._ALL_METHODS:
            raise ValueError("Method '%s' is not defined" % method)

        ndim = len(self.grid)
        xi = _ndim_coords_from_arrays(xi, ndim=ndim)
        if xi.shape[-1] != len(self.grid):
            raise ValueError("The requested sample points xi have dimension "
                             "%d, but this RegularGridInterpolator has "
                             "dimension %d" % (xi.shape[1], ndim))

        xi_shape = xi.shape
        xi = xi.reshape(-1, xi_shape[-1])

        # find nans in input
        nans = np.any(np.isnan(xi), axis=-1)

        if self.bounds_error:
            for i, p in enumerate(xi.T):
                if not np.logical_and(np.all(self.grid[i][0] <= p),
                                      np.all(p <= self.grid[i][-1])):
                    raise ValueError("One of the requested xi is out of bounds "
                                     "in dimension %d" % i)

        indices, norm_distances, out_of_bounds = self._find_indices(xi.T)
        if method == "linear":
            result = self._evaluate_linear(indices,
                                           norm_distances,
                                           out_of_bounds)
        elif method == "nearest":
            result = self._evaluate_nearest(indices,
                                            norm_distances,
                                            out_of_bounds)
        elif method in self._SPLINE_METHODS:
            if is_method_changed:
                self._validate_grid_dimensions(self.grid, method)
            result = self._evaluate_spline(self.values.T, xi,
                                           self._SPLINE_DEGREE_MAP[method])

        if not self.bounds_error and self.fill_value is not None:
            result[out_of_bounds] = self.fill_value

        # f(nan) = nan, if any
        if np.any(nans):
            result[nans] = np.nan
        return result.reshape(xi_shape[:-1] + self.values.shape[ndim:])

    def _evaluate_linear(self, indices, norm_distances, out_of_bounds):
        # slice for broadcasting over trailing dimensions in self.values
        vslice = (slice(None),) + (None,)*(self.values.ndim - len(indices))

        # find relevant values
        # each i and i+1 represents a edge
        edges = itertools.product(*[[i, i + 1] for i in indices])
        values = 0.
        for edge_indices in edges:
            weight = 1.
            for ei, i, yi in zip(edge_indices, indices, norm_distances):
                weight *= np.where(ei == i, 1 - yi, yi)
            values += np.asarray(self.values[edge_indices]) * weight[vslice]
        return values

    def _evaluate_nearest(self, indices, norm_distances, out_of_bounds):
        idx_res = [np.where(yi <= .5, i, i + 1)
                   for i, yi in zip(indices, norm_distances)]
        return self.values[tuple(idx_res)]

    def _validate_grid_dimensions(self, points, method):
        k = self._SPLINE_DEGREE_MAP[method]
        for i, point in enumerate(points):
            ndim = len(np.atleast_1d(point))
            if ndim <= k:
                raise ValueError(f"There are {ndim} points in dimension {i},"
                                 f" but method {method} requires at least "
                                 f" {k+1} points per dimension.")

    def _evaluate_spline(self, values, xi, spline_degree):
        # ensure xi is 2D list of points to evaluate
        if xi.ndim == 1:
            xi = xi.reshape((1, xi.size))
        m, n = xi.shape

        # Non-stationary procedure: difficult to vectorize this part entirely
        # into numpy-level operations. Unfortunately this requires explicit
        # looping over each point in xi.

        # can at least vectorize the first pass across all points in the
        # last variable of xi.
        last_dim = n - 1
        first_values = self._do_spline_fit(self.grid[last_dim],
                                           values,
                                           xi[:, last_dim],
                                           spline_degree)

        # the rest of the dimensions have to be on a per point-in-xi basis
        result = np.empty(m, dtype=self.values.dtype)
        for j in range(m):
            # Main process: Apply 1D interpolate in each dimension
            # sequentially, starting with the last dimension.
            # These are then "folded" into the next dimension in-place.
            folded_values = first_values[j]
            for i in range(last_dim-1, -1, -1):
                # Interpolate for each 1D from the last dimensions.
                # This collapses each 1D sequence into a scalar.
                folded_values = self._do_spline_fit(self.grid[i],
                                                    folded_values,
                                                    xi[j, i],
                                                    spline_degree)

            result[j] = folded_values

        return result

    @staticmethod
    def _do_spline_fit(x, y, pt, k):
        local_interp = make_interp_spline(x, y, k=k, axis=0)
        values = local_interp(pt)
        return values

    def _find_indices(self, xi):
        # find relevant edges between which xi are situated
        indices = []
        # compute distance to lower edge in unity units
        norm_distances = []
        # check for out of bounds xi
        out_of_bounds = np.zeros((xi.shape[1]), dtype=bool)
        # iterate through dimensions
        for x, grid in zip(xi, self.grid):
            i = np.searchsorted(grid, x) - 1
            i[i < 0] = 0
            i[i > grid.size - 2] = grid.size - 2
            indices.append(i)

            # compute norm_distances, incl length-1 grids,
            # where `grid[i+1] == grid[i]`
            denom = grid[i + 1] - grid[i]
            with np.errstate(divide='ignore', invalid='ignore'):
                norm_dist = np.where(denom != 0, (x - grid[i]) / denom, 0)
            norm_distances.append(norm_dist)

            if not self.bounds_error:
                out_of_bounds += x < grid[0]
                out_of_bounds += x > grid[-1]
        return indices, norm_distances, out_of_bounds


def interpn(points, values, xi, method="linear", bounds_error=True,
            fill_value=np.nan):
    """
    Multidimensional interpolation on regular or rectilinear grids.

    Strictly speaking, not all regular grids are supported - this function
    works on *rectilinear* grids, that is, a rectangular grid with even or
    uneven spacing.

    Parameters
    ----------
    points : tuple of ndarray of float, with shapes (m1, ), ..., (mn, )
        The points defining the regular grid in n dimensions. The points in
        each dimension (i.e. every elements of the points tuple) must be
        strictly ascending or descending.

    values : array_like, shape (m1, ..., mn, ...)
        The data on the regular grid in n dimensions. Complex data can be
        acceptable.

    xi : ndarray of shape (..., ndim)
        The coordinates to sample the gridded data at

    method : str, optional
        The method of interpolation to perform. Supported are "linear" and
        "nearest", and "splinef2d". "splinef2d" is only supported for
        2-dimensional data.

    bounds_error : bool, optional
        If True, when interpolated values are requested outside of the
        domain of the input data, a ValueError is raised.
        If False, then `fill_value` is used.

    fill_value : number, optional
        If provided, the value to use for points outside of the
        interpolation domain. If None, values outside
        the domain are extrapolated.  Extrapolation is not supported by method
        "splinef2d".

    Returns
    -------
    values_x : ndarray, shape xi.shape[:-1] + values.shape[ndim:]
        Interpolated values at input coordinates.

    Notes
    -----

    .. versionadded:: 0.14

    Examples
    --------
    Evaluate a simple example function on the points of a regular 3-D grid:

    >>> from scipy.interpolate import interpn
    >>> def value_func_3d(x, y, z):
    ...     return 2 * x + 3 * y - z
    >>> x = np.linspace(0, 4, 5)
    >>> y = np.linspace(0, 5, 6)
    >>> z = np.linspace(0, 6, 7)
    >>> points = (x, y, z)
    >>> values = value_func_3d(*np.meshgrid(*points, indexing='ij'))

    Evaluate the interpolating function at a point

    >>> point = np.array([2.21, 3.12, 1.15])
    >>> print(interpn(points, values, point))
    [12.63]

    See Also
    --------
    NearestNDInterpolator : Nearest neighbor interpolation on unstructured
                            data in N dimensions

    LinearNDInterpolator : Piecewise linear interpolant on unstructured data
                           in N dimensions

    RegularGridInterpolator : interpolation on a regular or rectilinear grid
                              in arbitrary dimensions (`interpn` wraps this
                              class).

    RectBivariateSpline : Bivariate spline approximation over a rectangular mesh

    scipy.ndimage.map_coordinates : interpolation on grids with equal spacing
                                    (suitable for e.g., N-D image resampling)

    """
    # sanity check 'method' kwarg
    if method not in ["linear", "nearest", "splinef2d"]:
        raise ValueError("interpn only understands the methods 'linear', "
                         "'nearest', and 'splinef2d'. You provided %s." %
                         method)

    if not hasattr(values, 'ndim'):
        values = np.asarray(values)

    ndim = values.ndim
    if ndim > 2 and method == "splinef2d":
        raise ValueError("The method splinef2d can only be used for "
                         "2-dimensional input data")
    if not bounds_error and fill_value is None and method == "splinef2d":
        raise ValueError("The method splinef2d does not support extrapolation.")

    # sanity check consistency of input dimensions
    if len(points) > ndim:
        raise ValueError("There are %d point arrays, but values has %d "
                         "dimensions" % (len(points), ndim))
    if len(points) != ndim and method == 'splinef2d':
        raise ValueError("The method splinef2d can only be used for "
                         "scalar data with one point per coordinate")

    # sanity check input grid
    for i, p in enumerate(points):
        diff_p = np.diff(p)
        if not np.all(diff_p > 0.):
            if np.all(diff_p < 0.):
                # input is descending, so make it ascending
                points, values = _make_points_and_values_ascending(points,
                                                                   values)
            else:
                raise ValueError("The points in dimension %d must be strictly "
                                 "ascending or descending" % i)
        if not np.asarray(p).ndim == 1:
            raise ValueError("The points in dimension %d must be "
                             "1-dimensional" % i)
        if not values.shape[i] == len(p):
            raise ValueError("There are %d points and %d values in "
                             "dimension %d" % (len(p), values.shape[i], i))
    grid = tuple([np.asarray(p) for p in points])

    # sanity check requested xi
    xi = _ndim_coords_from_arrays(xi, ndim=len(grid))
    if xi.shape[-1] != len(grid):
        raise ValueError("The requested sample points xi have dimension "
                         "%d, but this RegularGridInterpolator has "
                         "dimension %d" % (xi.shape[1], len(grid)))

    if bounds_error:
        for i, p in enumerate(xi.T):
            if not np.logical_and(np.all(grid[i][0] <= p),
                                                np.all(p <= grid[i][-1])):
                raise ValueError("One of the requested xi is out of bounds "
                                "in dimension %d" % i)

    # perform interpolation
    if method == "linear":
        interp = RegularGridInterpolator(points, values, method="linear",
                                         bounds_error=bounds_error,
                                         fill_value=fill_value)
        return interp(xi)
    elif method == "nearest":
        interp = RegularGridInterpolator(points, values, method="nearest",
                                         bounds_error=bounds_error,
                                         fill_value=fill_value)
        return interp(xi)
    elif method == "splinef2d":
        xi_shape = xi.shape
        xi = xi.reshape(-1, xi.shape[-1])

        # RectBivariateSpline doesn't support fill_value; we need to wrap here
        idx_valid = np.all((grid[0][0] <= xi[:, 0], xi[:, 0] <= grid[0][-1],
                            grid[1][0] <= xi[:, 1], xi[:, 1] <= grid[1][-1]),
                           axis=0)
        result = np.empty_like(xi[:, 0])

        # make a copy of values for RectBivariateSpline
        interp = RectBivariateSpline(points[0], points[1], values[:])
        result[idx_valid] = interp.ev(xi[idx_valid, 0], xi[idx_valid, 1])
        result[np.logical_not(idx_valid)] = fill_value

        return result.reshape(xi_shape[:-1])
