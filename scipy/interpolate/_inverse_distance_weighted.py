import numpy as np
from scipy.spatial import cKDTree, distance
from scipy.interpolate.interpnd import NDInterpolatorBase, _ndim_coords_from_arrays


class InverseDistanceWeightedNDInterpolator(NDInterpolatorBase):
    """
    Inverse Distance Weighted (IDW) interpolation.

    Parameters
    ----------
    x : (npoints, ndims) ndarray of floats
        Data point coordinates.
    y : (npoints, ) ndarray of floats
        Data values associated with the points.
    rescale : bool, optional
        Rescale points to a unit cube before performing interpolation.
        This is useful if some of the input dimensions have incommensurable units
        and differ by many orders of magnitude. Default is False.
    fill_value : float or nan, optional
        Value used to fill in for requested points outside of the convex hull
        of the input points. Default is nan.
    tree_options : dict, optional
        Options passed to the underlying `cKDTree`.
    local : bool, optional
        If True, perform local interpolation, else global. Default is True.

    Methods
    -------
    __call__(xi, p=2, k=3, distance_upper_bound=np.inf, **query_options)
        Evaluate the interpolator at given points.

    Examples
    --------
    >>> points = np.array([[0, 0], [1, 1], [2, 2]])
    >>> values = np.array([0, 1, 2])
    >>> interpolator = InverseDistanceWeightedNDInterpolator(points, values)
    >>> interpolator([1.5, 1.5])
    array([1.42105263])
    >>> def exp_weight(dists, power, decorrelation_length=1.5):
    ...     return np.exp(-(dists / decorrelation_length) ** p)
    >>> interpolator([1.5, 1.5], weight_func=exp_weight)
    array([1.38311737])

    Notes
    -----
    The interpolation uses an inverse distance weighting method, where the
    influence of each data point on the interpolated values is inversely
    proportional to the distance from the point to a power `p` or specified
    by the user.
    """

    def __init__(
        self, x, y, rescale=False, fill_value=np.nan, tree_options=None, local=True
    ):
        NDInterpolatorBase.__init__(
            self, x, y, rescale=rescale, need_contiguous=False, need_values=False
        )
        if tree_options is None:
            tree_options = dict()
        self.tree = cKDTree(self.points, **tree_options)
        self.values = np.asarray(y)

        self._local = local
        self._fill_value = fill_value

    def __call__(
        self,
        xi,
        weight_func=None,
        p=2,
        k=3,
        distance_upper_bound=np.inf,
        **query_options,
    ):
        """
        Evaluate the interpolator at given points.

        The interpolating weights are given by:

            w_i = weight_func(d_i, p)

        The default weighting function is the inverse function:

            w_i = 1/d_i^p

        Parameters
        ----------
        xi : array_like
            Points to interpolate the data at.
        weight_func : callable, optional
            Manual way to determine weights. Default is 1/dist^p.
        p : float, optional
            The power of the inverse distance weighting. Default is 2.
        k : int, optional
            The number of nearest neighbors to use for interpolation. Default is 3.
        distance_upper_bound : float, optional
            The maximum distance to consider for points to influence the
            interpolation at a given point. Default is np.inf.

        **query_options : dict, optional
            Additional keyword arguments to pass to the cKDTree query method.

        Returns
        -------
        ndarray
            Interpolated values at the input points ``xi``.
        """
        xi = _ndim_coords_from_arrays(xi, ndim=self.points.shape[1])
        xi = self._check_call_shape(xi)
        xi = self._scale_x(xi)

        if self._local:
            interp_values = self._local_idw_interpolation(
                xi,
                weight_func=weight_func,
                p=p,
                k=k,
                distance_upper_bound=distance_upper_bound,
                **query_options,
            )
        else:
            interp_values = self._global_idw_interpolation(
                xi, weight_func=weight_func, p=p
            )

        return interp_values

    def _local_idw_interpolation(
        self, xi, weight_func, p, k, distance_upper_bound, eps=1e-7, **query_options
    ):
        """
        Perform local Inverse Distance Weighting interpolation.

        In local IDW, only the k-nearest neighbors of each query point are considered
        for interpolation. This method tends to be faster, but requires some careful
        thought on hyperparameters chosen for the IDW interpolation.

        Parameters
        ----------
        xi : ndarray
            The query points where interpolation is performed.
        p : float
            The power parameter of IDW. Higher values assign greater weight to closer
            points.
        k : int
            The number of nearest neighbors to consider for each query point.
        distance_upper_bound : float
            Maximum distance for points to be considered as neighbors.
        eps : float, optional
            Small constant to prevent division by zero in weight calculation. This can
            be changed depending on the distance scale of different problems.
        **query_options : dict, optional
            Additional options to pass to the KDTree query.

        Returns
        -------
        ndarray
            The interpolated values at the query points.
        """
        xi_flat = xi.reshape(-1, xi.shape[-1])
        original_shape = xi.shape
        flattened_shape = xi_flat.shape

        dist, i = self.tree.query(
            xi_flat, k=k, distance_upper_bound=distance_upper_bound, **query_options
        )

        # for points that lie outside of distance_upper_bound, returned dist from
        # cKDTree will be inf
        valid_mask = np.isfinite(dist)
        to_interpolate_mask = np.any(valid_mask, axis=1)

        if weight_func is not None:
            weights = weight_func(dist, p)
        else:
            weights = 1.0 / np.maximum(dist, eps) ** p

        # create correct shape for the returned interpolated values and for the
        # intermediate relevant_values and weighted_values
        interp_shape = (
            flattened_shape[:-1] + self.values.shape[1:]
            if self.values.ndim > 1
            else flattened_shape[:-1]
        )

        if np.issubdtype(self.values.dtype, np.complexfloating):
            relevant_values = np.full(
                interp_shape + (k,), np.nan, dtype=self.values.dtype
            )
            weighted_values = np.full(interp_shape + (k,), 0.0, dtype=self.values.dtype)
            interp_values = np.full(
                interp_shape, self._fill_value, dtype=self.values.dtype
            )
        else:
            relevant_values = np.full(interp_shape + (k,), np.nan)
            weighted_values = np.full(interp_shape + (k,), 0.0)
            interp_values = np.full(interp_shape, self._fill_value)

        # fetch value of valid nearest neighbors and multiply by weights
        relevant_values[valid_mask] = self.values[i[valid_mask]]
        np.multiply(weights, relevant_values, where=valid_mask, out=weighted_values)

        # sum(x*w)/sum(w)
        numerator = np.sum(weighted_values, axis=-1)
        denominator = np.sum(weights, where=valid_mask, axis=-1)

        # only set values for those points that had real data points within
        # distance_upper_bound
        interp_values[to_interpolate_mask] = (
            numerator[to_interpolate_mask] / denominator[to_interpolate_mask]
        )

        # handle case where there are any/mutiple equidistant points

        zero_dist_mask_points = np.any(dist == 0, axis=1)
        if True in zero_dist_mask_points:
            zero_dist_mask_full = dist == 0
            values = np.mean(self.values[i], where=zero_dist_mask_full, axis=1)
            # clunky handle of case when numpy still wants to take the mean for points
            # with all False (i.e. no points with 0 dists)
            interp_values[zero_dist_mask_points] = values[~np.isnan(values)]

        new_shape = (
            original_shape[:-1] + self.values.shape[1:]
            if self.values.ndim > 1
            else original_shape[:-1]
        )
        interp_values = interp_values.reshape(new_shape)

        return interp_values

    def _global_idw_interpolation(self, xi, weight_func=None, p=2, eps=1e-7):
        """
        Perform global Inverse Distance Weighting interpolation.

        In global IDW, all points in the dataset contribute to the interpolation of each
        query point.

        Parameters
        ----------
        xi : ndarray
            The query points where interpolation is performed.
        p : float
            The power parameter of IDW.
        eps : float, optional
            Small constant to prevent division by zero in weight calculation.

        Returns
        -------
        ndarray
            The interpolated values at the query points.
        """
        xi_flat = xi.reshape(-1, xi.shape[-1])
        original_shape = xi.shape
        flattened_shape = xi_flat.shape

        interp_shape = (
            flattened_shape[:-1] + self.values.shape[1:]
            if self.values.ndim > 1
            else flattened_shape[:-1]
        )

        if np.issubdtype(self.values.dtype, np.complexfloating):
            interp_values = np.full(
                interp_shape, self._fill_value, dtype=self.values.dtype
            )
        else:
            interp_values = np.full(interp_shape, self._fill_value)

        for i in range(len(xi_flat)):
            dist = distance.cdist(xi_flat[i, np.newaxis], self.points).ravel()

            if weight_func is not None:
                weights = weight_func(dist, p)
            else:
                weights = 1.0 / np.maximum(dist, eps) ** p

            zero_dist_mask = dist == 0

            if True in zero_dist_mask:
                interp_values[i] = np.mean(self.values[zero_dist_mask])

            else:
                if np.any(np.isfinite(dist)):
                    numerator = np.sum(weights * self.values)
                    denominator = np.sum(weights)
                    interp_values[i] = numerator / denominator

        new_shape = (
            original_shape[:-1] + self.values.shape[1:]
            if self.values.ndim > 1
            else original_shape[:-1]
        )
        interp_values = interp_values.reshape(new_shape)

        return interp_values
