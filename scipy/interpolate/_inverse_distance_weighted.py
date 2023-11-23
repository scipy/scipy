import numpy as np
from scipy.spatial import cKDTree
from .interpnd import NDInterpolatorBase, _ndim_coords_from_arrays


class InverseDistanceWeightedNDInterpolator(NDInterpolatorBase):
    def __init__(self, x, y, rescale=False, fill_value=np.nan, tree_options=None, local=True):
        NDInterpolatorBase.__init__(self, x, y, rescale=rescale,
                                    need_contiguous=False,
                                    need_values=False)
        if tree_options is None:
            tree_options = dict()
        self.tree = cKDTree(self.points, **tree_options)
        self.values = np.asarray(y)

        self._local = local
        self._fill_value = fill_value

    def __call__(self, xi, p=2, k=6, distance_upper_bound=np.inf, **query_options):
        xi = _ndim_coords_from_arrays(xi, ndim=self.points.shape[1])
        xi = self._check_call_shape(xi)
        xi = self._scale_x(xi)

        if self._local:
            interp_values = self._local_idw_interpolation(xi, p=p, k=k, distance_upper_bound=distance_upper_bound,
                                                          **query_options)
        else:
            interp_values = self._gloabl_idw_interpolation(xi, p=p)

        return interp_values

    def _local_idw_interpolation(self, xi, p, k, distance_upper_bound, **query_options):
        xi_flat = xi.reshape(-1, xi.shape[-1])
        original_shape = xi.shape
        flattened_shape = xi_flat.shape

        dist, i = self.tree.query(xi_flat, k=k, distance_upper_bound=distance_upper_bound, **query_options)
        valid_mask = np.isfinite(dist)

        weights = 1.0 / dist ** p

        # create a holder interp_values array with shape (n*m*.., k, l, ...) and fill with nans.
        interp_shape = flattened_shape[:-1] + self.values.shape[1:] if self.values.ndim > 1 else flattened_shape[:-1]

        if np.issubdtype(self.values.dtype, np.complexfloating):
            interp_values = np.full(interp_shape, self._fill_value, dtype=self.values.dtype)
        else:
            interp_values = np.full(interp_shape, self._fill_value)

        if self.values.ndim == 1:
            interp_values[valid_mask] = np.sum(weights[valid_mask, ...] * self.values[i[valid_mask]], axis=-1) / np.sum(
                weights[valid_mask], axis=-1)
        else:
            interp_values[valid_mask] = np.sum(weights[valid_mask, ...] * self.values[i[valid_mask], ...], axis=-1) / np.sum(
                weights[valid_mask], axis=-1)

        new_shape = original_shape[:-1] + self.values.shape[1:] if self.values.ndim > 1 else original_shape[:-1]
        interp_values = interp_values.reshape(new_shape)

        return interp_values
