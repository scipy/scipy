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
        self._k = k
        self._local = local

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
        dist, i = self.tree.query(xi, k=self._k, distance_upper_bound=distance_upper_bound, **query_options)

        weights = 1.0 / dist ** p

        mask = np.isfinite(dist)
        ret_array = np.ones(i.shape)
        ret_array[~mask] = np.nan
        ret_array[mask] = weights[mask] * self.values[i[mask]]
        ret_array = np.sum(ret_array, axis=2) / np.sum(weights, axis=2)

        return ret_array
