"""
Functions copied from newer versions of NumPy.
"""

import numpy as np
from ._pep440 import Version


if Version(np.__version__) >= Version('1.15.0'):
    from numpy import take_along_axis
else:
    def _make_along_axis_idx(arr_shape, indices, axis):
        import numpy.core.numeric as _nx

        # compute dimensions to iterate over
        if not _nx.issubdtype(indices.dtype, _nx.integer):
            raise IndexError('`indices` must be an integer array')
        if len(arr_shape) != indices.ndim:
            raise ValueError(
                "`indices` and `arr` must have the same number of dimensions")
        shape_ones = (1,) * indices.ndim
        dest_dims = list(range(axis)) + [None] + list(range(axis+1,
                                                            indices.ndim))

        # build a fancy index, consisting of orthogonal aranges, with the
        # requested index inserted at the right location
        fancy_index = []
        for dim, n in zip(dest_dims, arr_shape):
            if dim is None:
                fancy_index.append(indices)
            else:
                ind_shape = shape_ones[:dim] + (-1,) + shape_ones[dim+1:]
                fancy_index.append(_nx.arange(n).reshape(ind_shape))

        return tuple(fancy_index)

    def take_along_axis(arr, indices, axis):
        from numpy.core.multiarray import normalize_axis_index

        # normalize inputs
        if axis is None:
            arr = arr.flat
            arr_shape = (len(arr),)  # flatiter has no .shape
            axis = 0
        else:
            axis = normalize_axis_index(axis, arr.ndim)
            arr_shape = arr.shape

        # use the fancy index
        return arr[_make_along_axis_idx(arr_shape, indices, axis)]
