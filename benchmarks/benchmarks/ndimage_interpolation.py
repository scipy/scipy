import numpy as np

from .common import Benchmark

try:
    from scipy.ndimage import (geometric_transform, affine_transform, rotate,
                               zoom, shift, map_coordinates)
except ImportError:
    pass


def shift_func_2d(c):
    return (c[0] - 0.5, c[1] - 0.5)


def shift_func_3d(c):
    return (c[0] - 0.5, c[1] - 0.5, c[2] - 0.5)


class NdimageInterpolation(Benchmark):
    param_names = ['shape', 'order', 'mode']
    params = [
        [(64, 64), (512, 512), (2048, 2048), (16, 16, 16), (128, 128, 128)],
        [0, 1, 3, 5],
        ['mirror', 'constant']
    ]

    def setup(self, shape, order, mode):
        rstate = np.random.RandomState(5)
        self.x = rstate.standard_normal(shape)
        self.matrix_2d = np.asarray([[0.8, 0, 1.5],
                                     [0, 1.2, -5.]])
        self.matrix_3d = np.asarray([[0.8, 0, 0, 1.5],
                                     [0, 1.2, 0, -5.],
                                     [0, 0, 1, 0]])

    def time_affine_transform(self, shape, order, mode):
        if self.x.ndim == 2:
            matrix = self.matrix_2d
        else:
            matrix = self.matrix_3d
        affine_transform(self.x, matrix, order=order, mode=mode)

    def time_rotate(self, shape, order, mode):
        rotate(self.x, 15, order=order, mode=mode)

    def time_shift(self, shape, order, mode):
        shift(self.x, (-2.5,) * self.x.ndim, order=order, mode=mode)

    def time_zoom(self, shape, order, mode):
        zoom(self.x, (1.3,) * self.x.ndim, order=order, mode=mode)

    def time_geometric_transform_mapping(self, shape, order, mode):
        if self.x.ndim == 2:
            mapping = shift_func_2d
        if self.x.ndim == 3:
            mapping = shift_func_3d
        geometric_transform(self.x, mapping, order=order, mode=mode)

    def time_map_coordinates(self, shape, order, mode):
        coords = np.meshgrid(*[np.arange(0, s, 2) + 0.3 for s in self.x.shape])
        map_coordinates(self.x, coords, order=order, mode=mode)

    def peakmem_rotate(self, shape, order, mode):
        rotate(self.x, 15, order=order, mode=mode)

    def peakmem_shift(self, shape, order, mode):
        shift(self.x, 3, order=order, mode=mode)
