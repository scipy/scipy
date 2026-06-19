import zlib

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


_ZOOM_CASES = [
    ('2d_cubic_up_sinusoid_f32',
     (160, 192), (1.7, 1.25), 3, 'float32', 'sinusoid'),
    ('2d_linear_aniso_random_f32',
     (384, 320), (1.0, 0.37), 1, 'float32', 'random'),
    ('3d_cubic_down_random_f32',
     (64, 56, 24), (0.5, 0.5, 0.5), 3, 'float32', 'random'),
    ('3d_linear_down_random_f32',
     (80, 72, 28), (0.5, 0.5, 0.5), 1, 'float32', 'random'),
    ('2d_linear_up_sinusoid_f64',
     (192, 224), (1.25, 1.25), 1, 'float64', 'sinusoid'),
    ('2d_linear_down_random_f32',
     (192, 160), (0.5, 0.5), 1, 'float32', 'random'),
    ('2d_cubic_down_random_f32',
     (160, 192), (0.5, 0.5), 3, 'float32', 'random'),
    ('2d_cubic_down_random_f64',
     (256, 256), (0.37, 0.37), 3, 'float64', 'random'),
]


def _zoom_seed(name, shape, zoom_factor, order, dtype_name, pattern):
    payload = f"{name}|{shape}|{zoom_factor}|{order}|{dtype_name}|{pattern}"
    return zlib.crc32(payload.encode("utf-8")) & 0xFFFFFFFF


def _make_zoom_input(name, shape, zoom_factor, order, dtype_name, pattern):
    dtype = np.dtype(dtype_name)
    rng = np.random.default_rng(
        _zoom_seed(name, shape, zoom_factor, order, dtype_name, pattern)
    )
    if pattern == 'random':
        return rng.random(shape, dtype=dtype)
    if pattern == 'sinusoid':
        grids = np.meshgrid(
            *[np.arange(n, dtype=np.float64) for n in shape],
            indexing='ij',
        )
        out = np.zeros(shape, dtype=np.float64)
        for axis, grid in enumerate(grids):
            out += np.sin(2.0 * np.pi * (0.07 + 0.05 * axis) * grid)
        return (out / float(len(grids))).astype(dtype, copy=False)
    raise NotImplementedError(pattern)


class NdimageZoom(Benchmark):
    param_names = ['case']
    params = [_ZOOM_CASES]

    def setup(self, case):
        name, shape, zoom_factor, order, dtype_name, pattern = case
        self.x = _make_zoom_input(
            name, shape, zoom_factor, order, dtype_name, pattern
        )
        self.zoom_factor = zoom_factor
        self.order = int(order)
        self.prefilter = self.order > 1

    def time_zoom_mirror_grid_false(self, case):
        zoom(
            self.x,
            self.zoom_factor,
            order=self.order,
            mode='mirror',
            prefilter=self.prefilter,
            grid_mode=False,
        )
