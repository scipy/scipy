import numpy as np

from .common import Benchmark

try:
    from scipy.ndimage import uniform_filter
except ImportError:
    pass


class NdimageUniformFilter(Benchmark):
    param_names = ['shape', 'size', 'mode']
    params = [
        [(512, 512), (2048, 2048), (4096, 4096), (128, 128, 128), (512, 512, 512)],
        [1, 3, 5, 7, 9],
        ['reflect', 'constant', 'nearest', 'mirror', 'wrap']
    ]

    def setup(self, shape, order, mode):
        rstate = np.random.RandomState(5)
        self.x = rstate.uniform(shape)

    def time_uniform_filter(self, shape, size, mode):
        uniform_filter(self.x, size=size, mode=mode)

    def time_uniform_filter_nan(self, shape, size, mode):
        uniform_filter(self.x, size=size, mode=mode, nan_policy='omit')

