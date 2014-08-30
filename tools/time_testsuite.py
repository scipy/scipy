import time

from scipy import (cluster, constants, fftpack, integrate, interpolate, io,
        linalg, misc, ndimage, odr, optimize, signal, sparse,
        spatial, special, stats, weave)


modules = [cluster, constants, fftpack, integrate, interpolate, io,
           linalg, misc, ndimage, odr, optimize, signal, sparse,
           spatial, special, stats, weave]


def test_module(module, mode='fast'):
    """Run tests for a module, then return how much time they took."""
    t_start = time.time()
    module.test(mode)
    t_stop = time.time()
    return t_stop - t_start


timings = dict()
for module in modules:
    timings[module.__name__] = test_module(module, mode='full')


for key in reversed(sorted(timings, key=lambda key: timings[key])):
    print("%25s: %1.3f" % (key, timings[key]))
