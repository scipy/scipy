import sys
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
    result = module.test(mode)
    t_stop = time.time()
    return (result, t_stop - t_start)


timings = dict()
tests_passed = True
for module in modules:
    result, total_time = test_module(module, mode='full')
    timings[module.__name__] = total_time
    tests_passed = tests_passed and result.wasSuccessful()


for key in reversed(sorted(timings, key=lambda key: timings[key])):
    print("%25s: %1.3f" % (key, timings[key]))


if tests_passed:
    sys.exit(0)
else:
    sys.exit(1)
