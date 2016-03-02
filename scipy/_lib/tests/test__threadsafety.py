from __future__ import division, print_function, absolute_import

import threading
import time

from numpy.testing import assert_equal, assert_raises

from scipy._lib._threadsafety import ReentrancyLock, non_reentrant, ReentrancyError


def test_parallel_threads():
    results = []

    lock = ReentrancyLock("failure")

    def worker(k):
        with lock:
            time.sleep(0.01*(3 - k))
            results.append(k)

    threads = [threading.Thread(target=lambda k=k: worker(k))
               for k in range(3)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert_equal(results, sorted(results))


def test_reentering():
    @non_reentrant()
    def func(x):
        return func(x)

    assert_raises(ReentrancyError, func, 0)
